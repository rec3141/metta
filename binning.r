library(data.table) #fread
library(seqinr) #read.fasta, uco
library(parallel) #mclapply
source("bin/FastPCA.R") #FastPCA
library(Rtsne) #Rtsne
library(Morpho)
library(igraph)
library(scales)
library(gplots)
library(dendextend) #cutree
library(rgl)
library(parallelDist)
library(plyr)

gcut <- function(data,split,range=y) {do.call(paste, data.frame(do.call(rbind, lapply(strsplit(data,split,fixed=T), "[", range))))}

runsetup = F
if(runsetup==T) {

#read in annotation files and depths
dt.in <- fread("annotation/kaiju/scaffolds_gt2000.node_orf_bin_ann_depths_tax.tsv",sep="\t",head=F,quote="", stringsAsFactors=F)
head.in <- fread("annotation/metabat/spades_all_norm_ASGARD_20181114_2018-12-04-15-09-06/ASGARD_20181114/scaffolds_gt2000.fasta.depth.tx",nrows=0)
colnames(dt.in) <- c("node","orf","metabat.bin","prokka.ann",colnames(head.in)[2:ncol(head.in)])
#rownames(dt.in) <- dt.in$orf

#make metadata table
dt.meta <- dt.in[,1:6,with=F]
df.meta <- as.data.frame(dt.meta)

#make numeric tables with coverages/depths and variances
dt.num <- dt.in[,7:ncol(dt.in),with=F]
dt.var <- dt.num[,grep(".var",colnames(dt.num)),with=F]
dt.num <- dt.num[,grep(".var",colnames(dt.num),invert=T),with=F]
df.num <- as.data.frame(dt.num)

# attempt a normal distribution
# dt.num <- log2(d.num+0.0001)
# dt.num <- dt.num - log2(0.0001)

rownames(df.meta) <- dt.meta$orf
rownames(dt.meta) <- dt.meta$orf
rownames(dt.var) <- dt.meta$orf
rownames(dt.num) <- dt.meta$orf
rownames(df.num) <- dt.meta$orf

rm(dt.in)
gc()

#read scaffold paths
paths.in <- read.table("binning/spades_all_2018-10-03-10-44-21.paths.tsv",sep="\t")
gpaths <- graph_from_data_frame(paths.in, directed = FALSE, vertices = NULL)

# read in scaffold taxonomy from kaiju for coloring

# read in scaffold taxonomy from kaiju for PCA
tax.in <- read.table("annotation/kaiju/scaffolds.1000bpp.kaiju_names_all.tsv",head=F,stringsAsFactors=F,sep="\t",quote="",fill=T)
tax.class <- tax.in[tax.in[,1] == "C",]
tax.class[tax.class[,4]=="",4] <- "unknown;"
tax.class[tax.class[,4]==" ",4] <- "unknown;"

zz <- textConnection(paste(as.character(tax.class[,4]),collapse="\n"))
tax.tax <- read.table(zz,head=F,stringsAsFactors=F,sep=";",quote="",fill=T,comment.char="",strip.white=T)
tax.all <- cbind(tax.class[,1:3],tax.tax)
colnames(tax.all) <- c("classified","orf","taxid","na","superkingdom","phylum","class","order","family","genus","species","na2")
#fill in missing
tax.out <- apply(tax.all,1,function(x) {
	saved <- "unknown";
	for(j in 5:11) {
		if(is.na(x[j])) {
			x[j] <- saved
		} else {
			saved <- x[j]
		}
	}
	x
})

tax.out <- as.data.frame(t(tax.out),stringsAsFactors=F)
rownames(tax.out) <- tax.out$orf

tax.out$node <- df.meta[match(tax.out$orf,df.meta$orf),"node"]

#calculate distances on unique partitions (columns), then add them together by membership

### dimensional reduction on separate matrices

# make graph PCA and taxonomy PCA
# balance among them = machine learning

# need to aggregate over orfs and just look at scaffolds
cov.node <- as.data.frame( df.num[order(df.meta$node,df.meta$orf),])
cov.meta <- as.data.frame(df.meta[order(df.meta$node,df.meta$orf),])

cov.node.uniq <- cov.node[!duplicated(cov.meta$node),]
cov.meta.uniq <- cov.meta[!duplicated(cov.meta$node),]
cov.meta.uniq$normAvgDepth <- rowMeans(prop.table(as.matrix(cov.node.uniq),2))

#cov.node.uniq <- unique(cov.node)
#cov.meta.uniq <- cov.meta[rownames(cov.node.uniq),]

rownames(cov.node.uniq) <- cov.meta.uniq$node
rownames(cov.meta.uniq) <- cov.meta.uniq$node

#first coverage matrix with Fast PCA
# cov.pca.list <- FastPCA(scale(log2(cov.node.uniq+1e-6)-log2(1e-6),center=F,scale=T),50)
cov.pca.list <- FastPCA(scale(cov.node.uniq,center=F,scale=T),50)
cov.pca <- cov.pca.list$x

plot(as.data.frame(cov.pca[,1:3]),pch='.')

#then codon usage with Fast PCA
#get codon usage table, takes a minute for 300Mb file
#ffn <- read.fasta("annotation/prokka/prokka_scaffolds.1000bpp/scaffolds.1000bpp.ffn")
#ffn <- saveRDS(ffn,file="ffn.rds")
#ffn <- readRDS(file="ffn.rds")
#ffn.uco.list <- mclapply(ffn,mc.cores=14, FUN=function(x) {uco(x, frame=0, as.data.frame=T, NA.rscu=T)}) #5s for 1k, 3 minutes for 400k on 12 threads
#ffn.uco <- do.call(cbind,ffn.uco.list)
#saveRDS(ffn.uco,file="ffn.uco.rds")
ffn.uco <- readRDS(file="ffn.uco.rds")
ffn.eff <- t(ffn.uco[,grepl(pattern="eff",colnames(ffn.uco))])
ffn.rscu <- t(ffn.uco[,grepl(pattern="RSCU",colnames(ffn.uco))])
ffn.freq <- t(ffn.uco[,grepl(pattern="freq",colnames(ffn.uco))])
rownames(ffn.eff) <- sub(".eff","",rownames(ffn.eff))
rownames(ffn.rscu) <- sub(".RSCU","",rownames(ffn.rscu))
rownames(ffn.freq) <- sub(".freq","",rownames(ffn.freq))

#only keep those rows with matches in df.in
ffn.eff <- ffn.eff[df.meta$orf,]
ffn.rscu <- ffn.rscu[df.meta$orf,]
ffn.freq <- ffn.freq[df.meta$orf,]

#these need to be averaged over scaffolds first
#aggby should be a data.frame as long as x with repeated indices where rows of x should be merged
# aggregate is very slow
# agg.eff <- aggregate(ffn.eff,aggby,"mean")
# agg.rscu <- aggregate(ffn.rscu,aggby,"mean")
# agg.freq <- aggregate(ffn.freq,aggby,"mean")

#use data.table because its so fast
ffn.eff.dt <- as.data.table(ffn.eff)
aggby <- as.data.frame(df.meta$node[match(rownames(ffn.eff), df.meta$orf)])
agg.eff <- as.data.frame(ffn.eff.dt[, lapply(.SD, mean), by = aggby])
rownames(agg.eff) <- agg.eff[,1]
agg.eff <- agg.eff[,-1]

ffn.rscu.dt <- as.data.table(ffn.rscu)
aggby <- as.data.frame(df.meta$node[match(rownames(ffn.rscu), df.meta$orf)])
agg.rscu <- as.data.frame(ffn.rscu.dt[, lapply(.SD, mean), by = aggby])
rownames(agg.rscu) <- agg.rscu[,1]
agg.rscu <- agg.rscu[,-1]

ffn.freq.dt <- as.data.table(ffn.freq)
aggby <- as.data.frame(df.meta$node[match(rownames(ffn.freq), df.meta$orf)])
agg.freq <- as.data.frame(ffn.freq.dt[, lapply(.SD, mean), by = aggby])
rownames(agg.freq) <- agg.freq[,1]
agg.freq <- agg.freq[,-1]

freq.pca.list <- FastPCA(scale(agg.freq,center=F,scale=T),50)
freq.pca <- freq.pca.list$x
plot(as.data.frame(freq.pca[,1:2]),pch='.')

rscu.pca.list <- FastPCA(scale(agg.rscu,center=F,scale=T),50)
rscu.pca <- rscu.pca.list$x
plot(as.data.frame(rscu.pca[,1:2]),pch='.')

eff.pca.list <- FastPCA(scale(agg.eff,center=F,scale=T),50)
eff.pca <- eff.pca.list$x
plot(as.data.frame(eff.pca[,1:2]),pch='.')

codon.pca.list <- FastPCA(scale(cbind(freq.pca,rscu.pca,eff.pca),center=F,scale=T),50)
codon.pca <- codon.pca.list$x
#rownames are nodes already
#rownames(codon.pca) <- rownames(df.meta)[match(rownames(agg.freq),df.meta$node)]

#center each of the columns before going into Rtsne
codon.pca <- scale(codon.pca, center=T, scale=T)
cov.pca <- scale(cov.pca, center=T, scale=T)

# should be fixed now
# cov is contigs
# codon is contigs

goodrows <- intersect(rownames(cov.pca),rownames(codon.pca))
cov.pca.good <- cov.pca[goodrows,]
codon.pca.good <- codon.pca[goodrows,]
cc.pca.good <- cbind(cov.pca.good,codon.pca.good)

#how many samples is the contig found in?
# pdf(file="rs.pdf",width=40,height=40)
# rs <- rowSums(cov.node.uniq>0)
# hist(aa,breaks=96)
# bin.list <- cov.meta.uniq[,"metabat.bin"]
# bin.colors <- rainbow(length(unique(cov.meta.uniq$metabat.bin)))[as.numeric(as.factor(bin.list))]
# plot(rs, log10(cov.meta.uniq$contigLen), pch=21,bg=alpha("grey",0.4),col=NULL,cex=2)
# points(rs, log10(cov.meta.uniq$contigLen), pch=21,bg=alpha(bin.colors,0.4),col=NULL,cex=1)
# dev.off()

### 
# depth-first traverse
# or sample-freq first?
# cluster highest-coverage sequences first, 
# remove them from consideration once binned in bins >200kbp
# or replace them with the median? of the bin
# repeat
# then re-assemble

# include a fork where bins are only built within (heatmap) clades of related samples 
# balance between sequencing depth and diversity of samples
# related bins can be merged after
# a heatmap of core genes in each genome with phylogeny on the row dendrograms

# with additional samples, use good bins as magnets to attract new scaffolds for faster binning


h <- hist(cov.meta$totalAvgDepth,breaks=200)
plot(log10(h$counts) ~ h$mids)


## use ribosomal protein sequence distribution to determine sample clusters
## then bin by sample cluster

#use as a prior to determine how many bins to look for?

#FIRST: cluster ribosomal proteins by tsne3d
ribo.orf <- unique(cov.meta[grep("ribosomal",cov.meta$prokka.ann,ignore.case=T),"orf"])
ribo.meta <- df.meta[rownames(df.num[ribo.orf,]),]
ribo.mat <- df.num[rownames(ribo.meta),]
ribo.mat <- ribo.mat[!duplicated(ribo.meta$node),]
ribo.meta <- ribo.meta[!duplicated(ribo.meta$node),]
rownames(ribo.mat) <- ribo.meta$node
ribo.pca.list <- FastPCA(ribo.mat,50)
ribo.pca <- ribo.pca.list$x
ribo.pca <- cbind(ribo.pca,ribo.pca[,50]+(1:nrow(ribo.pca))*1e-6)
ribo.tsne.list <- Rtsne(ribo.pca,dims=2, max_iter=1000, verbose=TRUE, theta=0.5, pca=FALSE)
ribo.tsne <- ribo.tsne.list$Y

rs <- rowSums(ribo.mat>0)[match(rownames(ribo.pca),ribo.meta$node)]
# tax.list <- tax.out[ribo.meta$orf[match(rownames(ribo.pca),ribo.meta$node)],"family"] #why does this take forever?
# bin.list <- ribo.meta[match(rownames(ribo.pca),ribo.meta$node),"metabat.bin"]
tax.list <- tax.out[ribo.meta$orf[match(rownames(ribo.pca),ribo.meta$node)],"family"] #why does this take forever?
bin.list <- ribo.meta[match(rownames(ribo.pca),ribo.meta$node),"metabat.bin"]
size.list <- ribo.meta$contigLen[match(rownames(ribo.pca),ribo.meta$node)]
plot(rs,log2(size.list),col=tax.colors,pch='.',cex=2)

tax.colors <- rainbow(max(as.numeric(as.factor(tax.list)),na.rm=T))[as.numeric(as.factor(tax.list))]
bin.colors <- rainbow(length(unique(cov.meta.uniq$metabat.bin)))[as.numeric(as.factor(bin.list))]

plot(as.data.frame(ribo.tsne),pch=21,col=NULL,bg=tax.colors,main="taxa")
plot(as.data.frame(ribo.tsne),pch=21,col=NULL,bg=bin.colors,main="bins")

# pdf(file="hm.pdf",width=24,height=24)
# hm <- heatmap.2(t(sqrt(sqrt(as.matrix(ribo.mat)))),trace="none")
hm <- heatmap.2(t((as.matrix(ribo.mat>0)*1)),trace="none",hclustfun = function(x) hclust(x,method = 'ward.D'))
# dev.off()

sample.clusters <- cutree(hm$rowDendrogram,k=5)
heatmap.2(t((as.matrix(ribo.mat>0)*1)),trace="none",hclustfun = function(x) hclust(x,method = 'ward.D'),RowSideColors=rainbow(5)[sample.clusters])

dev.off()


save.image(file="binning.Rdata")

} else {

#load.image("binning.Rdata")

}



#second: kmeans clustering from 50:100
#third: pull in related contigs by coverage/graph/composition
#fourth: assign bins by heirarchical clustering


# will have to figure out what to do when problem contigs get repeatedly pulled to the top because they're not getting incorporated into bins
# counter for how many times it's been seen?
# used at a frequency of 1/counter?

# remove singletons? or try to bin them up first
# codon usage gets more limited at shorter sequences, maybe increase the proportion of other vectors as a function of length

# might want to use longest to shortest rather than highest cov to lowest
# when shorter contigs get included there will be weird outliers with extreme coverage
# but... the bins are related to sequencing depth

# how about... make solid bins with len>5k
# then iterate over each bin, adding contigs in it's sphere
# FastPCA(cov,codon) on bins+contigs, then
# distance matrix from each bin to all others
# replicate distant contigs with coverages at 1/n* and n* to pull in duplications/repeats


#could parallelize this, but want to print. could use ggsave

sample.clusters <- cutree(hm$rowDendrogram,k=1)

output.list <- mclapply(sort(unique(sample.clusters)), mc.cores=12, 
function(clus) {
# for(clus in sample.clusters) {

	node.kmeans <- list()
	tsne.all <- list()
	iter <- 0

	print(clus)

	seenrows <- NULL

	seed <- 1
	
	ndim <- ifelse(sum(sample.clusters==clus)>=50,50,sum(sample.clusters==clus))

	minlen <- 8000
	
# 	pdf(file=paste0("binning",clus,".pdf"),width=20,height=20)

	#first bin those that are ONLY in this cluster? not for now
	#select samples in cluster, select contigs in samples
	cov.sample <- cov.node.uniq[,sample.clusters==clus]

	#only keep very large contigs
	cov.sample <- cov.sample[which(cov.meta.uniq$contigLen > minlen),]

	#removes singletons/should bin them first to get novel stuff?
	cov.sample <- cov.sample[which(rowSums(cov.sample>0)>0),]

	cov.sample.meta <- cov.meta.uniq[rownames(cov.sample),]

	cov.sample.uniq <- cov.sample[!duplicated(cov.sample.meta$node),]
	cov.sample.meta.uniq <- cov.sample.meta[!duplicated(cov.sample.meta$node),]
	rownames(cov.sample.uniq) <- cov.sample.meta.uniq$node
	rownames(cov.sample.meta.uniq) <- cov.sample.meta.uniq$node
	
	#calculate these on reduced dataset
	cov.pca.sample.list <- FastPCA(cov.sample.uniq,ndim)
	cov.pca.sample <- cov.pca.sample.list$x
# 	cov.pca.sample <- cov.pca.list$x[rownames(cov.sample.uniq),]

	codon.sample.uniq <- cbind(freq.pca[rownames(cov.sample.uniq),],rscu.pca[rownames(cov.sample.uniq),],eff.pca[rownames(cov.sample.uniq),])
	codon.pca.sample.list <- FastPCA(codon.sample.uniq,ndim)
	codon.pca.sample <- codon.pca.sample.list$x
#	codon.pca.sample <- codon.pca.list$x[rownames(cov.sample.uniq),]

	saveseed <- 0

	#could parallelize this? no because depends on previous iter
	while(length(seenrows) < nrow(cov.sample.meta.uniq)) {
		iter <- iter + 1

		print(iter)
		
		for(seed in seq(min(setdiff(1:nrow(cov.sample.meta.uniq),saveseed)),nrow(cov.sample.meta.uniq),200)) {

			#find the minimum coverage to consider as a seed
			sample.prop <- rowSums(prop.table(as.matrix(cov.sample.uniq),2))
			mincov <- rev(sort(sample.prop))[seed]
			seedrows <- which(sample.prop >= mincov)
			seedrows <- setdiff(seedrows,seenrows)
			maxhops <- 4

			#pull in N neighbors up to N=maxhops
			newrows <- seedrows
			for(i in maxhops) {
				vrows <- match(cov.sample.meta.uniq$node[newrows],V(gpaths)$name)
				vrows <- vrows[!is.na(vrows)]
				nbors <- names(unlist(ego(gpaths, i, nodes = vrows, mode = "all", mindist = i-1))) 
				newrows <- match(nbors[!is.na(nbors)], cov.sample.meta.uniq$node)
				newrows <- newrows[!is.na(newrows)]
				seedrows <- sort(unique(c(seedrows, newrows)))
				len <- length(seedrows)
				if(len > 100000) break
			}
			if(len > 0) cat(paste("for seed = ",seed, ", total contigs: ", len, "\n"))
			if(len > 100000) break
		}

		if(length(seedrows)<ndim) break

		saveseed <- seed
		seenrows <- unique(c(seenrows,seedrows))

		cov.pca.filt <- cov.pca.sample[seedrows,]
		codon.pca.filt <- codon.pca.sample[seedrows,]

		#make graph pca
		node.new <- setdiff(rownames(cov.pca.filt),V(gpaths)$name)
		cat("contigs with links: ", length(intersect(rownames(cov.pca.filt),V(gpaths)$name)),"\n")
		gpaths.sample <- add_vertices(gpaths, length(node.new), attr=list(name=node.new))
		graph.dist <- as.matrix(distances(gpaths.sample, rownames(cov.pca.filt), rownames(cov.pca.filt), mode="all"))
		graph.dist[is.infinite(graph.dist)] <- NA
		graph.dist[is.finite(graph.dist)] <- graph.dist[is.finite(graph.dist)]/max(graph.dist[is.finite(graph.dist)])^2
		graph.dist[is.na(graph.dist)] <- 1

		graph.pca.list <- FastPCA(scale(graph.dist),ndim)
		graph.pca.filt <- graph.pca.list$x[rownames(cov.pca.filt),]

		#this is where I could change the balance using contigLen, etc.
		cc.pca.list <- FastPCA(cbind(cov.pca.filt, codon.pca.filt, graph.pca.filt[,1:3]),ndim)
		cc.pca.filt <- cc.pca.list$x[rownames(cov.pca.filt),]
		
		#add dummy col to make them unique
		cov.pca.filt	<-   cbind(cov.pca.filt,     cov.pca.filt[,ndim]+(1:nrow(cov.pca.filt))/1e8)
		codon.pca.filt	<- cbind(codon.pca.filt, codon.pca.filt[,ndim]+(1:nrow(codon.pca.filt))/1e8)
		graph.pca.filt	<- cbind(graph.pca.filt, graph.pca.filt[,3]+(1:nrow(graph.pca.filt))/1e8)
		cc.pca.filt		<-    cbind(cc.pca.filt,       cc.pca.filt[,ndim]+(1:nrow(cc.pca.filt))/1e8)

		# then Barnes-Hut t-SNE in 3D on combined PCAs
		# job.out <- mclapply(c("cov.pca","codon.pca","cc.pca"),mc.cores=3, function(x) {Rtsne(unique(eval(parse(text=x))), dims=3, max_iter=1500, verbose=TRUE, theta=0.5, pca=FALSE)})
		job.out <- mclapply(c("cov.pca.filt","codon.pca.filt","cc.pca.filt","graph.pca.filt"),mc.cores=4, function(x) {
			Rtsne(eval(parse(text=x)), dims=2, max_iter=2000, verbose=TRUE, theta=0.3, pca=FALSE)
		})

		tsne.all[[iter]] <- job.out

		cov.tsne3d <- job.out[[1]]$Y
		codon.tsne3d <- job.out[[2]]$Y
		cc.tsne3d <- job.out[[3]]$Y
		graph.tsne3d <- job.out[[4]]$Y
		
		# do kmeans here on high k, then aggregate and add each of those as bins?
		# and remove component contigs from table
		fkm.cov <- fastKmeans(cov.tsne3d, 100, iter.max = 10000, project = TRUE, threads = 0)$class
		fkm.codon <- fastKmeans(codon.tsne3d, 100, iter.max = 10000, project = TRUE, threads = 0)$class
		fkm.cc <- fastKmeans(cc.tsne3d, 100, iter.max = 10000, project = TRUE, threads = 0)$class
		fkm.graph <- fastKmeans(graph.tsne3d, 100, iter.max = 10000, project = TRUE, threads = 0)$class

		fkm.cov.s <- paste("cov",iter,fkm.cov,sep="_")
		fkm.codon.s <- paste("codon",iter,fkm.codon,sep="_")
		fkm.cc.s <- paste("cc",iter,fkm.cc,sep="_")
		fkm.graph.s <- paste("cc",iter,fkm.graph,sep="_")

		node.kmeans[[iter]] <- cbind(fkm.cov.s,fkm.codon.s,fkm.cc.s,fkm.graph.s)
		rownames(node.kmeans[[iter]]) <- rownames(cov.pca.filt)
		
		# saveRDS(cov.tsne3d.list,file="cov.tsne3d.list.rds")
		# saveRDS(cov.tsne3d.list,file="codon.tsne3d.list.rds")
		# saveRDS(cov.tsne3d.list,file="cc.tsne3d.list.rds")

#		tax.list <- tax.out[cov.sample.meta.uniq[seedrows,"orf"],"family"] #why does this take forever?
		tax.list <- tax.out$family[match(cov.sample.meta.uniq[seedrows,"orf"],rownames(tax.out))] #why does this take forever?
		bin.list <- cov.sample.meta.uniq[seedrows,"metabat.bin"]

		#tax.colors <- viridis_pal()(max(as.numeric(as.factor(tax.list)),na.rm=T))[as.numeric(as.factor(tax.list))]
		#bin.colors <- viridis_pal()(length(unique(cov.meta.uniq$metabat.bin)))[as.numeric(as.factor(bin.list))]

		tax.colors <- rainbow(max(as.numeric(as.factor(tax.list)),na.rm=T))[as.numeric(as.factor(tax.list))]
		bin.colors <- rainbow(length(unique(cov.sample.meta.uniq$metabat.bin)))[as.numeric(as.factor(bin.list))]

		size.contigs <- sqrt(sqrt(cov.sample.meta.uniq$contigLen[seedrows]))/10
		size.depth <- sqrt(cov.sample.meta.uniq$totalAvgDepth[seedrows])/5

		pdf(file=paste0("binning",clus,"_",iter,".pdf"),width=20,height=20)
		plot(cov.tsne3d,  pch=21,col=alpha("grey",0.2),bg=tax.colors,cex=size.contigs,main="coverage + tax + contigLen")
		plot(cov.tsne3d,  pch=21,col=alpha("grey",0.2),bg=bin.colors,cex=size.contigs,main="coverage + bin + contigLen")

		plot(codon.tsne3d,pch=21,col=alpha("grey",0.2),bg=tax.colors,cex=size.contigs,main="codon + tax + contigLen")
		plot(codon.tsne3d,pch=21,col=alpha("grey",0.2),bg=bin.colors,cex=size.contigs,main="codon + bin + contigLen")

		plot(graph.tsne3d,   pch=21,col=alpha("grey",0.2),bg=tax.colors,cex=size.contigs,main="graph + tax + contigLen")
		plot(graph.tsne3d,   pch=21,col=alpha("grey",0.2),bg=bin.colors,cex=size.contigs,main="graph + bin + contigLen")

		plot(cc.tsne3d,   pch=21,col=alpha("grey",0.2),bg=tax.colors,cex=size.contigs,main="cc + tax + contigLen")
		plot(cc.tsne3d,   pch=21,col=alpha("grey",0.2),bg=bin.colors,cex=size.contigs,main="cc + bin + contigLen")

		plot(cc.tsne3d,   pch=21,col=alpha("grey",0.2),bg=tax.colors,cex=size.depth,main="cc + tax + depth")
		plot(cc.tsne3d,   pch=21,col=alpha("grey",0.2),bg=bin.colors,cex=size.depth,main="cc + bin + depth")
		dev.off()

	}
	
list("kmeans"=node.kmeans,"tsne"=tsne.all)

}
)

# node.all <- do.call(rbind(node.kmeans))

#saveRDS(cov.tsne3d.list,file="binning/tsne-embed.rds")
#cov.tsne3d.list <- readRDS(file="binning/tsne-embed.rds")
#cov.tsne3d <- cov.tsne3d.list$Y
#rownames(cov.tsne3d) <- rownames(covtab)

stop()




	

# going back to using very large contigs as seed bins
# now take consensus of kmeans as evidence for a bin
kmeans.out <- output.list[[1]]$kmeans[[1]]

#exclude graph-only
kmeans.mat <- as.data.frame(kmeans.out[,c("fkm.cov.s","fkm.codon.s","fkm.cc.s")],stringsAsFactors=F)
kmeans.bins <- unique(kmeans.mat)

# add bins to database
agg.node <- cov.node.uniq
codon.pca.agg <- codon.pca.list$x
agg.meta <- cov.meta.uniq

add.meta <- NULL
add.node <- NULL

for(bin in 1:nrow(kmeans.bins)) {
	this.bin <- match_df(kmeans.mat,kmeans.bins[bin,])

	if(nrow(this.bin)<2) next
	new.node <- t(as.data.frame(colMeans(agg.node[rownames(this.bin),])))
	rownames(new.node) <- paste0("bin",bin)

#	don't remove component nodes? because some bins might reasonably use same contig
#	agg.node <- agg.node[!(row.names(agg.node) %in% rownames(this.bin)),]
#	agg.node <- rbind(agg.node,new.node)
	add.node <- rbind(add.node,new.node)
	
	new.codon <- t(as.data.frame(colMeans(codon.pca.agg[rownames(this.bin),])))
	codon.pca.agg <- rbind(codon.pca.agg,new.codon)
	
	new.meta <- as.data.frame(agg.meta[rownames(this.bin),],stringsAsFactors=T)
	new.meta <- cbind(NA,NA,NA,NA,sum(new.meta$contigLen),mean(new.meta$totalAvgDepth), mean(new.meta$normAvgDepth))
	rownames(new.meta) <- paste0("bin",bin)
	colnames(new.meta) <- colnames(agg.meta)
# 	agg.meta <- rbind(agg.meta,new.meta)
	add.meta <- rbind(add.meta,new.meta)	
}

agg.node <- rbind(agg.node,add.node)
agg.meta <- rbind(agg.meta,add.meta)

heatmap.2(as.matrix(sqrt(sqrt(add.node))),trace="none",hclustfun = function(x) hclust(x,method = 'ward.D'))
heatmap.2(cor(t(as.matrix((sqrt(agg.node[grep("bin",rownames(agg.node)),]))))),trace="none",hclustfun = function(x) hclust(x,method = 'ward.D'))

#heatmap of ribosomal proteins per bin


# how many genomes can we expect?
aa <- df.meta[rownames(df.num[ribo.orf,]),]
bb <- rle(sort(aa$prokka.ann))
cc <- as.data.frame(cbind(bb$values[rev(order(bb$lengths))],rev(sort(bb$lengths))),stringsAsFactors=F)
cc[,2] <- as.numeric(cc[,2])
num.ribo <- cbind(sum(cc[grep("5S",cc[,1]),2]), sum(cc[grep("16S",cc[,1]),2]), sum(cc[grep("23S",cc[,1]),2]))
colnames(num.ribo) <- c("5S","16S","23S")

rrna.meta <- df.meta[rownames(df.num[ribo.orf,]),]

rrna.meta <- rrna.meta[grep("5S|16S|23S",rrna.meta$prokka.ann),]


cov.bins.pca <- FastPCA(agg.node,50)$x
fkm.bins <- fastKmeans(cov.bins.pca[,1:3], 50, iter.max = 10000, project = TRUE, threads = 0)$class
rownames(fkm.bins) <- rownames(agg.node)

# go fishing for additional contigs
for(bin in 1:nrow(kmeans.bin)) {

	#do rough kmeans on PCA, find which cluster bin is in, do hc on that cluster
	bin.clus <- fkm.bins[paste0("bin",bin),]
	bin.pick <- which(fkm.bins==bin.clus)
	bin.cov <- bins.pca[bin.pick,]
	bin.codon <- codon.pca.agg[bin.pick,]

	bin.tsne <- Rtsne(cbind(bin.cov,bin.codon), max_iter=1000, verbose=TRUE, theta=0.5, pca=FALSE)$Y

	tax.list <- tax.out$family[match(agg.meta[bin.pick,"orf"],rownames(tax.out))] #why does this take forever?
	bin.list <- agg.meta[bin.pick,"metabat.bin"]

	tax.colors <- rainbow(max(as.numeric(as.factor(tax.list)),na.rm=T))[as.numeric(as.factor(tax.list))]
	bin.colors <- rainbow(length(unique(agg.meta$metabat.bin)))[as.numeric(as.factor(bin.list))]

	size.contigs <- sqrt(sqrt(agg.meta$contigLen[seedrows]))/10
	size.depth <- sqrt(agg.meta$normAvgDepth[seedrows])*100

	plot(bin.tsne,  pch=21,col=alpha("grey",0.2),bg=tax.colors,cex=size.contigs,main="coverage + tax + contigLen")
	plot(bin.tsne,  pch=21,col=alpha("grey",0.2),bg=bin.colors,cex=size.contigs,main="coverage + bin + contigLen")

}


#	don't really think this is working
# 	
# 	# correlations instead of distances to find dups?
# 	# 
# 	agg.mat <- as.matrix(agg.node)
# 	bin.row <- which(rownames(agg.mat)==paste0("bin",bin))
# 	cor.agg.bin <- apply(agg.mat,1,function(x) cor(agg.mat[bin.row,],x))
# 
# 	# how far in PC/cor space
# 	cov.pca.agg <- FastPCA(agg.mat,50)$x
# 	dist.agg.bin <- sweep(cov.pca.agg, 2, cov.pca.agg[bin.row,])
# 	plot(cor.agg.bin,log(rowSums(dist.agg.bin[,c(1:3)]^2)),pch='.')
# 
# 	cor.agg.new <- agg.meta[names(which(cor.agg.bin>=0.9)),]
# 





#quickie method/3dkmeans on large number of tests
goodrows <- intersect(rownames(cov.node.uniq),rownames(codon.pca))	

tax.list <- tax.out$family[match(cov.meta.uniq[goodrows,"orf"],rownames(tax.out))] #why does this take forever?
bin.list <- cov.meta.uniq[goodrows,"metabat.bin"]

tax.colors <- rainbow(max(as.numeric(as.factor(tax.list)),na.rm=T))[as.numeric(as.factor(tax.list))]
bin.colors <- rainbow(length(unique(cov.meta.uniq[goodrows,"metabat.bin"])))[as.numeric(as.factor(bin.list))]

size.contigs <- sqrt(sqrt(cov.meta.uniq[goodrows,"contigLen"]))/10
size.depth <- sqrt(rowSums(prop.table(cov.node.uniq)))[goodrows]/5

#plot pairwise log(abundance) + codon.pca1, 3d kmeans, consensus?
pdf(file="allvsall.pdf")
# for(i in 1:ncol(cov.node.uniq)) {
for(i in match(names(rev(sort((colSums(cov.node.uniq))))[1:10]),colnames(cov.node.uniq))) {
	print(i)
for(j in match(names(rev(sort((colSums(cov.node.uniq))))[1:10]),colnames(cov.node.uniq))) {
# 	for(j in seq(i,ncol(cov.node.uniq),8)) {
#	plotmat <- cbind(log10(cov.node.uniq[goodrows,i]+1e-6)-log10(1e-6),log10(cov.node.uniq[goodrows,j]+1e-6)-log10(1e-6),codon.pca[goodrows,2])
for(k in match(names(rev(sort((colSums(cov.node.uniq))))[1:10]),colnames(cov.node.uniq))) {
	if(i==j | i==k | j==k) next
	plotmat <- cbind(sqrt(sqrt(cov.node.uniq[goodrows,i])),sqrt(sqrt(cov.node.uniq[goodrows,j])),sqrt(sqrt(cov.node.uniq[goodrows,k])))

# 	plotmat <- cbind(sqrt(sqrt(cov.node.uniq[goodrows,i])),sqrt(sqrt(cov.node.uniq[goodrows,j])),codon.pca[goodrows,1])
# 	plotmat <- cbind(sqrt(sqrt(cov.node.uniq[goodrows,i])),sqrt(sqrt(cov.node.uniq[goodrows,j])),codon.pca[goodrows,2])
	plotmat[plotmat==0] <- NA
# 	plot3d(plotmat,col=alpha("grey",0.4),cex=size.contigs*2)
	plot3d(plotmat,col=alpha(bin.colors,0.4),cex=size.contigs*5)
	}
}
}
dev.off()



#get lovejoy sequences
