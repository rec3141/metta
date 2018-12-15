# this program merges FOAM and KEGG output files with and recasts them as matrices with samples in the column and hmms in the rows

library(data.table)
library(gplots)

prefix="scaffolds_gt2000"

foamdir="annotation/foam"
kaijudir="annotation/kaiju"

kh <- read.table(file.path(foamdir,"kegg.headers.txt"),head=T)
fh <- read.table(file.path(foamdir,"foam.headers.txt"),head=T)

foamfile <- paste0(prefix,".foam.depths.tsv")
keggfile <- paste0(prefix,".kegg.depths.tsv")

foamDT <- fread(file.path(kaijudir,foamfile),sep="\t",head=F,stringsAsFactors=F)
keggDT <- fread(file.path(kaijudir,keggfile),sep="\t",head=F,stringsAsFactors=F)

colnames(foamDT) <- colnames(fh)
colnames(keggDT) <- colnames(kh)

rownames(foamDT) <- make.unique(foamDT$foam_hmm)
rownames(keggDT) <- make.unique(keggDT$foam_hmm)

foam.meta <- as.data.frame(foamDT[,1:11,with=F])
kegg.meta <- as.data.frame(keggDT[,1:13,with=F])

rownames(foam.meta) <- rownames(foamDT)
rownames(kegg.meta) <- rownames(keggDT)

foam.num <- as.data.frame(foamDT[,12:ncol(foamDT),with=F])
foam.num <- foam.num[,grep(".var",colnames(foam.num),invert=T)]

kegg.num <- as.data.frame(keggDT[,14:ncol(keggDT),with=F])
kegg.num <- kegg.num[,grep(".var",colnames(kegg.num),invert=T)]

#works but is not faster
#	agg_orf <- setDT(kegg.in)[, lapply(.SD, mean, na.rm=TRUE), by=list(contig,orf,ec) ]
#	agg_orf <- setDT(kegg.in)[, lapply(colnames(kegg.in)[14:ncol(kegg.in)], mean, na.rm=TRUE), by=list(contig,orf,ec) ]

gcut <- function(data,split,range=y) {do.call(paste, data.frame(do.call(rbind, lapply(strsplit(data,split,fixed=T), "[", range))))}


#recast foam matrix

#subset the dataset
# for(patt in sort(gcut(data=unique(foam.meta$foam_1),split=" ",range=1))) {
pattlist <- sort(gcut(data=unique(foam.meta$foam_1),split=" ",range=1))

foam.list <- mclapply(pattlist, function(patt) {
	print(patt)
	picked <- grepl(pattern=patt,foam.meta$foam_1)

	foam.out <- NULL

	for(col in c("foam_1","foam_2","foam_3","foam_4","ko","foam_hmm")) {
		print(col)
		if(all(is.na(foam.meta[picked,col]))) next
		agg_orf <- aggregate(foam.num[picked,], by=list(foam.meta$contig[picked],foam.meta$orf[picked],as.matrix(foam.meta[picked,col])), mean, na.rm=T, simplify = TRUE, drop = TRUE)
		colnames(agg_orf) <- c("contig","orf",col,colnames(agg_orf)[4:ncol(agg_orf)])
		agg_contig <- aggregate(agg_orf, by=list(agg_orf$contig,agg_orf[,col]), mean, na.rm=T, simplify=T,drop=T)
		agg_col <- aggregate(agg_contig, by=list(agg_contig$Group.2), mean, na.rm=T, simplify=T, drop=T)
		rownames(agg_col) <- agg_col$Group.1
		agg_out <- t(agg_col[,-c(1:6)])
		foam.out <- cbind(foam.out,agg_out)
	}

	foam.out

}, mc.cores=21)

# write out and plot
for(i in 1:length(foam.list)) {
	foam.out <- foam.list[[i]]
	patt <- pattlist[i]
 	write.table(foam.out,file=file.path(kaijudir,paste0("foam.",patt,".tsv")),sep="\t",quote=F)

	picked <- unique(sample(1:ncol(foam.out),10000, replace=T))
	pdf(file=file.path(kaijudir,paste0("foam.",patt,".pdf")),width=20,height=20)
	heatmap.2(sqrt(sqrt(foam.out[,picked])),trace="none",margins=c(12,12))
	dev.off()
}


#recast kegg matrix
# pattlist <- gcut(data=unique(kegg.meta$ko_2),split=" ",range=2:7)
# pattlist <- gsub(pattern=" NA",repl="",x=pattlist)
# pattlist <- sort(pattlist)[35:52]

#pattlist <- make.names(unlist(lapply(unique(kegg.meta$ko_2),substr,7,1000)))
pattlist <- unlist(lapply(unique(kegg.meta$ko_2),substr,7,1000))
kegg.list <- mclapply(pattlist, function(patt) {

#for(patt in sort(patterns)[35:52]) {
	print(patt)
	if(is.na(patt)) { return(NULL); next }
	picked <- grepl(pattern=patt,kegg.meta$ko_2)

	kegg.out <- NULL
	#for(col in c("ko","ko_1","ko_2","ko_3","gene","description","ec")) {
	for(col in c("ko_1","ko_2","ko_3","ko")) {

		print(col)
		if(all(is.na(kegg.meta[picked,col]))) next

		agg_orf <- aggregate(kegg.num[picked,], by=list(kegg.meta$contig[picked],kegg.meta$orf[picked],as.matrix(kegg.meta[picked,col])), mean, na.rm=T, simplify = TRUE, drop = TRUE)
		colnames(agg_orf) <- c("contig","orf",col,colnames(agg_orf)[4:ncol(agg_orf)])
		agg_contig <- aggregate(agg_orf, by=list(agg_orf$contig,agg_orf[,col]), mean, na.rm=T, simplify=T,drop=T)
		agg_col <- aggregate(agg_contig, by=list(agg_contig$Group.2), mean, na.rm=T, simplify=T, drop=T)
		rownames(agg_col) <- agg_col$Group.1
		agg_out <- t(agg_col[,-c(1:6)])
		kegg.out <- cbind(kegg.out,agg_out)
		}
		
		kegg.out

}, mc.cores=53)

#write out and plot
for(i in 1:length(kegg.list)) {
	kegg.out <- kegg.list[[i]]
	if(is.null(kegg.out)) next
	if(class(kegg.out)=="try-error") next
	patt <- pattlist[i]
	print(patt)
 	write.table(kegg.out,file=file.path(kaijudir,paste0("kegg.",make.names(patt),".tsv")),sep="\t",quote=F)

	picked <- unique(sample(1:ncol(kegg.out),10000, replace=T))
	pdf(file=file.path(kaijudir,paste0("kegg.",make.names(patt),".pdf")),width=20,height=20)
	heatmap.2(sqrt(sqrt(kegg.out)),trace="none",margins=c(12,12))
	dev.off()
}






Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#aggregate by bin
dt.df <- as.data.frame(dt.in)
dt.num <- dt.df[,seq(10,759,2)]
dt.meta <- dt.df[,1:9]

agg.num <- dt.num[!duplicated(dt.meta$V3),]
agg.meta <- dt.meta[!duplicated(dt.meta$V3),]

#each orf gets a vote..., could also weight by length...
agg_bin <- aggregate(agg.num, by=list(agg.meta$V4), mean, na.rm=T, simplify = TRUE, drop = TRUE)
agg_size <- aggregate(agg.meta$V6, by=list(agg.meta$V4), sum, na.rm=T, simplify = TRUE, drop = TRUE)
agg_met <- aggregate(agg.meta, by=list(agg.meta$V4), Mode)
agg_met[order(-agg_size$x),]

rownames(agg_bin) <- agg_bin[,1]
agg_bin <- agg_bin[,-1]

agg_big <- agg_bin[agg_size$x >3e3,]

# write out and plot
	pdf(file=file.path("annotation/bins",paste0("agg.bin.pdf")),width=20,height=20)
	heatmap.2(sqrt(sqrt(as.matrix(agg_big))),trace="none",margins=c(12,12))
	dev.off()

agg_out <- cbind(agg_size, agg_met)  #agg.meta$V2[match(agg_size$Group.1,agg.meta$V4)]
# colnames(agg_out) <- c("bin","size","taxonomy")
# head(agg_out[order(-agg_out$size,agg_out$taxonomy),])

#could also aggregate by taxonomy here...

#tsne on aggregated bins

agg.pca.samples <- FastPCA(t(as.matrix(agg_big)),50)
agg.pca.bins <- FastPCA(as.matrix(agg_big),50)

agg.tsne.samples <- Rtsne(agg.pca.samples$x,check_duplicates=F,verbose=T)
agg.tsne.bins <- Rtsne(agg.pca.bins$x,check_duplicates=F,verbose=T)

pca.samples <- FastPCA(t(as.matrix(agg.num)),50)
tsne.samples <- Rtsne(pca.samples$x,check_duplicates=F,verbose=T)

pdf(file=file.path("annotation/bins",paste0("agg.tsne.pdf")),width=20,height=20)
	plot(tsne.samples$Y,pch=21,bg=alpha("blue",0.3),main="samples")
	plot(agg.tsne.samples$Y,pch=21,bg=alpha("blue",0.3),main="samples agg")
	plot(agg.tsne.bins$Y,pch=21,bg=alpha("blue",0.3),main="bins agg")
dev.off()



