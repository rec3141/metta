# read recast matrices back in
library(gplots)
#devtools::install_github('talgalili/dendextend')
library(dendextend)
library(data.table)

find_trees <- function (dend, selected_labels) {
	if (all(get_leaves_attr(dend,"label") %in% selected_labels)) { return(dend)}
	if (any(get_leaves_attr(dend[[1]],"label") %in% selected_labels)) { return(find_trees(dend[[1]], selected_labels))}
	else {return(find_trees(dend[[2]], selected_labels))}
}

get_subtrees <- function (dend, k) {
	clusters <- cutree(dend, k)
	dend_list <- lapply(unique(clusters), function(cluster.id) { find_trees(dend, names(which(clusters == cluster.id))) })
	class(dend_list) <- "dendlist"
	dend_list
}

splitree <- function(dend,max.size=96) {
	saved <- list()
	ct <- get_subtrees(dend,2)
	for(subtree in ct) {
		lab <- get_leaves_attr(subtree,"label")
		if(length(lab) > max.size) {
			saved[[length(saved)+1]] <- splitree(subtree,max.size)
		} else {
# 			print(length(lab))
			saved[[length(saved)+1]] <- lab
		}
	}
	flatten <- function(x) {
	  if (!inherits(x, "list")) return(list(x))
	  else return(unlist(c(lapply(x, flatten)), recursive = FALSE))
	}
	flatten(saved)
}


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


gcut <- function(data,split,range=y) {do.call(paste, data.frame(do.call(rbind, lapply(strsplit(data,split,fixed=T), "[", range))))}

sample.names <- read.table("sample-names.txt",sep="\t",row.names=1,stringsAsFactors=F)
rownames(sample.names) <- gsub("_",".",rownames(sample.names))

#plot all FOAM

pdf(file=paste0(file.path(kaijudir,"foam.all.pdf")),width=20,height=20)

for(filename in list.files(path=kaijudir,pattern="foam.*_.*.tsv")) {

	print(filename)
	patt <- gcut(data=filename,split=".",range=2)
	foam.in <- read.table(file.path(kaijudir,filename),sep="\t",row.names=1)
	rownames(foam.in) <- gcut(data=rownames(foam.in),split="_",range=1)
	rownames(foam.in) <- make.unique(sample.names$V2[match(rownames(foam.in),rownames(sample.names))])
#	cbind(rownames(sample.names)[match(rownames(foam.in),rownames(sample.names))],rownames(foam.in))

	hm <- heatmap.2(sqrt(sqrt(as.matrix(foam.in))),trace="none",margins=c(12,12),main=patt)
	splits <- splitree(hm$colDendrogram,192)
	for(j in 1:length(splits)) {
		print(j)
		foam.min <- as.matrix(foam.in[,splits[[j]]])
		if(ncol(foam.min)==1) foam.min <- cbind(foam.min,foam.min)
		colnames(foam.min) <- do.call(paste,cbind(colnames(foam.min),foam.meta[colnames(foam.min),"ko",drop=F]))
		hm <- heatmap.2(sqrt(sqrt(foam.min)),trace="none",margins=c(12,12),main=paste0(patt,".",j))	
	}
}

dev.off()



#plot all KEGG
pdf(file=paste0("kegg.all.pdf"),width=20,height=20)

for(filename in list.files(pattern="^kegg.*.tsv")) {

	if(filename=="kegg.mat.tsv") next
	
	print(filename)
	patt <- gcut(data=filename,split=".",range=2)
	foam.in <- read.table(filename,sep="\t",row.names=1)
	rownames(foam.in) <- gcut(data=rownames(foam.in),split="_",range=1)
	rownames(foam.in) <- make.unique(sample.names$V2[match(rownames(foam.in),rownames(sample.names))])
#	cbind(rownames(sample.names)[match(rownames(foam.in),rownames(sample.names))],rownames(foam.in))

	hm <- heatmap.2(sqrt(sqrt(as.matrix(foam.in))),trace="none",margins=c(12,12),main=patt)
	max.size=192
	dend <- hm$colDendrogram
	if(nleaves(dend)<=max.size) {
		splits <- list(get_leaves_attr(dend,"label"))
	} else {
		splits <- splitree(dend,max.size)
	}
	for(j in 1:length(splits)) {
		print(j)
		foam.min <- as.matrix(foam.in[,splits[[j]]])
		if(ncol(foam.min)==1) foam.min <- cbind(foam.min,foam.min)
		colnames(foam.min) <- do.call(paste,cbind(colnames(foam.min),foam.meta[colnames(foam.min),"ko",drop=F]))
		hm <- heatmap.2(sqrt(sqrt(foam.min)),trace="none",margins=c(12,12),main=paste0(patt,".",j))	
	}
}

dev.off()
