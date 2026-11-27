# Computational Biomedicine - Tutorial 4
# Wright-Fisher simulation and GEA quality control
# Run from a directory containing data/ (or the three .txt files).
# IMPORTANT: this script makes QC plots; it does not perform differential expression testing.

options(repos = c(CRAN = 'https://cloud.r-project.org'))
if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager')
if (!requireNamespace('edgeR', quietly=TRUE)) BiocManager::install('edgeR', ask=FALSE, update=FALSE)
if (!requireNamespace('ggplot2', quietly=TRUE)) install.packages('ggplot2')
if (!requireNamespace('limma', quietly=TRUE)) BiocManager::install('limma', ask=FALSE, update=FALSE)
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(ggplot2))

out <- 'Week4_R_outputs'
dir.create(out, showWarnings=FALSE, recursive=TRUE)
find_data <- function(filename) {
  candidates <- c(file.path('data', filename), filename,
                  file.path('Week4_Complete','data',filename))
  found <- candidates[file.exists(candidates)]
  if (!length(found)) stop('Missing ', filename, '. Set working directory to the folder with data/.')
  found[1]
}

# 2. Load and align RNA-seq counts and metadata
seqdata <- read.delim(find_data('GSE60450_Lactation-GenewiseCounts.txt'),
                      check.names=FALSE)
sampleinfo <- read.delim(find_data('SampleInfo.txt'), stringsAsFactors=FALSE)
sampleinfo_corrected <- read.delim(find_data('SampleInfo_Corrected.txt'),stringsAsFactors=FALSE)
stopifnot(all(c('EntrezGeneID','Length') %in% names(seqdata)))
cat('Raw table dimensions:', dim(seqdata),'\n')
print(head(seqdata[,seq_len(min(5,ncol(seqdata)))]))
print(head(sampleinfo))
# First two columns contain gene annotation, not sequencing counts.
rownames(seqdata) <- as.character(seqdata$EntrezGeneID)
counts <- as.matrix(seqdata[, -(1:2),drop=FALSE])
storage.mode(counts) <- 'integer'
# Counts use MCL1-DG whereas metadata use MCL1.DG
colnames(counts) <- sub('^MCL1-([A-Z]{2}).*$', 'MCL1.\\1', colnames(counts))
stopifnot(!anyDuplicated(colnames(counts)), !anyDuplicated(rownames(counts)))
match_metadata <- function(info, wanted) {
  idx <- match(wanted,info$SampleName)
  if (anyNA(idx)) stop('Sample IDs missing in metadata.')
  info <- info[idx,,drop=FALSE]
  rownames(info) <- info$SampleName
  stopifnot(identical(info$SampleName,wanted))
  info
}
sampleinfo <- match_metadata(sampleinfo,colnames(counts))
sampleinfo_corrected <- match_metadata(sampleinfo_corrected,colnames(counts))
write.csv(sampleinfo,file.path(out,'04_original_metadata.csv'),row.names=FALSE)
write.csv(sampleinfo_corrected,file.path(out,'05_corrected_metadata.csv'),row.names=FALSE)
changes <- data.frame(SampleName=sampleinfo$SampleName,
  old_CellType=sampleinfo$CellType,new_CellType=sampleinfo_corrected$CellType)
changes <- changes[changes$old_CellType != changes$new_CellType,,drop=FALSE]
write.csv(changes,file.path(out,'06_annotation_changes.csv'),row.names=FALSE)

# 3. EdgeR raw-count data container and library QC
y <- DGEList(counts=counts)
cat('DGEList structure\n'); str(y,max.level=2)
cat('DGEList samples\n'); print(y$samples)
lib_sizes <- data.frame(sample=colnames(y),library_size=y$samples$lib.size)
write.csv(lib_sizes,file.path(out,'07_library_sizes.csv'),row.names=FALSE)
pdf(file.path(out,'08_library_sizes.pdf'),width=11,height=6)
par(mar=c(10,4,3,1)); barplot(y$samples$lib.size/1e6,
 names.arg=colnames(y),las=2,cex.names=.75,
 ylab='Library size (millions of counts)',col='#4C95A8',
 main='Original RNA-seq library sizes')
dev.off()

# Filter genes: CPM > 0.5 in at least three of the 12 samples
myCPM <- cpm(y)
keep <- rowSums(myCPM>0.5)>2
cat('Raw genes:',nrow(y),' Retained genes:',sum(keep),'\n')
y <- y[keep,,keep.lib.sizes=TRUE]
cat('Filtered count matrix:',paste(dim(y$counts),collapse=' x '),'\n')
write.csv(data.frame(raw_genes=length(keep),filtered_genes=sum(keep),
                     samples=ncol(y)), file.path(out,'09_filter_summary.csv'),row.names=FALSE)

# Log CPM QC boxplot
logcpm <- cpm(y, log=TRUE)
overall_median <- median(as.vector(logcpm))
pdf(file.path(out,'10_logCPM_boxplot.pdf'),width=12,height=6)
par(mar=c(10,4,3,1)); boxplot(as.data.frame(logcpm),las=2,cex.axis=.8,
 col='#96C7CD',ylab='Log2 counts per million',main='Log-CPM after filtering')
abline(h=overall_median,lty=2,col='firebrick',lwd=2)
dev.off()

# MDS uses the same count data in both cases; only label colors change.
plot_mds_qc <- function(info, variable, filename, title) {
  f <- factor(info[[variable]])
  colors <- setNames(grDevices::hcl.colors(nlevels(f),'Dark 3'),levels(f))
  pdf(file.path(out,filename),width=10,height=7)
  limma::plotMDS(y, labels=info$SampleName,col=colors[as.character(f)],
                 main=title,cex=.85)
  legend('topright',legend=levels(f),col=colors,pch=16,cex=.9,bty='n')
  dev.off()
}
plot_mds_qc(sampleinfo,'CellType','11_MDS_original_CellType.pdf','MDS: original CellType')
plot_mds_qc(sampleinfo,'Status','12_MDS_original_Status.pdf','MDS: original Status')
plot_mds_qc(sampleinfo_corrected,'CellType','13_MDS_corrected_CellType.pdf','MDS: corrected CellType')
plot_mds_qc(sampleinfo_corrected,'Status','14_MDS_corrected_Status.pdf','MDS: corrected Status')

saveRDS(y, file.path(out,'y_filtered.rds'))
saveRDS(sampleinfo_corrected,file.path(out,'sampleinfo_corrected.rds'))
cat('\nDONE: Outputs saved in ',normalizePath(out),'\n',sep='')
cat('Important: MDS clustering suggests anomalies but does not independently prove label errors.\n')
