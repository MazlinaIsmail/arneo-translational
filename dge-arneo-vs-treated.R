# based on TR229

library(edgeR)
library(statmod)
edgeRUsersGuide()
library(tximport)
library(biomaRt)
ensembl <- useEnsembl(biomart='genes', dataset='hsapiens_gene_ensembl')
ensembl
library(reshape2)
library(ggridges)
library(ggplot2)

# select samples to exclude known outliers
sampledets <- read.table('/path/to/mastersamplenamekey.txt', header=T, sep='\t')
sampledets_keep <- subset(sampledets, SAMPLE_INCLUSION == 'Y' & (PROJECT == 'ARNEO' | PROJECT == 'PWB'))
nrow(sampledets_keep)
table(sampledets_keep$PROJECT)

dirs_of_interest <- c('/path/to/gene-quant/MI000059_8')
input_files <- list.files(dirs_of_interest, full.names=T)
length(input_files)
keep_input_files <- unlist(lapply(1:nrow(sampledets_keep), function(x){grep(sampledets_keep[x,14], input_files, value=T)}))
length(keep_input_files)
sample_id <- gsub('_1P$', '', sapply(1:length(keep_input_files), function(x) tail(strsplit(strsplit(keep_input_files[x], '/')[[1]], '\\.'), n=1)[[1]][1]))
new_name <- unlist(lapply(1:length(sample_id), function(x){
	sampledets_keep[grep(sample_id[x], sampledets_keep$FILENAME_ID_V2),6]
}))
# check match
data.frame(sample_id, new_name)
names(keep_input_files) <- new_name

# edgeR
txi.rsem <- tximport(keep_input_files, type='rsem', txIn=F, txOut=F, importer=read.delim)
head(txi.rsem$counts)

cts <- txi.rsem$counts
normMat <- txi.rsem$length

table(is.na(cts))
table(is.na(normMat))
dim(cts)
dim(normMat)

txi.rsem$length[txi.rsem$length == 0] <- 1
normMat <- txi.rsem$length
table(is.na(normMat))

# obtaining per-observation scaling factors for length, adjusted to avoid changing the magnitude of the counts
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat
# computing effective library sizes from scaled counts, to account for composition biases between samples
eff.lib <- calcNormFactors(normCts) * colSums(normCts)
# combining effective library sizes with the length factors, and calculating offsets for a log-link GLM
normMat <- sweep(normMat, 2, eff.lib, '*')
normMat <- log(normMat)
# creating a DGEList object for use in edgeR
# get treatment group
arneo_unblind <- read.table('/path/to/arneo-unblinding-key.txt', header=T, sep='\t')
head(arneo_unblind)
treatment <- sapply(1:ncol(cts), function(x) arneo_unblind[grep(colnames(cts)[x], arneo_unblind$PATIENT_ID),4]) 
treatment <- factor(treatment, levels=c('untreated', 'adt', 'adt_apa'))
treatment[is.na(treatment)] <- 'untreated'
treatment
table(treatment)
# design matrix
cts_df <- data.frame(cts)
dim(cts_df)
design <- model.matrix(~0+treatment, data=cts_df)
colnames(design) <- levels(treatment)
y <- DGEList(cts, group=treatment)
y <- scaleOffset(y, normMat)
# filtering
keep <- filterByExpr(y, group=treatment, min.prop=0.9, min.count=50)
table(keep)
y <- y[keep,]
# includes annotation for entrez ID
annot <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id'), filters='ensembl_gene_id', values=rownames(y$counts), mart=ensembl)
table(duplicated(annot$ensembl_gene_id))
# remove duplicated ensembl gene id
d <- duplicated(annot$ensembl_gene_id)
annot <- annot[!d,]
nrow(annot)
nrow(y$counts)
# remove gene IDs not in annot table
keep_id <- rownames(y$counts) %in% annot$ensembl_gene_id
table(keep_id)
y <- y[keep_id,]
nrow(y$counts)
y$genes <- annot

# plot MDS and PCA
col <- ifelse(grepl('untreated', treatment) == T, '#999999', ifelse(grepl('adt_apa', treatment) == T, '#E69F00', '#56B4E9')) # grey, orange, blue
plotMDS(y, col=col, pch=19, main='treatment: grey = naive, orange = adt_apa, blue = adt')
plotMDS(y, col=col, pch=19, gene.selection='common', main='treatment: grey = naive, orange = adt_apa, blue = adt')
plotMDS(y, col=col, gene.selection='common', main='treatment: grey = naive, orange = adt_apa, blue = adt')

# logcpm data
logcpm <- cpm(y, log=T)
# count data
count_df <- y$count

y <- estimateDisp(y, design, robust=T)
y$common.dispersion
plotBCV(y)
fit <- glmQLFit(y, design, robust=T)
plotQLDisp(fit)

# testing for DE genes 

# (1) adt_apa vs untreated
Adt_apaVsUnt <- makeContrasts(adt_apa - untreated, levels=design)
qlf <- glmQLFTest(fit, contrast=Adt_apaVsUnt)
topTags(qlf)
summary(decideTests(qlf))
plotMD(qlf)
apa_res <- topTags(qlf, Inf, adjust.method='BH')$table
tr <- glmTreat(fit, contrast=Adt_apaVsUnt, lfc=log2(1.5))
summary(decideTests(tr))
apa_tr <- topTags(tr, Inf)$table
keg <- kegga(qlf, species='Hs', geneid=qlf$genes$entrezgene_id)
keg_res_down <- topKEGG(keg, sort='Down', number=Inf)
keg_res_down$P.Down_format <- format(keg_res_down$P.Down, scientific=F)
keg_res_down$P.Up_format <- format(keg_res_down$P.Up, scientific=F)
keg_res_up <- topKEGG(keg, sort='Up', number=Inf)
keg_res_up$P.Down_format <- format(keg_res_up$P.Down, scientific=F)
keg_res_up$P.Up_format <- format(keg_res_up$P.Up, scientific=F)

outname <- '/path/to/tables/EDGER/2022-09-19_topKEGG-down-apa_vs_naive.txt'
write.table(keg_res_down, outname, quote=F, sep='\t', col.names=T, row.names=F)
outname <- '/path/to/tables/EDGER/2022-09-19_topKEGG-up-apa_vs_naive.txt'
write.table(keg_res_up, outname, quote=F, sep='\t', col.names=T, row.names=F)

# (2) adt vs untreated
AdtVsUnt <- makeContrasts(adt - untreated, levels=design)
qlf <- glmQLFTest(fit, contrast=AdtVsUnt)
topTags(qlf)
summary(decideTests(qlf))
plotMD(qlf)
adt_res <- topTags(qlf, Inf, adjust.method='BH')$table
tr <- glmTreat(fit, contrast=AdtVsUnt, lfc=log2(1.5))
summary(decideTests(tr))
adt_tr <- topTags(tr, Inf)$table
keg <- kegga(qlf, species='Hs', geneid=qlf$genes$entrezgene_id)
keg_res_down <- topKEGG(keg, sort='Down', number=Inf)
keg_res_down$P.Down_format <- format(keg_res_down$P.Down, scientific=F)
keg_res_down$P.Up_format <- format(keg_res_down$P.Up, scientific=F)
keg_res_up <- topKEGG(keg, sort='Up', number=Inf)
keg_res_up$P.Down_format <- format(keg_res_up$P.Down, scientific=F)
keg_res_up$P.Up_format <- format(keg_res_up$P.Up, scientific=F)

outname <- '/path/to/tables/EDGER/2022-09-19_topKEGG-down-adt_vs_naive.txt'
write.table(keg_res_down, outname, quote=F, sep='\t', col.names=T, row.names=F)
outname <- '/path/to/tables/EDGER/2022-09-19_topKEGG-up-adt_vs_naive.txt'
write.table(keg_res_up, outname, quote=F, sep='\t', col.names=T, row.names=F)


# save objects in one file
# objects: count, logcpm, design, adt_res, apa_res
outname <- '/path/to/stable-data-table/arneo/TR229-arneo-edger-output_v2.RData'
save(count_df, logcpm, design, adt_res, apa_res, y, sampledets_keep, file=outname)

# normalize and voom transformation
y <- calcNormFactors(y)
v <- voom(y, design)

# plot dist
y_logcpm <- reshape2::melt(logcpm)
ggplot(y_logcpm, aes(y=Var2)) + geom_density_ridges(aes(x=value)) + theme_ridges()

voom_dat <- reshape2::melt(v$E)
ggplot(voom_dat, aes(y=Var2)) + geom_density_ridges(aes(x=value)) + theme_ridges()

outname <- '/path/to/stable-data-table/arneo/TR229-arneo-edger-output_v3.RData'
save(count_df, logcpm, design, adt_res, apa_res, y, sampledets_keep, v, voom_dat, file=outname)
