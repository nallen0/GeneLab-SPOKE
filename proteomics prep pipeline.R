{\rtf1\ansi\ansicpg1252\cocoartf2709
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 ## Load library\
\
library(DEP)\
library(dplyr)\
library(stringr)\
library(ggplot2)\
library(DESeq2)\
library(vsn)\
library(limma)\
library(ggplot2)\
library(ggrepel)\
library(enrichR)\
\
\
## Define functions\
get_Ratios <- function(a)\{\
  rownames(a) <- a$Accession\
  ratios <- a[, c(2:7,11:16)]\
  rownames(ratios) <- a$Accession\
  as.data.frame(ratios)\
\}\
\
## Set directory path\
\
setwd("/Users/noahallen/...")\
metadata_path <- "/Users/..."\
data_path <- "/Users/..."\
\
allfiles <- list.files(path = data_path)\
filenames <- str_split_fixed(allfiles, "_", Inf)\
filenames <- as.data.frame(filenames)\
tissue <- filenames$V7\
batch <- filenames$V8\
filepaths <- paste(data_path, allfiles, sep = "/")\
filepaths\
\
for(i in 1:length(filepaths)) \{   \
  assign(paste(tissue[i], batch[i]), read.delim(filepaths[i], sep = "\\t"))\
\}\
\
datasets <- as.data.frame(paste(tissue, batch))\
datasets\
\
\
##pick a data set\
\
my_tissue <- 'GST'\
x1 <- get(paste(my_tissue,'1'))\
x2 <- get(paste(my_tissue,'2'))\
all_protein_data <- merge(x=x1, y=x2, by="Accession")\
proteinColnames <- as.data.frame(list(colnames(all_protein_data)))\
my_data <- cbind(all_protein_data[c("Accession", "Ratios...126.....131..x","Ratios...127_N.....131..x","Ratios...127_C.....131..x","Ratios...128_N.....131..x",\
                                    "Ratios...128_C.....131..x","Ratios...129_N.....131..x","Ratios...129_C.....131..x","Ratios...130_N.....131..x","Ratios...130_C.....131..x",\
                                    "Ratios...126.....131..y","Ratios...127_N.....131..y","Ratios...127_C.....131..y","Ratios...128_N.....131..y","Ratios...128_C.....131..y",\
                                    "Ratios...129_N.....131..y","Ratios...129_C.....131..y","Ratios...130_N.....131..y","Ratios...130_C.....131..y")])\
\
\
### collect meta data regarding names and TMT labels\
protemoics_biobank <- read.delim(file = metadata_path, sep = "\\t")\
\
### Parse protein description and gene name\
description <- as.data.frame(all_protein_data$Description.x)\
description\
colnames(description) <- c("Description")\
a1 <- as.data.frame(str_split_fixed(description$Description, 'OS=',Inf))\
a2 <- as.data.frame(str_split_fixed(a1$V2,'GN=',Inf))\
a3 <- as.data.frame(str_split_fixed(a2$V2, 'PE=', Inf))\
names <- as.data.frame(gsub(" ","", a3$V1)) #remove whitespace after duplicate names\
colnames(names) <- c('V1')\
pdesc <- as.data.frame(cbind(a1$V1,names$`gsub(" ", "", pdesc$V2)`))\
\
###\
my_ratios <- get_Ratios(my_data)\
head(my_ratios, 5)\
\
### adjust column names\
cols <- as.data.frame(substr(colnames(my_ratios),10,14 ))\
cols1 <- as.data.frame(str_sub(colnames(my_ratios), -1,-1))\
cols\
cols2 <- as.data.frame(gsub("_","",cols$`substr(colnames(my_ratios), 10, 14)`))\
cols2\
cols2[1,1]=126\
cols2[7,1]=126\
cols2\
colnames(cols2) <- c("TMT_Label")\
batch_info <- c("1","1","1","1","1","1","2","2","2","2","2","2")\
cols2 <- paste(cols2$TMT_Label,batch_info, sep = "")\
cols2\
\
colnames(my_ratios) <- cols2\
\
\
##Get sample IDs and TMT labels\
IDS<- str_split_fixed(protemoics_biobank$Sample.Name, "_", Inf)\
IDS2 <- as.data.frame(IDS)\
Label<- as.data.frame(protemoics_biobank$Label)\
sample_info <- cbind(IDS2,Label)\
metadata_batchinfo <- c("1","1","1","2","2","2","1","1","1","2","2","2","na","na","na","na")\
sample_info <- cbind(sample_info, metadata_batchinfo)\
sample_info["ID"] <- cbind(paste(sample_info$`protemoics_biobank$Label`,sample_info$metadata_batchinfo, sep = ""))\
sample_info <- sample_info %>% filter(!row_number() %in% c(13:16))\
new_sampleID <- sample_info$ID\
\
#match experimental data to order of metadata , necessary for input into SummarizedExperiment\
my_ratios <- my_ratios[match(new_sampleID,colnames(my_ratios))]\
sample_info\
\
\
## Creating SE object\
expdesign <- data.frame(sample_info$ID,sample_info$V4,sample_info$V5)\
expdesign <- expdesign %>% filter(!row_number() %in% c(13:16))\
colnames(expdesign) <- c("label", "condition", "replicate")\
expdesign\
\
#add protein names and ID to ratio list, needed for SE generation\
my_ratios["name"] <- (names$V1)\
my_ratios["ID"] <- (rownames(my_ratios))\
head(my_ratios, 5)\
\
#remove duplicates and blank names\
my_clean_ratios <- my_ratios[!duplicated(my_ratios$name),]\
my_clean_ratios <- my_clean_ratios[!my_clean_ratios$name=="",]\
head(my_clean_ratios,5)\
\
#make SE object\
data_se <- make_se(my_clean_ratios, columns = c(1:12), expdesign = expdesign)\
data_se\
\
## Check frequencies of proteins\
protein_freq <- paste("Output/", my_tissue, " Protein Frequency.png", sep = "")\
png(protein_freq)\
plot_frequency(data_se)\
dev.off()\
\
## Filter for proteins that do show up in at least 5/6 of the replicates in at least one condition\
plot_missval(data_se)\
\
data_filt <- filter_missval(data_se, thr = 1)\
plot_missval(data_filt)\
\
MissingVal_prefilter <- paste("Output/", my_tissue, " missing values pre-filter.png", sep = "")\
png(MissingVal_prefilter)\
plot_missval(data_se)\
dev.off()\
\
MissingVal_postfilter <- paste("Output/", my_tissue, " missing values post-filter.png", sep = "")\
png(MissingVal_postfilter)\
plot_missval(data_filt)\
dev.off()\
\
## Normalization\
data_norm <- normalize_vsn(data_filt)\
\
# view SdPlots\
meanSdPlot(assay(data_filt))\
meanSdPlot(assay(data_norm))\
\
#view histogram of intensities\
hist(assay(data_filt), breaks = 1000, col = "#d95f0e")\
hist(assay(data_norm), breaks = 1000, col = "#d95f0e")\
\
#View and save un-normalized and normalized data\
normalized_plots <- paste("Output/", my_tissue, " Sample Normalization.png", sep = "")\
png(normalized_plots)\
plot_normalization(data_filt, data_norm)\
dev.off()\
\
## View and save PCA to observe batch effect before correction\
BE_pca_plot <- paste("Output/", my_tissue, " PCA before BE correction.png", sep = "")\
png(BE_pca_plot)\
plot_pca(data_norm)\
dev.off()\
\
## Remove batch effects between different mass-spec TMT runs\
norm_be_rm <- data.frame(assay(data_norm))\
colnames(norm_be_rm) <- colnames(my_clean_ratios)[1:(length(colnames(my_clean_ratios)) - 2)]\
batch_rm <- metadata_batchinfo[1:(length(colnames(norm_be_rm)))]\
norm_nobe = removeBatchEffect(norm_be_rm, batch = batch_rm)\
norm_nobe <- as.matrix(norm_nobe)\
assay(data_norm, withDimnames = F) <- norm_nobe #apply updated matrix to original SE object\
\
## View and save PCA to observe batch effect after correction\
rmBE_pca_plot <- paste("Output/", my_tissue, " PCA after BE correction.png", sep = "")\
png(rmBE_pca_plot)\
plot_pca(data_norm)\
dev.off()\
\
## View and save imputation plots, distribution of experimental groups\
batch_rm_plots <- paste("Output/", my_tissue, " Batch Effect Removal Plots.png", sep = "")\
png(batch_rm_plots)\
plot_imputation(data_se,data_norm)\
dev.off()\
\
## Shows data with/without missing values\
plot_detect(data_norm)\
\
\
## Imputation\
imputed_knn <- DEP::impute(data_norm, fun = "knn", k = 10, rowmax = 0.9)\
imputed_QRILC <- DEP::impute(data_norm, fun = "QRILC")\
\
plot_imputation(imputed_QRILC, imputed_knn)\
\
\
# DE analysis on imputated data\
\
differences <- test_diff(imputed_knn, type = "all", control = "GC")\
dep <- add_rejections(differences, alpha = 0.05, lfc = log2(1.5))\
get_results(dep)\
results <- get_results(dep)\
\
DEA_csv <- paste("Output/", my_tissue, "DEA.csv", sep = " ")\
write.csv(results, file = DEA_csv)\
\
plot_p_hist(dep)\
\
\
p <- ggplot(data =results , aes(x=FLT_vs_GC_ratio, y= -log10(FLT_vs_GC_p.adj),label= name)) +\
  geom_point() +\
  theme_DEP1() +\
  geom_vline(xintercept=c(-1, 1), col="red") +\
  geom_hline(yintercept=-log10(0.05), col="red")+\
  geom_text_repel() +\
  scale_color_manual()\
\
volcano_plot <- paste("Output/", my_tissue, " volcano plot.png", sep = "")\
png(volcano_plot)\
p\
dev.off()\
\
corr_plot <- paste("Output/", my_tissue, " correlation plot.png", sep = "")\
png(corr_plot)\
plot_cor(dep, significant = TRUE, lower = 0, upper = 1, pal = "Reds")\
dev.off()\
\
\
heatmap_knn <- plot_heatmap(dep, type = "centered", kmeans = TRUE, \
             k = 6, col_limit = 4, show_row_names = FALSE,\
             indicate = c("condition"), plot = F)\
heatmap_csv <- paste("Output/", my_tissue, "heatmap knn.csv", sep = " ")\
write.csv(heatmap_knn, file = heatmap_csv)\
\
heatmap_knn_plot <- paste("Output/", my_tissue, "heatmap knn.png", sep = "")\
png(heatmap_knn_plot)\
plot_heatmap(dep, type = "centered", kmeans = TRUE, \
             k = 6, col_limit = 4, show_row_names = FALSE,\
             indicate = c("condition"), plot = T)\
dev.off()\
\
\
#GSEA example using EnrichR\
\
dbs <- listEnrichrDbs() #list available databases\
\
ProteomicsDB <- test_gsea(dep, "GO_Biological_Process_2023")\
ProteomicsDB\
\
plotEnrich(ProteomicsDB, showTerms = 10)\
\
}