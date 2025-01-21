### This pipeline enables to predict eRNA-promoter interactions via "late" integration
### of three types of omics data: ChIP-Seq, Hi-C contact matrices, and RNA-Seq counts
### The goal is to predict eRNA-induced gene expression changes followed by binding
### of transcription factor of interest to cis-regulatory elements (CRE) located
### within intron sequences or intergenic regions

### List of Dependencies
## ChIPseeker
## dplyr
## GenomicRanges
## liftOver
## org.Hs.eg.db
## purrr
## readr
## rtracklayer
## Signac
## spatstat.core
## TxDb.Hsapiens.UCSC.hg19.knownGene

## function for annotating imported ChIP-seq peaks in .bed format
Peaks_annotation <- function(peaks, TSS_start = -1000, TSS_end = 1000, mode = "Introns") {
  require(ChIPseeker)
  require(dplyr)
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  
  peakAnno <- annotatePeak(peaks, tssRegion=c(TSS_start, TSS_end),
                           TxDb=txdb, annoDb="org.Hs.eg.db")
  ranges <- as.data.frame(peakAnno)
  
  ## stripping off intron number as this info not needed in this context
  for (i in 1:nrow(ranges)) {
    if (startsWith(ranges$annotation[i], "Intron") == TRUE) {
      ranges$annotation[i] <- "Intron" }
  }
  
  ## filtering Intronic or Intergenic peaks only
  if (mode == "Introns") {
    filtered_ranges <- ranges %>%
      filter(annotation == "Intron")
  }
  else if (mode == "Intergenic") {
    filtered_ranges <- ranges %>%
      filter(annotation == "Distal Intergenic")
  }
  else {
    print("Wrong annotation mode! Options: Introns or Intergenic") 
    return(NA)
  }

  filtered_ranges
}

## a series of functions for lifting over between different genome versions
## use UCSC chain files for lifting
lift_genome_hg38_to_hg19 <- function(original_version = "hg38", final_version = "hg19", ranges) {
  require(liftOver)
  ch = import.chain(paste0(getwd(), "/hg38ToHg19.over.chain"))
  ranges_lifted = unlist(liftOver(ranges, ch))
  genome(ranges_lifted) = "hg19"
  ranges_lifted
}

lift_genome_hg19_to_hg18 <- function(original_version = "hg19", final_version = "hg18", ranges) {
  require(liftOver)
  ch = import.chain(paste0(getwd(), "/hg19ToHg18.over.chain"))
  ranges_lifted = unlist(liftOver(ranges, ch))
  genome(ranges_lifted) = "hg18"
  ranges_lifted
}

lift_genome_hg18_to_hg19 <- function(original_version = "hg18", final_version = "hg19", ranges) {
  require(liftOver)
  ch = import.chain(paste0(getwd(), "/hg18ToHg19.over.chain"))
  ranges_lifted = unlist(liftOver(ranges, ch))
  genome(ranges_lifted) = "hg19"
  ranges_lifted
}

lift_genome_hg19_to_hg38 <- function(original_version = "hg19", final_version = "hg38", ranges) {
  require(liftOver)
  ch = import.chain(paste0(getwd(), "/hg19ToHg38.over.chain"))
  ranges_lifted = unlist(liftOver(ranges, ch))
  genome(ranges_lifted) = "hg38"
  ranges_lifted
}

## functions for converting gene names 
ENTREZtoSYMBOL <- function(genes = character()) {
  require("org.Hs.eg.db")
  mapIds <- mapIds(org.Hs.eg.db, keys = genes, keytype = "ENTREZID", column="SYMBOL")
  mapIds_df <- as.data.frame(mapIds)
  mapIds_df$mapIds
}

ENSEMBLtoSYMBOL <- function(genes = character()) {
  require("org.Hs.eg.db")
  mapIds <- mapIds(org.Hs.eg.db, keys = genes, keytype = "ENSEMBL", column="SYMBOL")
  mapIds_df <- as.data.frame(mapIds)
  mapIds_df$mapIds
}

Process_HiC <- function(HiC_data, first_bin = 10, last_bin = 25) {
  ## preparing GRanges from Hi-C data bins
  HiC_data <- HiC_data %>%
    mutate(ranges1 = paste0(Chr, "_", Region1)) %>%
    mutate(ranges2 = paste0(Chr, "_", Region2))
  
  ## getting top-positive scores (n = 3) 
  df <- split(HiC_data, HiC_data$ranges1)
  ## filtering 10-25 bins
  df_filt_dist <- list()
  for (i in 1:length(df)) {
    df_filt_dist[[i]] <- df[[i]][first_bin:last_bin,]
  }
  names(df_filt_dist) <- names(df)
  
  df_ordered_filt <- lapply(df_filt_dist, function(x) x[order(x[, "score"], decreasing = TRUE),]) 
  HiC_top3 <- lapply(df_ordered_filt, function(x) x[1:3,]) 
  
  HiC_top3

}


## function to integrate ChIP-Seq, Hi-C & RNA-Seq data
## Input1: annotated ChIP-Seq peaks (Peaks_annotation function);
## Input2: eRNA .bed regions from database of choice;
## Input3: processed binned Hi-C data from "Process_HiC" function;
## Input4: default Hi-C genome version = "hg18" based on published Hi-C matrix being used [GSM1081530]
## Note: all other regions have to be in hg19; please use provided chain functions to convert
## this function outputs predicted eRNA-gene interactions in a table format: eRNA's genomic coordinates (hg18) and corresponding differentially-expressed gene
Integrate_3DSeq <- function(anno_peaks, eRNA, HiC_top3, HiC_version = "hg18") {
  require(GenomicRanges)
  library(spatstat.core)
  library(Signac)
  
  ## find overlaps between ChIP-seq peaks and eRNA ranges
  anno_peaks_ranges <- unlist(makeGRangesListFromDataFrame(anno_peaks, keep.extra.columns = F))
  eRNA_ranges <- GRanges(eRNA)
  ol_peaks_eRNA <- subsetByOverlaps(anno_peaks_ranges, eRNA_ranges, minoverlap = 1, ignore.strand = TRUE)
  
  ## creating granges from ranges1 (=region1) of HiC data 
  gr1 <- StringToGRanges(names(HiC_top3), sep = c("_", "_"))
  
  if (HiC_version == "hg18") {
    require(liftOver)
    ch = import.chain("~/hg19ToHg18.over.chain")
    ol_peaks_eRNA_lifted = unlist(liftOver(ol_peaks_eRNA, ch))
    genome(ol_peaks_eRNA_lifted) = "hg18"
  }
  
  ## find overlaps between ChIP-seq/eRNA peaks and Hi-C ranges (in hg18)
  ol_HiC <- suppressWarnings(findOverlaps(ol_peaks_eRNA_lifted, gr1, type = "within", ignore.strand = TRUE)) ### query is wholly contained within the HiC region
  rownumber_gr1 <- ol_HiC@to ## saving the number of HiC region containing overlap
  number_eRNA <- ol_HiC@from
  
  ## connecting to Region2, merging to granges (gr2_final contains most probable regions of eRNA-interactions)
  gr2_high_score <- list()
  for (i in 1:length(rownumber_gr1)){
    gr2_high_score[[i]] <- HiC_top3[[rownumber_gr1[i]]]
  }
  
  library(purrr)
  library(dplyr)
  gr2 <- list_rbind(gr2_high_score)
  gr2[,"eRNA_number"] <- rep(number_eRNA, each = 3) ## saving corresponding eRNA number from "ol_peaks_eRNA_lifted" object

  gr2 <- gr2 %>% filter(score != 0) ## keeping only meaningful interactions (score > 0)
  gr2_regions <- StringToGRanges(gr2$ranges2, sep = c("_", "_"))
  
  ## lifting gr2_ranges to hg19 back to be able to overlap with DEGs
  if (HiC_version == "hg18") {
    require(liftOver)
    ch = import.chain("~/hg18ToHg19.over.chain")
    gr2_regions_hg19 = unlist(liftOver(gr2_regions, ch))
    genome(gr2_regions_hg19) = "hg19"
  }
  
  ## annotating final regions with gene promoters and/or gene 5'-UTRs as we expect eRNA affect gene promoters
  
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
  library(ChIPseeker)
  regions_Anno <- suppressWarnings(annotatePeak(gr2_regions_hg19, tssRegion=c(-1000, 1000), TxDb=txdb,
                               level = "gene",
                               genomicAnnotationPriority = c("Promoter", "5UTR"),
                               overlap = "TSS", verbose =T))
  regions_Anno_df <-as.data.frame(regions_Anno)
  
  extract_df <- regions_Anno_df %>%
    filter(!is.na(annotation))
  
  granges_extract_df <- makeGRangesFromDataFrame(extract_df, keep.extra.columns = T)
  if (HiC_version == "hg18") {
    require(liftOver)
    ch = import.chain("~/hg19ToHg18.over.chain")
    granges_extract_df_hg18 = unlist(liftOver(granges_extract_df, ch))
    genome(granges_extract_df_hg18) = "hg18"
  }
  
  ### cleaning data and inner joining to get a table of eRNA (hg18)-gene interactions 
  eRNA_result <- gr2 %>% dplyr::select(score, ranges1, ranges2, eRNA_number)
  genes_result <- as.data.frame(granges_extract_df_hg18, row.names = NULL) 
  genes_result <- genes_result %>% dplyr::select(seqnames, start, end, width, 
                                                 strand, annotation, geneId)
  
  genes_result <- genes_result %>%
    dplyr::mutate(ranges2 = paste0(seqnames, "_", start, "_", end)) %>% 
    dplyr::select(ranges2, geneId, annotation) %>% 
    dplyr::distinct(geneId, .keep_all = T)
  
  merged_data <- inner_join(genes_result, eRNA_result, by = "ranges2")
  merged_data <- merged_data %>% dplyr::select(ranges1, ranges2, score, eRNA_number, annotation, geneId)
  
  ol_peaks_eRNA_lifted <- as.data.frame(ol_peaks_eRNA_lifted, row.names = NULL)
  ol_peaks_eRNA_lifted <- ol_peaks_eRNA_lifted %>% 
    dplyr::mutate(eRNA_number = 1:nrow(ol_peaks_eRNA_lifted),
           ranges_eRNA = paste0(seqnames, "_", start, "_", end)) %>% 
    dplyr::select(ranges_eRNA, eRNA_number)
  
  merged_data2 <- inner_join(merged_data, ol_peaks_eRNA_lifted, by = "eRNA_number")
  
  merged_data2 <- merged_data2 %>% dplyr::select(ranges_eRNA,score, geneId)
  merged_data2
  
}


##### START #####

## Import
## Importing Chip-Seq peaks, annotating them and saving as GRanges object
require(rtracklayer)
require(readr)

extraCols_Peak <- c(FoldChange="numeric", pVal="numeric", qVal="numeric", peak = "numeric")
Chip_peaks <- import.bed("Chtop_rep1_peaks_hg38.narrowPeak", extraCols = extraCols_Peak)

### please make sure to convert peaks to hg19 version if they are hg38
Chip_peaks_hg19 <- lift_genome_hg38_to_hg19(ranges = Chip_peaks)

## peak annotation. There are two modes - peaks that locate within the introns or peaks that locate in intergenic regions
## Promoters are defined as (TSS_start, TSS_end); promoters are (-1000, 1000) by default
## In the output, only peaks mapped either to introns (mode = "Introns", default) or intergenic regions (mode = "Intergenic") are filtered
anno_peaks <- suppressWarnings(Peaks_annotation(Chip_peaks_hg19, mode = "Introns"))

## Importing eRNA regions in .bed format, for example, from Enhancer Atlas 2.0 database:
url <- "http://www.enhanceratlas.org/data/download/enhancer/hs/HEK293T.bed" ## these are hg19 regions 
download.file(url, destfile = "./eRNA_HEK293.bed",method = "auto")
eRNA <- import.bed("eRNA_HEK293.bed") 

## Importing Hi-C processed contact matrices corresponding to your experimental design and/or cell line as a tab file.
HiC_data <- as.data.frame(read_delim("Hi-C_HEK293_GSM1081530_CTRL_r1_cis.index_hg18.txt")) ## this is a published Hi-C matrix of 40 kb bins in HEK293 cells (GEO record #GSM1081530, hg18)
colnames(HiC_data) <- c("Chr", "Region1", "Region2", "score")
head(HiC_data)
dim(HiC_data)

## Importing DEGs (e.g., downregulated) lists of genes derived via analysis of RNA-Seq counts (using DESEq2/salmon etc)
Downreg_DEGs <- read_table("DEGs-down_Chtop_KD.tabular")

## function to process Hi-C data: converting to GRanges, then
## getting top-positive scores between regions that are separated by user-specified bins (by default, if 1 bin = 40 kb, it's bin10-bin25 distance, e.g. 400-1000 kb) 
## output is filtered Hi-C matrix of interactions that contains every bin + top-3 interacting regions for it
HiC_processed <- Process_HiC(HiC_data, first_bin = 10, last_bin = 25)


## Running core 3D-Seq function: please make sure Chip-seq peaks and eRNA regions are in hg19 version
predictions <- Integrate_3DSeq(anno_peaks = anno_peaks, eRNA = eRNA, 
                               HiC_top3 = HiC_processed, HiC_version = "hg18")

## the result of this function is a table of eRNA regions (hg18) and genes that can be potentially controlled
## checking if there are downregulated genes there (first converting gene names to SYMBOLs)

predictions[,"SYMBOL"] <- ENTREZtoSYMBOL(genes = predictions$geneId)
symbols_down_DEGs <- ENSEMBLtoSYMBOL(genes = Downreg_DEGs$Chtop_genes_ens_down)

final_pred <- predictions %>% 
  filter(SYMBOL %in% symbols_down_DEGs) %>% 
  arrange(desc(score)) ### these are the final genes can be worth following up

## converting ranges to hg19 -> hg38
final_pred_ranges <- StringToGRanges(final_pred$ranges_eRNA, sep = c("_", "_"))

final_pred_ranges_hg19 <- lift_genome_hg18_to_hg19(ranges = final_pred_ranges)
final_pred_ranges_hg38 <- lift_genome_hg19_to_hg38(ranges = final_pred_ranges_hg19)

final_pred_ranges_hg38_df <- as.data.frame(final_pred_ranges_hg38, row.names = NULL)

### final predictions: 
final_pred ### eRNA ranges are in hg18
final_pred_ranges_hg38_df ### here are corresponding hg38 ranges

### The most probable interactions will have high score of Hi-C interactions
### printing out top-5 predictions

final_pred[1:5,]
final_pred_ranges_hg38_df[1:5,]

