# 3D-Seq Integration
This workflow allows integration of three types of omics data to predict eRNA-gene promoter interactions:
* ChIP-Seq peaks;
* Hi-C contact matrices;
* RNA-Seq (lists of differentially-expressed genes).
The workflow is written in R and includes the list of dependencies at the top of the script. To test the pipeline, please refer to exemplary open datasets provided in 'data' (ChIP-Seq peaks of CHTOP, a component of TRanscription and EXport complex <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130992>, Hi-C data in HEK293 cells <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM1081530>, list of pre-determined differentially expressed genes in CHTOP-knockdown cells <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113952>).

# Rationale
Cis-regulatory elements (CREs) control gene expression, orchestrating tissue identity, developmental timing and stimulus responses, which collectively define the thousands of unique cell types in the body <https://www.nature.com/articles/s41586-024-08070-z>. Enhancers (eRNAs) are CREs in the genome that can loop over long genomic regions to engage distant promoters to control target gene transcription <https://pmc.ncbi.nlm.nih.gov/articles/PMC6971116/>. Thus, the question of specific eRNA-gene promoters interactions that act in a temporal and spatial manner attracts considerable attention. The workflow presented here aims to predict eRNA-gene interactions facilitated by a specific transcription factor (TF) of interest using this TF's binding locations, eRNA genomic locations, chromatin conformation data (Hi-C), and the changes in gene expression upon genetic or therapeutic inhibition of this TF (RNA-Seq). Basically, if TF can bind eRNA regions and promote their expression, we can predict that and distant target genes using the workflow.

# Usage
Typical steps of the workflow include:
* Peak annotation with two modes - "Introns" (default) or "Intergenic": *Peaks_annotation()* function. Only peaks annotated to Introns or Intergenic regions are outputed. **Note** This function expects peaks in hg19 version of human genome, so please use provided *lift_genome* functions to convert if peaks are in another version.
* Hi-C data processing: *Process_HiC()* function. This function uses Hi-C contact information in tabular format, and outputs top-3 interacting regions for every genomic bin based on proximity score. User can specify the range of distance between regions, e.g. 10-25 bins with 40 kb resolution = 400-1000 kb. **Note** Automated conversion of other Hi-C formats (e.g. cool) has not yet been implemented.
* Integration of ChIP-Seq annotated peaks, eRNA regions, and Hi-C top-interacting regions that outputs a table of predicted eRNA (genomic coordinates, hg18)-gene interactions: *Integrate_3DSeq()* function. After that, we intersect the predicted genes with actually observed as downregulated/upregulated upon TF's perturbation, and the overlapping genes are worth following up. 
* Final output prints out top-5 interactions (eRNA regions in hg38 and corresponding genes).

# Suggested validations and visualizations
Predicted eRNA-gene interactions can be further validated using nascent transcriptome sequencing technologies such as TT-Seq, PRO-Seq, or GRO-Seq. Additional visualizations include:
* IGV/IGB tracks: in a genome viewer of choice, ChIP-Seq of TF of interest and underlying TT-Seq bigwig coverage (Control vs knockdown/overexpression of TF) on eRNA regions predicted to be driving specific gene expression;
* Hi-C data can be shown using a heatmap of interactions between eRNAs and gene promoters.