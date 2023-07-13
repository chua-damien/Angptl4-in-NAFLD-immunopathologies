# Human transcriptomic-guided disease staging approach

## Step 1: Obtain transcriptomic features for NAFLD staging based on public human datasets

Scripts and data files are available in "human" folder. 

Human NAFLD liver transcriptomic datasets (GSE130970, GSE162694, GSE207310, GS225740; with NAS and fibrosis stage) were obtained from GEO.

The datasets were integrated by modelling the effects from batch and gender in DESeq2. For visualization, batch- and gender-adjusted data was obtained using limma::removeBatchEffect function.

To identify genes associated with disease severity, a differential gene expression analysis was performed using NAS score and fibrosis stage as covariates while controlling for the batch and gender effects (design=~batch+gender+ NAS score+Fibrosis stage). 

Genes with a log2 fold change > ±0.1 (per unit change of NAS score) or > ±0.2 (per unit change of fibrosis stage) with an adjusted p-value <0.05 were considered NAFLD-associated DEGs. Boruta feature selection algorithm was used to shortlist gene signatures.

## Step 2: Curate a unified human-mouse NAFLD gene signature

Scripts and data files are available in "mouse" folder. 

Manual curation was carried out to include only genes positively correlated with disease progression and with a mouse ortholog based on biomaRT. This analysis resulted in a human-mouse unified NAFLD transcriptomic signature comprised of 54 genes. 

## Step 3: Inter-study/diet-induced model comparison using ssGSEA and the unidied NAFLD gene signature

Scripts and data files are available in "ssgsea using unified features" folder. 

For comparative transcriptomic between different diet-induced NAFLD models, liver transcriptomes of mice treated with MCD (GSE156918, GSE144443), CDAHFD (GSE137449), HFD (GSE222171, GSE212294, GSE198358, GSE213355, GSE180353, GSE205021), ALIOS (GSE198358), AMLN (GSE195483), and GAN (GSE225616) for varying durations were obtained from GEO. The gene count data from these mouse datasets, LIDPAD model, and the liver transcriptomes of in-house NASH patient cohorts were subjected to variance stabilizing transformation (vsd). Genes with both the human and mouse orthologs were retained. ssGSEA was performed on all the samples using the unified NAFLD transcriptomic signature and the enrichment scores were ranked for disease stratification and mapping for mouse NAFLD severity to the human counterparts.
