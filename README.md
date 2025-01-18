# Angptl4-in-NAFLD-immunopathologies

The project contains the R scripts of the bulk and single-cell transcriptomes of our diet-induced non-alcoholic fatty liver disease (NAFLD) mouse model. We formulated a specialized diet (Liver Disease Progression Aggravation Diet; LIDPAD) that can accelerate NAFLD progression in C57BL/6J mice when they are housed at thermoneutral condition (~ 30 degree Celcius). Manuscript 1 is under review by Advanced Science.

## Manuscript
Manuscript 1: Low, Zun Siong, et al. "The LIDPAD mouse model captures the multisystem interactions and extrahepatic complications in MASLD." Advanced Science 11.35 (2024): 2404326.
Link: https://advanced.onlinelibrary.wiley.com/doi/full/10.1002/advs.202404326

Manuscript 2: Angiopoietin-like 4 shapes the intrahepatic T-cell landscape via eIF2α signaling during steatohepatitis in diet-induced NAFLD.
Link: Pending

## Data
Bulk hepatic transcriptomes of LIDPAD-fed mice at different time points: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159911
Single-cell intrahepatic immune landscape during liver fibrosis: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214172
Bulk hepatic transcriptomes of LIDPAD-fed mice with LysMCre-Angptl4-knockout or treated with cAngptl4 monoclonal antibodies: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE214504

## R scripts
"LIDPAD model" folder contains the Rscripts for the following analyses:
  1. Bulk liver transcriptomes of LIDPAD- (and control diet-) fed mice at Week 1, 4, 8, 12, 16, 32 and 48 post feeding.
  2. The establishment of a human transcriptomic-guided disease staging approach.
  3. Liver transcriptomes of MASH mice given diet reversion.

"Liver immune landscape" folder contains the scripts for the following analyses:
  1. "scRNAseq" folder contains the scripts to annotate the intrahepatic immune clusters, identify the differentially expressed genes, functional enrichment and psedotime trajectory.
  2. "mAb treatment liver RNAseq" folder contains the Rscripts to perform differential expression analysis between the liver transcriptomes of control diet-fed mice, untreated, cAngptl4 mAb treated and LysMCre-Angptl4-knockout LIDPAD-fed mice.  
