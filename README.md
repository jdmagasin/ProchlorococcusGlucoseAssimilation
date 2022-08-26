# ProchlorococcusGlucoseAssimilation

This repository has analysis scripts used in the article *Differential timing for glucose assimilation in* Prochlorococcus *and coexistent microbial populations in the North Pacific Subtropical Gyre* by Muñoz-Marín et al., accepted for publication in *Microbiology Spectrum* in August 2022.


Microarray data is processed by the BASH script runAllScripts.sh which runs scripts/microarrays.R and scripts/egsea.R. However, you will
first have to obtain the raw microarray data files from NCBI GEO [series GSE154594](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154594).
The 16S rRNA data in data/16S_STABvida is processed interactively in R using scripts/analyze16S.R.


## *Prochlorococcus* gene expression changes in response to glucose -- microarray
Custom microarrays were used to detect *Prochlorococcus* gene expression changes after 4, 12, and 24 hour incubations with glucose.
The core script for this analysis is microarrays.R, which handles all processing steps: quality control checks, background subtraction,
probe normalization by quantiles, probe to gene conversion by median polishing, gene detection, and the identification of significantly
differentially expressed genes in glucose treatments versus controls.  The script depends heavily upon a [software pipeline](https://www.jzehrlab.com/microtools)
that was created by Jonathan Zehr's lab at UC Santa Cruz for the MicroTOOLs microarray design 
([Shilova et al., 2014](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4069398/)).  To use the pipeline for the present study,
the [*Prochlorococcus* microarray description](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL28884)
was substituted for the [MicroTOOLs description](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL24371).  Instructions for installing
the MicroTOOLs software pipeline (R package) are available from the provided link.


## Ensemble of Gene Set Enrichment Analyses (EGSEA) -- microarray
In addition to single-gene differential expression analysis, we looked for pathways that responded to glucose treatments using
[EGSEA](https://f1000research.com/articles/6-2010/v1).  EGSEA uses the consensus of up to 12 gene set enrichment analysis algorithms
to identify sets of genes that in aggregate change significantly.  To use EGSEA we defined pathways to include genes that increase 
(or decrease) together as the pathway is up- (or down-) regulated in response to glucose.  The published results define pathways
by ecotype.  However, egsea.R also supports EGSEA for each *Prochlorococcus* strain.


## Microbial community structure changes in response to glucose -- 16S rRNA sequences
Raw MiSeq sequences for V3 and V4 regions of 16S rRNA genes were processed using a QIME2/DADA2-based pipeline by a commericial
bioinformatics company STABvida.  However, the analysis described in the paper is mainly based on analyze16S.R.
