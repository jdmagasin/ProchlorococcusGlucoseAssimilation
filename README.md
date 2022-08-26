# ProchlorococcusGlucoseAssimilation

This repository has analysis scripts used in the article *Differential timing for glucose assimilation in Prochlorococcus and coexistent microbial populations in the North Pacific Subtropical Gyre* accepted for publication in *Microbiology Spectrum* in August 2022.


### *Prochlorococcus* gene expression changes in response to glucose -- microarray
Custom microarrays were used to detect *Prochlorococcus* gene expression changes after 4, 12, and 24 hour incubations with glucose.
The core script for this analysis is microarrays.R, which handles all processing steps: quality control checks, background subtraction,
probe normalization by quantiles, probe to gene conversion by median polishing, gene detection, and the identification of significantly
differentially expressed genes in glucose treatments versus controls.  The script depends entirely upon a [software pipeline](https://www.jzehrlab.com/microtools)
that was created by Jonathan Zehr's lab at UC Santa Cruz for the MicroTOOLs microarray design.  To use the pipeline for the present study, 
the [*Prochlorococcus* microarray description](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL28884)
was substituted for the [MicroTOOLs description](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL24371).


### Microbial community structure changes in response to glucose -- 16S rRNA sequencing
Raw MiSeq sequnces were processed using a DADA2-based pipeline by a bioinformatics company. However, the analysis
presented in the paper is mainly based on the script analyze16S.R.


## Ensemble of Gene Set Enrichment Analyses (EGSEA) -- microarray
In addition to single-gene differential expression analysis, we looked for pathways that responded to glucose.  Pathways
were defined limited to genes on the microarray, and only used genes that tend to change in the same direction (increase
together, or decrease together) as a pathway is up- or down-regulated.
