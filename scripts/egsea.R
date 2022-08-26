#!/usr/bin/env Rscript

##
## Adapted EGSEA code from used for NEMO (Shilova et al. 2020;  doi: 10.1371/journal.pone.0231771)
##     /Volumes/JMagasin_SD1/Archive/Users/jmagasin/Desktop/Researchy/ZehrLabMicroarray/NEMO2/
##        Make_NEMO_Figures.Rmd
##        Manuscript_Drafting/NEMO_NE7_GeneSets.Rmd
##

cat("\n\n################## Starting the EGSEA script ##################\n\n")

library(MicroTOOLs)
library(EGSEA)
library(knitr)
library(kableExtra)
library(dplyr)       # For mutate()
load('microarrays_script.RData')  # This is creaed by microarrays.R


## These impact which EGSEA results we show (not how EGSEA is run).
PVALCUT = 0.01        # from NEMO
LOG2FC  = log2(1.5)   # used 1.2 in NEMO
cat("EGSEA will use adjusted p-value cut off <",PVALCUT,"to determine significant (DE) gene sets.\n")
cat("Additionally, this script will point out gene sets that are DE and also have an average fold",
    "change that is >",round(2^LOG2FC,2),".\n\n")

## For repeatable random sampling (if necessary) by EGSEA alorithms
cat("Setting random number seed.\n")
set.seed(123)


## Define gene sets.  Key point is that all these genes must change in the same
## direction in response to an environmental (or other) stimulus.  Mari reviewed
## my original set against the literature to define the following (see email
## April 11 at 9:45am).
##
## May 4:
## - Drop edd per Mari email on Apr 23 because it might add noise b/c it is part
##   of both PPP and Entner-Dudoroff.  This is follow-up on earlier email where
##   I propose this "rule" for defining gene sets (draft for how it would appear
##   in the ms:
##        Our goal with EGSEA was to compare ecotype responses with respect to
##        key pathways.  For each ecotype and pathway, we defined the gene set
##        to include the genes that are widely accepted (across strains...) to
##        function in the pathway. Genes that function in multiple, unrelated
##        (uncoupled? Independent?) pathways were usually excluded.
##
## April 21 changes in response to Mari email:
## - Moved edd from PentosePhosphate to EntnerDudoroff.
## - Dropped ugd from PentosePhosphate.  It "is not considered part of the hub"
##   of PPP and it functions in amino sugar and amino nucleotide metabolism.
## Earlier changes:
## - Dropped glk from Entner-Dudoroff because uncertain whether the glk on the array are
##   in Pro's E-D
## - Removed prk from CarbonFixation b/c observed opposite to rbcS/L
## - KaiB function in Pro is uncertain.  First separated kaiB from kaiC into its own
##   CircadianRhythmKaiB, but then dropped it since not signif.  [Note that Zinser et.
##   al (2009) saw both kaiB and kaiC to be signif periodic.]
##
Pathways <-
    list(CellDivision                = c('ftsZ','minC','minD','minE'), # Mari checked: all should go in same direction
         CircadianRhythmKaiC         = c('kaiC'),
        #CircadianRhythmKaiB         = c('kaiB'),
         EnergyMetabolism            = c('atpA','atpB','atpC','atpE','atpF','atpI'),
         PhotosystemI                = c('psaA','psaB'),
         CarbonFixation              = c('rbcL','rbcS','rbcS/L'),
         Respiration                 = c('ctaC/coxB','ctaE/coxC','cyoB/coxA','cyoC'),
         PentosePhosphate            = c('opcA','rbsK','tal','zwf'),
         Glycolysis                  = c('pykF','glk'),
         Glycolysis.KrebsCycle       = c('pdhA','pdhB','lpdA'),
         KrebsCycleAndRelated        = c('fumC','ppc','sdhB'),
         GlycogenMetabolism          = c('glgC'),
         Glycolysis.Gluconeogenesis  = c('eno','fda','pgi','cbbA','fda/cbbA'),
         EntnerDudoroff              = c('eda','gdh'),
         PyruvateMetabolism          = c('dld','ddA','lldD'),
         SugarTransporter            = c('glcH'),
         AminoAcidMetabolism         = c('argH','glmS'),
         NitrogenMetabolism          = c('glnA','glsF'))

## Make sure the Pathways agree with annotation recorded in Fig S4 (annotTab
## created by microarrays.R).  My names in Pathways are not necessarily as in
## annotTab, but we can at least make sure that each set of genes in Pathways is
## from just one pathway.
for (pnam in names(Pathways)) {
    pgenes <- Pathways[[pnam]]
    x <- table(subset(annotTab, Gene.symbol %in% pgenes)$Pathway)
    x <- x[x>=1] #
    if (length(x) != 1) {
        warning("WARNING:  ", pnam," in Pathways has genes from multiple pathways as annotated in Fig S4 ",
             "specifically these pathways: ", paste(names(x),collapse=', '))
    }
}
rm(pnam,pgenes,x)

## This snippet was used to get the genes for the original gene set, which Mari
## subsequently adjusted based on her expertise + literature review.
if (FALSE) {
    ## What pathways were detected.
    for (pw in names(sort(table(fData(eset)$Pathway, useNA='always'), decreasing=T))) {
        cat(paste0("\nPathway '",pw,"'"))
        print(table(subset(fData(eset), Pathway==pw)$Gene.symbol, useNA='always'))
    }
}


##
## Define taxa that will be used in combination with pathways to create the gene
## sets that EGSEA tests.
##

## Ecotype is an option, but clade is more interesting (to show whether
## consistent patterns within ecotypes).
##Taxa <- list(pro.HL          = list(Ecotype=c('HLunk','HLI','HLII')),
##             pro.LL          = list(Ecotype=c('LLI','LLII/III','LLIV')))
Taxa <- list(pro.HLI         = list(Ecotype='HLI'),
             pro.HLII        = list(Ecotype='HLII'),
             pro.HLunk       = list(Ecotype='HLunk'),  # defined by microarrays.R
             pro.LLI         = list(Ecotype='LLI'),
             pro.LLII.III    = list(Ecotype='LLII/III'),
             pro.LLIV        = list(Ecotype='LLIV'))

## If include these, then there are many more gene sets for EGSEA to test. This
## results in fewer DE pathways, almost certainly because more gene sets --->
## higher adjusted [for multiple hypotheses] p-values.  The strains we included
## only to see if consistent with the higher order (ecotype and clade) gene
## sets, but for the paper we can/should instead check DE genes from specific
## strains.  So I am not includeing Taxa.strains in the final 'Taxa' used to
## make GeneSets.
Taxa.strains <- list(
             ## HLI
             pro.MED4        = list(Strain='MED4'),
             pro.MIT9515     = list(Strain='MIT9515'),
             ## HLII
             pro.AS9606      = list(Strain='AS9606'),
             pro.MIT0604     = list(Strain='MIT0604'),
             pro.MIT9215     = list(Strain='MIT9215'),
             pro.MIT9301     = list(Strain='MIT9301'),
             pro.MIT9302     = list(Strain='MIT9302'),
             pro.MIT9312     = list(Strain='MIT9312'),
             ## HL other (likelyHLstrains defined in microarrays.R)
             pro.MIT9292     = list(Strain='MIT9292'),
             pro.HOT208_60m  = list(Strain=likelyHLstrains[grep('HOT208_60m',likelyHLstrains)]),
             
             ## LLI
             pro.NATL1       = list(Strain='NATL1'),
             pro.NATL2       = list(Strain='NATL2'),
             pro.PAC1        = list(Strain='PAC1'),
             ## LLII/III
             pro.MIT0601     = list(Strain='MIT0601'),
             pro.MIT0602     = list(Strain='MIT0602'),
             pro.MIT0603     = list(Strain='MIT0603'),
             pro.MIT9211     = list(Strain='MIT9211'),
             pro.SS120       = list(Strain='SS120'),
             ## LLIV
             pro.MIT9303     = list(Strain='MIT9303'),
             pro.MIT9313     = list(Strain='MIT9313'))


## Define the gene sets. Note that a set is for a specific taxon and pathway.
GeneSets <- list()
for (pw in names(Pathways)) {
    for (tx in names(Taxa)) {
      	targets = GenesInFeatureSet(eset, append(Taxa[[tx]], list(Gene.symbol=Pathways[[pw]])))
    	  if (length(targets) > 0) {
    	      ## E.g. Pro lacks dmdA so don't add that gene set.
    	      GeneSets[[paste0(tx,'.',pw)]] = targets
    	  }
    }
}
## Drop any sets that have just one gene.  CAMERA probably won't find them DE
## anyway, but I don't want them cluttering up the table that shows all the
## sets.
GeneSets <- GeneSets[names(which(sapply(GeneSets, length) > 1))]

## Describe all the gene sets (counts for each detected gene).
cat("The gene sets used for EGSEA are defined for each taxonomic group and pathway.\n",
    "Here are the definitions. For each taxon the count of _detected_ genes is given.\n\n")
print(lapply(GeneSets, function(x)  table(fData(eset)[x, 'Gene.symbol'])))
cat("\n\n")

##
## This is adapted from NEMO_NE7_GeneSets.Rmd.
##
BigFunctionToDoEGSEA = function(eset, geneSets, contrastMatrix, maxFDR=0.05)
{
    geneIdsInOrder = rownames(eset)
    geneAnnot = fData(eset)[geneIdsInOrder, 'Pathway', drop=F]
    geneAnnot$GeneSet = NA
    geneAnnot$ID = NA
    for (gs in names(geneSets)) {
      geneAnnot[ geneSets[[gs]], 'GeneSet' ] = gs
      geneAnnot[ geneSets[[gs]], 'ID' ]      = gs  # Need ID or crash during heatmaps
    }
    ## Reorder for report
    geneAnnot = geneAnnot[,c('GeneSet','Pathway','ID')]
    gs.annot = buildCustomIdx(geneIDs = geneIdsInOrder,
                              gsets = geneSets,
                              anno = geneAnnot,
                              name = "ProchlorococcusArrayGeneSets",
                              species = "Prochlorococcus",
                              label='ProchlorococcusGeneSets')

    ## Hack. Need the design from when we did the DE analysis.  Easy enough to
    ## recreate b/c 'array2condition' still exists.  But first check agreement
    ## with the eset (e.g. arrays could have been dropped).
    stopifnot(colnames(eset) == names(array2condition))
    fac <- factor(array2condition)
    design <- model.matrix(~0 + fac)
    colnames(design) <- levels(fac)  # design's cols ordered; this works
    
    ## I want to use egsea() rather than egsea.ma(), even though ".ma" is for microarrays, b/c
    ## ".ma" looks like it wants Entrez Gene IDs (see code) which I do not have.  Luckily there
    ## is vooma() which makes a voom object for microarray data (rather than RNA-seq data).
    v = vooma(eset, design)
    stopifnot(rownames(v$E) == geneIdsInOrder)
    results = egsea(v,
                    contrasts = contrastMatrix,
                    gs.annots = gs.annot,
                    display.top = 100,
                    fdr.cutoff = maxFDR, # 0.05 is the default.
                    report = F)  # See just below.
    return( results )
}


### Take camResults [sic] and a contrast name, get the data frame, and make a
### simple matrix representation of it.  stat is NA if you want CAMERA's
### Direction for the set.  Or you can specify one of the statistics I added
### such as 'meanLogFC'.  You can set bigLogFC to a log2 threshold.
### If absolute meanLogFC is bigger than this (and stat==NA) then
### the "Up" or "Down" will get a "*" appended.
### Also, if EGSEA results are available, check if they do not support the
### CAMERA result.  See Does.EGSEA.NOT.Support.CAMERA()
EGSEAResultsToMatrix = function(egseaResults, contrast, bigAvgLogfcDir=NA, pCut=1)
{
    contrast = gsub('\\.vs\\.', '-', contrast)  # Style change, if needed
    top = subset(topSets(egseaResults, contrast=contrast, number=Inf, names.only=F),
                 p.adj < pCut)
    rows = names(Taxa)
    cols = names(Pathways)
    m = matrix('', nrow = length(rows), ncol = length(cols), dimnames = list(rows,cols))
    for (tx in rows) {
        txIdx = grep(paste0(tx,'\\.'), rownames(top))   # which sets are for this taxon?
        ##if (length(txIdx) == 0) { next }              # No! Use inner loop for '-' (vs. '')
        for (sg in cols) {
            sgIdx = grep(paste0('\\.',sg), rownames(top)) # which sets are for this gene category?
            idx = intersect(txIdx,sgIdx)
            nam = paste0(tx,'.',sg)
            if (length(idx) > 0) {
                v = top[idx, 'direction']
                if (!is.na(bigAvgLogfcDir) &&
                    abs(top[idx, 'avg.logfc.dir']) >= bigAvgLogfcDir) {
                    v = paste0(v,"*")
                }
                m[tx,sg] = v
            } else if (nam %in% names(GeneSets)) {
                m[tx,sg] = '-'  # Set exists but is not DE
            }
        }
    }
    t(m)  # Did not transpose in NEMO. Better here b/c more pathways (with long names) than strains.
}

## Replace strings 'Up' and 'Down' with arrows, solid or thin depending on whether
## 'Up*' or 'Up' (and same for 'Down*' and 'Down').
## Also make row and column labels nicer.
HTMLifyUpDownMatrix <- function(m, caption=NULL)
{
    m[m=='Up*']   = '&#x2B06;'  # HTML unicode solid up arrow
    m[m=='Down*'] = '&#x2B07;'  #      unicode solid down arrow
    m[m=='Up']    = '&#x2191;'  #      unicode thin up arrow
    m[m=='Down']  = '&#x2193;'  #      unicode thin down arrow
    m[m=='-']     = '&#x25CB;'  #      unicode white circle

    x <- rownames(m)
    x <- gsub('EntnerDudoroff','Entner-Dudoroff',x)
    x <- gsub('PentosePhosphate','Pentose phosphate',x)
    x <- gsub('CellDivision','Cell division',x)
    x <- gsub('EnergyMetabolism','Energy metab.',x)
    ##x <- gsub('GlycolysisPyruvateMetabolism', 'Glycolysis / pyruvate metab.',x)
    x <- gsub('PyruvateMetabolism', 'Pyruvate metab.',x)
    x <- gsub('CarbonFixation','Carbon fixation',x)
    x <- gsub('KrebsCycleAndRelated','Krebs cycle and related',x)
    x <- gsub('CircadianRhythm','Circadian rhythm',x)
    x <- gsub('AminoAcidMetabolism','Amino acid metab.',x)
    x <- gsub('NitrogenMetabolism','Nitrogen metab.',x)
    x <- gsub('GlycogenMetabolism','Glycogen metab.',x)
    x <- gsub('SugarTransporter','Sugar transporter',x)
    x <- gsub('Glycolysis.Gluconeogenesis','Glycolysis / gluconeogenesis',x)
    x <- gsub('Glycolysis.KrebsCycle','Glycolysis / Krebs cycle',x)
    x <- gsub('PhotosystemI','Photosystem I',x)
    rownames(m) <- x
    pnams <- x

    x <- colnames(m)
    x <- gsub('pro.HLunk','pro.HL *',x)
    x <- gsub('pro.HL','HL Pro',x)
    x <- gsub('pro.LL','LL Pro',x)
    x <- gsub('pro\\.','',x)
    colnames(m) <- x

    ## Add some color to the arrows. (Could probably do with mutate().)
    MyCellSpec = function(vec) {
        get_color = function(vec) {
            upArw = c('&#x2B06;','&#x2191;')
            dnArw = c('&#x2B07;','&#x2193;')
            v = ifelse(vec %in% upArw, 'red', ifelse(vec %in% dnArw,'blue', 'black'))
            v
        }
        cell_spec(vec, color = get_color(vec), escape=F)
    }
    m <- apply(m,2, MyCellSpec)
    rownames(m) <- pnams

    ## Additional formatting.  FIXME: Update select() to use a regexp to get HL clades
    ## and then LL clades. Commented out for now b/c HL Pro and LL Pro are no longer
    ## in the table.
    df <- as.data.frame(m) %>%
        ##select(c('HL Pro','LL Pro',everything())) %>%  # selects and orders the columns
        kable(format='html', align='c', caption = caption, escape=F) %>%
        kable_styling('bordered', full_width=F, position='left') %>%
        column_spec(2:ncol(m), width='1cm')
    df
}


## Adapted from NEMO to use HTMLify
MakeNiceEGSEAResultsTable = function(contrast, bigAvgLogfcDir = LOG2FC, pCut=PVALCUT)
{
    m = EGSEAResultsToMatrix(egsea.results, contrast, bigAvgLogfcDir = bigAvgLogfcDir, pCut = pCut)
    if (nrow(m)==0 && ncol(m)==0) { return(NULL) }
    ##m = m[,apply(m,2, function(r) !all(r=='')),drop=F]  # drop strains with nada
    ##if (!all(dim(m) > 0)) { return(NULL) }
    ##m = m[apply(m,1, function(r) !all(r=='')),,drop=F]  # drop pathways with nada
    ##if (!all(dim(m) > 0)) { return(NULL) }
    ## The conventions take more text to explain than I want in the caption.
    df = HTMLifyUpDownMatrix(m, caption=paste("DE pathways in",contrast))
    ##print(df)
    df
}

cat("\n>>> Doing EGSEA <<<\n")
egsea.results = BigFunctionToDoEGSEA(eset, GeneSets, contrastMatrix)
cat(">>> Finished EGSEA <<<\n")

cat("\nSaving main objects to egsea_script.RData\n")
save(egsea.results, Taxa, Pathways, GeneSets, file='egsea_script.RData')

## Example
if (FALSE) {
    contrast <- 'X2D_12h_glucose-X2D_12h_control'
    topSets(egsea.results, contrast=contrast, names.only=F)
    EGSEAResultsToMatrix(egsea.results, contrast, bigAvgLogfcDir=LOG2FC, pCut=PVALCUT)
    MakeNiceEGSEAResultsTable(contrast)
}


for (contrast in colnames(contrastMatrix)) {
    if (FALSE) {
        ## Tables in console
        cat("### EGSEA summary for",contrast,"\n")
        print(t(EGSEAResultsToMatrix(egsea.results, contrast, LOG2FC, PVALCUT)))
        cat("\n")
    } else {
        ## Tables in browser
        df = MakeNiceEGSEAResultsTable(contrast)
        fnam <- paste0('egsea.',contrast,'.png')
        if (!is.null(df)) {
            cat("Saving EGSEA table",fnam,".\n")
            save_kable(df, file=fnam)
        } else {
            cat("No DE gene sets to save so not creating",fnam,"\n")
        }
    }
}

## Also nice to have the full/raw EGSEA results for DE gene sets.
x <- lapply(colnames(contrastMatrix),
            function(cnam) {
                ss <- subset(topSets(egsea.results, contrast=cnam, number=Inf, names.only=F), p.adj < PVALCUT)
                ss$Contrast <- cnam
                ss$GeneSet <- rownames(ss)
                firstTwoCols <- c('Contrast','GeneSet')
                ss <- ss[,c(firstTwoCols, setdiff(colnames(ss), firstTwoCols))]
                rownames(ss) <- NULL
                ss
            })
x <- do.call(rbind,x)
fnam <- paste0('egsea.pCut',PVALCUT,'.csv')
cat("\nSaving table of EGSEA results/statistics for DE pathways to",fnam,".\n",
    "This graphical tables saved earlier are simplified versions of",fnam,".\n")
write.csv(x, file=fnam)

quit(save='no')
