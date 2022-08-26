#!/usr/bin/env Rscript

##
## Process Mari's Prochlorococcus microarrays.  Mainly this was to see if I
## detect the same genes as Mari, after trying out a few approaches to defining
## the noise against which SNR >= 5 to detect genes.  I define per sample noise
## levels that are based on the ERCC (which amplified very well), specifically I
## pick the least concentrated ERCC that has an intensity (predicted from linear
## model) that is > 2x above the intensity of structural hairpins ("3xSLv1"
## probes) in the sample.  This approach detects a total of 619 genes, less than
## the 714 Mari got.
##


## Makes the log look nicer.
BigStep = function(desc)
 {
    cat("\n")
    cat("################################################################################\n")
    cat("##\n")
    cat("## ",desc,"\n")
    cat("##\n")
}



##------------------------------------------------------------------------------
## Libs
##------------------------------------------------------------------------------
BigStep("Loading libraries")
library(Biobase)
library(limma)
library(ggplot2)
library(gridExtra) # for grid.arrange
library(reshape2)
library(MicroTOOLs) # For MyTopTable and various plot funcs
wideScreen(90)



##------------------------------------------------------------------------------
## Params
##------------------------------------------------------------------------------
BigStep("Setting important parameters")
BG_LOWEST_ERCC_2x_ABOVE_HAIRPINS <- TRUE
DETGENECUT <- 5
## Lowering FC from 1.5 to 1.3 because the actual fold changes of the genes
## found when one uses topTreat() will be higher.  See topTreat() docs and
## also the paper on TREAT:  http://europepmc.org/article/PMC/2654802
FC <- 1.3
ADJPVAL <- 0.05
PRED_N_TRANSCRIPTS <- 100

cat("Parameters:\n",
    "\t- Detected genes:  Detected genes in a sample have a signal to noise ratio >=",
    DETGENECUT, ".  Signal is the gene intensity (non-log2) in a sample.  Noise in the",
    "sample is lowest ERCC intensity (from linear model) that is more than twice the",
    "(median) Agilent hairpin intensity.\n",
    "\t- Differentially expressed genes:  Genes with a fold change that is significantly",
    "(p<",ADJPVAL,") less than",FC,"\n")

## Probe annotation mainly used to annotate the *genes* (fData)
##annot60K.csv <- file.path('DropBoxImport','Transcriptomic analysis','UseThisFigS4forFeatureData.csv')
annot60K.csv <- 'UseThisFigS4forFeatureData.csv'
annotTab <- read.csv(annot60K.csv)
## Tweaks.  Probably not needed anymore
x <- annotTab$Pathway
x <- gsub('Fixation','fixation',x)  # 'Carbon Fixation' and 'Carbon fixation'
x <- gsub('Metabolism','metabolism',x)
##x <- gsub('minoacid','mino acid',x)  # Mari fixed
## x <- gsub('Rythm','rhythm',x)       # Mari fixed
x <- gsub('Cycle','cycle',x)
##x <- gsub('Krebs cycle and genes related','Krebs cycle and related',x)
annotTab$Pathway <- factor(x)
x <- c('ProbeID',
       'ID_Gene','Strain','Ecotype','Gene.symbol','Description',
       'Protein_ID2','Protein.ID','Pathway','Location','COG.Function')
stopifnot(x %in% colnames(annotTab))
annotTab <- annotTab[,x]
## Define "SystematicName" which MicroTOOLs functions (e.g. median polishing) require.
## Use the ID_Gene, just as Mari did in Glucose_Hawaii_20Nov.R.  Do *not* use the last
## part of the probe ID which seems to identify the strain, not the gene..
annotTab$SystematicName <- factor(paste0('GID.', annotTab$ID_Gene))
stopifnot( apply(annotTab, 2, function(v) sum(is.na(v))) == 0 ) # Table has "", not NA
cat("Number of SystematicNames imported from Table S4: ",
    length(levels(annotTab$SystematicName)),"\n") # 1326, as stated in manuscript

##
## The following strains have on ecotye annotation but are likely HL because
## they are from 60m.  However, external to this script I checked this as
## follows: All probes for these strains were locally blastn'd against all
## available RefSeq Prochlorococcus nucleotide sequences (n=7666).  Only hits
## that were at >=95%nid and >=95% coverage of the probe were examined (and then
## again with 98% for both cut offs).  The subect RefSeqs that the probes hit
## were tabulated by their ecotype (when known -- used Table S4).  Always there
## were 10's to 100's of hits to HL strains and also to RefSeqs with no ecotype
## -- but almost never (only 2 hits) were there hits to LL strains.  Therefore,
## I will assign these strains to "HL".
##
cat("Assigning lots of HOT208 strains (and MIT9292) to ecotype 'HLunk' based on",
    "external blastn to RefSeq Pro sequences.\n")
likelyHLstrains <- c("HOT208_60m",        "HOT208_60m_805A16", "HOT208_60m_808G21", "HOT208_60m_808M21",
                     "HOT208_60m_810B23", "HOT208_60m_810P02", "HOT208_60m_813102", "HOT208_60m_813B04",
                     "HOT208_60m_813E23", "HOT208_60m_813G15", "HOT208_60m_813L03", "HOT208_60m_813O14",
                     "HOT208_60m_824C06", "HOT208_60m_824E10", "HOT208_60m_826P21", "MIT9292")
x <- annotTab$Ecotype
x[which(annotTab$Strain %in% likelyHLstrains)] <- 'HLunk'
annotTab$Ecotype <- factor(x)
rm(x)

## Define this helper function now because it uses annotTab.  (Used to use
## fData(eset).)  Save matrix of genes (rows) by cols (samples or comparison),
## along with annotation.  Matrix row names must be the gene SystematicNames.
WriteMatrixCsv <- function(mat, filename)
{
    want <- c('Ecotype','Strain', 'Gene.symbol', 'ID_Gene', 'Pathway',
              'Protein.ID', 'Protein_ID2','COG.Function',
              'Description', 'Description2', 'Location')
    want <- intersect(want, colnames(annotTab))  # because no Description2 if did my own probe->gene
    idx <- match(rownames(mat), annotTab$SystematicName)  # match 1st converts factors to charvecs
    stopifnot(!is.na(idx))
    featData <- annotTab[idx,want]  # 'featData' b/c used to use fData(eset)
    df <- cbind(mat,featData)
    cat("\nSaving table of differentiallly expressed genes to",filename,"\n")
    write.csv(df[order(df$Ecotype,df$Strain,df$Pathway,df$Gene.symbol),], filename)
}


##------------------------------------------------------------------------------
## IF DOING MY OWN QUANTILE NORMALIZATION AND MEDIAN POLISHING:
##
## Load probes, bg correct, normalize. Argh, but there's ERCC too.  This is not
## a small effort.  Will wait to see Mari's reply to email sent Jan 31.  (Note
## that Irina and I in April 2019 discussed the Pro arrays with Mari, some
## background issues.  Sounds like she redesigned/ran the arrays so that probes
## for strains absent in her samples would not be included/used (I think
## assuming these strains were present lead to BG being set too low, but not
## sure and would have to study the thread more carefully).
##   -- Have since switched to BG defined based on lowest ERCC that is
##      2x above the hairpins. So strains (absent or present) cannot impact
##      on BG.
## ------------------------------------------------------------------------------

BigStep("Loading the the arrays")
cat("Will load each array and do some checks on the feature extraction software protocol and",
    "that background subtraction was not already done.\n")

## From MicroTOOLs, I think these are the only features I need.
requiredFeatures = c('ProbeName', 'gProcessedSignal', 'gProcessedSigError',
    'FeatureNum', 'Row', 'Col', 'SubTypeMask', 'ControlType',
    'SystematicName', 'PositionX', 'PositionY')
probeDataList <- list()
probeFEParamsList <- list() # To check how feature extraction software was run. (added 15Mar2021)
rawArrayDir <- file.path('/Volumes','JMagasin_SD1','Experiments.noindex',
                         'MMunozProchloroArrays','Samples','60Karrays')
expProto <- 'GE1_1200_Jun14 (Read Only)' # Most of the arrays had this protocol used, but a few of
                                         # them had GE1_1105_Oct12.  This should not matter w.r.t.
                                         # background subtraction, which is off for most protocols
                                         # based on Table 13 in the ref guide (link below).
for (af in list.files(rawArrayDir, '^GSM.+\\.txt\\.gz$', full.names=T)) {
    cat('Reading',af,'\n')
    lines <- readLines(af, n=50, encoding='UTF-8')
    ## First 3 lines have feature extraction software parameters.
    df <- data.frame(feparam = strsplit(lines[[2]],"\t")[[1]],
                     value   = strsplit(lines[[3]],"\t")[[1]],
                     type    = strsplit(lines[[1]],"\t")[[1]])
    probeFEParamsList[[af]] <- df
    ## Check which feature extraction sofware protocol was used, and that
    ## background correction was not done.  Probably it would not matter if a
    ## different protocol was used.  If BG sub were done it also probably would
    ## not matter. We do our own later, and if it were already done then there
    ## likely would just be less to subtract.
    ## The docs for Protocol is GE1_1200_Jun14 are here (see page 5):
    ##    https://www.agilent.com/cs/library/usermanuals/public/G4460-90052_FE_RefGuide.pdf
    if ((x <- as.character(subset(df, feparam=='Protocol_Name')$value)) != expProto) {
        warning('Feature extraction software protocol for ',af,' was ',x,' rather ',
                'than ',expProto,' as expected.\n')
    }
    if ((x <- subset(df, feparam=='BGSubtractor_BackgroundCorrectionOn')$value) != 0) {
        warning('The feature extraction software might already have done background',
                'subtraction for ',af,' because BGSubtractor_BackgroundCorrectionOn',
                'was ',x,' rather than 0.\n')
    }
    if ((x <- subset(df, feparam=='BGSubtractor_BGSubMethod')$value) != 7) {
        warning('For ',af,' expected BGSubMethod to be 7 ("none") but got ',x,'\n')
    }
    if ((x <- subset(df, feparam=='BGSubtractor_SpatialDetrendOn')$value) != 1) {
        warning('For ',af,' expected SpatialDetrendOn to be 1 but it was ',x,'\n')
    }

    ## The data must begin after the FEATURES line.
    stopifnot( length(featLine <- grep('^FEATURES', lines)) == 1 )
    stopifnot( grepl('^DATA', lines[featLine+1]) )
    ## I know that read.table() messed up with missing files (\t\t) but am
    ## hoping I won't see that in Mari's arrays.  If it happens, see how to
    ## handle in MicroTOOLs' ReadFilesIntoListOfDataFrames().
    probeDataList[[af]] <- read.table(af, sep='\t', skip=featLine-1, header=T)[,requiredFeatures]
}
rm(requiredFeatures, rawArrayDir, af, featLine, lines, expProto, x, df)
stopifnot(lapply(probeDataList, nrow) == 62976)  # better all have same probe count!
## Simplify names (so they are not file paths)
names(probeDataList) <- gsub('^.+60Karrays/(GSM.+)\\.txt\\.gz', '\\1', names(probeDataList))
names(probeFEParamsList) <- names(probeDataList)
rm(probeFEParamsList)  # Not used below. Created for interactive mode.


BigStep("Making probe ExpressionSet")

## From MicroTOOLs (not exported).  There is an analogous function for the
## errors but not sure I will need that.
MakeMatrixOfProbeValues <- function(dataList, wantedValue)
### Make a table of probe values: data.probes
###   rows = probes; named 1..numProbes
###   cols = samples (or "hybridizations"); named as in dataList items
### Probe intensities are from gProcessedSignal because that is corrected for
### background.  [NO!  Depends on how FE software was run.] Other columns of
### interest are "gBGMedianSignal","gBGSubSignal".
### Note: Because dataList is not modified by this function, it will not be
### copied, so there should be no memory thrashing. ("pass by promise")
{
    ## All microarrays must have probes like the first so get first's names.
    probeNames <- as.character(unlist(dataList[[1]]$ProbeName))
    numProbes <- length(probeNames)
    numSamples <- length(dataList)

    ## Note: Probe names are not unique so do not use as matrix rownames
    probeMatrix <- matrix(nrow=numProbes, ncol=numSamples)
    for (i in 1: numSamples) { probeMatrix[,i] <- dataList[[i]][,wantedValue] }
    rownames(probeMatrix) <- c(1:numProbes)    # numeric = row names (simple)
    colnames(probeMatrix) <- names(dataList)   # arrays = col names
    return(probeMatrix)
}
data.probes <- MakeMatrixOfProbeValues(probeDataList, 'gProcessedSignal')

featData <- probeDataList[[1]][,c('ProbeName','SubTypeMask','ControlType','SystematicName')]
featData$originalFeatureNumber <- 1:nrow(featData)  # Because merge() messes up row order
## IMPORTANT: At this point the SystematicName is equal to (almost always) the
## ProbeName (62K total), so it does not identify the gene (probe set).  We need
## it to, so that median polishing will work.  (Otherwise it will only median
## polish spots for a probe.)
## Get the SystematicNames from annotTab, looking them up via the ProbeName.
featData <- merge(x=featData, y=annotTab, by.x='ProbeName', by.y='ProbeID',
                  all.x=T, all.y=F, sort=FALSE)
## Control probes (e.g. (+)E1A_r60_1 and GE_BrightCorner) have SystematicName.x
## from original featData but NA for SystematicName.y.  Use the one from
## original featData.
sny <- as.character(featData$SystematicName.y)
idx <- which(is.na(sny))
cat("There are",length(idx),"probes that do not have a SystematicName.",
    "These are mainly control probes (e.g. (+)E1A_r60_1).\n")
sny[idx] <- as.character(featData$SystematicName.x)[idx]
stopifnot(!is.na(sny))
## Also standardize ERCC SystematicName's.  The ProbeName always looks good but
## the SystematicName.y's are weird, e.g. DQ459430 for ERCC-00002_129. So use
## the ERCC name, and also drop the _129 because erccMeta doesn't have that.
## FIXME:  Maybe these DQ459430 etc. are not ERCC after all and that is why
## my ERCC plots make no sense.
idx <- grep('^ERCC',featData$ProbeName)
sny[idx] <- gsub('^(ERCC-[0-9]+)_.*','\\1',featData[idx,'ProbeName'])
## Finally, replace the SystematicName.{x,y}
featData$SystematicName <- factor(sny)
## Now fix up the rows so they are in the original order. This is critical.
## Had a bug before realized that merge() changed the order --> probe values
## and features were not in correspondence.
featData <- featData[order(featData$originalFeatureNumber, decreasing=F),]
rownames(featData) <- NULL
featData <- featData[,!colnames(featData) %in% c(paste0('SystematicName.', c('x','y')), 'originalFeatureNumber')]
## (Various manual checks that all looks fine.)

## Finally, the eset.  Lacks pData.
eset.probes <- new("ExpressionSet", exprs = data.probes,
                   featureData = new("AnnotatedDataFrame", data=featData))
rm(featData, data.probes, sny, idx)
cat("Created ExpressionSet for the probes.  There are", nrow(eset.probes),"probes,",  # 62976
    length(levels(fData(eset.probes)$SystematicName)),"SystematicNames,",             # 1340
    length(grep('^GID',levels(fData(eset.probes)$SystematicName))),"SystematicNames for UCYN-A genes,",#1200
    length(grep('^ERCC',levels(fData(eset.probes)$SystematicName))),"ERCC.\n")        #   96

## Plot raw log2 distributions, and do array quality checks.
cat("Plotting raw probe intensities in probes_boxplot.png\n")
PlotSampleDistributions(eset.probes,"probes_boxplot.png")
cat("Running Array Quality Metrics. This takes a minute...\n")
arrayQualityMetrics(eset.probes, outdir='AQMforMari60arrays', do.logtransform=T)
##  Mostly looks good.  Only GSM4675253_2D_24h_control (array #15) is flagged
##  because it is more distant from other arrays.  It has a weaker overall
##  intensity.  Could remove it... Yes, let's remove it.
cat("Removing GSM4675253_2D_24h_control because it was flagged by AQM.\n")
eset.probes <- eset.probes[, colnames(eset.probes) != 'GSM4675253_2D_24h_control']


BigStep("Background correction")
## Background correction.
## 15Mar2021: Added code above that shows that FE software had global background
## correction disabled (but spatial detrending was on, to handle different
## brightness depending on array position).  Do background correction here.
##
## The R docs for backgroundCorrect() are ~confusing, but they point to
## normexp.fit() which is more clear on how it works. In particular, I wanted to
## check that passing just the gProcessedSignal is sufficient (and that
## gBGMedianSignal or other per feature measurements is not necessary).  I
## suppose I *could* pass in Eb = gBGMedianSignal and method=normexp, which
## would cause first the subtraction of gBGMedianSignal and then normexp.fit
## (see the code).  Not sure it's worth it.  Would expect that normexp.fit()
## will determine the error (the normal part) regardless of whether one first
## substracts gBGMedianSignal (and the normal would just have a smaller mean if
## one did pass in Eb).  So I'm not going to pass in Eb.  Changed from using
## backgroundCorrect() to backgroundCorrect.matrix() only because "RG" sounds
## like two color arrays, but the former ends up calling the latter anyway.
exprs(eset.probes) = backgroundCorrect.matrix(eset.probes, method = 'normexp')


BigStep("Quantile normalization")
## Normalize!  Similar to MTOOLs, remove all non-experimental probes before
## normalizing, e.g. so that bright corners don't influence.  Drops
## Bright/DarkCorner and hairpins.  Only keep experimental (0), including ERCC.
## 15Mar2021: Keep the hairpins for the moment too.
cat("There are", nrow(eset.probes),"probes.",  # 62976
    "Will only keep experimetal probes (SubTypeMask==0) and hairpins.\n")
idx <- which((fData(eset.probes)$SubTypeMask == 0) |
             (fData(eset.probes)$SystematicName == "3xSLv1"))
eset.probes <- eset.probes[idx,]
cat("Now there are", nrow(eset.probes),"probes\n")  # 62577. Wouuld be 62269 w/o hairpins which is
                                                    # exactly what Mari has in Glucose_Hawaii_20Nov.R
eset.probes.norm <- normalize.ExpressionSet.quantiles(eset.probes, transfn='none')
## The eset has *NON*-log2 values.  PlotSampleDistributions() will take the log2.
rm(idx)
cat('Normalization done! Boxplot in probes_boxplot_norm.png\n')
PlotSampleDistributions(eset.probes.norm, 'probes_boxplot_norm.png',
                        main='Normalized probes')  # Wahoo! How orderly!
## Interquartile ranges are 0-5.


##------------------------------------------------------------------------------
## Get the structural hairpin negative controls 3xSLv1 (SubTypeMask 66)
##------------------------------------------------------------------------------
BigStep("Getting structural hairpins")
cat("The structural hairpins probes ('3xSLv1') are effectively negative controls",
    "because they should not hybridize to anything. We will use them later to help",
    "define the noise in each sample. (Even though we already subtracted the background.)\n")

idx <- which(as.character(fData(eset.probes.norm)$SystematicName) == "3xSLv1")
hairpins <- log2(exprs(eset.probes.norm)[idx,])
## Now drop them.
eset.probes.norm <- eset.probes.norm[-idx,]
cat("Dropped hairpins from the ExpressionSet which now has",
    nrow(eset.probes.norm), "probes.\n")  # Now it is 62269
rm(idx)

## Look over hairpins. Like ERCC, the hairpin intensities are from after BG
## correction and normalization.
cat("Here are the hairpin median intensities (log2) for each sample, if you want to",
    "compare them to the normalized intensities in probes_boxplot_norm.png, which have",
    "medians that range from",
    paste(round(range(log2(apply(exprs(eset.probes.norm),2,median))),1),collapse=' to '),".\n")
round(log2(apply(2^hairpins,2,median)),1)  # Range -2.5 to -0.83
cat("\nHere are the hairpin maxima (log2). Some are high, so later we",
    "use the *medians* for determining the noise level in each sample.\n")
round(log2(apply(2^hairpins,2,max)),1)



##------------------------------------------------------------------------------
## ERCC
##------------------------------------------------------------------------------

BigStep("Getting ERCC metadata")
cat("We need ERCC copy numbers for making linear models later.\n")

erccMeta <- read.csv('ercc.meta.fromMari.csv', row.names=1)
erccMeta <- erccMeta[,c('Subgroup','Attomoles_ul','Copy_number','Copy_number_used')]
## Looks like Mari diluted ERCC by 200x.
stopifnot(abs(erccMeta$Copy_number / erccMeta$Copy_number_used - 200) < 0.0001)
## Sort in increasing order of concentration.
erccMeta <- erccMeta[order(erccMeta$Copy_number_used),]

## I already removed the weird trailing parts of ERCC names when I made the eset,
## so the eset's ERCC names should match what is in erccMeta.  (I did not fix the
## ProbeNames but that should be okay.)
x <- levels(fData(eset.probes.norm)$SystematicName)
idx <- grep('^ERCC',x)
stopifnot(rownames(erccMeta) %in% x[idx])
if (!all(x[idx] %in% rownames(erccMeta))) {
    ## Think I saw this in NEMO where the ERCC on the array were not
    ## quite in sync with the ERCC table. (Probes dropped from the kit
    ## b/c they were not reliable?)
    warning('These ERCC are in the array but not the ERCC metadata.  According to ',
            'https://www.jgenomics.com/v04p0019.htm they are excluded from the commercial ',
            'ERCC kit.  Good, then we can use them as negative controls:)  ',
            paste(setdiff(x[idx],rownames(erccMeta)),collapse=', '))
    ## These are the ones warned about: ERCC-00007, ERCC-00018, ERCC-00023, ERCC-00128
}

idx <- which(fData(eset.probes.norm)$SystematicName %in% rownames(erccMeta))
cat("Saving boxplot of ERCC levels to probes_boxplot_norm_ERCC.png.\n")
PlotSampleDistributions(eset.probes.norm[idx,], "probes_boxplot_norm_ERCC.png",
                        main='ERCC probes (normalized)')  # Again, this func takes log2 of the data
## Compare to plot of all probes.  The IQR's for the ERCC are ~4-10 with some
## variability, vs. for all probes seem 0-5 and extremely even.  Points to
## that gene detection is going to throw out a bunch of probes.
rm(x,idx)


##------------------------------------------------------------------------------
## Median polish (and make gene eset)
##------------------------------------------------------------------------------
BigStep("Median polishing to convert probe intensities to gene intensities")
## Note that the eset produced will have log2 levels!

## Requires 'SystematicName' be the gene name in the probe eset.
geneIntensitiesAndErrors <- ProbeIntensitiesToGeneIntensities(eset.probes.norm, 'medianpolish')
exprs <- geneIntensitiesAndErrors$exprs
idx <- match(rownames(exprs), annotTab$SystematicName)  # exprs row names are SystematicName
## These are the probes with no annotation.  Check manually that they are just ERCC and hairpins.
unmatched <- rownames(exprs[which(is.na(idx)),])
cat("The following targets had no annotation:\n")
print(unmatched)
if (FALSE) {
    ## DON'T NEED THIS any more because below I create a log2 ERCC matrix
    ## that takes the median of the ERCC in each sample.

    ## Grab the ERCC and order them by concentration.  These are
    ## post-normalization so they are log2.
    ERCC <- exprs[grep('^ERCC',rownames(exprs)),]
    ERCC <- ERCC[rownames(erccMeta)[rownames(erccMeta) %in% rownames(ERCC)], ]
}
## Drop unmatched from the exprs
exprs <- exprs[!is.na(idx),]
## Prepare feature data for the gene eset.  Recycle the probe fData.
idx <- match(rownames(exprs), fData(eset.probes.norm)$SystematicName)
stopifnot(!is.na(idx))
featData <- fData(eset.probes.norm)[idx,]
stopifnot(nrow(featData)==nrow(exprs))
rownames(featData) <- featData$SystematicName
rownames(exprs) <- featData$SystematicName
idx <- which(colnames(featData) %in% c('ProbeName','SubTypeMask','ControlType','ID_Gene','SystematicName'))
featData <- featData[,-idx]
eset <- new("ExpressionSet", exprs = exprs,
            featureData = new("AnnotatedDataFrame", data=featData))  # no pData
cat("The gene ExpressionSet has",nrow(eset),"genes and",ncol(eset),"samples.\n")  # 1200 x 15
rm(unmatched,idx,exprs,featData)

## Just to emphasize that the eset has log2 levels (vs. non log2 which will
## always be >0).  Iirc this is not documented but is simply what the underlying
## median polishing func does.
stopifnot(apply(exprs(eset),2,min) < 0)  # Does not *have* to be true but it is for Mari's data.


##------------------------------------------------------------------------------
## Gene detection: For each gene its SNR must be >= 5, with SNR calculated as
## SNR = (Si – BG)/BG.  Si is gene's intensity and BG is the "low ERCC signals
## non-amplified."  <--- NO
## 7 Feb 2021:  Change BG to be specific to sample and based on the lowest 10%
## of gene intensities.
## 15 Mar 2021:  Determine BG based on linear models for the ERCC and negative
##    controls (hairpins and ERCC that were not in the kits).
## 30 July 20201: The actual check below does not subtract BG in the numerator.
##    Search below for:  rownames(dat)[dat[,i] >= DETGENECUT*2^bg[i]]
## ------------------------------------------------------------------------------
BigStep("Determination of noise level to use for gene detection")

## Make ERCC matrix with the log2 of the median of each ERCC in a sample.
## Matrix later used for gene detection.
MakeLog2ERCC <- function()
{
    ## Get ERCC levels.  Look up in erccMeta rownames.  (Use SystematicName, not
    ## ProbeName because ProbeNames for ERCC have the weird suffixes, which are
    ## not in the rownames(erccMeta).)
    x <- levels(fData(eset.probes.norm)$SystematicName)
    x <- x[grep("^ERCC", x)]
    x <- x[x %in% rownames(erccMeta)]
    idx <- which(fData(eset.probes.norm)$SystematicName %in% x)    
    dat <- data.frame(exprs(eset.probes.norm)[idx,])
    dat$ERCC <- factor(as.character(fData(eset.probes.norm)$SystematicName[idx]))
    stopifnot(dat$ERCC %in% rownames(erccMeta))

    ## Within each sample, for each ERCC, take the median of the non-log2
    ## levels.
    dat <- aggregate(dat[,setdiff(colnames(dat),'ERCC')], by=list(ercc=dat$ERCC), median)
    rownames(dat) <- dat$ercc
    dat <- dat[,-1]  # drop 'ercc' column
    stopifnot(rownames(dat) %in% rownames(erccMeta))
    dat <- dat[rownames(erccMeta),]   # Get in concentration order
    log2(dat)
}
cat("Getting ERCC log2 intensities.\n")
ERCC <- MakeLog2ERCC()


PlotERCC <- function(eset, si)
{
    idx <- grep('^ERCC',fData(eset)$SystematicName)
    df <- data.frame(
        sysname   = fData(eset)$SystematicName[idx],
        intensity = exprs(eset)[idx,si])
    idx <- match(as.character(df$sysname), rownames(erccMeta))
    df <- df[-which(is.na(idx)),]
    idx <- match(as.character(df$sysname), rownames(erccMeta))
    stopifnot(!is.na(idx))
    df$copy = erccMeta[rownames(erccMeta)[idx],'Copy_number_used']
    df
}


## Quick plot. Note that ERCC is ordered by concentration so x axis will also be.
CheckERCC <- function(i)
{
    dat <- data.frame(x=log2(erccMeta[rownames(ERCC),'Copy_number_used']),
                      y=ERCC[,i])
    m <- lm(y~x,dat)
    plot(x = dat$x,
         y = dat$y,
         main=paste('ERCC for',colnames(ERCC)[i]),
         xlab='log2 of copy number used', ylab='log2 gene intensity')
    points(x = dat$x,
           y = predict(m),
           type='l', col='red')
    legend(x='topleft', legend=paste('adj r^2 =',round(summary(m)$adj.r.squared,3)),
           bty='n')
    cat(colnames(ERCC)[i],'linear model has adjusted r^2 = ', round(summary(m)$adj.r.squared,3),'\n')
}
## Checked them all.  Only the lowest ERCC is clearly noise.  The next ERCC is
## near 5 and all are ~linear with adj r2 of 0.85-0.89 for all.
cat("Saving panel plot of ERCC linear models to ERCC_models.png\n")
png("ERCC_models.png", height=12, width=12, units="in", res=144)
par(mfrow=c(4,4))
for (i in 1:ncol(eset)) { CheckERCC(i) }
dev.off()


## From inspecting each of the models (in individual plots, not the big panel
## just drawn), the third least concentrated ERCC (61,98,117 each with ~172
## transcripts added) mainly agree with the model prediction. And the predicted
## intensity is about 2.5.


## Model ERCC transcript copies (log2) as a function of itensities observed for ERCC (log2).
MakeModels <- function(response)  # response is 'copy' or 'intensity'
{
    mods <- apply(ERCC, 2, function(ev) {
        dat <- data.frame(copy      = log2(erccMeta[names(ev),'Copy_number_used']),
                          intensity = ev)
        if (response=='copy') {
            m <- lm(copy~intensity, dat)
        } else {
            m <- lm(intensity~copy, dat)
        }
        stopifnot(summary(m)$adj.r.squared > 0.79)
        m
    })
    mods
}
transcript.models <- MakeModels('copy')  # models predict *log2* transcripts

## Use each model to predict the number of transcripts for
## the hairpins' observed intensities.
x <- sapply(colnames(hairpins),
            function(samp) {
                mod   <- transcript.models[[samp]]
                log2inten <- hairpins[,samp]
                inten <- data.frame(intensity=log2inten)
                predict(mod, newdata=inten)
            })
## Each sample (columns) is comprised of all the predicted log2 copy numbers
## for the hairpin log2 intensities.
cat("Remember the hairpins (negative controls)? Let's see how many transcripts",
    "would be predicted for them using the ERCC linear models.  Keep in mind that",
    "the least concentrated ERCC has Copy_number_used =",round(erccMeta[1,'Copy_number_used'],1),"\n")
print(round(2^apply(x,2,median),1)) # Range of 13-68, which agrees quite well with
                                    # the least concentrated ERCC with 43 copies added.


## What intensity predicted to observe N transcripts?
intensity.models <- MakeModels('intensity')
bg.Ntranscripts <- sapply(intensity.models,
            function(mod) {
                predict(mod, newdata=data.frame(copy=log2(PRED_N_TRANSCRIPTS)))
            })
cat("\n\nWhat intensity (log2) do the models predict for",PRED_N_TRANSCRIPTS,"transcripts?\n")
print(round(bg.Ntranscripts,1))
cat('\nSummary of predicted intensities for',PRED_N_TRANSCRIPTS,'transcripts:\n')
print(summary(bg.Ntranscripts))  # log2 intensities: mean=2.9, range=[1.8,4.0]
            # These are much higher (for N=500) than the predicted log2 medians for hairpin intensities
            # which were all negative.


## Take the first ERCC with predicted intensity that is > twice the median of
## the hairpins in the sample.  Rationale is that hairpins are noise but ERCC
## were definitely present and the linear models look good.  (Plus we're going
## to take SNR > 5 so should be well above noise.)
preds <- sapply(intensity.models, function(m) predict(m)) 
x <- log2(2 * apply(2^hairpins,2,median))  # Log 2 of *twice* the hairpins' median intensities
stopifnot(rownames(preds) == rownames(erccMeta))  # Recall, sorted from LOW to HIGH concentration
bg.hairpin2x <- sapply(names(x), function(samp) {
    names(which(preds[,samp] > x[samp]))[1]  # Gets name of 1st ERCC that is >2x hairpin (b/c preds sorted by conc)
})
## Set aside for a nice table to print in a moment
df <- data.frame(ERCC=bg.hairpin2x)
df$Copy_number_used <- round(erccMeta[as.character(df$ERCC),'Copy_number_used'],2)
bg.hairpin2x <- sapply(names(bg.hairpin2x), function(samp) {
    preds[bg.hairpin2x[samp],samp]  # Get predicted value of the ERCC that is >2x hairpin
})
stopifnot(rownames(df) == names(bg.hairpin2x))
df$noiseLevel <- round(bg.hairpin2x,2)
cat("\n\nWill define the following noise levels (log2 intensities) for each sample.",
    "Method = for each sample use the first ERCC that has a predited intensity that is more than twice",
    "the median of the hairpins.\n")
print(df)
summary(bg.hairpin2x)  # mean=1.7, range=[1.4,2.0]
rm(x,preds,df)


## Had considered using the bottom 10% (etc.) of genes for the BG, but I don't
## think that's a great approach (mabye okay if you had a genome for one
## organism).
bg.lowPct <- log2(t(apply(2^exprs(eset), 2,
                          function(ev) quantile(ev,probs=c(0.1,0.25,0.35,0.40)))))

##
## Here are the options for BG.  log2 intensties
##
cat("\n\nThis table shows log2 intensities for N=",PRED_N_TRANSCRIPTS,"transcripts, 2x the hairpins (which is the",
    "noise level we will use), and various bottom % of genes.\n")
df = data.frame(transcripts.N  = bg.Ntranscripts,
                hairpin2x      = bg.hairpin2x,
                lowPct         = bg.lowPct)
print(round(df,2))
cat("\nNotice that hairpin2x is higher than then lowest 10%, and sometimes 25%, of genes.")



## Secondary negative control: For comparison, what values do we see for the
## ERCC that were on the array but apparently not in the kit (or at least not in
## erccMeta)?  These are probe intensities but should be a fair comparison to
## bg.ercc because still have BG correction and normalization.
cat("Yet another comparison:  These ERCC are on the array but are not in the commercial kit,",
    "so they act as negative controls. These are the log2 of the medians of the probes in each sample.",
    "Also showing the medians of the hairpins for comparison (not the 2x noise level, the hairpins!)\n")
df <- as.data.frame(sapply(c('ERCC-00007', 'ERCC-00018', 'ERCC-00023', 'ERCC-00128'),
                           function(x) {
                               log2(apply(exprs(eset.probes.norm)[grep(x,fData(eset.probes.norm)$SystematicName),],
                                          2, median))
                           }))
x <- log2(apply(2^hairpins,2,median))
stopifnot(names(x) == rownames(df))
df$hairPins <- x
print(round(df,2))
cat("\nI would expect these missing ERCC to have log2 intensities similar to the hairpins",
    "and that seems to be the case.  Good:)  Probably we could use either the hairpins or",
    "these absent ERCC to define the noise level.\n")
rm(df,x)

## Let's go with hairpin2x because (1) it's above the hairpins; (2) its above
## the absent ERCC just checked; (3) transcripts500 is probably too high since
## it often throws out >50% of genes.
if (BG_LOWEST_ERCC_2x_ABOVE_HAIRPINS) {
    cat("Background defined based on 2x the observed intensities for structural",
        "hairpins.\n")
    bg <- bg.hairpin2x
} else {
    stop("You should define noise level based on the hairpins!\n")
}



BigStep("Detecting genes")

## This is where we look in each sample for genes with signal >= 5xBG.
dat <- 2^exprs(eset)
detGenes <- lapply(1:ncol(dat),
                   function(i) {
                       rownames(dat)[dat[,i] >= DETGENECUT*2^bg[i]]
                   })
cat("Here are the numbers of detected genes in each sample.\n")
print(summary(sapply(detGenes,length)))  # min=416, max=538, mean=488 genes detected
detGenes <- unique(unlist(detGenes))
cat("Total num genes detected (in one or more samples):",
    length(detGenes),"\n") # 619 genes detected (in >=1 sample) if BG_LOWEST_ERCC_2x_ABOVE_HAIRPINS,
                           # versus 714 Mari got.
eset.beforeDetection <- eset # Only for density plots below
eset.notDet <- eset[ !rownames(eset) %in% detGenes,]
eset <- eset[detGenes,]
cat("Removed",nrow(eset.notDet),"undetected genes from the ExpressionSet. Now there are", nrow(eset),"genes.\n")


## For each strain, compare the number of genes that were detected (in >=1 sample) to
## the number that were not detected (but that are on the array).
df <- merge(x=data.frame(table(fData(eset)$Strain)),
            y=data.frame(table(fData(eset.notDet)$Strain)),
            by='Var1')
colnames(df) <- c('Strain','DetGenes','UndetGenes')
df$Total <- df$DetGenes + df$UndetGenes  # num genes on array for strain
df$PctDet <- round(100*df$DetGenes/df$Total,1)
df$Ecotype <- fData(eset)$Ecotype[match(df$Strain, fData(eset)$Strain)]
cat("\n\nThis table describes genes detected by strain.\n")
df <- df[order(df$Ecotype,df$DetGenes,df$PctDet, decreasing=T),]
print(df)
cat("Total genes on array (detected or not): ",sum(df$Total),"\n")
rm(df)
## Most HL strains have a high % (much more than 50%) of their genes detected
## (not MED4, MIT9515). LL have <10%.  (Strains have mainly 30-40 genes on the
## array, so the %'s not w.r.t. a tiny number.)


##
## Hack!  Change sample names so that same as if I loaded the normalized
## genes Mari provided, so that can use DE code without modification.
##
colnames(eset) <- gsub('^GSM[0-9]+_','X',colnames(eset))
colnames(eset.beforeDetection) <- gsub('^GSM[0-9]+_','X',colnames(eset.beforeDetection))

## Save table of detected genes.
WriteMatrixCsv(exprs(eset), paste0('detectedGenes.csv'))


##
## NMDS to see if replicates cluster, and if samples cluster by time of day.
## (Experiment 1 samples were collected at noon, with glucose addition at that
## time.  Experiment 2 samples were collected at 4pm, with glucose addition at
## that time.  See paragraph in Mari's transcriptomic results.)
##
library(ggrepel)
dat <- exprs(eset)  # all detected genes
dat <- dist(1-cor(dat))  # distances (Euclidean) based on correlation between
                         # sample expression profiles.
nmds <- isoMDS(dat, k=2)
stress <- round(nmds$stress/100, 2)
cat("NMDS of Euclidean distances between 1-correlation of samples had stress", stress, "\n")
dat = data.frame(nmds$points, rownames(nmds$points))
colnames(dat) = c('x','y','Sample')
g <- ggplot(dat, aes(x,y,color=Sample)) + geom_point(size=3) +
    ## coord_fixed(ratio=1) +  I've never seen such a 1D NMDS. Haha, no wonder the stress is so low:)
    labs(title = paste("Samples represented by transcript levels of",nrow(eset),"detected genes"),
         caption = paste0("NMDS on Euclidean distances between 1-correlation; stress = ", stress,"."),
         x=NULL, y=NULL) +
    geom_text_repel(aes(label=Sample), size=3) +
    theme(legend.position="none") + theme_bw()
ggsave('detetectedGenesNMDS.png') # 7"x7"


##
## Rest of this section is to make a panel plot of all the gene intensity
## distributions.
##

##
## Density plot of gene intensities. Usually you will pass an eset that has just
## one sample (and use exprs(eset)[,samp,drop=F]). You can also pass in noise
## levels for the sample(s) a vector and the detection levels.  Noise vertical
## lines will all be orange. Detection vertical lines will all be green.  (So
## you do not need to name the entries in these vectors with the sample.)
##
PlotGeneIntensities <- function(eset, noise=NA, detect=NA)
{
    gdat <- melt(exprs(eset))
    colnames(gdat) <- c('Gene','Sample','Intensity')
    g <- ggplot(gdat, aes(x=Intensity)) + geom_density() +
        theme_bw() +
        labs(title=paste(gsub('^X','',levels(gdat$Sample)), collapse=', '),
             x='Gene intensity (log2)')
    if (!is.na(noise)) {
        g <- g + geom_vline(xintercept=noise, color='orange')
    }
    if (!is.na(detect)) {
        g <- g + geom_vline(xintercept=detect, color='green')
    }
    g
}

glist <- lapply(colnames(eset.beforeDetection), function(samp) {
    ## Eset col names start with X and lack the GSM....  'bg' names have the GSM...
    samp2 <- match(gsub('^X','',samp), gsub('^[^_]+_','',names(bg)))
    stopifnot(!is.na(samp2))
    bg.s2 <- bg[samp2]
    PlotGeneIntensities(eset[,samp,drop=F], noise=bg.s2, detect=log2(DETGENECUT) + bg.s2)
})

g <- arrangeGrob(grobs=glist)
##grid.arrange(grobs=glist)  # to draw now
ggsave('gene_intensities.png', plot=g, width=14, height=12, units='in', dpi=72)
cat('\nSaved gene_intensities.png which shows the log2 intensities of all genes, before',
    'determingin which were detected, for each sample. Also shows the noise level (orange)',
    'and the detection threshold (green).\n')
## Fixme: Would be cool to show in each plot the number of genes that were detected,
## perhaps by adding a text parameter to PlotGeneIntensities().
rm(glist,eset.beforeDetection)



## ------------------------------------------------------------------------------
## ECC fold changes
## ------------------------------------------------------------------------------
BigStep('Extra analysis of ERCC fold changes')
cat('Plot the distribution of fold changes of ERCC intensities _within_ each sample,',
    'to see how much more intense the most concentrated ERCC are.  See script for',
    'more info. Ideally the distributions would be similar across samples.\n')

## Each ERCC theoretically has the same number of transcripts added to each
## sample.  If each sample had the same number of Prochlorococcus community
## transcripts, then the relative abundance of the i'th ERCC would be the same
## across samples.  Although neither of these assumptions is true, perhaps the
## most abundant ERCC will be less affected (i.e. it will still dominate the
## pool of transcripts hybridized to the array).  So we might expect them to
## have similar intensities across arrays.  Let's check.

## For each of the normalized ERCC intensities in 'ERCC', divide by the min to
## get the fold change for each ERCC across samples.  Theoretically it will be
## closer to 1 for more concentrated ERCC.
##  -- Intensities in 'ERCC' are log2 of the median values seen across all spots
##     in each sample.
dat <- data.frame(as.matrix(t(apply(2^ERCC, 1, function(rv) rv/min(rv)))))
dat$ERCC <- rownames(dat)
if (T) {
    ## Optionally limit to ERCC with high (or low) concentration.
    x <- rownames(subset(erccMeta, Copy_number_used > 1E3))
    dat <- subset(dat, ERCC %in% x)
}
## rowMeans(dat[,-16])   # On average what are the fold changes for each ERCC.
##                       # Even if limit to transcripts >1E4, often the fc is
##                       # >2 (mean=2.5, max=4.1)
dat <- melt(dat, id='ERCC')
colnames(dat) <- c('ERCC','Sample','foldchange')
## Sort the ERCC by intensity!
dat$ERCC <- factor(dat$ERCC, levels=rownames(erccMeta))
if (F) {
    g <- ggplot(dat, aes(x=Sample, y=ERCC, fill=foldchange)) +
        geom_raster() +
        labs(title='ERCC fold changes across samples', x='', y='') +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
} else {
    ## Boxplot maybe better
    g <- ggplot(dat, aes(x=ERCC, y=foldchange)) +
        geom_boxplot() +
        geom_hline(yintercept=1, color='red') + # ideally, fc=1
        labs(title='ERCC fold changes across samples', x='', y='Fold change',
             caption=paste('Fold change is relative to the lowest intensity ERCC in',
                           'each sample. Only using ERCC with >1K copies.')) +
        theme_bw() +  # do this before rotating x text
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}
plot(g)
cat('Plotting ERCC_fold_changes.png.\n')
ggsave('ERCC_fold_changes.png',width=10, height=6, units='in', dpi=144)

## All ERCC: Usually the fold changes (of ERCC i relative to the min for i
##   across samples) are < 5.  Seem to go as high as >25 though for ERCC-00123,
##   and ERCC-00016 is noisy too.  Most of the high fold changes seem to be in
##   "2D" samples with glucose.
##summary(g$data$foldchange)
##   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##  1.000   1.410   2.044   2.796   3.412  26.955

## ERCC with > 500 transcripts used:  Similar to if all used

## ERCC with > 1K transcripts used: Max fold change 11.3, mean=2.5 The 2D
##   glucose samples have higher fold changes esp. the 12h and 24h, even for
##   ERCC-00116 which is 1.4E6 copies. This means the relative abunance of
##   ERCC-00116 (and most of the other ERCC) was higher in those two samples
##   compared to others. So there must have been less community mRNA collected
##   for those two samples.

## ERCC with >10K transcripts: Very similar to the >1K heat map.  max fc = 11.2
##   and mean=2.5.

## Yes, the super high ERCC have fc that are closer to 1.  But even for
## ERCC-00060 (7E5) and more concentrated the fold changes (again, taken across
## samples of the *medians* within each sample, so this should not be noisy) are
## often > 1.



##------------------------------------------------------------------------------
## Differential Expression
##------------------------------------------------------------------------------

BigStep("Differential expression")

## Make the design matrix (similar to MicroTOOLs' MakeLinearModelsForGenes())
## Defintion presumes loading Mari's normd genes.  See fixup code just below.
array2condition <- c("X4h_1glucose"='X4h_glucose', "X4h_1glucose..4array."='X4h_glucose',
                     "X4h_1control"='X4h_control', "X4h_1control..4array."='X4h_control',
                     "X12h_1glucose"='X12h_glucose', "X12h_2glucose"='X12h_glucose',
                     "X12h_control"='X12h_control',
                     "X24h_1glucose"='X24h_glucose',
                     "X24h_1control"='X24h_control',
                     "X2D_4h_1glucose"='X2D_4h_glucose', "X2D_4h_2glucose"='X2D_4h_glucose',
                     "X2D_4h_control"='X2D_4h_control',
                     "X2D_12h_glucose"='X2D_12h_glucose',
                     "X2D_12h_control"='X2D_12h_control',
                     "X2D_24h_glucose"='X2D_24h_glucose',
                     "X2D_24h_control"='X2D_24h_control')
if (!all(names(array2condition) == colnames(eset))) {
    warning('Names in array2condition do not exactly match the eset.  ',
            'If you did your own probe->gene processing, do not worry, I will fix.\n')
    ## FIXME: Not sure about this remapping!!
    x <- names(array2condition)
    x[x=="X4h_1control..4array."] <- 'X4h_2control'
    x[x=="X4h_1glucose..4array."] <- 'X4h_2glucose'
    names(array2condition) <- x
    rm(x)
    stopifnot(colnames(eset) %in% names(array2condition))
    array2condition <- array2condition[colnames(eset)]
}
stopifnot(names(array2condition) == colnames(eset))

fac <- factor(array2condition)
design <- model.matrix(~0 + fac)
colnames(design) <- levels(fac)  # design's cols ordered; this works
## Make the contrasts matrix that describes comparisons of interest.  I was
## surprised to find that even comparisons without replicates are handled by the
## DE code below.  Perhaps b/c I have lots of arrays to make a background, even
## though only 4 sets of replicates.
comparisons <- c(exp1.4h   = 'X4h_glucose-X4h_control',             # has replicates
                 exp1.12h  = 'X12h_glucose-X12h_control',           # has replicates
                 exp1.24h  = 'X24h_glucose-X24h_control',
                 exp2.4h   = 'X2D_4h_glucose-X2D_4h_control',       # has replicates
                 exp2.12h  = 'X2D_12h_glucose-X2D_12h_control')
##              exp2.24h  = 'X2D_24h_glucose-X2D_24h_control')      # Dropped b/c of AQM
contrastMatrix <- makeContrasts(contrasts=comparisons, levels=design)
## Later use the contrastMatrix for EGSEA.

data <- exprs(eset)
## Determine array quality weights.  (optional, could set null)
weights <- arrayWeights(data, design)
## As in MicroTOOLs, do not center the array medians (default).
cat('lmFit:  Making least-squares linear models for genes...')
fit <- lmFit(data, design, weights=weights)
cat('contrasts.fit:  Estimating coefficients and errors based on the models...')
fit <- contrasts.fit(fit, contrastMatrix)  # Get coeffs, stderrs 
## FIXME:  Using the default params for eBayes.  Makes sense?
cat('eBayes:  Looking for differentially expressed genes...')
fit <- eBayes(fit)  # p-values for differential expression
rm(fac,design,comparisons,data,weights)
cat('DE genes were searched for in the following comparisons: ',
    paste(colnames(fit), collapse='; '), '\n')


## Wrapper fro MyTopTable().
## IMPORTANT: MTT() by default uses 0 for the logFC, so it looks for genes with
## any change >0 that is significant (a low bar).  If you pass statFilters$logFC, then
## it will look for genes with a change >logFC that is significant (a high bar).
MariCheckDE <- function(fit, coefs=colnames(fit), number=Inf,
                        statFilters=list(logFC=log2(FC), adj.P.Val=ADJPVAL))
{
    MyTopTable(fit,coefs,number,statFilters,useTopTreat=T)
}

## Separately check for DE genes in each comparison b/c MTT() can't show whether
## a gene was DE in multiple comparisons.  Also fill in a matrix of genes by
## comparisons with the fold changes for significant genes.
## (Multiple testing corrections adjust p for multiple tests in a contrast. There
## is not an additional correction for testing different contrasts. I'm not sure there
## would need to be one. In any case, the separate tests below produce the same 94 genes
## that I would get if I called MariCheckDE() with coefs=colnames(fit).)
deMTTs <- list()
mat <- matrix(NA, nrow=nrow(eset), ncol=length(colnames(fit)),
              dimnames = list(rownames(eset), colnames(fit)))
for (coef in colnames(fit)) {
    mtt <- MariCheckDE(fit, coef)
    deMTTs[[coef]] <- mtt
    if (nrow(mtt) > 0) {
        mat[ rownames(mtt), coef ] <- 2^mtt[,'logFC']  ## Note the 2^. Reporting *fold changes*
    }
}
rm(coef,mtt)

stopifnot(colnames(mat)==colnames(fit))  # cols are conditions compared; rows are genes
cat('\nStarting with a matrix of',nrow(mat),'genes by',ncol(mat),'conditions compared.',
    'Will remove genes that were never DE and comparisons that had no DE genes.\n')
## Remove genes that were never DE.
idx <- which(apply(mat, 1, function(v) !all(is.na(v))))
mat <- mat[idx,]
## Remove comparisons that had no DE genes.
idx <- which(apply(mat, 2, function(v) !all(is.na(v))))
mat <- mat[,idx]
rm(idx)
cat('Now the matrix has',nrow(mat),'genes which were DE in at least one comparison and',
    ncol(mat),'conditions compared: ',  paste(colnames(mat), collapse='; '), '\n')

## Save a csv with DE genes and the annotation.
WriteMatrixCsv(mat, paste0('deGenes.foldchange',FC,'.adjpval',ADJPVAL,'.csv'))

## Heat map of the DE genes
eset.de <- eset[rownames(mat),]
x <- as.character(fData(eset.de)$Strain)
x[grep('^HOT',x)] <- 'HOT clones'
fData(eset.de)$Strain <- factor(x)

cat("\n\nGenerating a heat map of the DE genes, deGenes_heatmap.png.",
    "Samples are bootstrapped 1000x and those with >90% support have a disc.",
    "Genes (rows) are log2 intensities after standardizing (centering and scaling so sd=1)",
    "and the colors indicate positive (red) to negative (blue) sd.",
    "Rowside annotation is shown for the 11 most Strains and Pathways.",
    "The legend shows the number of genes in parenthesis (if >=5, for non-Other).\n")
png("deGenes_heatmap.png",, height=8, width=12, units="in", res=144)
DEGeneHeatmap(eset.de, fit, num=Inf, deGenes=rownames(eset.de),
              colorBy=list(geneIds  = rownames(eset.de),
                           features = c('Ecotype','Strain','Pathway')),
              maxLegendEntries=11,
              bootstrap = c(1000, 90))
dev.off()
rm(x)


##
## Volcano plots for 2D_12h glucose vs. control
##
## The plot will show all detected genes.  The DE ones will be identifiable by
## the added dashed lines.
##
## ~Tricky part.  Important to pass statFilters (via MariCheckDE).  If not, you
## will get different p-values, because topTable()'s [did I mean topTreat?] signif check is w.r.t. the
## specified fold change, which is 0 if you don't specify.)  BUT we also want
## p-vals for all detected genes. So pass adj.P.Val=1.
##
## 28 July 2021: As described in DraftMS/writingUpTranscriptomicResults.txt
## (which has R code), I want to filter out 17 "DE" genes that fell below
## detection cut offs in both of the conditions compared (e.g. 2D_12h_glucose
## vs. control).  This function shows all detected genes, so the un-DE-ified
## ones would still appear above the reference lines to demarcate DE, which
## would be confusing. So I should change this function to show only the genes
## detected in the conditions compared.  Will *not* make that fix here. Copy
## this function to writingUpTranscriptomicResults.txt and rename it
## appropriately.  (Technical hiccup: MariCheckDE uses the 'fit', not deTab, and
## 'fit' cannot/should not be changed to reflect detection in specific
## conditions [I think not feasible to attach condidtions to specific cells in
## the eset passed to limma].).
MakeVolcanoPlot <- function(coef)
{
    df <- MariCheckDE(fit, coef,
                      statFilters=list(logFC=log2(FC), adj.P.Val=1))[,c('logFC','adj.P.Val')]
    stopifnot(setequal(rownames(df),detGenes))
    df <- cbind(df, fData(eset)[rownames(df),])
    df$GID <- rownames(df)
    df$Significance <- -log10(df$adj.P.Val)
    df$FoldChange <- log10(2^df$logFC)
    ## Next line removed because 'unknowns' are now all 'HLunk'.
    ## df$Ecotype <- factor(gsub('^$','Unknown ecotype', as.character(df$Ecotype)))
    df$Ecotype2 <- factor(gsub('^LL.*','Low-light adapted',
                               gsub('^HL.*','High-light adapted',as.character(df$Ecotype))))

    ## Reduce legend clutter. Should prioritize pathways of DE genes.
    numDE <- sort(table(as.character(df[rownames(deMTTs[[coef]]),'Pathway'])),
                  decreasing=T)
    numDE <- numDE[1:min(length(numDE),10)]
    newpaths <- as.character(df$Pathway)
    newpaths[!newpaths %in% names(numDE)] <- 'Other'
    df$Pathway <- factor(newpaths, levels=c(names(numDE),'Other'))
    ## Rename levels to include DE gene counts <-- moved to later
    ##levels(df$Pathway) <- c(paste0(names(numDE),' [',numDE,']'),'Other')

    ## FIXME: The pathways palette should be the same as in the heat map.
    if (FALSE) {
        ## Old way. Inconsistent across volcano plots.
        pal <- c('red','blue','black','purple','darkorange','green','brown','cyan','yellow','darkgreen') ##,'gray')
        pal <- pal[1:length(levels(df$Pathway))]
        names(pal) <- levels(df$Pathway)
    } else {
        ## Hack. Based on DE genes I know.
        pal <- c(## Found in 2D_12h incubations
                 'Respiration'='red',
                 'Glycolysis/Gluconeogenesis'='blue',
                 'Carbon fixation'='black',
                 'Pentose Phosphate'='purple',
                 'GlcH sugar transporter'='darkorange',
                 'Energy metabolism'='lightgreen',
                 'Glycolysis'='brown',
                 'Krebs cycle and related genes'='cyan',
                 'Circadian Rhythm'='yellow',
                 'Entner-Dudoroff'='darkgreen',
                 ##'Other'='gray'
                 ## Added for 24h incubations
                 'Glycolysis/Krebs cycle'='blue',  # recycle Glycolysis/Gluconeogenesis
                 'Cell Division'='pink',
                 'Pentose Phosphate/Entner-Dudoroff'='purple',  # recylce Pentose Phosphate
                 'Photosystem I'='green')
        stopifnot(names(numDE) %in% names(pal))  # can map everything
    }
    ## Make legend include the number of DE genes for each pathway.
    levels(df$Pathway) <- c(paste0(names(numDE),' [',numDE,']'),'Other')
    np <- names(pal)
    idx <- match(np, names(numDE))
    names(pal) <- paste0(np,' [',numDE[idx],']')
    pal <- c(pal, Other='gray')
    
    g <- ggplot(df, aes(x=FoldChange, y=Significance, color=Pathway)) +  #, shape=Ecotype)) +
        geom_point() +
        scale_color_manual(name='Pathway', values=pal) +
        geom_vline(xintercept=c(-log10(FC),+log10(FC)), linetype='dotted') +
        geom_vline(xintercept=0, linetype='solid') +
        geom_hline(yintercept=-log10(ADJPVAL), linetype='dotted') +
        geom_hline(yintercept=0, linetype='solid') +
        labs(y='-log10(adjusted p-value)', x='log10(fold change)') + theme_bw() + 
        facet_wrap(vars(Ecotype2))
}

g <- MakeVolcanoPlot('X2D_12h_glucose-X2D_12h_control')
plot(g)
ggsave('deGenes_volcano_2D_12h.png', width=10, height=4, units='in', dpi=144)
## Brief analysis
gids <- rownames(deMTTs[['X2D_12h_glucose-X2D_12h_control']])  # the 90 that are DE
sort(table(subset(fData(eset)[gids,], Pathway=='Respiration')$Gene.symbol))  # HL/HOT208 coxAB [24]
sort(table(subset(fData(eset)[gids,], Pathway=='Glycolysis')$Gene.symbol))   # HL/HOT208 cbbA [10], pgi [4]
sort(table(subset(fData(eset)[gids,], Pathway=='Pentose Phosphate')$Gene.symbol))  # tal [9] and opcA [2]

## There are 21 DE genes in X24h_glucose-X24h_control. Okay, plot them too.
g <- MakeVolcanoPlot('X24h_glucose-X24h_control')
plot(g)
ggsave('deGenes_volcano_24h.png', width=10, height=4, units='in', dpi=144)


##
## Extra stuff
##
cat('Determining the median log2 transcript levels for pathays in each sample.\n')
cat('Also determining the mean and sd (of the medians) across all samples.\n')
## Better to exponentiate first b/c it does not make sense to do (a+b)/2 for log2 levels.
df <- aggregate(2^exprs(eset), by=list(pathway=fData(eset)$Pathway), function (x) log2(median(x)))
df$mean <- apply(df[,-1],1, function(x) log2(mean(2^x)))
df$sd   <- apply(df[,-1],1, function(x) log2(sd(2^x)))
df$numGenes <- table(fData(eset)$Pathway)[df$pathway]
df <- df[order(df$mean, decreasing=T),]
rownames(df) <- NULL
filename <- 'pathwayMeanLevels.csv'
write.csv(df, filename)
cat('Saved table',filename,'\n')

cat('Determining the median log2 transcript levels for genes in each sample.\n')
cat('Also determining the mean and sd (of the medians) across all samples.\n')
## Better to exponentiate first b/c it does not make sense to do (a+b)/2 for log2 levels.
df <- aggregate(2^exprs(eset), by=list(gene=fData(eset)$Gene.symbol), function (x) log2(median(x)))
df$mean <- apply(df[,-1],1, function(x) log2(mean(2^x)))
df$sd   <- apply(df[,-1],1, function(x) log2(sd(2^x)))
df$numGenes <- table(fData(eset)$Gene.symbol)[df$gene]
df <- df[order(df$mean, decreasing=T),]
rownames(df) <- NULL
filename <- 'geneMeanLevels.csv'
write.csv(df, filename)
cat('Saved table',filename,'\n')


cat('Saving R workspace in microarrays_script.RData\n')
save.image('microarrays_script.RData')
quit(save='no')
