## DEAR READER, this script is run interactively from the R console (not by runAllScripts.sh).

## The denoised data has 1657 "features" (distinct sequences) that represent
## 2,299,435 total sequences.  Can see these in the Denoise/JM_esxtracted HTML
## report [from STABvida].

## This must be how STABvida used QIIME2:
##     https://docs.qiime2.org/2020.8/plugins/available/dada2/denoise-paired/
## which seems to say that the joined paired-end sequences are "features".
## QIIME2 denoising via DADA2 produces ASV's, as I read "Grouping similar sequences" here:
##   https://docs.qiime2.org/2020.8/tutorials/qiime2-for-experienced-microbiome-researchers/
## and QIIME2 should produce a table of feature (ASV) observations in each sample.
## That's what I need.
##

## Aha! This seems to be the QIIME2 feature table, produced after/via DADA2
## denoising.
datapath <- file.path('data','16S_STABvida','Statistics')
featureTableFile <- file.path(datapath, 'feature-table.txt')

## Load ASV table.  Sample names take from first line in the file (#commented)
asvTab <- read.table(featureTableFile)
colnames(asvTab) <- c('ASV',
                      "11-N0023", "15-N0025", "19-N0027",  "23-N0028", "27-N0029",
                      "28-N0030", "3-N0021",  "31-N0031", "35-N0033", "36-N0034",
                      "39-N0035", "40-N0036", "43-N0037",  "44-N0038", "47-N0039",
                      "48-N0040", "7-N0022",  "Croco-N0020")
stopifnot(nrow(asvTab)==1657)

## Samples have quite different total amplicons.
summary(colSums(asvTab[,-1]))
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   29532  123790  134616  127746  140333  156799 


## Sample metadata.
sampleMetadataFile <- file.path(datapath, 'sample-metadata.tsv')
sampMetaTab <- read.table(sampleMetadataFile, header=F, skip=2, sep="\t")
## The final column is really called "Description" and it always starts with the
## "Sample" and the ID (e.g. 3-N0021) which is redundant.  If I chop off that
## part then I can see which samples are technical replicates, as Mari indicated
## in Nov 16 email.  So Description --> SampType
colnames(sampMetaTab) <- c("SampleID","Time","DayCycle","Replicate","Treatment_Glu",
                           "DataSet","SampType")
rownames(sampMetaTab) <- sampMetaTab$SampleID
sampMetaTab$SampType <- factor(gsub('^Sample[^_]+_','',sampMetaTab$SampType))
## Also have SID that is just the leading number part (or "Croco") of the
## SampleID, since it uniquely identifies the SampleID.
sampMetaTab$SID <- gsub('-N0.*$','',sampMetaTab$SampleID)


library(vegan) # for Bray-Curtis, and metaMDS
library(ggplot2)
library(ggrepel) # for non-overlapping labels
library(reshape2)
library(RColorBrewer)

## Palette with distinguishable colors (mainly).
##    display.brewer.all()
MakeColorVec <- function(namVec)
{
    pcolors <- unique(c(brewer.pal(9,'Set1'),
                        rev(brewer.pal(7,'Set2')),
                        brewer.pal(8,'Set3')[-c(2,12,6,8,9,1)]))
    stopifnot(length(pcolors) >= length(namVec))
    pcolors <- pcolors[1:length(namVec)]
    names(pcolors) <- namVec
    pcolors
}


## Adapted from Ana's 18S project
## Uses global sampMetaTab
## NOTE that this is how Bray and Canberra are implemented (in vegdist):
##  ‘canberra’    d[jk] = (1/NZ) sum (abs(x[ij]-x[ik])/(abs(x[ij])+abs(x[ik])))
##                       where NZ is the number of non-zero entries.
##  ‘bray’        d[jk] = (sum abs(x[ij]-x[ik]))/(sum (x[ij]+x[ik]))
##
## For my count data the abs() in the denom of canberra has no affect, so the
## 1/NZ in canberra is the only difference.  It will boost the similarity of
## samples that have lots of species in common (at non-zero).
##
## Note that 'bray' numertor ssems not restricted to just the species
## in both samples, as Wikipedia defines Bray-Curtis.
##
DoNMDS <- function(dat, method='bray', sampScale=0, title=NULL, log2=FALSE,
                   rarefy=0)  # Don't use with sampScale!
{
    if (sampScale > 0) {
        dat <- scale(dat, center=F, scale=colSums(dat)/sampScale)
    }
    if (log2) {
        dat <- log2(dat+1)
    }
    if (rarefy > 0) {
        dat <- as.data.frame(t(rrarefy(t(dat), rarefy))) # should warn if rarefy is too big
    }
    ## Use isoMDS only b/c below I divide stress by 100.  If used monoMDS then I
    ## think I do not have to.
    nmds <- metaMDS(t(dat), distance=method, k=2, engine='isoMDS')
    rnams <- rownames(nmds$points)
    df = data.frame(nmds$points,
                    Sample  = factor(sampMetaTab[rnams,'SID']),
                    Glucose = factor(sampMetaTab[rnams,'Treatment_Glu'], levels=c('+','-')),  #~red,~blue
                    Time    = factor(sampMetaTab[rnams,'Time']),
                    Day     = factor(sampMetaTab[rnams,'DayCycle']),
                    DataSet = factor(sampMetaTab[rnams,'DataSet']),
                    SampType= factor(sampMetaTab[rnams,'SampType']))
    colnames(df) <- c('x','y','Sample','Glucose','Time','Day','Dataset','SampType')
    method <- c(bray='Bray-Curtis dissimilarity',
                canberra='Canberra distance',
                gower='Gower',
                cao='Cao',
                jaccard='Jaccard')[method]  # Jaccard is ~same as Bray in metaMDS()
    sampScale <- c('','normalized ')[(sampScale > 0)+1]
    gt <- ggtitle(paste0("NMDS: ",method," between ",sampScale,"sample ASV abundances\n(stress=",
                         round(nmds$stress/100, 2),")"))
    ##g <- ggplot(df, aes(x,y,color=Glucose,shape=Day)) + geom_point() +
    ##g <- ggplot(df, aes(x,y,color=Time,shape=Dataset)) + geom_point() +
    g <- ggplot(df, aes(x,y,color=SampType,shape=Glucose)) + geom_point(size=6) +
        xlab(NULL) + ylab(NULL) + gt + coord_fixed() +
        geom_text_repel(aes(label=Sample), size=3) +  # Yay!
        theme_bw() + theme(panel.grid.major = element_blank(),
                           panel.grid.minor = element_blank()) +
        scale_color_manual(values = MakeColorVec(levels(df$SampType))) +
        scale_shape_manual(values=c('+'='+','-'='-'))
    list(g=g, nmds=nmds)
}

## First checks
set.seed(12345)
dat <- asvTab[,-c(1,which(colnames(asvTab)=="Croco-N0020"))]
NMDS <- list()
NMDS$Bray <- DoNMDS(dat)
NMDS$Canberra <- DoNMDS(dat, 'canberra')
NMDS$Bray.norm100K <- DoNMDS(dat, sampScale=100000)
NMDS$Canberra.norm100K <- DoNMDS(dat, 'canberra', sampScale=100000)
NMDS$Bray.rarefy100K <- DoNMDS(dat, rarefy=100000)
NMDS$Canberra.rarefy100K <- DoNMDS(dat, 'canberra', rarefy=100000)


## Does normalizing the sample sizes make a difference?
##   Bray-Curtis:  No.  Get the same plot (maybe flipped of course). Stresses are very similar.
##   Canberra:     Same as above.

## Bray-Curtis vs. Canberra?: Compared normalized NMDS's.  Similar
## stresses. Same successful replicates (27,28; 39,40).  Same outliers
## (35,7,44). Same overall positions.


## Limit to organisms that are >0.1% of any sample, and that
## are have >10 reads in at least 5 of the 17 samples
dat.pct <- scale(dat, center=F, scale=colSums(dat)/100)
idx <- which(apply(dat.pct,1, function(rv) any(rv>0.1)))  # 212 if 0.1%
## Now filter for denovo's with >10 reads in >=5 samples.
idx <- intersect(idx,
                 which(apply(dat, 1, function(vec) { sum(vec > 10) >= 5 })))  # down to 202
NMDS.filt <- list()
NMDS.filt$Bray <- DoNMDS(dat[idx,])
NMDS.filt$Canberra <- DoNMDS(dat[idx,], 'canberra')
NMDS.filt$Bray.norm100K <- DoNMDS(dat[idx,], sampScale=100000)  # Very similar to first check
NMDS.filt$Canberra.norm100K <- DoNMDS(dat[idx,], 'canberra', sampScale=100000)

## Does normalizing the sample sizes make a difference?
##   Bray-Curtis:  No.  Tiny stress for both. Same outliers. Same plot, just flipped.
##   Canberra:     Plots look identical (no flip, lucky me). Same stress (0.08).

## Note that Canberra (in vegdist) weights distances by the numer of nonzero
## entries, i.e. samples that _have_ the same species will be more similar
## (vs. similarity due based on lots of 0's in both samples).

## Bray-Curtis vs. Canberra?: Compared normalized NMDS's.  Plots are ~identical
## and stresses are 0.07 and 0.08.


## Filtering vs. not.  Just checking Bray.norm10K's.
## - Stress is twice as much in non-filtered (0.14) vs. in filtered (0.07).
## - Outliers are the same (7,35,44).
## - Good replicates are the same (39,40; 27,28)
## - Overall positions pretty similar.  Better spread in the non-filtered.
##   Some changes in closest neighbors (e.g. see 43 and 47).

## Main conclusions on method of NMDS:
##   - dissimilarity measure doesn't matter for the 2 I tried
##   - normalizing samples to have 100K reads made no difference
##   - filtering out rare ASV's halved the stress

## Experimental observations:
##   (1) 12_2day_-Glu replicates are great (39,40)
##   (2) 4h_2day_+Glu replicates are great (27,28) and similar to (39,40)
##
##   (3) Bad replicates:  12h_2day_+Glu replicates are very different (35,36)
##                        24h_2day_+Glu replicates are quite different (43,44)
##                        [ *Not* yellowish 19 and 3, they are not replicates. ]
##   (4) Outliers with no replicates:  4h_1day_-Glu (7)
##
##   (5) Overall I don't see clusters by + vs. - Glu.
##       The main cluster is 39,40, 27,28 which are all {12,4}h_2day_{-,+}Glu.  But these
##         are not all the 2day's.
##       The {12,24}h_1day_{+,-}Glu samples form a loose cluster (19,23,15,11), but the 4h_1day
##       samples are ~outliers (7,3[light yellow])
##       ==> Some clustering of samples by the day, across hour and Glu status.


## A few more checks.
#NMDS.filt$Cao            <- DoNMDS(dat[idx,], 'cao')
NMDS.filt$Gower          <- DoNMDS(dat[idx,], 'gower')
#NMDS.filt$Cao.norm100K   <- DoNMDS(dat[idx,], 'cao', sampScale=100000)
NMDS.filt$Gower.norm100K <- DoNMDS(dat[idx,], 'gower', sampScale=100000)



Nor by DayCycle because
##       39,40 and 27,28 are from different DayCycle's (1 and 2 resp.).
##       MAYBE some time patterns since 23 and 15 are 24h and 12h from day 1 (with
##        different Glu's); and 
