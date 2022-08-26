#!/bin/bash

## Run microarrays.R. Then run egsea.R. Put files in the specified Output
## directory.

OUTDIR=$1
if [ -z "$OUTDIR" ] || [ -d "$OUTDIR" ] ; then
    echo "Pass an output directory that does exist."
    exit -1
fi

echo
echo "Will run microarrays.R and egsea.R.  Results will be in $OUTDIR."
echo "Separate logs will be created for each R script."
echo "I will be quiet now.  You should put me me in the background and"
echo "and periodically check the logs to see how things are going."

## OUTDIR must be simple because we will use ../scripts/.  mkdir will fail if
## OUTDIR is foo/bar/baz.
mkdir $OUTDIR
if [ "$?" -ne 0 ] ; then
    echo "You must pass a simple output directory name, not a path."
    exit -1
fi
CWD=`pwd`
cd $OUTDIR
mkdir Tables Plots Plots/EGSEA
## ~Hack.  Help the script find these files.
ln -s ../DropBoxImport/Transcriptomic\ analysis/UseThisFigS4forFeatureData.csv .
ln -s ../ercc.meta.fromMari.csv .
../scripts/microarrays.R &> log.microarrays.txt
rm UseThisFigS4forFeatureData.csv ercc.meta.fromMari.csv
mv *.png Plots
mv *.csv Tables

../scripts/egsea.R &> log.egsea.txt
mv egsea*.png Plots/EGSEA
mv egsea*.csv Tables

rm Rplots.pdf cachedProbeToGeneMap.RData
mkdir RData
mv microarrays_script.RData egsea_script.RData  RData
cd $CWD

exit 0
