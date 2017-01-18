# variables:
# (1) constrained species
# (2) total number of combinations
NC=7
TNC=116280

# count output directories
outroot="outs/prep/${NC}D/"
ND=$(find $outroot -mindepth 1 -maxdepth 1 -type d | wc -l)
ND=$(($ND-1))
echo "Number of CSI directories: " $ND

# count number of files
outdir="outs/prep/${NC}D/CSI0/"
NF=$(ls $outdir -l | wc -l)
NF=$(($NF-1))
echo "Number of p-rep files: " $NF

# create submission file
WORKDIR=$HOME/project-cse/redkinUQ/run/
echo "#!/bin/bash
#PBS -l walltime=00:05:00
#PBS -l nodes=1:ppn=1
#PBS -N PostProc.${NC}
#PBS -q test 
#PBS -j oe
#PBS -o $WORKDIR/rccepp.${NC}.out
#PBS -e $WORKDIR/rccepp.${NC}.err
cd $WORKDIR
./postproc.exe $NC $ND $NF $TNC
">>ppsub.sh

# run postprocessor
qsub ppsub.sh
rm ppsub.sh

