# working directory:
WORKDIR=$HOME/project-cse/redkinUQ/run

# variables:
# (1) constrained species
# (2) number of nodes
# (3) processors per node
# (4) combination file id
NC=$1
NN=$2
PPN=$3
JOB=$4

# remove file if it exists
if [ -f submit.sh ]
then
    rm submit.sh
fi

echo "#!/bin/bash
#PBS -l walltime=00:05:00
#PBS -l nodes=${NN}:ppn=${PPN}
#PBS -N RCCE.${NC}D.${JOB} 
#PBS -q xpacc 
#PBS -j oe
#PBS -o $WORKDIR/outs/prompt/${NC}D/rcce.CSI${JOB}.out
#PBS -e $WORKDIR/outs/prompt/${NC}D/rcce.CSI${JOB}.err
cd $WORKDIR
mpiexec ./redkin $NC $JOB
">>submit.sh
