# working directory
WORKDIR=$HOME/project-cse/redkinUQ/run/

# variables:
# (1) constrained species
# (2) number of nodes
# (3) processors per node
# (4) job id
NC=7
NN=10
PPN=24
JOB=476

# total number of processors
NP=$(($NN*$PPN))

# run the rest of the jobs:
sh createscript.sh $NC $NN $PPN $JOB
#sh checkoutdir.sh $NC $JOB
qsub submit.sh
rm submit.sh
