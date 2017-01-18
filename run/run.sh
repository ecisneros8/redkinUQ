# working directory
WORKDIR=$HOME/project-cse/redkinUQ/run/

# variables:
# (1) constrained species
# (2) number of nodes
# (2) processors per node
NC=7
NN=10
PPN=24

# total number of processors
NP=$(($NN*$PPN))

# preprocessor:
if [ "$(ls -A "inps/comb/${NC}D/")" ]
then
    rm "inps/comb/${NC}D/*"
    echo "emptied inps/comb/${NC}D/..."
else
    echo "inps/comb/${NC}D/ is empty..."
fi
./preproc $NC $NP

# get number of jobs:
#while read line
#do
#  NJ=$line
#  echo "The total number of jobs to submit is ${NJ}"
#done < jobs.sh
#rm jobs.sh


# run the rest of the jobs:
#JOB=0
#while [ $JOB -lt $NJ ]
#do
#    sh createscript.sh $NC $NN $PPN $JOB
#    sh checkoutdir.sh $NC $JOB
#    qsub submit.sh
#    echo "Submitted job ${JOB} out of ${NJ}"
#    let JOB=JOB+1
#done
