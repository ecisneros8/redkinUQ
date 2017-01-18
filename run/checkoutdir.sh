# variables:
# (1) constrained species
# (2) combination group id
NC=$1
JOB=$2

# logs
logdir="outs/logs/${NC}D/CSI${JOB}"
if [ ! -d "$logdir" ]
then
    echo "$logdir" " created"
    mkdir "$logdir"
else
    if [ "$(ls -A "$logdir")" ]
    then
        echo "$logdir" " already exists (nothing's been deleted)"
        #cd "$logdir"
        #rm *
	#echo "emptied $logdir"
        #cd ../../../../
    fi
fi

# prep
prepdir="outs/prep/${NC}D/CSI${JOB}"
if [ ! -d "$prepdir" ]
then
    echo "$prepdir" " created"
    mkdir "$prepdir"
else
    if [ "$(ls -A "$prepdir")" ]
    then
        echo "$prepdir" " already exists (nothing's been deleted)"
        #cd "$prepdir"
        #rm *
	#echo "emptied $prepdir"
        #cd ../../../../
    fi
fi

# mout (create only, don't delete contents)
moutdir="outs/mout/${NC}D/CSI${JOB}"
if [ ! -d "$moutdir" ]
then
    mkdir  "$moutdir"
fi
