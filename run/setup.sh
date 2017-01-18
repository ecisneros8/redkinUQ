# local variables
outdir="outs"
inpdir="inps"

# output directories 
if [ ! -d "$outdir" ]
then
    mkdir "$outdir"
    mkdir "${outdir}/logs"   
    mkdir "${outdir}/mout"    
    mkdir "${outdir}/prep"
fi

# input directory
if [ ! -d "$inpdir" ]
then
    mkdir "$inpdir"
    mkdir "${inpdir}/combs"
fi
cp ../data/* $inpdir

# copy executable
cp ../src/redkin .
co ../src/preproc .



