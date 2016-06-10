# inputs: 
# name of mechanism
# rcce dimension
mech=$1
dim=$2

# local variables
outdir="outs"
inpdir="inps"
ceqdir="ceqthermo"

# output directories 
if [ ! -d "$outdir" ]
then
    mkdir "$outdir"
    mkdir "${outdir}/logs"
    mkdir "${outdir}/${dim}D"
    mkdir "${outdir}/qois"
    mkdir "${outdir}/rtms"
    mkdir "${outdir}/qois/${dim}D"
    mkdir "${outdir}/rtms/${dim}D"
fi

# input directory
if [ ! -d "$inpdir" ]
then
    mkdir "$inpdir"
fi
cp "../data/${mech}."* "${inpdir}/"

# ceq directory
if [ ! -d "$ceqdir" ]
then
    mkdir "$ceqdir"
fi
cp "../mechs/${mech}_ceqthermo.dat" "${ceqdir}/thermo.dat"
