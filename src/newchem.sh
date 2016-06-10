FILES="mech_defs.h GasThermo.h GasKinetics.h"

# get header files from root/mechs/
for fid in $FILES
do
    echo "$1_$fid"
    cp ../mechs/"$1_$fid" "$fid"
done

# get thermo files for CEQ from root/mechs

