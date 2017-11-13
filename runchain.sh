#!/bin/bash

# First make the outgoing plots using getKin_mode
# Assumes the NEUT files are no deeper than 2 dirs down and don't have "outgoing" in name
#for i in $(find ./ -maxdepth 2 -name "*.root" ! -name "*outgoing*" -type f); do
  #( ./getKin_mode $i 1 3 & );
#done

# Should now have the files
# Make the ratios
#cd scripts
#./runme.sh

# And now make the validation plots
# Look for outgoing plots produced by getKin_mode
# But not the ratios, which will have "outgoing*outgoing" in name
for i in $(find ./ -maxdepth 2 -name "*NEUT533*" ! -name "*outgoing*" -type f); do
  # Strip out ./ and here's the NEUT533 file
  neutfile=${i#./*}

  # Now get the validation file
  neut=${neutfile/.root/}_CC1pip_1pi_modes_outgoing
  minoo=${neut/NEUT533\//NEUT_MIN\/}
  minoo=${minoo/\//_}
  minoo=${minoo/NEUT533/MINOO}
  neut=${neut/\//_}
  # Now we have the full validaiton path
  validpath=scripts/${minoo}_${neut}.root
  arr=( ""$i" "1" "3" "$validpath"" )
  #echo $i 
  #echo $validpath
  ( eval ./getKin_mode_validation "$i" "1" "3" "$validpath" & )
done
