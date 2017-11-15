#!/bin/bash

# First make the outgoing plots using getKin_mode
RunAnal () {

  # First produce the outgoing plots for NEUT533 and Minoo
  neut533=${1#./}
  echo "Running getKin_mode on NEUT533: $neut533..."
  sleep 3
  ./getKin_mode ${neut533} 1 3

  # Then run getKin on Minoo
  minoo=${neut533/NEUT533\//NEUT_MIN\/}
  minoo=${minoo/_NEUT533_/_MINOO_}
  echo "Running getKin_mode on MINOO: $minoo..."
  sleep 3
  ./getKin_mode $minoo 1 3

  # Now have plots called outgoing
  outgoing=CC1pip_1pi_modes_outgoing.root
  outgoing_533=${neut533}_${outgoing}
  outgoing_minoo=${minoo}_${outgoing}
  echo "Running makeratio on ${outgoing_minoo} ${outgoing_533}..."
  sleep 3
  root -b -q -l 'makeratio.C("'${outgoing_minoo}'", "'${outgoing_533}'")';

  # Construct the validation string
  minoosub=${minoo/\//_}
  minoosub=${minoosub%*.root}_${outgoing%*.root}
  neutsub=${neut533/\//_}
  neutsub=${neutsub%*.root}_$outgoing

  # Then we have the weightfile (contains ratio)
  weightfile=${minoosub}_${neutsub}
  echo "Running validation check on ${neut533} with weight file ${weightfile}..."
  sleep 3
  ( eval ./getKin_mode_validation "${neut533}" "1" "3" "${weightfile}" )

  # Finally make the ratio
  # Get the original Minoo outgoing
  # Get the reweighted NEUT533 validation file
  neut533_valid=${neut533}_${outgoing%.root}_validation.root
  
  echo "Making ratio of Minoo to reweighted NEUT533 using ${outgoing_minoo} and ${neut533_valid}"
  sleep 3
  root -b -q -l 'makeratio.C("'${outgoing_minoo}'", "'${neut533_valid}'")'

  echo -e "\n"
}

# Do analysis on NEUT533 files and then run on equivalent Minoo file in the same script
# Allows for multi-threading since only dependency is that NEUT533 and MINOO file has been produced so we can make the ratio
for i in $(find ./ -maxdepth 2 -name "*NEUT533*.root" ! -name "*outgoing*" -type f); do
  ( RunAnal $i & )
done
