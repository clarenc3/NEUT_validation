#!/bin/bash

for i in $(find ../NEUT_MIN -maxdepth 1 -name "*MINOO_SEP2017*outgoing*" -type f); do 
  j=${i/MINOO/NEUT533}
  j=${j/NEUT_MIN/NEUT533}
  root -b -q -l 'makeratio.C("'${i}'", "'${j}'")';
done
