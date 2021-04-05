#!/bin/bash -l
# The settings come from ~/.bash_profile

Usage()
{
  echo "Usage: ${0##*/} -n <number_of_events> -o <outputdir> "
  echo ""
  exit
}


while [ $# -gt 0 ] ; do
case $1 in
      -n)
      NEV="$2";
      shift 2
      ;;
      -h)
      Usage
      ;;
      *) echo "Wrong input"; Usage;

esac
done

NEV=${NEV:-"10"}

rm -rf *.root *.dat *.log fort* hlt hough raw*  GRP *.ps AliHLT* recraw/*.root recraw/*.log recraw/*.ps
 aliroot -b -q "sim.C(${NEV})"      2>&1 | tee sim.log
mv syswatch.log simwatch.log
 aliroot -b -q rec.C      2>&1 | tee rec.log
mv syswatch.log recwatch.log
# aliroot -b -q ${ALICE_ROOT}/STEER/macros/CheckESD.C 2>&1 | tee check.log
# aliroot -b -q aod.C 2>&1 | tee aod.log

#cd recraw
#ln -s ../raw.root
# aliroot -b -q rec.C      2>&1 | tee rec.log
# aliroot -b -q  2>&1 aod.C | tee aod.log
