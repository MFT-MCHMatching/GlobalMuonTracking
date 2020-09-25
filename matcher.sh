#!/bin/bash

Usage()
{
cat <<-END
  ${0##*/}: Tool to study MCH-MFT Track matching

  Usage:
  To generate common MCH & MFT Tracks:

    ${0##*/} --gen -n <number_of_events> -o <outputdir> -j <jobs>
    --npions <number_of_pions>
    --nmuons <number_of_muons>


  To run track matching:

    ${0##*/} --match -o <outputdir>

  To run checks:

     ${0##*/} --check -n <number_of_events> -o <outputdir> -j <jobs>

END
exit
}


generateMCHMFTTracks()
{

  if [ -d "$OUTDIR" ]; then
      echo -e "${OUTDIR} already exists. Aborting."
      exit
  fi

  mkdir -p ${OUTDIR}
  cp -r ${SCRIPTDIR}/${GENERATOR}/* ${SCRIPTDIR}/*.bin ${SCRIPTDIR}/include ${SCRIPTDIR}/*.C ${SCRIPTDIR}/*.h ${SCRIPTDIR}/*.cxx ${SCRIPTDIR}/macrohelpers ${OUTDIR}
  pushd ${OUTDIR} && \

  sed -i -e s/NPIONS/${NPIONS}/g Config.C
  sed -i -e s/NMUONS/${NMUONS}/g Config.C

  echo "Running on `pwd` ..." && \


  ## 1) aliroot generation of MCH Tracks
  alienv setenv AliPhysics/latest-master-next-root6 -c bash ./runtest.sh -n ${NEV} | tee aliroot_gen.log ;

  ## 2) aliroot conversion of MCH tracks to temporary format
  alienv setenv AliPhysics/latest-master-next-root6 -c aliroot -q -l "ConvertMCHESDTracks.C+(\".\")"

  ## 3) Generate MFT tracks using same Kinematics.root
  alienv setenv O2/latest-dev-o2 -c o2-sim -g extkin --extKinFile Kinematics.root -m PIPE ITS MFT ABS SHIL -e TGeant3 -n ${NEV} -j $JOBS | tee O2Sim.log
  alienv setenv O2/latest-dev-o2 -c o2-sim-digitizer-workflow -b --skipDet TPC,ITS,TOF,FT0,EMC,HMP,ZDC,TRD,MCH,MID,FDD,PHS,FV0,CPV >  O2Digitizer.log
  alienv setenv O2/latest-dev-o2 -c o2-mft-reco-workflow -b > O2Reco.log

  echo " Leaving  ${OUTDIR}"
  popd

}

runMatching()
{

  if [ -d "$OUTDIR" ]; then
    pushd ${OUTDIR} && \
    echo "Running on `pwd` ..." && \

    ## 4) MFT MCH track matching & global muon track fitting:
    alienv setenv O2/latest-dev-o2 -c root.exe -l -q -b runMatching.C+

    echo " Leaving  ${OUTDIR}"
    popd

  fi
  echo -e "${OUTDIR} not found..."

}


runChecks()
{
  if [ -d "$OUTDIR" ]; then
  pushd ${OUTDIR} && \
  echo "Running on `pwd` ..." && \

  ## 5) Check global muon Tracks
  alienv setenv O2/latest-dev-o2 -c root.exe -l -q -b GlobalMuonChecks.C+

  echo " Leaving  ${OUTDIR}"
  popd
fi
echo -e "${OUTDIR} not found..."

}


while [ $# -gt 0 ] ; do
case $1 in
      -n)
      NEV="$2";
      shift 2
      ;;
      -j)
      JOBS="$2";
      shift 2
      ;;
      -o)
      OUTDIR="$2";
      shift 2
      ;;
      --npions)
      NPIONS="$2";
      shift 2
      ;;
      --nmuons)
      NMUONS="$2";
      shift 2
      ;;
      --g)
      GENERATOR="$2";
      shift 2
      ;;
      --gen)
      GENERATE="1";
      shift 1
      ;;
      --match)
      MATCHING="1";
      shift 1
      ;;
      --check)
      CHECKS="1";
      shift 1
      ;;
      -h)
      Usage
      ;;
      *) echo "Wrong input"; Usage;

esac
done

if [ -z ${GENERATE+x} ] && [ -z ${MATCHING+x} ] && [ -z ${CHECKS+x} ]
  then
  echo "Missing use mode!"
  echo " "
  Usage
fi


if [ -z ${OUTDIR+x} ]; then echo "Missing output dir" ; Usage ; fi
NEV=${NEV:-"10"}
JOBS=${JOBS:-"4"}
GENERATOR=${GENERATOR:-"gen"}
NPIONS=${NPIONS:-"10"}
NMUONS=${NMUONS:-"2"}

export ALIROOT_OCDB_ROOT=~
SCRIPTDIR=`dirname "$0"`


if ! [ -z ${GENERATE+x} ]; then generateMCHMFTTracks ; fi
if ! [ -z ${MATCHING+x} ]; then runMatching ; fi
if ! [ -z ${CHECKS+x} ]; then runChecks ; fi
