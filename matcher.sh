#!/bin/bash

Usage()
{
  cat <<-END
  ${0##*/}: Tool to study MCH-MFT Track matching

  Usage:
  1) Generate common MCH & MFT Tracks:

  ${0##*/} --gen -n <number_of_events> -o <outputdir> -j <jobs>

  Will create outputdir, include a copy of the macros and generate
  MCH tracks with aliroot and MFT tracks with O2. Alienv is automatically loaded.
  Do not run this script with a loaded environment.

  Option to further configure the generator:
  --npions <number_of_pions>
  --nmuons <number_of_muons>

  -g <gen> # selects which generator to use (as defined by folder name in `dirname "$0"`)
     default generator: gun


  2) Run track matching:

  ${0##*/} --match -o <outputdir>

  3) Run checks:

  ${0##*/} --check -n <number_of_events> -o <outputdir> -j <jobs>

END
  exit
}

updatecode() {
  cp -r ${SCRIPTDIR}/${GENERATOR}/* ${SCRIPTDIR}/*.bin ${SCRIPTDIR}/include ${SCRIPTDIR}/*.C ${SCRIPTDIR}/*.h ${SCRIPTDIR}/*.cxx ${SCRIPTDIR}/macrohelpers ${OUTDIR}
}

generateMCHMFTTracks()
{

  mkdir -p ${OUTDIR}
  updatecode
  pushd ${OUTDIR} && \

  sed -i -e s/NPIONS/${NPIONS}/g Config.C
  sed -i -e s/NMUONS/${NMUONS}/g Config.C

  echo "Running on `pwd` ..." && \


  ## 1) aliroot generation of MCH Tracks
  alienv setenv ${ALIROOTENV} -c bash ./runtest.sh -n ${NEV} | tee aliroot_gen.log ;

  ## 2) Generate MFT tracks using same Kinematics.root
  alienv setenv ${O2ENV} -c o2-sim -g extkin --extKinFile Kinematics.root -m PIPE ITS MFT ABS SHIL -e TGeant3 -n ${NEV} -j $JOBS | tee O2Sim.log
  alienv setenv ${O2ENV} -c o2-sim-digitizer-workflow -b --skipDet TPC,ITS,TOF,FT0,EMC,HMP,ZDC,TRD,MCH,MID,FDD,PHS,FV0,CPV >  O2Digitizer.log
  alienv setenv ${O2ENV} -c o2-mft-reco-workflow -b > O2Reco.log

  ## 3) aliroot conversion of MCH tracks to temporary format
  alienv setenv ${ALIROOTENV} -c aliroot -q -l "ConvertMCHESDTracks.C+(\".\")"

  echo " Leaving  ${OUTDIR}"
  popd

}

runMatching()
{

  if [ -d "$OUTDIR" ]; then
    if ! [ -z ${UPDATECODE+x} ]; then updatecode ; fi

    pushd ${OUTDIR} && \
    echo "Running on `pwd` ..." && \


    ## 4) MFT MCH track matching & global muon track fitting:
    alienv setenv ${O2ENV} -c root.exe -l -q -b runMatching.C+ | tee matching.log

    echo " Leaving  ${OUTDIR}"
    popd

  fi
  echo -e "${OUTDIR} not found..."

}


runChecks()
{
  if [ -d "$OUTDIR" ]; then
    if ! [ -z ${UPDATECODE+x} ]; then updatecode ; fi
    pushd ${OUTDIR} && \
    echo "Running on `pwd` ..." && \

    ## 5) Check global muon Tracks
    alienv setenv ${O2ENV} -c root.exe -l -q -b GlobalMuonChecks.C+ | tee checks.log

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
    -g)
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
    --convert)
    CONVERT="1";
    shift 1
    ;;
    --check)
    CHECKS="1";
    shift 1
    ;;
    --updatecode)
    UPDATECODE="1";
    shift 1
    ;;
    -h|--help)
    Usage
    ;;
    *) echo "Wrong input"; Usage;

  esac
done

# Ensure no enviroment is loaded
if ! [[ -z "$LOADEDMODULES" ]]
 then
   echo "Do not run this script with alienv environment loaded. Aborting..."
   echo "Run '${0##*/} --help'"
   exit
 fi


if [ -z ${GENERATE+x} ] && [ -z ${MATCHING+x} ] && [ -z ${CHECKS+x} ]
then
  echo "Missing use mode!"
  echo " "
  Usage
fi


if [ -z ${OUTDIR+x} ]; then echo "Missing output dir" ; Usage ; fi
NEV=${NEV:-"10"}
JOBS=${JOBS:-"4"}
GENERATOR=${GENERATOR:-"gun"}
NPIONS=${NPIONS:-"10"}
NMUONS=${NMUONS:-"2"}

export ALIROOT_OCDB_ROOT=~
SCRIPTDIR=`dirname "$0"`

#ALIROOTENV=AliPhysics/latest-master-release
ALIROOTENV=AliPhysics/latest-master-next-root6
O2ENV=O2/latest-dev-o2


if ! [ -z ${GENERATE+x} ]; then
  if [ -d "$OUTDIR" ]; then
    echo " Warning! `realpath ${OUTDIR}` already exists."
    read -p " Delete output & proceed (y/N)? " choice
    case "$choice" in
      y|Y )
      echo " "
      #rm -rf ${OUTDIR}/*;
      ;;
      *) exit ;
    esac
  fi
  generateMCHMFTTracks ;


fi

if ! [ -d "$OUTDIR" ]; then
  echo "$OUTDIR not found. Nothing to do"
  exit
fi

if ! [ -z ${MATCHING+x} ]; then runMatching ; fi

if ! [ -z ${CHECKS+x} ]; then runChecks ; fi
