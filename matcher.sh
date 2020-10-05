#!/bin/bash

Usage()
{
  cat <<-END
  ${0##*/}: a tool to study MCH-MFT Track matching

  Aliroot and O2 environments are automatically loaded.

  Usage:
  1) Generate MCH Tracks:

  ${0##*/} --genMCH -n <number_of_events> -o <outputdir> --npions <number_of_pions> --nmuons <number_of_muons>

  Will create outputdir, include a copy of the macros and generate
  MCH tracks with aliroot.


  Option to set number of pions and muons on each events:
  --npions <number_of_pions>
  --nmuons <number_of_muons>

  2) Generate MFT Tracks
     ${0##*/} --genMFT -o <outputdir> -j <jobs>

  3) Run track matching:

  ${0##*/} --match -o <outputdir>

  4) Run checks:

  ${0##*/} --check -o <outputdir>

  ======================================================================================
  Macros on ${SCRIPTDIR} will copied to <outputdir> when --genMCH option is used.
  To replace the macros on --genMFT, --match and --check steps, use option --updatecode

END
  exit
}

updatecode() {
  cp -r ${SCRIPTDIR}/generators/* ${SCRIPTDIR}/*.bin ${SCRIPTDIR}/include ${SCRIPTDIR}/*.C ${SCRIPTDIR}/*.h ${SCRIPTDIR}/*.cxx ${SCRIPTDIR}/macrohelpers ${OUTDIR}
}

generateMCHTracks()
{

  mkdir -p ${OUTDIR}
  updatecode
  pushd ${OUTDIR}

  #sed -i -e s/NPIONS/${NPIONS}/g Config.C
  #sed -i -e s/NMUONS/${NMUONS}/g Config.C

  echo "Generating MCH tracks on `pwd` ..."
  #echo ${MCHGENERATOR}_${NPIONS}pi_${NMUONS}mu_${NEV_}evts  > GENCFG

  ## 1) aliroot generation of MCH Tracks
  echo ${NEV_} > nMCHEvents
  export NEV=${NEV_}
  alienv setenv ${ALIROOTENV} -c bash ./runtest.sh -n ${NEV_} | tee aliroot_MCHgen.log

  ## 2) aliroot conversion of MCH tracks to temporary format
  echo " Converting MCH Tracks to O2-compatible format"
  alienv setenv ${ALIROOTENV} -c aliroot -b -q -l "ConvertMCHESDTracks.C+(\".\")" | tee MCH-O2Conversion.log

  echo " Finished MCH track generation `realpath ${OUTDIR}`"
  popd

}


generateMFTTracks()
{

  if ! [ -f "${OUTDIR}/Kinematics.root" ]; then
    echo " ERROR! MCH Tracks Kinematics.root not found on `realpath ${OUTDIR}/Kinematics.root` ... exiting."
    exit
  fi

  if ! [ -z ${UPDATECODE+x} ]; then updatecode ; fi
  pushd ${OUTDIR}

  echo "Generating MFT Tracks `pwd` ..."

  NEV_=`cat nMCHEvents`
  ## O2 simulation and generation of MFT tracks using same Kinematics.root
  alienv setenv ${O2ENV} -c o2-sim -g extkin --extKinFile Kinematics.root -m PIPE ITS MFT ABS SHIL -e TGeant3 -n ${NEV_} -j $JOBS | tee O2Sim.log
  alienv setenv ${O2ENV} -c o2-sim-digitizer-workflow -b --skipDet TPC,ITS,TOF,FT0,EMC,HMP,ZDC,TRD,MCH,MID,FDD,PHS,FV0,CPV >  O2Digitizer.log
  alienv setenv ${O2ENV} -c o2-mft-reco-workflow -b > O2Reco.log

  echo " Leaving MFT Track generation on `realpath ${OUTDIR}`"
  popd

}

runMatching()
{

  if ! [ -f "${OUTDIR}/tempMCHTracks.root" ]; then
    echo " Nothing to Match... MCH Tracks not found on `realpath ${OUTDIR}/tempMCHTracks.root` ..."
    EXITERROR="1"
  fi

  if ! [ -f "${OUTDIR}/mfttracks.root" ]; then
    echo " Nothing to Match... MFT Tracks not found on `realpath ${OUTDIR}/mfttracks.root` ..."
    EXITERROR="1"
  fi

  if ! [ -z ${EXITERROR+x} ]; then exit ; fi

  if [ -d "${OUTDIR}" ]; then
    if ! [ -z ${UPDATECODE+x} ]; then updatecode ; fi

    pushd ${OUTDIR}
    echo "Matching MCH & MFT Tracks on `pwd` ..."

    ## MFT MCH track matching & global muon track fitting:
    alienv setenv ${O2ENV} -c root.exe -l -q -b runMatching.C+ | tee matching.log

    echo " Finished matching on `realpath ${OUTDIR}`"
    popd

  fi

}


runChecks()
{

  if ! [ -f "${OUTDIR}/GlobalMuonTracks.root" ]; then
    echo " Nothing to check... Global Muon Tracks not found on `realpath ${OUTDIR}/GlobalMuonChecks.root` ..."
    EXITERROR="1"
  fi

  if ! [ -f "${OUTDIR}/mfttracks.root" ]; then
    echo " Nothing to check... MFT Tracks not found on `realpath ${OUTDIR}/mfttracks.root` ..."
    EXITERROR="1"
  fi

  if ! [ -z ${EXITERROR+x} ]; then exit ; fi

  if ! [ -z ${UPDATECODE+x} ]; then updatecode ; fi
  pushd ${OUTDIR}
  echo "Checking global muon tracks on `pwd` ..." && \

  ## Check global muon Tracks
  alienv setenv ${O2ENV} -c root.exe -l -q -b GlobalMuonChecks.C+ | tee checks.log
  echo " Finished checking Global muon tracks on `realpath ${OUTDIR}`"
  popd


}

SCRIPTDIR=`dirname "$0"`


while [ $# -gt 0 ] ; do
  case $1 in
    -n)
    NEV_="$2";
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
    NPIONS_="$2";
    shift 2
    ;;
    --nmuons)
    NMUONS_="$2";
    shift 2
    ;;
    -g)
    GENERATOR="$2";
    shift 2
    ;;
    --genMCH)
    GENERATEMCH="1";
    shift 1
    ;;
    --genMFT)
    GENERATEMFT="1";
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


if [ -z ${GENERATEMCH+x} ] && [ -z ${GENERATEMFT+x} ] && [ -z ${MATCHING+x} ] && [ -z ${CHECKS+x} ]
then
  echo "Missing use mode!"
  echo " "
  Usage
fi


if [ -z ${OUTDIR+x} ]; then echo "Missing output dir" ; Usage ; fi
NEV_=${NEV_:-"10"}
JOBS=${JOBS:-"4"}
GENERATOR=${GENERATOR:-"gun0_100GeV"}
NPIONS_=${NPIONS_:-"10"}
NMUONS_=${NMUONS_:-"2"}

export NPIONS=${NPIONS_}
export NMUONS=${NMUONS_}
export MCHGENERATOR=${GENERATOR}
export ALIROOT_OCDB_ROOT=~/alice/OCDB

#ALIROOTENV=${ALIROOTENV:-"AliPhysics/latest-master-release"}
#ALIROOTENV=${ALIROOTENV:-"AliPhysics/latest-master-next-root6"}
ALIROOTENV=${ALIROOTENV:-"AliPhysics/latest"}
O2ENV=${O2ENV:-"O2/latest-dev-o2"}



if ! [ -z ${GENERATEMCH+x} ]; then
  if [ -d "${OUTDIR}" ]; then
    echo " Warning! `realpath ${OUTDIR}` already exists."
    read -p " Delete output & proceed (y/N)? " choice
    case "$choice" in
      y|Y )
      rm -rf ${OUTDIR}/*.root;
      ;;
      *) exit ;
    esac
  fi
  generateMCHTracks ;
fi

if ! [ -z ${GENERATEMFT+x} ]; then
  generateMFTTracks ;
fi

if ! [ -z ${MATCHING+x} ]; then runMatching ; fi

if ! [ -z ${CHECKS+x} ]; then runChecks ; fi
