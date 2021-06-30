#!/bin/bash

MATCHINGRESULTS="GlobalMuonTracks.root matching.log MatchingPlane_eV*.png ML_Evaluation*.root"
CHECKRESULTS="GlobalMuonChecks.root checks.log"

Usage()
{
  cat <<-END
  ${0##*/}: a tool to study MCH-MFT Track matching

  Aliroot and O2 environments are automatically loaded.

  Usage:
  1) Generate MCH Tracks:
     Will create outputdir, include a copy of the macros and generate MCH tracks with aliroot.

  ${0##*/} --genMCH -n <number_of_events> -o <outputdir> -g <generator> <generator options>

    -g
     Sets one or more generators for MCH and MFT simulations. Options:

      gun0_100GeV - Box generator for pions and muons with total momentum 0 to 100 GeV. (default)
                    Set number of pions and muons on each event:
          --npions <number_of_pions>
          --nmuons <number_of_muons>

      MuBoxGun - Box generator for muons with total momentum 0 to 100 GeV.
                       Set number of muons on each event:
           --nmuons <number_of_muons>

      PiMuParam - AliGenParam pions and muons generator with realistic parametrized distributions.
                  Set number of pions and muons on each event:
          --npions <number_of_pions>
          --nmuons <number_of_muons>

      PiParam - AliGenParam pions generator with realistic parametrized distributions.
                  Set number of pions and muons on each event:
          --npions <number_of_pions>

      hijing  - Mimic Hijing generator for particle background corresponding to the 0-10 most central PbPb collisions

      dimuon  - Generate dimuons decays from J/Psi and/or Upsilon. At least one kind of mother particle must be defined.
                 To set the number of mother particles on each event:
         --njpsis <number of J/Psis>
         --nupsilons <number of Upsilons>

    Examples:
    ${0##*/} --genMCH -g gun0_100GeV -n 20 --nmuons 4 --npions 20 -o sampletest
    ${0##*/} --genMCH -g hijing -g dimuon --njpsis 10 --nupsilons 10 -n 20 -o sampletest

  2) Generate MFT Tracks:
     ${0##*/} --genMFT -o <outputdir> -j <jobs>

      Options:
       --shm-segment-size  <size_in_bytes> (default 5000000000)
         sets shm-segment-size for o2 digitizer and reconstruction workflows. Default O2 value of 1GB breaks at ~2e6 MFT tracks. 5GB for 3e6 MFT tracks.

  3) Run track matching:
     ${0##*/} --match --matchFcn <matching_function> --cutFcn <cut_function> --cutParam0 <val0> -o <outputdir>

     --matchPlaneZ
       Sets the z position of the matching plane; MCH and MFT tracks are propagated to this plane for matching evaluation.

     --InitMFTTracksFromVertexingParameters
       Initializes MFT tracks from the vertexing parameters: tracks are propagated to the matching plane starting from the main vertexing parameters.
       MCS effects due to MFT disks are added to the covariances matrix. 

     --matchSaveAll
       Save all MCH/MFT track combinations, not only the best match

     --matchFcn
       Sets the function to calculate matching score for a MCH-MFT track pair.  Built-in options:

         matchALL -  matching chi2 calculated by all 5 parameters: X, Y, Phi, Tanl, q/Pt; (default)

         matchXYPhiTanl - matching chi2 calculated by position and angles: X, Y, Phi, Tanl

         matchXY - matching chi2 calculated track positions: X, Y

         matchHiroshima - Hiroshima matching function

         trainedML - matching using a trained neural network (see TMVA interface bellow)

     --cutFcn
       Sets the function that defines the search window for matching candidates. Built-in options:

        cutDisabled - all MFT tracks tested for every MCH track; (default)

        cutDistance - MFT candidates are limited by a distance on the XY plane as defined by cutParam0 in cm;

        cutDistanceSigma - MFT candidates are limited by a distance on the XY plane as defined by
                           mCutParams[0]*TMath::Sqrt(mchTrack.getSigma2X()+mchTrack.getSigma2Y());

        cutDistanceAndAngles - MFT candidates are cut by
                                1) distance on the XY plane as defined by cutParam0 in cm;
                                2) Delta Phi (direction of p_t) by cutParam1 in radians;
                                3) Delta Theta (track polar angle) by cutParam2 in radians;
                                ** Use --cutParamN (see bellow)

        cutDistanceAndAngles3Sigma - MFT TDR cut (Section 6.5)

        cut3Sigma - 3 Sigma cut for all parameters.

        cutNSigma - N Sigma cut for all parameters. N msut be set with `--cutParam0 <N>`.


        cutDistanceAndAnglesVar - Cut based on observed MCH residuals variances

     --cutParam0 <val0>
       Sets mCutParams[0]=val0; (double)

     --cutParam1 <val1>
       Sets mCutParams[1]=val1; (double)

     --cutParam2 <val2>
       Sets mCutParams[2]=val2; (double)

     --enableChargeMatchCut
        Enables charge match cut (which is disabled by default)

     Example:
     ${0##*/} --match --matchFcn matchXYPhiTanl --cutFcn cutDistance --cutParam0 2.0 -o sampletest

  4) Run checks:

  ${0##*/} --check -o <outputdir>

  Other options:
     --verbose
       Enable verbose matching and checking output

  ======================================================================================
  The following output files are copied to a results subdirectory:
   Step 3) Matching:
   $MATCHINGRESULTS

   Step 4) Checking:
   $CHECKRESULTS
  ======================================================================================
  Contents of ${SCRIPTDIR} are copied to <outputdir> when --genMCH option is used.
  To replace the macros on --genMFT, --match and --check steps, use option --updatecode

  ======================================================================================
  Machine learning interface - ROOT TMVA

  1) Generate training data file:
    ${0##*/} --exportTrainingData NMCH_Tracks -o <outputdir>
     Creates a traning data root file. Each entry o on the tree contains
       40 parameters and 1 truth.
     Track-pairs can be selected by the matching cut functions. See option --cutFcn.
     Correct matches can be forced into training data file with --CorrectMatchIgnoreCut


     Example:
     ${0##*/} --exportTrainingData 42 --cutFcn cutDistance --cutParam0 2.0 -o outputdir

  2) Train neural network:
    ${0##*/} --train <Training_Method> <config_alias/es> --trainingdata <training_data_file.root>

     Example:
     ${0##*/} --train DNN --layout DL4.2 --strategy ts1 --MLoptions oo1 --trainingdata MLTraining_1000_MCHTracks.root -o outputdir

      This creates the Trained Network file "Trained_ML_<config_alias>_<training_data_file>.weights.xml" in the folder "outputdir/trainedMLs/weights/"
      Existing aliases are in ML configuration file "MLConfigs.xml". One can use any number of the three available. If none are set, then default values from TMVA will be used. The only restrition is that within the options of a method, options cannot have the same name, e.g., a layout and a strategy called "example1".
      Training method is mandatory.

  3) Run track-matching using trained network:
    ${0##*/} --match --matchFcn trainedML --weightfile weightfilename.xml -o outputdir
     Runs track-matching as any matching function with a default machine learning
     score cut of 0.5 (see --MLScoreCut). Also generates ML_Evaluation.root with
     basic performance assessment of the method.

     --MLScoreCut <val0>
     Configures matching score cut to select the final GlobalMuonTracks.root. (default: 0.5)

     Example:
     ${0##*/} --match --matchFcn trainedML --weightfile trainedML/weights/Regression_DNN_DL4.2_ts1_oo1__MLTraining_1000_MCHTracks.weights.xml -o outputdir


 ======================================================================================
  Machine learning interface - Python

  1) Generate training data file:
    ${0##*/} --exportTrainingData NMCH_Tracks -o <outputdir>
     Creates a traning data root file. Each entry o on the tree contains
       40 parameters and 1 truth.
     Track-pairs can be selected by the matching cut functions. See option --cutFcn.
     TODO: Configurable training data format

     Example:
     ${0##*/} --exportTrainingData 42 --cutFcn cutDistance --cutParam0 2.0 -o outputdir

  2) Train XGBoost:
    ${0##*/} --train --onPythonML --trainingdata <training_data_file.root>

     Example:
     matcher.sh --train --onPythonML --trainingdata MLTraining_1000_MCHTracks.root -o outputdir
END
  exit
}

updatecode() {
  cp -r ${SCRIPTDIR}/generators/* ${SCRIPTDIR}/*.bin ${SCRIPTDIR}/*.xml ${SCRIPTDIR}/include ${SCRIPTDIR}/*.C ${SCRIPTDIR}/*.h ${SCRIPTDIR}/*.cxx ${SCRIPTDIR}/*.py ${SCRIPTDIR}/macrohelpers ${OUTDIR}
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
  export SEED=${SEED:-"123456"}
  echo ${NEV_} > nMCHEvents
  export NEV=${NEV_}
  rm -rf MatcherGenConfig.txt
  alienv setenv ${ALIROOTENV} -c bash ./runtest.sh -n ${NEV_} | tee aliroot_MCHgen.log

  ## 2) aliroot conversion of MCH tracks to temporary format
  echo " Converting MCH Tracks to O2-compatible format"
  alienv setenv ${ALIROOTENV} -c aliroot -e 'gSystem->Load("libpythia6_4_25")' -b -q -l "ConvertMCHESDTracks.C+(\".\")" | tee MCH-O2Conversion.log
  popd
  echo " Finished MCH track generation `realpath ${OUTDIR}`"

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
  alienv setenv ${O2ENV} -c o2-sim-digitizer-workflow ${CUSTOM_SHM} -b --skipDet TPC,ITS,TOF,FT0,EMC,HMP,ZDC,TRD,MCH,MID,FDD,PHS,FV0,CPV >  O2Digitizer.log
  alienv setenv ${O2ENV} -c o2-mft-reco-workflow ${CUSTOM_SHM} -b > O2Reco.log
  popd
  echo " Finished MFT Track generation on `realpath ${OUTDIR}`"

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
    alienv setenv ${O2ENV} -c root.exe -e 'gSystem->Load("libO2MCHTracking")' -l -q -b runMatching.C+ | tee matching.log
    RESULTSDIR="Results`cat MatchingConfig.txt`"
    mkdir -p ${RESULTSDIR}
    cp ${MATCHINGRESULTS} "${RESULTSDIR}"

    popd
    echo " Finished matching on `realpath ${OUTDIR}`"

  fi

}

exportMLTrainningData()
{

  if ! [ -f "${OUTDIR}/tempMCHTracks.root" ]; then
    echo " Nothing to export... MCH Tracks not found on `realpath ${OUTDIR}` ..."
    EXITERROR="1"
  fi

  if ! [ -f "${OUTDIR}/mfttracks.root" ]; then
    echo " Nothing to export... MFT Tracks not found on `realpath ${OUTDIR}` ..."
    EXITERROR="1"
  fi

  if ! [ -z ${EXITERROR+x} ]; then exit ; fi

  if [ -d "${OUTDIR}" ]; then
    if ! [ -z ${UPDATECODE+x} ]; then updatecode ; fi

    pushd ${OUTDIR}
    echo "Exporting ML Traning data file on `pwd` ..."
    ## MFT MCH track matching & global muon track fitting:
    alienv setenv ${O2ENV} -c root.exe -e 'gSystem->Load("libO2MCHTracking")' -l -q -b runMatching.C+ | tee training_data_gen.log
    RESULTSDIR="MLTraining`cat MatchingConfig.txt`"
    mkdir -p ${RESULTSDIR}
    cp training_data_gen.log MLTraining_*.root "${RESULTSDIR}"

    popd
    echo " Finished exporting ML Traning data. File saved on `realpath ${RESULTSDIR}`"
  fi


}


trainML()
{

  if ! [ -f "MLConfigs.xml" ]; then
    echo " Machine Learning configuration file absent..."
    cp MLConfigs.xml "${OUTDIR}"
  fi

  if ! [ -f "${ML_TRAINING_FILE}" ]; then
    echo " ERROR: could not open data file! "
    exit
  fi

  export ML_TYPE=${ML_TYPE:-"Regression"}
  if [ $ML_TEST ]; then
      export ML_NTEST=${ML_NTEST:-"0.1"}
  fi

  if [ -d "${OUTDIR}" ]; then
      if ! [ -z ${UPDATECODE+x} ]; then updatecode ; fi
      pushd ${OUTDIR}

      if ! [ -z ${ML_PYTHON+x} ]; then
	      alienv setenv ${O2ENV} -c python3 python_training.py | tee MLtraining.log
      else
	      alienv setenv ${O2ENV} -c root.exe -e 'gSystem->Load("libO2MCHTracking")' -l -q -b runMatching.C+ | tee MLtraining.log
      fi
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
  RESULTSDIR="Results`cat MatchingConfig.txt`"
  mv ${MATCHINGRESULTS} ${CHECKRESULTS} "${RESULTSDIR}"
  echo " Results moved to `realpath ${RESULTSDIR}`"
  popd
  echo " Finished checking Global muon tracks on `realpath ${OUTDIR}`"

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
    export NPIONS="$2";
    shift 2
    ;;
    --nmuons)
    export NMUONS="$2";
    shift 2
    ;;
    --njpsis)
    export NJPSI="$2";
    shift 2
    ;;
    --nupsilons)
    export NUPSILON="$2";
    shift 2
    ;;
    -g)
    if [ -z ${GENERATOR+x} ]
    then
      GENERATOR="$2";
    else
      GENERATOR="${GENERATOR}.$2";
    fi
    shift 2
    ;;
    --seed)
    export SEED="$2";
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
    --shm-segment-size)
    CUSTOM_SHM="--shm-segment-size $2";
    shift 1
    ;;
    --match)
    MATCHING="1";
    shift 1
    ;;
    --matchPlaneZ)
    export MATCHING_PLANEZ="$2";
    shift 2
    ;;
    --InitMFTTracksFromVertexingParameters)
    export INIT_MFT_FROM_VERTEXING="1";
    shift 1
    ;;
    --matchSaveAll)
    export MATCH_SAVE_ALL="1";
    shift 1
    ;;
    --matchFcn)
    export MATCHING_FCN="${2}_";
    shift 2
    ;;
    --cutFcn)
    export MATCHING_CUTFCN="${2}_";
    shift 2
    ;;
    --cutParam0)
    export MATCHING_CUTPARAM0="$2";
    shift 2
    ;;
    --cutParam1)
    export MATCHING_CUTPARAM1="$2";
    shift 2
    ;;
    --cutParam2)
    export MATCHING_CUTPARAM2="$2";
    shift 2
    ;;
    --enableChargeMatchCut)
    export ENABLECHARGEMATCHCUT="1";
    shift 1
    ;;
    --CorrectMatchIgnoreCut)
    export ML_CORRECTMATCHIGNORECUT="1";
    shift 2
    ;;
    --exportTrainingData)
    export ML_EXPORTTRAINDATA="$2";
    shift 2
    ;;
    --onPythonML)
    export ML_PYTHON="1";
    shift 1
    ;;
    --weightfile)
    export ML_WEIGHTFILE="`realpath $2`";
    shift 2
    ;;
    --MLScoreCut)
    export ML_SCORECUT="$2";
    shift 2
    ;;
    --train)
    export TRAIN_ML_METHOD="$2";
    shift 2
    ;;
    --mltest)
    export ML_TEST="1";
    shift 1
    ;;
    --ntest)
    export ML_NTEST="$2";
    shift 2
    ;;
    --type)
    export ML_TYPE="$2";
    shift 2
    ;;
    --layout)
    export ML_LAYOUT="$2";
    shift 2
    ;;
    --strategy)
    export ML_TRAINING_STRAT="$2";
    shift 2
    ;;
    --MLoptions)
    export ML_GENERAL_OPT="$2";
    shift 2
    ;;
    --trainingdata)
    export ML_TRAINING_FILE="`realpath $2`";
    shift 2
    ;;
    --bkg)
    export ML_BKG_FILE="`realpath $2`";
    shift 2
    ;;
    --testdata)
    export ML_TESTING_FILE="`realpath $2`";
    shift 2
    ;;
    --testbkg)
    export ML_TESTING_BKG="`realpath $2`";
    shift 2
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
    --verbose)
    export VERBOSEMATCHING="1";
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


if [ -z ${GENERATEMCH+x} ] && [ -z ${GENERATEMFT+x} ] && [ -z ${MATCHING+x} ] && [ -z ${CHECKS+x} ] && [ -z ${ML_EXPORTTRAINDATA+x} ] && [ -z ${TRAIN_ML_METHOD+x} ]
then
  echo "Missing use mode!"
  echo " "
  Usage
fi


if [ -z ${OUTDIR+x} ]; then echo "Missing output dir" ; Usage ; fi
NEV_=${NEV_:-"4"}
JOBS="1" # ${JOBS:-"1"} # Forcing O2 simulation with one worker: necessary to keep event ordering
GENERATOR=${GENERATOR:-"gun0_100GeV"}
CUSTOM_SHM="--shm-segment-size 5000000000"

export MCHGENERATOR=${GENERATOR}
export ALIROOT_OCDB_ROOT=${ALIROOT_OCDB_ROOT:-$HOME/alice/OCDB}

ALIROOTENV=${ALIROOTENV:-"AliRoot/latest-master-next-root6"}
O2ENV=${O2ENV:-"O2/latest-dev-o2"}
#O2ENV=${O2ENV:-"O2/latest-f754608ed4-o2"}

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
if ! [ -z ${ML_EXPORTTRAINDATA+x} ]; then exportMLTrainningData ; fi
if ! [ -z ${TRAIN_ML_METHOD+x} ]; then trainML ; fi
if ! [ -z ${CHECKS+x} ]; then runChecks ; fi
