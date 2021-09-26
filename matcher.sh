#!/bin/bash

MATCHINGRESULTS="globalfwdtracks.root MatchingConfig.txt matching.log MatchingPlane_eV*.png ML_Evaluation*.root"
CHECKRESULTS="GlobalMuonChecks.root checks.log images"

Usage()
{
  cat <<-END
  ${0##*/}: a tool to study MCH-MFT Track matching

  Aliroot and O2 environments are automatically loaded.

  Usage:

  1) Generate MCH and MFT Tracks:
     Will create outputdir, include a copy of the macros and generate MCH and MFT tracks with O2 simulation

  ${0##*/} --genMFTMCH -n <number_of_events> -o <outputdir> -g <background generator> -s <signal generator>

    -g  - O2 generator (pythia8pp, pythia8hi, ...). The output of this generator is used as 
             background events in case -s is called
    
    -s  - Generate a signal to embed over the generator configured by -g. Available signal generators:
         BoxMuons : A generator for forward muons with box distribution


    Examples:
    ${0##*/} --genMFTMCH -g pythia8hi -s BoxMuons --nmuons 10 -n 20 -o sampletest

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
  cp -r ${SCRIPTDIR}/O2Generators ${SCRIPTDIR}/*.bin ${SCRIPTDIR}/*.xml ${SCRIPTDIR}/include ${SCRIPTDIR}/*.C ${SCRIPTDIR}/*.h ${SCRIPTDIR}/*.cxx ${SCRIPTDIR}/macrohelpers ${OUTDIR}
}



embedBoxGenMuons () {
  NMUONS=${NMUONS:-"1"}
    # Generate muons with flat P distribution
    o2-sim -n ${NEV_}  $JOBS -g external -m ${MODULES} -o sgn \
       --configKeyValues "GeneratorExternal.fileName=O2Generators/MatcherGenerators.C;GeneratorExternal.funcName=matcherMuBoxGen(${NMUONS})" \
       --embedIntoFile o2sim_Kine.root > embedBoxGenMuons.log 2>&1
}

generateBackGroundEventsO2() {
  echo "Generating background events..."

  o2-sim -g $GENERATOR -m ${MODULES} -n ${NEV_} $JOBS | tee O2Sim.log
}

generateSignalEventsO2() {
  echo "Generating signal events..."

  if ! [ -z ${SIGNALGEN+x} ]; then
     case $SIGNALGEN in
       BoxMuons)
        embedBoxGenMuons;
        shift 1
       ;;
       *) echo "Wrong '-s' signal generator!"; echo "Run '${0##*/} --help'"; exit;
  esac
   SIMS=${SIMS},sgn
  fi
  }


loadO2 () {
echo "Loading O2 environment"
eval `/usr/local/bin/alienv load ${O2ENV} `
}

generateMFTMCHTracksO2()
{
  MODULES="PIPE ITS MFT MCH MID ABSO SHIL DIPO"
  SIMS=o2sim
  mkdir -p ${OUTDIR}
  updatecode
  if ! [ -z ${UPDATECODE+x} ]; then updatecode ; fi
  pushd ${OUTDIR}

  echo "Generating MFT and MFT Tracks with O2: `pwd` ..."

  # Load O2 environment if needed
  type o2-sim &> /dev/null ||  loadO2;


  time generateBackGroundEventsO2 ;
  time generateSignalEventsO2 ;

  echo "Running digitizer ..."
  time o2-sim-digitizer-workflow ${CUSTOM_SHM} -n ${NEV_} --sims ${SIMS} -b --onlyDet MFT,MCH,MID >  O2Digitizer.log
  
  echo "Running MFT reco workflow ..."

  time o2-mft-reco-workflow ${CUSTOM_SHM} -b > O2MFTReco.log
  
  echo "Running MCH reco workflow ..."
  time o2-mch-reco-workflow ${CUSTOM_SHM} -b > O2MCHReco.log
  
  popd
  echo " Finished MFT and MCH Track generation on `realpath ${OUTDIR}`"

}

runO2Match()
{
    # Load O2 environment if needed
  type o2-sim &> /dev/null ||  loadO2;
  echo "_default_o2-globalfwd-matcher-workflow" > MatchingConfig.txt
  echo "Running o2-globalfwd-matcher-workflow (O2) ... "
  o2-globalfwd-matcher-workflow ${CUSTOM_SHM} -b > matching.log
}

runMatching()
{

  if ! [ -f "${OUTDIR}/mchtracks.root" ]; then
    echo " Nothing to Match... MCH Tracks not found on `realpath ${OUTDIR}/mchtracks.root` ..."
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

    time runO2Match ;
    
    #alienv setenv ${O2ENV} -c root.exe -e 'gSystem->Load("libO2MCHTracking")' -l -q -b runMatching.C+ | tee matching.log
    RESULTSDIR="Results`cat MatchingConfig.txt`"
    mkdir -p ${RESULTSDIR}
    cp ${MATCHINGRESULTS} "${RESULTSDIR}" &> /dev/null

    popd
    echo " Finished matching on `realpath ${OUTDIR}`. Results copied to ${RESULTSDIR}"

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

  if ! [ -f "${OUTDIR}/globalfwdtracks.root" ]; then
    echo " Nothing to check... Global Muon Tracks not found on `realpath ${OUTDIR}/globalfwdtracks.root` ..."
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
  type root.exe &> /dev/null ||  loadO2;

  root.exe -l -q -b GlobalMuonChecks.C | tee checks.log
  RESULTSDIR="Results`cat MatchingConfig.txt`"
  mv ${MATCHINGRESULTS} ${CHECKRESULTS} "${RESULTSDIR}"  &> /dev/null
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
    JOBS="-j $2";
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
    -s)
    SIGNALGEN="$2";
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
    --genMFTMCH)
    GENERATEMFTMCH="1";
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


if [ -z ${GENERATEMCH+x} ] && [ -z ${GENERATEMFT+x} ] && [ -z ${GENERATEMFTMCH+x} ] && [ -z ${MATCHING+x} ] && [ -z ${CHECKS+x} ] && [ -z ${ML_EXPORTTRAINDATA+x} ] && [ -z ${TRAIN_ML_METHOD+x} ]
then
  echo "Missing use mode!"
  echo " "
  Usage
fi

if [ -z ${OUTDIR+x} ]; then echo "Missing output dir" ; Usage ; fi
NEV_=${NEV_:-"4"}
GENERATOR=${GENERATOR:-"pythia8hi"}
CUSTOM_SHM="--shm-segment-size 5000000000"

O2ENV=${O2ENV:-"O2/latest-dev-o2"}

if ! [ -z ${GENERATEMFTMCH+x} ]; then
  time generateMFTMCHTracksO2 ;
fi

if ! [ -z ${MATCHING+x} ]; then runMatching ; runChecks ;fi
if ! [ -z ${ML_EXPORTTRAINDATA+x} ]; then exportMLTrainningData ; fi
if ! [ -z ${TRAIN_ML_METHOD+x} ]; then trainML ; fi
if ! [ -z ${CHECKS+x} ]; then runChecks ; fi
