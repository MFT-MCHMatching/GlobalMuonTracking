#!/bin/bash

piScanAnalysis ()
{
    NEVTS=${NEVTS:-"5000"};
    NMUONS=${NMUONS:-"2"};
    CUTFNC=${CUTFNC:-"cutDistanceAndAngles3Sigma"};
    MATCHFCN=${MATCHFCN:-"matchXYPhiTanl"};
    FILELIST="list_of_checks_${NEVTS}Evts_${NMUONS}Mu_${MATCHFCN}_${CUTFNC}.txt";
    find | grep ${CUTFNC} | grep ${MATCHFCN}_ | grep ${NMUONS}Mu | grep Checks.root | grep ${NEVTS}E | parallel ls {} > ${FILELIST};
    cat ${FILELIST};
    root.exe "PionMultiplicityScanAnalysis.C(\"${FILELIST}\")"
}


piScanAnalysis
