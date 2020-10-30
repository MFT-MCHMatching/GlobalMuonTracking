#!/bin/bash
MATCHERCRIPT=$HOME/tools/MCHMFTMatching/matcher.sh

NEV=500
MATCHFCNLIST="matchXY matchXYPhiTanl matchALL"
CUTFCNLIST="cutDisabled cutDistance cutDistanceAndAngles3Sigma"
#MATCHFCNLIST="matchALL"
#CUTFCNLIST="cutDistance_"
NPIONLIST="2 10 20 100 200"
NMUONS=2
OUTPUTDIR=MUGun_PiBackground_{1}Pi_${NMUONS}Mu_${NEV}Ev
NJOBS=1 # Must set -j 1 to keep event ordering in O2 kinematics

generateMCH () {
parallel --delay 10 ${MATCHERCRIPT} --genMCH -g PiParam --npions {1} \
                          -g MuBoxGun --nmuons ${NMUONS} \
                                                                  -n ${NEV} -o $OUTPUTDIR ::: ${NPIONLIST}
}

generateMFT () {
parallel --delay 15 ${MATCHERCRIPT} --genMFT -j ${NJOBS} -o $OUTPUTDIR ::: ${NPIONLIST}

}

matchandcheck () {
# Using -j 1 here because match and check need exclusive access to outputdir.
parallel -j 1 -u ${MATCHERCRIPT}  --match --matchFcn {2} --cutFcn {3} --check --updatecode \
                  -o $OUTPUTDIR ::: ${NPIONLIST} ::: ${MATCHFCNLIST} ::: ${CUTFCNLIST}
}

generateMCH
generateMFT
matchandcheck
