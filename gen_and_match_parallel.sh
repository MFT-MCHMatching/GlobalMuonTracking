#!/bin/bash
MATCHERCRIPT=$HOME/tools/MCHMFTMatching/matcher.sh

NEV=5000
MATCHFCNLIST="matchXYPhiTanl"
CUTFCNLIST="cutDisabled cutDistance cutDistanceAndAngles3Sigma"
NPIONLIST="2 5 10 15 20 50 100 150 200 400" # 800 1200"
NMUONS=2
OUTPUTDIR=MUGun_PiBackground_{1}Pi_${NMUONS}Mu_${NEV}Ev
NJOBS=24

generateMCH () {
parallel --delay 10 ${MATCHERCRIPT} --genMCH -g PiParam --npions {1} \
                          -g MuBoxGun --nmuons ${NMUONS} \
		          					  -n ${NEV} -o $OUTPUTDIR ::: ${NPIONLIST}
}

generateMFT () {
# Using -j 1 here because o2-sim runs simulation on parallel
parallel -j 1 -u ${MATCHERCRIPT} --genMFT -j ${NJOBS} -o $OUTPUTDIR ::: ${NPIONLIST}

}

matchandcheck () {
# Using -j 1 here because match and check need exclusive access to outputdir.
parallel -j 1 -u ${MATCHERCRIPT}  --match --matchFcn {2} --cutFcn {3} --check --updatecode \
                  -o $OUTPUTDIR ::: ${NPIONLIST} ::: ${MATCHFCNLIST} ::: ${CUTFCNLIST}
}

generateMCH
generateMFT
matchandcheck
