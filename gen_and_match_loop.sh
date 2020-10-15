#!/bin/bash
shopt -s expand_aliases
alias matcher.sh=$HOME/tools/MCHMFTMatching/matcher.sh

NEV=5000
CUTFCNLIST="cutDisabled cutDistance"
NPIONLIST="10 20 50 100 150 200 "
NMUONS=2

generate () {
    ## Generate MCH and MFT Tracks
    for NPIONS in ${NPIONLIST}
    do
      matcher.sh --genMCH --genMFT \
							  -g PiParam --npions ${NPIONS} -g MuBoxGun --nmuons ${NMUONS} \
							  -n ${NEV} -j 24 -o MUGun_PiBackground_${NPIONS}Pi_${NMUONS}Mu_${NEV}
    done
}

matchandcheck () {
    ## Match & Check
    for NPIONS in ${NPIONLIST}
    do
	for CUTFCN in ${CUTFCNLIST}
	do
	    matcher.sh  --match --cutFcn ${CUTFCN} --check --updatecode \
                  -o MUGun_PiBackground_${NPIONS}Pi_${NMUONS}Mu_${NEV}
	done
    done
}

generate
matchandcheck
