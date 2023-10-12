#!/bin/bash

#outdir=/store/cmst3/group/hgcal/Geant4/
outdir=/store/group/dpg_hgcal/comm_hgcal/yumiao
#outdir=/store/user/yumiao/hgcal/tb/Geant4/
tag=2023Sep
for v in 120; do # 130 131 132 120 121 122 123 124 125 126 127 129; do
    #any eta>5 will generate a 0 angle incidence simulation
    python submitProd.py -t ${tag} -v ${v} -m 0 --etas 50 -d e- -n 5000 -e ${outdir} -g --enList 150 --nRuns 5 --wcuresol 0
    #python submitDigi.py -t ${tag} -v ${v} -m 0 --etas 50 -d e- -n 0 -e ${outdir} -b 0 -g --enList 150 --nRuns 5

    #python submitProd.py -t ${tag} -v ${v} -m 0 --etas 50 -d e+ -n 5000 -e ${outdir} -g --enList 80 100 --nRuns 5
    #python submitDigi.py -t ${tag} -v ${v} -m 0 --etas 50 -d e+ -n 0 -e ${outdir} -b 0 -g --enList 80 100 --nRuns 5
done
