#!/usr/bin/env bash
ARGS=`getopt -o "" -l ",energy:,eta:,run:,ic:,npuvtx:" -n "getopts_${0}" -- "$@"`
eval set -- "$ARGS"
while true; do
case "$1" in
--energy)
if [ -n "$2" ]; then
ENERGY="${2}";
echo "energy: ${ENERGY}";
fi
shift 2;;
--eta)
if [ -n "$2" ]; then
ETA=$(echo ${2} | sed 's/\.//');
echo "eta: ${ETA}";
fi
shift 2;;
--run)
if [ -n "$2" ]; then
Step="${2}";
echo "run: ${Step}";
fi
shift 2;;
--ic)
if [ -n "$2" ]; then
IC="${2}";
echo "ic: ${IC}";
fi
shift 2;;
--npuvtx)
if [ -n "$2" ]; then
NPUVTX="${2}";
echo "npuvtx: ${NPUVTX}";
fi
shift 2;;
--)
shift
break;;
esac
done

localdir=`pwd`
export HOME=/afs/cern.ch/user/y/yumiao
cd /afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE
source g4env.sh
echo $PATH
cd $localdir
/afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/userlib/bin/digitizer -c /afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/userlib/DigiConfig.cfg -n 0 -i /eos/cms/store/cmst3/group/hgcal/Geant4//git2023Sep/e-/HGcal_version120_model0_BOFF_en${ENERGY}_eta${ETA}_run${Step}.root -o $localdir/ --granulStr=0-1:1  --noiseStr=0-1:0.15 --threshStr=0-1:5 --interCalib=${IC} --nSiLayers=3 --nPU=${NPUVTX} --puPath=root://eoscms//eos/cms/store/cmst3/group/hgcal/Standalone/V12/MinBias/ -a false
echo "--Local directory is " $localdir
echo home=$HOME
echo path=$PATH
echo ldlibpath=$LD_LIBRARY_PATH
ls *
eos mkdir -p /eos/cms/store/cmst3/group/hgcal/Geant4//git2023Sep/e-
eos cp $localdir/DigiPFcal.root /eos/cms/store/cmst3/group/hgcal/Geant4//git2023Sep/e-/Digi_300uversion120_model0_BOFF_npuvtx${NPUVTX}_ic${IC}_en${ENERGY}_eta${ETA}_run${Step}.root
if (( "$?" != "0" )); then
echo " --- Problem with copy of file DigiPFcal.root to EOS. Keeping locally."
else
eossize=`eos ls -l /eos/cms/store/cmst3/group/hgcal/Geant4//git2023Sep/e-/Digi_300uversion120_model0_BOFF_npuvtx${NPUVTX}_ic${IC}_en${ENERGY}_eta${ETA}_run${Step}.root | awk '{print $5}'`
localsize=`ls -l DigiPFcal.root | awk '{print $5}'`
if [ ${eossize} != ${localsize} ]; then
echo " --- Copy of digi file to eos failed. Localsize = ${localsize}, eossize = ${eossize}. Keeping locally..."
else
echo " --- Size check done: Localsize = ${localsize}, eossize = ${eossize}"
echo " --- File DigiPFcal.root successfully copied to EOS: /eos/cms/store/cmst3/group/hgcal/Geant4//git2023Sep/e-/Digi_300uversion120_model0_BOFF_npuvtx${NPUVTX}_ic${IC}_en${ENERGY}_eta${ETA}_run${Step}.root"
rm DigiPFcal.root
fi
fi
