#!/usr/bin/env bash
ARGS=`getopt -o "" -l ",energy:,eta:,run:,granularity:" -n "getopts_${0}" -- "$@"`
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
ETAX10=$(echo ${2} | sed 's/\.//');
echo "etax10: ${ETAX10}";
ETA="${2}";
echo "eta: ${ETA}";
fi
shift 2;;
--run)
if [ -n "$2" ]; then
Step="${2}";
echo "run: ${Step}";
fi
shift 2;;
--granularity)
if [ -n "$2" ]; then
GRAN="${2}";
echo "granularity: ${GRAN}";
fi
shift 2;;
--)
shift
break;;
esac
done

localdir=`pwd`
echo "Job local dir: ${localdir}"
MACFILE="/afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/git2023Sep/version_120/model_0/e-/BOFF/200u/g4steer_en${ENERGY}_eta${ETAX10}_run${Step}.mac"
export HOME=/afs/cern.ch/user/y/yumiao
cd /afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/
source /afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/g4env.sh
cd $localdir
if [ "${GRAN}" -eq 0 ]; then
PFCalEE "${MACFILE}" --model 0 --version 120 --eta ${ETA} --shape 1 --wcuseed 42 --wcuresol 0.0 --fineGranularity
elif [ "${GRAN}" -eq -1 ]; then
PFCalEE "${MACFILE}" --model 0 --version 120 --eta ${ETA} --shape 1 --wcuseed 42 --wcuresol 0.0 --ultraFineGranularity
else
PFCalEE "${MACFILE}" --model 0 --version 120 --eta ${ETA} --shape 1 --wcuseed 42 --wcuresol 0.0
fi
mv PFcal.root HGcal_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.root
echo "--Local directory is $localdir" >> /afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/git2023Sep/version_120/model_0/e-/BOFF/200u/g4_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.log
echo home=$HOME >> /afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/git2023Sep/version_120/model_0/e-/BOFF/200u/g4_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.log
echo path=$PATH >> /afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/git2023Sep/version_120/model_0/e-/BOFF/200u/g4_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.log
echo ldlibpath=$LD_LIBRARY_PATH >> /afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/git2023Sep/version_120/model_0/e-/BOFF/200u/g4_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.log
ls -ltrh * >> /afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/git2023Sep/version_120/model_0/e-/BOFF/200u/g4_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.log
eos mkdir -p /eos/cms/store/group/dpg_hgcal/comm_hgcal/yumiao/git2023Sep/e-
eos cp HGcal_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.root /eos/cms/store/group/dpg_hgcal/comm_hgcal/yumiao/git2023Sep/e-/HGcal_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.root
if (( "$?" != "0" )); then
echo " --- Problem with copy of file PFcal.root to EOS. Keeping locally." >> /afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/git2023Sep/version_120/model_0/e-/BOFF/200u/g4_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.log
else
eossize=`eos ls -l /eos/cms/store/group/dpg_hgcal/comm_hgcal/yumiao/git2023Sep/e-/HGcal_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.root | awk '{print $5}'`
localsize=`ls -l HGcal_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.root | awk '{print $5}'`
if [ "${eossize}" != "${localsize}" ]; then
echo " --- Copy of sim file to eos failed. Localsize = ${localsize}, eossize = ${eossize}. Keeping locally..." >> /afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/git2023Sep/version_120/model_0/e-/BOFF/200u/g4_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.log
else
echo " --- Size check done: Localsize = ${localsize}, eossize = ${eossize}" >> /afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/git2023Sep/version_120/model_0/e-/BOFF/200u/g4_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.log
echo " --- File PFcal.root successfully copied to EOS: /eos/cms/store/group/dpg_hgcal/comm_hgcal/yumiao/git2023Sep/e-/HGcal_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.root" >> /afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/git2023Sep/version_120/model_0/e-/BOFF/200u/g4_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.log
rm HGcal_version120_model0_BOFF_en${ENERGY}_eta${ETAX10}_run${Step}.root
fi
fi
echo "--deleting core files and hepmc files: too heavy!!"
rm core.*
cp * /afs/cern.ch/user/y/yumiao/public/HGCAL_Raw_Data_Handling/CommonModeTB_submission_version/PFCal/PFCalEE/git2023Sep/version_120/model_0/e-/BOFF/200u/
echo "All done"
