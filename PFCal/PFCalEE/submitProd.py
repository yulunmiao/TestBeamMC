#!/usr/bin/env python

import os, sys
import argparse
import math
import numpy as np
from utils import create_dir, SubmitBase

git_tag=os.popen('git describe --tags --abbrev=0').read()

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-t', '--git-tag'     , dest='gittag'     , help='git tag version', default=git_tag)
parser.add_argument(      '--nRuns'       , dest='nRuns'      , type=int,   help='number of run, 0-indexed', default=-1)
parser.add_argument('-v', '--version'     , dest='version'    , type=int,   help='detector version', default=3)
parser.add_argument('-m', '--model'       , dest='model'      , type=int,   help='detector model', default=3)
parser.add_argument('-a', '--etas'        , dest='etas'       , type=float, help='incidence eta', nargs='+')
parser.add_argument('-p', '--phi'         , dest='phi'        , type=float, help='incidence phi angle in pi unit' , default=0.5)
parser.add_argument(      '--shape'       , dest='shape'      , type=int,   help='shape', default=1) # 1 = hexagons, 2=diamonds, 3=triangles, 4=squares
parser.add_argument('-b', '--Bfield'      , dest='Bfield'     , type=float, help='B field value in Tesla'       , default=0)
parser.add_argument('-d', '--datatype'    , dest='datatype'   ,             help='data type or particle to shoot', default='e-')
parser.add_argument('-f', '--datafile'    , dest='datafile'   ,             help='full path to HepMC input file', default='') #data/example_MyPythia.dat
parser.add_argument('-F', '--datafileeos' , dest='datafileeos',             help='EOS path to HepMC input file', default='') #/eos/cms/store/cmst3/group/hgcal/HGCalMinbias/Pythia8/
parser.add_argument('-n', '--nevts'       , dest='nevts'      , type=int,   help='number of events to generate' , default=1000)
parser.add_argument(    '--wcuseed'       , dest='wcuseed'    , type=int,   help='Seed to use when randomizing WCu density' , default=42)
parser.add_argument(    '--wcuresol'      , dest='wcuresol'   , type=float,  help='Relative resolution to use when randomizing WCu density' , default=-1)
parser.add_argument('-o', '--out'         , dest='out'        ,             help='output directory'             , default=os.getcwd() )
parser.add_argument('-e', '--eosOut'      , dest='eos'        ,             help='eos path to save root file to EOS',         default='')
parser.add_argument('-g', '--gun'         , dest='dogun'      ,             help='use particle gun.', action="store_true")
parser.add_argument(      '--enList'      , dest='enList'     , type=int,   help='E_T list to use with gun', nargs='+', default=[5,10,20,30,40,60,80,100,150,200])
parser.add_argument('-S', '--no-submit'   , dest='nosubmit'   ,             help='Do not submit batch job.', action="store_true")
opt, _ = parser.parse_known_args()

print(opt)

###################################################################################################
###################################################################################################
###################################################################################################
class SubmitProd(SubmitBase):
    def __init__(self, **kwargs):
        super(SubmitProd, self).__init__(**kwargs)

        self.etaint_tag  = '$(ETAX10)'
        
        #variables
        self.condor_submit_name = 'condorSubmitProd.sub'
        self.mac_var = '$(MACFILE)'
        self.mac_name = self._unique_name( (('prefix', 'g4steer'),
                                           ('en', self.shellify_tag(self.en_tag)),
                                           ('eta', self.shellify_tag(self.etaint_tag)),
                                           ('run', self.shellify_tag(self.run_tag)),
                                           ('ext', 'mac')) )

        self.tags = (self.en_tag, self.eta_tag, self.run_tag, self.gran_tag)
        self.labels = ('energy', 'eta', 'run', 'granularity')
        
    def gen_uniform_int_random_seeds_(self, low, high, size):
        np.random.seed()
        r = np.random.uniform(low=low, high=high, size=size)
        return [int(x) for x in r]
        
    def write_shell_script_file(self):
        with open('{}/runJob.sh'.format(self.outDir), 'w') as s:
            s.write('#!/usr/bin/env bash\n')

            #input arguments: energy, eta and run
            s.write('ARGS=`getopt -o "" -l ",energy:,eta:,run:,granularity:" -n "getopts_${0}" -- "$@"`\n')
            s.write('eval set -- "$ARGS"\n')
            s.write('while true; do\n')
            s.write('case "$1" in\n')
            for l,t in zip(self.labels, self.tags):
                s.write('--'+l+')\n')
                s.write('if [ -n "$2" ]; then\n')
                if l=='eta':
                    tmp = "$(echo ${2} | sed 's/\.//');"
                    s.write('{}='.format(self.clean_tag(self.etaint_tag))+tmp+'\n')
                    s.write('echo "'+l+'x10: {}";\n'.format(self.shellify_tag(self.etaint_tag)))
                s.write('{}="${{2}}";\n'.format(self.clean_tag(t)))
                s.write('echo "'+l+': {}";\n'.format(self.shellify_tag(t)))
                s.write('fi\n')
                s.write('shift 2;;\n')
            s.write('--)\n')
            s.write('shift\n')
            s.write('break;;\n')
            s.write('esac\n')
            s.write('done\n\n')
            
            s.write('localdir=`pwd`\n')
            s.write('echo "Job local dir: ${localdir}"\n')
            s.write('{}="{}/{}"\n'.format(self.clean_tag(self.mac_var),self.outDir,self.mac_name))
            s.write('export HOME={}\n'.format(os.environ['HOME']))
            s.write('cd {}/\n'.format(os.getcwd()))
            s.write('source {}/g4env.sh\n'.format(os.getcwd()))
            s.write('cd $localdir\n')
            if len(self.p.datafileeos)>0:
                s.write('eos cp {} {}\n'.format( os.path.join(self.p.datafileeos,self.p.datafile),self.p.datafile))

            cmd = ( 'PFCalEE "{}" --model {} --version {} --eta {} --shape {} --wcuseed {} --wcuresol {}'
                    .format(self.shellify_tag(self.mac_var), self.p.model, self.p.version, self.shellify_tag(self.eta_tag), self.p.shape,self.p.wcuseed,self.p.wcuresol) )
            s.write('if [ "${GRAN}" -eq 0 ]; then\n')
            s.write(cmd + ' --fineGranularity\n')
            s.write('elif [ "${GRAN}" -eq -1 ]; then\n')
            s.write(cmd + ' --ultraFineGranularity\n')
            s.write('else\n')
            s.write(cmd + '\n')
            s.write('fi\n')
            

            outTag = 'version{}_model{}_{}'.format(self.p.version, self.p.model, self.bfield)
            outTag += '_en{}_eta{}'.format(self.shellify_tag(self.en_tag),self.shellify_tag(self.etaint_tag)) 
            if self.p.phi != 0.5: outTag += '_phi{n:.{r}f}pi'.format(n=self.p.phi,r=3)
            outTag += '_run{}'.format(self.shellify_tag(self.run_tag))
            logfile = os.path.join(self.outDir, 'g4_'+outTag+'.log')
            s.write('mv PFcal.root HGcal_{}.root\n'.format(outTag))
            s.write('echo "--Local directory is $localdir" >> {}\n'.format(logfile))
            s.write('echo home=$HOME >> {}\n'.format(logfile))
            s.write('echo path=$PATH >> {}\n'.format(logfile))
            s.write('echo ldlibpath=$LD_LIBRARY_PATH >> {}\n'.format(logfile))
            s.write('ls -ltrh * >> {}\n'.format(logfile))
            if len(self.p.eos)>0:
                s.write('eos mkdir -p {}\n'.format(self.eosDirOut))
                s.write('eos cp HGcal_{}.root {}/HGcal_{}.root\n'.format(outTag,self.eosDirOut,outTag))
                s.write('if (( "$?" != "0" )); then\n')
                s.write('echo " --- Problem with copy of file PFcal.root to EOS. Keeping locally." >> {}\n'.format(logfile))
                s.write('else\n')
                s.write('eossize=`eos ls -l {}/HGcal_{}.root | awk \'{{print $5}}\'`\n'.format(self.eosDirOut,outTag))
                s.write('localsize=`ls -l HGcal_{}.root | awk \'{{print $5}}\'`\n'.format(outTag))
                s.write('if [ "${eossize}" != "${localsize}" ]; then\n')
                s.write('echo " --- Copy of sim file to eos failed. Localsize = ${{localsize}}, eossize = ${{eossize}}. Keeping locally..." >> {}\n'.format(logfile))
                s.write('else\n')
                s.write('echo " --- Size check done: Localsize = ${{localsize}}, eossize = ${{eossize}}" >> {}\n'.format(logfile))
                s.write('echo " --- File PFcal.root successfully copied to EOS: {}/HGcal_{}.root" >> {}\n'.format(self.eosDirOut,outTag,logfile))
                s.write('rm HGcal_{}.root\n'.format(outTag))
                s.write('fi\n')
                s.write('fi\n')

            s.write('echo "--deleting core files and hepmc files: too heavy!!"\n')
            s.write('rm core.*\n')
            if len(self.p.datafileeos)>0:
                s.write('rm {}\n'.format(self.p.datafile))
            s.write('cp * {}/\n'.format(self.outDir))
            s.write('echo "All done"\n')

    def write_geant4_files(self):
        """
        Writes all required geant4 input files, one
        for each run (different seed) and energy.
        """
        niters = self.p.nRuns*len(self.p.enList)*len(self.p.etas)
        gen_kwargs = dict(low=0, high=100000, size=niters)
        seeds1 = self.gen_uniform_int_random_seeds_(**gen_kwargs)
        seeds2 = self.gen_uniform_int_random_seeds_(**gen_kwargs)
        
        for run in range(self.p.nRuns):
            for iet,et in enumerate(self.p.enList):
                for ieta,eta in enumerate(self.p.etas):
                    gen_idx = ( ( run * len(self.p.enList) * len(self.p.etas) ) +
                                ( iet * len(self.p.etas) ) +
                                ( ieta ) )
                    assert(gen_idx < niters)
                    
                    t = ( ('prefix', 'g4steer'), ('en', et), ('eta', int(eta*10.)), ('run', run), ('ext', 'mac') )
                    this_mac_name = self._unique_name(t)
                    with open('{}/{}'.format(self.outDir, this_mac_name), 'w') as s:
                        s.write('/control/verbose 0\n')
                        s.write('/control/saveHistory\n')
                        s.write('/run/verbose 0\n')
                        s.write('/event/verbose 0\n')
                        s.write('/tracking/verbose 0\n')
                        s.write('/N03/det/setField {n:.{r}f} T\n'.format(n=self.p.Bfield,r=1))
                        s.write('/N03/det/setModel {}\n'.format(self.p.model))
                        s.write('/random/setSeeds {} {}\n'.format(seeds1[gen_idx], seeds2[gen_idx]) )
                        if self.p.dogun :
                            s.write('/generator/select particleGun\n')
                            s.write('/gun/particle {} \n'.format(self.p.datatype))
                            en = et*math.cosh(eta) if eta<5 else et
                            s.write('/gun/energy {n:.{r}f} GeV\n'.format(n=en, r=6))
                            if self.p.model != 2 and eta<5:
                                alpha = 2*math.atan(math.exp(-1.*eta));
                                s.write('/gun/direction {} {} {}\n'.format(math.cos(math.pi*self.p.phi)*math.sin(alpha),math.sin(math.pi*self.p.phi)*math.sin(alpha),math.cos(alpha)))
                            
                        else :
                            s.write('/generator/select hepmcAscii\n')
                            s.write('/generator/hepmcAscii/open {}\n'.format(self.p.datafile))
                            s.write('/generator/hepmcAscii/verbose 0\n')
                        s.write('/run/beamOn {}\n'.format(self.p.nevts))


    def write_condor_submission_file(self):
        """
        Writes one single condor submission file, which is expanded to multiple
        jobs for different energies, etas and runs.
        """
        with open('{}/{}'.format(self.outDir,self.condor_submit_name), 'w') as s:
            s.write('universe = vanilla\n')
            s.write('Executable = {}/runJob.sh\n'.format(self.outDir))
            s.write( ('Arguments = --energy {} --eta {} --run {} --granularity {}\n'
                      .format(self.en_tag, self.eta_tag, self.run_tag, self.gran_tag)) )
            #s.write('Requirements = (OpSysAndVer =?= "CentOS7")\n')
            s.write('MY.WantOS = "el7"\n')

            t = ( ('prefix', 'prod'), ('en', self.en_tag), ('eta', self.eta_tag),
                  ('run', self.run_tag) )
            out_name = self._unique_name( t + (('ext', 'out'),) )
            err_name = self._unique_name( t + (('ext', 'err'),))
            log_name = self._unique_name( t + (('ext', 'log'),))
            s.write('Output = {}/{}\n'.format(self.outDir,out_name))
            s.write('Error = {}/{}\n'.format(self.outDir,err_name))
            s.write('Log = {}/{}\n'.format(self.outDir,log_name))
            s.write('RequestMemory = 2GB\n')
            s.write('+JobFlavour = "testmatch"\n')            
            s.write('JobBatchName = prod_' + self.p.gittag + '_' + str(self.p.version) + '_' + self.p.datatype + '\n')
            s.write('Queue {nruns} {entag}, {etatag}, {gtag} from (\n'.format( nruns=self.p.nRuns, entag=self.clean_tag(self.en_tag),
                                                                               etatag=self.clean_tag(self.eta_tag),
                                                                               gtag=self.clean_tag(self.gran_tag)))
            for et in self.p.enList:
                for eta in self.p.etas:
                    gran = 1 if eta < 2.4 else 0
                    if (opt.model==3): gran = -1
                    s.write('{}, {}, {}\n'.format(et,str(eta),gran))
            s.write(')')

###################################################################################################
###################################################################################################
###################################################################################################

bval = 'BON' if opt.Bfield>0 else 'BOFF'
lab = '200u'
odir = '{}/git{}/version_{}/model_{}/{}/{}/{}'.format(opt.out,opt.gittag,opt.version,opt.model,opt.datatype,bval,lab)
if opt.phi != 0.5: odir='{out}/phi_{n:.{r}f}pi'.format(out=odir,n=opt.phi,r=3)
eos_partial = opt.eos[1:] if os.path.isabs(opt.eos) else opt.eos
edir = os.path.join('/eos', 'cms', eos_partial, 'git' + opt.gittag, opt.datatype)

subprod = SubmitProd(outDir=odir, eosDirOut=edir, bfield=bval, params=opt)
subprod.write_shell_script_file()
subprod.write_geant4_files()
subprod.write_condor_submission_file()

os.system('chmod u+rwx {}/runJob.sh'.format(odir))
if opt.nosubmit:
    os.system('echo condor_submit {}/{}'.format(odir, subprod.condor_submit_name)) 
else:
    os.system('condor_submit {}/{}'.format(odir, subprod.condor_submit_name))
