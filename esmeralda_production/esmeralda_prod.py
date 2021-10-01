from __future__ import print_function

import argparse
from time import sleep
import magic
from subprocess import call
import sys, os
from datetime import datetime
from os.path import dirname
from glob import glob

MC = 0

# ------------------------------------------------------------------------------
# esmeralda_prod.py
# Modified from irene_prod.py
#

tag='uncorrected'
#tag='20181011-20-g01f4a61'
#tag='20181011-19-g25d838b'
trigger_type = 2

def get_parser(args=None):
    parser = argparse.ArgumentParser(description='Script to produce HDF5 files')
    parser.add_argument('-j','--jobs',
                        action='store',
                        help='jobs',
                        required='True')
    parser.add_argument('-r','--run',
                        action='store',
                        help='run number',
                        required='True')
    parser.add_argument('-t','--type',
                        action='store',
                        help='run type',
                        required='True')
    return parser


#get options
args = get_parser().parse_args()
opts = vars(args) # dict
#print(args)

jobs = int(args.jobs)
run = args.run
runtype = args.type

#----------- check and make dirs
def checkmakedir( path ):
    if os.path.isdir( path ):
        print('hey, directory already exists!:\n' + path)
    else:
        os.makedirs( path )
        print('creating directory...\n' + path)

# IO dirs
if(MC == 1):
    PATHIN = '/analysis/MC/{}/hdf5/hdst/'.format(run)
    PATHOUT = '/analysis/MC/{}/hdf5/cdst/'.format(run)
    CONFIGSDIR = '/analysis/MC/{}/hdf5/configs/'.format(run)
    JOBSDIR = '/analysis/MC/{}/hdf5/jobs'.format(run)
else:
    PATHIN     = f'/analysis/{run}/hdf5/prod/v1.2.0-642-g282c69de/{tag}/hdst/trigger{trigger_type}/'
    PATHOUT    = f'/analysis/{run}/hdf5/prod/v1.2.0-642-g282c69de/{tag}/cdst/trigger{trigger_type}/'
    PATHLOG    = f'/analysis/{run}/hdf5/prod/v1.2.0-642-g282c69de/{tag}/logs/cdst/trigger{trigger_type}/'
    CONFIGSDIR = f'/analysis/{run}/hdf5/prod/v1.2.0-642-g282c69de/{tag}/configs/'
    JOBSDIR    = f'/analysis/{run}/hdf5/prod/v1.2.0-642-g282c69de/{tag}/jobs/'

checkmakedir(PATHIN)
checkmakedir(PATHOUT)
checkmakedir(PATHLOG)
checkmakedir(CONFIGSDIR)
checkmakedir(JOBSDIR)

#input files
files = glob(PATHIN + '*.h5')
if(MC):
    files = sorted(files, key=lambda s: int(s.split('.')[-3].split('_')[-1]))
else:
    files = sorted(files, key=lambda s: int(s.split('_')[-6]))  #6
    #files = sorted(files, key=lambda s: int(s.split('.')[-2].split('_')[-1]))

#open template
#TODO read from params

templates = {
'th'    : '/home/shifter/esmeralda_production/templates/esmeralda_th.conf',
'thdetsim'  : '/home/shifter/esmeralda_production/templates/esmeralda_detsim_th.conf'
}

if runtype in templates:
    template_file = templates[runtype]
else:
    template_file = runtype

template = open(template_file).read()
params = {'run': run}

exec_template_file = '/home/shifter/esmeralda_production/templates/esmeralda.sh'
exec_template = open(exec_template_file).read()
exec_params = {'jobsdir' : JOBSDIR,
               'run' : run,
               'filenum' : 0}

# build file list with outputs and filter it
to_process = []
#check hdf5 file are completely written
ftype = magic.Magic()
for f in files:
    #if file complete type would be: "Hierarchical Data Format (version 5) data"
    if ftype.from_file(f) == 'data':
        continue
    filename = f.split('/')[-1]
    filename_split = filename.split('_')
    #filename_split = filename.split('.')
    filename_split[0] = 'cdst'
    filename_out = '_'.join(filename_split)
    #filename_out = "{}.h5".format(filename_out)

    #if file already exists, skip
    fout = PATHOUT + '/' + filename_out
    if os.path.isfile(fout):
        print("skip ", fout)
        continue
    params['filein'] = PATHIN + filename
    params['fileout'] = PATHOUT + filename_out
    config_file = CONFIGSDIR + filename_out + '.conf'
    print(config_file)
    open(config_file, 'w').write(template.format(**params))
    to_process.append(config_file)

#remove old jobs
jobs_files = glob(JOBSDIR + 'esmeralda_*.sh')
map(os.remove, jobs_files)

jobfilename = JOBSDIR + 'esmeralda_0_trigger{}.sh'.format(trigger_type)
jobfile = open(jobfilename, 'w')
jobfile.write(exec_template.format(**exec_params))

nfiles = int(len(to_process) / jobs) + 1
print("Files = {}".format(nfiles))
count_jobs = 0
for i, config in enumerate(to_process, start=1):
    cmd = f'city esmeralda {config} 1>>{PATHLOG}esmeralda_{run}_{count_jobs}_trigger{trigger_type}.out 2>>{PATHLOG}esmeralda_{run}_{count_jobs}_trigger{trigger_type}.err\n'

    if i % nfiles == 0 or i==len(to_process):
        jobfile.write(cmd)
        jobfile.write('\n\necho date\ndate\n')
        jobfile.close()
        
        # Start the next file if we still have more to process.
        if(i < len(to_process)):
            count_jobs += 1
            jobfilename = JOBSDIR + 'esmeralda_{}_trigger{}.sh'.format(count_jobs,trigger_type)
            jobfile = open(jobfilename, 'w')
            exec_params['filenum'] = count_jobs
            jobfile.write(exec_template.format(**exec_params))

    else:
        jobfile.write(cmd)
#jobfile.write('\n\necho date\ndate\n')
#jobfile.close()


#send jobs
for i in range(0, count_jobs+1):
    cmd = 'qsub {}esmeralda_{}_trigger{}.sh'.format(JOBSDIR, i, trigger_type)
    print(cmd)
    #call(cmd, shell=True, executable='/bin/bash')
    os.system(cmd)
    sleep(0.5)

sys.exit()
