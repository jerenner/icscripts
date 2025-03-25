import argparse
import numpy    as np
import tables   as tb

from functools import partial
from invisible_cities.io.channel_param_io import subset_param_reader as spr


parser = argparse.ArgumentParser(description='Extract gain values from the param. hdf5 file.')
parser.add_argument('file_in', type=str, help='Name of the file') #positional argument
ns     = parser.parse_args()

data_file   = ns.file_in
run_number  = data_file[data_file.find('R')+1:data_file.find('R')+6] #to change based on number of digit

param_names = ['pedestal_sigma', 'poisson_mu', 'gain', 'gain_sigma', 'n_gaussians_chi2']

read_params = partial(spr, table_name='FIT_pmt_ngau_4',
                      param_names=param_names)

with open('pmtGain_R'+run_number+'.txt', 'w') as out_file:
    with tb.open_file(data_file) as df:
        for sens, (pars, errs) in read_params(df):
            out_file.write(run_number+', 100000, ' +str(sens)+', '
                           +str(pars['pedestal_sigma'  ])+', '+str(errs['pedestal_sigma'])+', '
                           +str(pars['poisson_mu'      ])+', '+str(errs['poisson_mu'    ])+', '
                           +str(pars['gain'            ])+', '+str(errs['gain'          ])+', '
                           +str(pars['gain_sigma'      ])+', '+str(errs['gain_sigma'    ])+', '
                           +str(errs['n_gaussians_chi2'])+'\n')
        
