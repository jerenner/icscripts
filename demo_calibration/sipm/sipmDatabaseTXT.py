import numpy    as np
import tables   as tb
import argparse

def num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

parser   = argparse.ArgumentParser(description='Upload the DB. Calculation of the mean value of the parameters for different runs.')
parser.add_argument('integer', metavar='Minimum run', type=int)
parser.add_argument(  'files', metavar='Name of the files', nargs='+') #positional argument
ns       = parser.parse_args()

min_run  = ns.integer
files    = ns.files
run_nos  = [f[f.find('R')+1:f.find('R')+5] for f in files]

sensor_number  = []
gain_list      = [[] for i in range(len(files))]
gain_err_list  = [[] for i in range(len(files))]
sigma_list     = [[] for i in range(len(files))]
sigma_err_list = [[] for i in range(len(files))]

for k, file in enumerate(files):
    with open(file, 'r') as df:
        lines  = df.read().split('\n')
        for line in lines:
            if line == '':
                continue
            _, _, sens_no, gain, gain_err, sigma, sigma_err, _, _, _ = [num(x) for x in line.split(',')]
            if k == 0:
                sensor_number.append(sens_no)
            gain_list     [k].append(     gain)
            gain_err_list [k].append( gain_err)
            sigma_list    [k].append(    sigma)
            sigma_err_list[k].append(sigma_err)

gain      = np.mean       (np.array([j for j in      gain_list]), axis=0)
gain_err  = np.linalg.norm(np.array([j for j in  gain_err_list]), axis=0)
sigma     = np.mean       (np.array([j for j in     sigma_list]), axis=0)
sigma_err = np.linalg.norm(np.array([j for j in sigma_err_list]), axis=0)
#last_run  = max(run_nos)


#param_names = ['MinRun', 'MaxRun', 'SensorID', 'Centroid', 'ErrorCentroid', 'Sigma', 'ErrorSigma']
with open('sipmDBvalues_R'+str(min_run)+'.txt', 'w') as out_file:
    for n,sens in enumerate(sensor_number):
        out_file.write(str(min_run)+',100000,'+str(sens)        +','
                        +str(gain[n]) +','+str(gain_err[n]) +','
                        +str(sigma[n])+','+str(sigma_err[n])+'\n')

