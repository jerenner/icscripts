# icscripts

## Esmeralda production scripts on frontend1demo

The scripts are in the `esmeralda_production` folder

```
[shifter@frontend1demo ~]$ cd esmeralda_production/
[shifter@frontend1demo esmeralda_production]$ ls
ana-code                          krmapgen      out_9928.dat  out_9969.dat
combine_cdst.py                   maps          out_9929.dat  t_evol.png
esmeralda_0_trigger2.sh.e9701696  notebooks     out_9931.dat  templates
esmeralda_0_trigger2.sh.o9701696  out_9925.dat  out_9932.dat
esmeralda_prod.py                 out_9926.dat  out_9966.dat
icsetup.sh                        out_9927      out_9967
[shifter@frontend1demo esmeralda_production]$ python esmeralda_prod.py
usage: esmeralda_prod.py [-h] -j JOBS -r RUN -t TYPE
esmeralda_prod.py: error: the following arguments are required: -j/--jobs, -r/--run, -t/--type
```
-j is the number of cores to run on (in parallel)<br>
-r is the run number<br>
-t is the type, specifying which parameter file to use (L88 of esmeralda_prod.py shows the possible options)

At the moment the "types" are:
```
templates = {
'th'    : '/home/shifter/esmeralda_production/templates/esmeralda_th.conf',
'thdetsim'  : '/home/shifter/esmeralda_production/templates/esmeralda_detsim_th.conf'
}
```

For data, use type `th`. Parameters for the Esmeralda city can be changed in `templates/esmeralda_th.conf`. For example, to launch run 9859 on 20 cores:

```
[shifter@frontend1demo esmeralda_production]$ python esmeralda_prod.py -j 20 -r 9859 -t th
```

The output CDSTs will appear in `/analysis/9859/hdf5/prod/v1.2.0-642-g282c69de/uncorrected/cdst` for run 9859.

One can check the number of running processes with `qstat`. When all processes are done, one can combine the CDSTs into a single file using `combine_cdst.py`:

```
[shifter@frontend1demo esmeralda_production]$ python combine_cdst.py
usage: combine_cdst.py [-h] -r RUN
combine_cdst.py: error: the following arguments are required: -r/--run
```

As this may take awhile, one can run the job in the background so that it will continue to run, even if the ssh session is closed:

```
python combine_cdst.py -r 9859 >& out_9859.dat &
```
