#!/bin/bash
#PBS -N {run}_esmeralda_{filenum}
#PBS -q short
#PBS -o /dev/null
#PBS -e /dev/null

echo date
date
source /software/miniconda3/etc/profile.d/conda.sh
export ICTDIR=/software/IC
export ICDIR=$ICTDIR/invisible_cities
export PATH="$ICTDIR/bin:$PATH"
export PYTHONPATH=$ICTDIR:$PYTHONPATH
source activate IC-3.7-2020-06-16
