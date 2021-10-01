from __future__ import print_function

import numpy as np
import argparse
import sys, os
from glob import glob
import tables

# ------------------------------------------------------------------------------
# combine_cdst.py
# Combine pytables from analyses over many cdst files.
#

MC = False

tag = 'uncorrected'
#tag='20181011-20-g01f4a61'
#tag='20181011-19-g25d838b'
trigger_type = 2
include_tracks = True
include_dst = True

def get_parser(args=None):
    parser = argparse.ArgumentParser(description='Script to produce HDF5 files')
    parser.add_argument('-r','--run',
                        action='store',
                        help='run number',
                        required='True')
    return parser

# get arguments
args = get_parser().parse_args()
opts = vars(args) # dict
run = args.run

if(MC):
    pyt_dir = "/analysis/MC/{}/hdf5/cdst".format(run)
else:
    pyt_dir = f'/analysis/{run}/hdf5/prod/v1.2.0/{tag}/cdst/trigger{trigger_type}/'

# input/output files
files = glob(pyt_dir + '/cdst_*.h5')
files = sorted(files, key=lambda s: int(s.split('_')[-6])) # 6
#files = sorted(files, key=lambda s: int(s.split('_')[-2].split('.')[-2][3:]))
out_file = "{}/cdst_combined_{}.h5".format(pyt_dir,run)

# Open the final tables file.
fcombined = tables.open_file("{0}/cdst_combined_{1}.h5".format(pyt_dir,run), "w", filters=tables.Filters(complib="blosc", complevel=9))
group_Summary = fcombined.create_group(fcombined.root, "Summary")
if(include_tracks):
    group_Tracking = fcombined.create_group(fcombined.root, "Tracking")
if(include_dst):
    group_DST = fcombined.create_group(fcombined.root, "DST")

# Open the first file.
f1 = tables.open_file(files[0], 'r')
summary_combined = f1.copy_node('/Summary', name='Events', newparent=group_Summary)
if(include_tracks):
    tracking_combined = f1.copy_node('/Tracking', name='Tracks', newparent=group_Tracking)
if(include_dst):
    dst_combined = f1.copy_node('/DST', name='Events', newparent=group_DST)
f1.close()

# Process nruns files.
for fname in files[1:]:

    print("-- Adding file {0}".format(fname));

    # Open the next file and extract the elements.
    fn = tables.open_file(fname, 'r')
    if("/Summary" in fn):
        summary_events = fn.root.Summary.Events
        summary_combined.append(summary_events.read())
    if(include_tracks and "/Tracking" in fn):
        tracking_tracks = fn.root.Tracking.Tracks
        tracking_combined.append(tracking_tracks.read())
    if(include_dst and "/DST" in fn):
        dst_events = fn.root.DST.Events
        dst_combined.append(dst_events.read())
    else:
        print(" --> Skipping file due to no table found")

    # Close the file.
    fn.close()

    # Flush the combined file.
    fcombined.flush()

# Close the combined hdf5 file.
print("Saving combined file {0}/cdst_combined_{1}.h5".format(pyt_dir,run))
fcombined.close()
