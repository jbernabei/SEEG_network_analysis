import sys
import glob
import json
import time

from util import *
from util_connectivity import *
from util_virtual_resection import *
#from util_plot import *

if __name__ == '__main__':
    # Get arguments
    try:
        patient_id = sys.argv[1]
    except TypeError:
        print 'Please enter an appropriate PATIENT_ID.'
        raise
    try:
        epoch_length = sys.argv[2] # seconds
    except IndexError:
        epoch_length = 1
    try:
        run_null = int(sys.argv[3])
    except IndexError:
        run_null = 0

    try:
        starting_null_id = int(sys.argv[4])-1
    except IndexError:
        starting_null_id = 0

    # Check if connectivity already exists
    count_conn = 0
    for k,v in data["PATIENTS"][patient_id]['Events']['Ictal'].items():
        try:
            if(v["STATUS"] == 'ALL_DROPOUT'):
                continue
            else:
                count_conn += 1
        except Exception:
            count_conn += 1

    # Check if all connectivity adjacency matrices have been computed
    if(len(glob.glob('%s/%s/*multiband*'%(os.path.expanduser(data['COMP_DIR']),patient_id))) != count_conn):
        # Compute multi band connectivity and store adjacency matricies
        # blah
        # pass
        compute_multiband_connectivity(patient_id, epoch_length)

    # Get synchronizability
        vr_sync(patient_id, data=data)


    if(len(glob.glob('%s/%s/*noderes*'%(os.path.expanduser(data['COMP_DIR']),patient_id))) != count_conn):
        # Compute node-level virtual resection
        nodal_virtual_resection(patient_id, data=data)


    
