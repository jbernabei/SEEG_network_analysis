import sys
import glob
import json
import time

from util_synch import *
from util_connectivity_synch import *

if __name__ == '__main__':
    # Get arguments
    try:
        patient_id = sys.argv[1]
    except TypeError:
        print 'Please enter an appropriate PATIENT_ID.'
        raise
    try:
        epoch_length = int(sys.argv[2]) # seconds
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

    # Get multiband connectivity 
        compute_multiband_connectivity(patient_id, epoch_length)

    # Get synchronizability
        vr_sync(patient_id, data=data)




    
