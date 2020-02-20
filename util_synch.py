#!/usr/bin/python
'''
Utility module to compute virtual resection.
'''
from util import *

np.random.seed(sum(map(ord, "aesthetics")))

with open('../data/DATA.json') as json_data_file:
    data = json.load(json_data_file)

warnings.filterwarnings('ignore')

def base_synchronizability(adj):
    """
    Function for computing base synchronizability of the entire network
    Parameters
    ----------
        adj: ndarray, shape (N, N)
            Undirected, symmetric adjacency matrix with N nodes

    Returns
    -------
        base_sync: float
            Base synchronizability of the network
    """

    # Standard param checks
    errors.check_type(adj, np.ndarray)
    errors.check_dims(adj, 2)
    if not (adj == adj.T).all():
        raise Exception('Adjacency matrix is not undirected and symmetric.')
    if(np.isnan(adj).any()):
        return np.nan

    # Get data attributes
    n_node = adj.shape[0]

    # Get the original synchronizability
    base_sync = synchronizability(adj)

    return base_sync



def vr_sync(patient_id, dilate_radius=0, data=data):
    """
    Function for computing c_resection(t).
    Parameters
    ----------
        patient_id: str
            Patient ID in DATA.json
        unique_id: str
            Unique UUID code for npz files to load adjacency matrices
        dilate_radius: int
            Radius of dilation (>0) or erosion (<2) of the resection estimate.
        data: dict
            Dictionary of data loaded from DATA.json (or TEST_DATA.json during unit tests). Default is DATA.json.
    Returns
    -------
        None
            Saves all adjacency matrices in different bands as npz files in comp_dir.
    """

    # Generate list of cartoon map labels
    labels = map(lambda x: x.split(',')[1].replace('\n',''), open(os.path.expanduser(
        data['PATIENTS'][patient_id]['ELECTRODE_LABELS']
        ),'r').readlines())

    # Get path
    comp_dir = os.path.expanduser(data['COMP_DIR'])
    data_dir = os.path.expanduser(data['DATA_DIR'])

    # Load ignored node labels
    ignored_node_labels = data['PATIENTS'][patient_id]['IGNORE_ELECTRODES']
    for ignored_node_label in ignored_node_labels:
        if(ignored_node_label not in labels):
            labels.append(ignored_node_label)

    # Create output UUID codx
    unique_idx = []

    # Load ictal clips and get data as T x N for T = epoch_length (seconds) * fs
    for event_type, events in data['PATIENTS'][patient_id]['Events'].items():
        unique_id = str(uuid.uuid4())
        for event_id in events.keys():
            try:
                if(events[event_id]['STATUS'] == 'ALL_DROPOUT'):
                        continue # unusable clip
            except KeyError:
                pass

            fn = os.path.join(data_dir, patient_id, events[event_id]['FILE'])
            channels = []

            # Get channels, ECoG Data, Fsx
            with h5py.File(fn) as f:
                evData = f['evData'].value
                Fs = f['Fs'].value
                for column in f['channels']:
                    row_data = []
                    for row_number in range(len(column)):
                        row_data.append(''.join(map(unichr, f[column[row_number]][:])))
                    channels.append(row_data)
            Fs = int(Fs[0][0])
            channels = channels[0]
            # evData = scipy.stats.zscore(evData,axis=1)
            T = evData.shape[0]

            # Correspond label names
            labels_dict = correspond_label_names(channels, labels)

            # Load electrodes to ignore
            ignored_node_idx  = map(lambda x: labels_dict[x][0], ignored_node_labels)
            for ii,node_id in enumerate(ignored_node_idx):
                print 'Ignoring node label: %s because label %s is in IGNORE_ELECTRODES'%(channels[node_id],ignored_node_labels[ii])
            channels = list(np.delete(np.array(channels),ignored_node_idx))

            # Recorrespond label names
            labels_dict = correspond_label_names(channels, labels)

            # For each clip, load up adjacency matrices
            adj_file = np.load(os.path.join(comp_dir,patient_id,'%s.%s.%s.multiband.npz'%(patient_id,event_type,event_id)))

            epoch_length = int(adj_file['epoch_length'])
            all_adj_alphatheta = adj_file['all_adj_alphatheta']
            all_adj_beta = adj_file['all_adj_beta']
            all_adj_lowgamma = adj_file['all_adj_lowgamma']
            all_adj_highgamma = adj_file['all_adj_highgamma']
            all_adj_broadband_CC = adj_file['all_adj_broadband_CC']
            epochs = int(T/(epoch_length*Fs))

            assert all_adj_alphatheta.shape[2] == epochs
            assert all_adj_beta.shape[2] == epochs
            assert all_adj_lowgamma.shape[2] == epochs
            assert all_adj_highgamma.shape[2] == epochs
            assert all_adj_broadband_CC.shape[2] == epochs

            # Compute base synchronizability of network
            base_sync_alphatheta = np.zeros((epochs,))
            base_sync_beta = np.zeros((epochs,))
            base_sync_lowgamma = np.zeros((epochs,))
            base_sync_highgamma = np.zeros((epochs,))
            base_sync_broadband_CC = np.zeros((epochs,))

            for epoch in range(epochs):
                if(np.isnan(all_adj_alphatheta[:,:,epoch]).any()):
                    base_sync_alphatheta[epoch] = np.nan
                else:
                    base_sync_alphatheta[epoch] = base_synchronizability(all_adj_alphatheta[:,:,epoch])
                if(np.isnan(all_adj_beta[:,:,epoch]).any()):
                    base_sync_beta[epoch] = np.nan
                else:
                    base_sync_beta[epoch] = base_synchronizability(all_adj_beta[:,:,epoch])
                if(np.isnan(all_adj_lowgamma[:,:,epoch]).any()):
                    base_sync_lowgamma[epoch] = np.nan
                else:
                    base_sync_lowgamma[epoch] = base_synchronizability(all_adj_lowgamma[:,:,epoch])
                if(np.isnan(all_adj_highgamma[:,:,epoch]).any()):
                    base_sync_highgamma[epoch] = np.nan
                else:
                    base_sync_highgamma[epoch] = base_synchronizability(all_adj_highgamma[:,:,epoch])
                if(np.isnan(all_adj_broadband_CC[:,:,epoch]).any()):
                    base_sync_broadband_CC[epoch] = np.nan
                else:
                    base_sync_broadband_CC[epoch] = base_synchronizability(all_adj_broadband_CC[:,:,epoch])

            # Save with appropriate name
            print 'Writing sync for patient %s event %s %s'%(patient_id,event_type,event_id)
            sync_fn = os.path.join(comp_dir,patient_id,'%s.%s.%s.sync.%s.npz'%(patient_id,event_type,event_id,unique_id))
            np.savez(open(sync_fn,'w'), base_sync_alphatheta=base_sync_alphatheta, base_sync_beta=base_sync_beta, base_sync_lowgamma=base_sync_lowgamma, base_sync_highgamma=base_sync_highgamma, base_sync_broadband_CC=base_sync_broadband_CC)
            pipeline_fn = os.path.join(comp_dir,patient_id,'%s.%s.%s.sync.%s.pipedef.json'%(patient_id,event_type,event_id,unique_id))
            timestamp = datetime.datetime.now()
            unique_idx.append((unique_id,event_type,event_id))
    return unique_idx

