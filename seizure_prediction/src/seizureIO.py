from scipy.io import loadmat
import pandas as pd
import csv
import os

ictal_class_dict = {'interictal': 0, 'preictal': 1}
def load_sample(patient, segment, ictal_class):
    """
    Load sample from training set 1
    
    df = load_sample(patient=1, segment=1, ictal_class='interictal')
    """
    fname = '../input/train_1/{0}_{1}_{2}.mat'.format(
        patient, segment, ictal_class_dict[ictal_class])
    return load_file(fname)

def load_file(fname):
    """
    Load sample from file
    """
    mat = loadmat(fname)
    mdata = mat['dataStruct']
    mtype = mdata.dtype
    ndata = {n: mdata[n][0,0] for n in mtype.names}
    data_headline = ndata['channelIndices'][0]
    data_raw = ndata['data']
    pdata = pd.DataFrame(data_raw,columns=data_headline)
    iEEGsamplingRate = ndata['iEEGsamplingRate']
    nSamplesSegment = ndata['nSamplesSegment']
    return pdata, iEEGsamplingRate, nSamplesSegment

def get_patient_seq_class(fname):
    fname = os.path.basename(fname)
    I, J, K = fname.split('.')[0].split('_')
    return I, J, K

preictal, r, n = load_sample(patient=1, segment=1, ictal_class='preictal')

outputfname = '../processed/train_1_f_std.csv'

feature_header = ["ch{0}_std".format(i) for i in range(1, 17)]
header = ["patient", "sequence", "class"] + feature_header

with open(outputfname, 'w') as f:
    csv_writer = csv.writer(f)
    csv_writer.writerow(header)

inputdir = '../input/train_1'
for fname in os.listdir(inputdir):
    I, J, K = get_patient_seq_class(fname)
    try:
        d, r, n = load_file(os.path.join(inputdir, fname))
        df = pd.DataFrame(d.std()).transpose()
        df.columns = feature_header
        df["patient"] = I
        df["sequence"] = J
        df["class"] = K
        df = df[header]
        df.to_csv(outputfname, header=False, mode="a")
    except Exception as exc:
        print((fname, exc))


if False:
    freqs = fftfreq(n, 1/r)

    pre1 = preictal[1.0]
    ft_pre1_values = fft(pre1)

    ft_pre1 = pd.Series(data=ft_pre1_values, index=freqs[0])

    preabs = np.abs(ft_pre1)
    presum = np.sqrt(np.sum(ft_pre1**2))
    prerel = preabs / presum
