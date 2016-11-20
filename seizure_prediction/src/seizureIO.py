from scipy.io import loadmat
import pandas as pd
import numpy as np
import csv
import os, sys

import lastwinning_features

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

def get_patient_seq_class(fname, istest=False):
    fname = os.path.basename(fname)
    if istest:
        prefix, I, J = fname.split('.')[0].split('_')
        K = 'not_available_for_testset'
    else:
        I, J, K = fname.split('.')[0].split('_')
    return I, J, K

def folder2featurefile(inputdir, outputfname, write_header=False, istest=False):
    feature_header = ["feat{0}".format(i) for i in range(0, 1610)]
    if istest:
        header = ["patient", "sequence"] + feature_header
    else:
        header = ["patient", "sequence", "class"] + feature_header

    if write_header:
        with open(outputfname, 'w') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerow([""] + header)

    for fname in os.listdir(inputdir):
        I, J, K = get_patient_seq_class(fname, istest=istest)
        try:
            d, r, n = load_file(os.path.join(inputdir, fname))

            #df = pd.DataFrame(d.std()).transpose()
            #df.columns = feature_header
            
            feature_vals = lastwinning_features.calculate_features(d.values, r, n)
            data = np.array([feature_vals])
            df = pd.DataFrame.from_records(data, columns=feature_header)

            df["patient"] = I
            df["sequence"] = J
            if not(K == 'not_available_for_testset'):
                df["class"] = K
            df = df[header]
            df.to_csv(outputfname, header=False, mode="a")
        except Exception as exc:
            print((fname, exc))

            
if __name__ == "__main__":
    outputfname = sys.argv[1]
    if outputfname.startswith('test'):
        outputfname = os.path.join('..', 'processed', outputfname)
        print('Generate test data: {0}'.format(outputfname))
        
        inputdir = '../input/test_1_new'
        folder2featurefile(inputdir, outputfname, write_header=True, istest=True)

        inputdir = '../input/test_2_new'
        folder2featurefile(inputdir, outputfname, write_header=False, istest=True)

        inputdir = '../input/test_3_new'
        folder2featurefile(inputdir, outputfname, write_header=False, istest=True)
    else:
        outputfname = os.path.join('..', 'processed', outputfname)
        print('Generate training data: {0}'.format(outputfname))

        inputdir = '../input/train_1'
        folder2featurefile(inputdir, outputfname, write_header=True)

        inputdir = '../input/train_2'
        folder2featurefile(inputdir, outputfname, write_header=False)

        inputdir = '../input/train_3'
        folder2featurefile(inputdir, outputfname, write_header=False)

    
    
if False:
    freqs = fftfreq(n, 1/r)

    pre1 = preictal[1.0]
    ft_pre1_values = fft(pre1)

    ft_pre1 = pd.Series(data=ft_pre1_values, index=freqs[0])

    preabs = np.abs(ft_pre1)
    presum = np.sqrt(np.sum(ft_pre1**2))
    prerel = preabs / presum
