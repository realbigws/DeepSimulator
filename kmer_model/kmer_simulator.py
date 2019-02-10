import numpy as np
from functools import partial
from multiprocessing import Pool
import multiprocessing
import scipy.stats as st
from random import *
import scipy.signal
import pickle
import argparse
from tqdm import *


#--------------- step 0: official kmer pore model ---------#
#-> load '.pkl' for official pore model
#   data structure is a dict with 4096 dims (6mer 4^6)
#   data.keys()   :   kmer
#   data.values() :   (mean, vari)
def load_official_poremodel(input_pkl):
    f=open(input_pkl, 'rb')
    data=pickle.load(f)
    f.close()
    return data

def sequence_official_poremodel(sequence, kmer_poremodel):
    k=len(kmer_poremodel.keys()[0])
    length=len(sequence)
    # check sequence length
    if length < k:
        # Assign mean and std value
        kmer_means=list()
        kmer_stdvs=list()
        [kmer_means.extend((float(90.2083),)*1) for i in range(length)]
        [kmer_stdvs.extend((float(2.0),)*1) for i in range(length)]
    else:
        # Divide sequence into kmers
        kmers = [sequence[i:i + k] for i in range(0, length - k + 1)]
        # Assign mean and std value
        kmer_means, kmer_stdvs = zip(*[kmer_poremodel[kmer] for kmer in kmers])
        kmer_means=list(kmer_means)
        kmer_stdvs=list(kmer_stdvs)
        # Append tail
        [kmer_means.extend((float(90.2083),)*1) for i in range(k-1)]
        [kmer_stdvs.extend((float(2.0),)*1) for i in range(k-1)]
    # return
    kmer_means = np.array(kmer_means)
    kmer_stdvs = np.array(kmer_stdvs)
    return kmer_means,kmer_stdvs


#--------------- step 1: load input sequence ---------------#
def get_seq_list(file_name):
    with open(file_name, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    seq_list = filter(lambda x: x!='', lines)
    seq_list = filter(lambda x: '>' not in x, seq_list)
    return seq_list

def get_id_list(file_name):
    with open(file_name, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    lines = filter(lambda x: '>' in x, lines)
    id_list = map(lambda x: x.split('|')[0][1:], lines)
    return id_list


#------ output functions ------#
def write_output(result, file_name):
    with open(file_name, 'w') as f:
        for i in result:
            temp = str(i)+'\n'
            f.write(temp)

def write_alignment(result, file_name):
    with open(file_name, 'w') as f:
        for i in result:
            temp = str(i[0]+1)+' '+str(i[1]+1)+'\n'
            f.write(temp)


#---------- step 2: repeat length sample -----------#
def rep_rvs(size,a):
    a = a*5
    array_1 = np.ones(int(size*(0.075-0.015*a))).astype(int)
    samples = st.alpha.rvs(3.3928495261646932+a,
        -7.6451557771999035+(2*a), 50.873948369526737,
        size=(size-int(size*(0.075-0.015*a)))).astype(int)
    samples = np.concatenate((samples, array_1), 0)
    addi = np.array(abs(np.random.normal(2,1,size))).astype(int)
    samples[samples<2] = 2
    samples[samples>40] = 40
    samples[samples<8] += addi[samples<8]
    np.random.shuffle(samples)
    samples[samples<8] += addi[samples<8]
    return samples

def repeat_n_time(a, result):
    rep_times = rep_rvs(len(result), a)
    out = list()
    ali = list()
    pos = 0
    for i in range(len(result)):
        k = rep_times[i]
        cur = [result[i]] * k
        out.extend(cur)
        for j in range(k):
            ali.append((pos,i))
            pos = pos + 1
    event_idx = np.repeat(np.arange(len(result)), rep_times)
    return out,ali,event_idx
    
def repeat_k_time(k, result):
    out = list()
    ali = list()
    pos = 0
    for i in range(len(result)):
        cur = [result[i]] * k
        out.extend(cur)
        for j in range(k):
            ali.append((pos,i))
            pos = pos + 1
    return out,ali


#------------- step 3: low pass filter for signal simulation -----#
#-> low pass filter
#   sampling_rate = 4000.0, cut_off_freq = 1750.0, bandwidth_freq = 40.0
def low_pass_filter(sampling_rate, cut_off_freq, bandwidth_freq):
    # Read input parameter
    fS = sampling_rate  # Sampling rate.
    fL = cut_off_freq   # Cutoff frequency.
    fb = bandwidth_freq # Bandwidth frequency

    # Generate frequency bin
    b = fb / fS
    N = int(np.ceil((4 / b)))
    if not N % 2: N += 1  # Make sure that N is odd.
    n = np.arange(N)

    # Compute sinc filter.
    h = np.sinc(2 * fL / fS * (n - (N - 1) / 2.))

    # Compute Blackman window.
    w = 0.42 - 0.5 * np.cos(2 * np.pi * n / (N - 1)) + \
        0.08 * np.cos(4 * np.pi * n / (N - 1))

    # Compute h and h_start
    h = h * w
    h /= np.sum(h)
    impulse = np.repeat(0., len(h))
    impulse[0] = 1.
    h_response = scipy.signal.lfilter(h, 1, impulse)
    h_start = np.argmax(h_response)

    # return
    return h,h_start,N


#---------- step 4: add Gaussian noise ----------#
def add_noise(std, l):
    noise = np.random.normal(0, std, l)
    return noise



#----------- main program: sequence to raw signal --------------#
# default parameters: 
#     repeat_alpha=0.1
#     event_std=1.0
#     filter_freq=850
#     noise_std=1.5
def sequence_to_true_signal(sequence, kmer_poremodel='null', perfect=0, p_len=1,
    repeat_alpha=0.1, independ=1, event_std=1.0, filter_freq=850, noise_std=1.5 ):
    #--- get kmer signal ----#
    mean_result, std_result = sequence_official_poremodel(sequence, kmer_poremodel)
    #--- add gauss noise ----#
    if perfect:
        final_result, final_ali = repeat_k_time(p_len, mean_result)
    else:
        #-> 1. repeat N times 
        indep_result, final_ali, event_idx = repeat_n_time(repeat_alpha, mean_result)
        event_std = np.random.uniform(-1*event_std*std_result[event_idx], event_std*std_result[event_idx])
        depend_result = mean_result[event_idx] + event_std
        #-> 2. independent kmer 
        if independ:
            final_result = indep_result
        else:
            final_result = depend_result
        #-> 3. low pass filter
        if filter_freq>0:
            h,h_start,N = low_pass_filter(4000.0, filter_freq, 40.0)
            final_result = np.convolve(final_result,h)[h_start+1:-(N-h_start-1)+1]
        #-> 4. add gauss noise
        if noise_std>0:
            final_result = final_result + add_noise(noise_std, len(final_result))
    #--- make integer -------#
    final_result = np.array(final_result)
    final_result = np.array(map(int, 5.7*final_result+14))
    return final_result, final_ali


#=================== main =======================#
if __name__ == '__main__':

    #--------- argument part -------------#
    parser = argparse.ArgumentParser(description='convert the \
        input sequence to nanopore signal')
    parser.add_argument('-i', action='store', dest='input', required=True, 
        help='the input file')
    parser.add_argument('-p', action='store', dest='output', required=True,
        help='prefix the output file')
    parser.add_argument('-l', action='store', dest='alignment', required=True,
        help='prefix the alignment file')
    parser.add_argument('-m', action='store', dest='poremodel', required=True,
        help='official kmer pore model')
    parser.add_argument('-t', action='store', dest='threads',
        type=int, help='the number of threads used. (default is 1)', default=1)
    parser.add_argument('-e', action='store', dest='event_std',
        type=float, help='set the std of the event. \
        The higher the value, the more variable the event. (default is 1.0)',
        default=1.0)
    parser.add_argument('-f', action='store', dest='filter_freq',
        type=float, help='change the cut frequency in the low pass filter. \
        The higher the value, the smoother the signal. (default is 850) \
        Set -1 to disable',
        default=850)
    parser.add_argument('-s', action='store', dest='noise_std',
        type=float, help='set the std of the noise. \
        The higher the value, the blurred the signal. (default is 1.5)',
        default=1.5)
    parser.add_argument('--perfect', action='store', dest='perfect',
        type=bool, help='Do you want a perfect signal and sequence',
        default=False)
    parser.add_argument('--independ', action='store', dest='independ',
        type=bool, help='Signal varies during each event or not',
        default=False)


    #---------- input list ---------------#
    arg = parser.parse_args()
    seq_list = get_seq_list(arg.input)
    id_list = get_id_list(arg.input)

    #---------- load pore model ----------#
    kmer_poremodel=load_official_poremodel(arg.poremodel)

    #---------- partial function ---------#
    func=partial(sequence_to_true_signal, \
    	kmer_poremodel=kmer_poremodel, perfect=arg.perfect, independ=arg.independ, \
    	event_std=arg.event_std, filter_freq=arg.filter_freq, noise_std=arg.noise_std)

    #---------- multi process ------------#
    p = Pool(arg.threads)
    result_list = []
    result_list = list(tqdm(
        p.imap(func, seq_list), 
        total=len(seq_list)))
    p.close()
    p.join()

    #----------- output part -------------#
    for i in range(len(result_list)):
        final_signal, final_ali = result_list[i]
        write_output(final_signal, arg.output+'_{}.txt'.format(id_list[i]))
        write_alignment(final_ali, arg.alignment+'_{}.ali'.format(id_list[i]))


