from data_pre import *
from con_reg_seq import *
import numpy as np
from multiprocessing import Pool
import multiprocessing
import sys
import os
import scipy.stats as st
from random import *
import scipy.signal


#----- use GPU for TensorFlow -----#
os.environ['TF_CPP_MIN_LOG_LEVEL']='2'

#--------------- step 1: generate chunk seq ---------------#
def generate_chunk_seq(sequence, chunk_size=100, stride=100):
    num_full_chunk = int((float(len(sequence)-4)/stride))
    sequence_list = list()
    # deal with extremely short sequence
    if len(sequence)==0:
        sys.exit('The sequence length is ZERO!!')
    if num_full_chunk==0:
        single_chunk = np.tile(sequence, (chunk_size+4, 1))
        sequence_list.append(single_chunk[:chunk_size+4])
        return np.array(sequence_list)

    for i in range(num_full_chunk):
        sequence_list.append(sequence[i*stride:i*stride+chunk_size+4])
    # deal with the small chunk, using mirrioning to get a full chunk
    if(len(sequence)!=(num_full_chunk*stride+4)):
        remainder = len(sequence)-(num_full_chunk*stride+4)
        add_chunk= sequence[-(remainder+2):]
        comp_chunk = np.flip(sequence[-(chunk_size-remainder+2):],0)
        add_chunk = np.vstack((add_chunk,comp_chunk))
        sequence_list.append(add_chunk)
    return np.array(sequence_list)

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

def convert_to_input(seq_list):
    p = Pool()
    encoding_list = map(sequence_encoding, seq_list)
    seq_chunk_list = p.map(generate_chunk_seq, encoding_list)
    p.close()
    p.join()
    return seq_chunk_list


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
    return out,ali
    
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


#---------- step 3: add Gaussian noise ----------#
def add_noise(std, l):
    noise = np.random.normal(0, std, l)
    return noise




#------ low pass filter for signal simulation
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



#----------- main program: sequence to raw signal --------------#
# default parameters: 
#     repeat_alpha=0.1
#     filter_freq=850
#     noise_std=1.5
def raw_to_true_signal(result_pred, sequence, repeat_alpha, filter_freq, noise_std, perfect, p_len):
    result_pred = np.array(result_pred)
    result_pred = result_pred.flatten()
    final_result = result_pred[:len(sequence)]   #-> this is Z-score
    final_result = np.array(final_result*12.868652 + 90.208199)
    #--- add gauss noise ----#
    if perfect:
        final_result, final_ali = repeat_k_time(p_len, final_result)
    else:
        #-> 1. repeat N times 
        final_result, final_ali = repeat_n_time(repeat_alpha, final_result)
        #-> 2. low pass filter
        if filter_freq>0:
            h,h_start,N = low_pass_filter(4000.0, filter_freq, 40.0)
            final_result = np.convolve(final_result,h)[h_start+1:-(N-h_start-1)+1]
        #-> 3. add gauss noise
        if noise_std>0:
            final_result = final_result + add_noise(noise_std, len(final_result))
    #--- make integer -------#
    final_result = np.array(final_result)
    final_result = np.array(map(int, 5.7*final_result+14))
    return final_result, final_ali

