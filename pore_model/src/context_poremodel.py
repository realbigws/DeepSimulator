from data_pre import *
import math
import numpy as np
from multiprocessing import Pool
import multiprocessing
from con_reg_seq import *
import cPickle
import gc
import os
import argparse
import sys
import tqdm
import scipy.stats as st

os.environ['TF_CPP_MIN_LOG_LEVEL']='2'

# the threshold for creat bins
threshold = 2.5
num_bins = 25 + 2
bin_size = threshold*2/(num_bins-2)

data_type =['sequence', 'fix_value', 'can_value', 'label','adp']


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
    samples[samples<2] = 2
    samples[samples>40] = 40
    np.random.shuffle(samples)
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


#----------- main program: sequence to raw signal --------------#
def raw_to_true_signal(result_pred, sequence, repeat_alpha, noise_std, perfect, p_len):
	result_pred = np.array(result_pred)
	result_pred = result_pred.flatten()
	final_result = result_pred[:len(sequence)]   #-> this is Z-score
	final_result = np.array(map(int, 5.7*(final_result*12.868652 + 90.208199) + 14 ))
	#final_result = np.array(map(int, 3.8*(final_result*13.239392 + 98.225867)))
	if perfect:
		final_result, final_ali = repeat_k_time(p_len, final_result)
	else:
		final_result, final_ali = repeat_n_time(repeat_alpha, final_result)
	if noise_std>0:
		final_result = final_result + add_noise(noise_std, len(final_result))
	return final_result, final_ali


#=================== main =======================#
if __name__ == '__main__':
	os.environ['CUDA_VISIBLE_DEVICES'] ='-1'
	parser = argparse.ArgumentParser(description='convert the \
		input sequence to context-dependent pore model')
	parser.add_argument('-i', action='store', dest='input', required=True, 
		help='the input file')
	parser.add_argument('-o', action='store', dest='output', required=True,
		help='the output file')
	parser.add_argument('-t', action='store', dest='threads',
		type=int, help='the number of threads used', default=1)

	arg = parser.parse_args()
	seq_list = get_seq_list(arg.input)
	id_list = get_id_list(arg.input)
	seq_chunk_list = convert_to_input(seq_list)


	result_list = []
	p = Pool(arg.threads)
	result_list = list(tqdm.tqdm(
		p.imap(model_whole_set_check, seq_chunk_list), 
		total=len(seq_chunk_list)))
	p.close()
	p.join()


	for i in range(len(result_list)):
		final_signal, final_ali = raw_to_true_signal(result_list[i], 
			seq_list[i], 0, 0, 1, 1)
		write_output(final_signal, arg.output+'_{}.txt'.format(id_list[i]))

