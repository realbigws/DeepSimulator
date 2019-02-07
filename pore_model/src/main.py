from poremodel_util import *
import numpy as np
from multiprocessing import Pool
import multiprocessing
import os
import argparse
import tqdm


#=================== main =======================#
if __name__ == '__main__':
	os.environ['CUDA_VISIBLE_DEVICES'] ='-1'

	#--------- argument part -------------#
	parser = argparse.ArgumentParser(description='convert the \
		input sequence to nanopore signal')
	parser.add_argument('-i', action='store', dest='input', required=True, 
		help='the input file')
	parser.add_argument('-p', action='store', dest='output', required=True,
		help='prefix the output file')
	parser.add_argument('-l', action='store', dest='alignment', required=True,
		help='prefix the alignment file')
	parser.add_argument('-t', action='store', dest='threads',
		type=int, help='the number of threads used', default=1)
	parser.add_argument('-a', action='store', dest='alpha',
		type=float, help='change the distribution of the signal repeat time, \
		value between 0 and 1, 0.1 (default) would give the distribution best \
		simulate the real case, 0 would give distribution whose basecalling \
		result is slightly worse than the real case, 1 would give the almost \
		perfect basecalling result using Albacore', default=0.1)
	parser.add_argument('-s', action='store', dest='std',
		type=float, help='set the std of the signal. \
		The higher the value, the blurred the signal. (default is 1.5)',
		default=1.5)
	parser.add_argument('-f', action='store', dest='freq',
		type=float, help='change the cut frequency in the low pass filter. \
		The higher the value, the smoother the signal. (default is 850)',
		default=850)
	parser.add_argument('--perfect', action='store', dest='perfect',
		type=bool, help='Do you want a perfect signal and sequence',
		default=False)
	parser.add_argument('--perflen', action='store', dest='perflen',
		type=int, help='repeat length for perfect mode',
		default=6)

	#---------- input list ---------------#
	arg = parser.parse_args()
	seq_list = get_seq_list(arg.input)
	id_list = get_id_list(arg.input)
	seq_chunk_list = convert_to_input(seq_list)

	#---------- deep simulator -----------#
	result_list = []
	p = Pool(arg.threads)
	result_list = list(tqdm.tqdm(
		p.imap(model_whole_set_check, seq_chunk_list), 
		total=len(seq_chunk_list)))
	p.close()
	p.join()

	#---------- output results -----------#
	for i in range(len(result_list)):
		final_signal, final_ali = raw_to_true_signal(result_list[i], 
			seq_list[i], arg.alpha, arg.freq, arg.std, arg.perfect, arg.perflen)
		write_output(final_signal, arg.output+'_{}.txt'.format(id_list[i]))
		write_alignment(final_ali, arg.alignment+'_{}.ali'.format(id_list[i]))

