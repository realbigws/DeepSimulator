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

