import argparse
import sys
import os
import string
import numpy as np
import math

def arg_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument("--bigram_file", default='./bigram.txt', type=str)
	parser.add_argument("--encode_file", default='./encode.txt', type=str)
	parser.add_argument("--test_file", default='./test.txt', type=str)
	parser.add_argument("--output_file", default='./decode.txt', type=str)
	args = parser.parse_args()

	return args

def load_prob(bigram_file, encode_file, c2i):

	size = len(c2i)
	bigram_table = np.zeros((size, size))
	encode_table = np.zeros((size, size))

	with open(bigram_file, 'r') as f:
		for line in f.readlines():
			info = line.strip().split(' ')
			if len(info) == 3:
				bigram_table[c2i[info[0]], c2i[info[1]]] = float(info[2])
			if len(info) == 2:
				if line.startswith(' '):
					bigram_table[c2i[' '], c2i[info[0]]] = float(info[1])
				else:
					bigram_table[c2i[info[0]], c2i[' ']] = float(info[1])

	# P(z|x) = P(x|z)P(z)/P(x) and P(x) P(x) are same
	with open(encode_file, 'r') as f:
		for line in f.readlines():
			info = line.strip().split(' ')
			if len(info) == 3:
				encode_table[c2i[info[0]], c2i[info[1]]] = float(info[2])
	encode_table[c2i[' '], c2i[' ']] = 1.

	return bigram_table, encode_table

def indexing():

	alphabet = [str(letter) for letter in string.ascii_lowercase]
	number = [str(num) for num in range(10)]
	total = alphabet + number
	total.append(' ')

	c2i = {}
	i2c = {}

	for i, char in enumerate(total):
		c2i[char] = i
		i2c[i] = char

	return c2i, i2c

def load_test(test_file):

	tokens = []

	with open(test_file, 'r') as f:
		for line in f.readlines():
			token = line.strip().split(' ')
			tokens += token

	return tokens

def count_prob(path, b_table, e_mat, p_mat, l_flag):
	prob = 0
	for i, pair in enumerate(path):
		if i == 0:
			continue
		else:
			pair_pre = path[i-1]
			if l_flag == 0:
				if b_table[pair_pre[1], pair[1]] == 0. or e_mat[pair[0], pair[1]] == 0.:
					prob = -float('Inf')
					break
				prob = p_mat[pair_pre[0], pair_pre[1]] + math.log(b_table[pair_pre[1], pair[1]]) + math.log(e_mat[pair[0], pair[1]])
			else:
				if b_table[pair_pre[1], pair[1]] == 0.:
					prob = -float('Inf')
					break
				prob = p_mat[pair_pre[0], pair_pre[1]] + math.log(b_table[pair_pre[1], pair[1]])


	return prob

def max_path(token ,b_table, e_table, c2i, i2c):
	
	size = len(c2i)
	max_path = []
	e_mat = []
	p_mat = []
	l_flag = 0
	
	token += ' '

	for c in token:
		index = c2i[c]
		e_mat.append(e_table[:, index])
		p_mat.append(np.zeros(size))
	e_mat = np.array(e_mat)
	p_mat = np.array(p_mat)

	for p in range(len(token)):
		if p == len(token)-1:
			l_flag = 1
		if p == 0:
			for i in range(size):
				p_mat[0, i] = np.log(e_mat[0, i] * b_table[c2i[' '], i]) if e_mat[0, i] * b_table[c2i[' '], i] != 0. else -float('Inf')
		else:
			max_pair = [(p-1, 0), (p, 0)]
			max_prob = count_prob(max_pair, b_table, e_mat, p_mat, l_flag)
			for pre in range(size):
				for cur in range(size):
					pair = [(p-1, pre), (p, cur)]
					tmp_prob = count_prob(pair, b_table, e_mat, p_mat, l_flag)
					if tmp_prob > max_prob:
						max_prob = tmp_prob
						max_pair = pair
			max_path.append(max_pair[0])

			for cur in range(size):
				if b_table[max_pair[0][1], cur] == 0. or e_mat[p, cur] == 0.:
					p_mat[p, cur] = -float('Inf')
					continue
				p_mat[p, cur] = p_mat[p, cur] + math.log(b_table[max_pair[0][1], cur]) + math.log(e_mat[p, cur])

	ans_str = ''
	for path in max_path:
		ans_str += i2c[path[1]]

	return ans_str


def main():
	
	args = arg_parser()

	c2i, i2c = indexing()

	b_table, e_table = load_prob(args.bigram_file, args.encode_file, c2i)

	tokens = load_test(args.test_file)

	with open(args.output_file, 'w') as f:
		ans_list = []
		print >> sys.stderr, 'total tokens: '+str(len(tokens))
		for i, token in enumerate(tokens):
			ans_str = max_path(token, b_table, e_table, c2i, i2c)
			ans_list.append(ans_str)
			print >> sys.stderr, '\rdone processing '+str(i)+str('/')+str(len(tokens))+' tokens',
		print >> sys.stderr, ''
		out = ' '.join(ans_list)+'\n'
		f.write(out)


if __name__ == '__main__':
	main()

