#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int c2i(char c) {
	if ((int) c == 32)
		return 36;
	if ((int) c <= 122 && (int) c >= 97)
		return ((int) c) - 97;
	if ((int) c <= 57 && (int) c >= 48)
		return ((int) c) - 48 + 26;
	return 0;
}

char i2c(int i) {
	if (i == 36)
		return ' ';
	if (i <= 25 && (int) i >= 0)
		return (char) (i + 97);
	if (i <= 35 && (int) i >= 26)
		return (char) (i - 26 + 48);
	return ' ';
}

void load_prob(double** table_ptr, char* path) {
	char s[100], prob_str[100];
	char c1, c2;
	double prob;
	FILE *fp = fopen(path, "r");
	while(fgets(s, 100,fp) != NULL) {
		s[strlen(s)-1] = '\0';
		c1 = s[0];
		c2 = s[2];
		strcpy(prob_str, &s[4]);
		prob = (double) atof(prob_str);
		table_ptr[c2i(c1)][c2i(c2)] = prob;
	}
}

double calculate_prob(double p_mat[37], double e_mat[37], double** bigram_table, int c_1, int c) {
	double ans = (bigram_table[c_1][c] && e_mat[c]) != 0.0 ? p_mat[c_1] + log(bigram_table[c_1][c]) + log(e_mat[c]) : -INFINITY;
	return ans;
}

char* viterbi(double** bigram_table, double** encode_table, char* token) {
	double e_mat[strlen(token)][37], p_mat[strlen(token)][37];
	int path_mat[strlen(token)][37];

	for (int t = 0; t < strlen(token); t++) {
		for (int c = 0; c < 37; c++) {
			e_mat[t][c] = encode_table[c][c2i(token[t])];
			path_mat[t][c] = -1;
		}
	}

	for (int t = 0; t < strlen(token); t++) {
		if (t == 0) {
			for (int c = 0; c < 37; c++) {
				p_mat[t][c] = (e_mat[t][c] && bigram_table[c2i(' ')][c]) != 0.0 ? log(e_mat[t][c]) + log(bigram_table[c2i(' ')][c]) : -INFINITY;
			}
		} else {
			for (int c = 0; c < 37; c++) {
				int max_pre = 0;
				double max_prob = -INFINITY;
				for (int c_1 = 0; c_1 < 37; c_1++) {
					double tmp_prob = calculate_prob(p_mat[t-1], e_mat[t], bigram_table, c_1, c);
					if (tmp_prob > max_prob) {
						max_prob = tmp_prob;
						max_pre = c_1;
					}
				}
				path_mat[t][c] = max_pre;

				int c_1 = path_mat[t][c];
				p_mat[t][c] = (bigram_table[c_1][c] && e_mat[t][c]) != 0.0 ? p_mat[t-1][c_1] + log(bigram_table[c_1][c]) + log(e_mat[t][c]) : -INFINITY;
			}
		}
	}
	int max_index;
	double max_prob = -INFINITY;
	for (int c = 0; c < 37; c++) {
		if (p_mat[strlen(token)-1][c] > max_prob) {
			max_prob = p_mat[strlen(token)-1][c];
			max_index = c;
		}
	}
	char* ans_str = malloc(sizeof(char)*(strlen(token)-1));
	for (int t = strlen(token)-1; t >= 1; t--) {
		int c_index = path_mat[t][max_index];
		ans_str[t-1] = i2c(c_index);
		max_index = c_index;
	}
	ans_str[strlen(token)-1] = '\0';
	return ans_str;
}

int main(int argc, char** argv) {

	char* bigram_path = "./src/bigram.txt";
	char* encode_path = "./src/encode.txt";
	char* test_path = "./src/test.txt";
	char* output_path = "./pre.txt";
	double **bigram_table, **encode_table;
	char token[50], inchar[2];
	inchar[1] = '\0';
	token[0] = '\0';
	bigram_table = (double**)malloc(sizeof(double*)*37);
	for (int i = 0; i < 37 ; i++) {
		bigram_table[i] = (double*)malloc(sizeof(double)*37);
	}
	encode_table = (double**)malloc(sizeof(double*)*37);
	for (int i = 0; i < 37 ; i++) {
		encode_table[i] = (double*)malloc(sizeof(double)*37);
	}

	load_prob(bigram_table, bigram_path);
	load_prob(encode_table, encode_path);

	FILE *fp = fopen(test_path, "r");
	FILE *fp_out = fopen(output_path, "w");
	while(fscanf(fp,"%c", &inchar[0]) != EOF) {
		strcat(token, inchar);
		if (inchar[0] == ' ') {
			char* ans_str = viterbi(bigram_table, encode_table, token);
			fprintf(fp_out, "%s ", ans_str);
			free(ans_str);
			token[0] = '\0';
		}
	}
	fprintf(fp_out, "\n");

	return 0;
}