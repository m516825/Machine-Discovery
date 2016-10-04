#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

void load_info(double** table_ptr, char* path) {
	char s[100], prob_str[100];
	int c1, c2;
	int info;
	FILE *fp = fopen(path, "r");
	while (fscanf(fp,"%d %d %d\n", &c1, &c2, &info) != EOF) {
		table_ptr[c1][c2] = (double) info;
	}
}

void initial(double** bigram_table, double** encode_table, double* prior_array) {
	char* info_path = "./valid/encode.bin";

	load_info(encode_table, info_path);

	//uniform init
	for (int i = 0; i < 37; i++) {
		double tmp_count = 0;
		for (int j = 0; j < 37; j++) {
			bigram_table[i][j] = 1./37.;
			tmp_count += encode_table[i][j];
		}
		for (int j = 0; j < 37; j++) {
			encode_table[i][j] /= tmp_count;
		}
		prior_array[i] = 1./37.;
	}
}

void E_process(double* _p_a, double** b_t,double** e_t,double* p_a, int* ob_dat, int T, double* b_g_sum, double** b_e_sum, double* e_g_sum, double** e_g_sum_only_j) {
	double **alpha, **beta, **gamma;
	double*** epsilon;
	// init space
	alpha = (double**)calloc(T, sizeof(double*));
	beta = (double**)calloc(T, sizeof(double*));
	gamma = (double**)calloc(T, sizeof(double*));
	epsilon = (double***)calloc(T, sizeof(double**));

	for (int t = 0; t < T; t++) {
		alpha[t] = (double*)calloc(37, sizeof(double));
		beta[t] = (double*)calloc(37, sizeof(double));
		gamma[t] = (double*)calloc(37, sizeof(double));
		epsilon[t] = (double**)calloc(37, sizeof(double*));
		for (int i = 0; i < 37; i++) {
			epsilon[t][i] = (double*)calloc(37, sizeof(double));
		}
	}

	// calculate alpha
	for (int t = 0; t < T; t++) {
		if (t == 0) {
			for (int i = 0; i < 37; i++) {
				alpha[t][i] = p_a[i] * e_t[i][ob_dat[t]];
			}
		} else {
			for (int j = 0; j < 37; j++) {
				double tmp_sum = 0;
				for (int i = 0; i < 37; i++) {
					tmp_sum += alpha[t-1][i] * b_t[i][j];
				}
				alpha[t][j] = tmp_sum * e_t[j][ob_dat[t]];
			}
		}
	}

	// calculate beta
	for (int t = T-1; t >= 0; t--) {
		if (t == T-1) {
			for (int i = 0; i < 37; i++) {
				beta[t][i] = 1;
			}
		} else {
			for (int i = 0; i < 37; i++) {
				for (int j = 0; j < 37; j++) {
					beta[t][i] += b_t[i][j] * e_t[j][ob_dat[t+1]] * beta[t+1][j];
				}
			}
		}
	}

	// calculate gamme
	for (int t = 0; t < T; t++) {
		double tmp_sum = 0;
		for (int i = 0; i < 37; i++) {
			tmp_sum += alpha[t][i] * beta[t][i];
		}
		for (int i = 0; i < 37; i++) {
			gamma[t][i] = tmp_sum != 0. ? alpha[t][i] * beta[t][i] / tmp_sum : 0.;
		}	
	}

	// calculate epsilon
	for (int t = 0; t < T-1; t++) {
		double tmp_sum = 0;
		for (int i = 0; i < 37; i++) {
			for (int j = 0; j < 37; j++) {
				tmp_sum += alpha[t][i] * b_t[i][j] * e_t[j][ob_dat[t+1]] * beta[t+1][j];
			}	
		}
		for (int i = 0; i < 37; i++) {
			for (int j = 0; j < 37; j++) {
				epsilon[t][i][j] = tmp_sum != 0. ? alpha[t][i] * b_t[i][j] * e_t[j][ob_dat[t+1]] * beta[t+1][j] / tmp_sum : 0.;
			}
		}
	}

	/* adding info to table' */

	// adding to p_a'
	for (int i = 0; i < 37; i++) {
		_p_a[i] += gamma[0][i];
	}
	// adding to b_t'
	for (int i = 0; i < 37; i++) {
		double g_sum = 0;
		for (int t = 0; t < T-1; t++) {
			g_sum += gamma[t][i];
		}
		for (int j = 0; j < 37; j++) {
			double e_sum = 0;
			for (int t = 0; t < T-1; t++) {
				e_sum += epsilon[t][i][j];
			}
			// _b_t[i][j] += g_sum != 0. ? e_sum / g_sum : 0.;
			b_e_sum[i][j] += (double) e_sum;
		}
		b_g_sum[i] += (double) g_sum;
	}
	// adding e_t'
	for (int i = 0; i < 37; i++) {
		double g_sum = 0;
		for (int t = 0; t < T; t++) {
			g_sum += gamma[t][i];
		}
		for (int j = 0; j < 37; j++) {
			double g_sum_only_j = 0;
			for (int t = 0; t < T; t++) {
				if (ob_dat[t] == j) {
					g_sum_only_j += gamma[t][i];
				} 
			}
			// _e_t[i][j] += g_sum != 0. ? g_sum_only_j / g_sum : 0.;
			e_g_sum_only_j[i][j] += (double) g_sum_only_j;
		}
		e_g_sum[i] += (double) g_sum;
	}

	// free space
	for (int t = 0; t < T; t++) {
		for (int i = 0; i < 37; i++) {
			free(epsilon[t][i]);
		}
		free(alpha[t]);
		free(beta[t]);
		free(gamma[t]);
		free(epsilon[t]);
	}
	free(alpha);
	free(beta);
	free(gamma);
	free(epsilon);
}

void restore_prob(double** bigram_table, double** encode_table, double* prior_array, double** _bigram_table, double** _encode_table, double* _prior_array) {
	for (int i = 0; i < 37; i++) {
		for (int j = 0; j < 37; j++) {
			bigram_table[i][j] = _bigram_table[i][j];
			encode_table[i][j] = _encode_table[i][j];
		}
		prior_array[i] = _prior_array[i];
		_prior_array[i] = 0.;
	}
}

bool stop_process(double** bigram_table, double** encode_table, double* prior_array, double** _bigram_table, double** _encode_table, double* _prior_array) {
	double b_diff = 0;
	double e_diff = 0;
	double p_diff = 0;
	for (int i = 0; i < 37; i++) {
		for (int j = 0; j < 37; j++) {
			b_diff += fabs(bigram_table[i][j] - _bigram_table[i][j]);
			e_diff += fabs(encode_table[i][j] - _encode_table[i][j]);
		}
		p_diff += fabs(prior_array[i] - _prior_array[i]);
	}
	if ((b_diff && e_diff && p_diff) < 0.00001) {
		return true;
	} else {
		printf("b_diff: %lf, e_diff: %lf, p_diff: %lf\n", b_diff, e_diff, p_diff);
		return false;
	}
}

int main(int argc, char** argv) {
	double **bigram_table, **encode_table;
	double *prior_array;

	double **_bigram_table, **_encode_table;
	double *_prior_array;
	int *observed_dat;
	int T = 50;
	int split_num = 0;
	int num;
	char* test_path = "./valid/test.num";

	observed_dat = (int*)calloc(T, sizeof(int));
	bigram_table = (double**)calloc(37, sizeof(double*));
	encode_table = (double**)calloc(37, sizeof(double*));
	prior_array = (double*)calloc(37, sizeof(double));

	_bigram_table = (double**)calloc(37, sizeof(double*));
	_encode_table = (double**)calloc(37, sizeof(double*));
	_prior_array = (double*)calloc(37, sizeof(double));

	for (int i = 0; i < 37; i++) {
		bigram_table[i] = (double*)calloc(37, sizeof(double));
		encode_table[i] = (double*)calloc(37, sizeof(double));

		_bigram_table[i] = (double*)calloc(37, sizeof(double));
		_encode_table[i] = (double*)calloc(37, sizeof(double));
	}

	double* b_gamma_sum = (double*)calloc(37, sizeof(double));
	double* e_gamma_sum = (double*)calloc(37, sizeof(double));
	double** b_epsilon_sum = (double**)calloc(37, sizeof(double*));
	double** e_gamma_sum_only_j = (double**)calloc(37, sizeof(double*));
	for (int i = 0; i < 37; i++) {
		b_epsilon_sum[i] = (double*)calloc(37, sizeof(double));
		e_gamma_sum_only_j[i] = (double*)calloc(37, sizeof(double));
	}

	initial(bigram_table, encode_table, prior_array);

	while(1) {
		FILE* fp = fopen(test_path, "r");
		int loaded = 0;
	
		for (int i = 0; i < 37; i++) {
			for (int j = 0; j < 37; j++) {
				b_epsilon_sum[i][j] = 0;
				e_gamma_sum_only_j[i][j] = 0;
			}
			b_gamma_sum[i] = 0;
			e_gamma_sum[i] = 0;
		}

		while(fscanf(fp, "%d ", &num) != EOF) {
			observed_dat[loaded] = num;
			loaded += 1;
			if (loaded == 50) {
				//
				E_process(_prior_array, bigram_table, encode_table, prior_array, observed_dat, loaded, b_gamma_sum, b_epsilon_sum, e_gamma_sum, e_gamma_sum_only_j);
				// next token
				loaded = 0;
				split_num += 1;
			}
		}
		if (loaded > 0) {
			E_process(_prior_array, bigram_table, encode_table, prior_array, observed_dat, loaded, b_gamma_sum, b_epsilon_sum, e_gamma_sum, e_gamma_sum_only_j);
			split_num += 1;
		}

		// calculate prob after E process
		for (int i = 0; i < 37; i++) {
			for (int j = 0; j < 37; j++) {
				_bigram_table[i][j] = b_epsilon_sum[i][j] / b_gamma_sum[i];
				_encode_table[i][j] = e_gamma_sum_only_j[i][j] / e_gamma_sum[i];
			}
			_prior_array[i] /= (double) split_num;
		}

		/*
			compare two prob table here!!!!
		*/
		if (stop_process(bigram_table, encode_table, prior_array, _bigram_table, _encode_table, _prior_array)) {
			break;
		}

		// restore new prob table
		restore_prob(bigram_table, encode_table, prior_array, _bigram_table, _encode_table, _prior_array);

	}

	return 0;
}

