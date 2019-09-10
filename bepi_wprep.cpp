#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <memory.h>
#include <limits.h>
#include <pthread.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include "src/mergesort.h"
#include "src/queue.h"
#include "src/bmpcreator.h"
#include "src/LU.h"

#define DELIMITER '\t'
#define MAX_NODES 50000000
#define MAX_EDGES 2000000000 
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define WIDTH  1024
#define HEIGHT 1024
#define RESTART 0.15
#define FILTER_THRESHOLD 100
#define SLASH_RATIO 100
#define BIG_BLOCK_THRESHOLD 1000

int nthreads = 16;

#include "paralution.hpp"
using namespace paralution;

/* pthread helper */
typedef struct {
	int *array;
	int *subarray;
	int p, q, dir;
} arg_struct;

typedef struct {
	size_t part;
	size_t block_count;
	size_t *block_idx;
	size_t *h11_idx;
	double **H11;
	int *h11_row, *h11_col;
	double *h11_data;
} arg_inv;

void readFile(char*);
void loadEdges(FILE*, int*, int*, float*, int*, int*, int*, size_t*);
void* merge_sort_t(void *arg);
void merge_sort_pthread(int*, int*, int, int, int);
void trp_push_back(size_t*, size_t*, int**, int**, double**, int, int, double);
void* H_inv_t(void *);
void LUinverse(double *, size_t);

/*************************************/
long long current_timestamp() {
	struct timeval te; 
	gettimeofday(&te, NULL); // get current time
	return te.tv_sec*1000LL + te.tv_usec/1000;
}
/************************************/

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: ./cogi_prep <filename.tsv>\n");
		return -1;
	}

	mkdir("prep", 0700);

	long long starttime, endtime;
	starttime = current_timestamp();
	readFile(argv[1]);

	endtime = current_timestamp();
	fprintf(stdout, ">>> Process terminated successfully. (took %.3f sec)\n",
			  (endtime-starttime)/1000.0);

	return 0;
}

void readFile(char* filename) {

	FILE *inFile = fopen(filename, "r");
	int max_index = 0;
	size_t edge_count = 0;
	int i, j, k, u, v, tmpval, tmp, top, bottom, topfront, topmid;
	int iteration, deadend_count, ptop, pbottom;
	int max_degree, max_idx;
	int cclabel_count, nlabel, visit, curtag;
	int init_length, count;
	int tmpidx;
	int SLASH_K;


	if (inFile == NULL){
		fprintf(stderr, "File %s not found!\n", filename);
		return;
	}
	
	fprintf(stdout, ">>> File opened successfully.\n");

	int *start_nf = (int*) malloc (sizeof(int) * MAX_EDGES);
	int *start_nt = (int*) malloc (sizeof(int) * MAX_EDGES);
	float *start_vv = (float*) malloc (sizeof(float) * MAX_EDGES);
	int *in_degree  = (int*) malloc (sizeof(int) * MAX_NODES);
	int *out_degree = (int*) malloc (sizeof(int) * MAX_NODES);

	/* Start reading file to add edge data on Memory */
	loadEdges(inFile, 
			  start_nf, 
			  start_nt,
			  start_vv,
			  in_degree, 
			  out_degree, 
			  &max_index, 
			  &edge_count);
	fclose(inFile);
	fprintf(stdout, ">>>   Max index: %d\n", max_index);
	fprintf(stdout, ">>>   Edges: %zu\n", edge_count);
	in_degree  = (int*) realloc (in_degree, sizeof(int) * (max_index + 1));
	out_degree = (int*) realloc (out_degree, sizeof(int) * (max_index + 1));

	//int *nf  = (int*) malloc (sizeof(int) * edge_count);
	int *nt  = (int*) malloc (sizeof(int) * edge_count);
	//int *nfR = (int*) malloc (sizeof(int) * edge_count);
	int *ntR = (int*) malloc (sizeof(int) * edge_count);
	int *nf_idx  = (int*) malloc (sizeof(int) * (max_index + 2));
	int *nfR_idx = (int*) malloc (sizeof(int) * (max_index + 2));
	int *off1 = (int*) calloc (1, sizeof(int) * (max_index + 1));
	int *off2 = (int*) calloc (1, sizeof(int) * (max_index + 1));

	float *ev = (float*) malloc (sizeof(float) * edge_count);
	//start_vv = (float*) realloc (start_vv, sizeof(float) * edge_count);

	nfR_idx[0] = nf_idx[0] = 0;
	for (i = 0; i <= max_index; i++){
		nf_idx[i+1] = out_degree[i] + nf_idx[i];
		nfR_idx[i+1] = in_degree[i] + nfR_idx[i];
	}

	for (i = 0; i < edge_count; i++){
		u = start_nf[i];
		v = start_nt[i];

		//nf[nf_idx[u] + off1[u]] = u;
		nt[nf_idx[u] + off1[u]] = v;
		ev[nf_idx[u] + off1[u]] = start_vv[i];

		//nfR[nfR_idx[v] + off2[v]] = v;
		ntR[nfR_idx[v] + off2[v]] = u;

		off1[u]++;
		off2[v]++;
	}

	free(start_nf);
	free(start_nt);
	free(start_vv);


	int *reordered  = (int*) malloc (sizeof(int) * (max_index + 1));
	int *new_index  = (int*) malloc (sizeof(int) * (max_index + 1));
	int *tag_order  = (int*) malloc (sizeof(int) * (max_index + 1));
	int *CCtag      = (int*) calloc (sizeof(int) , (max_index + 1));
	int *CCcount    = (int*) calloc (sizeof(int) , (max_index + 1));
	int *CClabel    = (int*) malloc (sizeof(int) * (max_index + 1));
	int *CCrank     = (int*) malloc (sizeof(int) * (max_index + 1));
	int *tot_degree = (int*) malloc (sizeof(int) * (max_index + 1));
	int *tmp_degree = (int*) malloc (sizeof(int) * (max_index + 1));

	queue *mQueue   = init_queue(edge_count);


	top = 0;
	bottom = max_index;
	deadend_count = 0;

	/* Deadends */
	for (i = 0; i <= max_index; i++){
		if (out_degree[i] == 0){
			CCtag[i] = -1;
			reordered[bottom] = i;
			new_index[i] = bottom;
			deadend_count++;
			bottom--;
		}
		else{
			reordered[top] = i;
			top++;
		}
		tmp_degree[i] = tot_degree[i] = out_degree[i] + in_degree[i];
	}

	for (i = 0; i <= max_index; i++){
		for (j = nf_idx[i], tmp = nf_idx[i+1] - 1; j <= tmp;){
			if (CCtag[nt[j]] == -1) { // Prevent revisit for efficiency
				tmpval = nt[tmp];
				nt[tmp] = nt[j];
				nt[j] = tmpval;
				tmp --;

				off1[i]--;
			}
			else {
				j++;
			}
		}
	}

	top = 0;

	SLASH_K = MAX(1, (max_index - deadend_count)/100);
	fprintf(stdout, ">>>   SLASH_K = %d\n", SLASH_K);

	assert(deadend_count == max_index - bottom);
	fprintf(stdout, ">>>   Deadend count: %d\n", deadend_count);

	top = 0;
	iteration = 0;
	nlabel = cclabel_count = 0;

	while(top <= bottom){
		iteration++;
		ptop = top;
		pbottom = bottom;
		fprintf(stdout, ">>>     Iteration: %2d [%.2f%%] (%d ~ %d)\n", iteration, 
			((double)(top + (max_index - bottom))*100.0/(double) max_index), top, bottom);
		
		merge_sort_pthread(tmp_degree, reordered, top, bottom, 1);

		for (i = 0; i < SLASH_K && top <= bottom; i++){
			CCtag[reordered[bottom]] = -1;
			new_index[reordered[bottom]] = bottom;
			bottom--;
		}

		if (top > bottom) break;

		/* this one is another bottleneck but hard to code in pthread */
		nlabel = cclabel_count;
		for (i = 0; i <= max_index; i++){
			if (CCtag[i]) continue;
			curtag = ++cclabel_count;
			CClabel[curtag] = curtag;
			queue_put(mQueue, i);
			while (!is_queue_empty(mQueue)){
				visit = queue_get(mQueue);
				if (CCtag[visit]) continue;
				CCtag[visit] = curtag;
				CCcount[curtag]++;
				
				for (j = nf_idx[visit], tmp = nf_idx[visit] + off1[visit] - 1; j <= tmp;){
					if (CCtag[nt[j]] == 0) {
						queue_put(mQueue, nt[j]);
						j++;
					}
					else if (CCtag[nt[j]] == -1) { // Prevent revisit for efficiency
						tmpval = nt[tmp];
						nt[tmp] = nt[j];
						nt[j] = tmpval;
						
						tmp --;
						off1[visit]--;

						// do not increase j
					}
					else {
						j++;
					}
				}
				for (j = nfR_idx[visit], tmp = nfR_idx[visit] + off2[visit] - 1; j <= tmp;){
					if (CCtag[ntR[j]] == 0) {
						queue_put(mQueue, ntR[j]);
						j++;
					}
					else if (CCtag[ntR[j]] == -1) { // Prevent revisit for efficiency
						tmpval = ntR[tmp];
						ntR[tmp] = ntR[j];
						ntR[j] = tmpval;
						
						tmp --;
						off2[visit]--;

						// do not increase j
					}
					else {
						j++;
					}
				}
			}
		}

		merge_sort_pthread(CCcount, CClabel, nlabel + 1, cclabel_count, -1);
		int label_big_index = nlabel + 1;
		for (i = nlabel+1; i <= cclabel_count; i++){
			CCrank[CClabel[i]] = (cclabel_count + 1 - i) + 1;
			if (CCcount[CClabel[i]] > BIG_BLOCK_THRESHOLD) label_big_index = i;
		}
		//printf("nlabel + 1 = %d, label_big_idx = %d\n", nlabel+1, label_big_index);
		//printf("blocks = %d\n", cclabel_count - nlabel);

		topfront = top;
		//topmid = bottom - CCcount[nlabel + 1] + 1;
		topmid = bottom;

		for (i = 0; i <= max_index; i++){
			if (CCtag[i] > nlabel){
				bool slash_again = false;
				for (k = nlabel+1; k <= label_big_index; k++){
					if (CCtag[i] == CClabel[k])
						slash_again = true;
				}
				//if(CCtag[i] == CClabel[nlabel + 1] || (true && nlabel + 2 <= cclabel_count && CCtag[i] == CClabel[nlabel + 2])){
				if(slash_again){
					CCtag[i] = 0;
					reordered[topmid] = i;
					//new_index[i] = topmid;
					tmp_degree[topmid] = tot_degree[i];
					//topmid++;
					topmid--;
				}
				else {
					//CCtag[i] = -1;
					reordered[topfront] = i;
					tag_order[topfront] = CCrank[CCtag[i]];
					topfront++;
				}
			}
		}

		//assert(topfront == bottom - CCcount[nlabel+1] + 1);
		//assert(topmid == bottom + 1);
		//printf("topmid = %d, topfront = %d\n", topmid, topfront);
		//assert(topmid + 1 == topfront);


		//merge_sort_pthread(tag_order, reordered, top, bottom - CCcount[nlabel+1], 1);
		merge_sort_pthread(tag_order, reordered, top, topfront - 1, 1);
		while (top < topfront){
			new_index[reordered[top]] = top;
			top++;
		}

		//top = bottom - CCcount[nlabel+1] + 1;
		//assert(top == bottom - CCcount[nlabel+1] + 1);
		assert(top == topfront);

	}

	fprintf(stdout, ">>>   Iteration over.\n");


	bottom = max_index + 1 - deadend_count;
	top = bottom/10*7 + 1; //(max_index - deadend_count) * 7 / 10; // do as you want
	while(CCtag[reordered[top]] == -1) top--;
	top ++;
	printf(">>> top = %d, bottom = %d, max_index = %d\n", top, bottom, max_index);

	/*
	FILE *CCtag_file = fopen("CCtag_list.txt", "w");
	fprintf(CCtag_file, "%d\t%d\t%d\t%d\n", i, reordered[i], CCtag[reordered[i]], CCcount[CCtag[reordered[i]]]);
	for (i = 0; i <= cclabel_count; i++){
		fprintf(CCtag_file, "%d\t%d\t%d\t%d\n", i, reordered[i], CCtag[reordered[i]], CCcount[CCtag[reordered[i]]]);
	}
	fclose(CCtag_file);

	int *CCcount_cpy   = (int*) malloc (sizeof(int) * (cclabel_count + 1));
	int *CCnumber_cpy  = (int*) malloc (sizeof(int) * (cclabel_count + 1));
	memcpy(CCcount_cpy, CCcount, sizeof(int)*(cclabel_count + 1));
	for (i = 0; i <= cclabel_count; i++) CCnumber_cpy[i] = i;
	merge_sort_pthread(CCcount_cpy, CCnumber_cpy, 0, cclabel_count+1, -1);
	FILE *CCcount_file = fopen("CCcount_list.txt", "w");
	for (i = 0; i <= cclabel_count; i++){
		//if (CCcount_cpy[i] == 0) break;
		fprintf(CCcount_file, "%d\t%d\n", CCnumber_cpy[i], CCcount_cpy[i]);
	}
	fclose(CCcount_file);
	*/

	FILE *CCtag_file = fopen("prep/CCtag_list.txt", "w");
	FILE *CCtag_file_big = fopen("prep/CCtag_list_big.txt", "w");
	FILE *CCcount_file = fopen("prep/CCcount_list.txt", "w");
	fprintf(CCtag_file, "Node\tBlockIdx\tCount\n");
	fprintf(CCtag_file_big, "Node\tBlockIdx\tCount\n");
	fprintf(CCcount_file, "BlockIdx\tCount\n");

	/* H11 block count and each size */
	size_t block_count = 0;
	size_t *block_idx = (size_t*) malloc (sizeof(size_t)*(max_index - deadend_count + 1));
	size_t *h11_idx = (size_t*) malloc(sizeof(size_t)*(max_index - deadend_count) + 1); // realloc needed
	size_t *h11_size = (size_t*) malloc(sizeof(size_t)*(max_index - deadend_count) + 1); // realloc needed
	double **H11 = (double**) malloc (sizeof(double*)*(max_index - deadend_count));
	size_t mtmp = 1;
	size_t h11_tot = 0;
	h11_idx[0] = block_idx[0] = 0;
	j = 0;
	for (i = 1; i <= top; i++){
		u = reordered[i];
		v = reordered[i-1];
		if (CCtag[u] != CCtag[v] || i == top) {
			block_idx[block_count+1] = block_idx[block_count] + mtmp;
			h11_idx[block_count+1] = h11_idx[block_count] + mtmp*mtmp;
			h11_size[block_count] = mtmp;
			//printf("%zu ", h11_idx[block_count+1]);
			H11[block_count] = (double*) calloc (sizeof(double), mtmp*mtmp);
			h11_tot += mtmp*mtmp;

			for (; j < i; j++){
				fprintf(CCtag_file, "%d\t%d\t%zu\n", reordered[j], block_count, mtmp);
				if (mtmp >= FILTER_THRESHOLD){
					fprintf(CCtag_file_big, "%d\t%d\t%zu\n", reordered[j], block_count, mtmp);
				}
			}
			fprintf(CCcount_file, "%d\t%zu\n", block_count, mtmp);

			block_count++;
			//printf("%zu ", mtmp);

			
			if (CCtag[u] == -1 || i == top) break;
			mtmp = 1;
		}
		else {
			mtmp ++;
		}
	}

	fclose(CCtag_file);
	fclose(CCtag_file_big);
	fclose(CCcount_file);

	/* Initialize matrixes into sparse triplets */
	size_t h11_nnz = 0, h11_max = h11_tot;
	//size_t h11_nnz = 0, h11_max = edge_count / 20; //h11_tot;
	size_t h12_nnz = 0, h12_max = edge_count / 5;
	size_t h21_nnz = 0, h21_max = edge_count / 5;
	size_t h22_nnz = 0, h22_max = edge_count / 5;
	size_t h31_nnz = 0, h31_max = edge_count / 5;
	size_t h32_nnz = 0, h32_max = edge_count / 5;

	int *h11_row = (int*) malloc(sizeof(int)*h11_max);
	int *h11_col = (int*) malloc(sizeof(int)*h11_max);
	double *h11_data = (double*) malloc(sizeof(double)*h11_max);

	int *h12_row = (int*) malloc(sizeof(int)*h12_max);
	int *h12_col = (int*) malloc(sizeof(int)*h12_max);
	double *h12_data = (double*) malloc(sizeof(double)*h12_max);

	int *h21_row = (int*) malloc(sizeof(int)*h21_max);
	int *h21_col = (int*) malloc(sizeof(int)*h21_max);
	double *h21_data = (double*) malloc(sizeof(double)*h21_max);

	int *h22_row = (int*) malloc(sizeof(int)*h22_max);
	int *h22_col = (int*) malloc(sizeof(int)*h22_max);
	double *h22_data = (double*) malloc(sizeof(double)*h22_max);

	int *h31_row = (int*) malloc(sizeof(int)*h31_max);
	int *h31_col = (int*) malloc(sizeof(int)*h31_max);
	double *h31_data = (double*) malloc(sizeof(double)*h31_max);

	int *h32_row = (int*) malloc(sizeof(int)*h32_max);
	int *h32_col = (int*) malloc(sizeof(int)*h32_max);
	double *h32_data = (double*) malloc(sizeof(double)*h32_max);

	/* Initialize completed */
	char *red   = (char*) malloc (sizeof(char)*WIDTH*HEIGHT);
	char *green = (char*) malloc (sizeof(char)*WIDTH*HEIGHT);
	char *blue  = (char*) malloc (sizeof(char)*WIDTH*HEIGHT);
	memset(red,   0xFF, sizeof(char)*WIDTH*HEIGHT);
	memset(green, 0xFF, sizeof(char)*WIDTH*HEIGHT);
	memset(blue,  0xFF, sizeof(char)*WIDTH*HEIGHT);
	//vector< T > trp11, trp12, trp21, trp22, trp31, trp32, trpH;

	size_t bidx = 0;
	double val;
	for (i = 0; i <= max_index && i < bottom; i++){ // else deadend
		u = reordered[i];

		if (i >= block_idx[bidx+1] && bidx < block_count){
			bidx ++;
		}

		if (i < top){
			H11[bidx][(i-block_idx[bidx])*(block_idx[bidx+1]-block_idx[bidx])
						+ (i-block_idx[bidx])] = 1.0;
		}
		else if (i < bottom){
			//trp22.push_back(T(i-block_idx[bidx], i-block_idx[bidx], 1.0));
			trp_push_back(&h22_nnz, &h22_max, &h22_row, &h22_col, &h22_data,
							i-block_idx[bidx], i-block_idx[bidx], 1.0);
		}

		for (j = nf_idx[u]; j < nf_idx[u+1]; j++){
			v = new_index[nt[j]];

			size_t offiv = (((size_t)i)*HEIGHT/max_index)*WIDTH + (((size_t)v)*WIDTH/max_index);
			*(red + offiv ) = *(green + offiv ) = 0;

			//val = -(1-RESTART)/(double)(out_degree[u]);
			val = -(1-RESTART)*ev[j];//esum[u];
			if (i < top) {
				if (v < top) {
					//printf("H11[%d][%d,%d]\n", bidx, v-block_idx[bidx], i-block_idx[bidx]);
					H11[bidx][(v-block_idx[bidx])*(block_idx[bidx+1]-block_idx[bidx])
								+ (i-block_idx[bidx])] = val;
				}
				else if (v < bottom) {
					//trp21.push_back(T(v - top, i, val));
					trp_push_back(&h21_nnz, &h21_max, &h21_row, &h21_col, &h21_data,
									v-top, i, val);
				}
				else {
					//trp31.push_back(T(v - bottom, i, val));
					//printf("3[%zu/%zu]", h31_nnz, h31_max);
					trp_push_back(&h31_nnz, &h31_max, &h31_row, &h31_col, &h31_data,
									v-bottom, i, val);
				}
			}
			else if (i < bottom) {
				if (v < top) {
					//trp12.push_back(T(v, i - top, val));
					trp_push_back(&h12_nnz, &h12_max, &h12_row, &h12_col, &h12_data,
									v, i-top, val);
				}
				else if (v < bottom) {
					//trp22.push_back(T(v - top, i - top, val));
					trp_push_back(&h22_nnz, &h22_max, &h22_row, &h22_col, &h22_data,
									v-top, i-top, val);
				}
				else {
					//trp32.push_back(T(v - bottom, i - top, val));
					trp_push_back(&h32_nnz, &h32_max, &h32_row, &h32_col, &h32_data,
									v-bottom, i-top, val);
				}
			}
			else {
				fprintf(stderr, "ERROR\n");
				return;
			}
		}
	}

	printf(">>> Inserting values into Paralution...\n");
	merge_sort_pthread(tmp_degree, reordered, top, bottom, 1);
	printf(">>>   H11: %zu\n"
		   "      H12: %zu\n"
		   "      H21: %zu\n"
		   "      H22: %zu\n"
		   "      H31: %zu\n"
		   "      H32: %zu\n",
			h11_tot, h12_nnz, h21_nnz, h22_nnz, h31_nnz, h32_nnz);

	char resultbmp[100] = "prep/adjmat.bmp";
	drawbmp (resultbmp, WIDTH, HEIGHT, red, green, blue);
	printf(">>> Plot [%s] created\n", resultbmp);
	free(red); free(green); free(blue);

	printf(">>> Inversing H11...\n");
	init_paralution();
	info_paralution();
	/*
	pthread_t threads_inv[nthreads];
	arg_inv argi[nthreads];
	for (i = 0; i < nthreads; i++){
		argi[i].part = (size_t)i;
		argi[i].block_count = block_count;
		argi[i].block_idx = block_idx;
		argi[i].h11_idx = h11_idx;
		argi[i].H11 = H11;
		argi[i].h11_row = h11_row;
		argi[i].h11_col = h11_col;
		argi[i].h11_data = h11_data;

		pthread_create(&threads_inv[i], NULL, H_inv_t, (void *)&argi[i]);
	}

	for (i = 0; i < nthreads; i++){
		pthread_join(threads_inv[i], NULL);
	}
	*/

	h11_nnz = 0;
    for(bidx = 0; bidx < block_count; bidx++){
		if (h11_size[bidx] < 100 || true) {
	        LUinverse(H11[bidx], h11_size[bidx]);
			// H11 is now inversed for small size blocks.
			for (j = 0; j < h11_size[bidx]; j++){
				for (k = 0; k < h11_size[bidx]; k++){
					val =  H11[bidx][j * h11_size[bidx] + k];
					if (val != 0){
						h11_row[h11_nnz] = block_idx[bidx] + j;
						h11_col[h11_nnz] = block_idx[bidx] + k;
						h11_data[h11_nnz] = val;
						h11_nnz++;
					}
				}
			}
		}
		else {
			printf("... matrix size = %d\n", h11_size[bidx]);
			size_t lin_max = h11_size[bidx] * h11_size[bidx];
			int *lin_row = (int*) malloc (sizeof(int)*lin_max);
			int *lin_col = (int*) malloc (sizeof(int)*lin_max);
			double *lin_data = (double*) malloc (sizeof(double)*lin_max);
			size_t lin_nnz = 0;

			// Don't confused!!!* not yet inversed.
			for (j = 0; j < h11_size[bidx]; j++){
				for (k = 0; k < h11_size[bidx]; k++){
					val =  H11[bidx][j * h11_size[bidx] + k];
					if (val != 0){
						lin_row[h11_nnz] = block_idx[bidx] + j;
						lin_col[h11_nnz] = block_idx[bidx] + k;
						lin_data[h11_nnz] = val;
						lin_nnz++;
					}
				}
			}

			LocalMatrix<double> lin_mat;
			lin_mat.SetDataPtrCOO(&lin_row, &lin_col, &lin_data, "lin",
									lin_nnz, h11_size[bidx], h11_size[bidx]);
			/*
			Chebyshev<LocalMatrix<double>, LocalVector<double>, double> ls;
			ls.SetOperator(lin_mat);
			ls.Set(-1, -1);
			ls.Build();
			*/
			GMRES<LocalMatrix<double>, LocalVector<double>, double> ls;
			ls.SetOperator(lin_mat);
			ls.Build();

			LocalVector<double> lin_rhs, lin_x;
			lin_rhs.Allocate("lin_rhs", h11_size[bidx]);
			lin_x.Allocate("lin_x", h11_size[bidx]);

			for (k = 0; k < h11_size[bidx]; k++){
				lin_rhs.Zeros();
				lin_rhs[k] = 1.0;
				ls.Solve(lin_rhs, &lin_x);
				
				for(j = 0; j < h11_size[bidx]; j++){
					val = lin_x[j];
					if (val != 0){
						h11_row[h11_nnz] = block_idx[bidx] + j;
						h11_col[h11_nnz] = block_idx[bidx] + k;
						h11_data[h11_nnz] = val;
						h11_nnz++;
						printf("%f\n", val);
					}
				}
			}
			lin_rhs.Clear();
			lin_x.Clear();
			lin_mat.Clear();
			
		}

    }
    //h11_nnz = h11_max;
	//assert(h11_nnz == h11_max);

	free(block_idx);
	for (i = 0; i < block_count; i++) free(H11[i]);
	free(H11);
	printf(">>> Inversing H11 completed.\n");


	LocalVector<int> save_new_index, save_reordered;
	save_new_index.SetDataPtr(&new_index, "new_index", max_index+1);
	save_reordered.SetDataPtr(&reordered, "reordered", max_index+1);
	save_new_index.WriteFileASCII("prep/new_index.ASCII");
	save_reordered.WriteFileASCII("prep/reordered.ASCII");
	save_new_index.Clear();
	save_reordered.Clear();

	
	/*
	int *h11_row_r = (int*) malloc(sizeof(int)*h11_max);
	int *h11_col_r = (int*) malloc(sizeof(int)*h11_max);
	double *h11_data_r = (double*) malloc(sizeof(double)*h11_max);

	h11_nnz = 0;
	for (size_t si = 0; si < h11_max; si++){
		if (h11_data[si] != 0){
			h11_row_r[h11_nnz] = h11_row[si];
			h11_col_r[h11_nnz] = h11_col[si];
			h11_data_r[h11_nnz] = h11_data[si];
			h11_nnz++;
		}
	}

	free(h11_row);
	free(h11_col);
	free(h11_data);
	*/

	LocalMatrix<double> Hv11, H12, H21, H22, H31, H32, S, St1, St2;
	Hv11.SetDataPtrCOO(&h11_row, &h11_col, &h11_data, "H11 inv",
						h11_nnz, top, top);
	H12.SetDataPtrCOO(&h12_row, &h12_col, &h12_data, "H12",
						h12_nnz, top, bottom-top);
	H21.SetDataPtrCOO(&h21_row, &h21_col, &h21_data, "H21",
						h21_nnz, bottom-top, top);
	H22.SetDataPtrCOO(&h22_row, &h22_col, &h22_data, "H22",
						h22_nnz, bottom-top, bottom-top);
	H31.SetDataPtrCOO(&h31_row, &h31_col, &h31_data, "H31",
						h31_nnz, max_index+1-bottom, top);
	H32.SetDataPtrCOO(&h32_row, &h32_col, &h32_data, "H32",
						h32_nnz, max_index+1-bottom, bottom-top);

	St1.MatrixMult(Hv11, H12);
	St2.MatrixMult(H21, St1); St1.Clear();

	S.CloneFrom(H22);
	S.MatrixAdd(St2, 1.0, -1.0, false); St2.Clear();

	Hv11.WriteFileMTX("prep/Hv11.mtx");
	H12.WriteFileMTX("prep/H12.mtx");
	H21.WriteFileMTX("prep/H21.mtx");
	H22.WriteFileMTX("prep/H22.mtx");
	H31.WriteFileMTX("prep/H31.mtx");
	H32.WriteFileMTX("prep/H32.mtx");
	S.WriteFileMTX("prep/S.mtx");
	printf(">>> Matrixes are successfully saved into file.\n");


	Hv11.Clear();
	H12.Clear();
	H21.Clear();
	H22.Clear();
	H31.Clear();
	H32.Clear();
	S.Clear();

	stop_paralution();


	free(off1);
	free(off2);

	free(in_degree);
	free(out_degree);
	//free(nf);
	free(nt);
	free(nf_idx);
	//free(nfR);
	free(ntR);
	free(nfR_idx);
	free(tmp_degree);
	free(tot_degree);
	free(CCcount);
	free(CClabel);
	free(CCrank);
	free(CCtag);
	free_queue(mQueue);

	free(reordered);
	free(new_index);
}

void loadEdges(FILE *inFile,
			   int *node_from, 
			   int *node_to, 
			   float *edge_value, 
			   int *in_degree,
			   int *out_degree,
			   int *max_index, 
			   size_t *edge_count){
	char * line = NULL;
	size_t len = 0;
	ssize_t read = 0;
	char *ptr, *uu, *vv, *val;
	int u, v;
	float value;

	*max_index = 0;
	*edge_count = 0;
	while (read = getline(&line, &len, inFile) != EOF)
	{
		ptr = line;
		uu = ptr;
		while(*ptr++ != DELIMITER);
		*(ptr-1) = '\0';
		vv = ptr;
		while(*ptr++ != DELIMITER);
		*(ptr-1) = '\0';
		val = ptr;

		u = atoi(uu);
		v = atoi(vv);
		value = atof(val);

		//printf("%d, %d, %f\n", u, v, value);


		if (*max_index < u) *max_index = u;
		if (*max_index < v) *max_index = v;

		*node_from++	= u;
		*node_to++		= v;
		*edge_value++	= value;

		out_degree[u]++;
		in_degree[v]++;

		(*edge_count)++;
	}
}

void* merge_sort_t(void *arg)
{
	arg_struct *a = (arg_struct*) arg;

	merge_sort(a->array, a->subarray, a->p, (a->q)-1, a->dir);
	return NULL;
}

void merge_sort_pthread(int *array, int *subarray, int top, int bottom, int dir){
	int i, j;
	int* p = (int*)malloc(sizeof(int)*nthreads);
	int* q = (int*)malloc(sizeof(int)*nthreads);
	arg_struct args[nthreads];

	pthread_t threads[nthreads];

	for (i = 0; i < nthreads; i++){
		p[i] = top + (bottom - top + 1) * i / nthreads;
		q[i] = top + (bottom - top + 1) * (i + 1) / nthreads;
		args[i].array = array;
		args[i].subarray = subarray;
		args[i].p = p[i];
		args[i].q = q[i];
		args[i].dir = dir;

		pthread_create(&threads[i], NULL, merge_sort_t, (void *)&args[i]);
	}


	for (i = 0; i < nthreads; i++){
		pthread_join(threads[i], NULL);
	}


	int* tmp    = (int*)malloc(sizeof(int)*(bottom-top+1));
	int* subtmp = (int*)malloc(sizeof(int)*(bottom-top+1));

	int marrval, msubval, mj;

	i = top;
	while(i <= bottom){
		if (dir < 0){
			marrval = INT_MIN;
		}
		else {
			marrval = INT_MAX;
		}
		msubval = -1;
		mj = -1;
		for(j = 0; j < nthreads; j++){
			if ( p[j] < q[j] ){
				if (dir < 0 && marrval < array[p[j]] ){
					marrval = array[p[j]];
					msubval = subarray[p[j]];
					mj = j;
				}
				if (dir > 0 && marrval > array[p[j]] ){
					marrval = array[p[j]];
					msubval = subarray[p[j]];
					mj = j;
				}
			}
		}

		tmp[i-top]    = marrval;
		subtmp[i-top] = msubval;
		p[mj]++;
		i++;
	}

	memcpy(array + top, tmp, sizeof(int)*(bottom-top+1));
	memcpy(subarray + top, subtmp, sizeof(int)*(bottom-top+1));

	free(tmp);
	free(subtmp);
}

void trp_push_back(size_t* nnz, size_t* nnz_max, int** row, int** col, double** data,
				   int i, int j, double v){
	if (*nnz == *nnz_max){
		*nnz_max = (*nnz_max) / 2 * 3;
		*row  = (int*) realloc (*row, sizeof(int) * (*nnz_max));
		*col  = (int*) realloc (*col, sizeof(int) * (*nnz_max));
		*data  = (double*) realloc (*data, sizeof(double) * (*nnz_max));
	}

	*((*row) + (*nnz)) = i;
	*((*col) + (*nnz)) = j;
	*((*data) + (*nnz)) = v;
	
	*nnz = *nnz + 1;
}

void* H_inv_t(void *arg)
{
	arg_inv *a = (arg_inv*) arg;
	size_t bidx, tidx;
	size_t block_count = a->block_count;
	size_t *block_idx = a->block_idx;
	size_t *h11_idx = a->h11_idx;
	double **H11 = a->H11;
	int *h11_row = a->h11_row;
	int *h11_col = a->h11_col;
	double *h11_data = a->h11_data;
	double val;

	size_t part = a->part;
	size_t i, j, k;
	for(bidx = part; bidx < block_count; bidx+=nthreads){
		//if ( (bidx+1) % ((block_count+1)/10) == 0)
		//	printf("[%d / %d]\n", bidx+1, block_count);
		LUinverse(H11[bidx], block_idx[bidx+1]-block_idx[bidx]);
		for (j = 0; j < block_idx[bidx+1] - block_idx[bidx]; j++){
			for (k = 0; k < block_idx[bidx+1] - block_idx[bidx]; k++){
				val =  H11[bidx][j * (block_idx[bidx+1] - block_idx[bidx]) + k];
				tidx = h11_idx[bidx] + j*(block_idx[bidx+1] - block_idx[bidx]) + k;
				//assert(tidx < h11_idx[bidx+1]);

				// bloating...
				//printf("|%zu", tidx);
				h11_row[tidx] = block_idx[bidx] + j;
				h11_col[tidx] = block_idx[bidx] + k;
				h11_data[tidx] = val;
				//*((*row) + (*nnz)) = i;
				//*((*col) + (*nnz)) = j;
				//*((*data) + (*nnz)) = v;
			}
		}
	}
	return NULL;
}


void LUinverse(double *A, size_t n){
	int i, j;

	mat M = (double**) malloc(sizeof(double*) * n);
	M[0] = A;
	for(size_t i = 1; i < n; i++) M[i] = M[0] + n * i;

	mat L, U, Lv, Uv, P, W1, W2;
	L = mat_new(n);
	U = mat_new(n);
	P = mat_new(n);

	mat_LU(M, L, U, P, n);
	Lv = mat_invL(L, n);
	Uv = mat_invU(U, n);
	W1 = mat_mul(Lv, P, n);
	W2 = mat_mul(Uv, W1, n);

	for (int i = 0; i < n; i++)
		for (int j= 0; j < n; j++)
			A[i*n + j] = W2[i][j];

	free(M);
	mat_del(L);
	mat_del(Lv);
	mat_del(U);
	mat_del(Uv);
	mat_del(P);
	mat_del(W1);
	mat_del(W2);
}
