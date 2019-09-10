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
void loadEdges(FILE*, int*, int*, int*, int*, int*, size_t*);
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

void loadQueryList(FILE *inFile, int *query_count, int *query_list){
	char * line = NULL;
	size_t len = 0;
	ssize_t read = 0;
	char *ptr, *uu;
	int u;

	while (read = getline(&line, &len, inFile) != EOF)
	{
		ptr = line;
		uu = ptr;

		u = atoi(uu);

		(*query_count)++;
		*query_list++ = u;
	}
}

void cogi_query(char *query_file){
	int i;
	long long starttime, endtime;
	starttime = current_timestamp();

	FILE *inFile = fopen(query_file, "r");
	if (inFile == NULL){
		fprintf(stderr, "File %s not found!\n", query_file);
		return;
	}

	int query_count;
	int *query_list = (int*) malloc (sizeof(int) * MAX_NODES);
	
	loadQueryList(inFile, &query_count, query_list);

	fprintf(stdout, "[QueryList]:", i, query_list[i]);
	for(i = 0; i < query_count; i++){
		fprintf(stdout, " %d", query_list[i]);
	}
	fprintf(stdout, "\n");

	fprintf(stdout, ">>> Loading Matrix...\n");

	init_paralution();
	info_paralution();

	LocalMatrix<double> Hv11, H12, H21, H22, H31, H32, S;
	Hv11.ReadFileMTX("prep/Hv11.mtx");
	 H12.ReadFileMTX("prep/H12.mtx");
	 H21.ReadFileMTX("prep/H21.mtx");
	 H22.ReadFileMTX("prep/H22.mtx");
	 H31.ReadFileMTX("prep/H31.mtx");
	 H32.ReadFileMTX("prep/H32.mtx");
	   S.ReadFileMTX("prep/S.mtx");

	LocalVector<int> new_index, reordered;
	new_index.ReadFileASCII("prep/new_index.ASCII");
	reordered.ReadFileASCII("prep/reordered.ASCII");

	fprintf(stdout, ">>> Loading Completed!\n");

	int top = Hv11.get_nrow(),
		bottom = Hv11.get_nrow() + H21.get_nrow(),
		max_index = Hv11.get_nrow() + H21.get_nrow() + H31.get_nrow() -1;

	LocalVector<double> q1_vec, q2_vec, q3_vec, r1_vec, r2_vec, r3_vec;
	q1_vec.Allocate("q1 vector", Hv11.get_nrow());
	q2_vec.Allocate("q2 vector", H22.get_nrow());
	q3_vec.Allocate("q3 vector", H31.get_nrow());
	q1_vec.Zeros();
	q2_vec.Zeros();
	q3_vec.Zeros();
	r1_vec.Allocate("r1 vector", q1_vec.get_size());
	r2_vec.Allocate("r2 vector", q2_vec.get_size());
	r3_vec.Allocate("r3 vector", q3_vec.get_size());

	LocalVector<double> T1, T2, T3, t2_vec;
	T1.Allocate("T1 vector", q1_vec.get_size());
	T2.Allocate("T2 vector", q2_vec.get_size());
	t2_vec.Allocate("t2 vector", q2_vec.get_size());

	GMRES<LocalMatrix<double>, LocalVector<double>, double> ls;
	ls.SetOperator(S);
	ls.Build();

	LocalVector<double> H12xr2, H31xr1, H32xr2;
	H12xr2.Allocate("H12xr2", top);
	H31xr1.Allocate("H31xr1", max_index+1-bottom);
	H32xr2.Allocate("H32xr2", max_index+1-bottom);

	for (i = 0; i < query_count; i++){
		int qnew_index = new_index[query_list[i]]; //new_index[query_number];
		if (qnew_index < top){
			q1_vec[qnew_index] += RESTART/query_count;
		}
		else if (qnew_index < bottom){
			q2_vec[qnew_index - top] += RESTART/query_count;
		}
		else {
			q3_vec[qnew_index - bottom] += RESTART/query_count;
		}
	}

	Hv11.Apply(q1_vec, &T1);
	H21.Apply(T1, &T2); T1.Clear();
	t2_vec.ScaleAdd2(0, q2_vec, 1, T2, -1.0); T2.Clear();

	ls.Solve(t2_vec, &r2_vec);

	H12.Apply(r2_vec, &H12xr2);
	T3.CopyFrom(H12xr2);
	T3.ScaleAdd(-1.0, q1_vec);
	Hv11.Apply(T3, &r1_vec); T3.Clear();

	H31.Apply(r1_vec, &H31xr1);
	H32.Apply(r2_vec, &H32xr2);

	r3_vec.CopyFrom(q3_vec);
	r3_vec.ScaleAdd2(1.0, H31xr1, -1.0, H32xr2, -1.0);
	

	FILE *rfp = fopen("result.txt", "w");
	float result_sum = 0.0;
	int idx, iidx;
	double *results = (double*) malloc (sizeof(double) * (max_index + 1));
	for (i = 0; i <= max_index; i++){
		idx = new_index[i];
		if (idx < top){
			iidx = idx;
			//fprintf(rfp, "%f\n", r1_vec[iidx]);
			result_sum += r1_vec[iidx];
			results[i] = r1_vec[iidx];
		}
		else if (idx < bottom){
			iidx = idx - top;
			//fprintf(rfp, "%f\n", r2_vec[iidx]);
			result_sum += r2_vec[iidx];
			results[i] = r2_vec[iidx];
		}
		else {
			iidx = idx - bottom;
			//fprintf(rfp, "%f\n", r3_vec[iidx]);
			result_sum += r3_vec[iidx];
			results[i] = r3_vec[iidx];
		}
	}
	for (i = 0; i <= max_index; i++){
		fprintf(rfp, "%e\n", results[i] / result_sum);
	}

	fprintf(stdout, "Result sum = %f\n", result_sum);
	//fprintf(rfp, "Result sum = %f\n", result_sum);
	fclose(rfp);

	endtime = current_timestamp();

	fprintf(stdout, "Query: %.3f sec\n", (endtime-starttime)/1000.0);

	stop_paralution();
	free(query_list);
}

int main(int argc, char* argv[])
{
	if (argc != 2)
	{
		fprintf(stderr, "Usage: ./cogi_server <input.txt>\n");
		return -1;
	}

	cogi_query(argv[1]);
	fprintf(stdout, "Process terminated successfully.\n");

	return 0;
}

