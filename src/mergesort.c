#include <memory.h>
#include <stdlib.h>

/* MergeSort - Reference: http://milvus.tistory.com/69 */
void merge(int* array, int* subarray, int start, int mid, int end, int dir){
	int* tmp=(int*)malloc(sizeof(int)*(end-start+1));
	int* subtmp=(int*)malloc(sizeof(int)*(end-start+1));
	int tmp_index=0;
	int p=start,q=mid+1;
	int i;

	for(i=tmp_index; i<=end-start; i++){
		while(p<=mid && q<=end){
			if(dir>=0){
				if(array[p]>array[q]){
					tmp[tmp_index]=array[q];
					subtmp[tmp_index]=subarray[q];
					q++;
				}else{
					tmp[tmp_index]=array[p];
					subtmp[tmp_index]=subarray[p];
					p++;
				}
			}
			else{ //if(dir<0){
				if(array[p]<array[q]){
					tmp[tmp_index]=array[q];
					subtmp[tmp_index]=subarray[q];
					q++;
				}else{
					tmp[tmp_index]=array[p];
					subtmp[tmp_index]=subarray[p];
					p++;
				}
			}
			tmp_index++;
		}

		if(p>mid){
			while(q<=end){
				tmp[tmp_index]=array[q];
				subtmp[tmp_index]=subarray[q];
				q++;
				tmp_index++;
			}
		}
		else{
			while(p<=mid){
				tmp[tmp_index]=array[p];
				subtmp[tmp_index]=subarray[p];
				p++;
				tmp_index++;
			}
		}
	}//end for

	for(i=start; i<=end; i++){
		array[i]=tmp[i-start];
		subarray[i]=subtmp[i-start];
	}
	/*
	memcpy(array+start, tmp, sizeof(int)*(end-start+1));
	memcpy(subarray+start, subtmp, sizeof(int)*(end-start+1));
	*/

	free(tmp);
	free(subtmp);
}

void merge_sort(int* array, int* subarray, int start, int end, int dir){
	if(start>=end) return;

	int mid=(start+end)/2;

	merge_sort(array, subarray, start, mid, dir);
	merge_sort(array, subarray, mid+1, end, dir);

	merge(array, subarray, start, mid, end, dir);
}

/* MergeSort3 - Reference: http://milvus.tistory.com/69 */
void merge3(int* array, int* subarray, double* subarray2, size_t start, size_t mid, size_t end, int dir){
	int* tmp=(int*)malloc(sizeof(int)*(end-start+1));
	int* subtmp=(int*)malloc(sizeof(int)*(end-start+1));
	int* subtmp2=(int*)malloc(sizeof(int)*(end-start+1));
	int tmp_index=0;
	int p=start,q=mid+1;
	int i;

	for(i=tmp_index; i<=end-start; i++){
		while(p<=mid && q<=end){
			if(dir>=0){
				if(array[p]>array[q]){
					tmp[tmp_index]=array[q];
					subtmp[tmp_index]=subarray[q];
					subtmp2[tmp_index]=subarray2[q];
					q++;
				}else{
					tmp[tmp_index]=array[p];
					subtmp[tmp_index]=subarray[p];
					subtmp2[tmp_index]=subarray2[p];
					p++;
				}
			}
			else{ //if(dir<0){
				if(array[p]<array[q]){
					tmp[tmp_index]=array[q];
					subtmp[tmp_index]=subarray[q];
					subtmp2[tmp_index]=subarray2[q];
					q++;
				}else{
					tmp[tmp_index]=array[p];
					subtmp[tmp_index]=subarray[p];
					subtmp2[tmp_index]=subarray2[p];
					p++;
				}
			}
			tmp_index++;
		}

		if(p>mid){
			while(q<=end){
				tmp[tmp_index]=array[q];
				subtmp[tmp_index]=subarray[q];
				subtmp2[tmp_index]=subarray2[q];
				q++;
				tmp_index++;
			}
		}
		else{
			while(p<=mid){
				tmp[tmp_index]=array[p];
				subtmp[tmp_index]=subarray[p];
				subtmp2[tmp_index]=subarray2[p];
				p++;
				tmp_index++;
			}
		}
	}//end for

	for(i=start; i<=end; i++){
		array[i]=tmp[i-start];
		subarray[i]=subtmp[i-start];
		subarray2[i]=subtmp2[i-start];
	}
	/*
	memcpy(array+start, tmp, sizeof(int)*(end-start+1));
	memcpy(subarray+start, subtmp, sizeof(int)*(end-start+1));
	memcpy(subarray2+start, subtmp2, sizeof(int)*(end-start+1));
	*/

	free(tmp);
	free(subtmp);
	free(subtmp2);
}

void merge_sort3(int* array, int* subarray, double* subarray2, size_t start, size_t end, int dir){
	if(start>=end) return;

	int mid=(start+end)/2;

	merge_sort3(array, subarray, subarray2, start, mid, dir);
	merge_sort3(array, subarray, subarray2, mid+1, end, dir);

	merge3(array, subarray, subarray2, start, mid, end, dir);
}
