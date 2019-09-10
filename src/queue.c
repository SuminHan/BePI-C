#include <stdio.h>
#include <stdlib.h>
#include "queue.h"

queue* init_queue(int init_size){
	queue* q = (queue*) malloc(sizeof(queue));
	q->MAX = init_size;
	q->front = q->rear = 0;
	q->array = (int*) malloc(sizeof(int)*init_size);
	return q;
}

void clear_queue(queue* q){
    q->front = q->rear;
}

void free_queue(queue* q){
	free(q->array);
	free(q);
}

int queue_put(queue* q, int v){
    if ((q->rear + 1) % q->MAX == q->front){
        return -1;
    }

    q->array[q->rear] = v;
    q->rear = ++q->rear % q->MAX;
    return v;
}

int queue_get(queue* q){
    int i;
    if (q->front == q->rear){
        return -1;
    }

    i = q->array[q->front];
    q->front = ++q->front % q->MAX;
    return i;
}

int is_queue_empty(queue* q){
	return q->front == q->rear;
}
