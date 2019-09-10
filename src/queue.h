typedef struct{
	int MAX;
	int* array;
	int front;
	int rear;
} queue;

queue* init_queue(int init_size);
void clear_queue(queue* q);
void free_queue(queue* q);
int queue_put(queue* q, int v);
int queue_get(queue* q);
int is_queue_empty(queue* q);
