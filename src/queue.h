#ifndef QUEUE_H
#define QUEUE_H
/* TODO: Document functions */

/* Enqueue and dequeue of longs (V2) */
/* Count - Number of 'characters' in long (length of long) */
void enqueue_l(unsigned long long *arr, int count);
int dequeue_l(unsigned long long *arr, int *count);
long long get_queue_size_l();

/* Helper function to store queue to file without dequeue */
void store_queue_helper(FILE *file, int max_longs_needed);

/* Enqueue and dequeue of strings (V1) */
void enqueue_s(char *item);
int dequeue_s(char *item);
long long get_queue_size();

#endif