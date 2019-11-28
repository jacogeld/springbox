#ifndef QUEUE_H
#define QUEUE_H

void enqueue(char *item);
int dequeue(char *item);
int get_queue_size();
#endif

/* TODO: Create code to enqueue and dequeue across processes to reduce amounts
 * of strcpy and strdup needed to execute. Special case if process ==
 * ROOT_PROCESS (include as param?)
void enqueue_d(int process_id, int length);
int dequeue_d(int process_id, int length);
*/
