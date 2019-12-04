/******************************************************************************
 *	Priority queue using O(1) long design. Designed to be utilised by search.c
 *	compile: see search.c
 *
 *	@author C. R. Zeeman (caleb.zeeman@gmail.com)
 *	@version 1.4
 *	@date 2019-12-04
 *****************************************************************************/

/* Includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "search.h"

/* Definitions and macros */
#define DEBUG 	(0)
#define SUCCESS (1)
#define FAILURE (0)

/* function prototypes */
void enqueue_s(char *item);
int dequeue_s(char *item);
int get_queue_size();
void enqueue_l(unsigned long long *arr, int count);
int dequeue_l(unsigned long long *arr, int *count);

typedef struct NodeStruct {
	char *item;
	int size;
	struct NodeStruct *next;
} Node;

typedef struct NodeStructL {
	unsigned long long *item;
	int size;
	struct NodeStructL *next;
} NodeL;

/* globals */
Node *queue_head = NULL;
Node *queue_tail = NULL;
int queue_len = 0;

/* Temporarily having two queues, one for V1 and one for V2.
	TODO replace string-based with long-based */
NodeL *queue_head_l = NULL;
NodeL *queue_tail_l = NULL;
int queue_len_l = 0;

/* functions */
int get_queue_size() {
	return queue_len;
}

void enqueue_s(char *w) {
	Node *new_node = (Node*) malloc(sizeof(Node));
	Node *x = queue_head;
	new_node->item = strdup(w);
	new_node->size = strlen(w);
	new_node->next = NULL;
	if (queue_head == NULL) { /* Empty queue */
		queue_head = new_node;
		queue_tail = new_node;
	}
	else if (queue_head->next == NULL) { /* Only one item currently in queue */
		if (new_node->size < queue_head->size) {
			queue_head->next = new_node;
			queue_tail = new_node;
		}
		else {
			new_node->next = queue_head;
			queue_head = new_node;
		}
	}	/* Multiple items already in queue */
	else {
		if (new_node->size >= x->size) {
			new_node->next = x;
			queue_head = new_node;
		} 
		else {
			while (new_node->size < x->next->size) {
				if (x->next == NULL) { /* If at end */
					x->next = new_node;
					queue_tail = new_node;
					goto end;
				}
				x = x->next;
			}
			new_node->next = x->next;
			x->next = new_node;
		}
	}
	end:
	queue_len++;
	if (DEBUG)  {
		printf("<<<<<<enqueued: %s>>>>>>>\n", new_node->item);
	}
}

int dequeue_s(char *w) {
	Node *old_head;
	if (queue_head == NULL) {
		return FAILURE;
	}
	old_head = queue_head;
	strcpy(w, queue_head->item);
	queue_head = queue_head->next;
	/*free(old_head->item); TODO QUERY: Uncomment? Once all working */ 
	free(old_head); 
	if (queue_head == NULL) {
		queue_tail = NULL;
	}
	queue_len--;
	return SUCCESS;
}

/*---------------------------------------------------------------*/
int get_queue_size_l() {
	return queue_len_l;
}
void enqueue_l(unsigned long long *arr, int count) {
	NodeL *new_node = (NodeL*) malloc(sizeof(NodeL));
	NodeL *x = queue_head_l;
	int i, size = (int) (count / chars_per_long);
	if (count % chars_per_long != 0) {
		size++;	/* Extra long if needed */
	}
	new_node->item = malloc(sizeof(long long) * size);
	for (i = 0; i < size; i++) { /* Copy longs */
		new_node->item[i] = arr[i];
	}
	new_node->size = count;
	new_node->next = NULL;
	if (queue_head_l == NULL) { /* Empty queue */
		queue_head_l = new_node;
		queue_tail_l = new_node;
	}
	else if (queue_head_l->next == NULL) { /* Only one item currently in queue */
		if (new_node->size < queue_head_l->size) {
			queue_head_l->next = new_node;
			queue_tail_l = new_node;
		}
		else {
			new_node->next = queue_head_l;
			queue_head_l = new_node;
		}
	}	/* Multiple items already in queue */
	else {
		if (new_node->size >= x->size) {
			new_node->next = x;
			queue_head_l = new_node;
		} 
		else {
			while (new_node->size < x->next->size) {
				if (x->next == NULL) { /* If at end */
					x->next = new_node;
					queue_tail_l = new_node;
					goto end;
				}
				x = x->next;
			}
			new_node->next = x->next;
			x->next = new_node;
		}
	}
	end:
	queue_len_l++;
}
/* User must free arr before use. Count is set to length */
int dequeue_l(unsigned long long *arr, int* count) {
	NodeL *old_head;
	int i, size;
	*count = queue_head_l->size;
	size = (int) (*count / chars_per_long);
	if (*count % chars_per_long != 0) {
		size++;	/* Extra long if needed */
	}
	arr = malloc(sizeof(long long) * size);

	for (i = 0; i < *count; i++) { /* Copy longs */
		arr[i] = queue_head_l->item[i];
	}
	if (queue_head == NULL) {
		return FAILURE;
	}
	old_head = queue_head_l;
	queue_head_l = queue_head_l->next;
	free(old_head->item);
	free(old_head); 
	if (queue_head_l == NULL) {
		queue_tail_l = NULL;
	}
	queue_len_l--;
	return SUCCESS;
}

/* TODO
void enqueue_lm(unsigned long *arr, int count, int my_id, int root_id) {
	return;
}

int dequeue_lm(unsigned long *arr, int count, int my_id, int root_id) {
	return 0;
}*/