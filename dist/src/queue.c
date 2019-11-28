/******************************************************************************
 *  TODO: Describe code
 *  TODO: Compile & run commands
 *****************************************************************************/
 /* TODO: Convert from char *w to Work *w and create Work struct*/
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
void enqueue(char *item);
int dequeue(char *item);
int get_queue_size();

typedef struct NodeStruct {
	char *item;
	struct NodeStruct *next;
} Node;

/* globals */
Node *queue_head = NULL;
Node *queue_tail = NULL;
int queue_len = 0;

/* functions */
int get_queue_size() {
	return queue_len;
}

void enqueue(char *w) {
	Node *new_node = (Node*) malloc(sizeof(Node));
	new_node->item = strdup(w);
	new_node->next = NULL;
	if (queue_head == NULL) {
		queue_head = new_node;
	}
	else {
		queue_tail->next = new_node;
	}
	queue_tail = new_node;
	queue_len++;
	if (DEBUG)  {
		printf("<<<<<<enqueued: %s>>>>>>>\n", queue_tail->item);
	}
}

int dequeue(char *w) {
	Node *old_head;
	if (queue_head == NULL) {
		return FAILURE;
	}
	old_head = queue_head;
	strcpy(w, queue_head->item);
	queue_head = queue_head->next;
	free(old_head); 
	if (queue_head == NULL) {
		queue_tail = NULL;
	}
	queue_len--;
	return SUCCESS;
	
}
