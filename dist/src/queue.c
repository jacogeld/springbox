/******************************************************************************
 *	Standard priority queue. Designed to be utilised by search.c
 *	compile:	gcc -Wall -W -pedantic -o queue.o queue.c
 *			 OR see search.c
 *
 *	@author C. R. Zeeman (caleb.zeeman@gmail.com)
 *	@version 1.1
 *	@date 2019-11-29
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
	int size;
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