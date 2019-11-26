/*==============================================================================
 * gcc -O3 -W -Wall -pedantic -c queue.c
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "queue.h"

#define DEBUG (0)

typedef struct NodeStruct {
	char* item;
	struct NodeStruct* pred;
} Node;

Node* queue_head = NULL;
Node* queue_tail = NULL;
int queue_max = 0;
int queue_len = 0;

void enqueue(char* work) {
	if (DEBUG) {
		printf("ENQUEUE: \"%s\"\n", work);
	}
	Node* new_node = (Node*) malloc(sizeof(Node));
	new_node->item = strdup(work);
	new_node->pred = NULL;
	if (queue_head == NULL) {
		queue_head = new_node;
	} else {
		queue_tail->pred = new_node;
	}
	queue_tail = new_node;
	queue_len++;
	if (queue_len > queue_max) {
		queue_max = queue_len;
	}
}

int dequeue(char* work) {
	if (queue_head == NULL) {
		if (DEBUG) {
			printf("EMPTY QUEUE!\n");
		}
		return 0;
	}
	strcpy(work, queue_head->item);
	if (DEBUG) {
		printf("DEQUEUE: \"%s\"\n", work);
	}
	Node* old_head = queue_head;
	free(old_head->item);
	queue_head = old_head->pred;
	free(old_head);
	if (queue_head == NULL) {
		queue_tail = NULL;
	}
	queue_len--;
	return 1;
}

int get_queue_max() {
	return queue_max;
}

