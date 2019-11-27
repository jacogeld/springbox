/*==============================================================================
 * gcc -O3 -W -Wall -pedantic -c queue.c
 * gcc -O3 -W -Wall -pedantic -c worker.c
 * gcc -O3 -W -Wall -pedantic queue.o worker.o -o worker
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "queue.h"

#define DEBUG (0)
#define SEARCH_DEPTH (8)
#define ALPHA (4)

#define LAST ('/' + ALPHA)

int max = 0;
char maxw[10000];

int box_ok(char* w, int i) {
	int j = i - 1;
	while ((j >= 0) && (i - j <= ALPHA) && (w[i] != w[j])) {
		j--;
	}
	if ((j < 0) || (i - j > ALPHA)) {
		return 1;
	}
	for (int k = 0, l = i - j; l <= j; k++, l++) {
		int q = 0, n = i - j;
		for (; q <= n; q++) {
			if (w[k + q] != w[j + q]) {
				break;
			}
		}
		if (q > n) {
			return 0;
		}
	}
	return 1;
}

int limit = 10;

void process(char* w) {
	int k = strlen(w);
	int i = k;
	if (DEBUG) {
		printf("process(w==\"%s\")\n", w);
		printf("k==%d\n", k);
		printf("i==%d\n", i);
		if (limit-- < 0) { exit(1); }
	}
	w[i] = '/';
	while (i >= k) {
		if (DEBUG) {
			printf("exploring position i==%d w==", i);
			for (int ii = 0; ii <= i; ii++) {
				printf("%c", w[ii]);
			}
			printf("\n");
		}
		if (i == k + SEARCH_DEPTH) {
			if (DEBUG) {
				printf("too deep\n");
			}
			w[i] = '\0';
			enqueue(w);
			i--;
			continue;
		}
		w[i]++;
		if (w[i] > LAST) {
			if (DEBUG) {
				printf("done with this branch\n");
			}
			i--;
			continue;
		}
		if (box_ok(w, i)) {
			if (DEBUG) {
				printf("no repeat box\n");
			}
			i++;
			if (i > max) {
				max = i;
				w[i] = '\0';
				strcpy(maxw, w);
			}
			w[i] = '/';
		} else if (DEBUG) {
			printf("repeated box\n");
		}
	}
}

int main(int argc, char* argv[]) {
	if (argc != 1) {
		printf("Usage: %s\n", argv[0]);
	}

	char w[10000];
	enqueue("0123");
	int counter = 100;
	while (1) {
		if (dequeue(w)) {
			process(w);
		} else {
			break;
		}
		if (--counter < 0) {
			counter = 100;
			printf("(q%d) %d: %s\n", get_queue_max(), max, maxw);
		}
	}
	printf("max == %d\n", max);
	printf("maxw == %s\n", maxw);
	printf("queue_max == %d\n", get_queue_max());

}

