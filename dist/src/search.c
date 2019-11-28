/******************************************************************************
 *	TODO: Describe code
 *	compile:	mpicc -o ../bin/search.o queue.c search.c
 *	run:		mpirun -np X ../search.o		where X is number processes
 *
 *	@author C. R. Zeeman (caleb.zeeman@gmail.com)
 *	@version 0.1
 *	@date 2019-11-28
 *****************************************************************************/
/* Includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "queue.h"

/************************** Definitions and macros **************************/

#define DEBUG (0)
#define STANDARD_SEARCH_DEPTH (4)
#define INITIAL_SEARCH_DEPTH (3)
#define STANDARD_ALPHABET_SIZE (3)

#define LAST ('/' + STANDARD_ALPHABET_SIZE)
#define VALID (1)
#define INVALID (0)

/* Constants for MPI signalling */
#define ROOT_PROCESS		(0)
#define STOP_TAG			(2000)
#define SEND_MAX_TAG		(2001)
#define NEW_WORK_TAG		(2002)
#define WORK_TAG			(2003)
#define SENDING_MAX_TAG		(2004)
#define MAX_STRING_TAG		(2005)
#define EXPECT_QUEUE_TAG	(2006)
#define QUEUE_TAG			(2007)

/*************************** global variables ******************************/

/* Alphabet used -> alpha_size subset of alpha[] */
/* TODO: Allow dynamic setting via cmd line */

/* TODO: Documentation */
/* TODO: Changelog. Leave todos in after working V1, then can do versions
 * afterwards (e.g. 1.1, 1.2) */

 /* BUG: Occasionally an invalid case gets through process, even though the
  * check tools pick it up as invalid if manually checking */
int alpha_size = STANDARD_ALPHABET_SIZE;
char alpha[] = {'0','1','2','3','4','5', '6', '7', '8', '9'};

/* Theoretical longest word per alphabet size. l(n) = n(l(n-1) +2) */
int max_word_size[] = {0, 2, 8, 30, 128, 650, 3912, 27398, 219200, 1972818};

int max_length = 0;
char *max_word;
int my_id;		/* Process ID; set in main */
/*********************** function prototypes *******************************/

/* TODO: rename box_ok to box_repeat_exists */
int box_repeat_exists(char *word, int len);	/* TODO: Convert to inline*/
void process (char *word, int depth);		/* TODO: convert to inline*/
int box_ok(char *word, int len); 	/* TODO: remove and use 1 */
int box_ok_1(char *word, int len);
void process_1(char *word, int depth);	/* "" */

void manual_check(char *string, int len);	/*Check if pattern is box-valid */

/**************************** functions ***********************************/

/* TODO: Make all cases of 'enqueue' send the signal and data to root process if
 * it is not the root process so it can be enqueued. Or merge this into the
 * actual enqueue function (so it recieves in that function) */

int main(int argc, char *argv[])
{
	/* Setup */
	int ierr = 0, nbr_procs = 0, rec_count, send_count, an_id, i = 0,
	queue_size = 0, desired_size, stop = 0;
	MPI_Status status;
	char *w; 
	/* Keeps track of how many processes have been successfully stopped */
	int active_processes;

	/*TODO: Convert desired_size and malloc to handle alpha size from cmd */
	desired_size = max_word_size[alpha_size];
	max_word = (char*) malloc(sizeof(char) * (desired_size + 1));  /* +1 for '\0' */
	max_word[0] = '\0';

	
	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &nbr_procs);

	if (nbr_procs <= 1) {
		printf("Usage: \"mpirun -np X %s\" where X is number processess > 1\n",
			argv[0]);
		goto end;
	}

	active_processes = nbr_procs - 1; /* Not including root */

	w = (char*) malloc(sizeof(char) * (max_word_size[alpha_size] + 1));

	/* Root/Master */
	if (my_id == ROOT_PROCESS) {
		/*DELME*/
		/*
		char str[] = "012001210122010201120210212012";
		manual_check(str, sizeof(str)/sizeof(char));
		char str2[] = "0220102012101210012011021202102";
		manual_check(str2, sizeof(str2)/sizeof(char));
		*/
		/* Initial run to create branches for workers */
		w[0] = '0';
		w[1] = '\0';
		process(w, INITIAL_SEARCH_DEPTH);

		MPI_Barrier(MPI_COMM_WORLD);

		/* TODO: Send initial instructions to each process */
		for (an_id = 1; an_id < nbr_procs; an_id++) {
			if (!dequeue(w)) {
				break;	/* Failsafe */	
			}
			send_count = strlen(w);
			if (DEBUG) {
				printf("(1)Dequeued %s with length %d\n", w, send_count);
			}
			ierr = MPI_Send(&send_count, 1, MPI_INT, an_id, NEW_WORK_TAG,
				MPI_COMM_WORLD);
			ierr = MPI_Send(w, send_count, MPI_CHAR, an_id, WORK_TAG, 
				MPI_COMM_WORLD);
		}

		/* Communcation cycle */
		/* Wait for size found from process; if size > max, ask process to send
		 * string. Else, if need to stop, send stop signal. else send next
		 * workload signal followed by workload */
		queue_size = get_queue_size();
		while (1) {
			if (queue_size <= 0 ||max_length >= desired_size) {
				stop = 1;
			}
			/* Recieve a reply */
			ierr = MPI_Recv(&rec_count, 1, MPI_INT, MPI_ANY_SOURCE,
				MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if (DEBUG) {
				printf("Main process recieved a command: %d\n", status.MPI_TAG);
			}
			an_id = status.MPI_SOURCE;

			/* Recieving max found from process */
			if (status.MPI_TAG == SENDING_MAX_TAG) {
				if (rec_count > max_length) {
					max_length = rec_count;
					if (DEBUG) printf("(n)New max found, requesting string\n");
					/* Request longest string*/
					MPI_Send(&rec_count, 1, MPI_INT, an_id, SEND_MAX_TAG,
						MPI_COMM_WORLD);
					MPI_Recv(max_word, max_length, MPI_CHAR, an_id,
						MAX_STRING_TAG, MPI_COMM_WORLD, &status);
				}	
				/* Send stop command if needed... */
				if (stop) {
					/* TODO: Keep track of which processes stopped and if all
					 * stopped, then main can also end */
					 MPI_Send(&send_count, 1, MPI_INT, an_id, STOP_TAG,
					 	MPI_COMM_WORLD);
					active_processes--;
					if (active_processes == 0) {
						break;
					}
				} /* or send the next work available in queue */
				else if (!dequeue(w)) {
					break;	/* Failsafe */
				}
				else {
					send_count = strlen(w);
					if (DEBUG) {
						printf("(n)Dequeued %s with length %d\n", w, send_count);
					}
					ierr = MPI_Send(&send_count, 1, MPI_INT, an_id,
						NEW_WORK_TAG, MPI_COMM_WORLD);
					ierr = MPI_Send(w, send_count, MPI_CHAR, an_id, WORK_TAG, 
						MPI_COMM_WORLD);
				}
			}
			/* Need to queue a string */
			else if (status.MPI_TAG == EXPECT_QUEUE_TAG) {
			/* TODO: Move to new enqueue code 
				enqueue_d(an_id, rec_count); */
				MPI_Recv(w, rec_count, MPI_CHAR, an_id, QUEUE_TAG,
					MPI_COMM_WORLD, &status);
				if (DEBUG) {
					printf("(n)About to enqueue %s of size %d\n", w, rec_count);
				}
				enqueue(w);
			}
			/*
			...
			update max size by getting result from workers */
			/* TODO: two way communication: process sends max length back.
			 * Then main sends command signal: 1 int stating either "sending
			 * next workload", or "send me the max string back", or "stop" 
			 */
			/*...*/
			queue_size = get_queue_size();
		} 
		/* qs == 0 => explored all || ml == ds => found a longest pattern */
		/* Send stop signal to all processes */

		/* print status */
		printf("max == %d\n", max_length);
		printf("max word == %s\n", max_word);
		printf("queue size == %d\n", get_queue_size());
	}
	
	/* Worker/Slave */
	else {
		MPI_Barrier(MPI_COMM_WORLD);

		/* TODO: exception if count higher than array size */

		/* Continue until a message is recieved stating you must stop */
		while (status.MPI_TAG != STOP_TAG) {
			/* Get next command */
			ierr = MPI_Recv(&rec_count, 1, MPI_INT, ROOT_PROCESS,
				MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			/*if command is to send max string back, do that */
			if (status.MPI_TAG == SEND_MAX_TAG) {
				ierr = MPI_Send(max_word, max_length, MPI_CHAR,
					ROOT_PROCESS, MAX_STRING_TAG, MPI_COMM_WORLD);
			}

			 /* else if command is to get new work, get new work string */
			else if (status.MPI_TAG == NEW_WORK_TAG) {
				if (DEBUG) printf("*:rec_count == %d\n", rec_count);
				ierr = MPI_Recv(w, rec_count, MPI_CHAR, ROOT_PROCESS,
					WORK_TAG, MPI_COMM_WORLD, &status);

				w[rec_count] = '\0';
				if (DEBUG) {
					printf("process %d recieved %s with %d\n",
						my_id, w, rec_count);
				}
				process(w, STANDARD_SEARCH_DEPTH);		

				/* Send local max length back */
				ierr = MPI_Send(&max_length, 1, MPI_INT, ROOT_PROCESS,
					SENDING_MAX_TAG, MPI_COMM_WORLD);
			}
		}
	}
	
	end: ierr = MPI_Finalize();
	printf("***** PROCESS %d STOPPED ******\n", my_id);
	return EXIT_SUCCESS;
}

void process_1(char *word, int depth) {
	int k = strlen(word), i = k, valid = 1, char_sel = 0;

	while (i >= k) {
		break;
		/* Insert */
		/* Check if valid */
		/* if valid, increment i */
			/* if i = k+depth, queue and i-- */
		/* increment w[i] */
		/* if w[i] > '0' +  alpha_size - 1*/
			/* i-- */
	}
}

/* Jaco's process code. TODO: Refactor */
int limit = 10000;

void process(char *word, int depth) {
	int k = strlen(word);
	int i = k;
	int valid = 1;
	if (DEBUG) {
		printf("process(w==\"%s\")\n", word);
		printf("k==%d\n", k);
		printf("i==%d\n", i);
		if (limit-- < 0) { exit(1); }
	}
	word[i] = '/';
	while (i >= k) {
		if (DEBUG) {
			printf("exploring position i==%d w==", i);
			for (int ii = 0; ii <= i; ii++) {
				printf("%c", word[ii]);
			}
			printf("\n");
		}
		if (i == k + depth) {
			if (DEBUG) {
				printf("too deep\n");
			}
			word[i] = '\0';
			if (valid) {
			if (my_id != ROOT_PROCESS) {
			/* TODO: Move to new enqueue_d function*/
				MPI_Send(&i, 1, MPI_INT, ROOT_PROCESS, EXPECT_QUEUE_TAG,
					MPI_COMM_WORLD);
				MPI_Send(word, i, MPI_CHAR, ROOT_PROCESS, QUEUE_TAG,
					MPI_COMM_WORLD);
			} else {
				if (DEBUG) printf("(1)about to enqueue %s\n", word);
				enqueue(word);
			}
			}
			i--;
			continue;
		}
		word[i]++;
		if (word[i] > LAST) {
			if (DEBUG) {
				printf("done with this branch\n");
			}
			word[i] = '\0'; /* Temp */
			i--;
			continue;
		}
		valid = box_ok_1(word,i);
		if (valid) {
			if (DEBUG) {
				printf("no repeat box\n");
			}
			i++;
			if (i > max_length) {
				max_length = i;
				word[i] = '\0';
				strcpy(max_word, word);
			}
			word[i] = '/';
		} else if (DEBUG) {
			printf("repeated box\n");
		}
	}
	if (DEBUG) printf("Done process\n");
}

/* Jaco's Box_ok*/
int box_ok(char* w, int i) {
	int j = i - 1;
	while ((j >= 0) && (i - j <= STANDARD_ALPHABET_SIZE) && (w[i] != w[j])) {
		j--;
	}
	if ((j < 0) || (i - j > STANDARD_ALPHABET_SIZE)) {
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


/* TODO: Rename & Refactor */
int box_ok_1(char* pattern, int len)
{
	char c, c_end;
	int count = 1, b1_end = len, b1_start = b1_end, b2_start, b2_end, i = 0;
	int lb1, lb2;
	/* Step back until find char at pattern[len-1]
	 * Remember this box. Find 3rd repeat, check box. find 4th repeat. Check 3-4
	 * box and compare to 1-2 box. Repeat all boxes 
	 * optimisation: skip if boxes too long (> alpha_size) TODO
	*/
	/* Find box */
	c_end = pattern[b1_start--];
	while (b1_start >= 0) {
		c = pattern[b1_start];
		if (c == c_end) {
			count++;
			break;
		}
		b1_start--;
	}
	/* No box */
	if (b1_start < 0) {
		return VALID;
	}

	/* Run through blocks with same start/end and compare */
	b2_end = b1_start;
	b2_start = b1_start-1;
	while (b2_start >= 0) {
		c = pattern[b2_start];
		if (c == c_end) {
			lb1 = b1_end - b1_start + 1;
			lb2 = b2_end - b2_start + 1;
			/* Compare length */
			if (lb1 != lb2) {
				b2_end = b2_start;
				goto end_iteration;
			}
			count++;
			/* Check for odd repeats (i.e. common edges) */
			if (count % 2 == 1) {
				/* Compare characters */
				for (i = 0; i < lb2; i++) {
					if (pattern[b2_start + i] != pattern[b1_start + i]) {
						b2_end = b2_start;
						goto end_iteration;
					}
				}
				return INVALID;
				
			}

			/* Check for even repeats (i.e. new box) */
			else {
				/* Compare characters */
				for (i = 0; i < lb2; i++) {
					if (pattern[b2_start + i] != pattern[b1_start + i]) {
						b2_end = b2_start;
						goto end_iteration;
					}
				}
				return INVALID;
			}
		}
		end_iteration: b2_start--;
	}
	return VALID;
}

void manual_check(char *str, int len) {
	char *c = str;
	char *w = malloc(sizeof(char) * (len + 1));
	int i;
	printf("local test\n");
	for (i = 0; i < len; i++) {
		w[i] = c[i];
		w[i+1] = '\0';
		printf("i:%d\tValid:%d\n", i, box_ok(w, i));
	}
	free(w);
}
