/******************************************************************************
 *	Finds the longest string that has no repeated boxes. Utilises MPI to run
 *	multiple workers at once. Relies on a queue-based system to handle jobs
 *	between workers.
 *
 *  original process and box_ok code provided by Jaco Geldenhuys
 * 
 *	compile:	mpicc -o ../bin/search.o queue.c search.c
 *	run:		mpirun -np X ./search.o	N
 		where X is number processes,
	 	N alphabet size between 2 and 9 (default if not specified)
 *
 *	@author C. R. Zeeman (caleb.zeeman@gmail.com)
 *	@version 1.2
 *	@date 2019-11-29
 *****************************************************************************/
/* Includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "queue.h"

/************************** Definitions and macros **************************/
/* Configurations */
#define DEBUG (0)						/* Enables debug output */
#define STATUS (1)						/* Prints status messages */
#define PRINT_MAX (1)					/* Prints max whenever it is updated */
#define STANDARD_SEARCH_DEPTH (5)		/* Depth to search per job */
#define INITIAL_SEARCH_DEPTH (2)		/* Depth for initial run by ROOT */
#define STANDARD_ALPHABET_SIZE (3)		/* Default alphabet size */

/* Do not modify */
#define VALID (1)						/* Signals for valid and invalid */
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

/* TODO: Allow dynamic setting via cmd line */

/*************************** global variables ******************************/

/* TODO: BUG: If not stopping at max string, only after queue is empty,
it will rarely segfault. Unknown cause. TODO bugfix this */

/* Alphabet used -> alpha_size subset of alpha[] */
int alpha_size = STANDARD_ALPHABET_SIZE;
char alpha[] = {'0','1','2','3','4','5', '6', '7', '8', '9'};
char last;	/* Last char to process in alphabet. Set in main */

/* Theoretical longest word per alphabet size. l(n) = n(l(n-1) +2) */
int max_word_size[] = {0, 2, 8, 30, 128, 650, 3912, 27398, 219200, 1972818};
/* Properties of max word found */
int max_length = 0;
char *max_word;

int my_id;		/* Process ID; set in main */
/* V2 */
int *BOT;	/* Box_occurance_table. Bit table of if box has occured */
long *suff_table;	/* Suffix table for strings */
int suffix_depth;
/*********************** function prototypes *******************************/

void process (char *word, int depth);		/* TODO: convert to inline*/
int box_ok(char *word, int len); 			/* TODO: refactor and use below */
int box_valid(char *word, int len);			
int deep_check(char *string, int len);		
void manual_check(char *string, int len);
void process_1(char *word, int depth);

/**************************** functions ***********************************/

/* ROOT_PROCESS handles the queue, while all other processes do the
 * calculations. Communication is done via tags to denote message types.
 * Communcation cycle as follows:
 * ROOT does one iteration and adds items to queue, then sends one job to each
 * workers. The workers process the job, and once reach specified depth, sends a
 * signal back to the ROOT to recieve the item to be queued, followed by the
 * string itself. Once complete, the worker sends the maximum size it found back
 * to ROOT. if this is longer than previous maximum found by all processes, ROOT
 * requests the string and stores it. It then sends the next workload to the
 * process, and the cycle repeats until a) the queue is empty, or b) a string of
 * theoretical maximum length is found. A stop signal is then sent to all
 * workers (at appropriate times), and once all workers have stopped, ROOT
 * reports the maximum and stops itself.
 *
 * @args none
*/
int main(int argc, char *argv[])
{
	/* Setup */
	int ierr = 0, nbr_procs = 0, rec_count, send_count, an_id, i = 0,
	queue_size = 0, desired_size, stop = 0;
	MPI_Status status;
	char *w; 
	/* Keeps track of how many processes have been successfully stopped */
	int active_processes, waiting_count = 0;
	int *waiting_for_work;	/* For when queue is empty and continuing */

	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &nbr_procs);

	if (argc == 2) {
		alpha_size = atoi(argv[1]);
		if (alpha_size < 2 || alpha_size > 9) {
			if (my_id == ROOT_PROCESS) {	/* prevents excess messages */
				printf("Invalid/unknown alphabet size. Using N=3\n");
			}
			alpha_size = STANDARD_ALPHABET_SIZE;
		}
		else {
			if (STATUS) {
				if (my_id == ROOT_PROCESS) printf("Alphabet size: %d\n", alpha_size);
			}
		}
	} else {
		if (STATUS) {
			if (my_id == ROOT_PROCESS) {
				printf("Using default alphabet size: %d\n", STANDARD_ALPHABET_SIZE);
			}
		}
	}

	if (nbr_procs <= 1) {
		if (my_id == ROOT_PROCESS) {
			printf("Usage: \"mpirun -np X %s ALPHA_SIZE\" where X is number processess > 1\n",
				argv[0]);
		}
		goto end;
	}
	/*TODO: Convert desired_size and malloc to handle alpha size from cmd */
	desired_size = max_word_size[alpha_size];
	suffix_depth = alpha_size + 1;
	last = '/' + alpha_size;
	max_word = (char*) malloc(sizeof(char) * (desired_size + 1));  /* +1 for '\0' */
	max_word[0] = '\0';

	waiting_for_work = malloc(sizeof(int) * nbr_procs);
	for (i = 0; i < nbr_procs; i++) {
		waiting_for_work[i] = 0;
	}

	active_processes = nbr_procs - 1; /* Not including root */

	w = (char*) malloc(sizeof(char) * (max_word_size[alpha_size] + 1));

	/* Root/Master */
	if (my_id == ROOT_PROCESS) {
		/* Initial run to create branches for workers */
		/* TODO: Create initial string from 1 of each alphabet character */
		w[0] = '0';
		w[1] = '\0';
		process(w, INITIAL_SEARCH_DEPTH);

		MPI_Barrier(MPI_COMM_WORLD);

		/* Send initial instructions to each process */
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
			if (queue_size <= 0 || max_length >= desired_size) {
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
					if (PRINT_MAX) {
						printf("New max: %d %s\n", max_length, max_word);
					}
				}	
				/* Send stop command if needed... */
				if (stop) {
					/* Keep track of number of processes still running.
					 * If all stopped, then ROOT_PROCESS can also end */
					 MPI_Send(&send_count, 1, MPI_INT, an_id, STOP_TAG,
					 	MPI_COMM_WORLD);
					active_processes--;
					if (active_processes == 0) {
						break;
					}
				} /* ...or send the next work available in queue */
				else if (!dequeue(w)) {
					/* Queue temporarily empty;
					 wait until next queue request to send data to process
					 TODO: Verify that this works correctly. Hopefully never
					 needs to run */
					 if (STATUS) printf("stalled process %d due to empty queue\n",
					 	an_id);
					 waiting_for_work[an_id] = 1;
					 waiting_count++;
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
				w[rec_count] = '\0';
				if (DEBUG) {
					printf("(n)About to enqueue %s of size %d\n", w, rec_count);
				}
				enqueue(w);
			}
			queue_size = get_queue_size();
			/* Send to stalled processes due to empty queue. Ideally never
			   needs to run */
			if (waiting_count > 0 && queue_size >= waiting_count) {
				if (STATUS) {
					printf("Recover from a stall due to empty queue\n");
				}
				for (i = 0; i < nbr_procs; i++) {
					if (waiting_for_work[i] == 1) {
						dequeue(w);
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
				waiting_count = 0;
				queue_size = get_queue_size();
			}
		} 
		/* qs == 0 => explored all || ml == ds => found a longest pattern */
		/* Send stop signal to all processes */

		/* print status */
		printf("\n*********** FINISHED *********\n");
		printf("max == %d\n", max_length);
		printf("max word == %s\n", max_word);
		printf("queue size == %d\n", get_queue_size());
		printf("max valid? %d\n", deep_check(max_word, max_length));
		printf("******************************\n\n");
	}
	
	/* Worker/Slave */
	else {
		MPI_Barrier(MPI_COMM_WORLD);

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
				if (DEBUG) {
					printf("*:rec_count == %d\n", rec_count);
				}
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
	free(waiting_for_work);
	if (STATUS) printf("***** PROCESS %d STOPPED ******\n", my_id);
	return EXIT_SUCCESS;
}
/*  Jaco's process code. TODO: Refactor 
 *	Adds characters to pattern until certain depth is reached, checking that
 *	each new addition is box-valid. Once certain depth is reached, the remaining
 *	possible branches are added to the queue for future iterations to handle.
 *
 * @param word	the string/pattern so far
 * @param depth	how far to travel along the branches before queueing
*/
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
	/* TODO: remove this once queueing issue is gone */
	/* if (deep_check(word,i) == INVALID) return; */
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
			if (DEBUG) printf("(%d)string too deep %s\n", my_id, word);
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
		if (word[i] > last) {
			if (DEBUG) {
				printf("done with this branch\n");
			}
			word[i] = '\0'; /* Temp */
			i--;
			continue;
		}
		valid = box_valid(word,i);
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
	if (DEBUG) {
		printf("Done process\n");
	}
}

/* Jaco's Box_ok function */
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

/* checks if the last character added to a pattern is box-valid (i.e. the box
 * ending in the len'th character is never repeated.
 *
 * @param pattern	pattern/string to be checked
 * @param len		length of pattern
 * @return		VALID or INVALID respectively
 */

/* TODO: Refactor */
int box_valid(char* pattern, int len) {
	char c, c_end;
	int count = 1, b1_end = len, b1_start = b1_end, b2_start, b2_end, i = 0;
	int lb1, lb2;
	/* Step back until find char at pattern[len-1]
	 * Remember this box. Find 3rd repeat, check box. find 4th repeat. Check 3-4
	 * box and compare to 1-2 box. Repeat this for all possible boxes 
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

/*	Iterates over an entire pattern to confirm if it box-valid
 *	
 *	@param str	string to be checked
 *	@param len	length of string
 *	@return VALID or INVALID respectively
 */
int deep_check(char* str, int len) {
	char *c = str;
	char *w = malloc(sizeof(char) * (len+1));
	int i, valid = 1;
	for (i = 0; i < len; i++) {
		w[i] = c[i];
		w[i+1] = '\0';
		valid = box_valid(w, i);
		if (valid == 0){
			free(w);
			return INVALID;
		}
	}
	free(w);
	return VALID;
}

/*	Iterates over the entire pattern, adding characters iterably and states
 *	whether the pattern is box-valid at each character added.
 *	(Like deepcheck, except reports via stdout at each iteration)
 *
 *	@param str 	string to be checked
 *	@param len	length of string
 */
void manual_check(char *str, int len) {
	char *c = str;
	char *w = malloc(sizeof(char) * (len + 1));
	int i, valid = 1;
	printf("local test\n");
	for (i = 0; i < len; i++) {
		w[i] = c[i];
		w[i+1] = '\0';
		valid = box_ok(w,i);
		printf("i:%d\tValid:%d\n", i, valid);
	}
	free(w);
}

/* Refactoring in progress. Do not use */
void process_1(char *word, int depth) {
	int k = strlen(word), i = k, valid = 1, diff = alpha[0] + alpha_size - 1;
	char nul = alpha[0] - 1;

	word[i] = nul;
	while (i >= k) {
		if (i == k + depth) {
			word[i] = '\0';
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
			i--;
			continue;
		}
		word[i]++;
			/* outside bounds*/
		if (word[i] >  diff) {
			i--;
			continue;
		}
		if (box_ok(word, i)) {
			i++;
			if (i > max_length) {
				max_length = i;
				word[i] = '\0';
				strcpy(max_word, word);
			}
			word[i] = nul;
		}
		/* Insert */
		/* Check if valid */
		/* if valid, increment i */
			/* if i = k+depth, queue and i-- */
		/* increment w[i] */
		/* if w[i] > '0' +  alpha_size - 1*/
			/* i-- */
	}
}


