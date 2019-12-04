/******************************************************************************
 *	Finds the longest string that has no repeated boxes. Utilises MPI to run
 *	multiple workers at once. Relies on a queue-based system to handle jobs
 *	between workers.
 *
 *  original process and box_ok code provided by Jaco Geldenhuys
 * 
 *	compile:	mpicc -o ../bin/search.o queue.c search.c -lm
 *	run:		mpirun -np X ./search.o	N
 		where X is number processes,
	 	N alphabet size between 2 and 9 (default if not specified)

 * WARNING: Since one process is used for controlling the rest, do not specify
 * more than the available number of cores. Else the program can lock (admin
 * process essentially waits for itself) 
 *	@author C. R. Zeeman (caleb.zeeman@gmail.com)
 *	@version 1.4
 *	@date 2019-12-04
 *****************************************************************************/
/* Includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <math.h>
#include <limits.h>
#include "queue.h"
#include "search.h"

/************************** Definitions and macros **************************/
/* Configurations */
#define BUGFIX (0)
#define DEBUG (0)						/* Enables debug output */
#define DEBUG_V2 (0)					/* ^ for bitwise method */
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

/*************************** global variables ******************************/

/* TODO: BUG: If not stopping at max string, only after queue is empty, and 
processes = 4, it will rarely segfault. Unknown cause. TODO bugfix this
Error in enqueue (failure at address 0x8). Hopefully will be fixed when
queue is rewritten to new method */

/* Alphabet used -> alpha_size subset of alpha[] */
int alpha_size = STANDARD_ALPHABET_SIZE;
char alpha[] = {'0','1','2','3','4','5', '6', '7', '8', '9'};
char last;	/* Last char to process in alphabet. Set in main */

/* Theoretical longest word per alphabet size. l(n) = n(l(n-1) +2) */
int max_word_size[] = {0, 2, 8, 30, 128, 650, 3912, 27398, 219200, 1972818};
/* Properties of max word found */
int max_length = 0;
char *max_word;
/* Number of 'characters' that can fit in a long long */
int chars_per_long;

int my_id;		/* Process ID; set in main */

/* V2 */
unsigned long *ST;	/* Suffix_table for current. Table 1 */
int suffix_depth;
int *LOT;	/* Last_occurance_table; where alpha char last occured. Table 2 */
void *BOT; /* Box_occurance_table. Bit table of if box has occured. Table 3 */
unsigned long long BOT_size;	/* Number of bits in BOT */

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
	int ierr = 0, nbr_procs = 0, rec_count, send_count, an_id, i = 0, j = 0,
	queue_size = 0, desired_size, stop = 0;
	long partial_sum = 0;
	MPI_Status status;
	char *w; 
	/* Keeps track of how many processes have been successfully stopped */
	int active_processes, waiting_count = 0;
	int *waiting_for_work;	/* For when queue is empty and continuing */

	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &nbr_procs);
	/* Command line arguement details */
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
	/* constants and variable setup */
	desired_size = max_word_size[alpha_size];
	last = '/' + alpha_size;
	max_word = (char*) malloc(sizeof(char) * (desired_size + 1));  /* +1 for '\0' */
	max_word[0] = '\0';
	chars_per_long = sizeof(long long int) * 8 / (log(alpha_size) / log(2));
	waiting_for_work = malloc(sizeof(int) * nbr_procs);
	for (i = 0; i < nbr_procs; i++) {
		waiting_for_work[i] = 0;
	}

	active_processes = nbr_procs - 1; /* Not including root */

	w = (char*) malloc(sizeof(char) * (max_word_size[alpha_size] + 1));

	/* V2 setup */
	ST = malloc(sizeof(long) * (alpha_size+1)); /* memset in process */
	LOT = malloc(sizeof(int) * alpha_size);

	suffix_depth = alpha_size + 1;
	/*bit size of BOT. e.g. n=4: 4^1 + 4^2 + ... + 4^5 */
	for (i = 1; i <= suffix_depth; i++) {
		partial_sum = 1;
		for (j = 1; j <= i; j++) {
			partial_sum *= alpha_size;
		}
		BOT_size += partial_sum;
	}
	if (DEBUG_V2 && my_id == ROOT_PROCESS) printf("BOT_size: %lld\n", BOT_size);
	/*Padding to next byte */
	BOT_size += BOT_size % 8;
	BOT = malloc(BOT_size);

	/* Root/Master */
	if (my_id == ROOT_PROCESS) {
		/* Initial run to create branches for workers */
		/* Initial start is 1 of each alphabet character */
		for (i = 0; i < alpha_size; i++) {
			w[i] = alpha[i];
		}
		w[i] = '\0';

		if (STATUS) printf("char per long: %i\n", chars_per_long);
		
		process(w, INITIAL_SEARCH_DEPTH);

		MPI_Barrier(MPI_COMM_WORLD);

		/* Send initial instructions to each process */
		for (an_id = 1; an_id < nbr_procs; an_id++) {
			if (!dequeue_s(w)) {
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
			if (queue_size <= 0/*|| max_length >= desired_size*/) {
				stop = 1;
			}
			/* Recieve a reply */
			ierr = MPI_Recv(&rec_count, 1, MPI_INT, MPI_ANY_SOURCE,
				MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if (DEBUG) {
				printf("Main process recieved a command: %d\n", status.MPI_TAG);
			}
			if (BUGFIX && stop > 0) printf("(%d)stop == 1\n", my_id);
			an_id = status.MPI_SOURCE;

			/* Recieving max found from process */
			if (status.MPI_TAG == SENDING_MAX_TAG) {
				if (stop && BUGFIX) printf("<>");
				if (rec_count > max_length) {
					max_length = rec_count;
					if (DEBUG) printf("(n)New max found, requesting string\n");
					/* Request longest string*/
					MPI_Send(&rec_count, 1, MPI_INT, an_id, SEND_MAX_TAG,
						MPI_COMM_WORLD);
					MPI_Recv(max_word, max_length, MPI_CHAR, an_id,
						MAX_STRING_TAG, MPI_COMM_WORLD, &status);
					max_word[rec_count] = '\0';
					if (PRINT_MAX) {
						printf("New max: %d %s\n", max_length, max_word);
					}
				}	
				/* Send stop command if needed... */
				if (stop) {
					if (BUGFIX){
						printf("Sending stop to process %d. remaining: %d\n", an_id, active_processes-1);
					}
					/* Keep track of number of processes still running.
					 * If all stopped, then ROOT_PROCESS can also end */
					 MPI_Send(&send_count, 1, MPI_INT, an_id, STOP_TAG,
					 	MPI_COMM_WORLD);
					active_processes--;
					if (active_processes == 0) {
						break;
					}
				} /* ...or send the next work available in queue */
				else if (!dequeue_s(w)) {
					/* Queue temporarily empty;
					 wait until next queue request to send data to process
					 TODO: Verify that this works correctly. Hopefully never
					 needs to run */
					if (STATUS || BUGFIX) printf("stalled process %d due to empty queue\n",
						an_id);
					assert(an_id < nbr_procs);
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
				if (stop) { /* New work  TODO rework */
					if (BUGFIX) printf("queueing while stop > 0 by %d\n", an_id);
					stop = 0;
				}
			/* TODO: Move to new enqueue code 
				enqueue_d(an_id, rec_count); */
				assert(rec_count <= max_word_size[alpha_size]);
				MPI_Recv(w, rec_count, MPI_CHAR, an_id, QUEUE_TAG,
					MPI_COMM_WORLD, &status);
				w[rec_count] = '\0';
				enqueue_s(w);
			}
			queue_size = get_queue_size();
			/* Send to stalled processes due to empty queue. Ideally never
			   needs to run */
			if (waiting_count > 0 && queue_size >= waiting_count) {
				if (STATUS || BUGFIX) {
					printf("Recover from a stall due to empty queue\n");
				}
				for (i = 0; i < nbr_procs; i++) {
					if (waiting_for_work[i] == 1) {
						dequeue_s(w);
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
				ierr = MPI_Recv(w, rec_count, MPI_CHAR, ROOT_PROCESS,
					WORK_TAG, MPI_COMM_WORLD, &status);
				assert(rec_count <= max_word_size[alpha_size]);
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
		if (BUGFIX) printf("process %d recieved stop command\n", my_id);
	}
	
	end: ierr = MPI_Finalize();

	free(waiting_for_work);
	free(ST);
	free(LOT);
	free(BOT);
	if (STATUS) printf("***** PROCESS %d STOPPED ******\n", my_id);
	return EXIT_SUCCESS;
}
/*  Jaco's process code merged with V2 additions
 *	Adds characters to pattern until certain depth is reached, checking that
 *	each new addition is box-valid. Once certain depth is reached, the remaining
 *	possible branches are added to the queue for future iterations to handle.
 *
 * @param word	the string/pattern so far
 * @param depth	how far to travel along the branches before queueing
*/

void process(char *word, int depth) {
	int k = strlen(word);
	int i = 0, j = 0, letter = 0, lc, c,d, counter;
	int valid = 1, next, o_valid, have_gone_back = 0; /*HGB: see BOT reverse */
	char *temp_ptr;
	unsigned long *new_ST;
	unsigned long long exponent = 0, exponent_2, exp_temp, box;
	unsigned long long word_v = 0, old_word_v = 0, temp = 0;
	new_ST = malloc(sizeof(long) * (alpha_size+1));
	/* Reset tables */
	memset(ST, 0, sizeof(long) * (alpha_size+1));
	memset(LOT, 0, alpha_size * sizeof(int));
	memset(BOT, 0, BOT_size/8);

	if (DEBUG) {
		printf("process(w==\"%s\")\n", word);
		printf("k==%d\n", k);
		printf("i==%d\n", i);
		printf("BOT:\n");
		temp_ptr = (char*) BOT;
		for (i = 0; i < BOT_size/8; i++) {
			printf("[%d]: %d\n", i, *(temp_ptr+i));
		}
	}
	/* Temp: Convert string to long for V2. Queue will do this later
	   stored in base (alpha_size) */
	   /* TODO: Figure out how to handle long string (multiple longs) */
	exponent = 1;
	for (j = k-1; j >= 0; j--) {	/* exponent = alpha_size ^ (k-1-j) */
		word_v += (word[j] - alpha[0]) * exponent;
		exponent *= alpha_size;
	}
	exponent /= alpha_size;
	if (DEBUG_V2) {
		printf("word = %s\t", word);
		printf("word_v = %lld\texp = %lld\n", word_v, exponent);
	}

	/* Reconstruct ST, LOT, BOT */
	old_word_v = word_v;
	exp_temp = exponent;
	for (i = 0; i < k; i++) {
		/* Obtain letter from 'encoded' number */
		for (j = 0; j < alpha_size; j++) {
			if (DEBUG_V2) printf("LOT[%d]: %d\n", j, LOT[j]);
		}
		if (exponent == 0) {
			letter = word_v;
		}
		else {
			letter = (word_v / exponent);
		}

		/* update BOT */
		assert(letter < alpha_size);
		lc = LOT[letter];
		if (DEBUG_V2) printf("i: %d, letter: %d; lc = %d\n", i, letter, lc);
		if (lc != 0) {	/* Box occured */
		/* Only adding most closest box as only one that will
		   appear when using suffix table */
			if (i - lc > alpha_size) {	/* +1 both sides cancel out */
			   continue;
			}
			if (DEBUG_V2) printf("exp_temp: %lld\n", exp_temp);
			exponent_2 = exp_temp;
			temp = old_word_v;
			if (DEBUG_V2) printf("lc: %d\n", lc);
			for (j = 1; j < lc; j++) {	/* remove left of box */
			   temp %= exponent_2;
			   exponent_2 /= alpha_size;
			}
			if (DEBUG_V2) printf("exponent: %lld\texponent_2: %lld\n", exponent, exponent_2);
			if (DEBUG_V2) printf("temp: %lld\tword_v: %lld\tword_v mod exp: %lld\n",
			temp, word_v, word_v % exponent);
			//box = (temp - (word_v % exponent)) / exponent;
			box = temp / exponent;	/* remove right of box */

			if (DEBUG_V2) printf("box found. Value: %lld\n", box);
			/* flip bit */
			if (box > BOT_size && STATUS) {
				/* Impossible situation; failsafe for testing. TODO: remove */
				printf("WARNING: BOX > BOT_SIZE. %lld vs %lld\n", box, BOT_size);
				printf("word: %s\n", word);
				printf("i == %d, lc == %d\n", i, lc);
				printf("temp: %lld\tword_v: %lld\tword_v mod exp: %lld\n",
					temp, word_v, word_v % exponent);
			}
			else {
				/* Move (box/8) bytes and change (box%8) bit */
				BOT += box / 8;
				temp_ptr = (char*) BOT;
				exponent_2 = 128;
				for (j = 0; j < box % 8; j++) {
					exponent_2 /= 2;
				}
				*temp_ptr ^= exponent_2;
				if (DEBUG_V2) printf("Update BOT entry %lld\n", box);
		   		BOT -= box / 8;
		   }
		}

		/* update LOT */
		assert(letter <alpha_size);
		LOT[letter] = i+1;

		/* Setup for next iteration */
		word_v %= exponent;
		exponent /= alpha_size;
	}
	word_v = old_word_v;

	/* Update suffix table (ST) */
	exponent = 1;
	for (i = 0; i < alpha_size; i++) {	/* exp ^ alpha_size*/
		exponent *= alpha_size;
	}
	if (DEBUG_V2) printf("exp: %lld\tword_v: %lld\n", exponent, word_v);
	
	for (i = alpha_size; i>=0; i--) {
		ST[i] = word_v;
		if (DEBUG_V2) printf("added suffix %lld at index %d\n", word_v, i);
		word_v %= exponent;
		exponent /= alpha_size;
	}

	/* TODO: remove this once queueing issue is gone */
	if (deep_check(word,i) == INVALID) return;

	/* Iterations */
	i = k;
	counter = i;
	word_v = old_word_v;
	word_v *= alpha_size;
	//word_v--;
	if (DEBUG_V2) printf("word_v == %lld\n", word_v);
	word[i] = '/';
	while (i >= k) {
		/*delme*/ 
		if (DEBUG_V2) printf("<%d>Bot now:\n", i);
		temp_ptr = BOT;
		for (c = 0; c < BOT_size/8; c++) {
			if (DEBUG_V2) printf("[%d]: [%d]\n", c, *(temp_ptr+c));
		}
		if (DEBUG) {
			printf("exploring position i==%d w==", i);
			for (int ii = 0; ii <= i; ii++) {
				printf("%c", word[ii]);
			} 
			printf("\n");
			printf("word_v at start: %lld\n", word_v);
		}
		if (i == k + depth) {
			if (DEBUG) {
				printf("too deep\n");
			}
			word[i] = '\0';
			if (DEBUG) {
				printf("(%d)string too deep %s\n", my_id, word);
			}
			if (valid) {	/* Failsafe check */
				if (my_id != ROOT_PROCESS) {
					if (DEBUG) {
						printf("(%d) asking to enqueue %s\n", my_id, word);
					}
					/* TODO: Move to new enqueue_d function*/
					MPI_Send(&i, 1, MPI_INT, ROOT_PROCESS, EXPECT_QUEUE_TAG,
						MPI_COMM_WORLD);
					MPI_Send(word, i, MPI_CHAR, ROOT_PROCESS, QUEUE_TAG,
						MPI_COMM_WORLD);
				} else {
					if (DEBUG) printf("(0)about to enqueue %s\n", word);
					enqueue_s(word);
				}
			}
			i--;
			word_v /= alpha_size;
			
			if (DEBUG_V2) printf("word_v / alpha_size. now %lld\n", word_v);
			/* Update LOT */
			j = word_v % alpha_size;
			c = i;
			temp = word_v / alpha_size;
			assert(j < alpha_size);
			LOT[j] = 0;
			if (DEBUG_V2) printf("Updating LOT\n");
			while (c >= 0) {
				d = temp % alpha_size;
				if (d == j) {
					LOT[j] = c;
					break;
				}
				temp /= alpha_size;
				c--;
			}
			continue;
		}
		assert(i <= max_word_size[alpha_size]);
		word[i]++;

		if (word[i] > last) {
			//if (word_v % alpha_size == alpha_size-1)
			if (DEBUG) {
				printf("done with this branch\n");
			}
			word[i] = '\0'; /* Safety */
			i--;
			/* Update (revert) BOT entry, if any */
			/*	Caution: If 'undoing' last character added, need to
			 first check if it is valid (else could remove an older box).
			 If not the last character to be added, you know it's valid
			 (but can't use o_valid), so need to do a different check.
			 Hence 'have_gone_back' variable */
			if (o_valid || have_gone_back) {
				exp_temp = alpha_size;
				letter = word_v % alpha_size;
				if (DEBUG_V2) printf("word_v: %lld\n", word_v);
				if (DEBUG_V2) printf("BOT_revert: letter: %d\n", letter);
				for (j = 1; j <= alpha_size; j++) { /* find same char */
					if (DEBUG_V2) printf("BOT_revert: searching: %lld\n", ((word_v % (exp_temp*alpha_size)) /exp_temp));
					if (((word_v % (exp_temp*alpha_size)) /exp_temp) == letter) {
						box = word_v % (exp_temp*alpha_size); /* box value */
						if (DEBUG_V2) printf("BOT_revert: box %lld at j: %d\n", box, j);
						BOT += box / 8;
						temp_ptr = (char*) BOT;
						exp_temp = 128;
						for (j = 0; j < box % 8; j++) {
							exp_temp /= 2;
						}
						*temp_ptr ^= exp_temp;
						BOT -= box / 8;
						break;
					}
					exp_temp *= alpha_size;
				}
			}
			have_gone_back = 1;
			/* Update (revert to old) suffix table */
			if (DEBUG_V2) printf("undo: new_st:\n");
			for (j = 0; j < alpha_size; j++) {
				if (DEBUG_V2) printf("ST[%d] = %ld\n",j+1, ST[j+1]);
				new_ST[j] = ST[j+1] / alpha_size;
				if (DEBUG_V2) printf("new_st: [%d] - %ld\n", j, new_ST[j]);
			}
			exp_temp = 1;
			for (j = 0; j <= alpha_size; j++) {
				exp_temp *= alpha_size;
			}
			word_v /= alpha_size;
			if (DEBUG_V2) printf("word_v / alpha_size. now %lld\n", word_v);
			new_ST[alpha_size] = word_v % exp_temp;
			if (DEBUG_V2) {
				printf("word_v == %lld\n", word_v);
				printf("new_st[%d] == %ld\n", alpha_size,new_ST[alpha_size]);
				printf("\n");
			}
			for (j = 0; j <= alpha_size; j++){
				ST[j] = new_ST[j];
			}
			/* Update LOT */
			j = word_v % alpha_size;
			c = i;
			temp = word_v / alpha_size;
			assert(j < alpha_size);
			if (DEBUG_V2) printf("updating LOT\n");
			LOT[j] = 0;
			while (c >= 0) {
				d = temp % alpha_size;
				if (d == j) {
					LOT[j] = c;
					break;
				}
				temp /= alpha_size;
				c--;
			}
			continue;
		}
		if (DEBUG_V2) {
			printf("LOT:\n");
			for (j = 0; j < alpha_size; j++) {
				printf("[%d]: %d\n",j, LOT[j]);
			}
		}
		// else word_v++
		if (DEBUG_V2) printf("word[i] == %c\n", word[i]);
		if (word[i] != alpha[0]) {
			/* Update (revert) BOT entry, if any */
			if (o_valid || have_gone_back) {
				exp_temp = alpha_size;
				letter = (word_v % alpha_size);
				if (DEBUG_V2) printf("word_v: %lld\n", word_v);
				if (DEBUG_V2) printf("BOT_revert: letter: %d\n", letter);
				for (j = 1; j <= alpha_size; j++) { /* find same char */
					if (DEBUG_V2) printf("BOT_revert: searching: %lld\n", ((word_v % (exp_temp*alpha_size)) /exp_temp));
					if (((word_v % (exp_temp*alpha_size)) /exp_temp) == letter) {
						box = word_v % (exp_temp*alpha_size); /* box value */
						if (DEBUG_V2) printf("BOT_revert: box %lld at j: %d\n", box, j);
						BOT += box / 8;
						temp_ptr = (char*) BOT;
						exp_temp = 128;
						for (j = 0; j < box % 8; j++) {
							exp_temp /= 2;
						}
						*temp_ptr ^= exp_temp;
						BOT -= box / 8;
						break;
					}
					exp_temp *= alpha_size;
				}
			}
			word_v++;
			have_gone_back = 0;
			if (DEBUG_V2) printf("word_v++. now %lld\n", word_v);
		}
		/* clear new_table and check if valid. If valid, make new table
			the main table */
		/* Create new suffix table*/
		memset(new_ST, 0, sizeof(long) * alpha_size+1);
		next = word[i] - alpha[0];	/* TODO: Change to new method */

		if (DEBUG_V2) printf("added: %c\n", word[i]);
		new_ST[0] = next;
		if (DEBUG_V2) {
			printf("i = %d\n", i);
			printf("Addition: new_st:\nnew_st: [%d] - %ld\n", 0, new_ST[0]);
		}
		/* Seperate logic for if adding a new character to string, or simply 
			incrementing an existing one (no length increase) */
		if (next != 0) {	/* Increment alphabet char only */
			for (j = 1; j <= alpha_size; j++) {
				new_ST[j] = ST[j] + 1;
				if (DEBUG_V2) printf("new_st: [%d] - %ld\t st: [%d] - %ld\n",
					j, new_ST[j], j, ST[j]);			
			}
		}
		else {
			for (j = 1; j <= alpha_size; j++) {
				new_ST[j] = ST[j-1] * alpha_size + next;
				if (DEBUG_V2) printf("new_st: [%d] - %ld\n", j, new_ST[j]);
				counter++;
			}
		}
		if (DEBUG_V2) printf("\n");
		for (j = 0; j <= alpha_size; j++) {
			ST[j] = new_ST[j];
		}
		/* ---- V2: check if valid ---- */
		/* find box */
		/* TODO: Replace i with len everywhere, keep track of length
		using it. Speeds up other functions too */
		box = -1;
		temp_ptr = BOT;
		if (DEBUG_V2) {
			printf("BOT:\n");
			for (j = 0; j < BOT_size/8; j++) {
				printf("[%d]: %d\n", j, *(temp_ptr+j));
			}
			printf("ST:\n");
			for (j = 0; j <= alpha_size; j++) {
				printf("[%d]: %ld\n", j, ST[j]);
			}
			printf("LOT:\n");
			for (j = 0; j < alpha_size; j++) {
				printf("[%d]: %d\n",j, LOT[j]);
			}
		}
		assert(next < alpha_size);
		for (j = i; j > i-alpha_size; j--) {
			if (LOT[next] == j) {	/* box occured */
				box = ST[i+1-j]; 
				if (DEBUG_V2) printf("box occured at letter %d, j=%d\n", next, j);
				goto e_loop;
			}
		}
		e_loop:
		if (box == -1) { /* No box */
			o_valid = 1;
		} else {
			/* check if box ticked & update BOT if valid */
			/* Move (box/8) bytes and change (box%8) bit */
			if (DEBUG_V2) printf("box value: %lld\n", box);
			BOT += box / 8;
			exp_temp = 128;
			for (j = 0; j < box % 8; j++) {
				exp_temp /= 2;
			}
			temp_ptr = (char*) BOT;
			c = (*temp_ptr) & (exp_temp);
			if (DEBUG_V2) printf(">c == %d\n", c);
			if (c == 0) {
				if (DEBUG_V2) printf(">updating box\n");
				o_valid = 1;
				*temp_ptr ^= exp_temp;
			} else {
				o_valid = 0;
			}
			if (DEBUG_V2) {
				printf("Bot now:\n");
				temp_ptr = BOT;
				for (c = 0; c < BOT_size/8; c++) {
					printf("[%d]: [%d]\n", c, *(temp_ptr+c));
				}
				BOT -= box/8;
			}
		}
		if (o_valid == 1) {
			have_gone_back = 0;
			/* Update LOT */
			if (DEBUG_V2) printf("LOT[%d]: %d -> %d\n", next, LOT[next], i+1);
			LOT[next] = i+1;
		}
		/*	------------------	*/
		valid = box_valid(word,i);
		if (DEBUG_V2) {
			printf("word_v == %lld\n", word_v);
			printf("o_valid == %d, valid == %d, 0_valid == valid: %d\n",
				o_valid, valid, o_valid==valid);
		}
		if (STATUS && o_valid != valid) {
			printf("Warning: o_valid == %d, valid == %d, 0_valid == valid: %d\n",
				o_valid, valid, o_valid==valid);
			return;
		}
		if (valid) {
			if (DEBUG) {
				printf("no repeat box\n");
			}
			for (j = 0; j <= alpha_size; j++){
				ST[j] = new_ST[j];
			}
			if (DEBUG_V2) printf("swapped: ST[] = {%ld, %ld, %ld}\n", ST[0], ST[1], ST[2]);

			i++;
			word_v *= alpha_size;
			if (DEBUG_V2) printf("word_v * alpha_size. now %lld\n", word_v);
			assert(i <= max_word_size[alpha_size]);
			if (i > max_length) {
				max_length = i;
				word[i] = '\0';
				strcpy(max_word, word);
			}
			/* Setup for next run */
			word[i] = '/';
		} else if (DEBUG && 0) {
			printf("repeated box\n");
		}
	}
	if (DEBUG) {
		printf("Done process\n");
	}
	free(new_ST);
}
/* TODO: Create seperate functions for handling the long long[]
	Convert to inline later */

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
	assert(b1_start <= max_word_size[alpha_size]);
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

/* TODO: Once O(1) search complete, convert deep_check to use the 3 tables
   to check for repeating boxes */
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
		valid = box_valid(w,i);
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
				if (DEBUG) {
					printf("(1)about to enqueue %s\n", word);
				}
				enqueue_s(word);
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
		if (box_valid(word, i)) {
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


