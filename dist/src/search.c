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
 * 
 * If attempting to use debugging flags, strongly recommend
 * piping output to file.
 * 
 *	@author C. R. Zeeman (caleb.zeeman@gmail.com)
 *	@version 1.9
 *	@date 2019-12-10
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
#define BUGFIX (0)			/* DELME */
#define STRING_TEST (0) 	/* DELME */
#define DEEP_TEST (0) 		/* DELME */
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

/* Macros to avoid function calls. Slight speed boost */
/* Revert BOT entry, if any */
/* TODO: Fix bot revert to handle (char_len+1). i=33 doesn't seem to be reset */
/* D keeps track of when to switch to previous long's values. */
#define REVERT_BOT_ENTRY(d, j, word_v, exp_temp, exponent_2, letter, temp_ptr,temp, box, word_arr) do { \
	unsigned long long preamble = 0; \
	exp_temp = alpha_size;	\
	d = count-1; \
	temp = word_v; \
	letter = (word_v % alpha_size);	\
	if (d % chars_per_long == 0 && count % chars_per_long != 0) box = letter; \
	else box = 0; \
	exponent_2 = 1;	/*TODO TEST NEW METHOD */\
	if (DEBUG_V2) printf("word_v: %llu\n", word_v);	\
	if (DEBUG_V2) printf("BOT_revert: letter: %d\n", letter);	\
	for (j = 1; j <= alpha_size; j++) { /* find same char */	\
		if (DEBUG_V2) printf("d: %d, count %d\n", d, count); \
		/*if (d % (chars_per_long+1) == (chars_per_long) && count % (chars_per_long+1) != chars_per_long) { */\
		if (d % chars_per_long == 0 && count % chars_per_long != 0) { \
			if (DEBUG_V2) printf("swapped to prev. long\n"); \
			preamble = box; \
			exponent_2 = exp_temp; \
			exp_temp = 1; /* Alpha_size? */ \
			temp = word_arr[count / chars_per_long -1]; \
		} \
		if (DEBUG_V2) printf("BOT_revert: searching: %llu\n", ((temp % (exp_temp*alpha_size)) /exp_temp));	\
		box = temp % (exp_temp*alpha_size) * exponent_2 + preamble; /* box value TODO CHECK*/	\
		if (((temp % (exp_temp*alpha_size)) /exp_temp) == letter) {	\
			if (DEBUG_V2) printf("BOT_revert: box value: %llu at j: %d\n", box, j);	\
			if (DEBUG_V2) printf("exp_temp: %llu, temp: %llu, exponent_2: %llu, preamble: %llu\n", exp_temp, temp, exponent_2, preamble); \
			BOT += box / 8;	\
			temp_ptr = (char*) BOT;	\
			if (DEBUG_V2) printf("bot before: %u\n", *temp_ptr); \
			flip = 128;	\
			for (j = 0; j < box % 8; j++) {	\
				flip /= 2;	\
			}	\
			/* TEST */ \
			/*PROBLEM: TODO: Used to be OR (setting not revert) (fixed?) */ \
			flip = 255 - flip; \
			*temp_ptr &= flip;	\
			if (DEBUG_V2) printf("bot now: %u, flip: %u\n", *temp_ptr, flip); \
			BOT -= box / 8;	\
			break;	\
		}	\
		exp_temp *= alpha_size;	\
		d--;\
	}	\
} while (0)

#define CALC_SPILLOVER_BOX(lc, temp, exponent, exponent_2, exponent_3, word_arr,i, box, it) do { \
	if (lc != 0 && lc % (chars_per_long) == 0) lc = chars_per_long; \
	else lc %= (chars_per_long); \
	/* calculate box from end of long to current position */ \
	if (DEBUG_V2) printf("spillover box. it: %d\n", it); \
	box = temp / exponent; /* remove right of box */ \
	/* calculate box segment in previous long */ \
	temp = word_arr[it-1]; \
	exponent_3 = pow(alpha_size, i+1); \
	exponent_2 = pow(alpha_size, chars_per_long - lc + 1); \
	if (DEBUG_V2) printf("exponent: %llu\texponent_2: %llu\n", exponent, exponent_2); \
	if (DEBUG_V2) printf("temp: %llu\n", temp); \
	temp %= exponent_2; /* Remove left of box */\
	/* add together */ \
	if (DEBUG_V2) printf("box before add: %llu\n", box); \
	if (DEBUG_V2) printf("temp: %llu\n", temp); \
	if (DEBUG_V2) printf("exponent_3: %llu\n", exponent_3); \
	box += temp * exponent_3; /* TODO check if correct */ \
	if (DEBUG_V2) printf("new box: %llu\n", box); \
	assert(box <= BOT_size); \
	if (DEBUG_V2) printf("end of spillover method\n"); \
} while(0)
/* Two bugs: 1: After revert, if LOT was at i=i-1, doesn't track
2: [0] often gets set to 32 (first of last box) */

#define REVERT_LOT(j, c, temp, it, word_arr, d, label_name) do { \
	c += i / chars_per_long; /* TODO Check. Fix for over-mods */ \
	c %= (chars_per_long+1); /* Fix for fix */ \
	/* TODO Fixes above still needed? */ \
	if (DEBUG_V2) printf("updating (reverting) LOT\n"); \
	if (DEBUG_V2) printf("j: %d, c: %d, temp: %llu\n", j, c, temp); \
	LOT[j] = -1; /*TODO: Rework equation below */ \
	while (c > 0) { \
		d = temp % alpha_size; \
		if (DEBUG_V2) printf("c: %d d: %d, temp: %llu\n",c, d, temp); \
		if (d == j) { \
			if (DEBUG_V2) printf("c==d at val %d\n", d); \
			LOT[j] = c + ((int) (i / chars_per_long)) * chars_per_long; \
			if (i % chars_per_long == 0 && i != 0) { \
				/* TODO: FIx needed? */ \
				LOT[j] -= chars_per_long; \
			} \
			break; \
		} \
		temp /= alpha_size; \
		c--; \
	} \
	/* TODO: Update to search previous longs too */ \
	if (LOT[j] == -1) { \
		if (DEBUG_V2) printf("searching previous longs for %d. ii: %d\n",j, (count / chars_per_long) - 1); \
		for (ii = (count / chars_per_long) -1; ii >= 0; ii--) { \
			if (DEBUG_V2) printf("word_arr[%d] == %llu\n", ii, word_arr[ii]); \
			temp = word_arr[ii]; \
			c = chars_per_long; /* TODO: -1? */ \
			while (c > 0) { \
				d = temp % alpha_size; \
				if (DEBUG_V2) printf("c: %d d: %d, d==j: %d, temp: %llu\n",c, d, d==j,temp); \
				if (d == j) { \
					if (DEBUG_V2) printf("d==j\n"); \
					LOT[j] = c + (ii * chars_per_long); \
					goto label_name; \
				} \
				temp /= alpha_size; \
				c--; \
			} \
			if (DEBUG_V2) printf("out of loop\n"); \
		} \
	} \
} while (0)

/* TODO: BUG: arr[] had 0 entries when string > chars_per_long
	Error in saving. Likely in long_count and getting set to 0 */
/*************************** global variables ******************************/

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

void process (char *word, int depth, unsigned long long *word_arr, int len);		/* TODO: convert to inline*/
int box_ok(char *word, int len); 			/* TODO: refactor and use below */
int box_valid(char *word, int len);			
int deep_check(char *string, int len);		
int deep_check_l(unsigned long long *word_arr, int len);		
char *longs_to_string(unsigned long long *arr, int len);
unsigned long long *string_to_longs(char *string, int len);


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
/* TODO: N=1 breaks chars_per_long, and freezes program */
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
				printf("Invalid/unknown alphabet size. Using default.\n");
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
	/* TODO: If initial_depth^alpha_size < nbr_procs, initial_depth++ or smth*/
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
	if ( my_id == ROOT_PROCESS) printf("BOT_size: %llu\n", BOT_size);
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
		
		process(w, INITIAL_SEARCH_DEPTH, NULL, 0);

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
		/* DELME */
		if (STRING_TEST) {
			printf("Test: L->S\n");
			unsigned long long *arr = malloc(sizeof(unsigned long long) * 2);
			char * str;
			arr[0] = 1; str = longs_to_string(arr, 2);
			printf("'01' -> '%s'\n", str);
			free(str); arr[0] = 5; str = longs_to_string(arr, 3);
			printf("'101' -> '%s'\n", str);
			free(str); arr[0] = 3; str = longs_to_string(arr, 2);
			printf("'11' -> '%s'\n", str);
			free(str);
			arr[0] = 0x8000000000000000;
			arr[1] = 3;
			str = longs_to_string(arr, 66);
			printf("'100.....11' len=66 -> '%s'\n", str);
			printf("Test: S->L\n");
			unsigned long long *result;
			free(str); str = longs_to_string(arr, 66);
			result = string_to_longs(str, 66);
			free(str);
			printf("'01' -> %llu\n", (string_to_longs("01",2))[0]);
			printf("'11' -> %llu\n", (string_to_longs("11",2))[0]);
			printf("'101' -> %llu\n", (string_to_longs("101",3))[0]);
			printf("arr[0] == result[0]: %d\n", result[0] == arr[0]);
			printf("arr[1] == result[1]: %d\n", result[1] == arr[1]);
			free(result); free(arr);
		}
		if (DEEP_TEST) {
			unsigned long long *arr2 = malloc(sizeof(unsigned long long) * 2);
			arr2[0] = 1; arr2[1] = 0;
			printf("Deep_check test: '' -> valid?\n");
			arr2[0] = 1; printf("01 -> %d\n", deep_check_l(arr2,2));
			arr2[0] = 0; printf("00 -> %d\n", deep_check_l(arr2,2));
			arr2[0] = 2; printf("10 -> %d\n", deep_check_l(arr2,2));
			arr2[0] = 5; printf("101 -> %d\n", deep_check_l(arr2,3));
			arr2[0] = 7; printf("111 -> %d\n", deep_check_l(arr2,3));


		}
		/* <DELME> */
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
					if (STATUS) printf("stalled process %d due to empty queue\n",
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
				if (STATUS) {
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
				process(w, STANDARD_SEARCH_DEPTH, NULL, 0);		

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

/* Generally: i is index in array, j, k are actual lengths
i.e. '0':  i=0, j=k=1 */
void process(char *word, int depth, unsigned long long *word_arr, int len) {
	int k = strlen(word), size, count, long_count;
	int i = 0, j = 0, letter = 0, lc, c,d, counter, it, ii, flip;
	int valid = 1, next, o_valid, have_gone_back = 0; /*HGB: see BOT reverse */
	char *temp_ptr;
	unsigned long *new_ST;
	unsigned long long exponent = 0, exponent_2,exponent_3, exp_temp, box;
	unsigned long long *arr, word_v = 0, old_word_v = 0, temp = 0;
	new_ST = malloc(sizeof(long) * (alpha_size+1));
	len = k;
	/* Reset tables */
	memset(ST, 0, sizeof(long) * (alpha_size+1));
	memset(LOT, -1, alpha_size * sizeof(int));
	memset(BOT, 0, BOT_size/8);

	if (DEBUG) {
		printf("process(w==\"%s\")\n", word);
		printf("k==%d\n", k);
		printf("i==%d\n", i);
	}
	size = (int) (len / chars_per_long);
	if (len % chars_per_long != 0) {
		size++;
	}
	/* Temp: Convert string to long for V2. Queue will do this later
	   stored in base (alpha_size) */
	arr = string_to_longs(word, len);
	word_arr = arr;
	if (DEBUG_V2) {
		for (ii = 0; ii < size; ii++) {
			printf("arr[%d]: %llu\n", ii, arr[ii]);
		} 

	}
	/* Reconstruct ST, LOT, BOT */
	for (it = 0; it < size; it++) {
		if (DEBUG_V2) printf("it: %d, size: %d\n", it, size);
		word_v = arr[it];
		if (it < size - 1) { /* Non-last -> full long */
			count = chars_per_long;
			exponent = pow(alpha_size, chars_per_long-1);
		}
		else {	/* last long, could be < chars_per_long */
			count = len - (it * chars_per_long);
			exponent = pow(alpha_size, count-1);
		}
		old_word_v = word_v;
		exp_temp = exponent;
		if (DEBUG_V2) printf("count: %d\n", count);
		for (i = 0; i < count; i++) {
			/* Obtain letter from 'encoded' number */
			for (j = 0; j < alpha_size; j++) {
				if (DEBUG_V2) printf("LOT[%d]: %d\n", j, LOT[j]);
			}
			if (exponent == 0) {
				letter = word_v;
				exponent = 1;
			}
			else {
				letter = (word_v / exponent);
			}

			/* update BOT */
			assert(letter < alpha_size);
			lc = LOT[letter];	/* Compare to LOT */
			if (DEBUG_V2) printf("i: %d, letter: %d; lc = %d\n", i, letter, lc);
			if (lc != -1) {	/* Box occured */
			/* Only adding most closest box as only one that will
			appear when using suffix table */
				if (i + (it * chars_per_long) - lc +1 > alpha_size) { /* internal box exists */
					goto eloop;	/* TODO used to be continue; correct? */
				}
				/* Calculate box */
				if (DEBUG_V2) printf("exp_temp: %llu,\told_word_v: %llu\n", exp_temp, old_word_v);
				exponent_2 = exp_temp;
				temp = old_word_v;
				/* TODO: Possible bug? check check */
				if (it == ((lc / (chars_per_long+1)) - (lc % chars_per_long == 0))) {
					if (DEBUG_V2) printf("box in this long\n");
					/* Handle boxes in this long */
					lc %= chars_per_long;
					for (j = 1; j < lc; j++) {	/* remove left of box */
						temp %= exponent_2;
						exponent_2 /= alpha_size;
					}
					if (DEBUG_V2) printf("exponent: %llu\texponent_2: %llu\n", exponent, exponent_2);
					if (DEBUG_V2) printf("temp: %llu\tword_v: %llu\n", temp, word_v);
					//box = (temp - (word_v % exponent)) / exponent;
					box = temp / exponent;	/* remove right of box */
					if (DEBUG_V2) printf("box: %llu\n", box);
					assert(box <= BOT_size);
				}
				else {
					/* TODO: Test */
					if (DEBUG_V2) printf("spillover method\n");
					if (DEBUG_V2) { \
						printf("word: "); \
						for (j = 0; j < k; j++) { \
							printf("%c", word[j]); \
						} \
						printf("\nbox: %llu, temp %llu\n", box, temp); \
				} \
					/* handle boxes that spills over to previous long */
					CALC_SPILLOVER_BOX(lc, temp, exponent, exponent_2, exponent_3, word_arr,i, box, it);
				}
				if (DEBUG_V2) printf("box found. Value: %llu\n", box);
				/* flip bit */
				/* Move (box/8) bytes and change (box%8) bit */
				/* TODO: Can actually check box 1st as failsafe.  Quick */
				BOT += box / 8;
				temp_ptr = (char*) BOT;
				exponent_2 = 128;
				for (j = 0; j < box % 8; j++) {
					exponent_2 /= 2;
				}
				*temp_ptr ^= exponent_2;
				if (DEBUG_V2) printf("Updated BOT entry %llu\n", box);
				BOT -= box / 8;
			}
			eloop:
			/* update LOT */
			assert(letter <alpha_size && letter >= 0);
			LOT[letter] = i+1 + chars_per_long * it;

			/* Setup for next iteration */
			word_v %= exponent;
			exponent /= alpha_size;
		}
		word_v = old_word_v;
	}
	it--;	/* pointing to last long in array */
	/* Update suffix table (ST) */
	exponent = 1;
	for (i = 0; i < alpha_size; i++) {	/* exp ^ alpha_size*/
		exponent *= alpha_size; /*TODO check what this does exactly */
	 }
	exponent_2 = exponent;
	if (it == 0 || count >= alpha_size) { /* first long or no underflow risk */

		if (DEBUG_V2) printf("exp: %llu\tword_v: %llu\n", exponent, word_v);
		for (i = alpha_size; i>=0; i--) {
			ST[i] = word_v;
			if (DEBUG_V2) printf("added suffix %llu at index %d\n", word_v, i);
			word_v %= exponent;
			exponent /= alpha_size;
		}
	} else {	/* Subsequent longs with underflow risk */
		exponent = pow(alpha_size, count);
		if (DEBUG_V2) printf("exponent: %llu\n", exponent);
		for (i = count-1; i >= 0; i--) { /* Current long */
			word_v %= exponent;
			exponent /= alpha_size;
			ST[i] = word_v;
			if (DEBUG_V2) printf("added suffix %llu at index %d\n", word_v, i);
		}
		/* TODO: Test */
		word_v = word_arr[it-1];
		exponent = pow(alpha_size, count);
		for (i = alpha_size - (count -1); i > count-1; i--) {	/* Previous long */
			ST[i + count-1] = ST[count-1] + (word_v % exponent) * exponent;
			exponent *= alpha_size;
			if (DEBUG_V2) printf("ST[%d] now %lu\n", i + count-1, ST[i+count-1]);
		}
	}
	if (DEBUG_V2) {
		for (j = 0; j < size; j++) {
			printf("arr[%d] == %llu\n", j, word_arr[j]);
		}
	}
	/* TODO: remove this once queueing issue is gone */
	if (deep_check(word,i) == INVALID) return;

	/* Iterations */
	i = k;
	count = i+1;
	counter = count;
	word_v = old_word_v;
	/* TODO: increase long count if needed */
	word_v *= alpha_size;
	//word_v--;
	if (DEBUG_V2) printf("word_v == %llu\n", word_v);
	word[i] = '/';
	while (i >= k) {
		if (DEBUG) {
			printf("exploring position i==%d w==", i);
			for (ii = 0; ii <= i; ii++) {
				printf("%c", word[ii]);
			} 
			printf("\n");
			printf("word_v at start: %llu\n", word_v);
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
			/* TODO: Revert back to previous long if needed */
			long_count = count / (chars_per_long);
			if (count % chars_per_long == 0 && count != 0) {
				long_count--; /* if last char, shouldn't increment count */
			}
			arr[long_count] = word_v;
			if (DEBUG_V2) printf("i--: long_count == %d, count == %d\n", long_count, count);
			if (count % (chars_per_long) != 1 || count == 1) { /* Todo: used to be > 0 */
				word_v /= alpha_size;
			} else {
				if (DEBUG_V2) printf("count: %d reverting to arr[%d]: %llu\n", count, long_count-1, arr[long_count-1]);
				word_v = arr[long_count-1];
			}
			count--;
			
			if (DEBUG_V2) printf("word_v / alpha_size. now %llu\n", word_v);
			/* Update LOT */
			j = word_v % alpha_size;
			c = i % (chars_per_long+1);
			temp = word_v / alpha_size;
			assert(j < alpha_size && j >= 0);
			REVERT_LOT(j, c, temp, it, word_arr, d, end_l);
			end_l:
			if (DEBUG_V2) {
				printf("LOT:\n");
				for (j = 0; j < alpha_size; j++) {
					printf("[%d]: %d\n",j, LOT[j]);
				}
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
				REVERT_BOT_ENTRY(d, j, word_v, exp_temp, exponent_2, letter, temp_ptr, temp, box, word_arr);
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
			/* Handling 'underflow' to multiple longs  in suffix table */
			word_v /= alpha_size;
			if (DEBUG_V2) printf("word_v / alpha_size. now %llu\n", word_v);
			new_ST[alpha_size] = word_v % exp_temp;
			j = (i+1) % chars_per_long;
			if (DEBUG_V2) printf("i+1 mod chars_per_long: %d\n", j);
			if (j < alpha_size+1 && i / chars_per_long > 0) { /* previous long needed */
				if (DEBUG_V2) printf("rollover revert\n");
				/* TODO: Check if below accesses correct entry */
				temp = word_arr[(i+1) / chars_per_long - 1];
				if (DEBUG_V2) printf("String form of temp: %s\n",longs_to_string(&temp,32));
				j = alpha_size+1 - j;
				exp_temp = pow(alpha_size, j);
				d = ((temp % exp_temp) / (exp_temp/alpha_size));
				new_ST[alpha_size] = new_ST[alpha_size-1] + d * pow(alpha_size, alpha_size);
				if (DEBUG_V2) printf("letter: %d\n", d);
				if (DEBUG_V2) printf("temp = %llu, j = %d, exp_temp = %llu, new_st[al] = %ld\n", temp, j , exp_temp, new_ST[alpha_size]);
				/* TEST */
 			}
			/* TODO: Need to switch to old word_v here? */
			long_count = count / (chars_per_long);
			if (count % chars_per_long == 0 && count != 0) {
				long_count--; /* if last char, shouldn't increment count */
			}
			arr[long_count] = word_v;
			if (DEBUG_V2) printf("<>long_count: %d, count : %d\n", long_count, count);
			if (count % (chars_per_long) == 1 && count != 1 && long_count > 0) { /* TODO: -1? */
				if (DEBUG_V2) printf("\n\nreverting to arr[%d] as count == %d\n\n\n", long_count-1, count);
				word_v = arr[long_count-1];
			}
			count--;

			if (DEBUG_V2) {
				printf("word_v == %llu\n", word_v);
				printf("new_st[%d] == %ld\n", alpha_size,new_ST[alpha_size]);
				printf("\n");
			}
			for (j = 0; j <= alpha_size; j++){
				ST[j] = new_ST[j];
			}

			/* Update LOT */
			j = word_v % alpha_size;
			c = i % (chars_per_long+1);
			temp = word_v / alpha_size;
			assert(j < alpha_size && j >= 0);
			REVERT_LOT(j, c, temp, it, word_arr, d, end_l2);
			end_l2:
			if (DEBUG_V2) {
				printf("LOT:\n");
				for (j = 0; j < alpha_size; j++) {
					printf("[%d]: %d\n",j, LOT[j]);
				}
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
				REVERT_BOT_ENTRY(d, j, word_v, exp_temp, exponent_2, letter, temp_ptr, temp, box, word_arr);
			}
			/* Update LOT entry TODO CHECK IF LOCATION RIGHT */
			j = word_v % alpha_size;
			c = i % (chars_per_long+1);
			temp = word_v / alpha_size;
			assert(j < alpha_size && j >= 0);
			REVERT_LOT(j, c, temp, it, word_arr, d, end_l3);
			end_l3:
			word_v++;
			have_gone_back = 0;
			if (DEBUG_V2) printf("word_v++. now %llu\n", word_v);
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
		if (DEBUG_V2) printf("next: %d\n", next);
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
			if (DEBUG_V2) printf("box value: %llu\n", box);
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
			BOT -= box/8;
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
			printf("word_v == %llu\n", word_v);
			printf("o_valid == %d, valid == %d, 0_valid == valid: %d\n",
				o_valid, valid, o_valid==valid);
		}
		if (STATUS && o_valid != valid) {
			printf("Warning: o_valid == %d, valid == %d, 0_valid == valid: %d\n",
				o_valid, valid, o_valid==valid);
			exit(0);
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
			/* increase long count if needed */
			long_count = count / (chars_per_long);
			if (count % chars_per_long == 0 && count != 0) {
				long_count--; /* if last char, shouldn't increment count */
			}
			/* TODO: Used to be long_count = count / (chars_per_long+1);*/
			if (DEBUG_V2) printf("arr[%d] updated to %llu\n", long_count, word_v);
			arr[long_count] = word_v;
			if (DEBUG_V2) printf("i == %d, it == %d, count == %d, long_count == %d\n", i, it, count, long_count);
			//if (count % chars_per_long < chars_per_long-1) { /* Todo: fix to remove -1 */
			if (count % chars_per_long != 0 || count == 0) {
				word_v *= alpha_size;
			} else {
				if (DEBUG_V2) printf("increasing long count & realloc\n");
				arr[long_count] = word_v;
				if (DEBUG_V2) printf("arr[%d] == %llu\n", long_count, word_v);
				word_v = 0;
				word_arr = realloc(word_arr, sizeof(unsigned long long) * (long_count + 1));
				arr[long_count + 1] = 0;
			}
			count++;
			if (DEBUG_V2) printf("word_v * alpha_size. now %llu\n", word_v);
			assert(i <= max_word_size[alpha_size]);
			if (i > max_length) {
				max_length = i;
				word[i] = '\0';
				strcpy(max_word, word);
			}
			/* Setup for next run */
			word[i] = '/';
		} else if (DEBUG) {
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

/*************************** Utility Functions ********************************/

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

/* len - length of 'string' (# characters) */
int deep_check_l(unsigned long long *word_arr, int len) {
	unsigned long long word_v, old_word_v, exponent, exp_temp, box, temp, exponent_2, exponent_3;
	int size = (int) (len / chars_per_long); /* # longs */
	int i, j, k, c, valid = VALID, count, letter, lc;
	char *temp_ptr;
	void *BOT2 = malloc(BOT_size);
	int *LOT2 = malloc(sizeof(int) * alpha_size);
	memset(BOT2, 0, BOT_size/8);
	memset(LOT2, 0, alpha_size * sizeof(int));

	if (len % chars_per_long != 0) {
		size++;	/* Extra long if needed */
	}
	for (i = 0; i < size; i++) {
		word_v = word_arr[i];
		if (i < size - 1) { /* Non-last -> full long */
			count = chars_per_long;
			exponent = pow(alpha_size, chars_per_long-1);
		}
		else {	/* last long, could be < chars_per_long */
			count = len - (i * chars_per_long);
			exponent = pow(alpha_size, count-1);
		}
		old_word_v = word_v;
		exp_temp = exponent;
		for (j = 0; j < count; j++) {
		/* Get letter */
			if (exponent == 0) {
				letter = word_v;
				exponent = 1;
			} 
			else {
				letter = (word_v / exponent);
			}
			assert(letter < alpha_size);
			/* Compare to LOT */
			lc = LOT2[letter];
			if (lc != 0) { /* box occured */
				if (i - lc > alpha_size) { /* internal box exists */
					goto eloop;
				}
			
				/* Calculate box */
				exponent_2 = exp_temp;
				temp = old_word_v;
				if (((i * chars_per_long + j) / chars_per_long)
						== (lc / chars_per_long)) {
					/* Handle boxes in this long */
					lc %= chars_per_long;
					for (k = 1; k < lc; k++) {	/* remove left of box */
						temp %= exponent_2;
						exponent_2 /= alpha_size;
					}
					box = temp / exponent;	/* remove right of box */
					assert(box <= BOT_size);
				}
				else {
					/* TODO: Test */
					/* handle boxes that spills over to previous long */
					CALC_SPILLOVER_BOX(lc, temp, exponent, exponent_2, exponent_3, word_arr,i, box,i);
				}
				/* Check if BOT entry != 0. If so, return 0. else update BOT */
					/* Move (box/8) bytes and change (box%8) bit */
				if (DEBUG_V2) printf("box value: %llu\n", box);
				BOT2 += box / 8;
				exponent_2 = 128;
				for (k = 0; k < box % 8; k++) {
					exponent_2 /= 2;
				}
				temp_ptr = (char*) BOT2;
				c = (*temp_ptr) & (exponent_2);
				if (DEBUG_V2) printf(">c == %d\n", c);
				if (c == 0) {
					if (DEBUG_V2) printf(">updating box\n");
					*temp_ptr ^= exponent_2;
				} else {
					valid = INVALID;
					goto end;
				}
				BOT2 -= box/8;

			}
			eloop:
			/* Update LOT */
			LOT2[letter] = j+1;
			/* Setup for next iteration*/	
			word_v %= exponent;
			exponent /= alpha_size;
		}
	}
	end:
	free(BOT2);
	return valid;
}
/* Converts an 'encoded string' (long format) to a string */
char *longs_to_string(unsigned long long *arr, int len) {
	char letter, *string = malloc(sizeof(char) * (len+1));
	int i, j, count;
	unsigned long long word_v, exponent;
	int size = (int) (len / chars_per_long);	/* # longs */
	if (len % chars_per_long != 0) {
		size++;
	}
	for (i = 0; i < size; i++) { /* For each long */
		word_v = arr[i];
		if (i < size - 1) { /* Non-last -> full long */
			count = chars_per_long;
			exponent = pow(alpha_size, chars_per_long -1);
			/* TODO: Check */
		}
		else {	/* last long, could be < chars_per_long */
			count = len - (i * chars_per_long);
			exponent = pow(alpha_size, count - 1);
			/* TODO: Check */
		}
		for (j = 0; j < count; j++) { /* Get chars from long */
			letter = alpha[0] + word_v / exponent;
			assert(letter >= alpha[0] && letter <= alpha[alpha_size]);
			string[i * chars_per_long + j] = letter;
			word_v %= exponent;
			exponent /= alpha_size;
		}
	}
	string[len] = '\0';
	return string;
}

unsigned long long *string_to_longs(char *string, int len) {
	int i, j, count, size = (int) (len / chars_per_long);
	unsigned long long word_v, exponent, *arr;

	if (len % chars_per_long != 0) {
		size++;
	}
	arr = malloc(sizeof(unsigned long long) * size);
	memset(arr, size, sizeof(unsigned long long));

	for (i = 0; i < size; i++) { /* For each long */
		word_v = 0;
		if (i < size - 1) { /* Non-last -> full long */
			count = chars_per_long;
			exponent = pow(alpha_size, chars_per_long -1);
			/* TODO: Check */
		}
		else {	/* last long, could be < chars_per_long */
			count = len - (i * chars_per_long);
			exponent = pow(alpha_size, count - 1);
			/* TODO: Check */
		}
		for (j = 0; j < count; j++) { /* Get chars from long */
			word_v += (string[i * chars_per_long + j] - alpha[0]) * exponent;
			exponent /= alpha_size;
		}
	arr[i] = word_v;
	}
	return arr;
}