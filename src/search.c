/******************************************************************************
 *	Finds the longest string that has no repeated boxes. Utilises MPI to run
 *	multiple workers at once. Relies on a queue-based system to handle jobs
 *	between workers.
 * 
 *  Setup: Create a 'bin' and 'saves' directory in root of repo
 * 
 *	compile:	mpicc -o ../bin/search.o queue.c search.c -lm
 *	run:		mpirun -np X ./search.o	N <filename>
 		where X is number processes,
	 	N alphabet size between 2 and 9 (default if not specified)
		OPTIONAL: <filename> is the name of the queue file to restore from

 * WARNING: Since one process is used for controlling the rest, do not specify
 * more than the available number of cores. Else the program can lock (ROOT
 * process essentially waits for itself) 
 * 
 * If attempting to use debugging flags, strongly recommend
 * piping output to file.
 * 
 * 	Original scaffolding of sequential process provided by Jaco Geldenhuys
 * 
 *	@author C. R. Zeeman (caleb.zeeman@gmail.com)
 *	@version 2.5
 *	@date 2019-12-17
 *****************************************************************************/
/* Includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <math.h>
#include <limits.h>
#include <time.h>
#include "queue.h"
#include "search.h"

/************************** Definitions and macros **************************/
/*--- CONFIGURATIONS  ---*/
/* Enables debug output. Primarily related to overhead/MPI calls */
#define DEBUG (0)

/* Enables debug output. Primarily for version 2 (O(1) search) in process */
#define DEBUG_V2 (0)

/* Enables extra failsafes to ensure answer is definitely valid.
 * Some of these slow down processing speed, so once it has been
 * tested a few times, these should be disabled.
 * Recommended mode: 1.
 * 0 - No failsafes
 * 1 - Basic assertions and surrounding ifs
 * 2 - Above and enqueue deep_check_l
 * 3 - Above and string deep_check & assertion built into deep_check_l
 * 
 * WARNING: There's still a bug from the original sequential code
 * where it would rarely enqueue an invalid pattern. There's two solutions
 * to this:
 * a) set 2 above and the deep_check will prevent invalid patterns
 * from being enqueued. This is quite a bit slower.
 * b) Set 0 or 1 above, and a second check will prevent an invalid pattern
 * from being processed. Faster, however the queue will then also contain
 * possibly invalid patterns. Will be picked up before processing though.
 */
#define FAILSAFE_MODE (1)

#if FAILSAFE_MODE == 0
/* Disables assert statements if mode == 0 */
#define NDEBUG
#endif

/* Prints status messages (alphabet size being used, etc) */
#define STATUS (1)

/* Prints status messages related to stalling processes
   (Not needed unless wanting to see when a process is idling) */
#define STATUS_2 (0)

/* Prints maximum whenever it is updated / a new max is found */
#define PRINT_MAX (1)

/* -- SAVE FILE -- */
/* Toggles whether queue is saved to files */
#define STORE_QUEUE (1)

/* How often to store the queue. Stores queue every X minutes */
/* Note: saving/loading can take a minute or two in the worst case,
 so generally don't have this time too rapid. 2 hours works well (120) */
#define QUEUE_TIME (120)

/* Number of different save files for queue*/
#define COPIES_KEPT	(9)

/* Prints whenever the queue is saved to or reverted from a file */
#define PRINT_SAVE_LOAD (1)

/* Whether or not to store the maximum found so far in the save file too */
#define STORE_MAX (1)

/*-- SEARCH DEPTH CONFIGURATIONS -- */
	/* recommended is (9, 3) */
/* Depth to search per job */
#define STANDARD_SEARCH_DEPTH (9)

/* Depth for initial run by ROOT */
/* Note: Ensure that this value is high enough to generate queue entries
   for all child processes. Likewise, setting this too high results
   in higher setup time initially. Anywhere between 3-5 seems ideal */
#define INITIAL_SEARCH_DEPTH (3)

/* Default alphabet size, if none was specified */
#define STANDARD_ALPHABET_SIZE (3)

/* Observations on config:
  (6,3) - (9,3) seems to have very little memory usage
  (10,3) memory usage starts increasing slowly
  (11,3) takes forever to find first max strings. Broad search? */

/*--- END OF CONFIG ---*/

/* Signals for valid and invalid */
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

/* Macros to avoid function calls where possible */

/* Revert BOT entry, if any */
#define REVERT_BOT_ENTRY(d, j, word_v, exp_temp, exponent_2, letter, temp_ptr,temp, box, word_arr) do { \
	unsigned long long preamble = 0; \
	exp_temp = alpha_size;	\
	d = count-1; /*keeps track of when to switch to previous long's values*/ \
	temp = word_v; \
	letter = (word_v % alpha_size);	\
	if (d % chars_per_long == 0 && count % chars_per_long != 0) box = letter; \
	else box = 0; \
	exponent_2 = 1; \
	if (DEBUG_V2) printf("word_v: %llu\n", word_v);	\
	if (DEBUG_V2) printf("BOT_revert: letter: %d\n", letter);	\
	for (j = 1; j <= alpha_size; j++) { /* find same char */	\
		if (DEBUG_V2) printf("d: %d, count %d\n", d, count); \
		if (d % chars_per_long == 0 && count % chars_per_long != 0) { \
			if (DEBUG_V2) printf("swapped to prev. long\n"); \
			preamble = box; \
			exponent_2 = exp_temp; \
			exp_temp = 1; /* Alpha_size? */ \
			temp = word_arr[count / chars_per_long -1]; \
		} \
		if (DEBUG_V2) printf("BOT_revert: searching: %llu\n", ((temp % (exp_temp*alpha_size)) /exp_temp));	\
		assert (exp_temp > 0); \
		box = (temp % (exp_temp*alpha_size)) * exponent_2 + preamble; /* box value TODO CHECK*/	\
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

/* Calculates the value of a box that has 'spilled over' to a previous long*/
/*e.g. 32 chars per long, and box exists from i=31 to i=33 */
#define CALC_SPILLOVER_BOX(lc, temp, exponent, exponent_2, exponent_3, word_arr,i, box, it) do { \
	int x = 0; \
	if (lc != 0 && lc % (chars_per_long) == 0) lc = chars_per_long; \
	else lc %= chars_per_long; \
	/* calculate box from end of long to current position */ \
	if (DEBUG_V2) printf("spillover box. it: %d\n", it); \
	box = temp / exponent; /* remove right of box */ \
	/* calculate box segment in previous long */ \
	temp = word_arr[it-1]; \
	exponent_2 = 1; \
	exponent_3 = 1; \
	for (x = 0; x < i+1; x++) { \
		exponent_3 *= alpha_size; /*pow(alpha_size, i+1)*/ \
	} \
	for (x = 0; x < (chars_per_long - lc + 1); x++) { \
		exponent_2 *= alpha_size; /*pow(alpha_size, chars_per_long - lc + 1)*/ \
	} \
	if (DEBUG_V2) printf("exponent: %llu\texponent_2: %llu\n", exponent, exponent_2); \
	if (DEBUG_V2) printf("temp: %llu\n", temp); \
	if (exponent_2 == 0) { /* Temp. Del this if */ \
		printf("lc == %d, exponent_2 == %llu\n", lc, exponent_2); \
		printf("count: %d\n", count); \
	} \
	assert(exponent_2 > 0); \
	temp %= exponent_2; /* Remove left of box */\
	/* add together */ \
	if (DEBUG_V2) printf("box before add: %llu\n", box); \
	if (DEBUG_V2) printf("temp: %llu\n", temp); \
	if (DEBUG_V2) printf("exponent_3: %llu\n", exponent_3); \
	box += temp * exponent_3; \
	if (DEBUG_V2) printf("new box: %llu\n", box); \
	if (box > BOT_size) { \
		/* If assert fails, print surrounding variables */ \
		printf("box: %llu, lc: %d\n", box, lc); \
		printf("word_v: %llu, temp: %llu\n", word_arr[i], temp); \
		assert(box <= BOT_size); \
	} \
	if (DEBUG_V2) printf("end of spillover method\n"); \
} while(0)

/* Revert LOT entry to it's previous occurance of the letter */
#define REVERT_LOT(j, c, temp, it, word_arr, d, label_name) do { \
	/* Both of these are needed */ \
	c += i / chars_per_long; /* Adjustment for over-mods */ \
	c %= (chars_per_long+1); /* Adjustment for adjustment */ \
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
	/* if still not found, Search previous longs */ \
	if (LOT[j] == -1) { \
		if (DEBUG_V2) printf("searching previous longs for %d. ii: %d\n",j, (count / chars_per_long) - 1); \
		for (x = (count / chars_per_long) -1; x >= 0; x--) { \
			if (DEBUG_V2) printf("word_arr[%d] == %llu\n", x, word_arr[x]); \
			temp = word_arr[x]; \
			c = chars_per_long; /* TODO: -1? */ \
			while (c > 0) { \
				d = temp % alpha_size; \
				if (DEBUG_V2) printf("c: %d d: %d, d==j: %d, temp: %llu\n",c, d, d==j,temp); \
				if (d == j) { \
					if (DEBUG_V2) printf("d==j\n"); \
					LOT[j] = c + (x * chars_per_long); \
					goto label_name; \
				} \
				temp /= alpha_size; \
				c--; \
			} \
			if (DEBUG_V2) printf("out of loop\n"); \
		} \
	} \
	if ( STATUS && LOT[j] == -1 && i > alpha_size) { \
		printf("warning: LOT[%d] == -1 at i == %d\n", j, i); \
	} \
} while (0)

/* 	Sets count to long_count's value, and long_count to the
	number of longs needed to fit in said count */
#define UPDATE_COUNT(count, long_count) do { \
	count = long_count;	\
	long_count = (int) (count / chars_per_long); \
	/* if count is midway through a long, increase the count */ \
	if (count % chars_per_long != 0) { \
		long_count++; \
	} \
	assert(long_count <= max_longs_needed); \
} while (0)
/*************************** global variables ******************************/

/* Alphabet used: the alpha_size'd subset of alpha[] */
int alpha_size = STANDARD_ALPHABET_SIZE;
char alpha[] = {'0','1','2','3','4','5', '6', '7', '8', '9'};
char last;	/* Last char to process in alphabet. Set in main */

/* Theoretical longest word per alphabet size. l(n) = n(l(n-1) +2) */
int max_word_size[] = {0, 2, 8, 30, 128, 650, 3912, 27398, 219200, 1972818};
/* Properties of max word found */
int max_length = 0, max_longs_needed;
char *max_word;
unsigned long long *word_l_arr;
unsigned long long *max_word_l;
/* Number of 'characters' that can fit in a long long */
int chars_per_long;

int my_id;		/* Process's ID; set in main */

/* V2 */
unsigned long *ST;	/* Suffix_table for current. Table 1 */
int suffix_depth;
int *LOT;	/* Last_occurance_table; where alpha char last occured. Table 2 */
void *BOT; /* Box_occurance_table. Bit table of if box has occured. Table 3 */
unsigned long long BOT_size;	/* Number of bits in BOT */

/* File-saving */
char file_name[] = "../saves/n^-_.save"; /* '^', '_' replaced by function */
int swap_index = 12;	/* Index of '_' above */
int n_index = 10; /* Index of '^' above */
int current_save_number = 0;	/* For naming save files */

/*********************** function prototypes *******************************/
/* Long variants - O(1) search */
void process (int depth, unsigned long long *word_arr, int len);
int deep_check_l(unsigned long long *word_arr, int len);		
char *longs_to_string(unsigned long long *arr, int len);
unsigned long long *string_to_longs(char *string, int len);

void store_queue_to_file(char *save_name);
void restore_queue_from_file(char *save_name);

/* String variants */
int box_valid(char *word, int len);			
int deep_check(char *string, int len);	

/**************************** functions ***********************************/

/* ROOT_PROCESS handles the queue, while all other processes do the
 * calculations. Communication is done via tags to denote message types.
 * Communcation cycle as follows:
 * ROOT does one iteration and adds items to queue, then sends one job to each
 * workers. The workers process the job, and once reach specified depth, sends a
 * signal back to the ROOT to recieve the item to be queued, followed by the
 * pattern itself. Once complete, the worker sends the maximum size it found back
 * to ROOT. if this is longer than previous maximum found by all processes, ROOT
 * requests the pattern and stores it. It then sends the next workload to the
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
	int desired_size, stop = 0, count, i = 0, j = 0;
	unsigned long long partial_sum = 0;	/* Calculating box size */
	unsigned long long queue_size = 0;	/* Number entries in queue */
	/* Variables used for MPI communication */
	MPI_Status status;
	int ierr = 0, nbr_procs = 0, rec_count, send_count, an_id;
	char *w; /* used for temporary conversions and intial setup */
	/* Keeps track of how many processes have been successfully stopped */
	int active_processes;
	/* For when queue is empty and continuing */
	int waiting_count = 0, *waiting_for_work;
	/* For saving */
	clock_t start_time, current_time;
	double elapsed_time;
	int save_needed = 0, stall_count = 0;

	/* MPI setup */
	ierr = MPI_Init(&argc, &argv);
	ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	ierr = MPI_Comm_size(MPI_COMM_WORLD, &nbr_procs);

	/* Command line arguement details */
	if (argc == 2 || argc == 3) {
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
	}
	else {
		if (STATUS && my_id == ROOT_PROCESS) {
			printf("Using default alphabet size: %d\n", STANDARD_ALPHABET_SIZE);
		}
	}
	if (nbr_procs <= 1) {
		if (my_id == ROOT_PROCESS) {
			printf("Usage: \"mpirun -np X %s ALPHA_SIZE\" where X is number processess > 1\n",
				argv[0]);
		}
		goto end; /* Skips unneeded frees */
	}
	/* constants and variable setup */
	desired_size = max_word_size[alpha_size];
	last = '/' + alpha_size;
	max_word = (char*) malloc(sizeof(char) * (desired_size + 1));  /* +1 for '\0' */
	max_word[0] = '\0';
	chars_per_long = sizeof(long long int) * 8 / (log(alpha_size) / log(2));
	file_name[n_index] = '0' + alpha_size;

	/* Number of longs required to store maximum length pattern */
	max_longs_needed = (int)(desired_size / chars_per_long);
	if (desired_size % chars_per_long != 0) max_longs_needed++;
	max_word_l = malloc(sizeof(long long) * max_longs_needed);
	
	/* Setup for empty queue situation */
	waiting_for_work = malloc(sizeof(int) * nbr_procs);
	for (i = 0; i < nbr_procs; i++) {
		waiting_for_work[i] = 0;
	}
	for (i = 0; i < max_longs_needed; i++) {
			max_word_l[i] = 0;
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
	/*Padding to next byte */
	BOT_size += BOT_size % 8;
	BOT = malloc(BOT_size);

	/* Root/Master */
	if (my_id == ROOT_PROCESS) {
		start_time = clock();
		if (STATUS) {
			printf("BOT_size: %llu\n", BOT_size);
			printf("char per long: %i\n", chars_per_long);
			printf("Longs needed for max: %d\n", max_longs_needed);
		}
		if (argc == 3) { /* Restore from file*/
			restore_queue_from_file(argv[2]);
			if (get_queue_size_l() == 0) {
				printf("Nothing was able to be restored from file\n");
				return EXIT_FAILURE;
			}
		}
		else { 
			/* Initial run to create branches for workers */
			/* Initial start is 1 of each alphabet character */
			for (i = 0; i < alpha_size; i++) {
				w[i] = alpha[i];
			}
			w[i] = '\0';
			word_l_arr = string_to_longs(w, alpha_size);
			free(w);

			process(INITIAL_SEARCH_DEPTH, word_l_arr, alpha_size);

			free(word_l_arr);
		}
		/*	TOOL: To manually confirm if a pattern is valid,
			uncomment this and replace str with the pattern
			and fill in size. Will do both long and str check */
		//char str[] = "0123301212303223101302312001313201120213012310321302132301321031203210232013203123020302103201230102301203102310213";
		//int len_of_str = 115;
		//printf("str: %s\nisvalid: %d and %d\n", str, deep_check_l(string_to_longs(str,len_of_str),len_of_str), deep_check(str, len_of_str));
		word_l_arr = malloc(sizeof(long long) * max_longs_needed);
		MPI_Barrier(MPI_COMM_WORLD);

		/* Send initial instructions to each process */
		for (an_id = 1; an_id < nbr_procs; an_id++) {
			if (!dequeue_l(word_l_arr, &send_count)) {
				printf("<%d> dequeue failed at %d\n", my_id, an_id);
				break;	/* Failsafe */	
			}
			if (DEBUG) {
				printf("<%d>Dequeued [0]: %llu with length %d\n", my_id, word_l_arr[0], send_count);
				printf("<%d>queuesize: %lld\n", my_id, get_queue_size_l());

			}
			ierr = MPI_Send(&send_count, 1, MPI_INT, an_id, NEW_WORK_TAG,
				MPI_COMM_WORLD);
			UPDATE_COUNT(count, send_count);
			if (DEBUG) printf("<0>send_count: %d, count: %d\n", send_count, count);
			ierr = MPI_Send(word_l_arr, send_count, MPI_LONG_LONG, an_id, WORK_TAG, 
				MPI_COMM_WORLD);
		}

		/* Communcation cycle */
		/* Wait for size found from process; if size > max, ask process to send
		 * string. Else, if need to stop, send stop signal. else send next
		 * workload signal followed by workload */
		queue_size = get_queue_size_l();
		while (1) {
			if (STORE_QUEUE) {
				current_time = clock();
				elapsed_time = ((double) (current_time - start_time)) / CLOCKS_PER_SEC / 60;
				if (elapsed_time > QUEUE_TIME) {
					/* If a save is needed, continue as normal. Whenever
					   a process sends back a maximum (i.e. done), stall 
					   the process. Once all processes stall (i.e. queue 
					   has all entries), save it, and then resume 
					   each process by giving it new work */
					if (STATUS_2 && save_needed == 0) {
						printf("save needed\n");
					}
					save_needed = 1;
				}
			}
			if (queue_size <= 0 || max_length >= desired_size) {
				/* Reached goal. Either explored all, or found longest */
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
					/* Request longest string*/
					MPI_Send(&rec_count, 1, MPI_INT, an_id, SEND_MAX_TAG,
						MPI_COMM_WORLD);
					UPDATE_COUNT(count, rec_count);
					if (DEBUG) {
						printf("<%d>New max found, requesting string.\t", my_id);
						printf("size %d, rec-count %d\n", max_length, rec_count);
					}
					MPI_Recv(max_word_l, rec_count, MPI_LONG_LONG, an_id,
						MAX_STRING_TAG, MPI_COMM_WORLD, &status);
					free(max_word);
					max_word = longs_to_string(max_word_l, max_length);
					if (PRINT_MAX) {
						printf("New max: %d -- %s\n", max_length, max_word);
						//printf("arr[0]: %llu\n", max_word_l[0]);
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
				} 
				else if (save_needed == 1) {
					/* Empty out processes into queue before saving */
					stall_count++;
					waiting_for_work[an_id] = 1;
					if (STATUS_2) printf("process %d stalled for save\n", an_id);
					if (stall_count == active_processes) {
						if (STATUS_2) printf("saving file\n");
						store_queue_to_file(file_name);
						restore_queue_from_file(file_name);
						save_needed = 0;
						waiting_count = stall_count;
						stall_count = 0;
						start_time = clock();
					}
				} /* ...or send the next work available in queue */
				else if (!dequeue_l(word_l_arr, &send_count)) {
					/* Queue temporarily empty;
					 wait until next enqueue request to send data to process
					 This should never run (unless queuesize==0 check earlier
					 is removed) */
					if (STATUS_2) printf("stalled process %d due to empty queue\n",
						an_id);
					assert(an_id < nbr_procs);
					waiting_for_work[an_id] = 1;
					waiting_count++;
				}
				else {
					if (DEBUG) {
						printf("(n)Dequeued [0]: %llu with length %d\n", word_l_arr[0], send_count);
						printf("queuesize: %lld\n", get_queue_size_l());
					}
				//printf("count: %d, queue_size: %llu. Dequeued: %llu %llu %llu %llu\n",send_count, queue_size,word_l_arr[0],word_l_arr[1],word_l_arr[2],word_l_arr[3]);
					ierr = MPI_Send(&send_count, 1, MPI_INT, an_id,
						NEW_WORK_TAG, MPI_COMM_WORLD);
					UPDATE_COUNT(count, send_count);
					ierr = MPI_Send(word_l_arr, send_count, MPI_LONG_LONG, an_id, WORK_TAG, 
						MPI_COMM_WORLD);
				}
			}
			/* Need to queue a pattern */
			else if (status.MPI_TAG == EXPECT_QUEUE_TAG) {
				if (stop) { /* New work  TODO rework */
					stop = 0;
				}
				count = rec_count;
				rec_count = (int) (count / chars_per_long);
				if (count % chars_per_long != 0) {
					rec_count++;
				}
				assert(rec_count <= max_longs_needed && rec_count > 0);
				MPI_Recv(word_l_arr, rec_count, MPI_LONG_LONG, an_id, QUEUE_TAG,
					MPI_COMM_WORLD, &status);
				if (DEBUG) printf("<%d> asked to enqueue. word_l_arr[0]: %llu, size: %d, reccount: %d\n",
					my_id, word_l_arr[0], count, rec_count);
				enqueue_l(word_l_arr, count);
			}
			queue_size = get_queue_size_l();
			/* Send work to any stalled processes if available
			   This either occurs due to a once-empty queue, or due to a save */
			if (waiting_count > 0 && queue_size >= waiting_count) {
				for (i = 0; i < nbr_procs; i++) {
					if (waiting_for_work[i] == 1) {
						dequeue_l(word_l_arr, &send_count);
						if (DEBUG) {
							printf("(n)Dequeued [0]: %llu with length %d\n", word_l_arr[0], send_count);
							printf("queuesize: %lld\n", get_queue_size_l());
						}
						an_id = i;
						if (STATUS_2) printf("process %d restarted\n", an_id);
						ierr = MPI_Send(&send_count, 1, MPI_INT, an_id,
							NEW_WORK_TAG, MPI_COMM_WORLD);
						UPDATE_COUNT(count, send_count);
						ierr = MPI_Send(word_l_arr, send_count, MPI_LONG_LONG, an_id, WORK_TAG, 
							MPI_COMM_WORLD);
					}
				}
				waiting_count = 0;
				queue_size = get_queue_size_l();
				/*TODO: If queuesize=0 and all processes waiting, send stop
				signal to all. This should only occur on a bug, so not urgent */
			}
		}
		//free(max_word); /* Needed? TODO*/
		max_word = longs_to_string(max_word_l,max_length);

		/* print status */
		printf("\n*********** FINISHED *********\n");
		printf("max == %d\n", max_length);
		printf("max word == %s\n", max_word);
		printf("queue size == %lld\n", get_queue_size_l());
		printf("max valid? %d and %d\n", deep_check(max_word, max_length),
		deep_check_l(max_word_l,max_length));
		printf("******************************\n\n");
	}
	
	/* Worker/Slave */
	else {
		MPI_Barrier(MPI_COMM_WORLD);
		word_l_arr = malloc(sizeof(long long) * max_longs_needed);
		for (i = 0; i < max_longs_needed; i++) {
			word_l_arr[i] = 0;
		}

		/* Continue until a message is recieved stating you must stop */
		while (status.MPI_TAG != STOP_TAG) {
			/* Get next command */
			ierr = MPI_Recv(&rec_count, 1, MPI_INT, ROOT_PROCESS,
				MPI_ANY_TAG, MPI_COMM_WORLD, &status);

			/*if command is to send max string back, do that */
			if (status.MPI_TAG == SEND_MAX_TAG) {
				send_count = max_length;
				UPDATE_COUNT(count, send_count);
				if (DEBUG) printf("<%d>: sending max of length %d, send_count %d\t",
					my_id,max_length,send_count);
				if (DEBUG) printf("value[0]: %llu\n", max_word_l[0]);
				if (DEBUG) printf("value[1]: %llu\n", max_word_l[1]);
				ierr = MPI_Send(max_word_l, send_count, MPI_LONG_LONG,
					ROOT_PROCESS, MAX_STRING_TAG, MPI_COMM_WORLD);
			}

			 /* else if command is to get new work, get new work string */
			else if (status.MPI_TAG == NEW_WORK_TAG) {
				UPDATE_COUNT(count, rec_count);
				ierr = MPI_Recv(word_l_arr, rec_count, MPI_LONG_LONG, ROOT_PROCESS,
					WORK_TAG, MPI_COMM_WORLD, &status);
				assert(rec_count <= max_longs_needed);
				if (DEBUG) {
					printf("process %d recieved [0]: %llu with length %d\n",
						my_id, word_l_arr[0], count);
					if (count > chars_per_long) {
						printf("[1] == %llu\n", word_l_arr[1]);
					}
				}
				if (DEBUG) {
					w = longs_to_string(word_l_arr, count);
					printf("beginning process. w: %s\n", w);
					free(w);
					printf("word_l_arr[0]: %llu, count: %d\n", word_l_arr[0], count);
					printf("word_l_arr[1]: %llu\n", word_l_arr[1]);
				}
				process(STANDARD_SEARCH_DEPTH, word_l_arr, count);		;
				//printf("<%d> finished processing, sending max back\n", my_id);

				/* Send local max length back */
				ierr = MPI_Send(&max_length, 1, MPI_INT, ROOT_PROCESS,
					SENDING_MAX_TAG, MPI_COMM_WORLD);
			}
		}
	}

	free(max_word);
	free(waiting_for_work);
	free(ST);
	free(LOT);
	free(BOT);
	free(word_l_arr);
	end: 
	ierr = MPI_Finalize();
	if (STATUS) printf("***** PROCESS %d STOPPED ******\n", my_id);
	return EXIT_SUCCESS;
}
/*  Adds characters to pattern until certain depth is reached, checking that
 *	each new addition is box-valid. Once certain depth is reached, the remaining
 *	possible branches are added to the queue for future iterations to handle.
 *
 *	Instead of needing to reallocate space every time an extra
 *	long is needed, allocate space for the maximum length word.
 *	When queueing/dequeueing, the array is 'shortened' as much
 *	as possible.
 *
 * @param depth		how far to travel along the branches before queueing
 * @param word_arr	long array containing pattern
 * @param len		length of pattern (number of 'characters')
*/

/* TODO: Document how the tables work in process (O(1) search) */
/* Generally: i is index in array, j, k are actual lengths
i.e. '0':  i=0, j=k=1 */
void process(int depth, unsigned long long *word_arr, int len) {
	int k = len, size, count, long_count;
	int x, i = 0, j = 0, letter = 0, lc, c,d, counter, it, flip;
	int valid = 1, next, o_valid, have_gone_back = 0; /*HGB: see BOT reverse */
	char *temp_ptr, *string, next_char = '/';
	unsigned long *new_ST;
	unsigned long long exponent = 0, exponent_2,exponent_3, exp_temp, box;
	unsigned long long *arr, word_v = 0, old_word_v = 0, temp = 0;
	new_ST = malloc(sizeof(long) * (alpha_size+1));
	/* Reset tables */
	memset(ST, 0, sizeof(long) * (alpha_size+1));
	memset(LOT, -1, alpha_size * sizeof(int));
	memset(BOT, 0, BOT_size/8);
	string = longs_to_string(word_arr,len);
	if (DEBUG) {
		printf("<%d>New process: len == %d,", my_id, len);
		printf("pattern == %s\n", string);
		printf("Starting setup\n");
	}
	//printf("<%d> is processing %llu %llu %llu %llu\n",my_id, word_arr[0],word_arr[1],word_arr[2],word_arr[3]);
	size = (int) (len / chars_per_long);
	if (len % chars_per_long != 0) {
		size++;
	}
	arr = word_arr;

	if (DEBUG_V2) {
		for (x = 0; x < size; x++) {
			printf("arr[%d]: %llu\n", x, arr[x]);
		}
	}
	/* Reconstruct ST, LOT, BOT */
	for (it = 0; it < size; it++) {
		if (DEBUG_V2) printf("it: %d, size: %d\n", it, size);
		word_v = arr[it];
		exponent = 1;
		if (it < size - 1) { /* Non-last -> full long */
			count = chars_per_long;
			for (x = 0; (x < chars_per_long-1); x++) {
				exponent *= alpha_size; /*pow(alpha_size, chars_per_long-1)*/
			}
		}
		else {	/* last long, could be < chars_per_long */
			count = len - (it * chars_per_long);
			for (x = 0; x < count-1; x++) {
				exponent *= alpha_size; /*pow(alpha_size, count-1)*/
			}
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
			if (letter >= alpha_size) {
				printf("letter: %d\n", letter);
			}
			assert(letter < alpha_size);
			lc = LOT[letter];	/* Compare to LOT */
			if (DEBUG_V2) printf("i: %d, letter: %d; lc = %d\n", i, letter, lc);
			if (lc != -1) {	/* Box occured */
			/* Only adding most closest box as only one that will
			appear when using suffix table */
				if (i + (it * chars_per_long) - lc +1 > alpha_size) { /* internal box exists */
					goto eloop;
				}
				/* Calculate box */
				if (DEBUG_V2) printf("exp_temp: %llu,\told_word_v: %llu\n", exp_temp, old_word_v);
				exponent_2 = exp_temp;
				temp = old_word_v;
				j = lc / chars_per_long;
				if (lc % chars_per_long == 0 && lc != 0) j--;
				if (it == j) {
					if (DEBUG_V2) printf("box in this long\n");
					/* Handle boxes in this long */
					lc %= chars_per_long;
					for (j = 1; j < lc; j++) {	/* remove left of box */
						temp %= exponent_2;
						exponent_2 /= alpha_size;
					}
					if (DEBUG_V2) printf("exponent: %llu\texponent_2: %llu\n", exponent, exponent_2);
					if (DEBUG_V2) printf("temp: %llu\tword_v: %llu\n", temp, word_v);
					assert(exponent > 0);
					box = temp / exponent;	/* remove right of box */
					if (DEBUG_V2) printf("box: %llu\n", box);
					assert(box <= BOT_size);
				}
				else {
					if (DEBUG_V2) {
						printf("spillover method\n");
						printf("word_v: %llu, box: %llu, temp %llu\n",
							word_v,box, temp);
					}
					CALC_SPILLOVER_BOX(lc, temp, exponent, exponent_2, exponent_3, word_arr,i, box, it);
				}
				if (DEBUG_V2) printf("box found. Value: %llu\n", box);
				/* flip bit */
				/* Move (box/8) bytes and change (box%8) bit */
				assert(box <= BOT_size);
				temp_ptr = (char*) (BOT + (int) (box/8));
				flip = 128;
				for (j = 0; j < box % 8; j++) {
					flip /= 2;
				}
				if (((*temp_ptr) &flip) != 0) {
					/* Hotfix mentioned in FAILSAFE_MODE. Should a bad
					item be queued, it will be detected here when the boxes
					are being recalculated. When updating boxes, if box
					had already occured, then invalid and simply get next
					from queue. */
					goto end_p; //HOTFIX
					/*
					string = longs_to_string(arr,len);
					printf("string: %s\n", string);
					printf("valid: %d\n", deep_check(string,len));
					printf("valid_l: %d\n", deep_check_l(arr,len));
					assert(((*temp_ptr) & flip) == 0);
					*/
				}
				*temp_ptr ^= flip;
				if (DEBUG_V2) printf("Updated BOT entry %llu\n", box);
			}
			eloop:
			/* update LOT */
			assert(letter <alpha_size && letter >= 0);
			LOT[letter] = i+1 + chars_per_long * it;

			/* Setup for next iteration */
			assert(exponent > 0);
			word_v %= exponent;
			exponent /= alpha_size;
		}
		word_v = old_word_v;
	}
	it--;	/* pointing to last long in array */
	/* Update suffix table (ST) */
	exponent = 1;
	for (x = 0; x < alpha_size; x++) {
		exponent *= alpha_size; /*pow(alpha_size, alpha_size)*/
	}
	exponent_2 = exponent;
	if (it == 0 || count > alpha_size) { /* first long or no underflow risk */
		if (DEBUG_V2) printf("exp: %llu\tword_v: %llu\n", exponent, word_v);
		for (i = alpha_size; i>=0; i--) {
			ST[i] = word_v;
			if (DEBUG_V2) printf("added suffix %llu at index %d\n", word_v, i);
			assert(exponent > 0);
			word_v %= exponent;
			exponent /= alpha_size;
		}
	} else {	/* Subsequent longs with underflow risk */
		exponent = 1;
		for (x = 0; x < count; x++) {
			exponent *= alpha_size; /*pow(alpha_size, count)*/
		}
		if (DEBUG_V2) printf("exponent: %llu\n", exponent);
		for (i = count-1; i >= 0; i--) { /* Current long */
			word_v %= exponent;
			exponent /= alpha_size;
			ST[i] = word_v;
			if (DEBUG_V2) printf("added suffix %llu at index %d\n", word_v, i);
		}
		word_v = word_arr[it-1];
		exponent = 1;
		for (x = 0; x < count; x++) {
			exponent *= alpha_size; /*pow(alpha_size, count)*/
		}
		exponent_2 = alpha_size;
		for (i = count; i <= alpha_size; i++) { /* Previous long */
			ST[i] = ST[count-1] + (word_v % exponent_2) * exponent;
			exponent_2 *= alpha_size;
			if (DEBUG_V2) printf("ST[%d] now %lu\n", i + count-1, ST[i+count-1]);
		}
	}
	if (DEBUG_V2) {
		for (j = 0; j < size; j++) {
			printf("arr[%d] == %llu\n", j, word_arr[j]);
		}
	}
	if (DEBUG) printf("end of setup\n");

	/* Iterations */
	i = k;
	count = i+1;
	counter = count;

	word_v = old_word_v;
	if (count % chars_per_long == 1 && count != 1) {
		arr[count/chars_per_long] = word_v;
		word_v = 0;
	}

	/* TODO: increase long count if needed */
	word_v *= alpha_size;
	if (DEBUG_V2) printf("word_v == %llu\n", word_v);
	next_char = '/';
	while (i >= k) {
		if (DEBUG) {
			printf("exploring position i: %d, count:%d\n", i, count);
			printf("word_v at start: %llu\n", word_v);
		}
		if (i == k + depth) {;
			if (DEBUG) {
				printf("(%d)string too deep\n", my_id);
			}
			/* TODO: Figure out why invalids are getting put here
			(suspect it's due to a multiplication after a check to
			'setup' next one, which it then queues without checking)*/
			#if FAILSAFE_MODE >= 2
			if (deep_check_l(arr,count-1)) {	/* Failsafe check */
			#else
			if (1) {
			#endif
				if (my_id != ROOT_PROCESS) {
					/* TODO: swap to count-1 */
					MPI_Send(&i, 1, MPI_INT, ROOT_PROCESS, EXPECT_QUEUE_TAG,
						MPI_COMM_WORLD);
					j = (int) ((count-1) / chars_per_long);
					if ((count-1) % chars_per_long != 0) {
						j++;
					}
					assert(j <= max_longs_needed);
					MPI_Send(word_arr, j, MPI_LONG_LONG, ROOT_PROCESS, QUEUE_TAG,
						MPI_COMM_WORLD);
					//printf("<%d>count: %d. Enqueued: %llu %llu %llu %llu\n", my_id,j ,word_l_arr[0],word_l_arr[1],word_l_arr[2], word_l_arr[4]);
				} else {
					if (DEBUG) {
						free(string);
						string = longs_to_string(word_arr, count-1);
						printf("<%d>about to enqueue %s; i: %d, count-1: %d\n",my_id, string, i ,count-1);
						printf("word_arr[0] == %llu\n", word_arr[0]);
					}
					enqueue_l(word_arr, count-1);
					if (DEBUG) printf("queue count: %lld\n", get_queue_size_l());
				}
			}
			i--;
			/* Number of longs needed */
			long_count = count / (chars_per_long);
			if (count % chars_per_long == 0 && count != 0) {
				long_count--; /* if last char, shouldn't increment count */
			}
			/* TODO: possible first error*/
			arr[long_count] = word_v;
			if (DEBUG_V2) printf("i--: long_count == %d, count == %d\n", long_count, count);
			if (count % (chars_per_long) != 1 || count == 1) {
				word_v /= alpha_size;
			} else {
				if (DEBUG_V2) printf("count: %d reverting to arr[%d]: %llu\n", count, long_count-1, arr[long_count-1]);
				word_v = arr[long_count-1];
			}
			count--;
			next_char = alpha[0] + (word_v % alpha_size);
			
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
		next_char++;

		if (word_v % alpha_size == alpha_size-1) {
			if (DEBUG) {
				printf("done with this branch\n");
			}
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
			for (x = 0; x <= alpha_size; x++) {
				exp_temp *= alpha_size; /*pow(alpha_size, alpha_size+1)*/
			}
			word_v /= alpha_size;

			if (DEBUG_V2) printf("word_v / alpha_size. now %llu\n", word_v);
			new_ST[alpha_size] = word_v % exp_temp;

			/* Handling 'underflow' to multiple longs  in suffix table */
			j = (i+1) % chars_per_long;
			if (j == 0) j = chars_per_long;	/*chars_per -> end of long */
			if (i % chars_per_long == chars_per_long-1) j = 0; /* -> reverted */

			if (j < alpha_size+1 && i / chars_per_long > 0) { /* previous long needed */
				if (DEBUG_V2) printf("rollover revert for i:%d, j:%d\n",i, j);
				/* TODO: Check if below accesses correct entry */
				temp = word_arr[(i+1) / chars_per_long - 1];
				j = alpha_size+1 - j;
				exp_temp = 1;
				for (x = 0; x < j; x++) {
					exp_temp *= alpha_size; /*pow(alpha_size, j)*/
				}
				d = ((temp % exp_temp) / (exp_temp/alpha_size));
				exponent_3 = 1;
				for (x = 0; x < alpha_size; x++) {
					exponent_3 *= alpha_size; /*pow(alpha_size, alpha_size)*/
				}
				new_ST[alpha_size] = new_ST[alpha_size-1] + d * exponent_3;

				if (DEBUG_V2) printf("letter: %d\n", d);
				if (DEBUG_V2) printf("temp = %llu, j = %d, exp_temp = %llu, new_st[al] = %ld\n", temp, j , exp_temp, new_ST[alpha_size]);
 			}
			else {
				if (DEBUG_V2) printf("non-rollover revert for i:%d, j:%d\n",i,j);
			}
			if (DEBUG_V2) printf("st[%d]: %lu\n", alpha_size, new_ST[alpha_size]);

			/* Update long_count and update array */
			long_count = count / (chars_per_long);
			if (count % chars_per_long == 0 && count != 0) {
				long_count--; /* if last char, shouldn't increment count */
			}
			/* TODO: possible bug here with long_count setting */
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
			next_char = alpha[0] + (word_v % alpha_size);
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
		// else increment word_v
		//assert(next_char == word[i]); TODO clean up
		//if (word[i] != alpha[0]) {
		if (next_char != alpha[0]) {
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
		next = word_v % alpha_size;

		if (DEBUG_V2) printf("added: %d\n", next);
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
		else { /* Adding a new character to sequence */
			for (j = 1; j <= alpha_size; j++) {
				new_ST[j] = ST[j-1] * alpha_size + next;
				if (DEBUG_V2) printf("new_st: [%d] - %ld\n", j, new_ST[j]);
				counter++;
			}
		}
		if (DEBUG_V2) printf("\n");
		for (j = 0; j <= alpha_size; j++) { /* Copy ST over */
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
			/* checks if box occured & update BOT if valid */
			/* Move (box/8) bytes and change (box%8) bit */
			if (DEBUG_V2) printf("box value: %llu\n", box);
			assert(box <= BOT_size);
			temp_ptr = (char*) (BOT + (int) (box/8));
			flip = 128;
			for (j = 0; j < box % 8; j++) {
				flip /= 2;
			}
			c = (*temp_ptr) & (flip);
			if (DEBUG_V2) printf(">c == %d\n", c);
			if (c == 0) {
				if (DEBUG_V2) printf(">updating box\n");
				o_valid = 1;
				if (((*temp_ptr) &flip) != 0) {
					printf("valid: %d\n", deep_check_l(arr,len));
				}
				/* If failsafe enabled, assert to verify working correctly.
				If not,	and problem occurs, just skip this queued item */
				#if FAILSAFE_MODE > 0
				assert(((*temp_ptr) & flip) == 0);
				#else
				if (((*temp_ptr) & flip) != 0) {
					goto end_p; /* Documented above */
				}
				#endif
				*temp_ptr ^= flip;
			} else {
				o_valid = 0;
			}
		}
		if (o_valid == 1) {
			have_gone_back = 0;
			/* Update LOT */
			if (DEBUG_V2) printf("LOT[%d]: %d -> %d\n", next, LOT[next], i+1);
			LOT[next] = i+1;
		}
		if (o_valid) {
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
				arr[long_count] = word_v;
				word_v = 0;
				//printf("{%d} realloc, number longs: %u\n", my_id, (long_count + 1));
				//word_arr = realloc(word_arr, sizeof(unsigned long long) * (long_count + 2));
				arr[long_count + 1] = 0;
				if (DEBUG_V2) {
					printf("[increasing long count\n");
					printf("arr[%d] = %llu\n",long_count, arr[long_count]);
					printf("arr[%d] = %llu\n", long_count+1, arr[long_count+1]);
					printf("word_v now: %llu]\n", word_v);
				}
			}
			count++;
			if (DEBUG_V2) printf("word_v * alpha_size. now %llu\n", word_v);
			assert(i <= max_word_size[alpha_size]);
			if (i > max_length) {
				max_length = i;
				free(max_word);
				max_word = longs_to_string(max_word_l, max_length);
				x = max_length / chars_per_long;
				if (max_length % chars_per_long > 0) {
					x++;
				}
				for (j = 0; j < x; j++) {
					max_word_l[j] = word_arr[j];
				}
			}
			/* Setup for next run */
			next_char = '/';
		} else if (DEBUG) {
			printf("repeated box\n");
		}
	}
	if (DEBUG) {
		printf("Done process\n");
	}
	end_p:
	free(new_ST);
	free(string);
}

/*  Checks if the last character added to a pattern is box-valid
 *	(i.e. the box ending in the len'th character is never repeated)
 *
 *	@param pattern	pattern to be checked, stored as a string
 *	@param len		length of pattern
 *	@return			VALID or INVALID respectively
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
	#if FAILSAFE_MODE >= 1	
		assert(b1_start <= max_word_size[alpha_size]);
	#endif
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

/*  Performs a deep check of a pattern stored as a string
*	in which all boxes inside the pattern are checked for repeats, not
*  	just the last box added.
*
*	It does this by reconstructing the pattern, character by character,
*	while checking boxes as it goes.
*
*	@param 	str		pattern stored as a string
*	@param	len		length of the pattern
*	@return			VALID or INVALID respectively
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

/* 	Performs a deep check of a pattern stored as an array of longs 
*	in which all boxes inside the pattern are checked for repeats,
*	not just the last box added.
*
*	It does this by reconstructing the pattern, character by character,
*	while checking boxes as it goes.
*
*	@param 	word_arr	pattern stored as a long array
*	@param	len			length of the pattern (number of 'characters')
*	@return				VALID or INVALID respectively
*/
int deep_check_l(unsigned long long *word_arr, int len) {
	unsigned long long word_v, old_word_v, exponent, exp_temp, box, temp, exponent_2, exponent_3;
	int size = (int) (len / chars_per_long); /* # longs */
	int flip, i, j, k, c, x, valid = VALID, count, letter, lc, temp_int;
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
			exponent = 1;
			for (x = 0; x < (chars_per_long-1); x++) {
				exponent *= alpha_size; /*alpha_size ^chars_per_long-1*/
			}
		}
		else {	/* last long, could be < chars_per_long */
			count = len - (i * chars_per_long);
			exponent = 1;
			for (x = 0; x < (count-1); x++) {
				exponent *= alpha_size; /*alpha_size ^count-1*/
			}
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
			#if FAILSAFE_MODE >= 1
			if (letter >= alpha_size) {
				/* If assert would fail, print nearby situation.
				(provides some debug information in a crash) */
				printf("Assertion failed: letter >= alpha_size\n");
				printf("\nletter: %d, i: %d, j: %d, size: %d\n",
					letter,i,j, size);
				printf("count: %d, word_v: %llu, exponent: %llu\n",
					count, word_v, exponent);
				char *string = longs_to_string(word_arr, len);
				printf("word: %s\n", string);
				free(string);
				assert(letter < alpha_size);
			}
			#endif
			/* Compare to LOT */
			lc = LOT2[letter];
			if (lc != 0) { /* box occured */
				if (j+1 + (i*chars_per_long) - lc > alpha_size) {
					goto eloop; /* internal box exists. Skip */
				}
				/* Calculate box */
				exponent_2 = exp_temp;
				temp = old_word_v;
				temp_int = lc / chars_per_long;
				if (lc % chars_per_long == 0 && lc != 0) temp_int--;
				if (temp_int == i) {
					/* Handle boxes in this long */
					temp_int = lc % chars_per_long;
					for (k = 1; k < temp_int; k++) { /* remove left of box */
						temp %= exponent_2;
						exponent_2 /= alpha_size;
					}
					box = temp / exponent;	/* remove right of box */
					#if FAILSAFE_MODE >= 1
					if (box > BOT_size) {
						/* If assert would fail, print nearby situation.
					   (provides some debug information in a crash) */
						printf("<Assertion failed\nbox == %llu\n", box);
						printf("<len: %d, i: %d, j: %d, word_v == %llu\n",len,i,j,word_arr[i]);
						printf("<lc: %d, k: %d, count: %d, letter: %d\n",lc,k,count, letter);
						printf("<string: %s\n", longs_to_string(word_arr, len));
						assert(box <= BOT_size);
					}
					#endif	
				}
				else {
					/* handle boxes that spills over to previous long */
					CALC_SPILLOVER_BOX(lc, temp, exponent, exponent_2, exponent_3, word_arr,j, box,i);
					assert(box <= BOT_size);
				}
				/* Check if BOT entry != 0. If so, return 0. else update BOT */
					/* Move (box/8) bytes and change (box%8) bit */
				if (DEBUG_V2) printf("box value: %llu\n", box);
				temp_ptr = (char*) (BOT2 + (int) (box/8));
				flip = 128;
				assert(box <= BOT_size);
				for (k = 0; k < box % 8; k++) {
					flip /= 2;
				}
				c = (*temp_ptr) & (flip);
				if (DEBUG_V2) printf(">c == %d\n", c);
				if (c == 0) {
					if (DEBUG_V2) printf(">updating box\n");
					/* If failsafe enabled, assert to verify working correctly.
					If not,	and problem occurs, just report INVALID */
					#if FAILSAFE_MODE > 0
					assert(((*temp_ptr) & flip) == 0);
					#else
					if (((*temp_ptr) & flip) != 0) {
						valid = INVALID;
						goto end;
					}
					#endif
					*temp_ptr ^= flip;
				} else {
					valid = INVALID;
					//printf("i: %d, j: %d, box: %llu\n",i, j, box);
					goto end;
				}
			}
			eloop:
			/* Update LOT */
			LOT2[letter] = j+1 + (i * chars_per_long);
			/* Setup for next iteration*/	
			assert(exponent > 0);
			word_v %= exponent;
			exponent /= alpha_size;
		}
	}
	end:
	free(BOT2);
	free(LOT2);

#if FAILSAFE_MODE >= 3
	/* Temporary assert to confirm all is working */
	temp_ptr = longs_to_string(word_arr, len);
	assert(valid == deep_check(temp_ptr, len));
	free(temp_ptr);
#endif
	return valid;
}

/*	Converts a pattern stored as an array of longs into a string
 *
 *	@param	 arr	array of longs storing the pattern
 * 	@param	 len	length of the pattern (number of 'characters')
 * 	@return			String containing the pattern
 */
char *longs_to_string(unsigned long long *arr, int len) {
	char letter, *string = malloc(sizeof(char) * (len+1));
	int i, j, count;
	unsigned long long word_v;
	int size = (int) (len / chars_per_long); /* # longs */
	if (len % chars_per_long != 0) {
		size++;
	}
	for (i = 0; i < size; i++) { /* For each long */
		word_v = arr[i];
		if (i < size - 1) { /* Non-last -> full long */
			count = chars_per_long;
		}
		else { /* last long, could be < chars_per_long */
			count = len - (i * chars_per_long);
		}
		for (j = count-1; j >= 0; j--) { /* Get chars from long */
			letter = alpha[0] + (word_v % alpha_size);
			if (!(letter >= alpha[0] && letter <= alpha[alpha_size])) {
				printf("<<letter: %d, word_v: %llu, len: %d\n", letter, word_v, len);
				printf("i: %d, j: %d>>\n",i,j);
			}
			assert(letter >= alpha[0] && letter <= alpha[alpha_size]);
			string[i * chars_per_long + j] = letter;
			word_v /= alpha_size;
		}
	}
	string[len] = '\0';
	return string;
}

/*	Converts a pattern stored as a string into an array of longs 
 *
 *	@param 	string	String containing pattern
 *	@param	len		length of the string
 *	@return			array of longs storing the pattern
*/
unsigned long long *string_to_longs(char *string, int len) {
	int i, j, count, size = (int) (len / chars_per_long), x;
	unsigned long long word_v, exponent, *arr;

	if (len % chars_per_long != 0) {
		size++;
	}
	/* Make as big as possibly needed */
	arr = malloc(sizeof(unsigned long long) * max_longs_needed);
	for (i = 0; i < max_longs_needed; i++) {
		arr[i] = 0;
	}

	for (i = 0; i < size; i++) { /* For each long */
		word_v = 0;
		exponent = 1;
		if (i < size - 1) { /* Non-last means full long */
			count = chars_per_long;
			for (x = 0; x < (chars_per_long-1); x++) {
				exponent *= alpha_size; /*pow(alpha_size, chars_per_long -1)*/
			}
			/* TODO: Check */
		}
		else {	/* last long, could be less than chars_per_long  characters */
			count = len - (i * chars_per_long);
			for (x = 0; x < count-1; x++) {
				exponent *= alpha_size; /*pow(alpha_size, count -1)*/
			}
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

/* 	Saves the queue to a given file.
*	Dequeues each entry and writes its size and length to the file.
*	If enabled, will also write the current maximum entry to the file.
*
*	Should the file name be the same as the global file_name (the default),
*	it will also rename the file_name as needed to ensure that the maximum
*	number of copies (see configuration) is reached.
*	Standard format is: "nX-Y.save", with X being alphabet size, Y being
* 	the current entry, 0 <= Y < COPIES_KEPT.
*
*	Note: It will overwrite files of the same name. So make backups
*	as needed.
*
*	Note: Currently it dequeues the entire queue, so if you wish to continue,
*	this must be followed by a restore_file call.
*
*	File layout: alphabet size, followed by entries stored as:
*	max_len	max[0]	max[1]	max[2]	max[3]	...
*	len(l1)	l1[0]	l1[1]	l1[2]	l1[3] etc
*
*	@param 	save_name	File name to be saved to
*/

void store_queue_to_file(char *save_name) {
	FILE *file;
	unsigned long long queue_size, *arr;
	int size, i;
	if (my_id != ROOT_PROCESS) {
		fprintf(stderr, "Only root process can process queue\n");
		return;
	}
	else if (save_name == NULL) {
		fprintf(stderr, "NULL file name\n");
		return;
	}
	if (strcmp(save_name, file_name) == 0) {
		/* Custom renaming of saves */
		if (current_save_number >= COPIES_KEPT) {
			current_save_number = 0;
		}
		save_name[swap_index] = '0' + current_save_number;
		current_save_number++;
	}
	file = fopen(save_name,"w");
	if (!file) {
		fprintf(stderr, "Unable to create file. ");
		fprintf(stderr,
			"Please create directory 'storage' in root of this repo\n");
		return;
	}
	arr = malloc(sizeof(unsigned long long) * max_longs_needed);
	queue_size = get_queue_size_l();
	fprintf(file, "%d\n", alpha_size);
	if (PRINT_SAVE_LOAD) {
		printf("Storing queue: file: '%s'. Queue at start: %llu\n", save_name, queue_size);
	}
	/* Store current max */
	if (STORE_MAX) {
		fprintf(file, "%d", max_length);
		for (i = 0; i < max_longs_needed; i++) {
			fprintf(file, " %llu", max_word_l[i]);
		}
		fprintf(file, "\n");
	}
	/* Store the queue */
	while (queue_size) {
		for (i = 0; i < max_longs_needed; i++) {
			arr[i] = 0;
		}
		if (!dequeue_l(arr, &size)) {
			fprintf(stderr, "error dequeueing\n");
			goto close_s;
			return;
		}
		fprintf(file,"%d", size);
		for (i = 0; i < max_longs_needed; i++) {
			fprintf(file, " %llu", arr[i]);
		}
		fprintf(file, "\n");
		queue_size = get_queue_size_l();
	}
	if (PRINT_SAVE_LOAD) {
		printf("Storing queue: Queue at end: %llu\n", queue_size);
	}
	close_s:
	free(arr);
	fclose(file);
	return;
}

/* 	Restores the queue from a given file. Reads each entry line by line,
*	saving it's size and long values before enqueueing it.
*	Will only restore files that have the same corresponding alphabet size.
*
*	@param 	save_name	File name be restored from
*/
/* TODO: Optimisation: Reverse direction in file. As it currently stands 
	it restores it in descending order, causing extra transversals in queue */
void restore_queue_from_file(char *save_name) {
	FILE *file;
	int i, rec;
	unsigned long long value, *arr;
	if (my_id != ROOT_PROCESS) {
		fprintf(stderr, "Only root process can process queue\n");
		return;
	}
	else if (save_name == NULL) {
		fprintf(stderr, "NULL file name\n");
		return;
	}

	file = fopen(save_name,"r");

	if (file == NULL) {
		fprintf(stderr, "Unable to open file %s\n", save_name);
		return;
	}

	arr = malloc(sizeof(unsigned long long) * max_longs_needed);
	fscanf(file, "%d", &rec);
	/* Only read from a file with same alpha_size */
	if (rec != alpha_size) {
		fprintf(stderr, "invalid file to read from.\n");
		fprintf(stderr,
			"alpha_size: %d, file_alpha_size: %d\n", alpha_size, rec);
		goto close_r;
	}
	/* Get entries: size followed by longs of array */
	if (PRINT_SAVE_LOAD) {
		printf("Restoring queue: file: '%s'\n", save_name);
	}
	while (fscanf(file, "%d", &rec) > 0) {
		/* get longs */
		for (i = 0; i < max_longs_needed; i++) {
			if (fscanf(file, "%llu", &value) <= 0) {
				/* If a long is missing */
				fprintf(stderr, "unable to retrieve queue entry. ");
				fprintf(stderr, "Insufficient long entries. i: %d\n", i);
				goto close_r;
			}
			arr[i] =  value;
		}
		enqueue_l(arr, rec);
		for (i = 0; i < max_longs_needed; i++) {
			arr[i] = 0;
		}
	}

	if (PRINT_SAVE_LOAD) {
		printf("Restoring queue: Queue at end: %llu\n", get_queue_size_l());
	}
	close_r:
	free(arr);
	fclose(file);
	return;

}