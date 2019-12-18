Finds the longest string that has no repeated boxes. Utilises MPI to run
multiple workers at once. Relies on a queue-based system to handle jobs
between workers.

**  For exact setup,compile and execution instructions  **
**              see the 'search.c' file                 **

There are various configurations available (search depth, defaults,
displaying of debugging information / status update, etc).
All of these options are explained in detail in the 'search.c' file.

Further things to note with the configurations (discussed in 'search.c'):
SEARCH DEPTH CONFIGURATIONS
1) Standard search depth: Too high results in long delays between calls,
and slow saving procedure. Too low and setup costs start
making a relative impact.
Recommended: 6 or 7

2) Initial search depth:
This initial run creates branches for children to run.
Ensure that enough items are queued for all processes being used
(e.g. 24 cores will need ideally 3-4 to start on N>=3).
Any higher and the initial setup will take longer than necessary.
Side note: Ensure there is enough remaining if using N=2.
Setting it to 7+ when using N=2 will break it.

Recommended: 3 to 7

SAVE FILE
This program creates save files at specified intervals to allow easier
continuation if the program needs to stop. 
Set STORE_QUEUE to 1 to enable,
QUEUE_TIME to the number of minutes between saves, and
COPIES_KEPT to specify number of saves to be kept.
They will be stored in a directory 'saves', numbered as:
"nX-Y.save", with X the alphabet size, Y the copy number.

OUTPUT
Many outputs can be controlled via definitions in the config file.
This includes:
Debug output    - for future additions if needed
Status updates  - general details of current run
Maximums        - Whether to print when a new max is found
Save/load       - Whether to print when a save/load occurs

FAILSAFE_MODE
Since this is a computationally-heavy code, various levels of
checks are available. This allows a higher-level (slower processing)
mode to be used when additions are added to the code - to ensure everything
is working properly - and a lower level (faster) mode to be used
once everything has been verified, for any possible speed improvements