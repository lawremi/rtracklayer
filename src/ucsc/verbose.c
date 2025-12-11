/* verbose.c - write out status messages according to the
 * current verbosity level.  These messages go to stderr. */

/* Copyright (C) 2011 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */

#include "common.h"
#include "portable.h"
#include "verbose.h"

#define logVerbosity 0	/* The level of log verbosity.  0 is silent. */

void verboseVa(int verbosity, char *format, va_list args)
/* Log with at given verbosity vprintf formatted args. */
{
if (verbosity <= logVerbosity)
    {
#if logVerbosity > 0 // to avoid R CMD check warnings
	vfprintf(stderr, format, args);
        fflush(stderr);
#endif
    }
}

void verbose(int verbosity, char *format, ...)
/* Write printf formatted message to log (which by
 * default is stderr) if global verbose variable
 * is set to verbosity or higher. */
{
va_list args;
va_start(args, format);
verboseVa(verbosity, format, args);
va_end(args);
}

static long lastTime = -1;  // previous call time.

void verboseTimeInit(void)
/* Initialize or reinitialize the previous time for use by verboseTime. */
{
lastTime = clock1000();
}

void verboseTime(int verbosity, char *label, ...)
/* Print label and how long it's been since last call.  Start time can be
 * initialized with verboseTimeInit, otherwise the elapsed time will be
 * zero. */
{
assert(label != NULL);  // original version allowed this, but breaks some GCCs
if (lastTime < 0)
    verboseTimeInit();
long time = clock1000();
va_list args;
va_start(args, label);
verboseVa(verbosity, label, args);
verbose(verbosity, ": %ld millis\n", time - lastTime);
lastTime = time;
va_end(args);
}


boolean verboseDotsEnabled()
/* check if outputting of happy dots are enabled.  They will be enabled if the
 * verbosity is > 0, stderr is a tty and we don't appear to be running an
 * emacs shell. */
{
    return FALSE;
}

int verboseLevel(void)
/* Get verbosity level. */
{
return logVerbosity;
}
