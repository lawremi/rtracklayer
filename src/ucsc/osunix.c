/* Some wrappers around operating-system specific stuff. 
 *
 * This file is copyright 2002 Jim Kent, but license is hereby
 * granted for all use - public, private or commercial. */

#include "common.h"
#include <dirent.h>
#include <sys/utsname.h>
#include <sys/time.h>
#include <sys/statvfs.h>
#include <pwd.h>
#include <termios.h>
#include "portable.h"
#include "portimpl.h"
#include "_portimpl.h"
#include <sys/wait.h>
#include <regex.h>
#include <utime.h>




off_t fileSize(char *pathname)
/* get file size for pathname. return -1 if not found */
{
struct stat mystat;
ZeroVar(&mystat);
if (stat(pathname,&mystat)==-1)
    {
    return -1;
    }
return mystat.st_size;
}

long clock1000()
/* A millisecond clock. */
{
struct timeval tv;
static long origSec;
gettimeofday(&tv, NULL);
if (origSec == 0)
    origSec = tv.tv_sec;
return (tv.tv_sec-origSec)*1000 + tv.tv_usec / 1000;
}

void sleep1000(int milli)
/* Sleep for given number of 1000ths of second */
{
if (milli > 0)
    {
    struct timeval tv;
    tv.tv_sec = milli/1000;
    tv.tv_usec = (milli%1000)*1000;
    select(0, NULL, NULL, NULL, &tv);
    }
}

long clock1()
/* A seconds clock. */
{
struct timeval tv;
gettimeofday(&tv, NULL);
return tv.tv_sec;
}

char *getCurrentDir()
/* Return current directory.  Abort if it fails. */
{
static char dir[PATH_LEN];

if (getcwd( dir, sizeof(dir) ) == NULL )
    errnoAbort("getCurrentDir: can't get current directory");
return dir;
}

void setCurrentDir(char *newDir)
/* Set current directory.  Abort if it fails. */
{
if (chdir(newDir) != 0)
    errnoAbort("setCurrentDir: can't to set current directory: %s", newDir);
}

boolean makeDir(char *dirName)
/* Make dir.  Returns TRUE on success.  Returns FALSE
 * if failed because directory exists.  Prints error
 * message and aborts on other error. */
{
int err;
if ((err = mkdir(dirName, 0777)) < 0)
    {
    if (errno != EEXIST)
	{
	perror("");
	errAbort("Couldn't make directory %s", dirName);
	}
    return FALSE;
    }
return TRUE;
}


struct fileInfo *listDirXExt(char *dir, char *pattern, boolean fullPath, boolean ignoreStatFailures)
/* Return list of files matching wildcard pattern with
 * extra info. If full path is true then the path will be
 * included in the name of each file. */
{
struct fileInfo *list = NULL, *el;
struct dirent *de;
DIR *d;
int dirNameSize = strlen(dir);
int fileNameOffset = dirNameSize+1;
char pathName[512];

if ((d = opendir(dir)) == NULL)
    return NULL;
memcpy(pathName, dir, dirNameSize);
pathName[dirNameSize] = '/';

while ((de = readdir(d)) != NULL)
    {
    char *fileName = de->d_name;
    if (differentString(fileName, ".") && differentString(fileName, ".."))
	{
	if (pattern == NULL || wildMatch(pattern, fileName))
	    {
	    struct stat st;
	    bool isDir = FALSE;
	    int statErrno = 0;
	    strcpy(pathName+fileNameOffset, fileName);
	    if (stat(pathName, &st) < 0)
		{
		if (ignoreStatFailures)
		    statErrno = errno;
		else
    		    errAbort("stat failed in listDirX");
		}
	    if (S_ISDIR(st.st_mode))
		isDir = TRUE;
	    if (fullPath)
		fileName = pathName;
	    el = newFileInfo(fileName, st.st_size, isDir, statErrno, st.st_atime);
	    slAddHead(&list, el);
	    }
	}
    }
closedir(d);
slSort(&list, cmpFileInfo);
return list;
}

struct fileInfo *listDirX(char *dir, char *pattern, boolean fullPath)
/* Return list of files matching wildcard pattern with
 * extra info. If full path is true then the path will be
 * included in the name of each file. */
{
return listDirXExt(dir, pattern, fullPath, FALSE);
}

time_t fileModTime(char *pathName)
/* Return file last modification time.  The units of
 * these may vary from OS to OS, but you can depend on
 * later files having a larger time. */
{
struct stat st;
if (stat(pathName, &st) < 0)
    errAbort("stat failed in fileModTime: %s", pathName);
return st.st_mtime;
}


char *getHost()
/* Return host name. */
{
static char *hostName = NULL;
static char buf[128];
if (hostName == NULL)
    {
    hostName = getenv("HTTP_HOST");
    if (hostName == NULL)
        {
	hostName = getenv("HOST");
	if (hostName == NULL)
	    {
	    if (hostName == NULL)
		{
		static struct utsname unamebuf;
		if (uname(&unamebuf) >= 0)
		    hostName = unamebuf.nodename;
		else
		    hostName = "unknown";
		}
	    }
        }
    strncpy(buf, hostName, sizeof(buf));
    chopSuffix(buf);
    hostName = buf;
    }
return hostName;
}

char *semiUniqName(char *base)
/* Figure out a name likely to be unique.
 * Name will have no periods.  Returns a static
 * buffer, so best to clone result unless using
 * immediately. */
{
int pid = getpid();
int num = time(NULL)&0xFFFFF;
char host[512];
strcpy(host, getHost());
char *s = strchr(host, '.');
if (s != NULL)
     *s = 0;
subChar(host, '-', '_');
subChar(host, ':', '_');
static char name[PATH_LEN];
safef(name, sizeof(name), "%s_%s_%x_%x",
	base, host, pid, num);
return name;
}

static void eatSlashSlashInPath(char *path)
/* Convert multiple // to single // */
{
char *s, *d;
s = d = path;
char c, lastC = 0;
while ((c = *s++) != 0)
    {
    if (c == '/' && lastC == c)
        continue;
    *d++ = c;
    lastC = c;
    }
*d = 0;
}

static void eatExcessDotDotInPath(char *path)
/* If there's a /.. in path take it out.  Turns 
 *      'this/long/../dir/file' to 'this/dir/file
 * and
 *      'this/../file' to 'file'  
 *
 * and
 *      'this/long/..' to 'this'
 * and
 *      'this/..' to  ''   
 * and
 *       /this/..' to '/' */
{
/* Take out each /../ individually */
for (;;)
    {
    /* Find first bit that needs to be taken out. */
    char *excess= strstr(path, "/../");
    char *excessEnd = excess+4;
    if (excess == NULL || excess == path)
        break;

    /* Look for a '/' before this */
    char *excessStart = matchingCharBeforeInLimits(path, excess, '/');
    if (excessStart == NULL) /* Preceding '/' not found */
         excessStart = path;
    else 
         excessStart += 1;
    strcpy(excessStart, excessEnd);
    }

/* Take out final /.. if any */
if (endsWith(path, "/.."))
    {
    if (!sameString(path, "/.."))  /* We don't want to turn this to blank. */
	{
	int len = strlen(path);
	char *excessStart = matchingCharBeforeInLimits(path, path+len-3, '/');
	if (excessStart == NULL) /* Preceding '/' not found */
	     excessStart = path;
	else 
	     excessStart += 1;
	*excessStart = 0;
	}
    }
}

char *simplifyPathToDir(char *path)
/* Return path with ~ and .. taken out.  Also any // or trailing /.   
 * freeMem result when done. */
{
/* Expand ~ if any with result in newPath */
char newPath[PATH_LEN];
int newLen = 0;
char *s = path;
if (*s == '~')
    {
    char *homeDir = getenv("HOME");
    if (homeDir == NULL)
        errAbort("No HOME environment var defined after ~ in simplifyPathToDir");
    ++s;
    if (*s == '/')  /*    ~/something      */
        {
	++s;
	safef(newPath, sizeof(newPath), "%s/", homeDir);
	}
    else            /*   ~something        */
	{
	safef(newPath, sizeof(newPath), "%s/../", homeDir);
	}
    newLen = strlen(newPath);
    }
int remainingLen  = strlen(s);
if (newLen + remainingLen >= sizeof(newPath))
    errAbort("path too big in simplifyPathToDir");
strcpy(newPath+newLen, s);

/* Remove //, .. and trailing / */
eatSlashSlashInPath(newPath);
eatExcessDotDotInPath(newPath);
int lastPos = strlen(newPath)-1;
if (lastPos > 0 && newPath[lastPos] == '/')
    newPath[lastPos] = 0;

return cloneString(newPath);
}

void childExecFailedExit(char *msg)
/* Child exec failed, so quit without atexit cleanup */
{
fprintf(stderr, "child exec failed: %s\n", msg);
fflush(stderr);
_exit(1);  // Let the parent know that the child failed by returning 1.

/* Explanation:
_exit() is not the normal exit().  
_exit() avoids the usual atexit() cleanup.
The MySQL library that we link to uses atexit() cleanup to close any open MySql connections.
However, because the child's mysql connections are shared by the parent,
this causes the parent MySQL connections to become invalid,
and causes the puzzling "MySQL has gone away" error in the parent
when it tries to use its now invalid MySQL connections.
*/

}

static void execPStack(pid_t ppid)
/* exec pstack on the specified pid */
{
char *cmd[3], pidStr[32];
safef(pidStr, sizeof(pidStr), "%ld", (long)ppid);
cmd[0] = "pstack";
cmd[1] = pidStr;
cmd[2] = NULL;

// redirect stdout to stderr
if (dup2(2, 1) < 0)
    errAbort("dup2 failed");

execvp(cmd[0], cmd);

childExecFailedExit(cmd[0]); // cannot use the normal errAbort.

}

void vaDumpStack(char *format, va_list args)
/* debugging function to run the pstack program on the current process. In
 * prints a message, following by a new line, and then the stack track.  Just
 * prints errors to stderr rather than aborts. For debugging purposes
 * only.  */
{
static boolean inDumpStack = FALSE;  // don't allow re-entry if called from error handler
if (inDumpStack)
    return;
inDumpStack = TRUE;

fflush(stdout);  // clear buffer before forking
vfprintf(stderr, format, args);
fputc('\n', stderr);
fflush(stderr);
pid_t ppid = getpid();
pid_t pid = fork();
if (pid < 0)
    {
    perror("can't fork pstack");
    return;
    }
if (pid == 0)
    execPStack(ppid);
int wstat;
if (waitpid(pid, &wstat, 0) < 0)
    perror("waitpid on pstack failed");
else
    {
    if (WIFEXITED(wstat))
        {
        if (WEXITSTATUS(wstat) != 0)
            fprintf(stderr, "pstack failed\n");
        }
    else if (WIFSIGNALED(wstat))
        fprintf(stderr, "pstack signaled %d\n", WTERMSIG(wstat));
    }
inDumpStack = FALSE;
}

void dumpStack(char *format, ...)
/* debugging function to run the pstack program on the current process. In
 * prints a message, following by a new line, and then the stack track.  Just
 * prints errors to stderr rather than aborts. For debugging purposes
 * only.  */
{
va_list args;
va_start(args, format);
vaDumpStack(format, args);
va_end(args);
}

boolean maybeTouchFile(char *fileName)
/* If file exists, set its access and mod times to now.  If it doesn't exist, create it.
 * Return FALSE if we have a problem doing so (e.g. when qateam is gdb'ing and code tries 
 * to touch some file owned by www). */
{
if (fileExists(fileName))
    {
    struct utimbuf ut;
    ut.actime = ut.modtime = clock1();
    int ret = utime(fileName, &ut);
    if (ret != 0)
	{
	warn("utime(%s) failed (ownership?)", fileName);
	return FALSE;
	}
    }
else
    {
    FILE *f = fopen(fileName, "w");
    if (f == NULL)
	return FALSE;
    else
	carefulClose(&f);
    }
return TRUE;
}

boolean isRegularFile(char *fileName)
/* Return TRUE if fileName is a regular file. */
{
struct stat st;

if (stat(fileName, &st) < 0)
    return FALSE;
if (S_ISREG(st.st_mode))
    return TRUE;
return FALSE;
}
