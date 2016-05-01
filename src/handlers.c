#include "ucsc/common.h"
#include "ucsc/errAbort.h"

#include "handlers.h"

#define WARN_BUF_SIZE 512
static void R_warnHandler(char *format, va_list args) {
  char warn_buf[WARN_BUF_SIZE];
  vsnprintf(warn_buf, WARN_BUF_SIZE, format, args);
  warning(warn_buf);
}

static void R_abortHandler() {
  error("UCSC library operation failed");
}

extern int R_ignore_SIGPIPE;

void pushRHandlers() {
  pushAbortHandler(R_abortHandler);
  pushWarnHandler(R_warnHandler);
#ifndef WIN32
  R_ignore_SIGPIPE = 1;
#endif
}

void popRHandlers() {
  popAbortHandler();
  popWarnHandler();
#ifndef WIN32
  R_ignore_SIGPIPE = 0;
#endif
}
