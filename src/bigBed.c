#include "ucsc/common.h"
#include "ucsc/linefile.h"
#include "ucsc/localmem.h"
#include "ucsc/bbiFile.h"
#include "ucsc/bigBed.h"

#include "bigBed.h"
#include "bbiHelper.h"
#include "handlers.h"

/* --- .Call ENTRY POINT --- */
SEXP BBDFile_seqlengths(SEXP r_filename) {
  pushRHandlers();
  SEXP seqlengths;
  struct bbiFile * file = bigBedFileOpen((char *)CHAR(asChar(r_filename)));
  PROTECT(seqlengths = bbiSeqLengths(file));
  bigBedFileClose(&file);
  popRHandlers();
  UNPROTECT(1);
  return seqlengths;
}
