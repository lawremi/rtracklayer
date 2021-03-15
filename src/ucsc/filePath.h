/* filePath - stuff to handle file name parsing. */

/* Copyright (C) 2011 The Regents of the University of California 
 * See README in this or parent directory for licensing information. */
#ifndef FILEPATH_H
#define FILEPATH_H

#include "common.h"

char *expandRelativePath(char *baseDir, char *relPath);
/* Expand relative path to more absolute one. */

void undosPath(char *path);
/* Convert '\' to '/' in path. (DOS/Windows is typically ok with
 * this actually.) */

#endif /* FILEPATH_H */
