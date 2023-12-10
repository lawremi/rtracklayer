/*****************************************************************************
 * Copyright (C) 2000 Jim Kent.  This source code may be freely used         *
 * for personal, academic, and non-profit purposes.  Commercial use          *
 * permitted only by explicit agreement with Jim Kent (jim_kent@pacbell.net) *
 *****************************************************************************/
/* hmmstats.h - Stuff for doing statistical analysis in general and 
 * hidden Markov models in particular. */
#ifndef HMMSTATS_H
#define HMMSTATS_H

#define logScaleFactor 1000
/* Amount we scale logs by. */

double calcVarianceFromSums(double sum, double sumSquares, bits64 n);
/* Calculate variance. */

double calcStdFromSums(double sum, double sumSquares, bits64 n);
/* Calculate standard deviation. */

#endif /* HMMSTATS_H */

