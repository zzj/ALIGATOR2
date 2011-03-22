/**
 * @file   common.hpp
 * @author Zhaojun Zhang <zzj@cs.unc.edu>
 * @date   Mon Apr 27 21:55:38 2009
 * 
 * @brief  define the functions which are commonly used.
 * 
 * 
 */
#ifndef _COMMON_H_
#define _COMMON_H_




#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define MAX_SIZE_PER_LINE 40000	/**< max size per line */

FILE * open_file(const char *name, char * mode);

#endif /* _COMMON_H_ */
