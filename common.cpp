#ifndef _COMMON_H_
#define _COMMON_H_

#include <assert.h>
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

FILE * open_file(const char *name, char * mode){
     FILE * ret=NULL;
     ret=fopen(name, mode);
     if (ret==NULL)  {
	  fprintf(stderr, "ERROR: can not open file %s  at %s:%d\n",
			     name, __FILE__, __LINE__);
	  exit(0);
     }
     return ret;
}


#endif /* _COMMON_H_ */

