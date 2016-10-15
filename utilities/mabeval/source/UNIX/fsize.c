#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "extsymbols.h"
#include <errno.h>

/* return the file size of a file name */
/* fortran automatically passes strings via string, strlen*/
#if defined(_CRAY)

/* fortran character strings are passed by descriptor */
#include <fortran.h>
 int EXTERNAL_FSIZE( fdesc )
_fcd fdesc ; /* fortran character descriptor */


#else

/* fortran string and length are expanded automatically */
#if defined(BIT64) 
long EXTERNAL_FSIZE( name,len )
#else
int EXTERNAL_FSIZE( name,len )
#endif
  char *name ; /* string address */
  int  len  ; /* string length  */

#endif /* _CRAY */


{ 
struct stat stbuf;
  char fnam[32];
#if defined(_CRAY)
  /* _fcdtocp() and _fcdlen() are macros that extract the string address
     bits and string length respectively from fdesc */
char *name ; /* string address */
int  len  ; /* string length  */
  name    = _fcdtocp( fdesc ) ;
  len = _fcdlen( fdesc ) ;
#endif

  if ( len > 31 ) { return -999;}
  strncpy(fnam,name,len);
  fnam[len]='\0';
  if (stat(fnam, &stbuf) == -1 ) { 
   printf("Errorno %d len: %d name: %s\n**\n",errno,len,fnam);
  ;return errno;}
#if defined(BIT64)
  return (long) stbuf.st_size;
#else
  return  stbuf.st_size;
#endif
}
