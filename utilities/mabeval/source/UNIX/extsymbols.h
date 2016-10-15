/* 04-May-92 written by Ron Shepard */

/* Define the external symbols.  These must be defined so that */
/* the Fortran interface is portable.  This machine-dependency */
/* is localized in this file in order to simplify the C code.  */


/* #elif is not recognized by all C preprocessors, so punt. -rls */

#if defined(_CRAY) || defined(ardent)

  /* Fortran symbols are converted to upper case with no trailing "_". */
#define EXTERNAL_FALLOC FALLOC
#define EXTERNAL_FDATE  FDATE
#define EXTERNAL_HOSTNM HOSTNM
#define EXTERNAL_FWTIME FWTIME
#define EXTERNAL_FSIZE  FSIZE

#else

#if !defined(EXTNAME)
  /* AIX requires underscores if compiled with -qextname */ 
  /* Fortran symbols are converted to lower case with no trailing "_". */
#define EXTERNAL_FALLOC falloc
#define EXTERNAL_FDATE  fdate
#define EXTERNAL_HOSTNM hostnm
#define EXTERNAL_FWTIME fwtime
#define EXTERNAL_FLUSHSTDOUT flushstdout
#define EXTERNAL_FSIZE  fsize

#else
#if defined(EXTNAME)
  /* AIX requires underscores if compiled with -qextname */ 
  /* Fortran symbols are converted to lower case with no trailing "_". */
#define EXTERNAL_FALLOC falloc_
#define EXTERNAL_FDATE  fdate_
#define EXTERNAL_HOSTNM hostnm_
#define EXTERNAL_FWTIME fwtime_
#define EXTERNAL_FLUSHSTDOUT flushstdout_
#define EXTERNAL_FSIZE  fsize_ 

#else

  /* default case: */
  /* Fortran symbols are converted to lower case, and a trailing "_" */
  /* is added automatically by the compiler.                         */
#define EXTERNAL_FALLOC falloc_
#define EXTERNAL_FDATE  fdate_
#define EXTERNAL_HOSTNM hostnm_
#define EXTERNAL_FWTIME fwtime_
#define EXTERNAL_FSIZE  fsize_ 

#endif /*   */
#endif /*  */
#endif /* _CRAY */
