/*------------------------------------------------------------------------------
          1         2         3         4         5         6         7
01234567890123456789012345678901234567890123456789012345678901234567890123456789

'util.h'	created  12/09/94,
                   B. Verdonck, T. Geraud, J.F. Mangin

			modifie le 28/05/1997 par cointepa@ima.enst.fr

			Copyright (c) 1994
			Dept. IMA  TELECOM Paris
			46, rue Barrault  F-75634 PARIS Cedex 13
			All Rights Reserved

------------------------------------------------------------------------------*/
/** @file util.h
 * @brief General utility functions used allover the tivoli library.
 */

/*
|
|  This file groups general utility functions used allover the tivoli
|  library.  Please contact B.V. before adding additional functions.
|
*/

#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#ifndef TIVOLI_UTIL_H
#define TIVOLI_UTIL_H

/* Some handy shortcut macro's */
/* T_ refers to the Tivoli trace functionalities */
#define T_IN(s)  printftrace(IN,s)
#define T_OUT    printftrace(OUT, "")

#define T_RETURN_PB   printftrace(OUT,"");return(PB)
#define T_RETURN_NULL printftrace(OUT,"");return(NULL)
#define T_RETURN_OK   printftrace(OUT,"");return(OK)

/*JFM*/
/*define bete et mechants*/

#define SHOW_LIM 1
#define HIDE_LIM 0

/* boolean values */

#define TRUE	1
#define FALSE	0

#define YES     1
#define NO      0

#define ONE     1
#define ZERO    0

/* value returned by testing functions */

#define OK	1
#define PB	0

/* first argument of tracing function 'printftrace' */

#define IN  1
#define MSG 2
#define OUT 3

/* flag modes */

#define ON	1
#define OFF	0

/* default dir for help and usage files */
#define T_HELPDIR  "T_HELPDIR"
/* this has changed: this constant now contains the name of the
 * environment-variable in which to look for the help-directory
 */

/* Bert:   (see also volmanip.h) */

#define WORLD  1
#define DATA   2

typedef struct {
  double x, y, z;
  int state;
} Coord;

typedef struct {
  float x, y, z; } CoordF;

typedef struct {
  Coord p1, p2; } LineCoord;

typedef struct {
  Coord p1, p2, p3; } PlaneCoord;

typedef struct {
  double phi, theta, cosphi, sinphi, costheta, sintheta;
  int state;   } Angle;
                                /* phi: in the xy-plane */
                                /* theta: perp. to xy-plane */
/* def of Coord, LineCoord, PlaneCoord, Angle */

#ifdef __cplusplus
extern "C" {
#endif


/* global flags */

extern int	TRACE;	/* trace mode by default */
extern int	VERBOSE;	/* verbose mode by default */
extern int	EXEC;	/* execution mode by default */
extern int      WARNINGS;   /* print warnings by default */
extern int      ERRORS;     /* print errors by default */


/***********************************************/
/*                                             */
/*    error and information printing macros    */
/*                                             */
/***********************************************/

extern const char *T_LINE;

#ifndef T_INFO
#define T_INFO \
        if ( VERBOSE ) \
        { \
        printf("%s", T_LINE); \
        printf("%s", "Tivoli info: "); fprintf
#define INFO_T \
        ; \
        printf("%s", T_LINE);			\
        fflush(stdout); }
#endif

#ifndef T_WARNING
#define T_WARNING \
        if ( WARNINGS ) \
        { \
        (void) fflush  ( stdout ); \
        (void) fprintf ( stdout, "%s", T_LINE );				\
        (void) fprintf ( stdout, "%s", "! Tivoli warning: " ); (void) fprintf
#define WARNING_T \
        ; \
    (void) fprintf( stdout, "%s", T_LINE );	\
        (void) fflush ( stdout ); }
#endif

#ifndef T_ERROR
#define T_ERROR \
        if ( ERRORS ) \
        { \
        (void) fflush  ( stdout ); \
        (void) fprintf ( stdout, "%s", T_LINE ); \
        (void) fprintf ( stdout, "%s", "!! Tivoli error: " ); (void) fprintf
#define ERROR_T \
        ; \
        (void) fprintf( stdout, "%s", T_LINE ); \
        (void) fflush ( stdout ); }
#endif

/***************************/
/*                         */
/*    handy MACROS         */
/*                         */
/***************************/

/*----------------------------------------------------------------------------*/

#define vinilum(a,b) ((a)<(b)?(a):(b))  /* indeed, the minimum of a and b */
#define vaxilum(a,b) ((a)>(b)?(a):(b))  /* indeed, the maximum of a and b */

/*----------------------------------------------------------------------------*/


/***************************/
/*                         */
/*    math CONSTANTS       */
/*                         */
/***************************/

#define E_LNBASIS      	2.718282	/* the venerable e */
#define SQRT2		1.414214	/* sqrt(2) */
#define SQRT3		1.732051	/* sqrt(3) */
#define GOLDEN		1.618034	/* the golden ratio */
#define DEG2RAD  	0.01745329251994  /* convert degrees to radians */
#define RAD2DEG		57.29577951308232 /* convert radians to degrees */

#define PI_2 1.57079632679490     /* pi/2 */
#ifndef PI
#define PI  3.14159265358979      /* pi   */
#endif
#define PIx3_2  4.71238898038469  /* 3*pi/2 */
#define PIx2 6.28318530717959     /* 2*pi */
#define PIx3 9.42477796076938     /* 3*pi */

/*
** LERP(t,a,b) = Linear interpolation between 'a' and 'b' using 't' (0<=t<=1)
** ==> LERP(0) = a, LERP(1) = b
**
** HERP(t,a,b) = Hermite interpolation between 'a' and 'b' using 't' (0<=t<=1)
** ==> HERP(0) = a, HERP'(0) = 0, HERP(1) = b, HERP'(1) = 0
**
** CERP(t,a,b) = Cardinal interpolation between 'a' and 'b' using 't' (0<=t<=1)
** ==> CERP(0) = a, CERP'(0) = 0.5*(b-aa), CERP(1) = b, CERP'(1) = 0.5*(bb-a)
*/
#define LERP(t,a,b)       ((a) + ((b)-(a))*(t))
#define HERP(t,a,b)       ((a) + ((b)-(a))*(t)*(t)*(3.0-(t)-(t)))
#define CERP(t,aa,a,b,bb) ((a) + 0.5*(t)*((b)-(aa)+(t)*(2.0*(aa)-5.0*(a)\
                          +4.0*(b)-(bb)+(t)*((bb)-3.0*(b)+3.0*(a)-(aa)))))

/*
** ABS(t) = Absolute value of 't'
** SGN(t) = Sign value of 't'
** FLOOR(t) = Map 't' to the default integer
** ROUND(t) = Map 't' to the nearest integer
*/
#define ABS(t)          ((t) < 0 ? -(t) : (t))
#define FABS(t)   ABS(t)
#define SGN(t)          ((t) < 0 ? -1 : (t) > 0 ? 1 : 0)
#define FLOOR(t)        ((t) < 0 ? (int) ((t)-1.0) : (int) (t))
#define ROUND(t)        ((t) < 0 ? (int) ((t)-0.5) : (int) ((t)+0.5))


/***************************/
/*                         */
/*    safe ALLOCATION      */
/*                         */
/***************************/


/*----------------------------------------------------------------------------*/
	extern
	void		    *mallocT
/*
	abort program execution if there is no more memory available ;
	otherwise do 'malloc'
*/
/*----------------------------------------------------------------------------*/
(
	size_t		size
);


/*----------------------------------------------------------------------------*/
	extern
	void		    *callocT
/*
	abort program execution if there is no more memory available ;
	otherwise do 'calloc'
*/
/*----------------------------------------------------------------------------*/
(
	size_t		nelem,
	size_t		elsize
);


/*-----------------------------------------------------------------------------*/
        extern
        void             **matAlloc

/*      allocate (if possible) a consecutive memory segment for a matrix by using
        callocT.
*/
/*-----------------------------------------------------------------------------*/
(
                 int lines,
                 int columns,
                 int elem_size
);

/*-----------------------------------------------------------------------------*/
        extern
        void             matFree

/*      Free memory allocated to a matrix.
*/
/*-----------------------------------------------------------------------------*/
(
                 void **   mat
);

/*----------------------------------------------------------------------------*/

char *charMalloc(size_t  length);
/*----------------------------------------------------------------------------*/


/*==============================*
 |                              |
 |  error processing functions  |
 |                              |
 *==============================*/

/*-------------------------------------------------------------------------*/
/*                     E X I T  and co.                                    */
/*-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/
extern void t_in ( char *name );
extern void t_out ( char *name );
extern void t_error_out ( char *name );
extern void t_cl_out ( char *name );
/*-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/
extern void
exitFailure ( char *from, char *reason );
/*-------------------------------------------------------------------------*/

/*-------------------------------------------------------------------------*/
extern void
exitSuccess ( char *from );
/*-------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
        extern
        void exitError ( char *errorText );

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
	extern
	void			printfexit
/*
	abort program execution after printing the name of the function
	where problems occure
*/
/*----------------------------------------------------------------------------*/
(
	char		    *fname
);


/*----------------------------------------------------------------------------*/
	extern
	void			printftrace
/*
	if 'TRACE' is 'ON', print in file "./trace.txt" a trace message ;
	otherwise, do nothing
*/
/*----------------------------------------------------------------------------*/
(
	int			mode,
	char		    *fmt, ...
);

/***************************/
/*                         */
/*    USER I/O- functions  */
/*                         */
/***************************/

/* Ask for a float */
/**********************/
float ask_float(float *def, float *mini, float *maxi, char *comment);

/*-----------------------------------------------*/
/*demande un entier en specifiant eventuellement les bornes (qui
sont toujours actives) selon la valeur de show (true,false)*/
extern int getMyInt( int liminf, int limsup, int show);
/*JFM*/

/*----------------------------------------------------*/
/*demande un float (meme convention que getMyInt)*/
extern float getMyFloat(  float liminf, float limsup, int show);
/*JFM*/

/*-------------------------------------------------*/
/*demande une chaine de caractere en posant une question*/
extern char *getMyName( char *question, int maxsize);
/*JFM*/

/*---------------------------------------------------------*/
/*demande une reponse oui ou non a une question*/
extern int getManicheanRep(char *question);
/*JFM*/


/***************************/
/*                         */
/*    DISK I/O- functions  */
/*                         */
/***************************/

/*----------------------------------------------------------------------------*/
        extern
        void            matPrint
/* Function for printing data arrays in MATLAB format.
*/
/*----------------------------------------------------------------------------*/
(
                char     name[],
                int        rows,
                int        cols,
                double  **array
);

extern int file_open_read(FILE **fp, const char *name);
extern int file_open_write(FILE **fp, const char *name);
extern int file_open_append(FILE **fp, const char *name);
extern int file_close(FILE **fp);

extern int file_write(void *ptr, size_t size, size_t length, FILE *fp);

extern int m_filename(char *base, char *name, int number);
extern int m_point_filename(char *base, char *name, int number);
extern int m_zero_filename(char *base, char *name, int number);

extern char *composeName ( char *baseName, char *ext1, char *ext2 );



/***************************/
/*                         */
/*    ROUND - functions    */
/*                         */
/***************************/

int roundD2I(double f);
long int roundFloat2Long(float f);
int roundFloat2Int(float f);
unsigned char roundFloat2UChar(float f);

/***************************/
/*                         */
/*    MATH - functions    */
/*                         */
/***************************/
void lsincos(double f, double *s, double *c);

/*********************************/
/*                               */
/*    HELP- USAGE - functions    */
/*                               */
/*********************************/
int printUsage(char *usageOf);
/* Searches usage information in a text file "<usageOf>_U.txt".
   path: ./
         or the path indicated in the T_HELPDIR constant
   Output goes to stdout.
   Lines beginning with a '%' are skipped.
*/

int printHelp(char *helpon);
/* Same behaviour as printUsage, but for a file "<helpon>_H.txt".
   After the helptext, "<helpon>_U.txt" is given also.
*/

/*********************************/
/*                               */
/*    Line command utilites      */
/*                               */
/*********************************/
#define testSeeOrWrite(t1,t2)\
    if ( t1 == NO && t2 == NO )\
    {\
        (void)fprintf ( stdout,\
               "You have to specify an output argument (-o and/or -s).\n" );\
         t_cl_out ( argv[0] );\
    }

#define testSeeU8(t1,t2,vol)\
    if ( t1 == NO && t2 == YES && typeOf(vol) != U8BIT )\
    {\
        (void)fprintf ( stdout,\
               "I can see only U8BIT volumes.\n" );\
         t_error_out ( argv[0] );\
    }

#define writeSeeVolume(vol,name,visu,yesifwrite,yesifsee)\
    if ( yesifwrite )\
    {\
        iwriteVolume ( vol, name, "output" );\
    }\
    if ( yesifsee )\
    {\
        if ( typeOf(vol) != U8BIT )\
        {\
            (void)fprintf ( stdout, "Volume is not U8BIT, I cannot show it !\n" );\
        }\
        else\
        {\
             iseeVolume(vol, visu);\
        }\
    }




#ifdef __cplusplus
}
#endif

#ifndef CLOCKS_PER_SEC
#define CLOCKS_PER_SEC  1000000l
#endif

int
fegal(const float	a,
      const float	b);

int
degal(const double	a,
      const double	b);

#endif /* TIVOLI_UTIL_H */
