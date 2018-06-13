/*------------------------------------------------------------------------------
          1         2         3         4         5         6         7
01234567890123456789012345678901234567890123456789012345678901234567890123456789

'volume.h'	cree le 27/04/94 par geraud@ima.enst.fr
			modifie le 17/03/95
			modifie le 28/05/97 par cointepa@ima.enst.fr
			modifie le 24/07/97 by Yann.Cointepas@enst.fr
			   Added GIS.1.0 file format support

			Copyright (c) 1994
			Dept. IMA  TELECOM Paris
			46, rue Barrault
			F-75634 PARIS Cedex 13
			All Rights Reserved

------------------------------------------------------------------------------*/
/** @file volume.h
 * @brief Volume structures (to read first)
 *
 * This file defines the Volume structure and all the functions to handle it.
 *
 * If you are new to TIVOLI you should probably first get familiar with the following basic functions :
 *
     @li I/O : readVolume, writeVolume
     @li Data accessing : data_U8BIT, tab_U8BIT
     @li Memory management : freeVolume, createVolume, duplicateVolumeStructure, allocateVolumeData
     @li Size management : getSize
     @li Tests : testType, testFilled

@warning
   If you compile tivoli for Linux on Intel or for DEC and if you wish
   to access directly to image files saved on different platforms
   you should define GIS_BYTE_ORDER correctly.<BR>
   Examples :
   <UL>
   <LI>SUN : #define GIS_BYTE_ORDER GIS_BYTE_ORDER_SUN
   <LI>Linux Intel (PC) : #define GIS_BYTE_ORDER GIS_BYTE_ORDER_DEC
   <LI>Linux Motorola (Mac) : #define GIS_BYTE_ORDER GIS_BYTE_ORDER_SUN
   <LI>DEC : #define GIS_BYTE_ORDER GIS_BYTE_ORDER_DEC
   <LI>Amiga : #define GIS_BYTE_ORDER GIS_BYTE_ORDER_SUN
   </UL>
*/

#ifndef TIVOLI_VOLUME_H
#define TIVOLI_VOLUME_H

/* 20020326 O.C. - Linux
 *----------------------
 * gcc -O6 -fPIC -I/usr/X11R6/include -I./lib -I./include -I../libstring/include -c -o lib/volume.o lib/volume.c
 * lib/volume.c: In function `readVolume_ima':
 * lib/volume.c:8989: `FLT_MAX' undeclared (first use in this function)
 *----------------------
 * + #include <float.h>
 */

#include <float.h>

#ifndef TIVOLI_UTIL_H
#include <util.h>
#endif

#include <nifti1.h>

/* volume limitation */

/*new: changed from 2048 to 4096 */
#define IMAGE_MAXSIZE	131072		/**< maximum of (Volume*)->size. */
#define BORDER_MAXWIDTH	16		/**< maximum of (Volume*)->borderWidth */
#define NAME_MAXLEN	1024U		/**< maximum of volume name */
#define BIGVOLUME       1000000         /**< for volumes with size3Db > BIGVOLUME
                                           progress can be printed */

/* Different file formats */
/**GIS file format*/
#define GIS_TIVOLI 0
/**GIS file format*/
#define GIS_1_0    1

/* Possible and current byte order */
/**Byte order : big/little endian (Sun sparc, Motorola)*/
#define GIS_BYTE_ORDER_SUN 0
/**Byte order : big/little endian (Dec Alpha, Intel)*/
#define GIS_BYTE_ORDER_DEC 1
/**Byte order : big/little endian (Dec Alpha, Intel)*/
#define GIS_BYTE_ORDER_PC  GIS_BYTE_ORDER_DEC

/* If you compile tivoli for Linux on Intel or for DEC and if you wish
   to access directly to image files saved on different platforms
   you should define GIS_BYTE_ORDER correctly.

   Examples :

   SUN : #define GIS_BYTE_ORDER GIS_BYTE_ORDER_SUN
   Linux Intel (PC) : #define GIS_BYTE_ORDER GIS_BYTE_ORDER_DEC
   Linux Motorola (Mac) : #define GIS_BYTE_ORDER GIS_BYTE_ORDER_SUN
   DEC : #define GIS_BYTE_ORDER GIS_BYTE_ORDER_DEC
   Amiga : #define GIS_BYTE_ORDER GIS_BYTE_ORDER_SUN
   */

#ifndef GIS_BYTE_ORDER
#define GIS_BYTE_ORDER GIS_BYTE_ORDER_PC
/* #define GIS_BYTE_ORDER GIS_BYTE_ORDER_SUN */
#endif


#ifdef __cplusplus
extern "C" {
#endif


/**
* @brief Region of interest structure
*/
typedef struct {
  int  x1, x2, y1, y2, z1, z2, t1, t2, state;
} Roi;

/**
 * @brief Volume structure.

 * This is the base structure of TIVOLI. Use it to declare a Volume.
 * Use functions to access members properly.
 */
typedef struct {
  /* generic attributes : */

  /** Volume id*/
  long id;
  /** Volume name*/
  char name[NAME_MAXLEN+1];

  struct {
    char name[NAME_MAXLEN+1];
    int x;
    int y;
    int z;
    int t;
  } ref; /**< reference image */

  /* intrinsic caracteristics */
  int type;

  struct {
    int x; /* number of voxels per line */
    int y; /* number of lines */
    int z; /* number of slices */
    int t; /* number of images in the sequence */
  } size;

  struct {
    float x; /**< in millimeters */
    float y; /**< in millimeters */
    float z; /**< in millimeters */
    float t; /**< in second */
  } voxelSize;

  /** extrinsic caracteristic */
  int borderWidth;

  /** numerical values of the image contents
   ** size = size3Db (see below)
   */
  char *data;
  /* tab allow to access to the field data as an 3d array tab[slice][row][col] */
  char ***tab;

  /** security flag */
  int state;
  /*
  ** size3Db = (size.x + 2 * borderWidth) * (size.y + 2 * borderWidth)
  **           * (size.z + 2 * borderWidth) * sizeofType(volume)
  */
  size_t size3Db;

#ifdef MEMORYMAP
  char *swapfile;
#endif

  /* GIS format flags */
  char gisFileFormat;      /**< GIS_TIVOLI or GIS_1_0 */
  size_t gisOffset;        /**< GIS format: number of bytes to skip in .ima file to
			      access the point of coordinate (0,0,0) */
  char gisByteOrder;       /**< GIS format: Byte order on file GIS_BYTE_ORDER_SUN or
			      GIS_BYTE_ORDER_DEC. The byte order of the
			      current platform should be correctly defined
			      in GIS_BYTE_ORDER */

  /* nifti/analyse file format add */
  nifti_1_header		*nifti_hdr;
  /* .img/.hdr or .nii */
  unsigned			nitfi_file_extension;
} Volume;


/* macro for volume declaration */

#define DECL_VOLUME(aVolumePtr) Volume *aVolumePtr = NULL

#ifdef MEMORYMAP
extern int memoryMapping;
#endif

/* volume input-output specification in functions */

#define SAME_VOLUME 100
#define NEW_VOLUME  101



/* constant for the field (Volume*)->type */

#define U8BIT  10		/**< 8 bit unsigned integer */
#define S8BIT  11		/**< 8 bit signed integer */
#define U16BIT 20		/**< 16 bit unsigned integer */
#define S16BIT 21		/**< 16 bit signed integer */
#define U32BIT 40		/**< 32 bit unsigned integer */
#define S32BIT 41		/**< 32 bit signed integer */
#define FLOAT  400		/**< single precision real */
#define DOUBLE 800		/**< double precision real */



/* corresponding C declaration type */

typedef unsigned char	U8BIT_t;/**< 8 bit unsigned integer */
typedef char			S8BIT_t;/**< 8 bit signed integer */
typedef unsigned short	U16BIT_t;/**< 16 bit unsigned integer */
typedef short			S16BIT_t;/**< 16 bit signed integer */
  /* typedef unsigned long	U32BIT_t; *//**< 32 bit unsigned integer */
  /* typedef long		S32BIT_t; *//**< 32 bit signed integer */
typedef unsigned int	U32BIT_t;/**< 32 bit unsigned integer */
typedef int			S32BIT_t;/**< 32 bit signed integer */
typedef float			F32BIT_t;
typedef float			FLOAT_t;/**< single precision real */
typedef double			F64BIT_t;
typedef double			DOUBLE_t;/**< double precision real */

typedef char			ANY_t;	/**< any type above */


/* limits for each type */

#define  U8BIT_MIN            0 /**< 8 bit unsigned integer lower limit */
#define  U8BIT_MAX          255 /**< 8 bit unsigned integer upper limit */
#define  S8BIT_MIN         -128 /**< 8 bit signed integer lower limit */
#define  S8BIT_MAX          127 /**< 8 bit signed integer upper limit */
#define U16BIT_MIN            0 /**< 16 bit unsigned integer lower limit */
#define U16BIT_MAX        65535/**< 16 bit unsigned integer upper limit */
#define S16BIT_MIN       -32768/**< 16 bit signed integer lower limit */
#define S16BIT_MAX        32767/**< 16 bit signed integer upper limit */
#define U32BIT_MIN            0/**< 32 bit unsigned integer lower limit */
#define U32BIT_MAX   4294967295U/**< 32 bit unsigned integer upper limit */
#define S32BIT_MIN (-2147483647-1)/**< 32 bit signed integer lower limit */
#define S32BIT_MAX   2147483647/**< 32 bit signed integer upper limit */

#define FLOAT_MIN   (-FLT_MAX)/**< single precision real lower limit*/
#define FLOAT_MAX     FLT_MAX/**< single precision real upper limit*/
#define FLOAT_EPS     FLT_EPSILON/**< single precision real precision*/
#define DOUBLE_MIN  (-DBL_MAX)/**< double precision real lower limit*/
#define DOUBLE_MAX    DBL_MAX/**< double precision real upper limit*/
#define DOUBLE_EPS    DBL_EPSILON/**< double precision real precision*/

/*============================================*
 |                                            |
 |  printing functions (print to stdout)      |
 |                                            |
 *============================================*/
/*----------------------------------------------------------------------------*/
  /*
    print the nxnxn voxel values nearest from the first volume point (mode 0)
    or from the last point (mode 1)
  */
  extern void
  printfvolume(Volume*	volume,
	       int	n,
	       int	mode);

int
floategal(const float       a,
          const float       b);


/*----------------------------------------------------------------------------*/
  void
  printfRoi(Roi		roi,
	    FILE	*fp);

/*----------------------------------------------------------------------------*/
  extern void
  printVolumeList();


/*=================================================*
 |                                                 |
 |  testing functions  (return 0 (PB) on failure)  |
 |                                                 |
 *================================================*/

/*----------------------------------------------------------------------------*/
/*
	test if volume contains significant information (caracteristics and/or
	data) ; this function is called by all the following ones
*/
  extern int
  testExist(Volume	    *volume);

  /* idem for n volumes */
  extern int
  testExistN(int	n,
	     ...);

  /**
   * @brief Tests if the volume type is the one specified

   * @see testEqualType, typeOf
   */

/*----------------------------------------------------------------------------*/
  extern
  int			testType
/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			type
);


/*----------------------------------------------------------------------------*/
	extern
	int			testTypeN
/*
	test if the volume type is one of the list of n
*/
/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			n, ...
);


  /**
   * @brief	Tests if the volume contains significant data
   */
/*----------------------------------------------------------------------------*/
	extern
	int			testFilled
/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


	/**
	 * @brief   Tests if both volumes have the same dimensions (size.x .y .z and .t)
         *
	 * @see getSize, sizexOf, sizeyOf, sizezOf
	 */
/*----------------------------------------------------------------------------*/
	extern
	int			testEqualSize
/*----------------------------------------------------------------------------*/
(
	Volume	    *volume1,
	Volume	    *volume2
);


/*----------------------------------------------------------------------------*/
	extern
	int			testEqual2DSize
/*
	test if both volumes have the same dimensions (size.x and .y)
*/
/*----------------------------------------------------------------------------*/
(
	Volume	    *volume1,
	Volume	    *volume2
);


  /**
   * @brief   Tests if both volumes have the same type

   * @see testType, typeOf
   */
/*----------------------------------------------------------------------------*/
	extern
	int			testEqualType
/*----------------------------------------------------------------------------*/
(
	Volume	    *volume1,
	Volume	    *volume2
);


  /**
   * @brief Tests if both volumes have the same border width
   */

/*----------------------------------------------------------------------------*/
	extern
	int			testEqualBorderWidth
/*----------------------------------------------------------------------------*/
(
	Volume	    *volume1,
	Volume	    *volume2
);


  /**
   * @brief Tests if the slice index is in the volume (# between -bW and z+bW-1)
   */
/*----------------------------------------------------------------------------*/
	extern
	int			testSlice
/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			slice
);


  /**
   *@brief Tests if the image contains the slice (# between 0 and size.z-1)
   */
/*----------------------------------------------------------------------------*/
	extern
	int			testImageSlice
/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			slice
);


  /**
   * @brief verify that the volume contains correct information
   */
/*----------------------------------------------------------------------------*/
	extern
	int			verifyAll
/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/*----------------------------------------------------------------------------*/
extern
int                        verifyRoi
(
        Roi *roi,
        Volume *volume
);
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
extern
int testLevelInType ( double level, int type );

/*----------------------------------------------------------------------------*/


/*=====================================================*
 |                                                     |
 |  getting functions : access to the volume structure |
 |                                                     |
 *=====================================================*/

/*----------------------------------------------------------------------------*/
	extern
	void			getSize2

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int		    *sizeX,
	int		    *sizeY
);
#define getSize3 getSize

  /**
   * @brief Get the volume size in all 3 directions.
   *
   * @param sizeX A pointer to a variable that will contain the size in X direction
   * @param sizeY A pointer to a variable that will contain the size in Y direction
   * @param sizeZ A pointer to a variable that will contain the size in Z direction

   * @see sizexOF, sizeyOf, sizezOf, testEqualSize
   */
/*----------------------------------------------------------------------------*/
	extern
	void			getSize

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int		    *sizeX,
	int		    *sizeY,
	int		    *sizeZ
);
/*----------------------------------------------------------------------------*/
	extern
	void			getSize4

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int		    *sizeX,
	int		    *sizeY,
	int		    *sizeZ,
	int		    *sizeT
);


  /**
   * @brief Get the volume size in X direction.
   *
   * @see getSize, testEqualSize
   */
/*----------------------------------------------------------------------------*/

	int			sizexOf

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);
  /**
   * @brief Get the volume size in Y direction.
   *
   * @see getSize, testEqualSize
   */
/*----------------------------------------------------------------------------*/

	int			sizeyOf

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);
  /**
   * @brief Get the volume size in Z direction.
   *
   * @see getSize, testEqualSize
   */
/*----------------------------------------------------------------------------*/

	int			sizezOf

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);
  /**
   * @brief Get the volume size in T direction.
   */
/*----------------------------------------------------------------------------*/

	int			sizetOf

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);



/*----------------------------------------------------------------------------*/
	extern
	void			getVoxelSize2

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	float	    *voxelSizeX,
	float	    *voxelSizeY
);
#define getVoxelSize3 getVoxelSize
  /**
   * @brief Get the voxel size in all 3 directions.
   *
   * @param voxelSizeX A pointer to a variable that will contain the voxel size in X direction
   * @param voxelSizeY A pointer to a variable that will contain the voxel size in Y direction
   * @param voxelSizeZ A pointer to a variable that will contain the voxel size in Z direction
   */
/*----------------------------------------------------------------------------*/
	extern
	void			getVoxelSize

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	float	    *voxelSizeX,
	float	    *voxelSizeY,
	float	    *voxelSizeZ
);
/*----------------------------------------------------------------------------*/
	extern
	void			getVoxelSize4

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	float	    *voxelSizeX,
	float	    *voxelSizeY,
	float	    *voxelSizeZ,
	float	    *voxelSizeT
);


/*----------------------------------------------------------------------------*/
	extern
	int			borderWidthOf

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);

  /**
   * @brief Get the type of the volume (U8, U16, S16...)

   */
  extern
  int			typeOf(Volume	    *volume);

/*----------------------------------------------------------------------------*/
extern
char		    *typeName

/*----------------------------------------------------------------------------*/
(
	int			type
 );
/*----------------------------------------------------------------------------*/
extern
char		    *type2string

/*----------------------------------------------------------------------------*/
(
	int			type
 );
/*----------------------------------------------------------------------------*/
extern
int
string2type(const char*	name);


/* for use in util.h */
/*----------------------------------------------------------------------------*/
extern
char		    *stateName
/*----------------------------------------------------------------------------*/
(
	int			state
);


/*----------------------------------------------------------------------------*/
        extern
	void	    getMinMax

/*----------------------------------------------------------------------------*/
(
        Volume     *volume,
        double     *minVal,
        double     *maxVal
 );



/*----------------------------------------------------------------------------*/
	extern
	int	    *boundingBox

/*----------------------------------------------------------------------------*/
(
	Volume	    *Volume
);


/*----------------------------------------------------------------------------*/
        extern
	int	    gisGetFileFormat

/*----------------------------------------------------------------------------*/
(
        Volume     *volume
);

/*----------------------------------------------------------------------------*/
        extern
	size_t	    gisGetOffset

/*----------------------------------------------------------------------------*/
(
        Volume     *volume
);

/*----------------------------------------------------------------------------*/
        extern
	char	    gisGetByteOrder

/*----------------------------------------------------------------------------*/
(
        Volume     *volume
);


/*============================================================*
 |                                                            |
 |  setting functions : change values in the volume structure |
 |                                                            |
 *============================================================*/



/*----------------------------------------------------------------------------*/
	extern
	void			setSize2
/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			nx,
	int			ny
);
#define setSize3 setSize
/*----------------------------------------------------------------------------*/
	extern
	void			setSize
/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			nx,
	int			ny,
	int			nz
);
/*----------------------------------------------------------------------------*/
	extern
	void			setSize4

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			nx,
	int			ny,
	int			nz,
	int			nt
);



/*----------------------------------------------------------------------------*/
	extern
	void			setSizex

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			nx
);
/*----------------------------------------------------------------------------*/
	extern
	void			setSizey

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			ny
);
/*----------------------------------------------------------------------------*/
	extern
	void			setSizez

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			nz
);
/*----------------------------------------------------------------------------*/
	extern
	void			setSizet

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			nt
);



/*----------------------------------------------------------------------------*/
	extern
	void			setVoxelSize2
/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	float		dx,
	float		dy
);
#define setVoxelSize3 setVoxelSize
/*----------------------------------------------------------------------------*/
	extern
	void			setVoxelSize
/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	float		dx,
	float		dy,
	float		dz
);
/*----------------------------------------------------------------------------*/
	extern
	void			setVoxelSize4

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	float		dx,
	float		dy,
	float		dz,
	float		dt
);


/*----------------------------------------------------------------------------*/
	extern
	void			setType

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			type
);


/*----------------------------------------------------------------------------*/
	extern
	void			setBorderWidth

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			borderWidth
);


/*----------------------------------------------------------------------------*/
	extern
	void			setVolumeLevel

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	double			level
);


/*----------------------------------------------------------------------------*/
	extern
	void			setImageLevel

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	double			level
);


/*----------------------------------------------------------------------------*/
	extern
	void			setBorderLevel

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	double			level
);



/*----------------------------------------------------------------------------*/
        extern
	void	    gisSetFileFormat

/*----------------------------------------------------------------------------*/
(
        Volume     *volume, int
);

/*----------------------------------------------------------------------------*/
        extern
	void	    gisSetOffset

/*----------------------------------------------------------------------------*/
(
        Volume     *volume, size_t
);

/*----------------------------------------------------------------------------*/
        extern
	void	    gisSetByteOrder

/*----------------------------------------------------------------------------*/
(
        Volume     *volume, char
);




/*====================*
 |                    |
 |  offset functions  |
 |                    |
 *====================*/


/*----------------------------------------------------------------------------*/
	extern
	int			sizeofType

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);
/* the following function should disappear in the future: */
/*----------------------------------------------------------------------------*/
	extern
	long			offset_forANYtype_BetweenPoint

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/*----------------------------------------------------------------------------*/
	extern
	long			offsetFirstPoint

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/*----------------------------------------------------------------------------*/
	extern
	long			offsetLine

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/*----------------------------------------------------------------------------*/
	extern
	long			offsetPointBetweenLine

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/*----------------------------------------------------------------------------*/
	extern
	long			offsetSlice

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/*----------------------------------------------------------------------------*/
	extern
	long			offsetLineBetweenSlice

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/*----------------------------------------------------------------------------*/
	extern
	long			offsetVolume

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/*----------------------------------------------------------------------------*/
	extern
	long			offsetImage

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/*----------------------------------------------------------------------------*/

	long int	    *offsetBox

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			connexity
);







/*============================*
 |                            |
 |  data accessing functions  |
 |                            |
 *============================*/

/*
	Rq : we can not prevent with 'lint' the six occurances of the message
	'pointer cast may result in improper alignment' caused by these functions
*/

/**
 * @brief Data accessing (first method)
 *
 * There are two different methods to access the volume data. The first one uses this function (the second uses tab_U8BIT).
 * It returns a 1D array containing the volume data. The following code example shows how to use this function :
 * @code
        int nx, ny, nz, ix, iy, iz;
        long oPbL, oLbS;
        U8BIT_t *data;
	Volume *volume
        if (!testType ( volume, U8BIT ) ||
                       !testFilled ( volume ) )
                        goto abort;

        getSize ( volume, &nx, &ny, &nz );
        oPbL = offsetPointBetweenLine ( volume );
        oLbS = offsetLineBetweenSlice ( volume );

        data = data_U8BIT ( volume ) + offsetFirstPoint ( volume );

	for(iz=0;iz<nz;iz++) {
		for(iy=0;iy<ny;iy++) {
			for(ix=0;ix<nx;ix++) {
			  ... Do something ...
			  data++;
			}
			data += oPbL;
		}
		data += oLbS;
	}
    @endcode

    * @see offsetPointBetweenLine, offsetLineBetweenSlice, offsetFirstPoint, getSize, tab_U8BIT
    */
/*----------------------------------------------------------------------------*/
	extern
	U8BIT_t	    *data_U8BIT

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/**
 * @brief Data accessing (first method)
 *
 * There are two different methods to access the volume data. The first one uses this function (the second uses tab_S8BIT).
 * It returns a 1D array containing the volume data. The following code example shows how to use this function :
 * @code
        int nx, ny, nz, ix, iy, iz;
        long oPbL, oLbS;
        S8BIT_t *data;
	Volume *volume
        if (!testType ( volume, S8BIT ) ||
                       !testFilled ( volume ) )
                        goto abort;

        getSize ( volume, &nx, &ny, &nz );
        oPbL = offsetPointBetweenLine ( volume );
        oLbS = offsetLineBetweenSlice ( volume );

        data = data_S8BIT ( volume ) + offsetFirstPoint ( volume );

	for(iz=0;iz<nz;iz++) {
		for(iy=0;iy<ny;iy++) {
			for(ix=0;ix<nx;ix++) {
			  ... Do something ...
			  data++;
			}
			data += oPbL;
		}
		data += oLbS;
	}
    @endcode

    * @see offsetPointBetweenLine, offsetLineBetweenSlice, offsetFirstPoint, getSize, tab_S8BIT
    */
/*----------------------------------------------------------------------------*/
	extern
	S8BIT_t	    *data_S8BIT

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);

/**
 * @brief Data accessing (first method)
 *
 * There are two different methods to access the volume data. The first one uses this function (the second uses tab_U16BIT).
 * It returns a 1D array containing the volume data. The following code example shows how to use this function :
 * @code
        int nx, ny, nz, ix, iy, iz;
        long oPbL, oLbS;
        U16BIT_t *data;
	Volume *volume
        if (!testType ( volume, U16BIT ) ||
                       !testFilled ( volume ) )
                        goto abort;

        getSize ( volume, &nx, &ny, &nz );
        oPbL = offsetPointBetweenLine ( volume );
        oLbS = offsetLineBetweenSlice ( volume );

        data = data_U16BIT ( volume ) + offsetFirstPoint ( volume );

	for(iz=0;iz<nz;iz++) {
		for(iy=0;iy<ny;iy++) {
			for(ix=0;ix<nx;ix++) {
			  ... Do something ...
			  data++;
			}
			data += oPbL;
		}
		data += oLbS;
	}
    @endcode

    * @see offsetPointBetweenLine, offsetLineBetweenSlice, offsetFirstPoint, getSize, tab_U16BIT
    */

/*----------------------------------------------------------------------------*/
	extern
	U16BIT_t	    *data_U16BIT

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/**
 * @brief Data accessing (first method)
 *
 * There are two different methods to access the volume data. The first one uses this function (the second uses tab_S16BIT).
 * It returns a 1D array containing the volume data. The following code example shows how to use this function :
 * @code
        int nx, ny, nz, ix, iy, iz;
        long oPbL, oLbS;
        S16BIT_t *data;
	Volume *volume
        if (!testType ( volume, S16BIT ) ||
                       !testFilled ( volume ) )
                        goto abort;

        getSize ( volume, &nx, &ny, &nz );
        oPbL = offsetPointBetweenLine ( volume );
        oLbS = offsetLineBetweenSlice ( volume );

        data = data_S16BIT ( volume ) + offsetFirstPoint ( volume );

	for(iz=0;iz<nz;iz++) {
		for(iy=0;iy<ny;iy++) {
			for(ix=0;ix<nx;ix++) {
			  ... Do something ...
			  data++;
			}
			data += oPbL;
		}
		data += oLbS;
	}
    @endcode

    * @see offsetPointBetweenLine, offsetLineBetweenSlice, offsetFirstPoint, getSize, tab_S16BIT
    */
/*----------------------------------------------------------------------------*/
	extern
	S16BIT_t	    *data_S16BIT

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/**
 * @brief Data accessing (first method)
 *
 * There are two different methods to access the volume data. The first one uses this function (the second uses tab_U32BIT).
 * It returns a 1D array containing the volume data. The following code example shows how to use this function :
 * @code
        int nx, ny, nz, ix, iy, iz;
        long oPbL, oLbS;
        U32BIT_t *data;
	Volume *volume
        if (!testType ( volume, U32BIT ) ||
                       !testFilled ( volume ) )
                        goto abort;

        getSize ( volume, &nx, &ny, &nz );
        oPbL = offsetPointBetweenLine ( volume );
        oLbS = offsetLineBetweenSlice ( volume );

        data = data_U32BIT ( volume ) + offsetFirstPoint ( volume );

	for(iz=0;iz<nz;iz++) {
		for(iy=0;iy<ny;iy++) {
			for(ix=0;ix<nx;ix++) {
			  ... Do something ...
			  data++;
			}
			data += oPbL;
		}
		data += oLbS;
	}
    @endcode

    * @see offsetPointBetweenLine, offsetLineBetweenSlice, offsetFirstPoint, getSize, tab_U32BIT
    */
/*----------------------------------------------------------------------------*/
	extern
	U32BIT_t	    *data_U32BIT

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/**
 * @brief Data accessing (first method)
 *
 * There are two different methods to access the volume data. The first one uses this function (the second uses tab_S32BIT).
 * It returns a 1D array containing the volume data. The following code example shows how to use this function :
 * @code
        int nx, ny, nz, ix, iy, iz;
        long oPbL, oLbS;
        S32BIT_t *data;
	Volume *volume
        if (!testType ( volume, S32BIT ) ||
                       !testFilled ( volume ) )
                        goto abort;

        getSize ( volume, &nx, &ny, &nz );
        oPbL = offsetPointBetweenLine ( volume );
        oLbS = offsetLineBetweenSlice ( volume );

        data = data_S32BIT ( volume ) + offsetFirstPoint ( volume );

	for(iz=0;iz<nz;iz++) {
		for(iy=0;iy<ny;iy++) {
			for(ix=0;ix<nx;ix++) {
			  ... Do something ...
			  data++;
			}
			data += oPbL;
		}
		data += oLbS;
	}
    @endcode

    * @see offsetPointBetweenLine, offsetLineBetweenSlice, offsetFirstPoint, getSize, tab_S32BIT
    */
/*----------------------------------------------------------------------------*/
	extern
	S32BIT_t	   *data_S32BIT

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/**
 * @brief Data accessing (first method)
 *
 * There are two different methods to access the volume data. The first one uses this function (the second uses tab_FLOAT).
 * It returns a 1D array containing the volume data. The following code example shows how to use this function :
 * @code
        int nx, ny, nz, ix, iy, iz;
        long oPbL, oLbS;
        FLOAT_t *data;
	Volume *volume
        if (!testType ( volume, FLOAT ) ||
                       !testFilled ( volume ) )
                        goto abort;

        getSize ( volume, &nx, &ny, &nz );
        oPbL = offsetPointBetweenLine ( volume );
        oLbS = offsetLineBetweenSlice ( volume );

        data = data_FLOAT ( volume ) + offsetFirstPoint ( volume );

	for(iz=0;iz<nz;iz++) {
		for(iy=0;iy<ny;iy++) {
			for(ix=0;ix<nx;ix++) {
			  ... Do something ...
			  data++;
			}
			data += oPbL;
		}
		data += oLbS;
	}
    @endcode

    * @see offsetPointBetweenLine, offsetLineBetweenSlice, offsetFirstPoint, getSize, tab_FLOAT
    */
/*----------------------------------------------------------------------------*/
	extern
	FLOAT_t	    *data_FLOAT

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/**
 * @brief Data accessing (first method)
 *
 * There are two different methods to access the volume data. The first one uses this function (the second uses tab_DOUBLE).
 * It returns a 1D array containing the volume data. The following code example shows how to use this function :
 * @code
        int nx, ny, nz, ix, iy, iz;
        long oPbL, oLbS;
        DOUBLE_t *data;
	Volume *volume
        if (!testType ( volume, DOUBLE ) ||
                       !testFilled ( volume ) )
                        goto abort;

        getSize ( volume, &nx, &ny, &nz );
        oPbL = offsetPointBetweenLine ( volume );
        oLbS = offsetLineBetweenSlice ( volume );

        data = data_DOUBLE ( volume ) + offsetFirstPoint ( volume );

	for(iz=0;iz<nz;iz++) {
		for(iy=0;iy<ny;iy++) {
			for(ix=0;ix<nx;ix++) {
			  ... Do something ...
			  data++;
			}
			data += oPbL;
		}
		data += oLbS;
	}
    @endcode

    * @see offsetPointBetweenLine, offsetLineBetweenSlice, offsetFirstPoint, getSize, tab_DOUBLE
    */
/*----------------------------------------------------------------------------*/
	extern
	DOUBLE_t	    *data_DOUBLE

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/*----------------------------------------------------------------------------*/
	extern
	ANY_t	    *data_ANY

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/**
 * @brief Data accessing (second method)
 *
 * There are two different methods to access the volume data. The second one uses this function (the first uses data_U8BIT).
 * It returns a 3D array containing the volume data. The following code example shows how to use this function :
 * @code
	Volume *volume;
        int ix, iy, iz;
        U8BIT_t ***data;

        if (!testType ( volume, U8BIT ) ||
                       !testFilled ( volume ) )
                        goto abort;

        data = tab_U8BIT ( volume );

	Do something with data[iz][iy][ix]...
    @endcode

    * @see data_U8BIT
    */

/*----------------------------------------------------------------------------*/
	extern
	U8BIT_t    ***tab_U8BIT

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/**
 * @brief Data accessing (second method)
 *
 * There are two different methods to access the volume data. The second one uses this function (the first uses data_S8BIT).
 * It returns a 3D array containing the volume data. The following code example shows how to use this function :
 * @code
	Volume *volume;
        int ix, iy, iz;
        S8BIT_t ***data;

        if (!testType ( volume, S8BIT ) ||
                       !testFilled ( volume ) )
                        goto abort;

        data = tab_S8BIT ( volume );

	Do something with data[iz][iy][ix]...
    @endcode

    * @see data_S8BIT
    */
/*----------------------------------------------------------------------------*/
	extern
	S8BIT_t    ***tab_S8BIT

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);

/**
 * @brief Data accessing (second method)
 *
 * There are two different methods to access the volume data. The second one uses this function (the first uses data_U16BIT).
 * It returns a 3D array containing the volume data. The following code example shows how to use this function :
 * @code
	Volume *volume;
        int ix, iy, iz;
        U16BIT_t ***data;

        if (!testType ( volume, U16BIT ) ||
                       !testFilled ( volume ) )
                        goto abort;

        data = tab_U16BIT ( volume );

	Do something with data[iz][iy][ix]...
    @endcode

    * @see data_U16BIT
    */
/*----------------------------------------------------------------------------*/
	extern
	U16BIT_t    ***tab_U16BIT

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);

/**
 * @brief Data accessing (second method)
 *
 * There are two different methods to access the volume data. The second one uses this function (the first uses data_S16BIT).
 * It returns a 3D array containing the volume data. The following code example shows how to use this function :
 * @code
	Volume *volume;
        int ix, iy, iz;
        S16BIT_t ***data;

        if (!testType ( volume, S16BIT ) ||
                       !testFilled ( volume ) )
                        goto abort;

        data = tab_S16BIT ( volume );

	Do something with data[iz][iy][ix]...
    @endcode

    * @see data_S16BIT
    */
/*----------------------------------------------------------------------------*/
	extern
	S16BIT_t    ***tab_S16BIT

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);

/**
 * @brief Data accessing (second method)
 *
 * There are two different methods to access the volume data. The second one uses this function (the first uses data_U32BIT).
 * It returns a 3D array containing the volume data. The following code example shows how to use this function :
 * @code
	Volume *volume;
        int ix, iy, iz;
        U32BIT_t ***data;

        if (!testType ( volume, U32BIT ) ||
                       !testFilled ( volume ) )
                        goto abort;

        data = tab_U32BIT ( volume );

	Do something with data[iz][iy][ix]...
    @endcode

    * @see data_U32BIT
    */
/*----------------------------------------------------------------------------*/
	extern
	U32BIT_t    ***tab_U32BIT

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);

/**
 * @brief Data accessing (second method)
 *
 * There are two different methods to access the volume data. The second one uses this function (the first uses data_S32BIT).
 * It returns a 3D array containing the volume data. The following code example shows how to use this function :
 * @code
	Volume *volume;
        int ix, iy, iz;
        S32BIT_t ***data;

        if (!testType ( volume, S32BIT ) ||
                       !testFilled ( volume ) )
                        goto abort;

        data = tab_S32BIT ( volume );

	Do something with data[iz][iy][ix]...
    @endcode

    * @see data_S32BIT
    */
/*----------------------------------------------------------------------------*/
	extern
	S32BIT_t   ***tab_S32BIT

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);

/**
 * @brief Data accessing (second method)
 *
 * There are two different methods to access the volume data. The second one uses this function (the first uses data_FLOAT).
 * It returns a 3D array containing the volume data. The following code example shows how to use this function :
 * @code
	Volume *volume;
        int ix, iy, iz;
        FLOAT_t ***data;

        if (!testType ( volume, FLOAT ) ||
                       !testFilled ( volume ) )
                        goto abort;

        data = tab_FLOAT ( volume );

	Do something with data[iz][iy][ix]...
    @endcode

    * @see data_FLOAT
    */
/*----------------------------------------------------------------------------*/
	extern
	FLOAT_t    ***tab_FLOAT

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);

/**
 * @brief Data accessing (second method)
 *
 * There are two different methods to access the volume data. The second one uses this function (the first uses data_DOUBLE).
 * It returns a 3D array containing the volume data. The following code example shows how to use this function :
 * @code
	Volume *volume;
        int ix, iy, iz;
        DOUBLE_t ***data;

        if (!testType ( volume, DOUBLE ) ||
                       !testFilled ( volume ) )
                        goto abort;

        data = tab_DOUBLE ( volume );

	Do something with data[iz][iy][ix]...
    @endcode

    * @see data_DOUBLE
    */
/*----------------------------------------------------------------------------*/
	extern
	DOUBLE_t    ***tab_DOUBLE

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);






/*==========================*
 |                          |
 |  input-output functions  |
 |                          |
 *==========================*/



/* Macros for volume reading and writing in a command line program
   (depends on argv[0]). */
#define ireadVolume(filename,borderWidth,vol,refname)\
   (void)fprintf ( stdout, "Reading %s volume '%s'...", refname, filename );\
   (void)fflush ( stdout );\
   vol = readVolume(filename, borderWidth);\
   if ( vol == NULL )\
   {\
       printf ( "\nProblem during reading volume '%s'!\n", filename );\
       t_error_out ( argv[0] );\
   }\
   else (void)fprintf ( stdout, "\b\b\b\nReading finished.\n" );\
   (void)fflush ( stdout )

#define iwriteVolume(vol,filename,refname)\
   (void)fprintf ( stdout, "Writing %s volume to '%s'...", refname, filename );\
   (void)fflush ( stdout );\
   writeVolume( vol, filename );\
   (void)fprintf ( stdout, "\b\b\b finished.\n" );\
   (void)fflush ( stdout )

#define iseeVolume(vol,visu)\
    (void)fprintf ( stdout,\
           "Seeing output volume (strike h in window to get help)...\n" );\
    (void)fflush ( stdout );\
    visu = openXTWindow ( vol, 256 );\
    while ( waitforXTEvent () );\
    (void)fprintf ( stdout, "\b\b\b achieved\n" );


/* the following function existed earlier with other parameters
   please change them !  The dimension info is now returned as
   a volume structure
   The old function is renamed to readVolumeDim_OLD but should
   disappear in the future
*/
/*----------------------------------------------------------------------------*/

	Volume *readVolumeDim

/*----------------------------------------------------------------------------*/
(
	const char    *filename,
        Roi         *roi
);

/*----------------------------------------------------------------------------*/

	Volume *readVolumeDim_return

/*----------------------------------------------------------------------------*/
(
	const char    *filename,
        Roi         *roi
);

/*----------------------------------------------------------------------------*/
	extern
	int			readVolumeDim_OLD

/*----------------------------------------------------------------------------*/
(
	const char    *filename,
	int		    *nx,
	int		    *ny,
	int		    *nz,
	int		    *type,
	float	    *dx,
	float	    *dy,
	float	    *dz
);


/*----------------------------------------------------------------------------*/
  extern Volume*
  readVolume_ima(const char	*filename,
		 int		borderWidth);

  unsigned
  testFileExist(const char	*filename);

  /**
   * @brief Read a volume from disk
   *
   * This function is able to read gzipped IMA files.
   *
   * @param filename File name (either without extension or with .dim, .ima or .ima.gz extension).
   * @param borderWidth The border width
   *
   * @return A pointer to the volume

   * @see writeVolume
   */
  extern
  Volume	    *readVolume(const char	*filename,
				int		borderWidth);

/*----------------------------------------------------------------------------*/
extern
Volume	        *readVolumeRoi				/*NEW*/

/*----------------------------------------------------------------------------*/
(
	const char    *filename,
	int			borderWidth,
        Roi     roi
 );


/*----------------------------------------------------------------------------*/
extern
void			writeVolumeDim		       	/*NEW*/

/*----------------------------------------------------------------------------*/
(
 Volume	    *volume,
 const char    *filename
 );

  /**
   * @brief Write a volume to disk
   *
   * @param volume The volume to save
   * @param filename File name (either without extension or with .dim, .ima, .nii, .hdr, .imge or .ima.gz extension).
   *
   * @see readVolume
   */
  extern void
  writeVolume(Volume	    *volume,
	      const char    *filename);

  /**
   * @brief Write a ima/dim volume to disk
   *
   * @param volume The volume to save
   * @param filename File name (either without extension or with .dim, .ima or .ima.gz extension).
   *
   * @see readVolume
   */
  extern void
  writeVolume_ima(Volume	*volume,
		  const char    *filename);

  /**
   * @brief Writes a series of slices
   *
   * @param volume The volume to save
   * @param firstSlice The slice number to start from.  Slices go from 0 to n-1 independent of the presence of a border.
   * @param nSlices The number of consecutive slices to write.
   * @param filename File name (either without extension or with .dim, .ima or .ima.gz extension).
   *
   */
/*----------------------------------------------------------------------------*/

	void			writeVolumeSlices

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			firstSlice, /* slices go from 0 to n-1 */
	int			nSlices, /* the number of consecutive slices */
	const char    *filename
 );

/*----------------------------------------------------------------------------*/
	extern
	void			writeVolumeSlice

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			slice,   /* slices go from 0 to nz-1 */
	const char    *filename
);


/*----------------------------------------------------------------------------*/
	extern
	void                    writeSliceDim3D

/*----------------------------------------------------------------------------*/
(
        Volume     *volume,
	const char    *filename
);



/*==========================*
 |                          |
 |  volume creators         |
 |                          |
 *==========================*/




#define declareVolumeStructure3 declareVolumeStructure

  /**
   * @brief Declare a new volume.
   *
   * Declare a new volume structure. This function does not allocate the volume data (use allocateVolumeData or createVolume).
   *
   * @see allocateVolumeData, createVolume, duplicateVolumeStructure, freeVolume
   */
/*----------------------------------------------------------------------------*/
	extern
	Volume	    *declareVolumeStructure
/*----------------------------------------------------------------------------*/
(
	int			xSize,
	int			ySize,
	int			zSize,
	float		xVoxelSize,
	float		yVoxelSize,
	float		zVoxelSize,
	int			type,
	const char    *name,
	int			borderWidth
);

  /**
   * @brief Declare a new volume (4D).
   *
   * Declare a new volume structure. This function does not allocate the volume data (use allocateVolumeData or createVolume4).
   *
   * @see allocateVolumeData, createVolume4, duplicateVolumeStructure, freeVolume
   */
/*----------------------------------------------------------------------------*/
	extern
	Volume	    *declareVolumeStructure4

/*----------------------------------------------------------------------------*/
(
	int			xSize,
	int			ySize,
	int			zSize,
	int			tSize,
	float		xVoxelSize,
	float		yVoxelSize,
	float		zVoxelSize,
	float		tVoxelSize,
	int			type,
	const char    *name,
	int			borderWidth
);

  /**
   * @brief Duplicate a volume structure.
   *
   * This function does not allocate nor copy the volume data (use allocateVolumeData or copyVolume).
   *
   * @param volume The volume to duplicate
   * @param name The name of the new volume
   *
   * @see allocateVolumeData, createVolume, freeVolume, copyVolume
   */

/*----------------------------------------------------------------------------*/
	extern
	Volume	    *duplicateVolumeStructure

/*----------------------------------------------------------------------------*/	(
	Volume	    *volume,
	const char    *name
);

  /**
   * @brief Allocate volume data
   *
   * You must create the volume structure first (with declareVolumeStructure or duplicateVolumeStructure)
   *
   * @param volume The volume to duplicate
   * @param name The name of the new volume
   *
   * @see declareVolumeStructure, duplicateVolumeStructure, createVolume, freeVolume, copyVolume
   */
/*----------------------------------------------------------------------------*/
	extern
	void			allocateVolumeData

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


  /**
   * @brief Create a new volume
   *
   * This function declares a volume structure and allocate the volume data
   *
   * @see declareVolumeStructure, duplicateVolumeStructure, allocateVolumeData, freeVolume, copyVolume
   */
/*----------------------------------------------------------------------------*/
	extern
	Volume	    *createVolume		/*  prefer new version  */

/*----------------------------------------------------------------------------*/
(
	int			xSize,
	int			ySize,
	int			zSize,
	float		xVoxelSize,
	float		yVoxelSize,
	float		zVoxelSize,
	int			type,
	const char    *name,
	int			borderWidth
);
/*----------------------------------------------------------------------------*/
	extern
	Volume	    *createVolume4		/*NEW*/

/*----------------------------------------------------------------------------*/
(
	int			xSize,
	int			ySize,
	int			zSize,
	int			tSize,
	float		xVoxelSize,
	float		yVoxelSize,
	float		zVoxelSize,
	float		tVoxelSize,
	int			type,
	const char    *name,
	int			borderWidth
);


  /**
   * @brief Copy volume data
   *
   * This function does not copy the volume structure (use copyVolume)
   *
   * @param volumeR The volume to be read
   * @param volumeW The volume to be written
   * @see copyVolume
   */
/*----------------------------------------------------------------------------*/
	extern
	void			transferVolumeData

/*----------------------------------------------------------------------------*/
(
	Volume	    *volumeR,
	Volume	    *volumeW
);


/*----------------------------------------------------------------------------*/
	extern
	void			appendVolume

/*----------------------------------------------------------------------------*/
(
	const char    *filename,
	Volume	    *volume
);


  /**
   * @brief Copy volume
   *
   * Copy both structure and data
   *
   * @param volume The volume to be copied
   * @param name The name of the new volume
   * @see transferVolumeData
   */
/*----------------------------------------------------------------------------*/
	extern
	Volume	    *copyVolume

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	const char    *name
);




/*----------------------------------------------------------------------------*/
	extern
	void			copySlice

/*----------------------------------------------------------------------------*/
(
	Volume	    *volumeSource,
	int			sliceSource,
	Volume	    *volumeDestination,
	int			sliceDestination
);



/*----------------------------------------------------------------------------*/

	void			addBorder

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


/*----------------------------------------------------------------------------*/
	extern
	Volume	    *extractSlice

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			slice
);

/*----------------------------------------------------------------------------*/
        extern
	Volume	    *extractRoi

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
        Roi      *roi
 );

/*----------------------------------------------------------------------------*/
        extern
	Volume	    *roiInv

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
        int x,
        int y,
        int z,
        int nx,
        int ny,
        int nz

 );

/*----------------------------------------------------------------------------*/

	Volume	    *vol3DtoVol2D

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume,
	int			first,
	int			step,
	int			number,
	int			width
);


  /**
   * @brief Free volume data
   *
   * But not the volume structure (use freeVolume)
   *
   * @see freeVolume, allocateVolumeData
   *
   */
/*----------------------------------------------------------------------------*/
	extern
	void			freeVolumeData

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);


  /**
   * @brief Free volume
   *
   * Free both volume and data
   *
   * @see freeVolumeData, creatVolume
   *
   */
/*----------------------------------------------------------------------------*/
	extern
	void			freeVolume

/*----------------------------------------------------------------------------*/
(
	Volume	    *volume
);

void
listLIST(char	*str);


/*----------------------------------------------------------------------------*/
	extern
	Volume                 *axial2coronal
/*----------------------------------------------------------------------------*/
(
	Volume *volume
);



/*----------------------------------------------------------------------------*/
	extern
	Volume                 *coronal2axial
/*----------------------------------------------------------------------------*/
(
	Volume *volume
);


/*----------------------------------------------------------------------------*/
	extern
	Volume                 *axial2sagittal
/*----------------------------------------------------------------------------*/
(
	Volume *volume
);


/*----------------------------------------------------------------------------*/
	extern
	Volume                 *sagittal2axial
/*----------------------------------------------------------------------------*/
(
	Volume *volume
);




/*----------------------------------------------------------------------------*/
extern
double             goodCast (  double d,  int type  );
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
extern
void               gisInverseByteOrder (  Volume *  );
/*----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------*/
extern
void               writeRoi ( char *filename, int nx, int ny, int nz, int x1, int x2, int y1, int y2, int z1, int z2 );
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
extern
void               readRoi ( char *filename, int *nx, int *ny, int *nz, int *x1, int *x2, int *y1, int *y2, int *z1, int *z2 );
/*----------------------------------------------------------------------------*/

void
freeTivoliImage(Volume  **img);

#ifdef __cplusplus
}
#endif

#endif /* TIVOLI_VOLUME_H */
