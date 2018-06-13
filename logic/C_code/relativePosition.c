#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <volume.h>

#include <relativePosition.h>

#define ffmin(a,b) ((a)<(b)?(a):(b))
#define ffmax(a,b) ((a)>(b)?(a):(b))



static void neighborhoodX(int* , int , int );
static void neighborhoodY(int* , int , int );
static void neighborhoodZ(int* , int , int );


static void initNeighborhood(int* , int* , int* , int , int , int , int );
static void incNeighborhood(int* , int , int );
static void decNeighborhood(int* , int , int );



/*	==================================================================	*/
/*   angleVisible.c: computation of the "fuzzy landscape" representing	*/
/*		the degrees of satisfaction of the relative position 					*/
/*    in some direction.																*/
/*		Method using tabulation of angles											*/
/*		et search over the only points belonging to the support.				*/
/*					Isabelle Bloch - ENST Paris

										*/
/*     Modified by Olivier Colliot, Oct. 2002
       Add parameter k : the biggest angle such that mu!=0
       mu=max(0,1-beta/k) instead of mu=max(0,1-2beta/pi) (i.e. k=pi/2) */
/*	==================================================================	*/


Volume *angleVisible
(
 Volume	   *volume1,
 double		alpha1,
 double		alpha2,
 double k
 )
{
  Volume	    *resuVolume, *tabAngle;
  U8BIT_t ***info1, ***dist;
  DOUBLE_t ***angle;
  int nx, ny, nz, ix, iy, iz, ixp, iyp, izp, i;
  /* 	int inv, n; */
  double cos_angle, cos1, sin1, cos2, sin2, angle_min;
  double diffx, diffy, diffz, tmp, dirX, dirY, dirZ ;
  /* 	double pi = acos(-1.) ; */
  int card, *X, *Y, *Z;

  printftrace ( IN, "angleVisible" );

  /*	Control tests		*/
  /*	=============		*/

  if (  !testType ( volume1, U8BIT ) ||
	!testFilled ( volume1 ) )
    goto abort;

  /*	Bounds calculation	*/
  /*	==================	*/

  getSize ( volume1, &nx, &ny, &nz );
  if ( nz == 1 && alpha2 != 0. )
    {
      fprintf ( stderr,"warning: angle not compatible with image dimension\n" );
    }

  /*	Tab initialisation	*/
  /*	==================	*/

  info1 = tab_U8BIT ( volume1 );
  resuVolume = duplicateVolumeStructure ( volume1, "resu" );
  setBorderWidth ( resuVolume, 0 );
  allocateVolumeData ( resuVolume );
  dist = tab_U8BIT ( resuVolume );

  tabAngle = declareVolumeStructure ( (2*nx+1), (2*ny+1), (2*nz+1), 1., 1., 1., DOUBLE, "angle", 0 );
  allocateVolumeData ( tabAngle );
  angle = tab_DOUBLE ( tabAngle );

  X = (int*)malloc((nx*ny*nz)*sizeof(int));
  Y = (int*)malloc((nx*ny*nz)*sizeof(int));
  Z = (int*)malloc((nx*ny*nz)*sizeof(int));
  card = 0;

  printf("Search of object points\n");

  /*	Listing of object points (in the fuzzy case, points of the support)	*/
  /*	===================================================================	*/

  for(iz=0;iz<nz;iz++)
    for(iy=0;iy<ny;iy++)
      for(ix=0;ix<nx;ix++)
	{
	  if ( info1[iz][iy][ix] != 0 )
	    {
	      X[card] = ix ;
	      Y[card] = iy ;
	      Z[card] = iz ;
	      card++;
	    }
	}

  fprintf ( stderr, "cardinality of the object : %d\n", card) ;

  printf("Tabulation of angles\n");

  /*	Tabulation of angles	*/
  /*	====================	*/

  if ( EXEC )
    {

      cos1 = cos(alpha1);
      sin1 = sin(alpha1);
      cos2 = cos(alpha2);
      sin2 = sin(alpha2);
      dirX = cos2 * cos1 ;
      dirY = cos2 * sin1 ;
      dirZ = sin2 ;
      fprintf ( stderr, "dirX, dirY, dirZ : %f %f %f\n", (float)(dirX), (float)(dirY), (float)(dirZ) );


      for(iz=0;iz<(2*nz+1);iz++)
	for(iy=0;iy<(2*ny+1);iy++)
	  for(ix=0;ix<(2*nx+2);ix++)
	    {
	      diffx = (double)(ix) - (double)(nx) ;
	      diffy = (double)(iy) - (double)(ny) ;
	      diffz = (double)(iz) - (double)(nz) ;
	      tmp = sqrt ( diffx * diffx + diffy * diffy + diffz * diffz ) ;
	      cos_angle = ( dirX * diffx + dirY * diffy + dirZ * diffz ) / tmp ;
	      angle[iz][iy][ix] = acos ( cos_angle ) ;
	    }

      angle[nz][ny][nx] = 0 ;
      printf("Beginning computation...\n");

      /*	Visibility angle computation	*/
      /*	============================	*/

      for(izp=0;izp<nz;izp++)
	{
	  fprintf ( stderr, "%d\n", izp );
	  for(iyp=0;iyp<ny;iyp++)
	    {
	      fprintf ( stderr, "\b\b\b%3d", iyp );
	      for(ixp=0;ixp<nx;ixp++)
		{

		  /*       if ( info1[izp][iyp][ixp] != 0 )
			   dist[izp][iyp][ixp] = 255 ;
			   else
		  */
		  {

		    /*			weighting angle by membership value	*/
		    /*			angle_min = pi ;
		      for(i=0;i<card;i++)
		      {
		      ix = ixp - X[i] + nx ;
		      iy = iyp - Y[i] + ny ;
		      iz = izp - Z[i] + nz ;
		      tmp = 255. * angle[iz][iy][ix] / (double)(info1[Z[i]][Y[i]][X[i]]) ;
		      angle_min = ( tmp < angle_min ? tmp : angle_min ) ;
		      }
		      tmp = 1. - ( 2. * angle_min / pi );
		      if ( tmp < 0. ) tmp = 0. ;
		      dist[izp][iyp][ixp] = (int)( 255. * tmp ) ;
		    */

		    /*			combination of f(angle) by membership value	*/
		    angle_min = 0. ;
		    for(i=0;i<card;i++)
		      {
			ix = ixp - X[i] + nx ;
			iy = iyp - Y[i] + ny ;
			iz = izp - Z[i] + nz ;
			tmp = angle[iz][iy][ix] ;
			tmp = 1. -  tmp / k;
			if ( tmp < 0. ) tmp = 0. ;
			tmp = tmp * (double)(info1[Z[i]][Y[i]][X[i]]) / 255. ;
			angle_min = ( tmp > angle_min ? tmp : angle_min ) ;
		      }
		    dist[izp][iyp][ixp] = (int)( 255. * angle_min ) ;
		  }
		}
	    }
	}

    }

  printftrace ( OUT, "" );
  return( resuVolume );


 abort :
  printfexit ( "angleVisible\n" );
  return ( NULL );
}



/*	==================================================================	*/
/*   angleVisiblePropag.c: computation of the "fuzzy landscape" 			*/
/*		representing the degrees of satisfaction of the relative				*/
/*		position in a given direction.												*/
/*		Method using tabulation of angles											*/
/*		and a propagation	algorithm (provides an approximation).				*/
/*								Isabelle Bloch - ENST Paris							*/
/*     Modified by Olivier Colliot, Oct. 2002
       Add parameter k : the biggest angle such that mu!=0
       mu=max(0,1-beta/k) instead of mu=max(0,1-2beta/pi) (i.e. k=pi/2) */
/*	==================================================================	*/


Volume *angleVisiblePropag
(
 Volume	   *volume1,
 double		alpha1,
 double		alpha2,
 double k
 )
{
  Volume	    *resuVolume, *tabAngle, *resuX, *resuY, *resuZ;
  U8BIT_t ***info1, ***dist;
  S16BIT_t ***origX, ***origY, ***origZ;
  DOUBLE_t ***angle;
  int nx, ny, nz, ix, iy, iz, ixp, iyp, izp, i, connexity;
  /* 	int	inv, n; */
  int nx2, ny2, nz2;
  int *ixV, *iyV, *izV, X, Y, Z;
  /* 	double cos_angle; */
  double	cos1, sin1, cos2, sin2, angle_min;
  double diffx, diffy, diffz, tmp, dirX, dirY, dirZ ;
  /* 	double pi = acos(-1.) ; */
  double deux_pi_255, un_255, tmp2;
  /*deux_pi_255 = 2 / pi / 255.;*/
  deux_pi_255= 1 / k / 255.;
  un_255 = 1. / 255.;

  printftrace ( IN, "angleVisiblePropag" );

  /*	Control tests		*/
  /*	=============		*/

  if (  !testType ( volume1, U8BIT ) ||
	!testFilled ( volume1 ) )
    goto abort;

  /*	Bounds calculation	*/
  /*	==================	*/

  getSize ( volume1, &nx, &ny, &nz );
  if ( nz == 1 && alpha2 != 0. )
    {
      fprintf ( stderr,"warning: angle not compatible with image dimension\n" );
    }

  if ( nz == 1 )
    connexity = 8 ;
  else
    connexity = 26 ;

  /*	Tab initialisation	*/
  /*	==================	*/

  fprintf(stderr,"init\n");

  fprintf(stderr,"%d %d %d \n", nx,ny,nz);

  setBorderWidth ( volume1, 1 );
  info1 = tab_U8BIT ( volume1 );
  setBorderLevel ( volume1, 0 );
  resuVolume = duplicateVolumeStructure ( volume1, "resu" );
  setBorderWidth ( resuVolume, 0 );
  dist = tab_U8BIT ( resuVolume );

  nx2 = 2*nx+1;
  ny2 = 2*ny+1;
  nz2 = 2*nz+1;
  tabAngle = declareVolumeStructure ( nx2, ny2, nz2, 1., 1., 1., DOUBLE, "angle", 0 );
  angle = tab_DOUBLE ( tabAngle );

  resuX = duplicateVolumeStructure ( volume1, "resuX" );
  setType ( resuX, S16BIT );
  origX = tab_S16BIT ( resuX );
  resuY = duplicateVolumeStructure ( resuX, "resuY" );
  origY = tab_S16BIT ( resuY );
  resuZ = duplicateVolumeStructure ( resuX, "resuZ" );
  origZ = tab_S16BIT ( resuZ );

  if ( EXEC )
    {
      ixV = (int*)malloc((connexity+1)*sizeof(int)) ;
      iyV = (int*)malloc((connexity+1)*sizeof(int)) ;
      izV = (int*)malloc((connexity+1)*sizeof(int)) ;

      printf("Tabulation of angles\n");

      /*	Tabulation of angles	*/
      /*	====================	*/

      cos1 = cos(alpha1);
      sin1 = sin(alpha1);
      cos2 = cos(alpha2);
      sin2 = sin(alpha2);
      dirX = cos2 * cos1 ;
      dirY = cos2 * sin1 ;
      dirZ = sin2 ;
      fprintf ( stderr, "dirX, dirY, dirZ : %f %f %f\n", (float)(dirX), (float)(dirY), (float)(dirZ) );


      for(iz=0;iz<nz2;++iz)
	for(iy=0;iy<ny2;++iy)
	  for(ix=0;ix<nx2;++ix)
	    {
	      diffx = (double)(ix) - (double)(nx) ;
	      diffy = (double)(iy) - (double)(ny) ;
	      diffz = (double)(iz) - (double)(nz) ;
	      tmp = sqrt ( diffx * diffx + diffy * diffy + diffz * diffz ) ;
	      if ( tmp == 0. )
		{
		  angle[iz][iy][ix] = 0 ;
		  if ( dirX * diffx + dirY * diffy + dirZ * diffz != 0 )
		    {
		      fprintf ( stderr, "erreur : division par 0\n" );
		    }
		}
	      else
		{
		  angle[iz][iy][ix] = acos ( ( dirX * diffx + dirY * diffy + dirZ * diffz ) / tmp );
		}
	    }

      angle[nz][ny][nx] = 0 ;

      /*	Initialization of origin points	*/
      /*	===============================	*/

      for(iz=-1;iz<=nz;++iz)
	for(iy=-1;iy<=ny;++iy)
	  for(ix=-1;ix<=nx;++ix)
	    {
	      if ( info1[iz][iy][ix] != 0 )
		{
		  origX[iz][iy][ix] = ix ;
		  origY[iz][iy][ix] = iy ;
		  origZ[iz][iy][ix] = iz ;
		}
	      else
		{
		  origX[iz][iy][ix] = -1 ;
		  origY[iz][iy][ix] = -1 ;
		  origZ[iz][iy][ix] = -1 ;
		}
	    }


      printf("Beginning the computation...\n");

      /*	Visibility angle computation	*/
      /*	============================	*/

      /*	for(iy=0;iy<ny;iy++)
	{
	for(ix=0;ix<nx;ix++)
	fprintf ( stderr, "%4d", info1[0][iy][ix] );
	fprintf ( stderr, "\n") ;
	}
	fprintf ( stderr, "\n") ;
	fprintf ( stderr, "\n") ;

	for(iy=0;iy<ny;iy++)
	{
	for(ix=0;ix<nx;ix++)
	fprintf ( stderr, "%4d", iy*ny + ix );
	fprintf ( stderr, "\n") ;
	}
	fprintf ( stderr, "\n") ;
	fprintf ( stderr, "\n") ;
      */

      for(izp=0;izp<nz;++izp)
	{
	  /*		fprintf ( stderr, "%d\n", izp );	*/
	  neighborhoodZ(izV,izp,connexity);
	  for(iyp=0;iyp<ny;++iyp)
	    {
	      /*		fprintf ( stderr, "\b\b\b%3d", iyp );	*/
	      neighborhoodY(iyV,iyp,connexity);
	      for(ixp=0;ixp<nx;++ixp)
		{
		  neighborhoodX(ixV,ixp,connexity);

		  /*			combination of f(angle) by membership value	*/
		  angle_min = 0. ;
		  for(i=0;i<=connexity;++i)
		    {
		      X = origX[izV[i]][iyV[i]][ixV[i]] ;
		      Y = origY[izV[i]][iyV[i]][ixV[i]] ;
		      Z = origZ[izV[i]][iyV[i]][ixV[i]] ;
		      if ( X != -1 && Y != -1 && Z != -1 )
			{
			  ix = ixp - X + nx ;
			  iy = iyp - Y + ny ;
			  iz = izp - Z + nz ;
			  tmp = ( un_255 - deux_pi_255 * angle[iz][iy][ix] )
			    * (double)(info1[Z][Y][X]);
			  if ( tmp >= angle_min )
			    {
			      if ( tmp < 0. ) tmp = 0. ;
			      angle_min = tmp ;
			      origX[izp][iyp][ixp] = X ;
			      origY[izp][iyp][ixp] = Y ;
			      origZ[izp][iyp][ixp] = Z ;
			    }
			}
		    }

		  tmp2 = 255. * angle_min;
		  if ( tmp2 >= 255.50001 )
		    {
		      fprintf ( stderr, "depassement sup \n" );
		      tmp2 = 255.;
		    }
		  else
		    if ( tmp2 < 0. )
		      {
			fprintf ( stderr, "depassement inf \n" );
			tmp2 = 0.;
		      }
		    else
		      {
			tmp2 = floor( 0.49999 + tmp2 );
		      }
		  dist[izp][iyp][ixp] = (U8BIT_t)( tmp2 ) ;
		  /*			fprintf ( stderr, "%3d %2d", origX[izp][iyp][ixp], origY[izp][iyp][ixp] );	*/
		}
	      /*		fprintf ( stderr, "\n" );	*/
	    }
	}


      fprintf ( stderr, "\n" );
      fprintf ( stderr, "\n" );

      for(izp=nz-1;izp>=0;--izp)
	{
	  /*		fprintf ( stderr, "%d\n", izp );	*/
	  neighborhoodZ(izV,izp,connexity);
	  for(iyp=ny-1;iyp>=0;--iyp)
	    {
	      /*		fprintf ( stderr, "\b\b\b%3d", iyp );	*/
	      neighborhoodY(iyV,iyp,connexity);
	      for(ixp=nx-1;ixp>=0;--ixp)
		{
		  neighborhoodX(ixV,ixp,connexity);

		  angle_min = (double)(dist[izp][iyp][ixp])/255. ;
		  for(i=0;i<=connexity;++i)
		    {
		      X = origX[izV[i]][iyV[i]][ixV[i]] ;
		      Y = origY[izV[i]][iyV[i]][ixV[i]] ;
		      Z = origZ[izV[i]][iyV[i]][ixV[i]] ;
		      if ( X != -1 && Y != -1 && Z != -1 )
			{
			  ix = ixp - X + nx ;
			  iy = iyp - Y + ny ;
			  iz = izp - Z + nz ;
			  tmp = ( un_255 - deux_pi_255 * angle[iz][iy][ix] )
			    * (double)(info1[Z][Y][X]);
			  if ( tmp >= angle_min )
			    {
			      if ( tmp < 0. ) tmp = 0. ;
			      angle_min = tmp ;
			      origX[izp][iyp][ixp] = X ;
			      origY[izp][iyp][ixp] = Y ;
			      origZ[izp][iyp][ixp] = Z ;
			    }
			}
		    }

		  tmp2 = 255. * angle_min;
		  if ( tmp2 >= 255.50001 )
		    {
		      fprintf ( stderr, "depassement sup \n" );
		      tmp2 = 255.;
		    }
		  else
		    if ( tmp2 < 0. )
		      {
			fprintf ( stderr, "depassement inf \n" );
			tmp2 = 0.;
		      }
		    else
		      {
			tmp2 = floor( 0.49999 + tmp2 );
		      }
		  dist[izp][iyp][ixp] = (U8BIT_t)( tmp2 ) ;
		  /*			fprintf ( stderr, "%4d", (origY[izp][iyp][ixp]*(ny-1) + origY[izp][iyp][ixp]) );	*/
		}
	      /*			fprintf ( stderr, "\n" );	*/
	    }
	}

      free(ixV);
      free(iyV);
      free(izV);
    }

  freeVolume ( resuX ) ;
  freeVolume ( resuY ) ;
  freeVolume ( resuZ ) ;
  freeVolume ( tabAngle ) ;

  printftrace ( OUT, "" );
  return( resuVolume );


 abort :
  printfexit ( "angleVisiblePropag\n" );
  return ( NULL );
}

/*	==================================================================	*/

static void neighborhoodX(int* ixV, int ixp, int connexity)
{
  ixV[0] = ixp ;
  ixV[1] = ixp - 1 ;
  ixV[2] = ixp ;
  ixV[3] = ixp + 1 ;
  ixV[4] = ixp - 1 ;
  ixV[5] = ixp + 1 ;
  ixV[6] = ixp - 1 ;
  ixV[7] = ixp ;
  ixV[8] = ixp + 1 ;

  if(connexity == 26 )
    {
      ixV[9] = ixp - 1 ;
      ixV[10] = ixp ;
      ixV[11] = ixp + 1 ;
      ixV[12] = ixp - 1 ;
      ixV[13] = ixp ;
      ixV[14] = ixp + 1 ;
      ixV[15] = ixp - 1 ;
      ixV[16] = ixp ;
      ixV[17] = ixp + 1 ;

      ixV[18] = ixp - 1 ;
      ixV[19] = ixp ;
      ixV[20] = ixp + 1 ;
      ixV[21] = ixp - 1 ;
      ixV[22] = ixp ;
      ixV[23] = ixp + 1 ;
      ixV[24] = ixp - 1 ;
      ixV[25] = ixp ;
      ixV[26] = ixp + 1 ;

    }
}

/*	==================================================================	*/

static void neighborhoodY(int* iyV, int iyp, int connexity)
{
  iyV[0] = iyp ;
  iyV[1] = iyp - 1 ;
  iyV[2] = iyp - 1 ;
  iyV[3] = iyp - 1 ;
  iyV[4] = iyp ;
  iyV[5] = iyp ;
  iyV[6] = iyp + 1 ;
  iyV[7] = iyp + 1 ;
  iyV[8] = iyp + 1 ;

  if(connexity == 26 )
    {
      iyV[9] = iyp - 1 ;
      iyV[10] = iyp - 1 ;
      iyV[11] = iyp - 1;
      iyV[12] = iyp ;
      iyV[13] = iyp ;
      iyV[14] = iyp ;
      iyV[15] = iyp + 1 ;
      iyV[16] = iyp + 1 ;
      iyV[17] = iyp + 1 ;

      iyV[18] = iyp - 1 ;
      iyV[19] = iyp - 1 ;
      iyV[20] = iyp - 1;
      iyV[21] = iyp ;
      iyV[22] = iyp ;
      iyV[23] = iyp ;
      iyV[24] = iyp + 1 ;
      iyV[25] = iyp + 1 ;
      iyV[26] = iyp + 1 ;
    }
}

/*	==================================================================	*/

static void neighborhoodZ(int* izV, int izp, int connexity)
{
  int i;

  for(i=0;i<=8;++i)
    izV[i] = izp ;

  if(connexity == 26 )
    {
      for(i=9;i<=17;++i)
	izV[i] = izp - 1 ;
      for(i=18;i<=26;++i)
	izV[i] = izp + 1 ;
    }
}



/*	==================================================================	*/
/*   angleVisiblePropag2.c: computation of the "fuzzy landscape" 			*/
/*		representing the degrees of satisfaction of the relative				*/
/*		position in a given direction.												*/
/*		Method using tabulation of angles											*/
/*		and a propagation	algorithm (provides an approximation).				*/
/*		This second version works with floats and not doubles					*/
/*		and uses a faster access to the neighborhood points.					*/
/*								Isabelle Bloch - ENST Paris							*/
/*     Modified by Olivier Colliot, Oct. 2002
       Add parameter k : the biggest angle such that mu!=0
       mu=max(0,1-beta/k) instead of mu=max(0,1-2beta/pi) (i.e. k=pi/2) */
/*	==================================================================	*/


Volume*
angleVisiblePropag2(Volume*	volume1,
		    float	alpha1,
		    float	alpha2,
		    float	k)
{
  Volume			*resuVolume, *tabAngle, *resuX, *resuY, *resuZ;
  U8BIT_t			***info1, ***dist;
  S16BIT_t			***origX, ***origY, ***origZ;
  FLOAT_t			***angle;
  int				nx, ny, nz, ix, iy, iz, ixp, iyp, izp, i, connexity;
  /* 	int	 inv, n; */
  int				nx2, ny2, nz2;
  int				oldBorderWidth = 0;
  int				*ixV, *iyV, *izV, X, Y, Z;
  /* 	float cos_angle; */
  float				cos1, sin1, cos2, sin2, angle_min;
  float				diffx, diffy, diffz, tmp, dirX, dirY, dirZ ;
  /* 	float pi = acos(-1.) ; */
  float				deux_pi_255, un_255;
  /*deux_pi_255 = 2 / pi / 255.;*/
  deux_pi_255 = 1 / k / 255.;
  un_255 = 1. / 255.;

  printftrace (IN, "angleVisiblePropag2");

  /*	Control tests		*/
  /*	=============		*/

  if (!testType(volume1, U8BIT) ||
      !testFilled(volume1))
    goto abort;

  /*	Bounds calculation	*/
  /*	==================	*/

  getSize (volume1, &nx, &ny, &nz);
  if (nz == 1 && alpha2 != 0.)
    fprintf (stderr,"warning: angle not compatible with image dimension\n");

  if (nz == 1)
    connexity = 8 ;
  else
    connexity = 26 ;

  /*	Tab initialisation	*/
  /*	==================	*/

  fprintf(stderr,"init\n");
  fprintf(stderr,"%d %d %d \n", nx,ny,nz);

  oldBorderWidth = volume1->borderWidth;

  setBorderWidth (volume1, 1);
  info1 = tab_U8BIT (volume1);
  setBorderLevel (volume1, 0);
  resuVolume = duplicateVolumeStructure (volume1, "resu");
  setBorderWidth (resuVolume, 0);
  dist = tab_U8BIT (resuVolume);

  nx2 = 2*nx+1;
  ny2 = 2*ny+1;
  nz2 = 2*nz+1;
  tabAngle = declareVolumeStructure (nx2, ny2, nz2, 1., 1., 1., FLOAT, "angle", 0);
  angle = tab_FLOAT (tabAngle);

  resuX = duplicateVolumeStructure (volume1, "resuX");
  setType (resuX, S16BIT);
  origX = tab_S16BIT (resuX);
  resuY = duplicateVolumeStructure (resuX, "resuY");
  origY = tab_S16BIT (resuY);
  resuZ = duplicateVolumeStructure (resuX, "resuZ");
  origZ = tab_S16BIT (resuZ);

  if (EXEC)
    {
      ixV = (int*)malloc((connexity+1)*sizeof(int)) ;
      iyV = (int*)malloc((connexity+1)*sizeof(int)) ;
      izV = (int*)malloc((connexity+1)*sizeof(int)) ;

      printf("Tabulation of angles\n");

      /*	Tabulation of angles	*/
      /*	====================	*/

      cos1 = cos(alpha1);
      sin1 = sin(alpha1);
      cos2 = cos(alpha2);
      sin2 = sin(alpha2);
      dirX = cos2 * cos1 ;
      dirY = cos2 * sin1 ;
      dirZ = sin2 ;
      fprintf (stderr, "dirX, dirY, dirZ : %f %f %f\n", (float)(dirX), (float)(dirY), (float)(dirZ));


      for(iz=0;iz<nz2;++iz)
	for(iy=0;iy<ny2;++iy)
	  for(ix=0;ix<nx2;++ix)
	    {
	      diffx = (float)(ix) - (float)(nx) ;
	      diffy = (float)(iy) - (float)(ny) ;
	      diffz = (float)(iz) - (float)(nz) ;
	      tmp = sqrt (diffx * diffx + diffy * diffy + diffz * diffz) ;
	      if (tmp == 0.)
		{
		  angle[iz][iy][ix] = 0 ;
		}
	      else
		{
		  angle[iz][iy][ix] = acos ((dirX * diffx + dirY * diffy + dirZ * diffz) / tmp);
		}
	    }

      angle[nz][ny][nx] = 0 ;

      /*	Initialization of origin points	*/
      /*	===============================	*/

      for(iz=-1;iz<=nz;++iz)
	for(iy=-1;iy<=ny;++iy)
	  for(ix=-1;ix<=nx;++ix)
	    {
	      if (info1[iz][iy][ix] != 0)
		{
		  origX[iz][iy][ix] = ix ;
		  origY[iz][iy][ix] = iy ;
		  origZ[iz][iy][ix] = iz ;
		}
	      else
		{
		  origX[iz][iy][ix] = -1 ;
		  origY[iz][iy][ix] = -1 ;
		  origZ[iz][iy][ix] = -1 ;
		}
	    }


      printf("Beginning of computation...\n");

      /*	Visibility angle computation	*/
      /*	============================	*/

      /*	for(iy=0;iy<ny;iy++)
	{
	for(ix=0;ix<nx;ix++)
	fprintf (stderr, "%4d", info1[0][iy][ix]);
	fprintf (stderr, "\n") ;
	}
	fprintf (stderr, "\n") ;
	fprintf (stderr, "\n") ;

	for(iy=0;iy<ny;iy++)
	{
	for(ix=0;ix<nx;ix++)
	fprintf (stderr, "%4d", iy*ny + ix);
	fprintf (stderr, "\n") ;
	}
	fprintf (stderr, "\n") ;
	fprintf (stderr, "\n") ;
      */

      initNeighborhood(ixV,iyV,izV,0,0,0,connexity);
      for(izp=0;izp<nz;++izp)
	{
	  /*		fprintf (stderr, "%d\n", izp);	*/
	  for(iyp=0;iyp<ny;++iyp)
	    {
	      /*		fprintf (stderr, "\b\b\b%3d", iyp);	*/
	      for(ixp=0;ixp<nx;++ixp)
		{

		  /*			combination of f(angle) by membership value	*/
		  angle_min = 0. ;
		  for(i=0;i<=connexity;++i)
		    {
		      X = origX[izV[i]][iyV[i]][ixV[i]] ;
		      Y = origY[izV[i]][iyV[i]][ixV[i]] ;
		      Z = origZ[izV[i]][iyV[i]][ixV[i]] ;
		      if (X != -1 && Y != -1 && Z != -1)
			{
			  ix = ixp - X + nx ;
			  iy = iyp - Y + ny ;
			  iz = izp - Z + nz ;
			  tmp = (un_255 - deux_pi_255 * angle[iz][iy][ix])
			    * (float)(info1[Z][Y][X]);
			  if (tmp >= angle_min)
			    {
			      if (tmp < 0.) tmp = 0. ;
			      angle_min = tmp ;
			      origX[izp][iyp][ixp] = X ;
			      origY[izp][iyp][ixp] = Y ;
			      origZ[izp][iyp][ixp] = Z ;
			    }
			}
		    }
		  dist[izp][iyp][ixp] = (U8BIT_t)(floor(0.49999 + 255. * angle_min)) ;
		  /*			fprintf (stderr, "%3d %2d", origX[izp][iyp][ixp], origY[izp][iyp][ixp]);	*/
		  incNeighborhood(ixV, nx, connexity);
		}
	      incNeighborhood(iyV, ny, connexity);
	      /*		fprintf (stderr, "\n");	*/
	    }
	  incNeighborhood(izV, nz, connexity);
	}


      fprintf (stderr, "\n");
      fprintf (stderr, "\n");

      initNeighborhood(ixV,iyV,izV,nx-1,ny-1,nz-1,connexity);
      for(izp=nz-1;izp>=0;--izp)
	{
	  /*		fprintf (stderr, "%d\n", izp);	*/
	  for(iyp=ny-1;iyp>=0;--iyp)
	    {
	      /*		fprintf (stderr, "\b\b\b%3d", iyp);	*/
	      for(ixp=nx-1;ixp>=0;--ixp)
		{
		  angle_min = (float)(dist[izp][iyp][ixp])/255. ;
		  for(i=0;i<=connexity;++i)
		    {
		      X = origX[izV[i]][iyV[i]][ixV[i]] ;
		      Y = origY[izV[i]][iyV[i]][ixV[i]] ;
		      Z = origZ[izV[i]][iyV[i]][ixV[i]] ;
		      if (X != -1 && Y != -1 && Z != -1)
			{
			  ix = ixp - X + nx ;
			  iy = iyp - Y + ny ;
			  iz = izp - Z + nz ;
			  tmp = (un_255 - (deux_pi_255 * angle[iz][iy][ix])) * (float)(info1[Z][Y][X]);
			  if (tmp >= angle_min)
			    {
			      if (tmp < 0.) tmp = 0. ;
			      angle_min = tmp ;
			      origX[izp][iyp][ixp] = X ;
			      origY[izp][iyp][ixp] = Y ;
			      origZ[izp][iyp][ixp] = Z ;
			    }
			}
		    }
		  dist[izp][iyp][ixp] = (U8BIT_t)(floor(0.49999 + 255. * angle_min)) ;
		  /*			fprintf (stderr, "%4d", (origY[izp][iyp][ixp]*(ny-1) + origY[izp][iyp][ixp]));	*/
		  decNeighborhood(ixV, nx, connexity);
		}
	      decNeighborhood(iyV, ny, connexity);
	      /*		fprintf (stderr, "\n");	*/
	    }
	  decNeighborhood(izV, nz, connexity);
	}

      free(ixV);
      free(iyV);
      free(izV);
    }

  freeVolume(resuX); free(resuX);
  freeVolume(resuY); free(resuY);
  freeVolume(resuZ); free(resuZ);
  freeVolume(tabAngle); free(tabAngle);

  // --gyf 19/09/2009 : re-set the border to the original width
  setBorderWidth(volume1, oldBorderWidth);

  printftrace (OUT, "");
  return(resuVolume);


 abort :
  printfexit ("angleVisiblePropag2\n");
  return (NULL);
}



/*	==================================================================	*/

static void initNeighborhood(int* ixV, int* iyV, int* izV, int ixp, int iyp, int izp, int connexity)
{
  int i;

  ixV[0] = ixp ;
  ixV[1] = ixp - 1 ;
  ixV[2] = ixp ;
  ixV[3] = ixp + 1 ;
  ixV[4] = ixp - 1 ;
  ixV[5] = ixp + 1 ;
  ixV[6] = ixp - 1 ;
  ixV[7] = ixp ;
  ixV[8] = ixp + 1 ;

  if(connexity == 26)
    {
      ixV[9] = ixp - 1 ;
      ixV[10] = ixp ;
      ixV[11] = ixp + 1 ;
      ixV[12] = ixp - 1 ;
      ixV[13] = ixp ;
      ixV[14] = ixp + 1 ;
      ixV[15] = ixp - 1 ;
      ixV[16] = ixp ;
      ixV[17] = ixp + 1 ;

      ixV[18] = ixp - 1 ;
      ixV[19] = ixp ;
      ixV[20] = ixp + 1 ;
      ixV[21] = ixp - 1 ;
      ixV[22] = ixp ;
      ixV[23] = ixp + 1 ;
      ixV[24] = ixp - 1 ;
      ixV[25] = ixp ;
      ixV[26] = ixp + 1 ;

    }

  iyV[0] = iyp ;
  iyV[1] = iyp - 1 ;
  iyV[2] = iyp - 1 ;
  iyV[3] = iyp - 1 ;
  iyV[4] = iyp ;
  iyV[5] = iyp ;
  iyV[6] = iyp + 1 ;
  iyV[7] = iyp + 1 ;
  iyV[8] = iyp + 1 ;

  if(connexity == 26)
    {
      iyV[9] = iyp - 1 ;
      iyV[10] = iyp - 1 ;
      iyV[11] = iyp - 1;
      iyV[12] = iyp ;
      iyV[13] = iyp ;
      iyV[14] = iyp ;
      iyV[15] = iyp + 1 ;
      iyV[16] = iyp + 1 ;
      iyV[17] = iyp + 1 ;

      iyV[18] = iyp - 1 ;
      iyV[19] = iyp - 1 ;
      iyV[20] = iyp - 1;
      iyV[21] = iyp ;
      iyV[22] = iyp ;
      iyV[23] = iyp ;
      iyV[24] = iyp + 1 ;
      iyV[25] = iyp + 1 ;
      iyV[26] = iyp + 1 ;
    }

  for(i=0;i<=8;++i)
    izV[i] = izp ;

  if(connexity == 26)
    {
      for(i=9;i<=17;++i)
	izV[i] = izp - 1 ;
      for(i=18;i<=26;++i)
	izV[i] = izp + 1 ;
    }
}


/*	==================================================================	*/


static void incNeighborhood(int* iV, int n, int connexity)
{
  int i, n_1;

  if (iV[0] == n - 1)
    {
      n_1 = n - 1;
      for(i=0;i<=connexity;++i)
	iV[i] -= n_1;
      return;
    }

  for(i=0;i<=connexity;++i)
    ++iV[i];
}


/*	==================================================================	*/


static void decNeighborhood(int* iV, int n, int connexity)
{
  int i, n_1;

  if (iV[0] == 0)
    {
      n_1 = n - 1;
      for(i=0;i<=connexity;++i)
	iV[i] += n_1;
      return;
    }

  for(i=0;i<=connexity;++i)
    --iV[i];
}



/*	==================================================================	*/
/*		fuzzyAverageMatching: computes an average value in the pattern		*/
/*		matching between two fuzzy objects.											*/
/*							Isabelle Bloch - ENST Paris								*/
/*	==================================================================	*/


int fuzzyAverageMatching(Volume*	volume1,
			 Volume*	volume2,
			 FILE*		resu)
{
  double *info1, *info2, average, cardinal;
  float f_average, f_cardinal;
  int nx, ny, nz, ix, iy, iz, test;
  long oPbL1, oLbS1, oPbL2, oLbS2;

  printftrace (IN, "fuzzyAverageMatching");

  /*	Control tests		*/
  /*	=============		*/

  if ( !testEqualSize (volume1, volume2) ||
	!testType (volume1, DOUBLE) ||
	!testType (volume2, DOUBLE) ||
	!testFilled (volume1) ||
	!testFilled (volume2))
    goto abort;

  /*	Bounds calculation	*/
  /*	==================	*/

  getSize (volume1, &nx, &ny, &nz);
  oPbL1 = offsetPointBetweenLine (volume1);
  oLbS1 = offsetLineBetweenSlice (volume1);
  oPbL2 = offsetPointBetweenLine (volume2);
  oLbS2 = offsetLineBetweenSlice (volume2);

  /*	Tab initialisation	*/
  /*	==================	*/

  info1 = data_DOUBLE (volume1) + offsetFirstPoint (volume1);
  info2 = data_DOUBLE (volume2) + offsetFirstPoint (volume2);


  printf("Beginning of matching\n");

  test = 0;

  /*	Fuzzy pattern matching process	*/
  /*	==============================	*/

  if (EXEC)
    {
      average = 0.;
      cardinal = 0.;
      for(iz=0;iz<nz;iz++)
	{
	  for(iy=0;iy<ny;iy++)
	    {
	      for(ix=0;ix<nx;ix++)
		{
		  average = average + *info1 * *info2 ;
		  cardinal = cardinal + *info1 ;
		  info1++;
		  info2++;
		}
	      info1 += oPbL1;
	      info2 += oPbL2;
	    }
	  info1 += oLbS1;
	  info2 += oLbS2;
	}

      f_average = (float)(average/cardinal);
      f_cardinal = (float)(cardinal);

      printf("cardinal, average : %f %f\n", f_cardinal,f_average);
      fprintf(resu,"%f", f_average);
      test = 1;

      /*			fprintf(resu,"average : %f, cardinal : %f\n", f_average,f_cardinal);
       */
      return(test);
    }

  printftrace (OUT, "");


 abort :
  printfexit("fuzzyAverageMatching\n");
  /* --geoffroy
     return type is not a pointer : I replace NULL by 0
  */
  return (0);
}



/*	==================================================================	*/
/*		fuzzyPatternMatching: computes necessity and possibility degrees	*/
/*		in the pattern matching between two fuzzy objects.						*/
/*							Isabelle Bloch - ENST Paris								*/
/*	==================================================================	*/


int
fuzzyPatternMatching(Volume*	volume1,
		     Volume*	volume2,
		     FILE*	resu)
{
  int				test;
  double			*info1, *info2, possibility, necessity;
  float				f_possibility, f_necessity;
  int				nx, ny, nz, ix, iy, iz;
  long				oPbL1, oLbS1, oPbL2, oLbS2;

  printftrace(IN, "fuzzyPatternMatching");

  /*	Control tests		*/
  /*	=============		*/

  if (!testEqualSize(volume1, volume2) ||
      !testType(volume1, DOUBLE) ||
      !testType(volume2, DOUBLE) ||
      !testFilled(volume1) ||
      !testFilled(volume2))
    goto abort;

  /*	Bounds calculation	*/
  /*	==================	*/

  getSize (volume1, &nx, &ny, &nz);
  oPbL1 = offsetPointBetweenLine (volume1);
  oLbS1 = offsetLineBetweenSlice (volume1);
  oPbL2 = offsetPointBetweenLine (volume2);
  oLbS2 = offsetLineBetweenSlice (volume2);

  /*	Tab initialisation	*/
  /*	==================	*/

  info1 = data_DOUBLE (volume1) + offsetFirstPoint (volume1);
  info2 = data_DOUBLE (volume2) + offsetFirstPoint (volume2);

  test = 0;

  printf("Beginning of matching\n");

  /*	Fuzzy pattern matching process	*/
  /*	==============================	*/

  if (EXEC)
    {
      possibility = 0.;
      necessity = 1.;
      for(iz=0;iz<nz;iz++)
	{
	  for(iy=0;iy<ny;iy++)
	    {
	      for(ix=0;ix<nx;ix++)
		{
		  possibility =
		    ffmax (possibility, ffmin(*info1, *info2));
		  /*
		    possibility =
		    ffmax (possibility,  *info1 * *info2);
		  */
		  necessity =
		    ffmin (necessity, ffmax((1. - *info1), *info2));
		  info1++;
		  info2++;
		}
	      info1 += oPbL1;
	      info2 += oPbL2;
	    }
	  info1 += oLbS1;
	  info2 += oLbS2;
	}

      f_possibility = (float)(possibility);
      f_necessity = (float)(necessity);

      printf("%f %f\n", f_possibility,f_necessity);

      /*			fprintf(resu,"necessity : %f, possibility : %f\n",
	f_necessity,f_possibility);
      */
      fprintf(resu,"%f %f ", f_necessity,f_possibility);
      test = 1;
      return(test);
    }

  printftrace (OUT, "");


 abort :
  printfexit ("fuzzyPatternMatching\n");
  /* --geoffroy
     return type is not a pointer : I replace NULL by 0
  */
  return (0);
}
