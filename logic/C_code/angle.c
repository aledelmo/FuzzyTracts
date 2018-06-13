#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <volume.h>
#include <angle.h>
#include <relativePosition.h>

#define ffmin(a,b) ((a)<(b)?(a):(b))
#define ffmax(a,b) ((a)>(b)?(a):(b))

static void initNeighborhood(int* , int* , int* , int , int , int , int );
static void incNeighborhood(int* , int , int );
static void decNeighborhood(int* , int , int );

Volume*
angle(Volume*		volume1,
      float		alpha1,
      float		alpha2)
{
  Volume		*resuVolume, *tabAngle, *resuX, *resuY, *resuZ;
  U8BIT_t		***info1;
  FLOAT_t		***dist;
  S16BIT_t		***origX, ***origY, ***origZ;
  FLOAT_t		***angle;
  int			nx, ny, nz, ix, iy, iz, ixp, iyp, izp, i, connexity;
  /*	int n, inv; */
  int			nx2, ny2, nz2;
  int			*ixV, *iyV, *izV, X, Y, Z;
  /*	float cos_angle; */
  float			cos1, sin1, cos2, sin2, angle_min;
  float			diffx, diffy, diffz, tmp, dirX, dirY, dirZ ;
  /* float pi = acos(-1.) ; */

  printftrace( IN, "angleVisiblePropag2");
  printf("Hi\n");

  /*	Control tests		*/
  /*	=============		*/

  if (!testType(volume1, U8BIT) ||
	!testFilled(volume1))
    goto abort;

  /*	Bounds calculation	*/
  /*	==================	*/

  getSize (volume1, &nx, &ny, &nz);
  //if (nz == 1 && alpha2 != 0.)
  if (nz == 1 && fegal(alpha2, 0.) == 0)
    fprintf (stderr,"warning: angle not compatible with image dimension\n");

  if (nz == 1)
    connexity = 8 ;
  else
    connexity = 26 ;

  /*	Tab initialisation	*/
  /*	==================	*/
  fprintf(stderr,"init\n");
  fprintf(stderr,"%d %d %d \n", nx,ny,nz);

  setBorderWidth (volume1, 1);
  info1 = tab_U8BIT (volume1);
  setBorderLevel (volume1, 0);
  /*	resuVolume = duplicateVolumeStructure (volume1, "resu");*/
  resuVolume = declareVolumeStructure (nx, ny, nz, 1., 1., 1., FLOAT, "resu", 0);
  /*setType (resuVolume, FLOAT);*/
  setBorderWidth (resuVolume, 0);
  dist = tab_FLOAT (resuVolume);

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
  printf("Hi\n");

  if (EXEC)
    {
      ixV = (int*)malloc((connexity+1)*sizeof(int)) ;
      iyV = (int*)malloc((connexity+1)*sizeof(int)) ;
      izV = (int*)malloc((connexity+1)*sizeof(int)) ;

      printf("Tabulation of angles\n");

      /*	Tabulation of angles	*/
      /*	====================	*/

      cos1 = cosf(alpha1);
      sin1 = sinf(alpha1);
      cos2 = cosf(alpha2);
      sin2 = sinf(alpha2);
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
	      tmp = sqrtf(diffx * diffx + diffy * diffy + diffz * diffz) ;
	      if (fegal(tmp, 0.))
		{
		  angle[iz][iy][ix] = 0 ;
		}
	      else
		{
		  angle[iz][iy][ix] = acosf((dirX * diffx + dirY * diffy + dirZ * diffz) / tmp);
		}
	    }

      angle[nz][ny][nx] = 0 ;

      writeVolume(tabAngle,"angle");

      /*	Initialization of origin points	*/
      /*	===============================	*/

      for(iz=-1;iz<=nz;++iz)
	for(iy=-1;iy<=ny;++iy)
	  for(ix=-1;ix<=nx;++ix)
	    {
	      if (info1[iz][iy][ix] != 0)
		{
		  origX[iz][iy][ix] = (S16BIT_t)ix ;
		  origY[iz][iy][ix] = (S16BIT_t)iy ;
		  origZ[iz][iy][ix] = (S16BIT_t)iz ;
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
		  angle_min = 10000. ;
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
			  if (info1[Z][Y][X]!=0)
			    tmp =  angle[iz][iy][ix];
			  else
			    tmp = 10000.;
			  if (tmp <= angle_min)
			    {
			      angle_min = tmp ;
			      origX[izp][iyp][ixp] = (S16BIT_t)X ;
			      origY[izp][iyp][ixp] = (S16BIT_t)Y ;
			      origZ[izp][iyp][ixp] = (S16BIT_t)Z ;
			    }
			}
		    }
		  dist[izp][iyp][ixp] =  angle_min ;
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
		  angle_min = (float)(dist[izp][iyp][ixp]) ;
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
			  if (info1[Z][Y][X]!=0)
			    tmp =  angle[iz][iy][ix];
			  else
			    tmp = 10000;
			  if (tmp <= angle_min)
			    {
			      angle_min = tmp ;
			      origX[izp][iyp][ixp] = (S16BIT_t)X ;
			      origY[izp][iyp][ixp] = (S16BIT_t)Y ;
			      origZ[izp][iyp][ixp] = (S16BIT_t)Z ;
			    }
			}
		    }
		  dist[izp][iyp][ixp] =  angle_min  ;
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

  printftrace(OUT, "");
  return (resuVolume);


 abort :
  printfexit ("angle\n");
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
