/** @file relativePosition.h
 * @brief Fuzzy relative position between objects
 */



#include <volume.h>
#ifdef __cplusplus
extern "C" {
#endif


/*	==================================================================	*/
/*   angleVisible.c: computation of the "fuzzy landscape" representing	*/
/*		the degrees of satisfaction of the relative position 					*/
/*    in some direction.																*/
/*		Method using tabulation of angles											*/
/*		et search over the only points belonging to the support.				*/
/*					Isabelle Bloch - ENST Paris										*/
/*     Modified by Olivier Colliot, Oct. 2002
       Add parameter k : the biggest angle such that mu!=0
       mu=max(0,1-beta/k) instead of mu=max(0,1-2beta/pi) (i.e. k=pi/2) */
/*	==================================================================	*/

extern Volume *angleVisible
(
	Volume	   *volume1,
	double		alpha1,
	double		alpha2,
	double k
);


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

extern Volume *angleVisiblePropag
(
	Volume	   *volume1,
	double		alpha1,
	double		alpha2,
	double k
);


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

extern Volume *angleVisiblePropag2
(
	Volume	   *volume1,
	float		alpha1,
	float		alpha2,
	float k
);


/*	==================================================================	*/
/*		fuzzyAverageMatching																*/
/*	==================================================================	*/

extern int fuzzyAverageMatching
(
	Volume	   *volume1,
	Volume		*volume2,
	FILE			*resu
);


/*	==================================================================	*/
/*		fuzzyPatternMatching																*/
/*	==================================================================	*/

extern int fuzzyPatternMatching
(
	Volume	   *volume1,
	Volume		*volume2,
	FILE			*resu
);
#ifdef __cplusplus
}
#endif
