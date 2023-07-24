/***************************************************************************\
*          PREDICT: A satellite tracking/orbital prediction program         *
*          Project started 26-May-1991 by John A. Magliacane, KD2BD         *
*                        Last update: 04-May-2018                           *
*****************************************************************************
*                                                                           *
* This program is free software; you can redistribute it and/or modify it   *
* under the terms of the GNU General Public License as published by the     *
* Free Software Foundation; either version 2 of the License or any later    *
* version.                                                                  *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU         *
* General Public License for more details.                                  *
*                                                                           *
*****************************************************************************
*          See the "CREDITS" file for the names of those who have           *
*   generously contributed their time, talent, and effort to this project.  *
\***************************************************************************/

#ifndef DEEP_H_
#define DEEP_H_

/* Two-line-element satellite orbital data
   structure used directly by the SGP4/SDP4 code. */

typedef struct tle_t {
	double  epoch;
	double  xndt2o;
	double  xndd6o;
	double  bstar;
	double  xincl;
	double  xnodeo;
	double  eo;
	double  omegao;
	double  xmo;
	double  xno;
	long    catnr;
	int     elset;
	long    revnum;
	char    sat_name[25];
	char    idesg[9];
} tle_t;


/* Common arguments between deep-space functions used by SGP4/SDP4 code. */

typedef struct deep_arg_t{
	/* Used by dpinit part of Deep() */
	double  eosq, sinio, cosio, betao, aodp, theta2,
	sing, cosg, betao2, xmdot, omgdot, xnodot, xnodp;
	
	/* Used by dpsec and dpper parts of Deep() */
	double  xll, omgadf, xnode, em, xinc, xn, t;
	
	/* Used by thetg and Deep() */
	double  ds50;
} deep_arg_t;

struct deep_state_st
{
	double thgr, xnq, xqncl, omegaq, zmol, zmos, savtsn, ee2, e3,
	xi2, xl2, xl3, xl4, xgh2, xgh3, xgh4, xh2, xh3, sse, ssi, ssg, xi3,
	se2, si2, sl2, sgh2, sh2, se3, si3, sl3, sgh3, sh3, sl4, sgh4, ssl,
	ssh, d3210, d3222, d4410, d4422, d5220, d5232, d5421, d5433, del1,
	del2, del3, fasx2, fasx4, fasx6, xlamo, xfact, xni, atime, stepp,
	stepn, step2, preep, pl, sghs, xli, d2201, d2211, sghl, sh1, pinc,
	pe, shs, zsingl, zcosgl, zsinhl, zcoshl, zsinil, zcosil;
};
void Deep_init(struct deep_state_st *state, const tle_t * tle, deep_arg_t * deep_arg);
void Deep_sec(struct deep_state_st *state, const tle_t * tle, deep_arg_t * deep_arg); /* deep space secular effects */
void Deep_per(struct deep_state_st *state, const tle_t * tle, deep_arg_t * deep_arg); /* lunar-solar periodics */

#endif /* DEEP_H_ */
