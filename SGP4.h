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

#ifndef SGP4_H_
#define SGP4_H_

struct vector_t;
struct tle_t;

struct SGP4_state_st
{
	int simple;
	double aodp;
	double aycof;
	double c1;
	double c4;
	double c5;
	double cosio;
	double d2;
	double d3;
	double d4;
	double delmo;
	double omgcof;
	double eta;
	double omgdot;
	double sinio;
	double xnodp;
	double sinmo;
	double t2cof;
	double t3cof;
	double t4cof;
	double t5cof;
	double x1mth2;
	double x3thm1;
	double x7thm1;
	double xmcof;
	double xmdot;
	double xnodcf;
	double xnodot;
	double xlcof;
};

void SGP4_init(struct SGP4_state_st *state, const struct tle_t * tle);
void SGP4(struct SGP4_state_st *state, double tsince, const struct tle_t * tle, struct vector_t * pos, struct vector_t * vel);

#endif /* SGP4_H_ */
