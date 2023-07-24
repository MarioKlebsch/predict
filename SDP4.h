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

#ifndef SDP4_H_
#define SDP4_H_

#include "deep.h"

struct vector_t;
struct tle_t;


struct SDP4_state_st
{
	double x3thm1;
	double c1;
	double x1mth2;
	double c4;
	double xnodcf;
	double t2cof;
	double xlcof;
	double aycof;
	double x7thm1;

	struct deep_state_st deep_state;
	struct deep_arg_t    deep_arg;
};

void SDP4_init(struct SDP4_state_st *state, const tle_t * tle);
void SDP4(struct SDP4_state_st *state, double tsince, const tle_t * tle, struct vector_t * pos, struct vector_t * vel);

#endif /* SDP4_H_ */
