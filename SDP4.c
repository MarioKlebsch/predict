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

#include "SDP4.h"
#include "more_math.h"
#include "constants.h"
#include <math.h>


void SDP4_init(struct SDP4_state_st *state, const tle_t * tle)
{
	/* This function is used to calculate the position and velocity */
	/* of deep-space (period > 225 minutes) satellites. tle is a pointer to a tle_t     */
	/* structure with Keplerian orbital elements. Use Convert_Sat_State() to convert to km and km/s. */
	
	/* Initialization */
	
	/* Recover original mean motion (xnodp) and   */
	/* semimajor axis (aodp) from input elements. */
	
	const double a1=pow(xke/tle->xno,tothrd);
	state->deep_arg.cosio=cos(tle->xincl);
	state->deep_arg.theta2=state->deep_arg.cosio*state->deep_arg.cosio;
	state->x3thm1=3*state->deep_arg.theta2-1;
	state->deep_arg.eosq=tle->eo*tle->eo;
	state->deep_arg.betao2=1-state->deep_arg.eosq;
	state->deep_arg.betao=sqrt(state->deep_arg.betao2);
	const double del1=1.5*ck2*state->x3thm1/(a1*a1*state->deep_arg.betao*state->deep_arg.betao2);
	const double ao=a1*(1-del1*(0.5*tothrd+del1*(1+134/81*del1)));
	const double delo=1.5*ck2*state->x3thm1/(ao*ao*state->deep_arg.betao*state->deep_arg.betao2);
	state->deep_arg.xnodp=tle->xno/(1+delo);
	state->deep_arg.aodp=ao/(1-delo);
	
	/* For perigee below 156 km, the values */
	/* of s and qoms2t are altered.         */
	
	double s4=s;
	double qoms24=qoms2t;
	const double perigee=(state->deep_arg.aodp*(1-tle->eo)-ae)*xkmper;
	
	if (perigee<156.0)
	{
		if (perigee<=98.0)
			s4=20.0;
		else
			s4=perigee-78.0;
		
		qoms24=pow((120-s4)*ae/xkmper,4);
		s4=s4/xkmper+ae;
	}
	
	const double pinvsq=1/(state->deep_arg.aodp*state->deep_arg.aodp*state->deep_arg.betao2*state->deep_arg.betao2);
	state->deep_arg.sing=sin(tle->omegao);
	state->deep_arg.cosg=cos(tle->omegao);
	const double tsi=1/(state->deep_arg.aodp-s4);
	const double eta=state->deep_arg.aodp*tle->eo*tsi;
	const double etasq=eta*eta;
	const double eeta=tle->eo*eta;
	const double psisq=fabs(1-etasq);
	const double coef=qoms24*pow(tsi,4);
	const double coef1=coef/pow(psisq,3.5);
	const double c2=coef1*state->deep_arg.xnodp*(state->deep_arg.aodp*(1+1.5*etasq+eeta*(4+etasq))+0.75*ck2*tsi/psisq*state->x3thm1*(8+3*etasq*(8+etasq)));
	state->c1=tle->bstar*c2;
	state->deep_arg.sinio=sin(tle->xincl);
	const double a3ovk2=-xj3/ck2*pow(ae,3);
	state->x1mth2=1-state->deep_arg.theta2;
	state->c4=2*state->deep_arg.xnodp*coef1*state->deep_arg.aodp*state->deep_arg.betao2*(eta*(2+0.5*etasq)+tle->eo*(0.5+2*etasq)-2*ck2*tsi/(state->deep_arg.aodp*psisq)*(-3*state->x3thm1*(1-2*eeta+etasq*(1.5-0.5*eeta))+0.75*state->x1mth2*(2*etasq-eeta*(1+etasq))*cos(2*tle->omegao)));
	const double theta4=state->deep_arg.theta2*state->deep_arg.theta2;
	const double temp1=3*ck2*pinvsq*state->deep_arg.xnodp;
	const double temp2=temp1*ck2*pinvsq;
	const double temp3=1.25*ck4*pinvsq*pinvsq*state->deep_arg.xnodp;
	state->deep_arg.xmdot=state->deep_arg.xnodp+0.5*temp1*state->deep_arg.betao*state->x3thm1+0.0625*temp2*state->deep_arg.betao*(13-78*state->deep_arg.theta2+137*theta4);
	const double x1m5th=1-5*state->deep_arg.theta2;
	state->deep_arg.omgdot=-0.5*temp1*x1m5th+0.0625*temp2*(7-114*state->deep_arg.theta2+395*theta4)+temp3*(3-36*state->deep_arg.theta2+49*theta4);
	const double xhdot1=-temp1*state->deep_arg.cosio;
	state->deep_arg.xnodot=xhdot1+(0.5*temp2*(4-19*state->deep_arg.theta2)+2*temp3*(3-7*state->deep_arg.theta2))*state->deep_arg.cosio;
	state->xnodcf=3.5*state->deep_arg.betao2*xhdot1*state->c1;
	state->t2cof=1.5*state->c1;
	state->xlcof=0.125*a3ovk2*state->deep_arg.sinio*(3+5*state->deep_arg.cosio)/(1+state->deep_arg.cosio);
	state->aycof=0.25*a3ovk2*state->deep_arg.sinio;
	state->x7thm1=7*state->deep_arg.theta2-1;
	
	/* initialize Deep() */
	Deep_init(&state->deep_state,tle,&state->deep_arg);
}

void SDP4(struct SDP4_state_st *state, double tsince, const tle_t * tle, struct vector_t * pos, struct vector_t * vel)
{
	/* This function is used to calculate the position and velocity */
	/* of deep-space (period > 225 minutes) satellites. tsince is   */
	/* time since epoch in minutes, tle is a pointer to a tle_t     */
	/* structure with Keplerian orbital elements and pos and vel    */
	/* are vector_t structures returning ECI satellite position and */
	/* velocity. Use Convert_Sat_State() to convert to km and km/s. */

	/* Update for secular gravity and atmospheric drag */
	double xmdf=tle->xmo+state->deep_arg.xmdot*tsince;
	state->deep_arg.omgadf=tle->omegao+state->deep_arg.omgdot*tsince;
	const double xnoddf=tle->xnodeo+state->deep_arg.xnodot*tsince;
	const double tsq=tsince*tsince;
	state->deep_arg.xnode=xnoddf+state->xnodcf*tsq;
	const double tempa=1-state->c1*tsince;
	const double tempe=tle->bstar*state->c4*tsince;
	const double templ=state->t2cof*tsq;
	state->deep_arg.xn=state->deep_arg.xnodp;

	/* Update for deep-space secular effects */
	state->deep_arg.xll=xmdf;
	state->deep_arg.t=tsince;
	Deep_sec(&state->deep_state, tle, &state->deep_arg);

	xmdf=state->deep_arg.xll;
	const double a=pow(xke/state->deep_arg.xn,tothrd)*tempa*tempa;
	state->deep_arg.em=state->deep_arg.em-tempe;
	double xmam=xmdf+state->deep_arg.xnodp*templ;

	/* Update for deep-space periodic effects */
	state->deep_arg.xll=xmam;
	Deep_per(&state->deep_state, tle, &state->deep_arg);

	xmam=state->deep_arg.xll;
	const double xl=xmam+state->deep_arg.omgadf+state->deep_arg.xnode;
	const double beta=sqrt(1-state->deep_arg.em*state->deep_arg.em);
	state->deep_arg.xn=xke/pow(a,1.5);

	/* Long period periodics */
	const double axn=state->deep_arg.em*cos(state->deep_arg.omgadf);
	double temp=1/(a*beta*beta);
	const double xll=temp*state->xlcof*axn;
	const double aynl=temp*state->aycof;
	const double xlt=xl+xll;
	const double ayn=state->deep_arg.em*sin(state->deep_arg.omgadf)+aynl;

	/* Solve Kepler's Equation */
	const double capu=FMod2p(xlt-state->deep_arg.xnode);

	double temp2=capu;
	double sinepw, cosepw;
	double temp3, temp4, temp5, temp6;
	int i=0;
	do
	{
		sinepw=sin(temp2);
		cosepw=cos(temp2);
		temp3=axn*sinepw;
		temp4=ayn*cosepw;
		temp5=axn*cosepw;
		temp6=ayn*sinepw;
		const double epw=(capu-temp4+temp3-temp2)/(1-temp5-temp6)+temp2;
	  
		if (fabs(epw-temp2)<=e6a)
			break;

		temp2=epw;
	} while (i++<10);

	/* Short period preliminary quantities */
	const double ecose=temp5+temp6;
	const double esine=temp3-temp4;
	const double elsq=axn*axn+ayn*ayn;
	temp=1-elsq;
	const double pl=a*temp;
	const double r=a*(1-ecose);
	double temp1=1/r;
	const double rdot=xke*sqrt(a)*esine*temp1;
	const double rfdot=xke*sqrt(pl)*temp1;
	temp2=a*temp1;
	const double betal=sqrt(temp);
	temp3=1/(1+betal);
	const double cosu=temp2*(cosepw-axn+ayn*esine*temp3);
	const double sinu=temp2*(sinepw-ayn-axn*esine*temp3);
	const double u=AcTan(sinu,cosu);
	const double sin2u=2*sinu*cosu;
	const double cos2u=2*cosu*cosu-1;
	temp=1/pl;
	temp1=ck2*temp;
	temp2=temp1*temp;

	/* Update for short periodics */
	const double rk=r*(1-1.5*temp2*betal*state->x3thm1)+0.5*temp1*state->x1mth2*cos2u;
	const double uk=u-0.25*temp2*state->x7thm1*sin2u;
	const double xnodek=state->deep_arg.xnode+1.5*temp2*state->deep_arg.cosio*sin2u;
	const double xinck=state->deep_arg.xinc+1.5*temp2*state->deep_arg.cosio*state->deep_arg.sinio*cos2u;
	const double rdotk=rdot-state->deep_arg.xn*temp1*state->x1mth2*sin2u;
	const double rfdotk=rfdot+state->deep_arg.xn*temp1*(state->x1mth2*cos2u+1.5*state->x3thm1);

	/* Orientation vectors */
	const double sinuk=sin(uk);
	const double cosuk=cos(uk);
	const double sinik=sin(xinck);
	const double cosik=cos(xinck);
	const double sinnok=sin(xnodek);
	const double cosnok=cos(xnodek);
	const double xmx=-sinnok*cosik;
	const double xmy=cosnok*cosik;
	const double ux=xmx*sinuk+cosnok*cosuk;
	const double uy=xmy*sinuk+sinnok*cosuk;
	const double uz=sinik*sinuk;
	const double vx=xmx*cosuk-cosnok*sinuk;
	const double vy=xmy*cosuk-sinnok*sinuk;
	const double vz=sinik*cosuk;

	/* Position and velocity */
	pos->x=rk*ux;
	pos->y=rk*uy;
	pos->z=rk*uz;
	vel->x=rdotk*ux+rfdotk*vx;
	vel->y=rdotk*uy+rfdotk*vy;
	vel->z=rdotk*uz+rfdotk*vz;

	// 	TODO: Pass via the interface
	extern int calc_squint;
	extern double alat;
	extern double alon;
	extern double ax;
	extern double ay;
	extern double az;
	extern double phase;

	/* Calculations for squint angle begin here... */
	if (calc_squint)
	{
		const double bx=cos(alat)*cos(alon+state->deep_arg.omgadf);
		const double by=cos(alat)*sin(alon+state->deep_arg.omgadf);
		const double bz=sin(alat);
		const double cx=bx;
		const double cy=by*cos(xinck)-bz*sin(xinck);
		const double cz=by*sin(xinck)+bz*cos(xinck);
		ax=cx*cos(xnodek)-cy*sin(xnodek);
		ay=cx*sin(xnodek)+cy*cos(xnodek);
		az=cz;
	}
	
	/* Phase in radians */
	phase=xlt-state->deep_arg.xnode-state->deep_arg.omgadf+twopi;
    
	if (phase<0.0)
		phase+=twopi;

	phase=FMod2p(phase);
}

