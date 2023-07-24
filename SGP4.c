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

#include "SGP4.h"
#include "more_math.h"
#include "constants.h"
#include "deep.h"
#include <math.h>


void SGP4_init(struct SGP4_state_st *state, const struct tle_t * tle)
{
	/* This function is used to calculate the position and velocity */
	/* of near-earth (period < 225 minutes) satellites. tsince is   */
	/* time since epoch in minutes, tle is a pointer to a tle_t     */
	/* structure with Keplerian orbital elements and pos and vel    */
	/* are vector_t structures returning ECI satellite position and */
	/* velocity. Use Convert_Sat_State() to convert to km and km/s. */
	
	/* Initialization */
	
	/* Recover original mean motion (xnodp) and   */
	/* semimajor axis (aodp) from input elements. */
	
	const double a1=pow(xke/tle->xno,tothrd);
	state->cosio=cos(tle->xincl);
	const double theta2=state->cosio*state->cosio;
	state->x3thm1=3*theta2-1.0;
	const double eosq=tle->eo*tle->eo;
	const double betao2=1.0-eosq;
	const double betao=sqrt(betao2);
	const double del1=1.5*ck2*state->x3thm1/(a1*a1*betao*betao2);
	const double ao=a1*(1.0-del1*(0.5*tothrd+del1*(1.0+134.0/81.0*del1)));
	const double delo=1.5*ck2*state->x3thm1/(ao*ao*betao*betao2);
	state->xnodp=tle->xno/(1.0+delo);
	state->aodp=ao/(1.0-delo);
	
	/* For perigee less than 220 kilometers, the "simple"     */
	/* flag is set and the equations are truncated to linear  */
	/* variation in sqrt a and quadratic variation in mean    */
	/* anomaly.  Also, the c3 term, the delta omega term, and */
	/* the delta m term are dropped.                          */
	
	state->simple = (state->aodp*(1-tle->eo)/ae)<(220/xkmper+ae);
	
	/* For perigees below 156 km, the      */
	/* values of s and qoms2t are altered. */
	
	double s4=s;
	double qoms24=qoms2t;
	const double perigee=(state->aodp*(1-tle->eo)-ae)*xkmper;
	
	if (perigee<156.0)
	{
		if (perigee<=98.0)
			s4=20;
		else
			s4=perigee-78.0;
		
		qoms24=pow((120-s4)*ae/xkmper,4);
		s4=s4/xkmper+ae;
	}
	
	const double pinvsq=1/(state->aodp*state->aodp*betao2*betao2);
	const double tsi=1/(state->aodp-s4);
	state->eta=state->aodp*tle->eo*tsi;
	const double etasq=state->eta*state->eta;
	const double eeta=tle->eo*state->eta;
	const double psisq=fabs(1-etasq);
	const double coef=qoms24*pow(tsi,4);
	const double coef1=coef/pow(psisq,3.5);
	const double c2=coef1*state->xnodp*(state->aodp*(1+1.5*etasq+eeta*(4+etasq))+0.75*ck2*tsi/psisq*state->x3thm1*(8+3*etasq*(8+etasq)));
	state->c1=tle->bstar*c2;
	state->sinio=sin(tle->xincl);
	const double a3ovk2=-xj3/ck2*pow(ae,3);
	const double c3=coef*tsi*a3ovk2*state->xnodp*ae*state->sinio/tle->eo;
	state->x1mth2=1-theta2;
	
	state->c4=2*state->xnodp*coef1*state->aodp*betao2*(state->eta*(2+0.5*etasq)+tle->eo*(0.5+2*etasq)-2*ck2*tsi/(state->aodp*psisq)*(-3*state->x3thm1*(1-2*eeta+etasq*(1.5-0.5*eeta))+0.75*state->x1mth2*(2*etasq-eeta*(1+etasq))*cos(2*tle->omegao)));
	state->c5=2*coef1*state->aodp*betao2*(1+2.75*(etasq+eeta)+eeta*etasq);
	
	const double theta4=theta2*theta2;
	const double temp1=3*ck2*pinvsq*state->xnodp;
	const double temp2=temp1*ck2*pinvsq;
	const double temp3=1.25*ck4*pinvsq*pinvsq*state->xnodp;
	state->xmdot=state->xnodp+0.5*temp1*betao*state->x3thm1+0.0625*temp2*betao*(13-78*theta2+137*theta4);
	const double x1m5th=1-5*theta2;
	state->omgdot=-0.5*temp1*x1m5th+0.0625*temp2*(7-114*theta2+395*theta4)+temp3*(3-36*theta2+49*theta4);
	const double xhdot1=-temp1*state->cosio;
	state->xnodot=xhdot1+(0.5*temp2*(4-19*theta2)+2*temp3*(3-7*theta2))*state->cosio;
	state->omgcof=tle->bstar*c3*cos(tle->omegao);
	state->xmcof=-tothrd*coef*tle->bstar*ae/eeta;
	state->xnodcf=3.5*betao2*xhdot1*state->c1;
	state->t2cof=1.5*state->c1;
	state->xlcof=0.125*a3ovk2*state->sinio*(3+5*state->cosio)/(1+state->cosio);
	state->aycof=0.25*a3ovk2*state->sinio;
	state->delmo=pow(1+state->eta*cos(tle->xmo),3);
	state->sinmo=sin(tle->xmo);
	state->x7thm1=7*theta2-1;
	
	if (!state->simple)
	{
		const double c1sq=state->c1*state->c1;
		state->d2=4*state->aodp*tsi*c1sq;
		const double temp=state->d2*tsi*state->c1/3;
		state->d3=(17*state->aodp+s4)*temp;
		state->d4=0.5*temp*state->aodp*tsi*(221*state->aodp+31*s4)*state->c1;
		state->t3cof=state->d2+2*c1sq;
		state->t4cof=0.25*(3*state->d3+state->c1*(12*state->d2+10*c1sq));
		state->t5cof=0.2*(3*state->d4+12*state->c1*state->d3+6*state->d2*state->d2+15*c1sq*(2*state->d2+c1sq));
	}
}

void SGP4(struct SGP4_state_st *state, double tsince, const struct tle_t * tle, struct vector_t * pos, struct vector_t * vel)
{
	/* This function is used to calculate the position and velocity */
	/* of near-earth (period < 225 minutes) satellites. tsince is   */
	/* time since epoch in minutes, tle is a pointer to a tle_t     */
	/* structure with Keplerian orbital elements and pos and vel    */
	/* are vector_t structures returning ECI satellite position and */ 
	/* velocity. Use Convert_Sat_State() to convert to km and km/s. */

	const double xmdf=tle->xmo+state->xmdot*tsince;
	const double omgadf=tle->omegao+state->omgdot*tsince;
	const double tsq=tsince*tsince;
	const double xnoddf=tle->xnodeo+state->xnodot*tsince;
	const double xnode=xnoddf+state->xnodcf*tsq;
	double xmp=xmdf;;
	double omega=omgadf;
	double tempa=1-state->c1*tsince;
	double tempe=tle->bstar*state->c4*tsince;
	double templ=state->t2cof*tsq;

	/* Update for secular gravity and atmospheric drag. */
	if (!state->simple)
	{
		const double delomg=state->omgcof*tsince;
		const double delm=state->xmcof*(pow(1+state->eta*cos(xmdf),3)-state->delmo);
		const double temp=delomg+delm;
		xmp=xmdf+temp;
		omega=omgadf-temp;
		const double tcube=tsq*tsince;
		const double tfour=tsince*tcube;
		tempa += -state->d2*tsq-state->d3*tcube-state->d4*tfour;
		tempe += tle->bstar*state->c5*(sin(xmp)-state->sinmo);
		templ += state->t3cof*tcube+tfour*(state->t4cof+tsince*state->t5cof);
	}

	const double a=state->aodp*pow(tempa,2);
	const double e=tle->eo-tempe;
	const double xl=xmp+omega+xnode+state->xnodp*templ;
	const double beta=sqrt(1-e*e);
	const double xn=xke/pow(a,1.5);

	/* Long period periodics */
	const double axn=e*cos(omega);
	double temp=1/(a*beta*beta);
	const double xll=temp*state->xlcof*axn;
	const double aynl=temp*state->aycof;
	const double xlt=xl+xll;
	const double ayn=e*sin(omega)+aynl;

	/* Solve Kepler's Equation */
	const double capu=FMod2p(xlt-xnode);
	double temp2=capu;
	double sinepw,cosepw;
	double temp3,temp4,temp5,temp6;
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

		if (fabs(epw-temp2)<= e6a)
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
	const double xnodek=xnode+1.5*temp2*state->cosio*sin2u;
	const double xinck=tle->xincl+1.5*temp2*state->cosio*state->sinio*cos2u;
	const double rdotk=rdot-xn*temp1*state->x1mth2*sin2u;
	const double rfdotk=rfdot+xn*temp1*(state->x1mth2*cos2u+1.5*state->x3thm1);

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

	/* Phase in radians */
	extern double phase;	// TODO: Pass phase via the interface
	phase=xlt-xnode-omgadf+twopi;

	if (phase<0.0)
		phase+=twopi;

	phase=FMod2p(phase);
}
