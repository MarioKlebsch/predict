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
#include "deep.h"

#include <math.h>
#ifdef OLD
#include <time.h>
#include <sys/time.h>
#include <curses.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>
#include <unistd.h>
#include <fcntl.h>
#include <termios.h>

#include "predict.h"
//const char *version="2.2.5";
//char *predictpath="/usr/local/lib/predict-2.5.5/";
//int soundcard=0;
#endif


/* Constants used by SGP4/SDP4 code */

#define	km2mi		0.621371		/* km to miles */
#define deg2rad		1.745329251994330E-2	/* Degrees to radians */
#define pi		3.14159265358979323846	/* Pi */
#define pio2		1.57079632679489656	/* Pi/2 */
#define x3pio2		4.71238898038468967	/* 3*Pi/2 */
#define twopi		6.28318530717958623	/* 2*Pi  */
#define e6a		1.0E-6
#define tothrd		6.6666666666666666E-1	/* 2/3 */
#define xj2		1.0826158E-3		/* J2 Harmonic (WGS '72) */
#define xj3		-2.53881E-6		/* J3 Harmonic (WGS '72) */   
#define xj4		-1.65597E-6		/* J4 Harmonic (WGS '72) */
#define xke		7.43669161E-2
#define xkmper		6.378137E3		/* WGS 84 Earth radius km */
#define xmnpda		1.44E3			/* Minutes per day */
#define ae		1.0
#define ck2		5.413079E-4
#define ck4		6.209887E-7
#define f		3.35281066474748E-3	/* Flattening factor */
#define ge		3.986008E5 	/* Earth gravitational constant (WGS '72) */
#define s		1.012229
#define qoms2t		1.880279E-09
#define secday		8.6400E4	/* Seconds per day */
#define omega_E		1.00273790934	/* Earth rotations/siderial day */
#define omega_ER	6.3003879	/* Earth rotations, rads/siderial day */
#define zns		1.19459E-5
#define c1ss		2.9864797E-6
#define zes		1.675E-2
#define znl		1.5835218E-4
#define c1l		4.7968065E-7
#define zel		5.490E-2
#define zcosis		9.1744867E-1
#define zsinis		3.9785416E-1
#define zsings		-9.8088458E-1
#define zcosgs		1.945905E-1
#define zcoshs		1
#define zsinhs		0
#define q22		1.7891679E-6
#define q31		2.1460748E-6
#define q33		2.2123015E-7
#define g22		5.7686396
#define g32		9.5240898E-1
#define g44		1.8014998
#define g52		1.0508330
#define g54		4.4108898
#define root22		1.7891679E-6
#define root32		3.7393792E-7
#define root44		7.3636953E-9
#define root52		1.1428639E-7
#define root54		2.1765803E-9
#define thdt		4.3752691E-3
#define rho		1.5696615E-1
#define mfactor		7.292115E-5
#define sr		6.96000E5	/* Solar radius - km (IAU 76) */
#define AU		1.49597870691E8	/* Astronomical unit - km (IAU 76) */

/* Flow control flag definitions */

#define ALL_FLAGS              -1
#define SGP_INITIALIZED_FLAG   0x000001	/* not used */
#define SGP4_INITIALIZED_FLAG  0x000002
#define SDP4_INITIALIZED_FLAG  0x000004
#define SGP8_INITIALIZED_FLAG  0x000008	/* not used */
#define SDP8_INITIALIZED_FLAG  0x000010	/* not used */
#define DEEP_SPACE_EPHEM_FLAG  0x000040
#define LUNAR_TERMS_DONE_FLAG  0x000080
#define NEW_EPHEMERIS_FLAG     0x000100	/* not used */
#define DO_LOOP_FLAG           0x000200
#define RESONANCE_FLAG         0x000400
#define SYNCHRONOUS_FLAG       0x000800
#define EPOCH_RESTART_FLAG     0x001000
#define VISIBLE_FLAG           0x002000
#define SAT_ECLIPSED_FLAG      0x004000

/* Functions for testing and setting/clearing flags used in SGP4/SDP4 code */

static int isFlagSet(int flag)
{
	extern int Flags;
	return (Flags&flag);
}

static int isFlagClear(int flag)
{
	extern int Flags;
	return (~Flags&flag);
}

static void SetFlag(int flag)
{
	extern int Flags;
	Flags|=flag;
}

static void ClearFlag(int flag)
{
	extern int Flags;
	Flags&=~flag;
}

extern double AcTan(double sinx, double cosx);
extern double FMod2p(double x);
extern double Modulus(double arg1, double arg2);
extern double Julian_Date_of_Year(double year);

static double ThetaG(double epoch, deep_arg_t *deep_arg)
{
	/* The function ThetaG calculates the Greenwich Mean Sidereal Time */
	/* for an epoch specified in the format used in the NORAD two-line */
	/* element sets. It has now been adapted for dates beyond the year */
	/* 1999, as described above. The function ThetaG_JD provides the   */
	/* same calculation except that it is based on an input in the     */
	/* form of a Julian Date. */

	/* Reference:  The 1992 Astronomical Almanac, page B6. */


	/* Modification to support Y2K */
	/* Valid 1957 through 2056     */

	double year;
	double day=modf(epoch*1E-3,&year)*1E3;

	if (year<57)
		year+=2000;
	else
		year+=1900;

	const double UT=modf(day,&day);
	const double jd=Julian_Date_of_Year(year)+day;
	const double TU=(jd-2451545.0)/36525;
	double GMST=24110.54841+TU*(8640184.812866+TU*(0.093104-TU*6.2E-6));
	GMST=Modulus(GMST+secday*omega_E*UT,secday);
	deep_arg->ds50=jd-2433281.5+UT;
	return FMod2p(6.3003880987*deep_arg->ds50+1.72944494);
}


void Deep_init(struct deep_state_st *state, const tle_t * tle, deep_arg_t * deep_arg)
{
	state->thgr=ThetaG(tle->epoch,deep_arg);
	const double eq=tle->eo;
	state->xnq=deep_arg->xnodp;
	const double aqnv=1/deep_arg->aodp;
	state->xqncl=tle->xincl;
	const double xmao=tle->xmo;
	const double xpidot=deep_arg->omgdot+deep_arg->xnodot;
	const double sinq=sin(tle->xnodeo);
	const double cosq=cos(tle->xnodeo);
	state->omegaq=tle->omegao;
	
	/* Initialize lunar solar terms */
	const double day=deep_arg->ds50+18261.5;  /* Days since 1900 Jan 0.5 */
	
	if (day!=state->preep)
	{
		state->preep=day;
		const double xnodce=4.5236020-9.2422029E-4*day;
		const double stem=sin(xnodce);
		const double ctem=cos(xnodce);
		state->zcosil=0.91375164-0.03568096*ctem;
		state->zsinil=sqrt(1-state->zcosil*state->zcosil);
		state->zsinhl=0.089683511*stem/state->zsinil;
		state->zcoshl=sqrt(1-state->zsinhl*state->zsinhl);
		const double c=4.7199672+0.22997150*day;
		const double gam=5.8351514+0.0019443680*day;
		state->zmol=FMod2p(c-gam);
		double zx=0.39785416*stem/state->zsinil;
		const double zy=state->zcoshl*ctem+0.91744867*state->zsinhl*stem;
		zx=AcTan(zx,zy);
		zx=gam+zx-xnodce;
		state->zcosgl=cos(zx);
		state->zsingl=sin(zx);
		state->zmos=6.2565837+0.017201977*day;
		state->zmos=FMod2p(state->zmos);
	}
	
	/* Do solar terms */
	state->savtsn=1E20;
	double zcosg=zcosgs;
	double zsing=zsings;
	double zcosi=zcosis;
	double zsini=zsinis;
	double zcosh=cosq;
	double zsinh= sinq;
	double cc=c1ss;
	double zn=zns;
	double ze=zes;
	const double xnoi=1/state->xnq;
	
	
	/* Loop breaks when Solar terms are done a second */
	/* time, after Lunar terms are initialized        */
	double se,si,sl,sgh,sh;
	for (;;)
	{
		/* Solar terms done again after Lunar terms are done */
		const double a1=zcosg*zcosh+zsing*zcosi*zsinh;
		const double a3=-zsing*zcosh+zcosg*zcosi*zsinh;
		const double a7=-zcosg*zsinh+zsing*zcosi*zcosh;
		const double a8=zsing*zsini;
		const double a9=zsing*zsinh+zcosg*zcosi*zcosh;
		const double a10=zcosg*zsini;
		const double a2=deep_arg->cosio*a7+deep_arg->sinio*a8;
		const double a4=deep_arg->cosio*a9+deep_arg->sinio*a10;
		const double a5=-deep_arg->sinio*a7+deep_arg->cosio*a8;
		const double a6=-deep_arg->sinio*a9+deep_arg->cosio*a10;
		const double x1=a1*deep_arg->cosg+a2*deep_arg->sing;
		const double x2=a3*deep_arg->cosg+a4*deep_arg->sing;
		const double x3=-a1*deep_arg->sing+a2*deep_arg->cosg;
		const double x4=-a3*deep_arg->sing+a4*deep_arg->cosg;
		const double x5=a5*deep_arg->sing;
		const double x6=a6*deep_arg->sing;
		const double x7=a5*deep_arg->cosg;
		const double x8=a6*deep_arg->cosg;
		const double z31=12*x1*x1-3*x3*x3;
		const double z32=24*x1*x2-6*x3*x4;
		const double z33=12*x2*x2-3*x4*x4;
		double z1=3*(a1*a1+a2*a2)+z31*deep_arg->eosq;
		double z2=6*(a1*a3+a2*a4)+z32*deep_arg->eosq;
		double z3=3*(a3*a3+a4*a4)+z33*deep_arg->eosq;
		const double z11=-6*a1*a5+deep_arg->eosq*(-24*x1*x7-6*x3*x5);
		const double z12=-6*(a1*a6+a3*a5)+deep_arg->eosq*(-24*(x2*x7+x1*x8)-6*(x3*x6+x4*x5));
		const double z13=-6*a3*a6+deep_arg->eosq*(-24*x2*x8-6*x4*x6);
		const double z21=6*a2*a5+deep_arg->eosq*(24*x1*x5-6*x3*x7);
		const double z22=6*(a4*a5+a2*a6)+deep_arg->eosq*(24*(x2*x5+x1*x6)-6*(x4*x7+x3*x8));
		const double z23=6*a4*a6+deep_arg->eosq*(24*x2*x6-6*x4*x8);
		z1=z1+z1+deep_arg->betao2*z31;
		z2=z2+z2+deep_arg->betao2*z32;
		z3=z3+z3+deep_arg->betao2*z33;
		const double s3=cc*xnoi;
		const double s2=-0.5*s3/deep_arg->betao;
		const double s4=s3*deep_arg->betao;
		const double s1=-15*eq*s4;
		const double s5=x1*x3+x2*x4;
		const double s6=x2*x3+x1*x4;
		const double s7=x2*x4-x1*x3;
		se=s1*zn*s5;
		si=s2*zn*(z11+z13);
		sl=-zn*s3*(z1+z3-14-6*deep_arg->eosq);
		sgh=s4*zn*(z31+z33-6);
		sh=-zn*s2*(z21+z23);
		
		if (state->xqncl<5.2359877E-2)
			sh=0;
		
		state->ee2=2*s1*s6;
		state->e3=2*s1*s7;
		state->xi2=2*s2*z12;
		state->xi3=2*s2*(z13-z11);
		state->xl2=-2*s3*z2;
		state->xl3=-2*s3*(z3-z1);
		state->xl4=-2*s3*(-21-9*deep_arg->eosq)*ze;
		state->xgh2=2*s4*z32;
		state->xgh3=2*s4*(z33-z31);
		state->xgh4=-18*s4*ze;
		state->xh2=-2*s2*z22;
		state->xh3=-2*s2*(z23-z21);
		
		if (isFlagSet(LUNAR_TERMS_DONE_FLAG))
			break;
		
		/* Do lunar terms */
		state->sse=se;
		state->ssi=si;
		state->ssl=sl;
		state->ssh=sh/deep_arg->sinio;
		state->ssg=sgh-deep_arg->cosio*	state->ssh;
		state->se2=state->ee2;
		state->si2=state->xi2;
		state->sl2=state->xl2;
		state->sgh2=state->xgh2;
		state->sh2=state->xh2;
		state->se3=state->e3;
		state->si3=state->xi3;
		state->sl3=state->xl3;
		state->sgh3=state->xgh3;
		state->sh3=state->xh3;
		state->sl4=state->xl4;
		state->sgh4=state->xgh4;
		zcosg=state->zcosgl;
		zsing=state->zsingl;
		zcosi=state->zcosil;
		zsini=state->zsinil;
		zcosh=state->zcoshl*cosq+state->zsinhl*sinq;
		zsinh=sinq*state->zcoshl-cosq*state->zsinhl;
		zn=znl;
		cc=c1l;
		ze=zel;
		SetFlag(LUNAR_TERMS_DONE_FLAG);
	}
	
	state->sse=state->sse+se;
	state->ssi=state->ssi+si;
	state->ssl=state->ssl+sl;
	state->ssg=state->ssg+sgh-deep_arg->cosio/deep_arg->sinio*sh;
	state->ssh=	state->ssh+sh/deep_arg->sinio;
	
	/* Geopotential resonance initialization for 12 hour orbits */
	ClearFlag(RESONANCE_FLAG);
	ClearFlag(SYNCHRONOUS_FLAG);
	
	double bfact;
	if (!((state->xnq<0.0052359877) && (state->xnq>0.0034906585)))
	{
		if ((state->xnq<0.00826) || (state->xnq>0.00924))
			return;
		
		if (eq<0.5)
			return;
		
		SetFlag(RESONANCE_FLAG);
		const double eoc=eq*deep_arg->eosq;
		const double g201=-0.306-(eq-0.64)*0.440;
		
		double g211, g310, g322, g410, g422, g520;
		if (eq<=0.65)
		{
			g211=3.616-13.247*eq+16.290*deep_arg->eosq;
			g310=-19.302+117.390*eq-228.419*deep_arg->eosq+156.591*eoc;
			g322=-18.9068+109.7927*eq-214.6334*deep_arg->eosq+146.5816*eoc;
			g410=-41.122+242.694*eq-471.094*deep_arg->eosq+313.953*eoc;
			g422=-146.407+841.880*eq-1629.014*deep_arg->eosq+1083.435 * eoc;
			g520=-532.114+3017.977*eq-5740*deep_arg->eosq+3708.276*eoc;
		}
		
		else
		{
			g211=-72.099+331.819*eq-508.738*deep_arg->eosq+266.724*eoc;
			g310=-346.844+1582.851*eq-2415.925*deep_arg->eosq+1246.113*eoc;
			g322=-342.585+1554.908*eq-2366.899*deep_arg->eosq+1215.972*eoc;
			g410=-1052.797+4758.686*eq-7193.992*deep_arg->eosq+3651.957*eoc;
			g422=-3581.69+16178.11*eq-24462.77*deep_arg->eosq+12422.52*eoc;
			
			if (eq<=0.715)
				g520=1464.74-4664.75*eq+3763.64*deep_arg->eosq;
			
			else
				g520=-5149.66+29936.92*eq-54087.36*deep_arg->eosq+31324.56*eoc;
		}
		
		double g533, g521, g532;
		if (eq<0.7)
		{
			g533=-919.2277+4988.61*eq-9064.77*deep_arg->eosq+5542.21*eoc;
			g521=-822.71072+4568.6173*eq-8491.4146*deep_arg->eosq+5337.524*eoc;
			g532=-853.666+4690.25*eq-8624.77*deep_arg->eosq+5341.4*eoc;
		}
		
		else
		{
			g533=-37995.78+161616.52*eq-229838.2*deep_arg->eosq+109377.94*eoc;
			g521 =-51752.104+218913.95*eq-309468.16*deep_arg->eosq+146349.42*eoc;
			g532 =-40023.88+170470.89*eq-242699.48*deep_arg->eosq+115605.82*eoc;
		}
		
		const double sini2=deep_arg->sinio*deep_arg->sinio;
		const double f220=0.75*(1+2*deep_arg->cosio+deep_arg->theta2);
		const double f221=1.5*sini2;
		const double f321=1.875*deep_arg->sinio*(1-2*deep_arg->cosio-3*deep_arg->theta2);
		const double f322=-1.875*deep_arg->sinio*(1+2*deep_arg->cosio-3*deep_arg->theta2);
		const double f441=35*sini2*f220;
		const double f442=39.3750*sini2*sini2;
		const double f522=9.84375*deep_arg->sinio*(sini2*(1-2*deep_arg->cosio-5*deep_arg->theta2)+0.33333333*(-2+4*deep_arg->cosio+6*deep_arg->theta2));
		const double f523=deep_arg->sinio*(4.92187512*sini2*(-2-4*deep_arg->cosio+10*deep_arg->theta2)+6.56250012*(1+2*deep_arg->cosio-3*deep_arg->theta2));
		const double f542=29.53125*deep_arg->sinio*(2-8*deep_arg->cosio+deep_arg->theta2*(-12+8*deep_arg->cosio+10*deep_arg->theta2));
		const double f543=29.53125*deep_arg->sinio*(-2-8*deep_arg->cosio+deep_arg->theta2*(12+8*deep_arg->cosio-10*deep_arg->theta2));
		const double xno2=state->xnq*state->xnq;
		const double ainv2=aqnv*aqnv;
		double temp1=3*xno2*ainv2;
		double temp=temp1*root22;
		state->d2201=temp*f220*g201;
		state->d2211=temp*f221*g211;
		temp1=temp1*aqnv;
		temp=temp1*root32;
		state->d3210=temp*f321*g310;
		state->d3222=temp*f322*g322;
		temp1=temp1*aqnv;
		temp=2*temp1*root44;
		state->d4410=temp*f441*g410;
		state->d4422=temp*f442*g422;
		temp1=temp1*aqnv;
		temp=temp1*root52;
		state->d5220=temp*f522*g520;
		state->d5232=temp*f523*g532;
		temp=2*temp1*root54;
		state->d5421=temp*f542*g521;
		state->d5433=temp*f543*g533;
		state->xlamo=xmao+tle->xnodeo+tle->xnodeo-state->thgr-state->thgr;
		bfact=deep_arg->xmdot+deep_arg->xnodot+deep_arg->xnodot-thdt-thdt;
		bfact=bfact+state->ssl+	state->ssh+	state->ssh;
	}
	
	else
	{
		SetFlag(RESONANCE_FLAG);
		SetFlag(SYNCHRONOUS_FLAG);
		
		/* Synchronous resonance terms initialization */
		const double g200=1+deep_arg->eosq*(-2.5+0.8125*deep_arg->eosq);
		const double g310=1+2*deep_arg->eosq;
		const double g300=1+deep_arg->eosq*(-6+6.60937*deep_arg->eosq);
		const double f220=0.75*(1+deep_arg->cosio)*(1+deep_arg->cosio);
		const double f311=0.9375*deep_arg->sinio*deep_arg->sinio*(1+3*deep_arg->cosio)-0.75*(1+deep_arg->cosio);
		double f330=1+deep_arg->cosio;
		f330=1.875*f330*f330*f330;
		state->del1=3*state->xnq*state->xnq*aqnv*aqnv;
		state->del2=2*state->del1*f220*g200*q22;
		state->del3=3*state->del1*f330*g300*q33*aqnv;
		state->del1=state->del1*f311*g310*q31*aqnv;
		state->fasx2=0.13130908;
		state->fasx4=2.8843198;
		state->fasx6=0.37448087;
		state->xlamo=xmao+tle->xnodeo+tle->omegao-state->thgr;
		bfact=deep_arg->xmdot+xpidot-thdt;
		bfact=bfact+state->ssl+state->ssg+	state->ssh;
	}
	
	state->xfact=bfact-state->xnq;
	
	/* Initialize integrator */
	state->xli=state->xlamo;
	state->xni=state->xnq;
	state->atime=0;
	state->stepp=720;
	state->stepn=-720;
	state->step2=259200;
	
	return;
}


void Deep_sec(struct deep_state_st *state, const tle_t * tle, deep_arg_t * deep_arg) /* Entrance for deep space secular effects */
{
	deep_arg->xll=deep_arg->xll+state->ssl*deep_arg->t;
	deep_arg->omgadf=deep_arg->omgadf+state->ssg*deep_arg->t;
	deep_arg->xnode=deep_arg->xnode+	state->ssh*deep_arg->t;
	deep_arg->em=tle->eo+state->sse*deep_arg->t;
	deep_arg->xinc=tle->xincl+state->ssi*deep_arg->t;
	
	if (deep_arg->xinc<0)
	{
		deep_arg->xinc=-deep_arg->xinc;
		deep_arg->xnode=deep_arg->xnode+pi;
		deep_arg->omgadf=deep_arg->omgadf-pi;
	}
	
	if (isFlagClear(RESONANCE_FLAG))
		return;
	
	double delt, ft, xndot, xnddt, xldot;
	do
	{
		if ((state->atime==0) || ((deep_arg->t>=0) && (state->atime<0)) || ((deep_arg->t<0) && (state->atime>=0)))
		{
			/* Epoch restart */
			if (deep_arg->t>=0)
				delt=state->stepp;
			else
				delt=state->stepn;
			
			state->atime=0;
			state->xni=state->xnq;
			state->xli=state->xlamo;
		}
		else
		{
			if (fabs(deep_arg->t)>=fabs(state->atime))
			{
				if (deep_arg->t>0)
					delt=state->stepp;
				else
					delt=state->stepn;
			}
		}
		
		do
		{
			if (fabs(deep_arg->t-state->atime)>=state->stepp)
			{
				SetFlag(DO_LOOP_FLAG);
				ClearFlag(EPOCH_RESTART_FLAG);
			}
			else
			{
				ft=deep_arg->t-state->atime;
				ClearFlag(DO_LOOP_FLAG);
			}
			
			if (fabs(deep_arg->t)<fabs(state->atime))
			{
				if (deep_arg->t>=0)
					delt=state->stepn;
				else
					delt=state->stepp;
				
				SetFlag(DO_LOOP_FLAG | EPOCH_RESTART_FLAG);
			}
			
			/* Dot terms calculated */
			if (isFlagSet(SYNCHRONOUS_FLAG))
			{
				xndot=state->del1*sin(state->xli-state->fasx2)+state->del2*sin(2*(state->xli-state->fasx4))+state->del3*sin(3*(state->xli-state->fasx6));
				xnddt=state->del1*cos(state->xli-state->fasx2)+2*state->del2*cos(2*(state->xli-state->fasx4))+3*state->del3*cos(3*(state->xli-state->fasx6));
			}
			else
			{
				const double xomi=state->omegaq+deep_arg->omgdot*state->atime;
				const double x2omi=xomi+xomi;
				const double x2li=state->xli+state->xli;
				xndot=state->d2201*sin(x2omi+state->xli-g22)+state->d2211*sin(state->xli-g22)+state->d3210*sin(xomi+state->xli-g32)+state->d3222*sin(-xomi+state->xli-g32)+state->d4410*sin(x2omi+x2li-g44)+state->d4422*sin(x2li-g44)+state->d5220*sin(xomi+state->xli-g52)+state->d5232*sin(-xomi+state->xli-g52)+state->d5421*sin(xomi+x2li-g54)+state->d5433*sin(-xomi+x2li-g54);
				xnddt=state->d2201*cos(x2omi+state->xli-g22)+state->d2211*cos(state->xli-g22)+state->d3210*cos(xomi+state->xli-g32)+state->d3222*cos(-xomi+state->xli-g32)+state->d5220*cos(xomi+state->xli-g52)+state->d5232*cos(-xomi+state->xli-g52)+2*(state->d4410*cos(x2omi+x2li-g44)+state->d4422*cos(x2li-g44)+state->d5421*cos(xomi+x2li-g54)+state->d5433*cos(-xomi+x2li-g54));
			}
			
			xldot=state->xni+state->xfact;
			xnddt=xnddt*xldot;
			
			if (isFlagSet(DO_LOOP_FLAG))
			{
				state->xli=state->xli+xldot*delt+xndot*state->step2;
				state->xni=state->xni+xndot*delt+xnddt*state->step2;
				state->atime=state->atime+delt;
			}
		} while (isFlagSet(DO_LOOP_FLAG) && isFlagClear(EPOCH_RESTART_FLAG));
	} while (isFlagSet(DO_LOOP_FLAG) && isFlagSet(EPOCH_RESTART_FLAG));
	
	deep_arg->xn=state->xni+xndot*ft+xnddt*ft*ft*0.5;
	const double xl=state->xli+xldot*ft+xndot*ft*ft*0.5;
	double temp=-deep_arg->xnode+state->thgr+deep_arg->t*thdt;
	
	if (isFlagClear(SYNCHRONOUS_FLAG))
		deep_arg->xll=xl+temp+temp;
	else
		deep_arg->xll=xl-deep_arg->omgadf+temp;
	
	return;
}


void Deep_per(struct deep_state_st *state, const tle_t * tle, deep_arg_t * deep_arg)/* Entrance for lunar-solar periodics */
{
	const double sinis=sin(deep_arg->xinc);
	const double cosis=cos(deep_arg->xinc);
	
	if (fabs(state->savtsn-deep_arg->t)>=30)
	{
		state->savtsn=deep_arg->t;
		double zm=state->zmos+zns*deep_arg->t;
		double zf=zm+2*zes*sin(zm);
		double sinzf=sin(zf);
		double f2=0.5*sinzf*sinzf-0.25;
		double f3=-0.5*sinzf*cos(zf);
		const double ses=state->se2*f2+state->se3*f3;
		const double sis=state->si2*f2+state->si3*f3;
		const double sls=state->sl2*f2+state->sl3*f3+state->sl4*sinzf;
		state->sghs=state->sgh2*f2+state->sgh3*f3+state->sgh4*sinzf;
		state->shs=state->sh2*f2+state->sh3*f3;
		zm=state->zmol+znl*deep_arg->t;
		zf=zm+2*zel*sin(zm);
		sinzf=sin(zf);
		f2=0.5*sinzf*sinzf-0.25;
		f3=-0.5*sinzf*cos(zf);
		const double sel=state->ee2*f2+state->e3*f3;
		const double sil=state->xi2*f2+state->xi3*f3;
		const double sll=state->xl2*f2+state->xl3*f3+state->xl4*sinzf;
		state->sghl=state->xgh2*f2+state->xgh3*f3+state->xgh4*sinzf;
		state->sh1=state->xh2*f2+state->xh3*f3;
		state->pe=ses+sel;
		state->pinc=sis+sil;
		state->pl=sls+sll;
	}
	
	double pgh=state->sghs+state->sghl;
	double ph=state->shs+state->sh1;
	deep_arg->xinc=deep_arg->xinc+state->pinc;
	deep_arg->em=deep_arg->em+state->pe;
	
	if (state->xqncl>=0.2)
	{
		/* Apply periodics directly */
		ph=ph/deep_arg->sinio;
		pgh=pgh-deep_arg->cosio*ph;
		deep_arg->omgadf=deep_arg->omgadf+pgh;
		deep_arg->xnode=deep_arg->xnode+ph;
		deep_arg->xll=deep_arg->xll+state->pl;
	}
	
	else
	{
		/* Apply periodics with Lyddane modification */
		const double sinok=sin(deep_arg->xnode);
		const double cosok=cos(deep_arg->xnode);
		double alfdp=sinis*sinok;
		double betdp=sinis*cosok;
		const double dalf=ph*cosok+state->pinc*cosis*sinok;
		const double dbet=-ph*sinok+state->pinc*cosis*cosok;
		alfdp=alfdp+dalf;
		betdp=betdp+dbet;
		deep_arg->xnode=FMod2p(deep_arg->xnode);
		double xls=deep_arg->xll+deep_arg->omgadf+cosis*deep_arg->xnode;
		const double dls=state->pl+pgh-state->pinc*deep_arg->xnode*sinis;
		xls=xls+dls;
		const double xnoh=deep_arg->xnode;
		deep_arg->xnode=AcTan(alfdp,betdp);
		
		/* This is a patch to Lyddane modification */
		/* suggested by Rob Matson. */
		
		if (fabs(xnoh-deep_arg->xnode)>pi)
		{
			if (deep_arg->xnode<xnoh)
				deep_arg->xnode+=twopi;
			else
				deep_arg->xnode-=twopi;
		}
		
		deep_arg->xll=deep_arg->xll+state->pl;
		deep_arg->omgadf=xls-deep_arg->xll-cos(deep_arg->xinc)*deep_arg->xnode;
	}
	return;
}
