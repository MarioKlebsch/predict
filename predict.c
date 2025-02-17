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

#include <math.h>
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

#include "more_math.h"
#include "deep.h"
#include "SGP4.h"
#include "SDP4.h"
#include "predict.h"
#include "constants.h"
//const char *version="2.2.5";
//char *predictpath="/usr/local/lib/predict-2.5.5/";
//int soundcard=0;


struct sat_db_st;

struct sat_st {
	char line1[70];
	char line2[70];
	char name[25];
	// Line 1
	long catnum;         // Satellite catalog number
	char designator[10]; // International Designator:
	//    * 2 digits: last two digits of launch year,
	//    * 3 digits: launch number of the year
	//    * 3 chars:  piece of the launch
	int year;            // Epoch year (last two digits of year)
	double refepoch;     // Epoch (day of the year and fractional portion of the day)
	double drag;         // First derivative of mean motion; the ballistic coefficient
	double nddot6;       // Second derivative of mean motion
	double bstar;        // B*, the drag term, or radiation pressure coefficient
	long setnum;         // Element set number
	// Line 2
	double incl;         // nclination (degrees)
	double raan;         // Right ascension of the ascending node (degrees)
	double eccn;         // Eccentricity
	double argper;       // Argument of perigee (degrees)
	double meanan;       // Mean anomaly (degrees)
	double meanmo;       // Mean motion (revolutions per day)
	long orbitnum;       // Revolution number at epoch (revolutions)

	const struct sat_db_st *db;
}  sat_table[24];

struct qth_st{
	char callsign[17];
	double stnlat;
	double stnlong;
	int stnalt;
} qth;


struct transponder_st {
	char          name[80];
	double        uplink_start;
	double        uplink_end;
	double        downlink_start;
	double        downlink_end;
	unsigned char dayofweek;
	int           phase_start;
	int           phase_end;
};

struct sat_db_st {
	char name[25];
	long catnum;
	char squintflag;
	double alat;
	double alon;
	unsigned char transponders;
	struct transponder_st transponder[10];
}  sat_db[24];

/* Global variables for sharing data among functions... */

double	eclipse_depth=0,
	sat_azi, sat_ele, sat_range, sat_range_rate,
	sat_lat, sat_lon, sat_alt, sat_vel, phase,
	sun_azi, sun_ele, daynum, sat_footprint, aostime,
	lostime, ax, ay, az, rx, ry, rz, squint, alat, alon,
	sun_ra, sun_dec, sun_lat, sun_lon, sun_range, sun_range_rate,
	moon_az, moon_el, moon_dx, moon_ra, moon_dec, moon_gha, moon_dv;

char	qthfile[50], tlefile[50], dbfile[50], temp[80], output[25],
	serial_port[15], resave=0, reload_tle=0, netport[7],
	once_per_second=0, sat_sun_status, findsun,
	database=0, xterm, io_lat='N', io_lon='W';

int calc_squint;

int	antfd, ma256, socket_flag=0;
int Flags=0;

long	rv;

unsigned char val[256];

/* The following variables are used by the socket server.  They
   are updated in the MultiTrack() and SingleTrack() functions. */

char	tracking_mode[30];

struct server_data_st {
	char	visibility;

	float	az;
	float	el;
	float	longitude;
	float	lattitude;
	float	footprint;
	float	range;
	float	altitude;
	float	velocity;
	float	eclipse_depth;
	float	phase;
	float	squint;

	double	doppler;
	double	nextevent;

	long	aos;
	long	orbitnum;
};
struct server_data_st server_data[24];


unsigned short portbase=0;

/** Type definitions **/

/* Geodetic position structure used by SGP4/SDP4 code. */

typedef struct	{
		   double lat, lon, alt, theta;
		}  geodetic_t;


/* Common arguments between deep-space functions used by SGP4/SDP4 code. */

/* Global structure used by SGP4/SDP4 code. */

geodetic_t obs_geodetic;

/* Functions for testing and setting/clearing flags used in SGP4/SDP4 code */

int isFlagSet(int flag)
{
	return (Flags&flag);
}

int isFlagClear(int flag)
{
	return (~Flags&flag);
}

void SetFlag(int flag)
{
	Flags|=flag;
}

void ClearFlag(int flag)
{
	Flags&=~flag;
}

/* Remaining SGP4/SDP4 code follows... */


void Calculate_Solar_Position(double time, vector_t *solar_vector)
{
	/* Calculates solar position vector */

	double mjd, year, T, M, L, e, C, O, Lsa, nu, R, eps;

	mjd=time-2415020.0;
	year=1900+mjd/365.25;
	T=(mjd+Delta_ET(year)/secday)/36525.0;
	M=Radians(Modulus(358.47583+Modulus(35999.04975*T,360.0)-(0.000150+0.0000033*T)*Sqr(T),360.0));
	L=Radians(Modulus(279.69668+Modulus(36000.76892*T,360.0)+0.0003025*Sqr(T),360.0));
	e=0.01675104-(0.0000418+0.000000126*T)*T;
	C=Radians((1.919460-(0.004789+0.000014*T)*T)*sin(M)+(0.020094-0.000100*T)*sin(2*M)+0.000293*sin(3*M));
	O=Radians(Modulus(259.18-1934.142*T,360.0));
	Lsa=Modulus(L+C-Radians(0.00569-0.00479*sin(O)),twopi);
	nu=Modulus(M+C,twopi);
	R=1.0000002*(1.0-Sqr(e))/(1.0+e*cos(nu));
	eps=Radians(23.452294-(0.0130125+(0.00000164-0.000000503*T)*T)*T+0.00256*cos(O));
	R=AU*R;
	solar_vector->x=R*cos(Lsa);
	solar_vector->y=R*sin(Lsa)*cos(eps);
	solar_vector->z=R*sin(Lsa)*sin(eps);
	solar_vector->w=R;
}

int Sat_Eclipsed(const vector_t *pos, vector_t *sol, double *depth)
{
	/* Calculates satellite's eclipse status and depth */

	double sd_sun, sd_earth, delta;
	vector_t Rho, earth;

	/* Determine partial eclipse */

	sd_earth=ArcSin(xkmper/pos->w);
	Vec_Sub(sol,pos,&Rho);
	sd_sun=ArcSin(sr/Rho.w);
	Scalar_Multiply(-1,pos,&earth);
	delta=Angle(sol,&earth);
	*depth=sd_earth-sd_sun-delta;

	if (sd_earth<sd_sun)
		return 0;
	else
		if (*depth>=0)
			return 1;
	else
		return 0;
}

int select_ephemeris(const tle_t *tle)
{
	/* Selects the apropriate ephemeris type to be used */
	/* for predictions according to the data in the TLE */
	/* It also processes values in the tle set so that  */
	/* they are apropriate for the sgp4/sdp4 routines   */

	/* Period > 225 minutes is deep space */
	double dd1=(xke/tle->xno);
	const double dd2=tothrd;
	const double a1=pow(dd1,dd2);
	const double r1=cos(tle->xincl);
	dd1=(1.0-tle->eo*tle->eo);
	const double temp=ck2*1.5f*(r1*r1*3.0-1.0)/pow(dd1,1.5);
	const double del1=temp/(a1*a1);
	const double ao=a1*(1.0-del1*(tothrd*.5+del1*(del1*1.654320987654321+1.0)));
	const double delo=temp/(ao*ao);
	const double xnodp=tle->xno/(delo+1.0);

	/* Select a deep-space/near-earth ephemeris */

	return twopi/xnodp/xmnpda>=0.15625;
}


void Calculate_User_PosVel(double time, geodetic_t *geodetic, vector_t *obs_pos, vector_t *obs_vel)
{
	/* Calculate_User_PosVel() passes the user's geodetic position
	   and the time of interest and returns the ECI position and
	   velocity of the observer.  The velocity calculation assumes
	   the geodetic position is stationary relative to the earth's
	   surface. */

	/* Reference:  The 1992 Astronomical Almanac, page K11. */

	double c, sq, achcp;

	geodetic->theta=FMod2p(ThetaG_JD(time)+geodetic->lon); /* LMST */
	c=1/sqrt(1+f*(f-2)*Sqr(sin(geodetic->lat)));
	sq=Sqr(1-f)*c;
	achcp=(xkmper*c+geodetic->alt)*cos(geodetic->lat);
	obs_pos->x=achcp*cos(geodetic->theta); /* kilometers */
	obs_pos->y=achcp*sin(geodetic->theta);
	obs_pos->z=(xkmper*sq+geodetic->alt)*sin(geodetic->lat);
	obs_vel->x=-mfactor*obs_pos->y; /* kilometers/second */
	obs_vel->y=mfactor*obs_pos->x;
	obs_vel->z=0;
	Magnitude(obs_pos);
	Magnitude(obs_vel);
}

void Calculate_LatLonAlt(double time, vector_t *pos,  geodetic_t *geodetic)
{
	/* Procedure Calculate_LatLonAlt will calculate the geodetic  */
	/* position of an object given its ECI position pos and time. */
	/* It is intended to be used to determine the ground track of */
	/* a satellite.  The calculations  assume the earth to be an  */
	/* oblate spheroid as defined in WGS '72.                     */

	/* Reference:  The 1992 Astronomical Almanac, page K12. */

	double r, e2, phi, c;

	geodetic->theta=AcTan(pos->y,pos->x); /* radians */
	geodetic->lon=FMod2p(geodetic->theta-ThetaG_JD(time)); /* radians */
	r=sqrt(Sqr(pos->x)+Sqr(pos->y));
	e2=f*(2-f);
	geodetic->lat=AcTan(pos->z,r); /* radians */

	do
	{
		phi=geodetic->lat;
		c=1/sqrt(1-e2*Sqr(sin(phi)));
		geodetic->lat=AcTan(pos->z+xkmper*c*e2*sin(phi),r);

	} while (fabs(geodetic->lat-phi)>=1E-10);

	geodetic->alt=r/cos(geodetic->lat)-xkmper*c; /* kilometers */

	if (geodetic->lat>pio2)
		geodetic->lat-=twopi;
}

void Calculate_Obs(double time, const vector_t *pos, const vector_t *vel, geodetic_t *geodetic, vector_t *obs_set)
{
	/* The procedures Calculate_Obs and Calculate_RADec calculate         */
	/* the *topocentric* coordinates of the object with ECI position,     */
	/* {pos}, and velocity, {vel}, from location {geodetic} at {time}.    */
	/* The {obs_set} returned for Calculate_Obs consists of azimuth,      */
	/* elevation, range, and range rate (in that order) with units of     */
	/* radians, radians, kilometers, and kilometers/second, respectively. */
	/* The WGS '72 geoid is used and the effect of atmospheric refraction */
	/* (under standard temperature and pressure) is incorporated into the */
	/* elevation calculation; the effect of atmospheric refraction on     */
	/* range and range rate has not yet been quantified.                  */

	/* The {obs_set} for Calculate_RADec consists of right ascension and  */
	/* declination (in that order) in radians.  Again, calculations are   */
	/* based on *topocentric* position using the WGS '72 geoid and        */
	/* incorporating atmospheric refraction.                              */

	double sin_lat, cos_lat, sin_theta, cos_theta, el, azim, top_s, top_e, top_z;

	vector_t obs_pos, obs_vel, range, rgvel;

	Calculate_User_PosVel(time, geodetic, &obs_pos, &obs_vel);

	range.x=pos->x-obs_pos.x;
	range.y=pos->y-obs_pos.y;
	range.z=pos->z-obs_pos.z;

	/* Save these values globally for calculating squint angles later... */

	rx=range.x;
	ry=range.y;
	rz=range.z;

	rgvel.x=vel->x-obs_vel.x;
	rgvel.y=vel->y-obs_vel.y;
	rgvel.z=vel->z-obs_vel.z;

	Magnitude(&range);

	sin_lat=sin(geodetic->lat);
	cos_lat=cos(geodetic->lat);
	sin_theta=sin(geodetic->theta);
	cos_theta=cos(geodetic->theta);
	top_s=sin_lat*cos_theta*range.x+sin_lat*sin_theta*range.y-cos_lat*range.z;
	top_e=-sin_theta*range.x+cos_theta*range.y;
	top_z=cos_lat*cos_theta*range.x+cos_lat*sin_theta*range.y+sin_lat*range.z;
	azim=atan(-top_e/top_s); /* Azimuth */

	if (top_s>0.0) 
		azim=azim+pi;

	if (azim<0.0)
		azim=azim+twopi;

	el=ArcSin(top_z/range.w);
	obs_set->x=azim;	/* Azimuth (radians)   */
	obs_set->y=el;		/* Elevation (radians) */
	obs_set->z=range.w;	/* Range (kilometers)  */

	/* Range Rate (kilometers/second) */

	obs_set->w=Dot(&range,&rgvel)/range.w;

	/* Corrections for atmospheric refraction */
	/* Reference:  Astronomical Algorithms by Jean Meeus, pp. 101-104    */
	/* Correction is meaningless when apparent elevation is below horizon */

	/*** The following adjustment for
		 atmospheric refraction is bypassed ***/

	/* obs_set->y=obs_set->y+Radians((1.02/tan(Radians(Degrees(el)+10.3/(Degrees(el)+5.11))))/60); */

	obs_set->y=el;

	/**** End bypass ****/

	if (obs_set->y>=0.0)
		SetFlag(VISIBLE_FLAG);
	else
	{
		obs_set->y=el;  /* Reset to true elevation */
		ClearFlag(VISIBLE_FLAG);
	}
}

void Calculate_RADec(double time, const vector_t *pos, const vector_t *vel, geodetic_t *geodetic, vector_t *obs_set)
{
	/* Reference:  Methods of Orbit Determination by  */
	/*             Pedro Ramon Escobal, pp. 401-402   */

	double	phi, theta, sin_theta, cos_theta, sin_phi, cos_phi, az, el,
		Lxh, Lyh, Lzh, Sx, Ex, Zx, Sy, Ey, Zy, Sz, Ez, Zz, Lx, Ly,
		Lz, cos_delta, sin_alpha, cos_alpha;

	Calculate_Obs(time,pos,vel,geodetic,obs_set);

	az=obs_set->x;
	el=obs_set->y;
	phi=geodetic->lat;
	theta=FMod2p(ThetaG_JD(time)+geodetic->lon);
	sin_theta=sin(theta);
	cos_theta=cos(theta);
	sin_phi=sin(phi);
	cos_phi=cos(phi);
	Lxh=-cos(az)*cos(el);
	Lyh=sin(az)*cos(el);
	Lzh=sin(el);
	Sx=sin_phi*cos_theta;
	Ex=-sin_theta;
	Zx=cos_theta*cos_phi;
	Sy=sin_phi*sin_theta;
	Ey=cos_theta;
	Zy=sin_theta*cos_phi;
	Sz=-cos_phi;
	Ez=0.0;
	Zz=sin_phi;
	Lx=Sx*Lxh+Ex*Lyh+Zx*Lzh;
	Ly=Sy*Lxh+Ey*Lyh+Zy*Lzh;
	Lz=Sz*Lxh+Ez*Lyh+Zz*Lzh;
	obs_set->y=ArcSin(Lz);  /* Declination (radians) */
	cos_delta=sqrt(1.0-Sqr(Lz));
	sin_alpha=Ly/cos_delta;
	cos_alpha=Lx/cos_delta;
	obs_set->x=AcTan(sin_alpha,cos_alpha); /* Right Ascension (radians) */
	obs_set->x=FMod2p(obs_set->x);
}

/* .... SGP4/SDP4 functions end .... */

void bailout(const char * string)
{
	/* This function quits ncurses, resets and "beeps"
	   the terminal, and displays an error message (string)
	   when we need to bail out of the program in a hurry. */

	beep();	
	curs_set(1);
	bkgdset(COLOR_PAIR(1));
	clear();
	refresh();
	endwin();
	fprintf(stderr,"*** predict: %s!\n",string);
}

void TrackDataOut(int antfd, double elevation, double azimuth)
{
	/* This function sends Azimuth and Elevation data
	   to an antenna tracker connected to the serial port */

	size_t n;
	int port;
	char message[30]="\n";

	port=antfd;

	sprintf(message, "AZ%3.1f EL%3.1f \x0D\x0A", azimuth,elevation);
	n=write(port,message,strlen(message));

	if (n<0)
	{
		bailout("Error Writing To Antenna Port");
		exit(-1);
	}
}

int passivesock(const char *service, const char *protocol, int qlen)
{
	/* This function opens the socket port */

	struct servent *pse;
	struct protoent *ppe;
	struct sockaddr_in sin;
	int sd, type;
	
	memset((char *)&sin, 0, sizeof(struct sockaddr_in));
	sin.sin_family=AF_INET;
	sin.sin_addr.s_addr=INADDR_ANY;
	
	if ((pse=getservbyname(service,protocol)))
		sin.sin_port=htons(ntohs((unsigned short)pse->s_port)+portbase);

	else if ((sin.sin_port=htons((unsigned short)atoi(service)))==0)
	{
		bailout("Can't get service");
		exit(-1);
	}
	
	if ((ppe=getprotobyname(protocol))==0)
	{
		bailout("Can't get protocol");
		exit(-1);
	}
	
	if (strcmp(protocol,"udp")==0)
		type=SOCK_DGRAM;
	else
		type=SOCK_STREAM;
	
	sd=socket(PF_INET,type, ppe->p_proto);

	if (sd<0)
	{
		bailout("Can't open socket");
		exit(-1);
	}
	
	if (bind(sd,(struct sockaddr *)&sin,sizeof(sin))<0)
	{
		bailout("Can't bind");
		exit(-1);
	}
	
	if ((type==SOCK_STREAM) && listen(s,qlen)<0)
	{
		bailout("Listen fail");
		exit(-1);
	}

	return sd;
}

void socket_server(const char * predict_name)
{
	/* This is the socket server code */

	int i, j, sock;
	size_t n;
	socklen_t alen;
	struct sockaddr_in fsin;
	char buf[80], buff[1000], satname[50], tempname[30], ok;
	time_t t;
	FILE *fd=NULL;

	/* Open a socket port at "predict" or netport if defined */

	if (netport[0]==0)
		strncpy(netport,"predict",7);

	sock=passivesock(netport,"udp",10);
 	alen=sizeof(fsin);
	
	/* This is the main loop for monitoring the socket
	   port and sending back replies to clients */

	while (1)
	{
		/* Get datagram from socket port */
		
		if ((n=recvfrom(sock,buf,sizeof(buf),0,(struct sockaddr *)&fsin,&alen)) < 0)
			exit (-1);

		buf[n]=0;
		ok=0;

		/* Parse the command in the datagram */
		if ((strncmp("GET_SAT",buf,7)==0) && (strncmp("GET_SAT_POS",buf,11)!=0))
		{
			/* Parse "buf" for satellite name */
			for (i=0; buf[i]!=32 && buf[i]!=0 && i<39; i++);

			for (j=++i; buf[j]!='\n' && buf[j]!=0 && (j-i)<25; j++)
				satname[j-i]=buf[j];

			satname[j-i]=0;

			/* Do a simple search for the matching satellite name */

			for (i=0; i<24; i++)
			{
				if ((strncmp(satname,sat_table[i].name,25)==0) || (atol(satname)==sat_table[i].catnum))
				{
					struct server_data_st * const server = &server_data[i];
					long nxtevt=(long)rint(86400.0*(server->nextevent+3651.0));

					/* Build text buffer with satellite data */
					sprintf(buff,"%s\n%-7.2f\n%+-6.2f\n%-7.2f\n%+-6.2f\n%ld\n%-7.2f\n%-7.2f\n%-7.2f\n%-7.2f\n%ld\n%c\n%-7.2f\n%-7.2f\n%-7.2f\n",
							sat_table[i].name,
							server->longitude,     server->lattitude,
							server->az,            server->el,
							nxtevt,                server->footprint,
							server->range,         server->altitude,
							server->velocity,      server->orbitnum,
							server->visibility,    server->phase,
							server->eclipse_depth, server->squint);

					/* Send buffer back to the client that sent the request */
					sendto(sock,buff,strlen(buff),0,(struct sockaddr*)&fsin,sizeof(fsin));
					ok=1;
					break;
				}
			}
		}

		if (strncmp("GET_TLE",buf,7)==0)
		{
			/* Parse "buf" for satellite name */
			for (i=0; buf[i]!=32 && buf[i]!=0 && i<39; i++);

			for (j=++i; buf[j]!='\n' && buf[j]!=0 && (j-i)<25; j++)
				satname[j-i]=buf[j];

			satname[j-i]=0;

			/* Do a simple search for the matching satellite name */

			for (i=0; i<24; i++)
			{
				if ((strncmp(satname,sat_table[i].name,25)==0) || (atol(satname)==sat_table[i].catnum))
				{
					/* Build text buffer with satellite data */

					sprintf(buff,"%s\n%s\n%s\n",sat_table[i].name,sat_table[i].line1, sat_table[i].line2);
					/* Send buffer back to the client that sent the request */
					sendto(sock,buff,strlen(buff),0,(struct sockaddr*)&fsin,sizeof(fsin));
					ok=1;
					break;
				}
			}
		}

		if (strncmp("GET_DOPPLER",buf,11)==0)
		{
			/* Parse "buf" for satellite name */
			for (i=0; buf[i]!=32 && buf[i]!=0 && i<39; i++);

			for (j=++i; buf[j]!='\n' && buf[j]!=0 && (j-i)<25; j++)
				satname[j-i]=buf[j];

			satname[j-i]=0;

			/* Do a simple search for the matching satellite name */

			for (i=0; i<24; i++)
			{
				if ((strncmp(satname,sat_table[i].name,25)==0) || (atol(satname)==sat_table[i].catnum))
				{
					/* Get Normalized (100 MHz)
					   Doppler shift for sat[i] */

					sprintf(buff,"%f\n",server_data[i].doppler);

					/* Send buffer back to client who sent request */
					sendto(sock,buff,strlen(buff),0,(struct sockaddr*)&fsin,sizeof(fsin));
					ok=1;
					break;
				}
			}
		}

		if (strncmp("GET_LIST",buf,8)==0)
		{
			buff[0]=0;

			for (i=0; i<24; i++)
			{
				if (sat_table[i].name[0]!=0)
					strcat(buff,sat_table[i].name);

				strcat(buff,"\n");
			}

			sendto(sock,buff,strlen(buff),0,(struct sockaddr *)&fsin,sizeof(fsin));
			ok=1;
		}

		if (strncmp("RELOAD_TLE",buf,10)==0)
		{
			buff[0]=0;
			sendto(sock,buff,strlen(buff),0,(struct sockaddr *)&fsin,sizeof(fsin));
			reload_tle=1;
			ok=1;
		}

		if (strncmp("GET_SUN",buf,7)==0)
		{
			buff[0]=0;
			sprintf(buff,"%-7.2f\n%+-6.2f\n%-7.2f\n%-7.2f\n%-7.2f\n",sun_azi, sun_ele, sun_lat, sun_lon, sun_ra);
			sendto(sock,buff,strlen(buff),0,(struct sockaddr *)&fsin,sizeof(fsin));
			ok=1;
		}

		if (strncmp("GET_MOON",buf,8)==0)
		{
			buff[0]=0;
			sprintf(buff,"%-7.2f\n%+-6.2f\n%-7.2f\n%-7.2f\n%-7.2f\n",moon_az, moon_el, moon_dec, moon_gha, moon_ra);
			sendto(sock,buff,strlen(buff),0,(struct sockaddr *)&fsin,sizeof(fsin));
			ok=1;
		}

		if (strncmp("GET_MODE",buf,8)==0)
		{
			sendto(sock,tracking_mode,strlen(tracking_mode),0,(struct sockaddr *)&fsin,sizeof(fsin));
			ok=1;
		}

		if (strncmp("GET_VERSION",buf,11)==0)
		{
			buff[0]=0;
			sprintf(buff,"%s\n",version);
			sendto(sock,buff,strlen(buff),0,(struct sockaddr *)&fsin,sizeof(fsin));
			ok=1;
		}

		if (strncmp("GET_QTH",buf,7)==0)
		{
			buff[0]=0;
			sprintf(buff,"%s\n%g\n%g\n%d\n",qth.callsign, qth.stnlat, qth.stnlong, qth.stnalt);
			sendto(sock,buff,strlen(buff),0,(struct sockaddr *)&fsin,sizeof(fsin));
			ok=1;
		}

		if (strncmp("GET_TIME$",buf,9)==0)
		{
			buff[0]=0;
			t=time(NULL);
			sprintf(buff,"%s",asctime(gmtime(&t)));

			if (buff[8]==32)
				buff[8]='0';

			sendto(sock,buff,strlen(buff),0,(struct sockaddr *)&fsin,sizeof(fsin));
			buf[0]=0;
			ok=1;
		}

		if (strncmp("GET_TIME",buf,8)==0)
		{
			buff[0]=0;
			t=time(NULL);
			sprintf(buff,"%lu\n",(unsigned long)t);
			sendto(sock,buff,strlen(buff),0,(struct sockaddr *)&fsin,sizeof(fsin));
			ok=1;
		}

		if (strncmp("GET_SAT_POS",buf,11)==0)
		{
			/* Parse "buf" for satellite name and arguments */
			for (i=0; buf[i]!=32 && buf[i]!=0 && i<39; i++);

			for (j=++i; buf[j]!='\n' && buf[j]!=0 && (j-i)<49; j++)
				satname[j-i]=buf[j];

			satname[j-i]=0;

			/* Send request to predict with output
			   directed to a temporary file under /tmp */

			strcpy(tempname,"/tmp/XXXXXX\0");
			i=mkstemp(tempname);

			sprintf(buff,"%s -f %s -t %s -q %s -o %s\n",predict_name,satname,tlefile,qthfile,tempname);
			system(buff);

			/* Append an EOF marker (CNTRL-Z) to the end of file */

			fd=fopen(tempname,"a");
			fprintf(fd,"%c\n",26);  /* Control-Z */
			fclose(fd);

			buff[0]=0;

			/* Send the file to the client */

			fd=fopen(tempname,"rb");

			fgets(buff,80,fd);

			do
			{
				sendto(sock,buff,strlen(buff),0,(struct sockaddr *)&fsin,sizeof(fsin));
				fgets(buff,80,fd);
				/* usleep(2);  if needed (for flow-control) */

			} while (feof(fd)==0);

			fclose(fd);
			unlink(tempname);
			close(i);
			ok=1;
		}

		if (strncmp("PREDICT",buf,7)==0)
		{
			/* Parse "buf" for satellite name and arguments */
			for (i=0; buf[i]!=32 && buf[i]!=0 && i<39; i++);

			for (j=++i; buf[j]!='\n' && buf[j]!=0 && (j-i)<49; j++)
				satname[j-i]=buf[j];

			satname[j-i]=0;

			/* Send request to predict with output
			   directed to a temporary file under /tmp */

			strcpy(tempname,"/tmp/XXXXXX\0");
			i=mkstemp(tempname);

			sprintf(buff,"%s -p %s -t %s -q %s -o %s\n",predict_name, satname,tlefile,qthfile,tempname);
			system(buff);

			/* Append an EOF marker (CNTRL-Z) to the end of file */

			fd=fopen(tempname,"a");
			fprintf(fd,"%c\n",26);  /* Control-Z */
			fclose(fd);

			buff[0]=0;

			/* Send the file to the client */

			fd=fopen(tempname,"rb");

			fgets(buff,80,fd);

			do
			{
				sendto(sock,buff,strlen(buff),0,(struct sockaddr *)&fsin,sizeof(fsin));
				fgets(buff,80,fd);
				/* usleep(2);  if needed (for flow-control) */

			} while (feof(fd)==0);

			fclose(fd);
			unlink(tempname);
			close(i);
			ok=1;
		}

		if (ok==0)
			sendto(sock,"Huh?\n",5,0,(struct sockaddr *)&fsin,sizeof(fsin));
	} 	
}

void Banner(void)
{
	curs_set(0);
	bkgdset(COLOR_PAIR(3));
	clear();
	refresh();

	attrset(COLOR_PAIR(6)|A_REVERSE|A_BOLD);
	mvprintw(3,18,"                                           ");
	mvprintw(4,18,"         --== PREDICT  v%s ==--         ",version);
	mvprintw(5,18,"   Released by John A. Magliacane, KD2BD   ");
	mvprintw(6,18,"                  May 2018                 ");
	mvprintw(7,18,"                                           ");
}

void AnyKey(void)
{
	mvprintw(23,24,"<< Press Any Key To Continue >>");
	refresh();
	getch();
}

double FixAngle(double x)
{
	/* This function reduces angles greater than
	   two pi by subtracting two pi from the angle */

	while (x>twopi)
		x-=twopi;

	return x;
}

double PrimeAngle(double x)
{
	/* This function is used in the FindMoon() function. */

	x=x-360.0*floor(x/360.0);
	return x;
}

char *SubString(const char *string, unsigned start, unsigned end)
{
	/* This function returns a substring based on the starting
	   and ending positions provided.  It is used heavily in
	   the AutoUpdate function when parsing 2-line element data. */

	static char temp[80];
	unsigned x, y;
	if (end>=start)
	{
		for (x=start, y=0; x<=end && string[x]!=0; x++)
			if (string[x]!=' ')
			{
				temp[y]=string[x];
				y++;
			}

		temp[y]=0;
		return temp;
	}
	else
		return NULL;
}

void CopyString(const char *source, char *destination, unsigned start, unsigned end)
{
	/* This function copies elements of the string "source"
	   bounded by "start" and "end" into the string "destination". */

	unsigned j, k=0;

	for (j=start; j<=end; j++)
		if (source[k]!=0)
		{
			destination[j]=source[k];
			k++;
		}
}

char *Abbreviate(const char * string,int n)
{
	/* This function returns an abbreviated substring of the original,
	   including a '~' character if a non-blank character is chopped
	   out of the generated substring.  n is the length of the desired
	   substring.  It is used for abbreviating satellite names. */
	static char temp[80];

	strncpy(temp,string,79);

	if (temp[n]!=0 && temp[n]!=32)
	{
		temp[n-2]='~';
		temp[n-1]=temp[strlen(temp)-1];
	}

	temp[n]=0;

	return temp;
}

char KepCheck(const char *line1,const char *line2)
{
	/* This function scans line 1 and line 2 of a NASA 2-Line element
	   set and returns a 1 if the element set appears to be valid or
	   a 0 if it does not.  If the data survives this torture test,
	   it's a pretty safe bet we're looking at a valid 2-line
	   element set and not just some random text that might pass
	   as orbital data based on a simple checksum calculation alone. */

	int x;
	unsigned sum1, sum2;

	/* Compute checksum for each line */

	for (x=0, sum1=0, sum2=0; x<=67; sum1+=val[(int)line1[x]], sum2+=val[(int)line2[x]], x++);

	/* Perform a "torture test" on the data */

	x=(val[(int)line1[68]]^(sum1%10)) | (val[(int)line2[68]]^(sum2%10)) |
	  (line1[0]^'1')  | (line1[1]^' ')  | (line1[7]^'U')  |
	  (line1[8]^' ')  | (line1[17]^' ') | (line1[23]^'.') |
	  (line1[32]^' ') | (line1[34]^'.') | (line1[43]^' ') |
	  (line1[52]^' ') | (line1[61]^' ') | (line1[62]^'0') |
	  (line1[63]^' ') | (line2[0]^'2')  | (line2[1]^' ')  |
	  (line2[7]^' ')  | (line2[11]^'.') | (line2[16]^' ') |
	  (line2[20]^'.') | (line2[25]^' ') | (line2[33]^' ') |
	  (line2[37]^'.') | (line2[42]^' ') | (line2[46]^'.') |
	  (line2[51]^' ') | (line2[54]^'.') | (line1[2]^line2[2]) |
	  (line1[3]^line2[3]) | (line1[4]^line2[4]) |
	  (line1[5]^line2[5]) | (line1[6]^line2[6]) |
	  (isdigit(line1[68]) ? 0 : 1) | (isdigit(line2[68]) ? 0 : 1) |
	  (isdigit(line1[18]) ? 0 : 1) | (isdigit(line1[19]) ? 0 : 1) |
	  (isdigit(line2[31]) ? 0 : 1) | (isdigit(line2[32]) ? 0 : 1);

	return (x ? 0 : 1);
}

void InternalUpdate(struct sat_st *sat)
{
	/* Updates data in TLE structure based on
	   line1 and line2 stored in structure. */

	double tempnum;

	sat->catnum=atol(SubString(sat->line1,2,6));
	strncpy(sat->designator,SubString(sat->line1,9,16),8);
	sat->designator[9]=0;
	sat->year=atoi(SubString(sat->line1,18,19));
	sat->refepoch=atof(SubString(sat->line1,20,31));
	sat->drag=atof(SubString(sat->line1,33,42));
	tempnum=1.0e-5*atof(SubString(sat->line1,44,49));
	sat->nddot6=tempnum/pow(10.0,(sat->line1[51]-'0'));
	tempnum=1.0e-5*atof(SubString(sat->line1,53,58));
	sat->bstar=tempnum/pow(10.0,(sat->line1[60]-'0'));
	sat->setnum=atol(SubString(sat->line1,64,67));
	sat->incl=atof(SubString(sat->line2,8,15));
	sat->raan=atof(SubString(sat->line2,17,24));
	sat->eccn=1.0e-07*atof(SubString(sat->line2,26,32));
	sat->argper=atof(SubString(sat->line2,34,41));
	sat->meanan=atof(SubString(sat->line2,43,50));
	sat->meanmo=atof(SubString(sat->line2,52,62));
	sat->orbitnum=atof(SubString(sat->line2,63,67));
}

const char *noradEvalue(double value)
{
	/* Converts numeric values to E notation used in NORAD TLEs */
	static char output[25];
	char string[15];

	sprintf(string,"%11.4e",value*10.0);

	output[0]=string[0];
	output[1]=string[1];
	output[2]=string[3];
	output[3]=string[4];
	output[4]=string[5];
	output[5]=string[6];
	output[6]='-';
	output[7]=string[10];
	output[8]=0;

	return output;
}

void Data2TLE(const struct sat_st *sat, char *line1_out, char *line2_out)
{
	/* This function converts orbital data held in the numeric
	   portion of the sat tle structure to ASCII TLE format,
	   and places the result in ASCII portion of the structure. */
 
	int i;
	char string[15];
	unsigned sum;

	/* Fill lines with blanks */

	for (i=0; i<70; line1_out[i]=' ', line2_out[i]=' ', i++);

	line1_out[69]=0;
	line2_out[69]=0;

	/* Insert static characters */

	line1_out[0]='1';
	line1_out[7]='U'; /* Unclassified */
	line2_out[0]='2';

	line1_out[62]='0'; /* For publically released TLEs */

	/* Insert orbital data */

	sprintf(string,"%05ld",sat->catnum);
	CopyString(string,line1_out,2,6);
	CopyString(string,line2_out,2,6);

	CopyString(sat->designator,line1_out,9,16);

	sprintf(string,"%02d",sat->year);
	CopyString(string,line1_out,18,19);

	sprintf(string,"%12.8f",sat->refepoch);
	CopyString(string,line1_out,20,32);

	sprintf(string,"%.9f",fabs(sat->drag));

	CopyString(string,line1_out,33,42);

	if (sat->drag<0.0)
		line1_out[33]='-';
	else
		line1_out[33]=' ';

	CopyString(noradEvalue(sat->nddot6),line1_out,44,51);
	CopyString(noradEvalue(sat->bstar),line1_out,53,60);

	sprintf(string,"%4lu",sat->setnum);
	CopyString(string,line1_out,64,67);

	sprintf(string,"%9.4f",sat->incl);
	CopyString(string,line2_out,7,15);
				
	sprintf(string,"%9.4f",sat->raan);
	CopyString(string,line2_out,16,24);

	sprintf(string,"%13.12f",sat->eccn);
	
	/* Erase eccentricity's decimal point */

	for (i=2; i<=9; string[i-2]=string[i], i++);

	CopyString(string,line2_out,26,32);

	sprintf(string,"%9.4f",sat->argper);
	CopyString(string,line2_out,33,41);

	sprintf(string,"%9.5f",sat->meanan);
	CopyString(string,line2_out,43,50);

	sprintf(string,"%12.9f",sat->meanmo);
	CopyString(string,line2_out,52,62);

	sprintf(string,"%5lu",sat->orbitnum);
	CopyString(string,line2_out,63,67);

	/* Compute and insert checksum for line 1 and line 2 */

	for (i=0, sum=0; i<=67; sum+=val[(int)line1_out[i]], i++);

	line1_out[68]=(sum%10)+'0';

	for (i=0, sum=0; i<=67; sum+=val[(int)line2_out[i]], i++);

	line2_out[68]=(sum%10)+'0';

	line1_out[69]=0;
	line2_out[69]=0;
}

double ReadBearing(const char *input)
{
	/* This function takes numeric input in the form of a character
	   string, and returns an equivalent bearing in degrees as a
	   decimal number (double).  The input may either be expressed
	   in decimal format (74.2467) or degree, minute, second
	   format (74 14 48).  This function also safely handles
	   extra spaces found either leading, trailing, or
	   embedded within the numbers expressed in the
	   input string.  Decimal seconds are permitted. */
 
	char string[20];
	double bearing=0.0, seconds;
	int a, b, degrees, minutes;
	size_t length;

	/* Copy "input" to "string", and ignore any extra
	   spaces that might be present in the process. */

	string[0]=0;
	length=strlen(input);

	for (a=0, b=0; a<length && a<18; a++)
	{
		if ((input[a]!=32 && input[a]!='\n') || (input[a]==32 && input[a+1]!=32 && b!=0))
		{
			string[b]=input[a];
			b++;
		}	 
	}

	string[b]=0;

	/* Count number of spaces in the clean string. */

	length=strlen(string);

	for (a=0, b=0; a<length; a++)
		if (string[a]==32)
			b++;

	if (b==0)  /* Decimal Format (74.2467) */
		sscanf(string,"%lf",&bearing);

	if (b==2)  /* Degree, Minute, Second Format (74 14 48) */
	{
		sscanf(string,"%d %d %lf",&degrees, &minutes, &seconds);

		if (degrees<0.0)
		{
			minutes=-minutes;
			seconds=-seconds;
		}

		bearing=(double)degrees+((double)minutes/60)+(seconds/3600);
	}

	/* Bizarre results return a 0.0 */

	if (bearing>360.0 || bearing<-360.0)
		bearing=0.0;

	return bearing;
}

char ReadDataFiles(void)
{
	/* This function reads "predict.qth", "predict.tle",
	   and (optionally) "predict.db" files into memory.
  	   Return values are as follows:

	   0: Neither the qth nor the tle files were loaded
	   1: Only the qth file was loaded
	   2: Only the tle file was loaded
	   3: The qth and tle files were loaded successfully */

	char flag=0;
	FILE *fd;
	
	if ((fd=fopen(qthfile,"r")))
	{
		fgets(qth.callsign,16,fd);
		qth.callsign[strlen(qth.callsign)-1]=0;
		fscanf(fd,"%lf", &qth.stnlat);
		fscanf(fd,"%lf", &qth.stnlong);
		fscanf(fd,"%d",  &qth.stnalt);
		fclose(fd);

		obs_geodetic.lat=qth.stnlat*deg2rad;
		obs_geodetic.lon=-qth.stnlong*deg2rad;
		obs_geodetic.alt=((double)qth.stnalt)/1000.0;
		obs_geodetic.theta=0.0;

		flag=1;
	}

	memset(sat_table, 0, sizeof(sat_table));
	if ((fd=fopen(tlefile,"r")))
	{
		int x=0;
		while (x<24 && feof(fd)==0)
		{
			char name[80];
			char line1[80];
			char line2[80];
			/* Initialize variables */

			name[0]=0;
			line1[0]=0;
			line2[0]=0;

			/* Read element set */

			fgets(name,75,fd);
			fgets(line1,75,fd);
			fgets(line2,75,fd);

			if (KepCheck(line1,line2) && (feof(fd)==0))
			{
				/* We found a valid TLE! */

				/* Some TLE sources left justify the sat
				   name in a 24-byte field that is padded
				   with blanks.  The following lines cut
				   out the blanks as well as the line feed
				   character read by the fgets() function. */
 
				ssize_t y=(ssize_t)strlen(name);

				while (name[y]==' ' || name[y]==0 || name[y]=='\n' || name[y]=='\r' || y==0)
				{
					name[y]=0;
					y--;
				}
				
				/* Copy TLE data into the sat data structure */

				strncpy(sat_table[x].name,name,24);
				strncpy(sat_table[x].line1,line1,69);
				strncpy(sat_table[x].line2,line2,69);

				/* Update individual parameters */

				InternalUpdate(&sat_table[x]);

				x++;
			}
		}

		fclose(fd);
		flag+=2;
		resave=0;

		/* Load satellite database file */
		if ((fd=fopen(dbfile,"r")))
		{
			database=1;
			memset(sat_db,0, sizeof(sat_db));

			char line[80];
			fgets(line,sizeof(line),fd);

			while (strncmp(line,"end",3)!=0 && line[0]!='\n' && !feof(fd))
			{
				/* The first line is the satellite
				   name which is ignored here. */

				fgets(line,sizeof(line),fd);
				long catnum;
				sscanf(line,"%ld",&catnum);

				/* Search for match */
				struct sat_db_st      *current_sat         = NULL;
				struct transponder_st *current_transponder = NULL;
				for (size_t y=0; y<24; y++)
				{
					if (catnum==sat_table[y].catnum)
					{
						current_sat=&sat_db[y];
						current_sat->transponders = 0;
						current_transponder = current_sat->transponder;
						sat_table[y].db = current_sat;
						break;
					}
				}

				fgets(line,40,fd);
				if (current_sat)
				{
					if (strncmp(line,"No",2)!=0)
					{
						sscanf(line,"%lf, %lf",&current_sat->alat, &current_sat->alon);
						current_sat->squintflag=1;
					}
					else
						current_sat->squintflag=0;
				}

				fgets(line,sizeof(line),fd);
				while (strncmp(line,"end",3)!=0 && line[0]!='\n' && !feof(fd))
				{
					if (current_transponder)
					{
						if (strncmp(line,"No",2)!=0)
						{
							line[strlen(line)-1]=0;
							strcpy(current_transponder->name,line);
						}
						else
							current_transponder->name[0]=0;
					}

					fgets(line,sizeof(line),fd);
					if (current_transponder)
						sscanf(line,"%lf, %lf", &current_transponder->uplink_start, &current_transponder->uplink_end);

					fgets(line,sizeof(line),fd);
					if (current_transponder)
						sscanf(line,"%lf, %lf", &current_transponder->downlink_start, &current_transponder->downlink_end);

					fgets(line,sizeof(line),fd);
					if (current_transponder)
					{
						if (strncmp(line,"No",2)!=0)
						{
							unsigned char dayofweek=(unsigned char)atoi(line);
							current_transponder->dayofweek=dayofweek;
						}
						else
							current_transponder->dayofweek=0;
					}

					fgets(line,sizeof(line),fd);
					if (current_transponder)
					{
						if (strncmp(line,"No",2)!=0)
							sscanf(line,"%d, %d",&current_transponder->phase_start, &current_transponder->phase_end);
						else
						{
							current_transponder->phase_start=0;
							current_transponder->phase_end=0;
						}

						if (current_transponder->uplink_start  != 0.0 || current_transponder->downlink_start!= 0.0)
						{
							if (++current_sat->transponders >= 10)
								current_transponder = NULL;
							else
								current_transponder++;
						}
					}
					fgets(line,sizeof(line),fd);
				}
				fgets(line,sizeof(line),fd);
			}

			fclose(fd);
		}
	}

	return flag;
}

char CopyFile(const char *source, const char *destination)
{
	/* This function copies file "source" to file "destination"
	   in 64k chunks.  The permissions on the destination file
	   are set to rw-r--r--  (0644).  A 0 is returned if no
	   errors are encountered.  A 1 indicates a problem writing
	   to the destination file.  A 2 indicates a problem reading
	   the source file.  */

	char error=0;
	int sd=open(source,O_RDONLY);

	if (sd!=-1)
	{
		int dd=open(destination,O_WRONLY | O_CREAT| O_TRUNC, 0644);

		if (dd!=-1)
		{
			char buffer[65536];
			size_t x=read(sd,&buffer,65536);

			while (x)
			{
				write(dd,&buffer,x);
				x=read(sd,&buffer,65536);
			}

			close(dd);
		}
		else
			error=1;

		close(sd);
	}
	else
		error+=2;

	return error;
}

void SaveQTH(void)
{
	/* This function saves QTH data to the QTH data file. */

	FILE *fd=fopen(qthfile,"w");

	fprintf(fd,"%s\n",qth.callsign);
	fprintf(fd," %g\n",qth.stnlat);
	fprintf(fd," %g\n",qth.stnlong);
	fprintf(fd," %d\n",qth.stnalt);

	fclose(fd);
}

void SaveTLE(void)
{
	int x;
 	/* Save orbital data to tlefile */
	FILE *fd=fopen(tlefile,"w");
	if (!fd)
		return;
	for (x=0; x<24; x++)
	{
		/* Convert numeric orbital data to ASCII TLE format */

		char line1[70], line2[70];
		Data2TLE(&sat_table[x], line1, line2);

		/* Write name, line1, line2 to predict.tle */

		fprintf(fd,"%s\n", sat_table[x].name);
		fprintf(fd,"%s\n", line1);
		fprintf(fd,"%s\n", line2);
	}

	fclose(fd);
}

int AutoUpdate(const char * string)
{
	/* This function updates PREDICT's orbital datafile from a NASA
	   2-line element file either through a menu (interactive mode)
	   or via the command line.  string==filename of 2-line element
	   set if this function is invoked via the command line. */

	char filename[50], saveflag=0, interactive=0, savecount=0;

	float database_epoch=0.0, tle_epoch=0.0;
	int success=0, kepcount=0;
	FILE *fd;

	do
	{
		if (!*string)
		{
			int i;
			interactive=1;
			curs_set(1);
			bkgdset(COLOR_PAIR(3));
			refresh();
			clear();
			echo();

			for (i=5; i<8; i+=2)
				mvprintw(i,19,"------------------------------------------");

			mvprintw(6,19,"* Keplerian Database Auto Update Utility *");
			bkgdset(COLOR_PAIR(2));
			mvprintw(19,18,"Enter NASA Two-Line Element Source File Name");
			mvprintw(13,18,"-=> ");
			refresh();
			wgetnstr(stdscr,filename,49);
			clear();
			curs_set(0);
		}
		else
			strcpy(filename,string);

		/* Prevent "." and ".." from being used as a
		   filename, otherwise strange things happen. */

		if (!filename[0] || strncmp(filename,".",1)==0 || strncmp(filename,"..",2)==0)
			return 0;

		fd=fopen(filename,"r");

		if (interactive && !fd)
		{
			int i;
			bkgdset(COLOR_PAIR(5));
			clear();
			move(12,0);

			for (i=47; i>strlen(filename); i-=2)
				printw(" ");

			printw("*** ERROR: File \"%s\" not found! ***\n",filename);
			beep();
			attrset(COLOR_PAIR(7)|A_BOLD);
			AnyKey();
		}

		if (fd)
		{
			char str0[80], str1[80], str2[80];
			success=1;

			fgets(str0,75,fd);
			fgets(str1,75,fd);
			fgets(str2,75,fd);
		
			do
			{
				if (KepCheck(str1,str2))
				{
					int i;
					char line1[80], line2[80];
					/* We found a valid TLE!
					   Copy strings str1 and
					   str2 into line1 and line2 */

					strncpy(line1,str1,75);
					strncpy(line2,str2,75);
					kepcount++;

					/* Scan for object number in datafile to see
					   if this is something we're interested in */

					for (i=0; (i<24 && sat_table[i].catnum!=atol(SubString(line1,2,6))); i++);

					if (i!=24)
					{
						float database_year;
						float tle_year;
						/* We found it!  Check to see if it's more
						   recent than the data we already have. */

						if (sat_table[i].year<57)
							database_year=365.25*(100.0+(float)sat_table[i].year);
						else
							database_year=365.25*(float)sat_table[i].year;

						database_epoch=(float)sat_table[i].refepoch+database_year;

						tle_year=(float)atof(SubString(line1,18,19));

						if (tle_year<57.0)
							tle_year+=100.0;

						tle_epoch=(float)atof(SubString(line1,20,31))+(tle_year*365.25);

						/* Update only if TLE epoch >= epoch in data file
						   so we don't overwrite current data with older
						   data. */

						if (tle_epoch>=database_epoch)
						{
							if (saveflag==0)
							{
								if (interactive)
								{
									clear();
									bkgdset(COLOR_PAIR(2));
									mvprintw(3,35,"Updating.....");
									refresh();
									move(7,0);
								}
								saveflag=1;
							}

							if (interactive)
							{
								bkgdset(COLOR_PAIR(3));
								printw("     %-15s",sat_table[i].name);
							}

							savecount++;

							/* Copy TLE data into the sat data structure */

							strncpy(sat_table[i].line1,line1,69);
							strncpy(sat_table[i].line2,line2,69);
							InternalUpdate(&sat_table[i]);
						}
					}

					 fgets(str0,75,fd);     	
					 fgets(str1,75,fd);
					 fgets(str2,75,fd);
				}

				else
				{
					strcpy(str0,str1);   
					strcpy(str1,str2);   
					fgets(str2,75,fd);
				}
			
			} while (feof(fd)==0);

			fclose(fd);

			if (interactive)
			{
				bkgdset(COLOR_PAIR(2));

				if (kepcount==1)
					mvprintw(18,21,"  Only 1 NASA Two Line Element was found.");
				else
					mvprintw(18,21,"%3u NASA Two Line Elements were read.",kepcount);

				if (saveflag)
				{
					if (savecount==1)
						mvprintw(19,21,"  Only 1 satellite was updated.");
					else
					{
						if (savecount==24)
							mvprintw(19,21,"  All satellites were updated!");
						else
							mvprintw(19,21,"%3u out of 24 satellites were updated.",savecount);
					}
				}

				refresh();
			}
		}

		if (interactive)
		{
			noecho();

			if (strlen(filename) && fd!=NULL) 
			{
				attrset(COLOR_PAIR(4)|A_BOLD);
				AnyKey();
			}
		}

		if (saveflag)
			SaveTLE();
	}
	while (success==0 && interactive);

	return (saveflag ? 0 : -1);
}

struct sat_st * Select(void)
{
	/* This function displays the names of satellites contained
	   within the program's database and returns an index that
	   corresponds to the satellite selected by the user.  An
	   ESC or CR returns a -1. */

	int x, y, z, key=0;

	clear();

	bkgdset(COLOR_PAIR(2)|A_BOLD);
	printw("\n\n\t\t\t      Select a Satellite:\n\n");

	attrset(COLOR_PAIR(3)|A_BOLD);

	for (x=0, y=8, z=16; y<16; ++x, ++y, ++z)
	{
		printw("\n\t[%c]: %-15s", x+'A', Abbreviate(sat_table[x].name,15));
		printw("\t[%c]: %-15s", y+'A', Abbreviate(sat_table[y].name,15));
		printw("\t[%c]: %-15s\n", z+'A', Abbreviate(sat_table[z].name,15));
	}

	attrset(COLOR_PAIR(4)|A_BOLD);

	printw("\n\n\t\t<< Enter Selection  -  Press [ESC] To Exit >>");
	refresh();

	do
	{
		key=toupper(getch());

		if (key==27 || key=='\n')
			return NULL;

	} while (key<'A' || key>'X');

	return &sat_table[key-'A'];
}

long DayNum(int m,int d,int y)
{
	/* This function calculates the day number from m/d/y. */

	long dn;
	double mm, yy;

	if (m<3)
	{ 
		y--; 
		m+=12; 
	}

	if (y<57)
		y+=100;

	yy=(double)y;
	mm=(double)m;
	dn=(long)(floor(365.25*(yy-80.0))-floor(19.0+yy/100.0)+floor(4.75+yy/400.0)-16.0);
	dn+=d+30*m+(long)floor(0.6*mm-0.3);
	return dn;
}

double CurrentDaynum(void)
{
	/* Read the system clock and return the number
	   of days since 31Dec79 00:00:00 UTC (daynum 0) */

	/* int x; */
	struct timeval tptr;
	double usecs, seconds;

	/* x=gettimeofday(&tptr,NULL); */
	(void)gettimeofday(&tptr,NULL);

	usecs=0.000001*(double)tptr.tv_usec;
	seconds=usecs+(double)tptr.tv_sec;

	return ((seconds/86400.0)-3651.0);
}

const char *Daynum2String(double daynum)
{
	/* This function takes the given epoch as a fractional number of
	   days since 31Dec79 00:00:00 UTC and returns the corresponding
	   date as a string of the form "Tue 12Oct99 17:22:37". */


	/* Convert daynum to Unix time (seconds since 01-Jan-70) */
	const time_t t=(time_t)(86400.0*(daynum+3651.0));

	char timestr[26];
	sprintf(timestr,"%s",asctime(gmtime(&t)));

	if (timestr[8]==' ')
		timestr[8]='0';

	static char output[25];
	output[0]=timestr[0];
	output[1]=timestr[1];
	output[2]=timestr[2];
	output[3]=timestr[3];
	output[4]=timestr[8];
	output[5]=timestr[9];
	output[6]=timestr[4];
	output[7]=timestr[5];
	output[8]=timestr[6];
	output[9]=timestr[22];
	output[10]=timestr[23];
	output[11]=' ';

	for (int x=12; x<20; x++) output[x]=timestr[x-1];

	output[20]=0;
	return output;
}

double GetStartTime(const char *object)
{
	/* This function prompts the user for the time and date
	   the user wishes to begin prediction calculations,
	   and returns the corresponding fractional day number.
	   31Dec79 00:00:00 returns 0.  Default is NOW. */

	static const char *month[12]= {"Jan", "Feb", "Mar", "Apr", "May",
		"Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"};

	char string[30];
	int	 bozo_count=0;
	do
	{
		bkgdset(COLOR_PAIR(2)|A_BOLD);
		clear();

		printw("\n\n\n\t     Starting UTC Date and Time for Predictions of %s\n\n", object);

		bozo_count++;

		strcpy(string,Daynum2String(CurrentDaynum()));

		for (int x=4; x<24; x++)
			string[x-4]=string[x];

		attrset(COLOR_PAIR(4)|A_BOLD);
		printw("\t\t    Format: %s -or- ",string);
		string[7]=0;
		printw("%s",string);
		attrset(COLOR_PAIR(2)|A_BOLD);
		mvprintw(21,30,"Default is `NOW'");
		attrset(COLOR_PAIR(3)|A_BOLD);
		mvprintw(13,1,"Enter Start Date & Time >> ");
		curs_set(1);
		refresh();
		echo();
		string[0]=0;
		wgetnstr(stdscr,string,29);
		curs_set(0);
		noecho();

		if (!strlen(string)) /* Select `NOW' */
			return(CurrentDaynum());

		if (strlen(string)==7)
			strcat(string, " 00:00:00");

		const int good = isdigit(string[ 0]) && isdigit(string[ 1]) && /* Check Day */
		                 isalpha(string[ 2]) && isalpha(string[ 3]) && isalpha(string[4]) && /* Month */
		                 isdigit(string[ 5]) && isdigit(string[ 6]) && (string[ 7]==' ') && /* Year */
		                 isdigit(string[ 8]) && isdigit(string[ 9]) && (string[10]==':') && /* Hour */
		                 isdigit(string[11]) && isdigit(string[12]) && (string[13]==':') && /* Minute */
		                 isdigit(string[14]) && isdigit(string[15]); /* Seconds */

		/* Decode Month Number */
		string[2]=toupper(string[2]);
		string[3]=tolower(string[3]);
		string[4]=tolower(string[4]);

		int	mm=0;
		for (mm=1; mm<=12; mm++)
			if (strncmp(&string[2],month[mm-1],3)==0)
				break;

		if (good && mm <= 12)
		{
			const int dd  = 10*(string[ 0]-'0')+string[ 1]-'0';
			const int yy  = 10*(string[ 5]-'0')+string[ 6]-'0';
			const int hr  = 10*(string[ 8]-'0')+string[ 9]-'0';
			const int min = 10*(string[11]-'0')+string[12]-'0';
			const int sec = 10*(string[14]-'0')+string[15]-'0';

			return ((double)DayNum(mm,dd,yy)+((hr/24.0)+(min/1440.0)+(sec/86400.0)));
		}

		beep();
	} while (bozo_count<6);

	/* If the user can't enter the starting date/time
	   correctly after several attempts, then the user
	   is a "bozo" and obviously can't follow directions. */

	bailout("Too Many Errors");
	exit(-1);
}

void FindMoon(double daynum)
{
	/* This function determines the position of the moon, including
	   the azimuth and elevation headings, relative to the latitude
	   and longitude of the tracking station.  This code was derived
	   from a Javascript implementation of the Meeus method for
	   determining the exact position of the Moon found at:
	   http://www.geocities.com/s_perona/ingles/poslun.htm. */

	double	jd, ss, t, t1, t2, t3, d, ff, l1, m, m1, ex, om, l,
		b, w1, w2, bt, p, lm, h, ra, dec, z, ob, n, el, az,
		teg, th, mm, dv;

	jd=daynum+2444238.5;

	t=(jd-2415020.0)/36525.0;
	t2=t*t;
	t3=t2*t;
	l1=270.434164+481267.8831*t-0.001133*t2+0.0000019*t3;
	m=358.475833+35999.0498*t-0.00015*t2-0.0000033*t3;
	m1=296.104608+477198.8491*t+0.009192*t2+0.0000144*t3;
	d=350.737486+445267.1142*t-0.001436*t2+0.0000019*t3;
	ff=11.250889+483202.0251*t-0.003211*t2-0.0000003*t3;
	om=259.183275-1934.142*t+0.002078*t2+0.0000022*t3;
	om=om*deg2rad;
	
	/* Additive terms */

	l1=l1+0.000233*sin((51.2+20.2*t)*deg2rad);
	ss=0.003964*sin((346.56+132.87*t-0.0091731*t2)*deg2rad);
	l1=l1+ss+0.001964*sin(om);
	m=m-0.001778*sin((51.2+20.2*t)*deg2rad);
	m1=m1+0.000817*sin((51.2+20.2*t)*deg2rad);
	m1=m1+ss+0.002541*sin(om);
	d=d+0.002011*sin((51.2+20.2*t)*deg2rad);
	d=d+ss+0.001964*sin(om);
	ff=ff+ss-0.024691*sin(om);
	ff=ff-0.004328*sin(om+(275.05-2.3*t)*deg2rad);
	ex=1.0-0.002495*t-0.00000752*t2;
	om=om*deg2rad;

	l1=PrimeAngle(l1);
	m=PrimeAngle(m);
	m1=PrimeAngle(m1);
	d=PrimeAngle(d);
	ff=PrimeAngle(ff);
	om=PrimeAngle(om);

	m=m*deg2rad;
	m1=m1*deg2rad;
	d=d*deg2rad;
	ff=ff*deg2rad;

	/* Ecliptic Longitude */

	l=l1+6.28875*sin(m1)+1.274018*sin(2.0*d-m1)+0.658309*sin(2.0*d);
	l=l+0.213616*sin(2.0*m1)-ex*0.185596*sin(m)-0.114336*sin(2.0*ff);
	l=l+0.058793*sin(2.0*d-2.0*m1)+ex*0.057212*sin(2.0*d-m-m1)+0.05332*sin(2.0*d+m1);
	l=l+ex*0.045874*sin(2.0*d-m)+ex*0.041024*sin(m1-m)-0.034718*sin(d);
	l=l-ex*0.030465*sin(m+m1)+0.015326*sin(2.0*d-2.0*ff)-0.012528*sin(2.0*ff+m1);
	
	l=l-0.01098*sin(2.0*ff-m1)+0.010674*sin(4.0*d-m1)+0.010034*sin(3.0*m1);
	l=l+0.008548*sin(4.0*d-2.0*m1)-ex*0.00791*sin(m-m1+2.0*d)-ex*0.006783*sin(2.0*d+m);
	
	l=l+0.005162*sin(m1-d)+ex*0.005*sin(m+d)+ex*0.004049*sin(m1-m+2.0*d);
	l=l+0.003996*sin(2.0*m1+2.0*d)+0.003862*sin(4.0*d)+0.003665*sin(2.0*d-3.0*m1);

	l=l+ex*0.002695*sin(2.0*m1-m)+0.002602*sin(m1-2.0*ff-2.0*d)+ex*0.002396*sin(2.0*d-m-2.0*m1);

	l=l-0.002349*sin(m1+d)+ex*ex*0.002249*sin(2.0*d-2.0*m)-ex*0.002125*sin(2.0*m1+m);

	l=l-ex*ex*0.002079*sin(2.0*m)+ex*ex*0.002059*sin(2.0*d-m1-2.0*m)-0.001773*sin(m1+2.0*d-2.0*ff);

	l=l+ex*0.00122*sin(4.0*d-m-m1)-0.00111*sin(2.0*m1+2.0*ff)+0.000892*sin(m1-3.0*d);

	l=l-ex*0.000811*sin(m+m1+2.0*d)+ex*0.000761*sin(4.0*d-m-2.0*m1)+ex*ex*.000717*sin(m1-2.0*m);

	l=l+ex*ex*0.000704*sin(m1-2.0*m-2.0*d)+ex*0.000693*sin(m-2.0*m1+2.0*d)+ex*0.000598*sin(2.0*d-m-2.0*ff)+0.00055*sin(m1+4.0*d);

	l=l+0.000538*sin(4.0*m1)+ex*0.000521*sin(4.0*d-m)+0.000486*sin(2.0*m1-d);

	l=l-0.001595*sin(2.0*ff+2.0*d);

	/* Ecliptic latitude */

	b=5.128189*sin(ff)+0.280606*sin(m1+ff)+0.277693*sin(m1-ff)+0.173238*sin(2.0*d-ff);
	b=b+0.055413*sin(2.0*d+ff-m1)+0.046272*sin(2.0*d-ff-m1)+0.032573*sin(2.0*d+ff);

	b=b+0.017198*sin(2.0*m1+ff)+9.266999e-03*sin(2.0*d+m1-ff)+0.008823*sin(2.0*m1-ff);
	b=b+ex*0.008247*sin(2.0*d-m-ff)+0.004323*sin(2.0*d-ff-2.0*m1)+0.0042*sin(2.0*d+ff+m1);

	b=b+ex*0.003372*sin(ff-m-2.0*d)+ex*0.002472*sin(2.0*d+ff-m-m1)+ex*0.002222*sin(2.0*d+ff-m);

	b=b+0.002072*sin(2.0*d-ff-m-m1)+ex*0.001877*sin(ff-m+m1)+0.001828*sin(4.0*d-ff-m1);

	b=b-ex*0.001803*sin(ff+m)-0.00175*sin(3.0*ff)+ex*0.00157*sin(m1-m-ff)-0.001487*sin(ff+d)-ex*0.001481*sin(ff+m+m1)+ex*0.001417*sin(ff-m-m1)+ex*0.00135*sin(ff-m)+0.00133*sin(ff-d);

	b=b+0.001106*sin(ff+3.0*m1)+0.00102*sin(4.0*d-ff)+0.000833*sin(ff+4.0*d-m1);

	b=b+0.000781*sin(m1-3.0*ff)+0.00067*sin(ff+4.0*d-2.0*m1)+0.000606*sin(2.0*d-3.0*ff);

	b=b+0.000597*sin(2.0*d+2.0*m1-ff)+ex*0.000492*sin(2.0*d+m1-m-ff)+0.00045*sin(2.0*m1-ff-2.0*d);

	b=b+0.000439*sin(3.0*m1-ff)+0.000423*sin(ff+2.0*d+2.0*m1)+0.000422*sin(2.0*d-ff-3.0*m1);

	b=b-ex*0.000367*sin(m+ff+2.0*d-m1)-ex*0.000353*sin(m+ff+2.0*d)+0.000331*sin(ff+4.0*d);

	b=b+ex*0.000317*sin(2.0*d+ff-m+m1)+ex*ex*0.000306*sin(2.0*d-2.0*m-ff)-0.000283*sin(m1+3.0*ff);

	w1=0.0004664*cos(om*deg2rad);
	w2=0.0000754*cos((om+275.05-2.3*t)*deg2rad);
	bt=b*(1.0-w1-w2);

	/* Parallax calculations */

	p=0.950724+0.051818*cos(m1)+0.009531*cos(2.0*d-m1)+0.007843*cos(2.0*d)+0.002824*cos(2.0*m1)+0.000857*cos(2.0*d+m1)+ex*0.000533*cos(2.0*d-m)+ex*0.000401*cos(2.0*d-m-m1);

	p=p+0.000173*cos(3.0*m1)+0.000167*cos(4.0*d-m1)-ex*0.000111*cos(m)+0.000103*cos(4.0*d-2.0*m1)-0.000084*cos(2.0*m1-2.0*d)-ex*0.000083*cos(2.0*d+m)+0.000079*cos(2.0*d+2.0*m1);

	p=p+0.000072*cos(4.0*d)+ex*0.000064*cos(2.0*d-m+m1)-ex*0.000063*cos(2.0*d+m-m1);

	p=p+ex*0.000041*cos(m+d)+ex*0.000035*cos(2.0*m1-m)-0.000033*cos(3.0*m1-2.0*d);

	p=p-0.00003*cos(m1+d)-0.000029*cos(2.0*ff-2.0*d)-ex*0.000029*cos(2.0*m1+m);

	p=p+ex*ex*0.000026*cos(2.0*d-2.0*m)-0.000023*cos(2.0*ff-2.0*d+m1)+ex*0.000019*cos(4.0*d-m-m1);

	b=bt*deg2rad;
	lm=l*deg2rad;
	moon_dx=3.0/(pi*p);

	/* Semi-diameter calculation */
	/* sem=10800.0*asin(0.272488*p*deg2rad)/pi; */

	/* Convert ecliptic coordinates to equatorial coordinates */

	z=(jd-2415020.5)/365.2422;
	ob=23.452294-(0.46845*z+5.9e-07*z*z)/3600.0;
	ob=ob*deg2rad;
	dec=asin(sin(b)*cos(ob)+cos(b)*sin(ob)*sin(lm));
	ra=acos(cos(b)*cos(lm)/cos(dec));
	
	if (lm>pi)
		ra=twopi-ra;

	/* ra = right ascension */
	/* dec = declination */

	n=qth.stnlat*deg2rad;    /* North latitude of tracking station */

	/* Find siderial time in radians */

	t=(jd-2451545.0)/36525.0;
	teg=280.46061837+360.98564736629*(jd-2451545.0)+(0.000387933*t-t*t/38710000.0)*t;

	while (teg>360.0)
		teg-=360.0;

	th=FixAngle((teg-qth.stnlong)*deg2rad);
	h=th-ra;

	az=atan2(sin(h),cos(h)*sin(n)-tan(dec)*cos(n))+pi;
	el=asin(sin(n)*sin(dec)+cos(n)*cos(dec)*cos(h));

	moon_az=az/deg2rad;
	moon_el=el/deg2rad;

	/* Radial velocity approximation.  This code was derived
	   from "Amateur Radio Software", by John Morris, GM4ANB,
	   published by the RSGB in 1985. */

	mm=FixAngle(1.319238+daynum*0.228027135);  /* mean moon position */
	t2=0.10976;
	t1=mm+t2*sin(mm);
	dv=0.01255*moon_dx*moon_dx*sin(t1)*(1.0+t2*cos(mm));
	dv=dv*4449.0;
	t1=6378.0;
	t2=384401.0;
	t3=t1*t2*(cos(dec)*cos(n)*sin(h));
	t3=t3/sqrt(t2*t2-t2*t1*sin(el));
	moon_dv=dv+t3*0.0753125;

	moon_dec=dec/deg2rad;
	moon_ra=ra/deg2rad;
	moon_gha=teg-moon_ra;

	if (moon_gha<0.0)
		moon_gha+=360.0;
}

void FindSun(double daynum)
{
	/* This function finds the position of the Sun */

	/* Zero vector for initializations */
	static const vector_t zero_vector={0,0,0,0};

	/* Solar ECI position vector  */
	vector_t solar_vector=zero_vector;

	/* Solar observed azi and ele vector  */
	vector_t solar_set=zero_vector;

	/* Solar right ascension and declination vector */
	vector_t solar_rad=zero_vector;

	/* Solar lat, long, alt vector */
	geodetic_t solar_latlonalt;

	const double jul_utc=daynum+2444238.5;

	Calculate_Solar_Position(jul_utc, &solar_vector);
	Calculate_Obs(jul_utc, &solar_vector, &zero_vector, &obs_geodetic, &solar_set);
	sun_azi=Degrees(solar_set.x); 
	sun_ele=Degrees(solar_set.y);
	sun_range=1.0+((solar_set.z-AU)/AU);
	sun_range_rate=1000.0*solar_set.w;

	Calculate_LatLonAlt(jul_utc, &solar_vector, &solar_latlonalt);

	sun_lat=Degrees(solar_latlonalt.lat);
	sun_lon=360.0-Degrees(solar_latlonalt.lon);

	Calculate_RADec(jul_utc, &solar_vector, &zero_vector, &obs_geodetic, &solar_rad);

	sun_ra=Degrees(solar_rad.x);
	sun_dec=Degrees(solar_rad.y);
}

void PreCalc(const struct sat_st *sat, tle_t *tle_out)
{
	/* This function copies TLE data from PREDICT's sat structure
	   to the SGP4/SDP4's single dimensioned tle structure, and
	   prepares the tracking code for the update. */
	static const double temp=twopi/xmnpda/xmnpda;

	strcpy(tle_out->sat_name,sat->name);
	strcpy(tle_out->idesg,sat->designator);
	tle_out->catnr  = sat->catnum;
	tle_out->epoch  = (1000.0*(double)sat->year)+sat->refepoch;
	tle_out->xndt2o = sat->drag * temp;
	tle_out->xndd6o = sat->nddot6 * temp/xmnpda;
	tle_out->bstar  = sat->bstar / ae;
	tle_out->xincl  = sat->incl * deg2rad;
	tle_out->xnodeo = sat->raan * deg2rad;
	tle_out->eo     = sat->eccn;
	tle_out->omegao = sat->argper * deg2rad;
	tle_out->xmo    = sat->meanan * deg2rad;
	tle_out->xno    = sat->meanmo * temp*xmnpda;
	tle_out->revnum = sat->orbitnum;

	if (sat->db && sat->db->squintflag)
	{
		calc_squint=1;
		alat=deg2rad*sat->db->alat;
		alon=deg2rad*sat->db->alon;
	}
	else
		calc_squint=0;
 
	/* Clear all flags */

	ClearFlag(ALL_FLAGS);

	/* Select ephemeris type.  This function will set or clear the
	   DEEP_SPACE_EPHEM_FLAG depending on the TLE parameters of the
	   satellite.  It will also pre-process tle members for the
	   ephemeris functions SGP4 or SDP4, so this function must
	   be called each time a new tle set is used. */

	if (select_ephemeris(tle_out))
		SetFlag(DEEP_SPACE_EPHEM_FLAG);
	else
		ClearFlag(DEEP_SPACE_EPHEM_FLAG);
}

void Calc(const tle_t *tle)
{
	/* This is the stuff we need to do repetitively while tracking. */

	/* Zero vector for initializations */
	vector_t zero_vector={0,0,0,0};

	/* Satellite position and velocity vectors */
	vector_t vel=zero_vector;
	vector_t pos=zero_vector;

	/* Satellite Az, El, Range, Range rate */
	vector_t obs_set;

	/* Solar ECI position vector  */
	vector_t solar_vector=zero_vector;

	/* Solar observed azi and ele vector  */
	vector_t solar_set;

	/* Satellite's predicted geodetic position */
	geodetic_t sat_geodetic;

	const double jul_utc=daynum+2444238.5;

	/* Convert satellite's epoch time to Julian  */
	/* and calculate time since epoch in minutes */

	const double jul_epoch=Julian_Date_of_Epoch(tle->epoch);
	const double tsince=(jul_utc-jul_epoch)*xmnpda;
	const double age=jul_utc-jul_epoch;

	/* Call NORAD routines according to deep-space flag. */



	if (isFlagSet(DEEP_SPACE_EPHEM_FLAG))
	{
		static struct SDP4_state_st spd4_state;
		if (isFlagClear(SDP4_INITIALIZED_FLAG))
		{
			SetFlag(SDP4_INITIALIZED_FLAG);
			SDP4_init(&spd4_state, tle);
		}

		SDP4(&spd4_state, tsince, tle, &pos, &vel);
	}
	else
	{
		static struct SGP4_state_st spg4_state;
		if (isFlagClear(SGP4_INITIALIZED_FLAG))
		{
			SetFlag(SGP4_INITIALIZED_FLAG);
			SGP4_init(&spg4_state, tle);
		}
		SGP4(&spg4_state, tsince, tle, &pos, &vel);
	}

	/* Converts the satellite's position and velocity  */
	/* vectors from normalized values to km and km/sec */
	Scale_Vector(xkmper, &pos);
	Scale_Vector(xkmper*xmnpda/secday, &vel);

	/* Calculate velocity of satellite */

	Magnitude(&vel);
	sat_vel=vel.w;

	/** All angles in rads. Distance in km. Velocity in km/s **/
	/* Calculate satellite Azi, Ele, Range and Range-rate */

	Calculate_Obs(jul_utc, &pos, &vel, &obs_geodetic, &obs_set);

	/* Calculate satellite Lat North, Lon East and Alt. */

	Calculate_LatLonAlt(jul_utc, &pos, &sat_geodetic);

	/* Calculate squint angle */

	if (calc_squint)
		squint=(acos(-(ax*rx+ay*ry+az*rz)/obs_set.z))/deg2rad;

	/* Calculate solar position and satellite eclipse depth. */
	/* Also set or clear the satellite eclipsed flag accordingly. */

	Calculate_Solar_Position(jul_utc, &solar_vector);
	Calculate_Obs(jul_utc, &solar_vector, &zero_vector, &obs_geodetic, &solar_set);

	if (Sat_Eclipsed(&pos, &solar_vector, &eclipse_depth))
		SetFlag(SAT_ECLIPSED_FLAG);
	else
		ClearFlag(SAT_ECLIPSED_FLAG);

	if (isFlagSet(SAT_ECLIPSED_FLAG))
		sat_sun_status=0;  /* Eclipse */
	else
		sat_sun_status=1; /* In sunlight */

	/* Convert satellite and solar data */
	sat_azi=Degrees(obs_set.x);
	sat_ele=Degrees(obs_set.y);
	sat_range=obs_set.z;
	sat_range_rate=obs_set.w;
	sat_lat=Degrees(sat_geodetic.lat);
	sat_lon=Degrees(sat_geodetic.lon);
	sat_alt=sat_geodetic.alt;

	sat_footprint=12756.33*acos(xkmper/(xkmper+sat_alt));

	rv=(long)floor((tle->xno*xmnpda/twopi+age*tle->bstar*ae)*age+tle->xmo/twopi)+tle->revnum;

	sun_azi=Degrees(solar_set.x); 
	sun_ele=Degrees(solar_set.y);

	ma256=(int)rint(256.0*(phase/twopi));

	if (sat_sun_status)
	{
		if (sun_ele<=-12.0 && rint(sat_ele)>=0.0)
			findsun='+';
		else
			findsun='*';
	}
	else
		findsun=' ';
}

int AosHappens(struct sat_st const * sat)
{
	/* This function returns a 1 if the satellite pointed to by
	   "x" can ever rise above the horizon of the ground station. */

	double lin, sma, apogee;

	if (sat->meanmo==0.0)
		return 0;
	else
	{
		lin=sat->incl;

		if (lin>=90.0)
			lin=180.0-lin;

		sma=331.25*exp(log(1440.0/sat->meanmo)*(2.0/3.0));
		apogee=sma*(1.0+sat->eccn)-xkmper;

		if ((acos(xkmper/(apogee+xkmper))+(lin*deg2rad)) > fabs(qth.stnlat*deg2rad))
			return 1;
		else
			return 0;
	}
}

int Decayed(const struct sat_st *sat, double time)
{
	/* This function returns a 1 if it appears that the
	   satellite pointed to by 'x' has decayed at the
	   time of 'time'.  If 'time' is 0.0, then the
	   current date/time is used. */

	double satepoch;

	if (time==0.0)
		time=CurrentDaynum();

	satepoch=DayNum(1,0,sat->year)+sat->refepoch;

	if (satepoch+((16.666666-sat->meanmo)/(10.0*fabs(sat->drag))) < time)
		return 1;
	else
		return 0;
}

int Geostationary(struct sat_st const * sat)
{
	/* This function returns a 1 if the satellite pointed
	   to by "x" appears to be in a geostationary orbit */

	if (fabs(sat->meanmo-1.0027)<0.0002)

		return 1;
	else
		return 0;
}

double FindAOS(struct sat_st *sat, const tle_t *tle)
{
	/* This function finds and returns the time of AOS (aostime). */

	aostime=0.0;

	if (AosHappens(sat) && !Geostationary(sat) && !Decayed(sat,daynum))
	{
		Calc(tle);

		/* Get the satellite in range */

		while (sat_ele<-1.0)
		{
			daynum-=0.00035*(sat_ele*((sat_alt/8400.0)+0.46)-2.0);
			Calc(tle);
		}

		/* Find AOS */

		while (aostime==0.0)
		{
			if (fabs(sat_ele)<0.03)
				aostime=daynum;
			else
			{
				daynum-=sat_ele*sqrt(sat_alt)/530000.0;
				Calc(tle);
			}
		}
	}

	return aostime;
}

double FindLOS(struct sat_st *sat, const tle_t *tle)
{
	lostime=0.0;

	if (!Geostationary(sat) && AosHappens(sat) && !Decayed(sat,daynum))
	{
		Calc(tle);

		do
		{
			daynum+=sat_ele*sqrt(sat_alt)/502500.0;
			Calc(tle);

			if (fabs(sat_ele) < 0.03)
				lostime=daynum;

		} while (lostime==0.0);
	}

	return lostime;
}

double FindLOS2(struct sat_st *sat, const tle_t *tle)
{
	/* This function steps through the pass to find LOS.
	   FindLOS() is called to "fine tune" and return the result. */

	do
	{
		daynum+=cos((sat_ele-1.0)*deg2rad)*sqrt(sat_alt)/25000.0;
		Calc(tle);

	} while (sat_ele>=0.0);

	return(FindLOS(sat, tle));
}

double NextAOS(struct sat_st *sat, const tle_t *tle)
{
	/* This function finds and returns the time of the next
	   AOS for a satellite that is currently in range. */

	aostime=0.0;

	if (AosHappens(sat) && !Geostationary(sat) && !Decayed(sat,daynum))
		daynum=FindLOS2(sat, tle)+0.014;  /* Move to LOS + 20 minutes */

	return (FindAOS(sat, tle));
}

int Print(struct sat_st *sat, const char *string,char mode)
{
	/* This function buffers and displays orbital predictions
	   and allows screens to be saved to a disk file. */

	char type[20], spaces[80], head1[160], head2[70],
	     head3[72], satellite_name[25];
	int key, ans=0, t;
	size_t l, x;
	static char buffer[1450], lines, quit;
	static FILE *fd;

	/* Pass a NULL string to initialize the buffer, counter, and flags */

	if (string[0]==0)
	{
		lines=0;
		quit=0;
		buffer[0]=0;
		fd=NULL;
	}

	else
	{
		if (mode=='m')
		{
			sprintf(head1,"\n                    %s's Orbit Calendar for the Moon",qth.callsign);
			strncpy(satellite_name,"Moon\0",5);
		}

		if (mode=='o')
		{
			sprintf(head1,"\n                    %s's Orbit Calendar for the Sun",qth.callsign);
			strncpy(satellite_name,"Sun\0",4);
		}
	
		if (mode=='m' || mode=='o')

			sprintf(head2,"\n\t   Date     Time    El   Az   RA     Dec    GHA     Vel   Range\n");


		if (mode=='p')
			strcpy(type,"Orbit");

		if (mode=='v')
			strcpy(type,"Visual");

		if (mode=='s')
			strcpy(type,"Solar Illumination");

		if (mode!='m' && mode!='o')
		{
			strcpy(satellite_name,sat->name);

			l=strlen(qth.callsign)+strlen(satellite_name)+strlen(type);

			spaces[0]=0;

			for (x=l; x<60; x+=2)
				strcat(spaces," ");

			sprintf(head1,"\n%s%s's %s Calendar For %s", spaces, qth.callsign, type, satellite_name);

			if (mode=='s')
				sprintf(head2,"\n\t  Date     Mins/Day    Sun%c          Date      Mins/Day    Sun%c\n",37,37);
			else
			{
				if (calc_squint)

					sprintf(head2,"\n\t   Date     Time    El   Az  Phase  %s   %s    Range  Squint\n",(io_lat=='N'?"LatN":"LatS"),(io_lon=='W'?"LonW":"LonE"));
				else

					sprintf(head2,"\n\t   Date     Time    El   Az  Phase  %s   %s    Range  Orbit\n",(io_lat=='N'?"LatN":"LatS"),(io_lon=='W'?"LonW":"LonE"));
			}
		}

		sprintf(head3,"      -----------------------------------------------------------------");

		strcat(buffer,string);
		lines++;

		if (lines==18)
		{
			bkgdset(COLOR_PAIR(2)|A_BOLD);
			clear();
			addstr(head1);
			attrset(COLOR_PAIR(4)|A_BOLD);
			addstr(head2);
			addstr(head3);
			attrset(COLOR_PAIR(3)|A_BOLD);

			if (buffer[0]!='\n')
				printw("\n");

			addstr(buffer);
			attrset(COLOR_PAIR(4)|A_BOLD);

			if (buffer[0]=='\n')
				printw("\n");

			if (fd==NULL)
				mvprintw(23,63,"        ");
			else
				mvprintw(23,63,"Log = ON");

			mvprintw(23,6,"More? [y/n] >> ");
			curs_set(1);
			refresh();

			while (ans==0)
			{
				key=toupper(getch());

				if (key=='Y' || key=='\n' || key==' ')
				{
					key='Y';
					ans=1;
					quit=0;
				}
			
				if (key=='N' || key=='Q' || key==27)
				{
					key='N';
					ans=1;
					quit=1;
				}

				/* 'L' logs output to "satname.txt" */

				if (key=='L' && fd==NULL && buffer[0])
				{
					sprintf(temp,"%s.txt",satellite_name);

					l=strlen(temp)-4;

					for (x=0; x<l; x++)
					{
						t=temp[x];

						if (t==32 || t==17 || t==92 || t==42 || t==46 || t==47)
							t='_';

						temp[x]=t;
					}

					fd=fopen(temp,"a");
					fprintf(fd,"%s%s%s\n",head1,head2,head3);
					fprintf(fd,"%s",buffer);
					mvprintw(23,63,"Log = ON");
					move(23,21);
					refresh();
				}

				else if (fd!=NULL)
				{
					if (key=='L' || key=='N')
					{
						fprintf(fd,"%s\n\n",buffer);
						fclose(fd);
						fd=NULL;
						mvprintw(23,63,"        ");
						move(23,21);
						refresh();
					}

					else
						fprintf(fd,"%s",buffer);
				}
				buffer[0]=0;
			}

			lines=0;
			curs_set(0);
		}
	}
	return (quit);
}

int PrintVisible(struct sat_st *sat, const char *string)
{
	/* This function acts as a filter to display passes that could
	   possibly be optically visible to the ground station.  It works
	   by first buffering prediction data generated by the Predict()
	   function and then checking it to see if at least a part of the
	   pass is visible.  If it is, then the buffered prediction data
	   is sent to the Print() function so it can be displayed
	   to the user and optionally logged to a file. */

	static char buffer[10000];
	char line[80], plus, asterisk, visible;
	int x, y, quit=0;

	if (string[0]==0)
		buffer[0]=0;
	else
	{
		strcat(buffer,string);

		if (string[0]=='\n')
		{
			plus=0;
			visible=0;
			asterisk=0;

			for (x=0; buffer[x]!=0 && visible==0; x++)
			{
				if (buffer[x]=='+')
					plus++;

				if (buffer[x]=='*')
					asterisk++;

				/* At least 3 +'s or at least 2 +'s
				   combined with at least 2 *'s is
				   worth displaying as a visible pass. */

				if ((plus>3) || (plus>2 && asterisk>2))
					visible=1;
			}

			if (visible)
			{
				/* Dump buffer to Print() line by line */

				for (x=0, y=0; buffer[x]!=0 && quit==0; x++)
				{	
					line[y]=buffer[x];

					if (line[y]=='\n')
					{
						line[y+1]=0;
						quit=Print(sat,line,'v');
						line[0]=0;
						y=0;
					}

					else
						y++;
				}
			}

			buffer[0]=0;
		}
	}

	return quit;
}

void Predict(struct sat_st *sat, char mode)
{
	/* This function predicts satellite passes.  It displays
	   output through the Print() function if mode=='p' (show
	   all passes), or through the PrintVisible() function if
	   mode=='v' (show optically visible passes only). */

	int quit=0, breakout=0;
	double lastel=0;
	char string[80];

	tle_t tle;
	PreCalc(sat, &tle);
	daynum=GetStartTime(sat->name);
	clear();

	/* Trap geostationary orbits and passes that cannot occur. */

	if (AosHappens(sat) && !Geostationary(sat) && !Decayed(sat,daynum))
	{
		if (xterm)
		{
			char type[10];
			strcpy(type,"Orbit");  /* Default */

			if (mode=='v')
				strcpy(type,"Visual");

			fprintf(stderr,"\033]0;PREDICT: %s's %s Calendar For %s\007",qth.callsign, type, sat->name);
		}

		do
		{
			daynum=FindAOS(sat, &tle);
		
			/* Display the pass */

			while ( rint(sat_ele)>=0 && !quit)
			{
				if (calc_squint)

					sprintf(string,"      %s%4.0f %4.0f  %4d  %4.0f   %4.0f   %6.0f  %4.0f %c\n",
							Daynum2String(daynum),
							sat_ele,sat_azi,ma256,
							(io_lat=='N'?+1:-1)*sat_lat,(io_lon=='W'?360.0-sat_lon:sat_lon),
							sat_range,squint,findsun);

				else
					sprintf(string,"      %s%4.0f %4.0f  %4d  %4.0f   %4.0f   %6.0f  %6ld %c\n",
							Daynum2String(daynum),
							sat_ele,sat_azi,ma256,
							(io_lat=='N'?+1:-1)*sat_lat,(io_lon=='W'?360.0-sat_lon:sat_lon),
							sat_range,rv,findsun);

				lastel=sat_ele;

				if (mode=='p')
					quit=Print(sat, string,'p');

				if (mode=='v')
				{
					nodelay(stdscr,TRUE);
					attrset(COLOR_PAIR(4));
					mvprintw(23,6,"                 Calculating... Press [ESC] To Quit");

					/* Allow a way out if this
					   should continue forever... */

					if (getch()==27)
						breakout=1;

					nodelay(stdscr,FALSE);

					quit=PrintVisible(sat, string);
				}

				daynum+=cos((sat_ele-1.0)*deg2rad)*sqrt(sat_alt)/25000.0;
				Calc(&tle);
			}

			if (rint(lastel)!=0)
			{
				daynum=FindLOS(sat, &tle);
				Calc(&tle);

				if (calc_squint)
					sprintf(string,"      %s%4.0f %4.0f  %4d  %4.0f   %4.0f   %6.0f  %4.0f %c\n",
							Daynum2String(daynum),
							sat_ele,sat_azi,ma256,
							(io_lat=='N'?+1:-1)*sat_lat,(io_lon=='W'?360.0-sat_lon:sat_lon),
							sat_range,squint,findsun);

				else
					sprintf(string,"      %s%4.0f %4.0f  %4d  %4.0f   %4.0f   %6.0f  %6ld %c\n",
							Daynum2String(daynum),
							sat_ele,sat_azi,ma256,
							(io_lat=='N'?+1:-1)*sat_lat,(io_lon=='W'?360.0-sat_lon:sat_lon),
							sat_range,rv,findsun);

				if (mode=='p')
					quit=Print(sat, string,'p');

				if (mode=='v')
					quit=PrintVisible(sat, string);
			}

			if (mode=='p')
				quit=Print(sat, "\n",'p');

			if (mode=='v')
				quit=PrintVisible(sat, "\n");

			/* Move to next orbit */
			daynum=NextAOS(sat, &tle);

		}  while (!quit && !breakout && AosHappens(sat) && !Decayed(sat,daynum));
	}

	else
	{
		bkgdset(COLOR_PAIR(5)|A_BOLD);
		clear();

		if (!AosHappens(sat) || Decayed(sat,daynum))
			mvprintw(12,5,"*** Passes for %s cannot occur for your ground station! ***\n",sat->name);

		if (Geostationary(sat))
			mvprintw(12,3,"*** Orbital predictions cannot be made for a geostationary satellite! ***\n");

		beep();
		bkgdset(COLOR_PAIR(7)|A_BOLD);
		AnyKey();
		refresh();
	}
}

void PredictMoon(void)
{
	/* This function predicts "passes" of the Moon */

	int iaz, iel, lastel=0;
	char string[80], quit=0;
	double daynum, lastdaynum, moonrise=0.0;

	daynum=GetStartTime("the Moon");
	clear();

	if (xterm)
		fprintf(stderr,"\033]0;PREDICT: %s's Orbit Calendar for the Moon\007",qth.callsign);

	do
	{
		/* Determine moonrise */

		FindMoon(daynum);

		while (moonrise==0.0)
		{
			if (fabs(moon_el)<0.03)
				moonrise=daynum;
			else
			{
				daynum-=(0.004*moon_el);
				FindMoon(daynum);
			}
		}

		FindMoon(moonrise);
		daynum=moonrise;
		iaz=(int)rint(moon_az);
		iel=(int)rint(moon_el);

		do
		{
			/* Display pass of the moon */

			sprintf(string,"      %s%4d %4d  %5.1f  %5.1f  %5.1f  %6.1f%7.3f\n",Daynum2String(daynum), iel, iaz, moon_ra, moon_dec, moon_gha, moon_dv, moon_dx);
			quit=Print(NULL, string,'m');
			lastel=iel;
			lastdaynum=daynum;
			daynum+=0.04*(cos(deg2rad*(moon_el+0.5)));
			FindMoon(daynum);
			iaz=(int)rint(moon_az);
			iel=(int)rint(moon_el);

		} while (iel>3 && quit==0);

		while (lastel!=0 && quit==0)
		{
			daynum=lastdaynum;

			do
			{
				/* Determine setting time */

				daynum+=0.004*(sin(deg2rad*(moon_el+0.5)));
				FindMoon(daynum);
				iaz=(int)rint(moon_az);
				iel=(int)rint(moon_el);

			} while (iel>0);

			/* Print moonset */

			sprintf(string,"      %s%4d %4d  %5.1f  %5.1f  %5.1f  %6.1f%7.3f\n",Daynum2String(daynum), iel, iaz, moon_ra, moon_dec, moon_gha, moon_dv, moon_dx);
			quit=Print(NULL, string,'m');
			lastel=iel;
		}

		quit=Print(NULL, "\n",'m');
		daynum+=0.4;
		moonrise=0.0;

	} while (quit==0);
}

void PredictSun(void)
{
	/* This function predicts "passes" of the Sun. */

	int iaz, iel, lastel=0;
	char string[80], quit=0;
	double daynum, lastdaynum, sunrise=0.0;

	daynum=GetStartTime("the Sun");
	clear();

	if (xterm)
		fprintf(stderr,"\033]0;PREDICT: %s's Orbit Calendar for the Sun\007",qth.callsign);

	do
	{
		/* Determine sunrise */

		FindSun(daynum);

		while (sunrise==0.0)
		{
			if (fabs(sun_ele)<0.03)
				sunrise=daynum;
			else
			{
				daynum-=(0.004*sun_ele);
				FindSun(daynum);
			}
		}

		FindSun(sunrise);
		daynum=sunrise;
		iaz=(int)rint(sun_azi);
		iel=(int)rint(sun_ele);

		/* Print time of sunrise */

		do
		{
			/* Display pass of the sun */

			sprintf(string,"      %s%4d %4d  %5.1f  %5.1f  %5.1f  %6.1f%7.3f\n",Daynum2String(daynum), iel, iaz, sun_ra, sun_dec, sun_lon, sun_range_rate, sun_range);
			quit=Print(NULL, string,'o');
			lastel=iel;
			lastdaynum=daynum;
			daynum+=0.04*(cos(deg2rad*(sun_ele+0.5)));
			FindSun(daynum);
			iaz=(int)rint(sun_azi);
			iel=(int)rint(sun_ele);

		} while (iel>3 && quit==0);

		while (lastel!=0 && quit==0)
		{
			daynum=lastdaynum;

			do
			{
				/* Find sun set */

				daynum+=0.004*(sin(deg2rad*(sun_ele+0.5)));
				FindSun(daynum);
				iaz=(int)rint(sun_azi);
				iel=(int)rint(sun_ele);

			} while (iel>0);

			/* Print time of sunset */

			sprintf(string,"      %s%4d %4d  %5.1f  %5.1f  %5.1f  %6.1f%7.3f\n",Daynum2String(daynum), iel, iaz, sun_ra, sun_dec, sun_lon, sun_range_rate, sun_range);
			quit=Print(NULL, string,'o');
			lastel=iel;
		}

		quit=Print(NULL, "\n",'o');
		daynum+=0.4;
		sunrise=0.0;

	} while (quit==0);
}

char KbEdit(int x, int y)
{
	/* This function is used when editing QTH
	   and orbital data via the keyboard. */

	char need2save=0, input[25];

	echo();
	move(y-1,x-1);
	wgetnstr(stdscr,input,24);

	if (strlen(input)!=0)
	{
		need2save=1;  /* Save new data to variables */
		resave=1;     /* Save new data to disk files */
		strncpy(temp,input,24);
	}

	mvprintw(y-1,x-1,"%-25s",temp);

	refresh();
	noecho();

	return need2save;
}

void ShowOrbitData(void)
{
	/* This function permits displays a satellite's orbital
	   data.  The age of the satellite data is also provided. */

	int c, age;
	size_t namelength;
	double an_period, no_period, sma, c1, e2, satepoch;
	char days[5];

	for (;;)
	{
		struct sat_st * sat=Select();
		if (!sat)
			break;
		if (sat->meanmo==0.0)
			continue;

		bkgdset(COLOR_PAIR(2)|A_BOLD);
		clear();
		sma=331.25*exp(log(1440.0/sat->meanmo)*(2.0/3.0));
		an_period=1440.0/sat->meanmo;
		c1=cos(sat->incl*deg2rad);
		e2=1.0-(sat->eccn*sat->eccn);
		no_period=(an_period*360.0)/(360.0+(4.97*pow((xkmper/sma),3.5)*((5.0*c1*c1)-1.0)/(e2*e2))/sat->meanmo);
		satepoch=DayNum(1,0,sat->year)+sat->refepoch;
		age=(int)rint(CurrentDaynum()-satepoch);

		if (age==1)
			strcpy(days,"day");
		else
			strcpy(days,"days");

		namelength=strlen(sat->name);

		printw("\n");

		for (c=41; c>namelength; c-=2)
			printw(" ");

		printw("Orbital Data For %s / Catalog Number %ld\n",sat->name,sat->catnum);
		attrset(COLOR_PAIR(3)|A_BOLD);
		printw("\n\t\t  Issued %d %s ago on %s UTC\n\n",age,days,Daynum2String(satepoch));

		attrset(COLOR_PAIR(4)|A_BOLD);
		mvprintw(5,21,"Reference Epoch");
		mvprintw(6,21,"Inclination");
		mvprintw(7,21,"RAAN");
		mvprintw(8,21,"Eccentricity");
		mvprintw(9,21,"Arg of Perigee");
		mvprintw(10,21,"Mean Anomaly");
		mvprintw(11,21,"Mean Motion");
		mvprintw(12,21,"Decay Rate");
		mvprintw(13,21,"Nddot/6 Drag");
		mvprintw(14,21,"Bstar Drag Factor");
		mvprintw(15,21,"Semi-Major Axis");
		mvprintw(16,21,"Apogee Altitude");
		mvprintw(17,21,"Perigee Altitude");
		mvprintw(18,21,"Anomalistic Period");
		mvprintw(19,21,"Nodal Period");
		mvprintw(20,21,"Orbit Number");
		mvprintw(21,21,"Element Set Number");

		attrset(COLOR_PAIR(2)|A_BOLD);
		mvprintw(5,40,": %02d %.8f",sat->year,sat->refepoch);
		mvprintw(6,40,": %.4f deg",sat->incl);
		mvprintw(7,40,": %.4f deg",sat->raan);
		mvprintw(8,40,": %g",sat->eccn);
		mvprintw(9,40,": %.4f deg",sat->argper);
		mvprintw(10,40,": %.4f deg",sat->meanan);
		mvprintw(11,40,": %.8f rev/day",sat->meanmo);
		mvprintw(12,40,": %g rev/day/day",sat->drag);
		mvprintw(13,40,": %g rev/day/day/day",sat->nddot6);
		mvprintw(14,40,": %g 1/earth radii",sat->bstar);
		mvprintw(15,40,": %.4f km",sma);
		mvprintw(16,40,": %.4f km",sma*(1.0+sat->eccn)-xkmper);
		mvprintw(17,40,": %.4f km",sma*(1.0-sat->eccn)-xkmper);
		mvprintw(18,40,": %.4f mins",an_period);
		mvprintw(19,40,": %.4f mins",no_period);
		mvprintw(20,40,": %ld",sat->orbitnum);
		mvprintw(21,40,": %ld",sat->setnum);

		attrset(COLOR_PAIR(3)|A_BOLD);
		refresh();
		AnyKey();
	}
}	

void KepEdit(void)
{
	/* This function permits keyboard editing of the orbital database. */

	for (;;)
	{
		struct sat_st *sat=Select();
		
		if (!sat)
			break;

		bkgdset(COLOR_PAIR(3)|A_BOLD);
		clear();
		mvprintw(3,1,"\t\t   *  Orbital Database Editing Utility  *\n\n\n");
		attrset(COLOR_PAIR(4)|A_BOLD);

		printw("\n\t\t\tSpacecraft Name :");
		printw("\n\t\t\tCatalog Number  :");
		printw("\n\t\t\tDesignator      :");
		printw("\n\t\t\tReference Epoch :");
		printw("\n\t\t\tInclination     :");
		printw("\n\t\t\tRAAN            :");
		printw("\n\t\t\tEccentricity    :");
		printw("\n\t\t\tArg of Perigee  :");
		printw("\n\t\t\tMean Anomaly    :");
		printw("\n\t\t\tMean Motion     :");
		printw("\n\t\t\tDecay Rate      :");
		printw("\n\t\t\tNddot/6         :");
		printw("\n\t\t\tBstar Drag Term :");
		printw("\n\t\t\tOrbit Number    :");
		printw("\n\t\t\tElement Set No. :");

		attrset(COLOR_PAIR(2)|A_BOLD);

		mvprintw(7,42,"%s",sat->name);
		mvprintw(8,42,"%ld",sat->catnum);
		mvprintw(9,42,"%s",sat->designator);
		mvprintw(10,42,"%02d %.8f",sat->year,sat->refepoch);
		mvprintw(11,42,"%.4f",sat->incl);
		mvprintw(12,42,"%.4f",sat->raan);
		mvprintw(13,42,"%g",sat->eccn);
		mvprintw(14,42,"%.4f",sat->argper);
		mvprintw(15,42,"%.4f",sat->meanan);
		mvprintw(16,42,"%.8f",sat->meanmo);
		mvprintw(17,42,"%g",sat->drag);
		mvprintw(18,42,"%g",sat->nddot6);
		mvprintw(19,42,"%g",sat->bstar);
		mvprintw(20,42,"%ld",sat->orbitnum);
		mvprintw(21,42,"%ld",sat->setnum);

		curs_set(1);
		refresh();

		sprintf(temp,"%s",sat->name);

		if (KbEdit(43,8))
			strncpy(sat->name,temp,24);

		sprintf(temp,"%ld",sat->catnum);

		if (KbEdit(43,9))
			sscanf(temp,"%ld",&sat->catnum);

		sprintf(temp,"%s",sat->designator);

		if (KbEdit(43,10))
			sscanf(temp,"%s",sat->designator);

		sprintf(temp,"%02d %4.8f",sat->year,sat->refepoch);

		if (KbEdit(43,11))
			sscanf(temp,"%d %lf",&sat->year,&sat->refepoch);

		sprintf(temp,"%4.4f",sat->incl);

		if (KbEdit(43,12))
			sscanf(temp,"%lf",&sat->incl);

		sprintf(temp,"%4.4f",sat->raan);

		if (KbEdit(43,13))
			sscanf(temp,"%lf",&sat->raan);

		sprintf(temp,"%g",sat->eccn);

		if (KbEdit(43,14))
			sscanf(temp,"%lf",&sat->eccn);

		sprintf(temp,"%4.4f",sat->argper);

		if (KbEdit(43,15))
			sscanf(temp,"%lf",&sat->argper);

		sprintf(temp,"%4.4f",sat->meanan);

		if (KbEdit(43,16))
			sscanf(temp,"%lf",&sat->meanan);

		sprintf(temp,"%4.8f",sat->meanmo);

		if (KbEdit(43,17))
			sscanf(temp,"%lf",&sat->meanmo);

		sprintf(temp,"%g",sat->drag);

		if (KbEdit(43,18))
			sscanf(temp,"%lf",&sat->drag);

		sprintf(temp,"%g",sat->nddot6);

		if (KbEdit(43,19))
			sscanf(temp,"%lf",&sat->nddot6);

		sprintf(temp,"%g",sat->bstar);

		if (KbEdit(43,20))
			sscanf(temp,"%lf",&sat->bstar);

		sprintf(temp,"%ld",sat->orbitnum);

		if (KbEdit(43,21))
			sscanf(temp,"%ld",&sat->orbitnum);

		sprintf(temp,"%ld",sat->setnum);

		if (KbEdit(43,22))
			sscanf(temp,"%ld",&sat->setnum);

		curs_set(0);
	}

	if (resave)
	{
		SaveTLE();
		resave=0;
	}
}	

void QthEdit(void)
{
	/* This function permits keyboard editing of
	   the ground station's location information. */

	bkgdset(COLOR_PAIR(3)|A_BOLD);
	clear();
	curs_set(1);
	mvprintw(7,0,"\t\t *  Ground Station Location Editing Utility  *\n\n\n");

	attrset(COLOR_PAIR(4)|A_BOLD);
	printw("\n\t\t\tStation Callsign  :");
	printw("\n\t\t\tStation Latitude  :");
	printw("\n\t\t\tStation Longitude :");
	printw("\n\t\t\tStation Altitude  :");

	attrset(COLOR_PAIR(2)|A_BOLD);
	mvprintw(11,44,"%s",qth.callsign);

	if (io_lat=='N')
		mvprintw(12,44,"%g [DegN]",+qth.stnlat);
	else
		mvprintw(12,44,"%g [DegS]",-qth.stnlat);

	if (io_lon=='W')
		mvprintw(13,44,"%g [DegW]",+qth.stnlong);
	else
		mvprintw(13,44,"%g [DegE]",-qth.stnlong);

	mvprintw(14,44,"%d [m]",qth.stnalt);

	refresh();

	sprintf(temp,"%s",qth.callsign);

	mvprintw(18,12,"Enter the callsign or identifier of your ground station");

	if (KbEdit(45,12))
		strncpy(qth.callsign,temp,16);

	if (io_lat=='N')
		sprintf(temp,"%g [DegN]",+qth.stnlat);
	else
		sprintf(temp,"%g [DegS]",-qth.stnlat);

	if (io_lat=='N')
		mvprintw(18,12,"Enter your latitude in degrees NORTH  (south=negative) ");
	else
		mvprintw(18,12,"Enter your latitude in degrees SOUTH  (north=negative) ");
 
	mvprintw(19,12,"  Decimal (74.2467) or DMS (74 14 48) format allowed");

	if (KbEdit(45,13))
	{
		if (io_lat=='N')
			qth.stnlat=+ReadBearing(temp);
		else
			qth.stnlat=-ReadBearing(temp);
	}
 
	if (io_lon=='W')
		sprintf(temp,"%g [DegW]",+qth.stnlong);
	else
		sprintf(temp,"%g [DegE]",-qth.stnlong);
 
	if (io_lon=='W')
		mvprintw(18,12,"Enter your longitude in degrees WEST   (east=negative) ");
	else
		mvprintw(18,12,"Enter your longitude in degrees EAST   (west=negative) ");
 
	if (KbEdit(45,14))
	{
		if (io_lon=='W')
			qth.stnlong=+ReadBearing(temp);
		else
			qth.stnlong=-ReadBearing(temp);
	}
 
	move(19,12);
	clrtoeol();
	mvprintw(18,12,"    Enter your altitude above sea level (in meters)   ");

	sprintf(temp,"%d",qth.stnalt);

	if (KbEdit(45,15))
		sscanf(temp,"%d",&qth.stnalt);

	if (resave)
	{
		SaveQTH();
		resave=0;
	}
}

void SingleTrack(struct sat_st *sat, char speak)
{
	/* This function tracks a single satellite in real-time
	   until 'Q' or ESC is pressed.  x represents the index
	   of the satellite being tracked.  If speak=='T', then
	   the speech routines are enabled. */

	int	xponder=0, polarity=0;
	double oldaz=0, oldel=0;
	char	approaching=0, command[80], aos_alarm=0,
		eclipse_alarm=0, old_visibility=0;
	double	oldtime=0.0, nextaos=0.0, lostime=0.0, aoslos=0.0,
		downlink=0.0, uplink=0.0, downlink_start=0.0,
	downlink_end=0.0, uplink_start=0.0, uplink_end=0.0;
		//shift;
	long	newtime, lasttime=0;

	tle_t tle;
	PreCalc(sat, &tle);

	const struct sat_db_st *comsat = (sat->db && sat->db->transponders>0) ? sat->db :  NULL;
	if (comsat)
	{
		downlink_start = comsat->transponder[xponder].downlink_start;
		downlink_end   = comsat->transponder[xponder].downlink_end;
		uplink_start   = comsat->transponder[xponder].uplink_start;
		uplink_end     = comsat->transponder[xponder].uplink_end;

		if (downlink_start>downlink_end)
			polarity=-1;

		if (downlink_start<downlink_end)
			polarity=1;

		if (downlink_start==downlink_end)
			polarity=0;

		downlink=0.5*(downlink_start+downlink_end);
		uplink=0.5*(uplink_start+uplink_end);
	}

	daynum=CurrentDaynum();
	const int aoshappens=AosHappens(sat);
	const int geostationary=Geostationary(sat);
	const int decayed=Decayed(sat,0.0);

	if (xterm)
		fprintf(stderr,"\033]0;PREDICT: Tracking %-10s\007",sat->name);
	halfdelay(2);
	curs_set(0);
	bkgdset(COLOR_PAIR(3));
	clear();

	attrset(COLOR_PAIR(6)|A_REVERSE|A_BOLD);

	mvprintw(0,0,"                                                                                 ");
	mvprintw(1,0,"                     PREDICT Real-Time Satellite Tracking                        ");
	mvprintw(2,0,"                 Tracking: %-10sOn                                          ",Abbreviate(sat->name,9));
	mvprintw(3,0,"                                                                                 ");

	attrset(COLOR_PAIR(4)|A_BOLD);

	const int tshift = comsat ? 0 : 2;
	const int bshift = comsat ? 0 : -2;

	mvprintw(5+tshift,1,"Satellite     Direction     Velocity     Footprint    Altitude     Slant Range");
	mvprintw(6+tshift,1,"---------     ---------     --------     ---------    --------     -----------");
	mvprintw(7+tshift,1,"        .            Az           mi            mi          mi              mi");
	mvprintw(8+tshift,1,"        .            El           km            km          km              km");
	if (comsat)
	{
		mvprintw(12,1,"Uplink   :");
		mvprintw(13,1,"Downlink :");
		mvprintw(14,1,"Delay    :");
		mvprintw(14,55,"Echo      :");
		mvprintw(13,29,"RX:");
		mvprintw(13,55,"Path loss :");
		mvprintw(12,29,"TX:");
		mvprintw(12,55,"Path loss :");
	}
	mvprintw(16+bshift,1,"Eclipse Depth   Orbital Phase   Orbital Model   Squint Angle      AutoTracking");
	mvprintw(17+bshift,1,"-------------   -------------   -------------   ------------      ------------");

	for (;;)
	{
		attrset(COLOR_PAIR(6)|A_REVERSE|A_BOLD);
		daynum=CurrentDaynum();
		mvprintw(2,41,"%s",Daynum2String(daynum));
		attrset(COLOR_PAIR(2)|A_BOLD);
		Calc(&tle);

		attrset(COLOR_PAIR(4)|A_BOLD);
		mvprintw(7+tshift,8,(io_lat=='N'?"N":"S"));
		mvprintw(8+tshift,8,(io_lon=='W'?"W":"E"));

		sat_footprint=12756.33*acos(xkmper/(xkmper+sat_alt));

		attrset(COLOR_PAIR(2)|A_BOLD);

		mvprintw(7+tshift,1,"%6.2f",(io_lat=='N'?+1:-1)*sat_lat);
		mvprintw(8+tshift,1,"%6.2f",(io_lon=='W'?360.0-sat_lon:sat_lon));

		mvprintw(7+tshift,15,"%6.2f", sat_azi);
		mvprintw(8+tshift,15,"%+6.2f",sat_ele);

		mvprintw(7+tshift,29,"%5.0f ",(3600.0*sat_vel)*km2mi);
		mvprintw(8+tshift,29,"%5.0f ",3600.0*sat_vel);

		mvprintw(7+tshift,42,"%6.0f ",sat_footprint*km2mi);
		mvprintw(8+tshift,42,"%6.0f ",sat_footprint);

		mvprintw(7+tshift,55,"%5.0f ",sat_alt*km2mi);
		mvprintw(8+tshift,55,"%5.0f ",sat_alt);

		mvprintw(7+tshift,68,"%8.0f",sat_range*km2mi);
		mvprintw(8+tshift,68,"%8.0f",sat_range);

		mvprintw(18+bshift,3,"%+6.2f%c  ",eclipse_depth/deg2rad,176);
		mvprintw(18+bshift,20,"%5.1f",256.0*(phase/twopi));
		mvprintw(18+bshift,37,"%s",isFlagSet(DEEP_SPACE_EPHEM_FLAG) ? "SDP4" : "SGP4");

		char visibility;
		if (!sat_sun_status)
			visibility='N';
		else if (sun_ele<=-12.0 && sat_ele>=0.0)
			visibility='V';
		else
			visibility='D';
		server_data[sat-&sat[0]].visibility = visibility;

		if (comsat)
		{
			if (downlink!=0.0)
				mvprintw(13,11,"%11.5f MHz",downlink);

			else
				mvprintw(13,11,"               ");

			if (uplink!=0.0)
				mvprintw(12,11,"%11.5f MHz",uplink);

			else
				mvprintw(12,11,"               ");
		}

		if (antfd!=-1)
		{
			if (sat_ele>=0.0)
				mvprintw(18+bshift,67,"   Active   ");
			else
				mvprintw(18+bshift,67,"Standing  By");
		}
		else
			mvprintw(18+bshift,67,"Not  Enabled");

		if (calc_squint)
			mvprintw(18+bshift,52,"%+6.2f",squint);
		else
			mvprintw(18+bshift,54,"N/A");

		const double doppler100=-100.0e06*((sat_range_rate*1000.0)/299792458.0);

		if (comsat)
		{
			attrset(COLOR_PAIR(4)|A_BOLD);
			const size_t length=strlen(comsat->transponder[xponder].name)/2;
			mvprintw(10,(int)(40-length),"%s",comsat->transponder[xponder].name);
		}

		if (sat_ele>=0.0)
		{
			if (!aos_alarm)
			{
				beep();
				aos_alarm=1;
			}

			if (comsat)
			{
				attrset(COLOR_PAIR(4)|A_BOLD);

				mvprintw(14,34, (sat_range_rate <= -0.1) ? "Approaching" :
				                (sat_range_rate >=  0.1) ? "  Receding " :
				                                           "    TCA    ");

				attrset(COLOR_PAIR(2)|A_BOLD);

				const double delay=1000.0*((1000.0*sat_range)/299792458.0);

				if (downlink!=0.0)
				{
					const double dopp=1.0e-08*(doppler100*downlink);
					mvprintw(13,32,"%11.5f MHz",downlink+dopp);
					const double loss=32.4+(20.0*log10(downlink))+(20.0*log10(sat_range));
					mvprintw(13,67,"%7.3f dB",loss);
					mvprintw(14,13,"%7.3f   ms",delay);
				}

				else
				{
					mvprintw(13,32,"                ");
					mvprintw(13,67,"          ");
					mvprintw(14,13,"            ");
				}

				if (uplink!=0.0)
				{
					const double dopp=1.0e-08*(doppler100*uplink);
					mvprintw(12,32,"%11.5f MHz",uplink-dopp);
					const double loss=32.4+(20.0*log10(uplink))+(20.0*log10(sat_range));
					mvprintw(12,67,"%7.3f dB",loss);
				}

				else
				{
					mvprintw(12,32,"                ");
					mvprintw(12,67,"          ");
				}

				if (uplink!=0.0 && downlink!=0.0)
					mvprintw(14,67,"%7.3f ms",2.0*delay);

				else
					mvprintw(14,67,"              ");
			}

			if (speak=='T' && soundcard)
			{
				if (!eclipse_alarm && fabs(eclipse_depth)<0.015) /* ~1 deg */
				{
					/* Hold off regular announcements if
					   satellite is within about 2 degrees
					   of entering into or out of an
					   eclipse. */

					oldtime=CurrentDaynum();

					if ((old_visibility=='V' || old_visibility=='D') && visibility=='N')
					{
						sprintf(command,"%svocalizer/vocalizer eclipse &",predictpath);
						system(command);
						eclipse_alarm=1;
						oldtime-=0.000015*sqrt(sat_alt);
					}

					if (old_visibility=='N' && (visibility=='V' || visibility=='D'))
					{
						sprintf(command,"%svocalizer/vocalizer sunlight &",predictpath);
						system(command);
						eclipse_alarm=1;
						oldtime-=0.000015*sqrt(sat_alt);
					}
				}

				if ((CurrentDaynum()-oldtime)>(0.00003*sqrt(sat_alt)))
				{
					if (sat_range_rate<0.0)
						approaching='+';

					if (sat_range_rate>0.0)
						approaching='-';

					sprintf(command,"%svocalizer/vocalizer %.0f %.0f %c %c &",
							predictpath,sat_azi,sat_ele,approaching,visibility);
					system(command);
  					oldtime=CurrentDaynum();
					old_visibility=visibility;
				}

				if (sat_ele<=1.0 && approaching=='-')
				{
					/* Suspend regular announcements
					   as we approach LOS. */

					oldtime=CurrentDaynum();
				}
			}
		}
		else
		{
			lostime=0.0;
			aos_alarm=0;
			eclipse_alarm=0;

			if (comsat)
			{
				mvprintw(12,32,"                ");
				mvprintw(12,67,"          ");
				mvprintw(13,32,"                ");
				mvprintw(13,67,"          ");
				mvprintw(14,13,"            ");
				mvprintw(14,34,"           ");
				mvprintw(14,67,"          ");
			}
		}

		attrset(COLOR_PAIR(3)|A_BOLD);

		mvprintw(21,22,"Orbit Number: %ld",rv);

		/* Send data to serial port antenna tracker
		   either as needed (when it changes), or
		   once per second. */

		if (sat_ele>=0.0 && antfd!=-1)
		{
			newtime=(long)time(NULL);

			if ( fabs(oldel-sat_ele)>=1.0 || fabs(oldaz-sat_azi)>=1.0 || (once_per_second && newtime>lasttime))
			{
				TrackDataOut(antfd,sat_ele,sat_azi);
				oldel=sat_ele;
				oldaz=sat_azi;
				lasttime=newtime;
			}
		}

		mvprintw(23,22,"Spacecraft is currently ");

		if (visibility=='V')
			mvprintw(23,46,"visible    ");

		if (visibility=='D')
			mvprintw(23,46,"in sunlight");

		if (visibility=='N')
			mvprintw(23,46,"in eclipse ");

		attrset(COLOR_PAIR(4)|A_BOLD);
		mvprintw(20,5,"   Sun   ");
		mvprintw(21,5,"---------");
		attrset(COLOR_PAIR(3)|A_BOLD);
		mvprintw(22,5,"%-7.2fAz",sun_azi);
		mvprintw(23,4,"%+-6.2f  El",sun_ele);

		FindMoon(daynum);

		attrset(COLOR_PAIR(4)|A_BOLD);
		mvprintw(20,65,"  Moon  ");
		mvprintw(21,65,"---------");
		attrset(COLOR_PAIR(3)|A_BOLD);
		mvprintw(22,65,"%-7.2fAz",moon_az);
		mvprintw(23,64,"%+-6.2f  El",moon_el);

		if (geostationary)
		{
			if (sat_ele>=0.0)
				mvprintw(22,22,"Satellite orbit is geostationary");
			else
				mvprintw(22,22,"This satellite never reaches AOS");
			aoslos=-3651.0;
		}

		else if (!aoshappens || decayed)
		{
			mvprintw(22,22,"This satellite never reaches AOS");
			aoslos=-3651.0;
		}

		else if (sat_ele>=0.0 && daynum>lostime)
		{
			lostime=FindLOS2(sat, &tle);
			mvprintw(22,22,"LOS at: %s UTC  ",Daynum2String(lostime));
			aoslos=lostime;
		}

		else if (sat_ele<0.0 && aoshappens && daynum>aoslos)
		{
			daynum+=0.003;  /* Move ahead slightly... */
			nextaos=FindAOS(sat, &tle);
			mvprintw(22,22,"Next AOS: %s UTC",Daynum2String(nextaos));
			aoslos=nextaos;

			if (oldtime!=0.0 && speak=='T' && soundcard)
			{
				/* Announce LOS */

				sprintf(command,"%svocalizer/vocalizer los &",predictpath);
				system(command);
			}
		}

		/* This is where the variables for the socket server are updated. */

		if (socket_flag)
		{
			struct server_data_st *const server  = &server_data[sat-&sat[0]];
			server->az            = sat_azi;
			server->el            = sat_ele;
			server->lattitude     = sat_lat;
			server->longitude     = 360.0-sat_lon;
			server->footprint     = sat_footprint;
			server->range         = sat_range;
			server->altitude      = sat_alt;
			server->velocity      = sat_vel;
			server->orbitnum      = rv;
			server->doppler       = doppler100;
			server->nextevent     = aoslos;
			server->eclipse_depth = eclipse_depth/deg2rad;
			server->phase         = 360.0*(phase/twopi);
			server->squint        = calc_squint ? squint : 360.0;

			FindSun(daynum);

			strcpy(tracking_mode, sat->name);
		}

		/* Get input from keyboard */

		const int ans=tolower(getch());
		if (ans=='q' || ans==27)
			break;

		/* We can force PREDICT to speak by pressing 'T' */

		if (ans=='t')
			oldtime=0.0;

		/* If we receive a RELOAD_TLE command through the
		   socket connection or an 'r' through the keyboard,
		   reload the TLE file.  */

		if (reload_tle || ans=='r')
		{
			ReadDataFiles();
			reload_tle=0;
		}

		if (comsat)
		{
			if (ans==' ' && comsat->transponders>1)
			{
				xponder++;

				if (xponder>=comsat->transponders)
					xponder=0;

				move(10,1);
				clrtoeol();

				downlink_start = comsat->transponder[xponder].downlink_start;
				downlink_end   = comsat->transponder[xponder].downlink_end;
				uplink_start   = comsat->transponder[xponder].uplink_start;
				uplink_end     = comsat->transponder[xponder].uplink_end;

				if (downlink_start>downlink_end)
					polarity=-1;

				if (downlink_start<downlink_end)
					polarity=1;

				if (downlink_start==downlink_end)
					polarity=0;

				downlink=0.5*(downlink_start+downlink_end);
				uplink=0.5*(uplink_start+uplink_end);
			}

			if (ans=='>' || ans=='.')	/* Raise uplink frequency */
			{
				const double shift = (ans=='>') ? 0.001 /* 1 kHz */ : 0.0001 /* 100 Hz */;
				uplink   += shift*(double)abs(polarity);
				downlink += shift*(double)polarity;

				if (uplink>uplink_end)
				{
					uplink=uplink_start;
					downlink=downlink_start;
				}
			}

			if (ans=='<' || ans== ',')	/* Lower uplink frequency */
			{
				const double shift = (ans=='<') ? 0.001 /* 1 kHz */ : 0.0001 /* 100 Hz */;
				uplink   -= shift*(double)abs(polarity);
				downlink -= shift*(double)polarity;

				if (uplink<uplink_start)
				{
					uplink=uplink_end;
					downlink=downlink_end;
				}
			}
		}

		refresh();

		halfdelay(2);
	}

	cbreak();
	strcpy(tracking_mode, "NONE\n");
}

void MultiTrack(void)
{
	/* This function tracks all satellites in the program's
	   database simultaneously until 'Q' or ESC is pressed.
	   Satellites in range are HIGHLIGHTED.  Coordinates
	   for the Sun and Moon are also displayed. */

	struct data_st
	{
		struct sat_st *sat;
		int inrange;
		int ok2predict;
		double aos;
		double aos2;
		double los;
		double aoslos;
	} data_table[24];

	unsigned char sunstat=0;

	double nextcalctime=0.0;

	if (xterm)
		fprintf(stderr,"\033]0;PREDICT: Multi-Satellite Tracking Mode\007");

	curs_set(0);
	attrset(COLOR_PAIR(6)|A_REVERSE|A_BOLD);
	clear();

	mvprintw(0,0,"                                                                                ");
	mvprintw(1,0,"                     PREDICT Real-Time Multi-Tracking Mode                      ");
	mvprintw(2,0,"                    Current Date/Time:                                          ");
	mvprintw(3,0,"                                                                                ");

	attrset(COLOR_PAIR(2)|A_REVERSE);

	mvprintw(4,0," Satellite  Az   El %s  %s  Range  | Satellite  Az   El %s  %s  Range   ",(io_lat=='N'?"LatN":"LatS"),(io_lon=='W'?"LonW":"LonE"),(io_lat=='N'?"LatN":"LatS"),(io_lon=='W'?"LonW":"LonE"));

	for (int x=0; x<24; x++)
	{
		if (!Geostationary(&sat_table[x]) && AosHappens(&sat_table[x]) && !Decayed(&sat_table[x],0.0))
			data_table[x].ok2predict=1;
		else
			data_table[x].ok2predict=0;

		data_table[x].aoslos=0.0;
		data_table[x].los=0.0;
		data_table[x].aos=0.0;
		data_table[x].aos2=0.0;
		data_table[x].sat = &sat_table[x];
	}

	for (;;)
	{
		for (int z=0; z<24; z++)
		{
			const int y=z/2;

			int sat_index;
			int x;
			if (z%2)
			{
				sat_index=y+12;
				x=41;
			}
			else
			{
				sat_index=y;
				x=1;
			}
			struct sat_st  *sat  = &sat_table[sat_index];
			struct data_st *data = &data_table[sat_index];

			if (sat->meanmo!=0.0 && Decayed(sat,0.0)!=1)
			{
				daynum=CurrentDaynum();
				tle_t tle;
				PreCalc(sat, &tle);
				Calc(&tle);

				if (sat_ele>=0.0)
				{
					attrset(COLOR_PAIR(2)|A_BOLD);
					data->inrange=1;
				}

				else
				{
					attrset(COLOR_PAIR(2));
					data->inrange=0;
				}

				if (sat_sun_status)
				{
					if (sun_ele<=-12.0 && sat_ele>=0.0)
						sunstat='V';
					else
						sunstat='D';
				}

				else
					sunstat='N';

				mvprintw(y+6,x,"%-10s%3.0f  %+3.0f  %3.0f   %3.0f %6.0f %c", Abbreviate(sat->name,9),sat_azi,sat_ele,(io_lat=='N'?+1:-1)*sat_lat,(io_lon=='W'?360.0-sat_lon:sat_lon),sat_range,sunstat);

				if (socket_flag)
				{
					struct server_data_st *const server  = &server_data[sat_index];
					server->az            = sat_azi;
					server->el            = sat_ele;
					server->lattitude     = sat_lat;
					server->longitude     = 360.0-sat_lon;
					server->footprint     = sat_footprint;
					server->range         = sat_range;
					server->altitude      = sat_alt;
					server->velocity      = sat_vel;
					server->orbitnum      = rv;
					server->visibility    = sunstat;
					server->eclipse_depth = eclipse_depth/deg2rad;
					server->phase         = 360.0*(phase/twopi);

					server->doppler       = -100e06*((sat_range_rate*1000.0)/299792458.0);

					if (calc_squint)
						server->squint=squint;
					else
						server->squint=360.0;

					FindSun(daynum);
					strcpy(tracking_mode,"MULTI\n");
				}

				attrset(COLOR_PAIR(4)|A_BOLD);
				mvprintw(20,5,"   Sun   ");
				mvprintw(21,5,"---------");
				attrset(COLOR_PAIR(3)|A_BOLD);
				mvprintw(22,5,"%-7.2fAz",sun_azi);
				mvprintw(23,4,"%+-6.2f  El",sun_ele);

				FindMoon(daynum);

				attrset(COLOR_PAIR(4)|A_BOLD);
				mvprintw(20,65,"  Moon  ");
				mvprintw(21,65,"---------");
				attrset(COLOR_PAIR(3)|A_BOLD);
				mvprintw(22,65,"%-7.2fAz",moon_az);
				mvprintw(23,64,"%+-6.2f  El",moon_el);

				/* Calculate Next Event (AOS/LOS) Times */

				if (data->ok2predict && daynum>data->los && data->inrange)
					data->los=FindLOS2(sat, &tle);

				if (data->ok2predict && daynum>data->aos)
				{
					if (data->inrange)
						data->aos=NextAOS(sat, &tle);
					else
						data->aos=FindAOS(sat, &tle);
				}

				if (data->inrange)
					data->aoslos=data->los;
				else
					data->aoslos=data->aos;

				if (socket_flag)
				{
					if (data->ok2predict)
						server_data[sat_index].nextevent=data->aoslos;

					else
						server_data[sat_index].nextevent=-3651.0;
				}

				data->aos2=data->aos;
			}

			if (Decayed(sat,0.0))
			{
				attrset(COLOR_PAIR(2));
				mvprintw(y+6,x,"%-10s---------- Decayed ---------", Abbreviate(sat->name,9));

				if (socket_flag)
				{
					struct server_data_st *const server  = &server_data[sat_index];
					server->az            = 0.0;
					server->el            = 0.0;
					server->lattitude     = 0.0;
					server->longitude     = 0.0;
					server->footprint     = 0.0;
					server->range         = 0.0;
					server->altitude      = 0.0;
					server->velocity      = 0.0;
					server->orbitnum      = 0L;
					server->visibility    = 'N';
					server->eclipse_depth = 0.0;
					server->phase         = 0.0;
					server->doppler       = 0.0;
					server->squint        = 0.0;
					server->nextevent     = -3651.0;
				}
			}
 		}

		attrset(COLOR_PAIR(6)|A_REVERSE|A_BOLD);

		daynum=CurrentDaynum();
		mvprintw(2,39,"%s",Daynum2String(daynum));

		if (daynum>nextcalctime)
		{
			/* Bubble sort the AOS times */

			for (int z=22; z>=0; z--)
				for (int y=0; y<=z; y++)
					if (data_table[y].aos2 >= data_table[y+1].aos2)
					{
						const struct data_st tmp = data_table[y];
						data_table[y] = data_table[y+1];
						data_table[y+1] = tmp;
					}

			/* Display list of upcoming passes */

			attrset(COLOR_PAIR(4)|A_BOLD);
			mvprintw(19,31,"Upcoming Passes");
			mvprintw(20,31,"---------------");
			attrset(COLOR_PAIR(3)|A_BOLD);

			int z=-1;
			for (int x=0, y=0; x<21 && y!=3; x++)
			{
				if (data_table[x].ok2predict && data_table[x].aos2!=0.0)
				{
					mvprintw(y+21,19,"%10s on %s UTC",Abbreviate(data_table[x].sat->name,9),Daynum2String(data_table[x].aos2));

					if (z==-1)
						z=x;
					y++;
				}
			}

			if (z!=-1)
				nextcalctime=data_table[z].aos2;
		}

		refresh();
		halfdelay(2);  /* Increase if CPU load is too high */
		int ans=tolower(getch());
		if (ans=='q' || ans==27)
			break;

		/* If we receive a RELOAD_TLE command through the
		   socket connection, or an 'r' through the keyboard,
		   reload the TLE file.  */

		if (reload_tle || ans=='r')
		{
			ReadDataFiles();
			reload_tle=0;
			nextcalctime=0.0;
		}
	}

	cbreak();
	strcpy(tracking_mode, "NONE\n");
}

void Illumination(struct sat_st *sat)
{
	double startday, oneminute, sunpercent;
	int eclipses, minutes, quit, breakout=0;
	char string1[365], string[725], datestring[25], count;

	oneminute=1.0/(24.0*60.0);

	tle_t tle;
	PreCalc(sat, &tle);
	daynum=floor(GetStartTime(sat->name));
	startday=daynum;
	count=0;

	curs_set(0);
	clear();

	if (xterm)
		fprintf(stderr,"\033]0;PREDICT: %s's Solar Illumination Calendar For %s\007",qth.callsign, sat->name);

	do
	{
		attrset(COLOR_PAIR(4));
		mvprintw(23,6,"                 Calculating... Press [ESC] To Quit");
		refresh();

		count++;
		daynum=startday;

		for (minutes=0, eclipses=0; minutes<1440; minutes++)
		{
			Calc(&tle);

			if (sat_sun_status==0)
				eclipses++;

			daynum=startday+(oneminute*(double)minutes);
		}

		sunpercent=((double)eclipses)/((double)minutes);
		sunpercent=100.0-(sunpercent*100.0);

		strcpy(datestring,Daynum2String(startday));
		datestring[11]=0;
		sprintf(string1,"      %s    %4d    %6.2f%c",datestring,1440-eclipses,sunpercent,37);

		/* Allow a quick way out */

		nodelay(stdscr,TRUE);

		if (getch()==27)
			breakout=1;

		nodelay(stdscr,FALSE);

		startday+=18.0;

		daynum=startday;

		for (minutes=0, eclipses=0; minutes<1440; minutes++)
		{
			Calc(&tle);

			if (sat_sun_status==0)
				eclipses++;

			daynum=startday+(oneminute*(double)minutes);
		}

		sunpercent=((double)eclipses)/((double)minutes);
		sunpercent=100.0-(sunpercent*100.0);

		strcpy(datestring,Daynum2String(startday));
		datestring[11]=0;
		sprintf(string,"%s\t %s    %4d    %6.2f%c\n",string1,datestring,1440-eclipses,sunpercent,37);
		quit=Print(sat, string,'s');

		/* Allow a quick way out */

		nodelay(stdscr,TRUE);

		if (getch()==27)
			breakout=1;

		nodelay(stdscr,FALSE);

		if (count<18)
			startday-=17.0;
		else
		{
			count=0;
			startday+=1.0;
		}
	}
	while (quit!=1 && breakout!=1 && Decayed(sat,daynum)==0);
}

void MainMenu(void)
{
	/* Start-up menu.  Your wish is my command. :-) */

	Banner();
	attrset(COLOR_PAIR(4)|A_BOLD);
	mvprintw(10,28,"--==[ Main Menu ]==--");

	attrset(COLOR_PAIR(3)|A_BOLD);
	mvprintw(13,1,"[P]: Predict Satellite Passes");
	mvprintw(14,1,"[V]: Predict Visible Passes");
	mvprintw(15,1,"[S]: Solar Illumination Predictions");
	mvprintw(16,1,"[L]: Lunar Predictions");
	mvprintw(17,1,"[O]: Solar Predictions");
	mvprintw(18,1,"[T]: Single Satellite Tracking Mode");
	mvprintw(19,1,"[M]: Multi-Satellite Tracking Mode");

	mvprintw(13,40,"[I]: Program Information");
	mvprintw(14,40,"[G]: Edit Ground Station Information");
	mvprintw(15,40,"[D]: Display Satellite Orbital Data");
	mvprintw(16,40,"[U]: Update Sat Elements From File");
	mvprintw(17,40,"[E]: Manually Edit Orbital Elements");
	mvprintw(18,40,"[B]: Edit Transponder Database");
	mvprintw(19,40,"[Q]: Exit PREDICT");

	if (socket_flag)
	{
		attrset(COLOR_PAIR(4)|A_BOLD);
		mvprintw(22,33,"Server Mode");
	}

	refresh();

	if (xterm)
		fprintf(stderr,"\033]0;PREDICT: Version %s\007",version); 
}

void ProgramInfo(void)
{
	Banner();
	attrset(COLOR_PAIR(3)|A_BOLD);

	printw("\n\n\n\n\n\t\tPREDICT version : %s\n",version);
	printw("\t\tQTH file loaded : %s\n",qthfile);
	printw("\t\tTLE file loaded : %s\n",tlefile);
	printw("\t\tDatabase file   : ");

	if (database)
		printw("Loaded\n");
	else
		printw("Not loaded\n");

	if (antfd!=-1)
	{
		printw("\t\tAutoTracking    : Sending data to %s",serial_port);

		if (once_per_second)
			printw(" every second");

		printw("\n");
	}

	else
		printw("\t\tAutoTracking    : Not enabled\n");

	printw("\t\tRunning Mode    : ");

	if (socket_flag)
		printw("Network server on port \"%s\"\n",netport);
	else
		printw("Standalone\n");

	printw("\t\tVocalizer       : ");

	if (soundcard)
		printw("Soundcard present");
	else
		printw("No soundcard available");

	refresh();
	attrset(COLOR_PAIR(4)|A_BOLD);
	AnyKey();
}

void NewUser(void)
{

	Banner();
	attrset(COLOR_PAIR(3)|A_BOLD);

	mvprintw(12,2,"WELCOME to PREDICT!  Since you are a new user to the program, default\n");
	printw("  orbital data and ground station location information was copied into\n");
	printw("  your home directory to get you going.  Please select option [G] from\n");
	printw("  PREDICT's main menu to edit your ground station information, and update\n");
	printw("  your orbital database using option [U] or [E].  Enjoy the program!  :-)");
	refresh();

	/* Make "~/.predict" subdirectory */

	sprintf(temp,"%s/.predict",getenv("HOME"));
	mkdir(temp,0777);

	/* Copy default files into ~/.predict directory */

	sprintf(temp,"%sdefault/predict.tle",predictpath);

	CopyFile(temp,tlefile);

	sprintf(temp,"%sdefault/predict.db",predictpath);

	CopyFile(temp,dbfile);

	sprintf(temp,"%sdefault/predict.qth",predictpath);

	CopyFile(temp,qthfile);

	attrset(COLOR_PAIR(4)|A_BOLD);
	AnyKey();
}

void db_edit(void)
{
	clear();
	attrset(COLOR_PAIR(3)|A_BOLD);
	mvprintw(2,15,"* PREDICT Transponder Database Editing Utility *");
	attrset(COLOR_PAIR(2)|A_BOLD);
	mvprintw(13,33,"Coming Soon!");
	attrset(COLOR_PAIR(4)|A_BOLD);
	refresh();
	AnyKey();
}

int QuickFind(const char *string, const char *outputfile)
{
	int x, y, z, step=1;
	long start, now, end, count;
	char satname[50], startstr[20], endstr[20];
	time_t t;
	FILE *fd;

	if (outputfile[0])
		fd=fopen(outputfile,"w");
	else
		fd=stdout;

	startstr[0]=0;
	endstr[0]=0;

	ReadDataFiles();

	for (x=0; x<48 && string[x]!=0 && string[x]!='\n'; x++)
		satname[x]=string[x];

	satname[x]=0;
	x++;

	for (y=0; string[x+y]!=0 && string[x+y]!='\n'; y++)
		startstr[y]=string[x+y];

	startstr[y]=0;
	y++;

	for (z=0; string[x+y+z]!=0 && string[x+y+z]!='\n'; z++)
		endstr[z]=string[x+y+z];

	endstr[z]=0;
 
	/* Do a simple search for the matching satellite name */

	for (z=0; z<24; z++)
	{
		if ((strcmp(sat_table[z].name,satname)==0) || (atol(satname)==sat_table[z].catnum))
		{
			start=atol(startstr);

			if (endstr[strlen(endstr)-1]=='m')
			{
				step=60;
				endstr[strlen(endstr)-1]=0;
			}
			
			if (endstr[0]=='+')
				end=start+((long)step)*atol(endstr);
			else
				end=atol(endstr);

			struct sat_st *sat=&sat_table[z];

			t=time(NULL);
			now=(long)t;

			if (start==0)
				start=now;

			if (startstr[0]=='+')
			{
				start=now;

				if (startstr[strlen(startstr)-1]=='m')
				{
					step=60;
					startstr[strlen(startstr)-1]=0;
				}

				end=start+((long)step)*atol(startstr);

				/* Prevent a list greater than
				   24 hours from being produced */

				if ((end-start)>86400)
				{
					start=now;
					end=now-1;
				}
			}

			if ((start>=now-31557600) && (start<=now+31557600) && end==0)
			{
				/* Start must be one year from now */
				/* Display a single position */
				daynum=((start/86400.0)-3651.0);
				tle_t tle;
				PreCalc(sat, &tle);
				Calc(&tle);

				if (Decayed(sat,daynum)==0)
					fprintf(fd,"%ld %s %4.0f %4.0f %4d %4.0f %4.0f %6.0f %6ld %c\n",
							start,Daynum2String(daynum),
							sat_ele,sat_azi,ma256,
							sat_lat,360.0-sat_lon,
							sat_range,rv,findsun);
				break;
			}

			else
			{
				/* Display a whole list */
				for (count=start; count<=end; count+=step)
				{
					daynum=((count/86400.0)-3651.0);
					tle_t tle;
					PreCalc(sat, &tle);
					Calc(&tle);

					if (Decayed(sat,daynum)==0)
						fprintf(fd,"%ld %s %4.0f %4.0f %4d %4.0f %4.0f %6.0f %6ld %c\n",
								count,Daynum2String(daynum),
								sat_ele,sat_azi,ma256,
								sat_lat,360.0-sat_lon,
								sat_range,rv,findsun);
				}
				break;
			}
		}
	}

	if (outputfile[0])
		fclose(fd);

	return 0;
}

int QuickPredict(const char *string, const char *outputfile)
{
	int x, y, z, lastel=0;
	long start, now;
	double doppler100=0.0;
	char satname[50], startstr[20];
	time_t t;
	FILE *fd;

	if (outputfile[0])
		fd=fopen(outputfile,"w");
	else
		fd=stdout;

	startstr[0]=0;

	ReadDataFiles();

	for (x=0; x<48 && string[x]!=0 && string[x]!='\n'; x++)
		satname[x]=string[x];

	satname[x]=0;
	x++;

	for (y=0; string[x+y]!=0 && string[x+y]!='\n'; y++)
		startstr[y]=string[x+y];

	startstr[y]=0;
	y++;

	/* Do a simple search for the matching satellite name */

	for (z=0; z<24; z++)
	{
		if ((strcmp(sat_table[z].name,satname)==0) || (atol(satname)==sat_table[z].catnum))
		{
			start=atol(startstr);
			struct sat_st *sat=&sat_table[z];

			t=time(NULL);
			now=(long)t;

			if (start==0)
				start=now;

			if ((start>=now-31557600) && (start<=now+31557600))
			{
				/* Start must within one year of now */
				daynum=((start/86400.0)-3651.0);
				tle_t tle;
				PreCalc(sat, &tle);
				Calc(&tle);

				if (AosHappens(sat) && !Geostationary(sat) && !Decayed(sat,daynum))
				{
					/* Make Predictions */
					daynum=FindAOS(sat, &tle);

					/* Display the pass */

					while ( rint(sat_ele)>=0)
					{
						fprintf(fd,"%.0f %s %4.0f %4.0f %4d %4.0f %4.0f %6.0f %6ld %c %f\n",
								floor(86400.0*(3651.0+daynum)),Daynum2String(daynum),
								sat_ele,sat_azi,ma256,
								sat_lat,360.0-sat_lon,
								sat_range,rv,findsun,doppler100);
						lastel=sat_ele;
						daynum+=cos((sat_ele-1.0)*deg2rad)*sqrt(sat_alt)/25000.0;
						Calc(&tle);
					}

					if (lastel!=0)
					{
						daynum=FindLOS(sat, &tle);
						Calc(&tle);
						fprintf(fd,"%.0f %s %4.0f %4.0f %4d %4.0f %4.0f %6.0f %6ld %c %f\n",
								floor(86400.0*(3651.0+daynum)),Daynum2String(daynum),
								sat_ele,sat_azi,ma256,
								sat_lat,360.0-sat_lon,
								sat_range,rv,findsun,doppler100);
					}
				}
				break;
			}
		}
	}

	if (outputfile[0])
		fclose(fd);

	return 0;
}

int QuickDoppler100(const char *string, const char *outputfile)
{

	/* Do a quick predict of the doppler for non-geo sattelites, returns UTC epoch seconds, 
	   UTC time and doppler normalized to 100MHz for every 5 seconds of satellite-pass as a CSV*/

	int x, y, z, lastel=0;
	long start, now;
	double doppler100;
	char satname[50], startstr[20];
	time_t t;
	FILE *fd;

	if (outputfile[0])
		fd=fopen(outputfile,"w");
	else
		fd=stdout;

	startstr[0]=0;

	ReadDataFiles();

	for (x=0; x<48 && string[x]!=0 && string[x]!='\n'; x++)
		satname[x]=string[x];

	satname[x]=0;
	x++;

	for (y=0; string[x+y]!=0 && string[x+y]!='\n'; y++)
		startstr[y]=string[x+y];

	startstr[y]=0;
	y++;

	/* Do a simple search for the matching satellite name */

	for (z=0; z<24; z++)
	{
		if ((strcmp(sat_table[z].name,satname)==0) || (atol(satname)==sat_table[z].catnum))
		{
			start=atol(startstr);
			struct sat_st *sat=&sat_table[z];

			t=time(NULL);
			now=(long)t;

			if (start==0)
				start=now;

			if ((start>=now-31557600) && (start<=now+31557600))
			{
				/* Start must within one year of now */
				daynum=((start/86400.0)-3651.0);
				tle_t tle;
				PreCalc(sat, &tle);
				Calc(&tle);

				if (AosHappens(sat) && !Geostationary(sat) && !Decayed(sat,daynum))
				{
					/* Make Predictions */
					daynum=FindAOS(sat, &tle);

					/* Display the pass */

					while ( rint(sat_ele)>=0)
					{
						doppler100=-100.0e06*((sat_range_rate*1000.0)/299792458.0);
						fprintf(fd,"%.0f,%s,%f\n",floor(86400.0*(3651.0+daynum)),Daynum2String(daynum),doppler100);
						lastel=sat_ele;
						daynum+=cos((sat_ele-1.0)*deg2rad)*sqrt(sat_alt)/500000.0;
						Calc(&tle);
					}

					if (lastel!=0)
					{
						doppler100=-100.0e06*((sat_range_rate*1000.0)/299792458.0);
						daynum=FindLOS(sat, &tle);
						Calc(&tle);
						fprintf(fd,"%.0f,%s,%f\n",floor(86400.0*(3651.0+daynum)),Daynum2String(daynum),doppler100);
					}
				}
				break;
			}
		}
	}

	if (outputfile[0])
		fclose(fd);

	return 0;
}


int main(int argc,char *argv[])
{
	size_t x,z;
	int y, key=0;
	char updatefile[80], quickfind=0, quickpredict=0,
	     quickstring[40], outputfile[42], quickdoppler100=0,
	     tle_cli[50], qth_cli[50], interactive=0;
	struct termios oldtty, newtty;
	pthread_t thread;
	char *env=NULL;
	FILE *db;

	/* Set up translation table for computing TLE checksums */

	for (x=0; x<=255; val[x]=0, x++);
	for (x='0'; x<='9'; val[x]=x-'0', x++);

	val['-']=1;

	updatefile[0]=0;
	outputfile[0]=0;
	temp[0]=0;
	tle_cli[0]=0;
	qth_cli[0]=0;
	dbfile[0]=0;
	netport[0]=0;
	serial_port[0]=0;
	once_per_second=0;
		
	y=argc-1;
	antfd=-1;

	/* Make sure entire "quickstring" array is initialized before use */

	for (x=0; x<40; quickstring[x]=0, x++);

	/* Scan command-line arguments */

	for (x=1; x<=y; x++)
	{
		if (strcmp(argv[x],"-f")==0)
		{
			quickfind=1;
			z=x+1;
			while (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				if ((strlen(quickstring)+strlen(argv[z]))<37)
				{
					strncat(quickstring,argv[z],15);
					strcat(quickstring,"\n");
					z++;
				}
			}
			z--;
		}

		if (strcmp(argv[x],"-p")==0)
		{
			quickpredict=1;
			z=x+1;

			while (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				if ((strlen(quickstring)+strlen(argv[z]))<37)
				{
					strncat(quickstring,argv[z],15);
					strcat(quickstring,"\n");
					z++;
				}
			}
			z--;
		}

		if (strcmp(argv[x],"-dp")==0)
		{
			quickdoppler100=1;
			z=x+1;

			while (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				if ((strlen(quickstring)+strlen(argv[z]))<37)
				{
					strncat(quickstring,argv[z],15);
					strcat(quickstring,"\n");
					z++;
				}
			}
			z--;
		}

		if (strcmp(argv[x],"-u")==0)
		{
			z=x+1;
			while (z<=y && argv[z][0] && argv[z][0]!='-')
			{
				if ((strlen(updatefile)+strlen(argv[z]))<75)
				{
					strncat(updatefile,argv[z],75);
					strcat(updatefile,"\n");
					z++;
				}
			}
			z--;	
		}


		if (strcmp(argv[x],"-t")==0)
		{
			z=x+1;
			if (z<=y && argv[z][0] && argv[z][0]!='-')
				strncpy(tle_cli,argv[z],48);
		}

		if (strcmp(argv[x],"-q")==0)
		{
			z=x+1;
			if (z<=y && argv[z][0] && argv[z][0]!='-')
				strncpy(qth_cli,argv[z],48);
		}

		if (strcmp(argv[x],"-a")==0)
		{
			z=x+1;
			if (z<=y && argv[z][0] && argv[z][0]!='-')
				strncpy(serial_port,argv[z],13);
		}

		if (strcmp(argv[x],"-a1")==0)
		{
			z=x+1;
			if (z<=y && argv[z][0] && argv[z][0]!='-')
				strncpy(serial_port,argv[z],13);
			once_per_second=1;
		}

		if (strcmp(argv[x],"-o")==0)
		{
			z=x+1;
			if (z<=y && argv[z][0] && argv[z][0]!='-')
				strncpy(outputfile,argv[z],40);
		}

		if (strcmp(argv[x],"-n")==0)
		{
			z=x+1;
			if (z<=y && argv[z][0] && argv[z][0]!='-')
				strncpy(netport,argv[z],5);
		}

		if (strcmp(argv[x],"-s")==0)
			socket_flag=1;

		if (strcmp(argv[x],"-north")==0) /* Default */
			io_lat='N';

		if (strcmp(argv[x],"-south")==0)
			io_lat='S';

		if (strcmp(argv[x],"-west")==0)  /* Default */
			io_lon='W';

		if (strcmp(argv[x],"-east")==0)
			io_lon='E';
	}

	/* We're done scanning command-line arguments */

	/* If no command-line (-t or -q) arguments have been passed
	   to PREDICT, create qth and tle filenames based on the
	   default ($HOME) directory. */

	env=getenv("HOME");

	if (qth_cli[0]==0)
		sprintf(qthfile,"%s/.predict/predict.qth",env);
	else
		/* sprintf(qthfile,"%s%c",qth_cli,0); */
		sprintf(qthfile,"%s",qth_cli);

	if (tle_cli[0]==0)
		sprintf(tlefile,"%s/.predict/predict.tle",env);
	else
		/* sprintf(tlefile,"%s%c",tle_cli,0); */
		sprintf(tlefile,"%s",tle_cli);

	/* Test for interactive/non-interactive mode of operation
	   based on command-line arguments given to PREDICT. */

	if (updatefile[0] || quickfind || quickpredict || quickdoppler100)
		interactive=0;
	else
		interactive=1;

	if (interactive)
	{
		sprintf(dbfile,"%s/.predict/predict.db",env);

		/* If the transponder database file doesn't already
		   exist under $HOME/.predict, and a working environment
		   is available, place a default copy from the PREDICT
		   distribution under $HOME/.predict. */

		db=fopen(dbfile,"r");
		if (!db)
		{
			sprintf(temp,"%sdefault/predict.db",predictpath);
			CopyFile(temp,dbfile);
		}

		else
			fclose(db);
	}

	x=ReadDataFiles();

	if (x>1)  /* TLE file was loaded successfully */
	{
		if (updatefile[0])  /* -u was passed to PREDICT */
		{
			y=0;
			z=0;
			temp[0]=0;

			while (updatefile[y]!=0)
			{
				while (updatefile[y]!='\n' && updatefile[y]!=0 && y<79)
				{
					temp[z]=updatefile[y];
					z++;
					y++;
				}

				temp[z]=0;

				if (temp[0])
				{
					AutoUpdate(temp);
					temp[0]=0;
					z=0;
					y++;
				}
			}

			exit(0);
		}
	}

	if (x==3)  /* Both TLE and QTH files were loaded successfully */
	{
		if (quickfind)  /* -f was passed to PREDICT */
			exit(QuickFind(quickstring,outputfile));

		if (quickpredict)  /* -p was passed to PREDICT */
			exit(QuickPredict(quickstring,outputfile));

		if (quickdoppler100)  /* -dp was passed to PREDICT */
			exit(QuickDoppler100(quickstring,outputfile));
	}

	else
	{
		if (tle_cli[0] || qth_cli[0])
		{
			/* "Houston, we have a problem..." */

			printf("\n%c",7);

			if (x^1)
				printf("*** ERROR!  Your QTH file \"%s\" could not be loaded!\n",qthfile);

			if (x^2)
				printf("*** ERROR!  Your TLE file \"%s\" could not be loaded!\n",tlefile);

			printf("\n");

			exit(-1);
		}
	}

	if (interactive)
	{
		/* We're in interactive mode.  Prepare the screen */

		/* Are we running under an xterm or equivalent? */

		env=getenv("TERM");

		if (env!=NULL && strncmp(env,"xterm",5)==0)
			xterm=1;
		else
			xterm=0; 

		/* Start ncurses */

		initscr();
		start_color();
		cbreak();
		noecho();
		scrollok(stdscr,TRUE);
		curs_set(0);

		init_pair(1,COLOR_WHITE,COLOR_BLACK);
		init_pair(2,COLOR_WHITE,COLOR_BLUE);
		init_pair(3,COLOR_YELLOW,COLOR_BLUE);
		init_pair(4,COLOR_CYAN,COLOR_BLUE);
		init_pair(5,COLOR_WHITE,COLOR_RED);
		init_pair(6,COLOR_RED,COLOR_WHITE);
		init_pair(7,COLOR_CYAN,COLOR_RED);

		if (x<3)
		{
			/* A problem occurred reading the
			   default QTH and TLE files, and
			   no -t or -q options were
			   provided on the command-line. */

			NewUser();
			x=ReadDataFiles();
			QthEdit();
		}
	}

	if (x==3)
	{
		/* Open serial port to send data to
		   the antenna tracker if present. */

		if (serial_port[0])
		{
			/* Make sure there's no trailing '/' */

			x=strlen(serial_port);

			if (serial_port[x-1]=='/')
				serial_port[x-1]=0;

			antfd=open(serial_port, O_WRONLY|O_NOCTTY);

			if (antfd!=-1)
			{
				/* Set up serial port */

				tcgetattr(antfd, &oldtty);
				memset(&newtty, 0, sizeof(newtty));

				/* 9600 baud, 8-bits, no parity,
				   1-stop bit, no handshaking */

				newtty.c_cflag=B9600|CS8|CLOCAL;
				newtty.c_iflag=IGNPAR;
				newtty.c_oflag=0;
				newtty.c_lflag=0;

				tcflush(antfd, TCIFLUSH);
				tcsetattr(antfd, TCSANOW, &newtty);
			}

			else
			{
				bailout("Unable To Open Antenna Port");
				exit(-1);
			}
		}
	
		/* Socket activated here.  Remember that
		   the socket data is updated only when
		   running in the real-time tracking modes. */

		if (socket_flag)
		{
			pthread_create(&thread,NULL,(void *)socket_server,(void *)argv[0]);
			bkgdset(COLOR_PAIR(3));
			MultiTrack();
		}

		MainMenu();

		do
		{	
			key=getch();

			if (key!='T')
				key=tolower(key);

			struct sat_st *sat;
			switch (key)
			{
				case 'p':
				case 'v':
					Print(NULL, "",0);
					PrintVisible(NULL, "");
					sat=Select();

					if (sat && sat->meanmo!=0.0 && !Decayed(sat,0.0))
						Predict(sat, key);

					MainMenu();
					break;

				case 'l':
					Print(NULL, "",0);
					PredictMoon();
					MainMenu();
					break;

				case 'o':
					Print(NULL, "",0);
					PredictSun();
					MainMenu();
					break;

				case 'u':
					AutoUpdate("");
					MainMenu();
					break;

				case 'e':
					KepEdit();
					MainMenu();
					break;

				case 'd':
					ShowOrbitData();
					MainMenu();
					break;

				case 'g':
					QthEdit();
					MainMenu();
					break;

				case 't':
				case 'T':
					sat=Select();

					if (sat && sat->meanmo!=0.0 && !Decayed(sat,0.0))
						SingleTrack(sat, key);

					MainMenu();
					break;

				case 'm':
					MultiTrack();
					MainMenu();
					break;

				case 'i':
					ProgramInfo();
					MainMenu();
					break;

				case 'b':
					db_edit();
					MainMenu();
					break;

				case 's':
					sat=Select();
					if (sat && sat->meanmo!=0.0 && !Decayed(sat,0.0))
					{
						Print(sat, "",0);
						Illumination(sat);
					}
					MainMenu();
					break;
			}

		} while (key!='q' && key!=27);

		if (antfd!=-1)
		{
			tcsetattr(antfd,TCSANOW,&oldtty);
			close(antfd);
		}

		curs_set(1);	
		bkgdset(COLOR_PAIR(1));
		clear();
		refresh();
		endwin();
	}

	exit(0);
}

