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

#include "more_math.h"
#include "constants.h"

#include <math.h>
#include <time.h>


int Sign(double arg)
{
	/* Returns sign of a double */

	if (arg>0)
		return 1;
	else if (arg<0)
		return -1;
	else
		return 0;
}

double Sqr(double arg)
{
	/* Returns square of a double */
	return (arg*arg);
}

double Cube(double arg)
{
	/* Returns cube of a double */
	return (arg*arg*arg);
}

double Radians(double arg)
{
	/* Returns angle in radians from argument in degrees */
	return arg*deg2rad;
}

double Degrees(double arg)
{
	/* Returns angle in degrees from argument in radians */
	return arg/deg2rad;
}

double ArcSin(double arg)
{
	/* Returns the arcsine of the argument */

	if (fabs(arg)>=1.0)
		return(Sign(arg)*pio2);

	return(atan(arg/sqrt(1.0-arg*arg)));
}

double ArcCos(double arg)
{
	/* Returns arccosine of argument */
	return pio2-ArcSin(arg);
}

void Magnitude(vector_t *v)
{
	/* Calculates scalar magnitude of a vector_t argument */
	v->w=sqrt(Sqr(v->x)+Sqr(v->y)+Sqr(v->z));
}

void Vec_Add(const vector_t *v1, const vector_t *v2, vector_t *v3)
{
	/* Adds vectors v1 and v2 together to produce v3 */
	v3->x=v1->x+v2->x;
	v3->y=v1->y+v2->y;
	v3->z=v1->z+v2->z;
	Magnitude(v3);
}

void Vec_Sub(const vector_t *v1, const vector_t *v2, vector_t *v3)
{
	/* Subtracts vector v2 from v1 to produce v3 */
	v3->x=v1->x-v2->x;
	v3->y=v1->y-v2->y;
	v3->z=v1->z-v2->z;
	Magnitude(v3);
}

void Scalar_Multiply(double k, const vector_t *v1, vector_t *v2)
{
	/* Multiplies the vector v1 by the scalar k to produce the vector v2 */
	v2->x=k*v1->x;
	v2->y=k*v1->y;
	v2->z=k*v1->z;
	v2->w=fabs(k)*v1->w;
}

void Scale_Vector(double k, vector_t *v)
{ 
	/* Multiplies the vector v1 by the scalar k */
	v->x*=k;
	v->y*=k;
	v->z*=k;
	Magnitude(v);
}

double Dot(const vector_t *v1, const vector_t *v2)
{
	/* Returns the dot product of two vectors */
	return (v1->x*v2->x+v1->y*v2->y+v1->z*v2->z);
}

double Angle(vector_t *v1, vector_t *v2)
{
	/* Calculates the angle between vectors v1 and v2 */
	Magnitude(v1);
	Magnitude(v2);
	return(ArcCos(Dot(v1,v2)/(v1->w*v2->w)));
}

void Cross(const vector_t *v1, const vector_t *v2 ,vector_t *v3)
{
	/* Produces cross product of v1 and v2, and returns in v3 */
	v3->x=v1->y*v2->z-v1->z*v2->y;
	v3->y=v1->z*v2->x-v1->x*v2->z;
	v3->z=v1->x*v2->y-v1->y*v2->x;
	Magnitude(v3);
}

void Normalize(vector_t *v)
{
	/* Normalizes a vector */
	v->x/=v->w;
	v->y/=v->w;
	v->z/=v->w;
	v->w=1.0;
}

double AcTan(double sinx, double cosx)
{
	/* Four-quadrant arctan function */

	if (cosx==0.0)
	{
		if (sinx>0.0)
			return (pio2);
		else
			return (x3pio2);
	}

	else
	{
		if (cosx>0.0)
		{
			if (sinx>0.0)
				return (atan(sinx/cosx));
			else
				return (twopi+atan(sinx/cosx));
		}

		else
			return (pi+atan(sinx/cosx));
	}
}

double FMod2p(double x)
{
	/* Returns mod 2PI of argument */

	double ret_val=x;
	const int i=ret_val/twopi;
	ret_val-=i*twopi;

	if (ret_val<0.0)
		ret_val+=twopi;

	return ret_val;
}

double Modulus(double arg1, double arg2)
{
	/* Returns arg1 mod arg2 */

	double ret_val=arg1;
	const int i=ret_val/arg2;
	ret_val-=i*arg2;

	if (ret_val<0.0)
		ret_val+=arg2;

	return ret_val;
}

double Frac(double arg)
{
	/* Returns fractional part of double argument */
	return(arg-floor(arg));
}

int Round(double arg)
{
	/* Returns argument rounded up to nearest integer */
	return((int)floor(arg+0.5));
}

double Int(double arg)
{
	/* Returns the floor integer of a double arguement, as double */
	return(floor(arg));
}

double Julian_Date_of_Year(double year)
{
	/* The function Julian_Date_of_Year calculates the Julian Date  */
	/* of Day 0.0 of {year}. This function is used to calculate the */
	/* Julian Date of any date by using Julian_Date_of_Year, DOY,   */
	/* and Fraction_of_Day. */

	/* Astronomical Formulae for Calculators, Jean Meeus, */
	/* pages 23-25. Calculate Julian Date of 0.0 Jan year */

	year=year-1;
	const long A=year/100;
	long i=A/4;
	const long B=2-A+i;
	i=365.25*year;
	i+=30.6001*14;
	return i+1720994.5+B;
}

double Julian_Date_of_Epoch(double epoch)
{ 
	/* The function Julian_Date_of_Epoch returns the Julian Date of     */
	/* an epoch specified in the format used in the NORAD two-line      */
	/* element sets. It has been modified to support dates beyond       */
	/* the year 1999 assuming that two-digit years in the range 00-56   */
	/* correspond to 2000-2056. Until the two-line element set format   */
	/* is changed, it is only valid for dates through 2056 December 31. */


	/* Modification to support Y2K */
	/* Valid 1957 through 2056     */

	double year;
	const double day=modf(epoch*1E-3, &year)*1E3;

	if (year<57)
		year=year+2000;
	else
		year=year+1900;

	return (Julian_Date_of_Year(year)+day);
}

static int DayOfYear (int yr, int mo, int dy)
{
	/* The function DOY calculates the day of the year for the specified */
	/* date. The calculation uses the rules for the Gregorian calendar   */
	/* and is valid from the inception of that calendar system.          */

	static const int days[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

	int day=0;
	for (int i=0; i<mo-1; i++ )
		day+=days[i];

	day=day+dy;

	/* Leap year correction */

	if ((yr%4==0) && ((yr%100!=0) || (yr%400==0)) && (mo>2))
		day++;

	return day;
}

static double Fraction_of_Day(int hr, int mi, double se)
{
	/* Fraction_of_Day calculates the fraction of */
	/* a day passed at the specified input time.  */

	const double dhr=(double)hr;
	const double dmi=(double)mi;

	return ((dhr+(dmi+se/60.0)/60.0)/24.0);
}

double Julian_Date(struct tm *cdate)
{
	/* The function Julian_Date converts a standard calendar   */
	/* date and time to a Julian Date. The procedure Date_Time */
	/* performs the inverse of this function. */

	return  Julian_Date_of_Year(cdate->tm_year)+
	        DayOfYear(cdate->tm_year,cdate->tm_mon,cdate->tm_mday)+
	        Fraction_of_Day(cdate->tm_hour,cdate->tm_min,cdate->tm_sec)+5.787037e-06; /* Round up to nearest 1 sec */
}

void Date_Time(double julian_date, struct tm *cdate)
{
	/* The function Date_Time() converts a Julian Date to
	standard calendar date and time. The function
	Julian_Date() performs the inverse of this function. */

	time_t jtime=(julian_date-2440587.5)*86400.0;
	*cdate=*gmtime(&jtime);
}

double Delta_ET(double year)
{
	/* The function Delta_ET has been added to allow calculations on   */
	/* the position of the sun.  It provides the difference between UT */
	/* (approximately the same as UTC) and ET (now referred to as TDT).*/
	/* This function is based on a least squares fit of data from 1950 */
	/* to 1991 and will need to be updated periodically. */

	/* Values determined using data from 1950-1991 in the 1990 
	Astronomical Almanac.  See DELTA_ET.WQ1 for details. */

	return 26.465+0.747622*(year-1950)+1.886913*sin(twopi*(year-1975)/33);
}

double ThetaG_JD(double jd)
{
	/* Reference:  The 1992 Astronomical Almanac, page B6. */

	const double UT=Frac(jd+0.5);
	jd=jd-UT;
	const double TU=(jd-2451545.0)/36525;
	double GMST=24110.54841+TU*(8640184.812866+TU*(0.093104-TU*6.2E-6));
	GMST=Modulus(GMST+secday*omega_E*UT,secday);

	return (twopi*GMST/secday);
}

