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
/* General three-dimensional vector structure used by SGP4/SDP4 code. */

#ifndef MORE_MATH_H_
#define MORE_MATH_H_


int Sign(double arg);	/* Returns sign of a double */
double Sqr(double arg);	/* Returns square of a double */
double Cube(double arg);/* Returns cube of a double */
double Radians(double arg);/* Returns angle in radians from argument in degrees */
double Degrees(double arg); /* Returns angle in degrees from argument in radians */
double ArcSin(double arg);/* Returns the arcsine of the argument */
double ArcCos(double arg); /* Returns arccosine of argument */

typedef struct vector_t {
	double x, y, z, w;
}  vector_t;


void Magnitude(vector_t *v);/* Calculates scalar magnitude of a vector_t argument */
void Vec_Add(const vector_t *v1, const vector_t *v2, vector_t *v3); /* Adds vectors v1 and v2 together to produce v3 */
void Vec_Sub(const vector_t *v1, const vector_t *v2, vector_t *v3); /* Subtracts vector v2 from v1 to produce v3 */
void Scalar_Multiply(double k, const vector_t *v1, vector_t *v2); /* Multiplies the vector v1 by the scalar k to produce the vector v2 */
void Scale_Vector(double k, vector_t *v);	/* Multiplies the vector v1 by the scalar k */
double Dot(const vector_t *v1, const vector_t *v2);	/* Returns the dot product of two vectors */
double Angle(vector_t *v1, vector_t *v2);	/* Calculates the angle between vectors v1 and v2 */
void Cross(const vector_t *v1, const vector_t *v2 ,vector_t *v3);	/* Produces cross product of v1 and v2, and returns in v3 */
void Normalize(vector_t *v);	/* Normalizes a vector */
double AcTan(double sinx, double cosx);	/* Four-quadrant arctan function */
double FMod2p(double x);	/* Returns mod 2PI of argument */
double Modulus(double arg1, double arg2);	/* Returns arg1 mod arg2 */
double Frac(double arg);	/* Returns fractional part of double argument */
int Round(double arg);	/* Returns argument rounded up to nearest integer */
double Int(double arg);	/* Returns the floor integer of a double arguement, as double */

double Julian_Date_of_Year(double year);	/* The function Julian_Date_of_Year calculates the Julian Date  */
double Julian_Date_of_Epoch(double epoch);	/* The function Julian_Date_of_Epoch returns the Julian Date of an epoch */

/* The function Julian_Date converts a standard calendar   */
/* date and time to a Julian Date. The procedure Date_Time */
/* performs the inverse of this function. */
struct tm;
double Julian_Date(struct tm *cdate);

/* The function Date_Time() converts a Julian Date to
standard calendar date and time. The function
Julian_Date() performs the inverse of this function. */
void Date_Time(double julian_date, struct tm *cdate);

/* The function Delta_ET has been added to allow calculations on   */
/* the position of the sun.  It provides the difference between UT */
/* (approximately the same as UTC) and ET (now referred to as TDT).*/
/* This function is based on a least squares fit of data from 1950 */
/* to 1991 and will need to be updated periodically. */
double Delta_ET(double year);
double ThetaG_JD(double jd);

#endif /* MORE_MATH_H_ */
