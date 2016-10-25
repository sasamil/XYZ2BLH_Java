/***************************************************************************
 *   Copyright (C) 2016 by Саша Миленковић                                 *
 *   sasa.milenkovic.xyz@gmail.com                                         *
 *   									   *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *   ( http://www.gnu.org/licenses/gpl-3.0.en.html )                       *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
 
package xyz2blh;

import java.lang.*;

/**
 *
 * @author Саша Миленковић
 */
public class Main {
  static public double tolerance  =  1.e-12;
  static public double mtolerance = -tolerance;
  static public double dtolerance =  tolerance * tolerance;

  static public double PI = 3.141592653589793238463;
  static public double TWOPI = PI * 2.0;
  static public double HALFPI = PI * 0.5;
  static public double ro = 3600.0 * 180.0 / PI; // rad2sec
  static public double ONETHIRD = 1./3.;

  //proj constants
  static public double  a_bessel	 = 6377397.155;
  static public double  a_bessel_2	 = a_bessel*a_bessel;
  static public double  b_bessel	 = 6356078.962818189;
  static public double  b_bessel_2	 = b_bessel*b_bessel;
  static public double  f_bessel	 = 1.0 - b_bessel/a_bessel; //1.0/299.1528128L;
  static public double  f_bessel2	 = a_bessel/b_bessel - 1.0;
  static public double  j_e2_bessel      = b_bessel_2/a_bessel_2;
  static public double  e2_bessel	 = 1.0 - j_e2_bessel; //0.0816968312225269L;
  static public double  e2_bessel2	 = (a_bessel_2 - b_bessel_2) / b_bessel;

  static public double  a_wgs	  = 6378137.0;
  static public double  a_wgs_2  = a_wgs*a_wgs;
  static public double  b_wgs	  = 6356752.314245179;
  static public double  b_wgs_2  = b_wgs*b_wgs;
  static public double  f_wgs	  = (a_wgs - b_wgs) / a_wgs; //1.0L/298.257223563L;
  static public double  f_wgs2	  = (a_wgs - b_wgs) / b_wgs;
  static public double  j_e2_wgs = b_wgs_2/a_wgs_2;
  static public double  e2_wgs	  = 1.0 - j_e2_wgs; //0.08181919084262149L;
  static public double  e2_wgs2 = (a_wgs_2 - b_wgs_2) / b_wgs;

//---------------------------------------------------------------------------
// equation to solve:  x^3 + a*x^2 + b*x + c
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] ± i*x[2], return 1
public static int solveP3(double[] x, double a, double b, double c)
{
    double a2 = a*a;
    double q  = (a2 - 3.*b)/9.;
    double r  = (a*(2.*a2-9.*b) + 27.*c)/54.;
    double r2 = r*r;
    double q3 = q*q*q;
    double A,B;
    if(r2<q3) {
        double t=r/Math.sqrt(q3);
        if( t<-1.) t=-1.;
        if( t> 1.) t= 1.;
        t=Math.acos(t);
        a/=3.; q=-2.*Math.sqrt(q);
        x[0]=q*Math.cos(t/3.)-a;
        x[1]=q*Math.cos((t+TWOPI)/3.)-a;
        x[2]=q*Math.cos((t-TWOPI)/3.)-a;
        return(3);
    } 
    else {
        //A =-pow(Math.abs(r)+sqrt(r2-q3),1./3.);
        A = -Math.exp(Math.log(Math.abs(r)+Math.sqrt(r2-q3)) * ONETHIRD);
        if( r<0. ) A=-A;
        //B = A==0? 0 : B=q/A;
        B = (0.==A ? 0. : q/A);

        a/=3.;
        x[0] =(A+B)-a;
        x[1] =-0.5*(A+B)-a;
        x[2] = 0.5*Math.sqrt(3.)*(A-B);
	if(Math.abs(x[2])<tolerance) { x[2]=x[1]; return(2); }
        return(1);
    }
}// solveP3

//===========================================================================
// equation to solve:  x^4 + a*x^3 + b*x62 + c*x + d
// re - array of size 4, im - array of size 2
// returns number of
// In case 4 real roots:  re[0], re[1], re[2], re[3]                                              return 4
//         2 real roots:  re[0], re[1], re[2] = re[3]  or   x[0] ± i*x[1], x[2], x[3]             return 2
//         0 real roots : x[0] ± i*x[1], x[2] ± i*x[3]                                            return 0
// In case x[1] or x[3] are imaginnary -
public static int solveP4(double[] re, double[] im, double a, double b, double c, double d)
{
    double q1, q2, p1, p2, D, sqd, y;
    double a3 = -b;
    double b3 =  a*c -4.*d;
    double c3 = -a*a*d - c*c + 4.*b*d;
    double x3[] = new double[3];

    int iRetval=4, iZeroes = solveP3(x3, a3, b3, c3);

    y = x3[0];
    if(iZeroes != 1.)
    {
        if(Math.abs(x3[1]) > Math.abs(y)) y = x3[1];
        if(Math.abs(x3[2]) > Math.abs(y)) y = x3[2];
    }

    D = y*y - 4.*d;
    if(D<tolerance /*Math.abs(D) < eps*/)
    {
        q1 = q2 = y * 0.5;

        D = a*a - 4.*(b-y);
        if(D<tolerance /*Math.abs(D) < eps*/)
            p1 = p2 = a * 0.5;
        else
        {
            sqd = Math.sqrt(D);
            p1 = (a + sqd/*sqrt(D)*/) * 0.5;
            p2 = (a - sqd/*sqrt(D)*/) * 0.5;
        }
    }
    else
    {
        sqd = Math.sqrt(D);
        q1 = (y + sqd/*sqrt(D)*/) * 0.5;
        q2 = (y - sqd/*sqrt(D)*/) * 0.5;

        p1 = (a*q1-c)/(q1-q2);
        p2 = (c-a*q2)/(q1-q2);
    }

    im[0] = im[1] = 0.0;

    // x^2 + p1*x + q1 = 0
    D = p1*p1 - 4.*q1;
    if(D < 0.0)
    {
        iRetval -= 2.;
        re[0] = re[1] = -p1 * 0.5;
        im[0] = Math.sqrt(-D) * 0.5;
    }
    else
    {
        sqd = Math.sqrt(D);
        re[0] = (-p1 + sqd/*sqrt(D)*/) * 0.5;
        re[1] = (-p1 - sqd/*sqrt(D)*/) * 0.5;
    }

    // x^2 + p2*x + q2 = 0
    D = p2*p2 - 4.*q2;
    if(D < 0.0)
    {
        iRetval -= 2.;
        re[2] = re[3] = -p2 * 0.5;
        im[1] = Math.sqrt(-D) * 0.5;
    }
    else
    {
        sqd = Math.sqrt(D);
        re[2] = (-p2 + sqd/*sqrt(D)*/) * 0.5;
        re[3] = (-p2 - sqd/*sqrt(D)*/) * 0.5;
    }

    return iRetval;
}// solveP4

//===========================================================================
//input: fi, lambda (in radians) and h (m)
//output: geocentric coordinates (m)
public static XYZtriplet BLh2XYZ(double a, double b, XYZtriplet blh)
{
	double fi = blh.x;
	double lambda = blh.y;
	double h = blh.z;

	double e_2 = 1.0 - b*b/a/a;
	double sn = Math.sin(fi);
	double cs = Math.cos(fi);
	double N = a/Math.sqrt(1.0 - e_2*sn*sn);

	double X = (N + h) * cs*Math.cos(lambda);
	double Y = (N + h) * cs*Math.sin(lambda);
	double Z = (N*(1-e_2) + h) * sn;

	XYZtriplet retval = new XYZtriplet();
	retval.x = X;
	retval.y = Y;
	retval.z = Z;

	return retval;
}

//===========================================================================
//input: geocentric coordinates (m)
//output: fi, lambda (in radians) and h (m)
public static XYZtriplet geocentric_to_geodetic(double a, double b, XYZtriplet p)
{
/* local defintions and variables */
/* end-criterium of loop, accuracy of sin(Latitude) */
    int maxiter = 30;
    double genau = tolerance; //1.E-12;
    double genau2 = (genau*genau);
    double e_2 = 1.0 - b*b/a/a;

    boolean At_Pole;    /* indicates location is in polar region */
    int iter;        /* # of continous iteration, max. 30 is always enough (s.a.) */
    double P;        /* distance between semi-minor axis and location */
    double RR;       /* distance between center and location */
    double CT;       /* sin of geocentric latitude */
    double ST;       /* cos of geocentric latitude */
    double RX;
    double RK;
    double RN;       /* Earth radius at location */
    double CPHI0;    /* cos of start or old geodetic latitude in iterations */
    double SPHI0;    /* sin of start or old geodetic latitude in iterations */
    double CPHI;     /* cos of searched geodetic latitude */
    double SPHI;     /* sin of searched geodetic latitude */
    double SDPHI;    /* end-criterium: addition-theorem of sin(Latitude(iter)-Latitude(iter-1)) */

    double X = p.x;
    double Y = p.y;
    double Z = (0.0!=p.z ? p.z : 0.0);   //Z value not always supplied
    double Longitude;
    double Latitude;
    double Height;
    XYZtriplet retval = new XYZtriplet();
    At_Pole = false;
    P = Math.sqrt(X*X+Y*Y);
    RR = Math.sqrt(X*X+Y*Y+Z*Z);

/*  special cases for latitude and longitude */
    if (P/a < genau) {

/*  special case, if P=0. (X=0., Y=0.) */
        At_Pole = true;
        Longitude = 0.0;

/*  if (X,Y,Z)=(0.,0.,0.) then Height becomes semi-minor axis
 *  of ellipsoid (=center of mass), Latitude becomes PI/2 */
        if (RR/a < genau) {
            Latitude = HALFPI;
            Height   = -b;

            retval.x = Latitude;
            retval.y = Longitude;
            retval.z = Height;
            return retval;
        }
    }
    else {
/*  ellipsoidal (geodetic) longitude
 *  interval: -PI < Longitude <= +PI */
        Longitude=Math.atan2(Y,X);
    }
/* --------------------------------------------------------------
 * Following iterative algorithm was developped by
 * "Institut für Erdmessung", University of Hannover, July 1988.
 * Internet: www.ife.uni-hannover.de
 * Iterative computation of CPHI,SPHI and Height.
 * Iteration of CPHI and SPHI to 10**-12 radian resp.
 * 2*10**-7 arcsec.
 * --------------------------------------------------------------
 */
    CT = Z/RR;
    ST = P/RR;
    RX = 1.0/Math.sqrt(1.0-e_2*(2.0-e_2)*ST*ST);
    CPHI0 = ST*(1.0-e_2)*RX;
    SPHI0 = CT*RX;
    iter = 0;

/* loop to find sin(Latitude) resp. Latitude
 * until |sin(Latitude(iter)-Latitude(iter-1))| < genau */
    do
    {
        iter++;
        RN = a/Math.sqrt(1.0-e_2*SPHI0*SPHI0);

/*  ellipsoidal (geodetic) height */
        Height = P*CPHI0+Z*SPHI0-RN*(1.0-e_2*SPHI0*SPHI0);

        RK = e_2*RN/(RN+Height);
        RX = 1.0/Math.sqrt(1.0-RK*(2.0-RK)*ST*ST);
        CPHI = ST*(1.0-RK)*RX;
        SPHI = CT*RX;
        SDPHI = SPHI*CPHI0-CPHI*SPHI0;
        CPHI0 = CPHI;
        SPHI0 = SPHI;
    }
    while (SDPHI*SDPHI > genau2 && iter < maxiter);

/*      ellipsoidal (geodetic) latitude */
    Latitude=Math.atan(SPHI/Math.abs(CPHI));

    retval.x = Latitude;
    retval.y = Longitude;
    retval.z = Height;
    
    return retval;
  } // geocentric_to_geodetic()

//---------------------------------------------------------------------------
// The method used here is derived from 'An Improved Algorithm for
// Geocentric to Geodetic Coordinate Conversion, by Ralph Toms, Feb 1996
// Note: Variable names follow the notation used in Toms, Feb 1996
public static XYZtriplet geocentric_to_geodetic_bowring (double Geocent_a, double Geocent_b, XYZtriplet p)
{ 
    double PI_OVER_2 = 3.14159265358979323e0 / 2.0;
    int FALSE = 0;
    int TRUE = 1;
    double COS_67P5 = 0.38268343236508977;  /* cosine of 67.5 degrees */
    double AD_C = 1.0026000;            /* Toms region 1 constant */

    double W;        /* distance from Z axis */
    double W2;       /* square of distance from Z axis */
    double T0;       /* initial estimate of vertical component */
    double T1;       /* corrected estimate of vertical component */
    double S0;       /* initial estimate of horizontal component */
    double S1;       /* corrected estimate of horizontal component */
    double Sin_B0;   /* sin(B0), B0 is estimate of Bowring aux variable */
    double Sin3_B0;  /* cube of sin(B0) */
    double Cos_B0;   /* cos(B0) */
    double Sin_p1;   /* sin(phi1), phi1 is estimated latitude */
    double Cos_p1;   /* cos(phi1) */
    double Rn;       /* Earth radius at location */
    double Sum;      /* numerator of cos(phi1) */
    int At_Pole;     /* indicates location is in polar region */

    /* sm */
    double X = p.x;
    double Y = p.y;
    double Z = p.z;
    double Geocent_a2 = Geocent_a * Geocent_a;
    double Geocent_b2 = Geocent_b * Geocent_b;
    double Geocent_e2 = (Geocent_a2 - Geocent_b2) / Geocent_a2;
    double Geocent_ep2 = (Geocent_a2 - Geocent_b2) / Geocent_b2;

    XYZtriplet blh = new XYZtriplet();

    At_Pole = FALSE;
    if (X != 0.0)
    {
        //*Longitude = atan2(Y,X);
        blh.y = Math.atan2(Y,X);
    }
    else
    {
        if (Y > 0)
        {
            //*Longitude = PI_OVER_2;
            blh.y = PI_OVER_2;
        }
        else if (Y < 0)
        {
            //*Longitude = -PI_OVER_2;
            blh.y = -PI_OVER_2;
        }
        else
        {
            At_Pole = TRUE;
            //*Longitude = 0.0;
            blh.y = 0.0;

            if (Z > 0.0)
            {  /* north pole */
                //*Latitude = PI_OVER_2;
                blh.x = PI_OVER_2;
            }
            else if (Z < 0.0)
            {  /* south pole */
                //*Latitude = -PI_OVER_2;
                blh.x = -PI_OVER_2;
            }
            else
            {  /* center of earth */
                //*Latitude = PI_OVER_2;
                //*Height = -Geocent_b;
                blh.x = PI_OVER_2;
                blh.z = -Geocent_b;
                return blh;
            }
        }
    }
    W2 = X*X + Y*Y;
    W = Math.sqrt(W2);
    T0 = Z * AD_C;
    S0 = Math.sqrt(T0 * T0 + W2);
    Sin_B0 = T0 / S0;
    Cos_B0 = W / S0;
    Sin3_B0 = Sin_B0 * Sin_B0 * Sin_B0;
    T1 = Z + Geocent_b * Geocent_ep2 * Sin3_B0;
    Sum = W - Geocent_a * Geocent_e2 * Cos_B0 * Cos_B0 * Cos_B0;
    S1 = Math.sqrt(T1*T1 + Sum * Sum);
    Sin_p1 = T1 / S1;
    Cos_p1 = Sum / S1;
    Rn = Geocent_a / Math.sqrt(1.0 - Geocent_e2 * Sin_p1 * Sin_p1);

    blh.x = Math.atan(Sin_p1 / Cos_p1);

    if (Cos_p1 >= COS_67P5)
    {
        //*Height = W / Cos_p1 - Rn;
        blh.z = W / Cos_p1 - Rn;
    }
    else if (Cos_p1 <= -COS_67P5)
    {
        //*Height = W / -Cos_p1 - Rn;
        blh.z = W / -Cos_p1 - Rn;
    }
    else
    {
        //*Height = Z / Sin_p1 + Rn * (Geocent_e2 - 1.0);
        blh.z = Z / Sin_p1 + Rn * (Geocent_e2 - 1.0);
    }
    if (At_Pole != FALSE)
    {
        //*Latitude = atan(Sin_p1 / Cos_p1);
        blh.z = Math.atan(Sin_p1 / Cos_p1);
    }

    return blh;
} // geocentric_to_geodetic_bowring

//---------------------------------------------------------------------------
// Bowring: /home/sasamil/Geodesy/BLH2XYZ/COMPARISON OF DIFFERENT ALGORITHMS.pdf  or
// /home/sasamil/Geodesy/Books/GPS_Hofmann-Wellenhof (p.232)   or
// /home/sasamil/Geodesy/BLH2XYZ/Comparision2!! (p.8)
public static XYZtriplet bowring2(double a, double b, XYZtriplet p)
{
	double f2  = a/b;
	double epr2 = f2*f2 - 1.0;
	double e2 = (a*a - b*b)/(a*a);
	double X = p.x;
	double Y = p.y;
	double Z = p.z;
	double P = Math.sqrt(X*X + Y*Y);

	double tgtheta, cstheta, sntheta, snfi, N, tgfi, tgfiold, delta;
	XYZtriplet blh = new XYZtriplet();

    if(0.0==X && 0.0==Y)
    {
        blh.x = HALFPI;
        blh.y = 0.0;
        blh.z = Z - b;
        return blh;
    }
    else if(0.0==Z)
    {
        blh.x = 0.0;
        blh.y = Math.atan(Y/X);
        blh.z = P - a;
        return blh;
    }

    else
    {
        tgfi = f2*f2*Z/P;
        do
        {
            tgfiold = tgfi;
            tgtheta = tgfi/f2;
            //double theta = atan(f2*Z/P);
            cstheta = 1./Math.sqrt(1. + tgtheta*tgtheta); //cos(theta);
            sntheta = tgtheta*cstheta; //tg/sqrt(1. + tg*tg); //sin(theta);

            tgfi = (Z + epr2*b*sntheta*sntheta*sntheta)/(P - e2*a*cstheta*cstheta*cstheta);
            //fi = atan(tg2);
            //snfi = sin(fi);
            delta = tgfi-tgfiold;
        }
        while(delta*delta > dtolerance);

        snfi = tgfi / Math.sqrt(1. + tgfi*tgfi);
        N = a / Math.sqrt(1. - e2*snfi*snfi);

        blh.x = Math.atan(tgfi);  //E = atan( (1-f_bessel)*tan(B) );
        blh.y = Math.atan(Y/X);
        blh.z = P*tgfi/snfi - N;

        return blh;
    }
} // bowring

//---------------------------------------------------------------------------
// Bowring: http://gis.stackexchange.com/questions/28446/..
// ..computational-most-efficient-way-to-convert-cartesian-to-geodetic-coordinates
public static XYZtriplet bowring3(double a, double b, XYZtriplet pin)
{
    double X = pin.x;
    double Y = pin.y;
    double Z = pin.z;
    double a2 = a*a;
    double f2  = a/b;
    double eb2 = f2*f2 - 1.0;
    double e2 = (a2 - b*b)/a2;
    double d  = eb2*b;
    //double ome2 = 1.0 - e2;
    double p, tu, tu2, su3, cu, cu3, tp, cp, sp, delta, tpold;

    XYZtriplet blh = new XYZtriplet();

    if(0.0==X && 0.0==Y)
    {
        blh.x = HALFPI;
        blh.y = 0.0;
        blh.z = Z - b;
        return blh;
    }
    else if(0.0==Z)
    {
        blh.x = 0.0;
        blh.y = Math.atan(Y/X);
        p = Math.sqrt(X*X + Y*Y);
        blh.z = p - a;
        return blh;
    }
    else
    {
        p = Math.sqrt(X*X + Y*Y);
        tp = f2*f2*Z/p;
        do
        {
            tpold = tp;
            tu  = tp/f2; //b*Z*(1.0 + d/r)/(a*p);
            tu2 = tu*tu;
            cu  = (1.0/Math.sqrt(1.0 + tu2));
            cu3 = cu*cu*cu;
            su3 = cu3*tu2*tu;
            tp  = (Z + d*su3)/(p - e2*a*cu3);
            delta = tp-tpold;
        }
        while(delta*delta > dtolerance);

        cp  = 1.0/Math.sqrt(1.0 + tp*tp);
        sp  = cp*tp;

        blh.x = Math.atan(tp);  //E = atan( (1-f_bessel)*tan(B) );
        blh.y = Math.atan(Y/X);
        blh.z = p*cp + Z*sp - a*Math.sqrt(1.0 - e2*sp*sp);

        return blh;
    }
} // bowring2

//---------------------------------------------------------------------------
// the basis is: (1+f')*X*tgE - b*e'^2*sinE - Y = 0
// Solving system y = x^4 + A*t^3 + B*t - 1 by Newton-Raphson ( t=tg(E/2) )
public static XYZtriplet nr2(double a, double b, XYZtriplet p)
{
    //double tolerance = 1.E-12;
    double step, x1, x1p2, y1, y1pr,
        f2, e22, f2R, X, Y, Z, R, A, B, C, tgEpocetno,
        tg, cs;

    XYZtriplet blh = new XYZtriplet();

    X = p.x;
    Y = p.y;
    Z = p.z;
    R = Math.sqrt(X*X + Y*Y);

    if(0.0==X && 0.0==Y)
    {
        blh.x = HALFPI;
        blh.y = 0.0;
        blh.z = Z - b;
        return blh;
    }
    else if(0.0==Z)
    {
        blh.x = 0.0;
        blh.y = Math.atan(Y/X);
        blh.z = R - a;
        return blh;
    }

    else
    {
        f2  = a/b;
        e22 = a*f2 - b;
        f2R = f2*R;

        A = 2.*(f2R + e22) / Z;
        B = 2.*(f2R - e22) / Z;
        C = 3.*A;

        tgEpocetno = f2*Z/R;
        x1 = tgEpocetno / (1.+Math.sqrt(1.+tgEpocetno*tgEpocetno));  // x=t=tan(E/2);
        x1p2 = x1*x1; // Die erste Iteration ist ausgeschlossen, weil Funktion könnte negativ sein (ausschließlich).

	y1 = x1*(B + x1p2*(A + x1)) - 1.;
        y1pr = x1p2*(C + 4.*x1) + B;
        step = y1/y1pr;

	x1 -= step;

        do
        {
            x1p2 = x1*x1; // Die erste Iteration ist ausgeschlossen, weil Funktion könnte negativ sein (ausschließlich).

            y1 = x1*(B + x1p2*(A + x1)) - 1.;
            y1pr = x1p2*(C + 4.*x1) + B;
            step = y1/y1pr;

            x1 -= step;
        }
        while (step > tolerance);

        tg = (2.*x1)/(1.-x1*x1) * f2;
        cs = 1. / Math.sqrt(1.+tg*tg);
        //t = (1.-t)/(1.+t); // t = tg(pi/4 - E/2) -->  t = tg(E/2)
        blh.x = Math.atan(tg);
        blh.y = Math.atan(Y/X);
        blh.z = (R - a*(1.-x1)/(1.+x1))*cs + (Z-b)*tg*cs;

        return blh;
    }
} // nr2

//---------------------------------------------------------------------------
// Solving system y = x^4 + A*t^3 + B*t - 1 by quadratic interpolation
// t = tg(E/2)
public static XYZtriplet sq2(double a, double b,  XYZtriplet p)
{
    double f2, e22, X, Y, Z, R, A, B, 
            x1, x1p2, y1, x2, x2p2, y2, oldx2, oldy,
            tg, cs;


    XYZtriplet blh = new XYZtriplet();

    f2  = a/b;
    e22 = a*f2 - b;
    X = p.x;
    Y = p.y;
    Z = p.z;
    R = Math.sqrt(X*X + Y*Y);

    A = 2.*(f2*R + e22) / Z;
    B = 2.*(f2*R - e22) / Z;

    x1 = Z/R * f2;  x1 = x1 / (1. + Math.sqrt(1.+x1*x1));
    x1p2 = x1*x1;   y1 = x1*(B + x1p2*(A + x1)) - 1.;

    x2 = x1 - y1/(x1p2*(3.*A + 4.*x1) + B);
    x2p2 = x2*x2;   y2 = x2*(B + x2p2*(A + x2)) - 1.;

    oldx2=x1p2;
    oldy=y1;
    
    do
    {
	x1p2 = (x2p2*y1 - x1p2*y2) / (y1-y2);
	x1   = Math.sqrt(x1p2);
	y1=x1*(B + x1p2*(A + x1)) - 1.;

	x2p2 = oldx2;
        y2 = oldy;
	oldx2 = x1p2;
        oldy = y1;
    }
    while( (y1>tolerance) | (y1<mtolerance) );

    tg = (2.*x1)/(1.-x1*x1) * f2;
    cs = 1. / Math.sqrt(1.+tg*tg);
    //t = (1.-t)/(1.+t); // t = tg(pi/4 - E/2) -->  t = tg(E/2)
    blh.x = Math.atan(tg);
    blh.y = Math.atan(Y/X);
    blh.z = (R - a*(1.-x1)/(1.+x1))*cs + (Z-b)*tg*cs;

    return blh;
} // sq2

//---------------------------------------------------------------------------
// the basis is: (1+f')*X*tgE - b*e'^2*sinE - Y = 0   <=>   a/b * X * tgE  - e22*sinE - Y
// Solving system y = x^4 + A*t^3 + B*t^2 + C*t + D by direct method  ( t=tg(E/2) )
public static XYZtriplet direct_solution(double a, double b,  XYZtriplet p)
{
    //double f2  = a/b;
    //double e22 = a*f2 - b;
    double X = p.x;
    double Y = p.y;
    double Z = p.z;
    double R = Math.sqrt(X*X + Y*Y);

    //double A = 2*(f2*R + e22) / Z;
    //double B = 2*(f2*R - e22) / Z;

    double f  = b/a;
    double sqd = -f * Z/R;
    double pom = (f*b-a)/R;

    double A = 2. * sqd;
    double C = A;
    double D = sqd*sqd;
    double B = 1 - pom*pom + D;

    double re[] = new double[4];
    double im[] = new double[2];

    solveP4(re, im, A, B, C, D);

    double t;
    XYZtriplet blh = new XYZtriplet();
    if(0.0==im[0])
        t = re[0]>0.0 ? re[0] : re[1];
    else
        t = re[2]>0.0 ? re[2] : re[3];

    ///*
    double sn = t/Math.sqrt(1. + t*t);
    //double cs = sn/x1;		x1=tg(E)

    double xe = a*sn/t;
    double ye = b*sn;

    double tgphi = t/f;
    double cosphi = 1.  / Math.sqrt(1.+tgphi*tgphi);

    blh.x = Math.atan(tgphi);
    blh.y = Math.atan(Y/X);
    blh.z = (R-xe)*cosphi + (Z-ye)*tgphi*cosphi;
    //*/

    return blh;
} //direct_solution

//*---------------------------------------------------------------------------
/// the basis is: (1+f')*X*tgE - b*e'^2*sinE - Y = 0
// Solving system y = t^4 + A*t^3 + B*t - 1 by direct method
public static XYZtriplet direct_solution_opt(double a, double b,  XYZtriplet p)
{
	double  f2, e22, X, Y, Z, R,
            aa, c, q1, q2, p1, p2,
            D, sqd, y, q, r, A, B,
            t, t2, tg, sn, xe, ye, tgphi, cosphi;

    f2  = a/b;
    e22 = a*f2 - b;
    X = p.x;
    Y = p.y;
    Z = p.z;
    R = Math.sqrt(X*X + Y*Y);

    aa = 2*(f2*R + e22) / Z;
    c  = 2*(f2*R - e22) / Z;

    q  = -(aa*c + 4.)/3.;
    r  =  (aa*aa - c*c)*0.5;
    A  = -Math.pow(Math.abs(r)+Math.sqrt(r*r-q*q*q),1./3.);
    //A = -exp(log(Math.abs(r)+sqrt(r*r-q*q*q)) * ONETHIRD);
    if( r<0 ) A=-A;
    B = (0==A ? 0 : q/A);
    y = A+B;

    D = y*y + 4.;
    sqd = Math.sqrt(D);
    q1 = (y + sqd) * 0.5;
    q2 = (y - sqd) * 0.5;

    p1 = (aa*q1-c)/(q1-q2);

	D = p1*p1 - 4*q1;
	if(D >= 0.0)
        t = (-p1 + Math.sqrt(D)) * 0.5;
	else
	{
        p2 = (c-aa*q2)/(q1-q2);
        t = (-p2 + Math.sqrt(p2*p2 - 4*q2)) * 0.5; // D2 = p2*p2 - 4*q2;
    }

    XYZtriplet blh = new XYZtriplet();
    t2 = t*t;
    t *= 2;
    tg = t/(1-t2);
    sn = t/(1+t2);

    xe = a*sn/tg;
    ye = b*sn;

    tgphi = tg*f2;
    blh.x = Math.atan(tgphi);			//E = atan( (1-f_bessel)*tan(B) );
    blh.y = Math.atan(Y/X);
    cosphi = 1 / Math.sqrt(1+tgphi*tgphi);
    blh.z = (R - xe)*cosphi + (Z-ye)*tgphi*cosphi;

    return blh;
}


//===========================================================================
//===========================================================================
    public static void main(String[] args) 
    {
	XYZtriplet blh    = new XYZtriplet();
        XYZtriplet xyz    = new XYZtriplet();
	XYZtriplet blhret = new XYZtriplet();

	/// Ртањ
	double B = 43.7761 * PI/180.0;
	double L = 21.8933 * PI/180.0;
	//double H = 2170.0;
	//double H = 400000.0; // Dove
	double H = 36000000.0; // BeiDou-1
	//double H = -5157000.0;

	blh.x = B;
	blh.y = L;
	blh.z = H;

	xyz = BLh2XYZ(a_bessel, b_bessel, blh);
        //System.out.println("x=" + xyz.x + "  y=" + xyz.y + "  z=" + xyz.z);
        
	System.out.println("Bowring-proj test:");
	blhret = geocentric_to_geodetic_bowring(a_bessel, b_bessel, xyz);
	System.out.println(blhret.x - B);
	System.out.println(blhret.y - L);
	System.out.println(blhret.z - H);
	System.out.println();

        System.out.println("Moritz-Heiskenen test:");
	blhret = geocentric_to_geodetic(a_bessel, b_bessel, xyz);
	System.out.println(blhret.x - B);
	System.out.println(blhret.y - L);
	System.out.println(blhret.z - H);
	System.out.println();

        System.out.println("Bowring2 test:");
	blhret = bowring2(a_bessel, b_bessel, xyz);
	System.out.println(blhret.x - B);
	System.out.println(blhret.y - L);
	System.out.println(blhret.z - H);
	System.out.println();

        System.out.println("Bowring3 test:");
	blhret = bowring3(a_bessel, b_bessel, xyz);
	System.out.println(blhret.x - B);
	System.out.println(blhret.y - L);
	System.out.println(blhret.z - H);
	System.out.println();

        System.out.println("Newton-Raphson test:");
	blhret = nr2(a_bessel, b_bessel, xyz);
	System.out.println(blhret.x - B);
	System.out.println(blhret.y - L);
	System.out.println(blhret.z - H);
	System.out.println();
        
        System.out.println("Semiquadratic test:");
	blhret = sq2(a_bessel, b_bessel, xyz);
	System.out.println(blhret.x - B);
	System.out.println(blhret.y - L);
	System.out.println(blhret.z - H);
	System.out.println();
        
        System.out.println("Direct solution test:");
	blhret = direct_solution(a_bessel, b_bessel, xyz);
	System.out.println(blhret.x - B);
	System.out.println(blhret.y - L);
	System.out.println(blhret.z - H);
	System.out.println();
        
        System.out.println("Direct solution optimized test:");
	blhret = direct_solution_opt(a_bessel, b_bessel, xyz);
	System.out.println(blhret.x - B);
	System.out.println(blhret.y - L);
	System.out.println(blhret.z - H);
	System.out.println();
        
        //------------------------------
        //------------------------------
        long t = System.nanoTime();
	for(double fi=0.1; fi<90.0; fi += 0.1)
        {
            //for(double h=0.0; h<10000.0; h+=1.0)
            //for(double h=395000; h<405000; h+=1.0)
            for(double h=35779000; h<35789000; h+=1.0)
            {
            	blh.x = fi * PI/180.0;
               	blh.z = h;
		xyz = BLh2XYZ(a_bessel, b_bessel, blh);

		blhret = geocentric_to_geodetic_bowring(a_bessel, b_bessel, xyz);
            }
        }
	t = System.nanoTime() - t;
        System.out.println("Bowring proj:  " + t*1e-9 + " sec (9*10^5 executions)");

        //------------------------------
        t = System.nanoTime();
	for(double fi=0.1; fi<90.0; fi += 0.1)
        {
            //for(double h=0.0; h<10000.0; h+=1.0)
            //for(double h=395000; h<405000; h+=1.0)
            for(double h=35779000; h<35789000; h+=1.0)
            {
            	blh.x = fi * PI/180.0;
               	blh.z = h;
		xyz = BLh2XYZ(a_bessel, b_bessel, blh);

		blhret = geocentric_to_geodetic(a_bessel, b_bessel, xyz);
            }
        }
	t = System.nanoTime() - t;
        System.out.println("Moritz-Heiskenen:  " + t*1e-9);

        //------------------------------
        t = System.nanoTime();
	for(double fi=0.1; fi<90.0; fi += 0.1)
        {
            //for(double h=0.0; h<10000.0; h+=1.0)
            //for(double h=395000; h<405000; h+=1.0)
            for(double h=35779000; h<35789000; h+=1.0)
            {
            	blh.x = fi * PI/180.0;
               	blh.z = h;
		xyz = BLh2XYZ(a_bessel, b_bessel, blh);

		blhret = bowring2(a_bessel, b_bessel, xyz);
            }
        }
	t = System.nanoTime() - t;
        System.out.println("bowring2:  " + t*1e-9);

        //------------------------------
        t = System.nanoTime();
	for(double fi=0.1; fi<90.0; fi += 0.1)
        {
            //for(double h=0.0; h<10000.0; h+=1.0)
            //for(double h=395000; h<405000; h+=1.0)
            for(double h=35779000; h<35789000; h+=1.0)
            {
            	blh.x = fi * PI/180.0;
               	blh.z = h;
		xyz = BLh2XYZ(a_bessel, b_bessel, blh);

		blhret = bowring3(a_bessel, b_bessel, xyz);
            }
        }
	t = System.nanoTime() - t;
        System.out.println("bowring3:  " + t*1e-9);

        //------------------------------
        t = System.nanoTime();
	for(double fi=0.1; fi<90.0; fi += 0.1)
        {
            //for(double h=0.0; h<10000.0; h+=1.0)
            //for(double h=395000; h<405000; h+=1.0)
            for(double h=35779000; h<35789000; h+=1.0)
            {
            	blh.x = fi * PI/180.0;
               	blh.z = h;
		xyz = BLh2XYZ(a_bessel, b_bessel, blh);

		blhret = nr2(a_bessel, b_bessel, xyz);
            }
        }
	t = System.nanoTime() - t;
        System.out.println("nr2:  " + t*1e-9);

        //------------------------------
        t = System.nanoTime();
	for(double fi=0.1; fi<90.0; fi += 0.1)
        {
            //for(double h=0.0; h<10000.0; h+=1.0)
            //for(double h=395000; h<405000; h+=1.0)
            for(double h=35779000; h<35789000; h+=1.0)
            {
            	blh.x = fi * PI/180.0;
               	blh.z = h;
		xyz = BLh2XYZ(a_bessel, b_bessel, blh);

		blhret = sq2(a_bessel, b_bessel, xyz);
            }
        }
	t = System.nanoTime() - t;
        System.out.println("sq2:  " + t*1e-9);

        //------------------------------
        t = System.nanoTime();
	for(double fi=0.1; fi<90.0; fi += 0.1)
        {
            //for(double h=0.0; h<10000.0; h+=1.0)
            //for(double h=395000; h<405000; h+=1.0)
            for(double h=35779000; h<35789000; h+=1.0)
            {
            	blh.x = fi * PI/180.0;
               	blh.z = h;
		xyz = BLh2XYZ(a_bessel, b_bessel, blh);

		blhret = direct_solution(a_bessel, b_bessel, xyz);
            }
        }
	t = System.nanoTime() - t;
        System.out.println("direct_solution:  " + t*1e-9);

        //------------------------------
        t = System.nanoTime();
	for(double fi=0.1; fi<90.0; fi += 0.1)
        {
            //for(double h=0.0; h<10000.0; h+=1.0)
            //for(double h=395000; h<405000; h+=1.0)
            for(double h=35779000; h<35789000; h+=1.0)
            {
            	blh.x = fi * PI/180.0;
               	blh.z = h;
		xyz = BLh2XYZ(a_bessel, b_bessel, blh);

		blhret = direct_solution_opt(a_bessel, b_bessel, xyz);
            }
        }
	t = System.nanoTime() - t;
        System.out.println("direct_solution_opt:  " + t*1e-9);
    }  
}
