/*
 * $Id: Clausen.cpp,v 1.2 2008/07/16 07:05:09 ps Exp $
 * $Source: /cs/research/multires/cvsroot/src/ConfEquiv/Clausen.cpp,v $
 */
#include "Clausen.h"
#include <math.h>

/*
 * C version of r1mach and d1mach from core (netlib)
 * specialized for IEEE arithmetic
 *
 * MACHINE CONSTANTS (s for single, d for double)
 * {S|D}1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
 * {S|D}1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
 * {S|D}1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
 * {S|D}1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
 * {S|D}1MACH(5) = LOG10(B)
 */
static double mach[5] = {2.2250738585072014e-308, 1.7976931348623157e+308,
                         1.1102230246251565e-16, 2.2204460492503131e-16,
                         3.0102999566398120e-01};

/*
 * chebyshev expansion of `cl(t)/t + log(t)' around Pi/6
 * accurate to 20 decimal places
 */
static double clpi6[14] = {
    2 * 1.0057346496467363858, .0076523796971586786263,
    .0019223823523180480014,   .53333368801173950429e-5,
    .68684944849366102659e-6,  .63769755654413855855e-8,
    .57069363812137970721e-9,  .87936343137236194448e-11,
    .62365831120408524691e-12, .12996625954032513221e-13,
    .78762044080566097484e-15, .20080243561666612900e-16,
    .10916495826127475499e-17, .32027217200949691956e-19};

/*
 * chebyshev expansion of cl(t/2)/t + log(t)' around Pi/2
 * accurate to 20 decimal places
 */
static double clpi2[19] = {
    2 * .017492908851746863924 + 2 * 1.0057346496467363858,
    .023421240075284860656 + .0076523796971586786263,
    .0060025281630108248332 + .0019223823523180480014,
    .000085934211448718844330 + .53333368801173950429e-5,
    .000012155033501044820317 + .68684944849366102659e-6,
    .46587486310623464413e-6 + .63769755654413855855e-8,
    .50732554559130493329e-7 + .57069363812137970721e-9,
    .28794458754760053792e-8 + .87936343137236194448e-11,
    .27792370776596244150e-9 + .62365831120408524691e-12,
    .19340423475636663004e-10 + .12996625954032513221e-13,
    .17726134256574610202e-11 + .78762044080566097484e-15,
    .13811355237660945692e-12 + .20080243561666612900e-16,
    .12433074161771699487e-13 + .10916495826127475499e-17,
    .10342683357723940535e-14 + .32027217200949691956e-19,
    .92910354101990447850e-16,
    .80428334724548559541e-17,
    .72598441354406482972e-18,
    .64475701884829384587e-19,
    .58630185185185185187e-20};

/*
 * chebyshev expansion of `-cl(Pi-t)/(Pi-t) + log(2)' around 5Pi/6
 * accurate to 20 decimal places
 */
static double cl5pi6[19] = {
    2 * .017492908851746863924, .023421240075284860656,
    .0060025281630108248332,    .000085934211448718844330,
    .000012155033501044820317,  .46587486310623464413e-6,
    .50732554559130493329e-7,   .28794458754760053792e-8,
    .27792370776596244150e-9,   .19340423475636663004e-10,
    .17726134256574610202e-11,  .13811355237660945692e-12,
    .12433074161771699487e-13,  .10342683357723940535e-14,
    .92910354101990447850e-16,  .80428334724548559541e-17,
    .72598441354406482972e-18,  .64475701884829384587e-19,
    .58630185185185185187e-20};

/*
 * evaluate a chebyshev series
 * adapted from fortran csevl
 */
double csevl(double x, double* cs, int n) {
    double b2, b1 = 0, b0 = 0, twox = 2 * x;

    while (n--) {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + cs[n];
    }

    return (b0 - b2) / 2;
}

/*
 * from the original fortran inits
 * april 1977 version.  w. fullerton, c3, los alamos scientific lab.
 *
 * initialize the orthogonal series so that inits is the number of terms
 * needed to insure the error is no larger than eta.  ordinarily, eta
 * will be chosen to be one-tenth machine precision.
 */
static int inits(double* series, int n, double eta) {
    double err = 0;

    while (err <= eta && n--) {
        err += fabs(series[n]);
    }

    return n++;
}

double clausen(double x) {
    static int nclpi6 = 0, nclpi2 = 0, ncl5pi6 = 0;
    /*
     * right half (Pi <= x < 2 Pi)
     */
    int rh = 0;
    double f;

    if (!nclpi6) {
        nclpi6  = inits(clpi6, sizeof clpi6 / sizeof *clpi6, mach[2] / 10);
        nclpi2  = inits(clpi2, sizeof clpi2 / sizeof *clpi2, mach[2] / 10);
        ncl5pi6 = inits(cl5pi6, sizeof cl5pi6 / sizeof *cl5pi6, mach[2] / 10);
    }

    /*
     * get to canonical interval
     */
    if ((x = fmod(x, 2 * M_PI)) < 0) {
        x += (double)(2 * M_PI);
    }
    if (x > (double)M_PI) {
        rh = 1;
        x  = (double)(2 * M_PI) - x;
    }

    if (x == 0) {
        f = x;
    } else if (x <= (double)(M_PI / 3)) {
        f = csevl(x * (double)(6 / M_PI) - 1, clpi6, nclpi6) * x - x * log(x);
    } else if (x <= (double)(2 * M_PI / 3)) {
        f = csevl(x * (double)(3 / M_PI) - 1, clpi2, nclpi2) * x - x * log(x);
    } else { /* x <= Pi */
        f = ((double)M_LN2 -
             csevl(5 - x * (double)(6 / M_PI), cl5pi6, ncl5pi6)) *
            ((double)M_PI - x);
    }

    return rh ? -f : f;
}
