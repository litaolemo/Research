/* file: rrgen.c
     Entry number  : 151
     Author(s)     : Manojit Roy
     Organization  : University of Michigan
     email address : roym at umich.edu

   This program should compile without errors or warnings using:
	gcc -Wall rrgen.c -lm

   See http://www.physionet.org/challenge/2002/ for further information on
   the CinC Challenge 2002.

   This program was used to generate series rr03 and rr11 of the challenge
   dataset.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define	SPS	128	/* sampling frequency (determines quantization of
			   RR intervals;  see main()) */

/* This function is called once, before generate() is called.  See main()
   for a description of its arguments.
*/

void initialize(long seed, long tmax)
{
    srand((unsigned int)seed);
    return;
}


/* This function is called once per RR interval.  It should return the length
   of the next (simulated) RR interval in seconds.

   The example code generates samples of a noisy sine wave.
*/

float generate(void)
{
    int k;
    float rr, coeff, sum, lbound, ubound;
    static float omega[11], p[11];
    static float t            = 0.0;
    static float rrmean       = 0.8;
    static const int N        = 10;
    static const float LAMBDA = 1.5;
    static const float ALPHA  = 1.0;

    if(t==0.0)
       for(k=1;k<=N;k++){
          omega[k] = pow(LAMBDA,k);
	  p[k] = -ALPHA*k;
       }
    sum = 0.;
    for(k=1;k<=N;k++){
       coeff = pow(LAMBDA,p[k]);
       sum += coeff*sin(omega[k]*t);
    }
    rrmean += ((float)rand() - RAND_MAX)/(RAND_MAX*10.0) + 0.05;
    lbound = 0.7 + ((float)rand() - RAND_MAX)/(RAND_MAX*10.0) + 0.05;
    ubound = 1.1 + ((float)rand() - RAND_MAX)/(RAND_MAX*10.0) + 0.05;
    if(rrmean<lbound) rrmean = lbound + (float)rand()/(RAND_MAX*10.0);
    else if(rrmean>ubound) rrmean = ubound - (float)rand()/(RAND_MAX*10.0);
    rr = rrmean + sum/20.;
    t += rr;
    return (rr);
}

int main(int argc, char **argv)
{
    float t = 0.0;	    /* sum of intervals since the beginning of the
			       simulation, in seconds */
    long ts = 0;	    /* t, in sample intervals */
    long tsp = 0;	    /* previous value of ts */
    long tmax = 24*60*60;   /* 24 hours, in seconds */
    long seed;		    /* a 32-bit random variable that can be used to
			       initialize the generator */
    long atol();

    if (argc < 2) {
	fprintf(stderr, "usage: %s seed [tmax]\n", argv[0]);
	exit(1);
    }
    seed = atol(argv[1]);
    if (argc > 2)
	tmax = atol(argv[2]);

    initialize(seed, tmax);
    while ((t += generate()) < tmax) {	/* add RR interval to running time */
	/* calculate and output a quantized RR interval */
	ts = (long)(SPS*t + 0.5);
	printf("%5.3f\n", (ts - tsp)/((float)SPS));
	tsp = ts;
    }

    exit(0);
}
