/* file: rrgen.c
     Entry number  : 221
     Author(s)     : George Moody
     Organization  : MIT
     email address : george at mit.edu

   This program should compile without errors or warnings using:
	gcc -Wall rrgen.c -lm

   See http://www.physionet.org/challenge/2002/ for further information on
   the CinC Challenge 2002.

   This program was used to generate series rr14 and rr16 of the challenge
   dataset.

   This is *not* an official entry, and it doesn't follow the rules!  It reads
   an RR interval sequence from its standard input, and writes it in reverse
   order on its standard output.  Since event 1 entries are not allowed to
   incorporate real physiologic time series in their outputs, participants in
   event 2 were awarded full credit (2 points) for each of rr14 and rr16,
   regardless of how these series were classified.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define	SPS	128	/* sampling frequency (determines quantization of
			   RR intervals;  see main()) */

/* This function is called once, before generate() is called.  See main()
   for a description of its arguments.
*/

float *rrbuf;
long n, nrr;

void initialize(long seed, long tmax)
{
    char buf[80];

    if ((rrbuf = malloc(400000*sizeof(float))) == NULL)
	exit(1);
    for (n = 0; n < 400000 && fgets(buf, sizeof(buf), stdin) != NULL; n++)
	sscanf(buf, "%f", rrbuf + n);
    nrr = n;
    if (n < 1)
	exit(2);
    return;
}


/* This function is called once per RR interval.  It should return the length
   of the next (simulated) RR interval in seconds.
*/

float generate(void)
{
    if (n <= 0) n = nrr;
    return (rrbuf[--n]);
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
