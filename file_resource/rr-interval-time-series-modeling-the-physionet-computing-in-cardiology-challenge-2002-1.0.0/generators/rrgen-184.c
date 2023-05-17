/* file: rrgen.c
     Entry number  : 184
     Author(s)     : Phil Langley, Emma Bowers, Michael Drinnan, John Allan,
		      Alan Murray
     Organization  : Medical Physics, Freeman Hospital, Newcastle upon Tyne, UK
     email address : philip.langley at ncl.ac.uk,
		      emma.bowers at nuth.northy.nhs.uk

   This program should compile without errors or warnings using:
	gcc -Wall rrgen.c -lm

   See http://www.physionet.org/challenge/2002/ for further information on
   the CinC Challenge 2002.

   This program was used to generate series rr17 and rr25 of the challenge
   dataset.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define	SPS	128	/* sampling frequency (determines quantization of
			   RR intervals;  see main()) */



#define cSampleRate 5
#define cDT (1.0 / (double) cSampleRate)
#define cSamplesPerDay (24 * 60 * 60 * cSampleRate)
#define cMaxDays 2
#define cSamplesTotal (cSamplesPerDay * cMaxDays)
#define cMaxAlloc (cSamplesTotal + 1000)
#define cMaxCycles 1000
#define cPI (3.141592654)

#define cDayActiveAmpMaxPos 0.15
#define cDayActiveAmpMinPos 0.05
#define cDayActiveAmpVarPos 0.25
#define cDayActiveFrqMaxPos 0.3
#define cDayActiveFrqMinPos 0.1666
#define cDayActiveFrqVarPos 0.15

#define cDayActiveAmpMaxNeg 0.15
#define cDayActiveAmpMinNeg 0.05
#define cDayActiveAmpVarNeg 0.25
#define cDayActiveFrqMaxNeg 0.3
#define cDayActiveFrqMinNeg 0.1666
#define cDayActiveFrqVarNeg 0.15

#define cDayPassiveAmpMaxPos 0.2
#define cDayPassiveAmpMinPos 0.05
#define cDayPassiveAmpVarPos 0.05
#define cDayPassiveFrqMaxPos 0.2666
#define cDayPassiveFrqMinPos 0.0666
#define cDayPassiveFrqVarPos 1.0

#define cDayPassiveAmpMaxNeg 0.3
#define cDayPassiveAmpMinNeg 0.05
#define cDayPassiveAmpVarNeg 0.15
#define cDayPassiveFrqMaxNeg 0.2666
#define cDayPassiveFrqMinNeg 0.0666
#define cDayPassiveFrqVarNeg 1.0

#define cNightActiveAmpMaxPos 0.3
#define cNightActiveAmpMinPos 0.05
#define cNightActiveAmpVarPos 0.05
#define cNightActiveFrqMaxPos 0.2666
#define cNightActiveFrqMinPos 0.0666
#define cNightActiveFrqVarPos 1.0

#define cNightActiveAmpMaxNeg 0.3
#define cNightActiveAmpMinNeg 0.05
#define cNightActiveAmpVarNeg 0.05
#define cNightActiveFrqMaxNeg 0.2666
#define cNightActiveFrqMinNeg 0.0666
#define cNightActiveFrqVarNeg 1.0

#define cNightPassiveAmpMaxPos 0.2
#define cNightPassiveAmpMinPos 0.05
#define cNightPassiveAmpVarPos 0.1
#define cNightPassiveFrqMaxPos 0.2
#define cNightPassiveFrqMinPos 0.0666
#define cNightPassiveFrqVarPos 0.15

#define cNightPassiveAmpMaxNeg 0.2
#define cNightPassiveAmpMinNeg 0.05
#define cNightPassiveAmpVarNeg 0.1
#define cNightPassiveFrqMaxNeg 0.2
#define cNightPassiveFrqMinNeg 0.0666
#define cNightPassiveFrqVarNeg 0.15

double *nT = NULL ;
double *nRR = NULL ;
double *nLFeffect = NULL ;
double *nWhitenoise = NULL ;
double *nPinknoise = NULL ;
double *nResp = NULL ;


double Min( double A, double B )
	{
   return A < B ? A : B ;
   }


double Max( double A, double B )
	{
   return A > B ? A : B ;
   }


/* Set part of the array (IndexFrom to IndexTo) to equal Value */
void SetArray( double *Array, int IndexFrom, int IndexTo, double Value )
	{
   int Index ;
   if( ( IndexFrom < 0 ) || ( IndexTo >= cSamplesPerDay ) || ( IndexFrom > IndexTo ) )
   	fprintf( stderr, "Can't write to array elements %d to %d\n", IndexFrom, IndexTo ) ;
   else
   	for( Index = IndexFrom ; Index <= IndexTo ; Index ++ )
      	Array[ Index ] = Value ;
   }


/* Pick a random number between 0 and 1 */
/* rand( ) picks a number between 1 and RAND_MAX */
/* This is scaled down appropriately */
/* Numbers are from a flat distribution */
double Rand0To1( void )
	{
   double Val ;
   Val = rand( ) ;
   return Val / (double) RAND_MAX ;
   }



/* Pick a random number from a normal distribution */
/* Approximation to normal distribution works by adding 100 random values from a flat distribution */
/* The value is then scaled for mean of 0 and SD of 1 */
/* The factor 3.464 is square root of 12 */
double RandGauss( void )
	{
   int Count ;
   double Total ;

   Total = 0.0 ;
   for( Count = 0 ; Count < 100 ; Count ++ )
      Total += (double) rand( ) ;
   Total -= 50.0 * (double) RAND_MAX ;
   Total /= 10.0 * (double) RAND_MAX ;
   Total *= 3.464 ;
   return Total ;
   }


/* Convert from Time to the corresponding Index into the array */
int TimeToIndex( double Time )
	{
   double Index ;
   Index = (double ) Time * (double) cSampleRate ;
   return (int) Index ;
   }

/* Convert from Index to the corresponding Time */
double IndexToTime( int Index )
	{
   double Time ;
   Time = (double ) Index / (double) cSampleRate ;
   return Time ;
   }

/*
%   a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
%                         - a(2)*y(n-1) - ... - a(na+1)*y(n-na)
*/
void Filter( double *X, int N, double *a, double *b, int Order, double *Y )
	{
   int SampleNo ;
   int SourceSampleNo ;
   int TapNo ;
   double Output ;
	for( SampleNo = 0 ; SampleNo < N ; SampleNo ++ )
   	{
      Output = b[ 0 ] * X[ SampleNo ] ;
      for( TapNo = 1 ; TapNo <= Order ; TapNo ++ )
      	{
         SourceSampleNo = SampleNo - TapNo ;
         if( SourceSampleNo < 0 )
         	break ;
         Output += b[ TapNo ] * X[ SourceSampleNo ] - a[ TapNo ] * Y[ SourceSampleNo ] ;
      	}
      Y[ SampleNo ] = Output ;
      }
   }



/* Implement a first-order low-pass filter for 1/f rolloff */
/* Idea is first to calculate the corresponding time constant */
/* This leads to multiplication factor for exponential smoothing */
/* This is equivalent to a 1st order RC filter */
/* General idea is: NewOutput  =  ( Z * OldOutput )  +  ( ( 1 - Z ) * CurrentInput ) */
/* If Z is nearly 1, the time constant is long */
void FirstOrderLowPassFilter( double *X, int N, double CutoffFreq, double *Y )
	{
   double TimeConstant ;
	double DecayPerSec ;
	double DecayPerSample ;
	double OneMinusDecayPerSample ;
	double Output ;
   int SampleNo ;
   TimeConstant = 1.0 / ( 2.0 * cPI * CutoffFreq ) ;
	DecayPerSec = 1.0 / TimeConstant ;
   DecayPerSample = DecayPerSec / (double) cSampleRate ;
   OneMinusDecayPerSample = 1.0 - DecayPerSample ;
   Output = 0.0 ;
	for( SampleNo = 1 ; SampleNo <= N ; SampleNo ++ )
   	{
      Output = ( OneMinusDecayPerSample * Output ) + ( DecayPerSample * X[ SampleNo ] ) ;
      Y[ SampleNo ] = Output ;
      }
   }



/* Sort the values in the Raw array */
/* Number of values is given by N */
/* Results are put in Cooked array */
/* Uses simple insertion sort algorithm */
/* Go thru Raw array N times */
/* Each time thru, find the smallest element */
/* Add this element to the end of the Cooked array */
/* Then set the element to a very big number so it doesn't get found again */
void Sort( double *Raw, int N, double *Cooked )
	{
   int Iteration ;
   int Smallest ;
   int Test ;
   for( Iteration = 0 ; Iteration <= N - 1 ; Iteration ++ )
   	{
      Smallest = 1 ;
      for( Test = 1 ; Test <= N ; Test ++ )
      	if( Raw[ Test ] < Raw[ Smallest ] )
         	Smallest = Test ;
      Cooked[ Iteration ] = Raw[ Smallest ] ;
      Raw[ Smallest ] = 1e250 ;
      }
   }



/* This section adapted from MATLAB code by Phil & Emma */
/*
%development code for Computers in Cardiology Challenge 2002
%24 hr RR simulator for normal subjects

%CVPER group
%16th April 2002
*/




/*Calculating the frequency and amplitude of the respiratory component*/
void choosefreqamp( double maxamp, double minamp, double percentagechangeamp, double previousamp, double maxfreq, double minfreq, double variance, double tlength, double i, double *numberoftimeelement, double *freq, double *amp )
	{
   double meanfreq ;

	/*choosing the amplitude randomly within a percentage range of the previous value*/
   *amp = ( Rand0To1( ) * 2.0 * percentagechangeamp * previousamp ) + ( ( 1.0 - percentagechangeamp ) * previousamp ) ;
	/*if the above value does not fall in the given range guess again */
	*amp = Min( *amp, maxamp + Rand0To1( ) * 0.02 - 0.01 ) ;
	*amp = Max( *amp, minamp + Rand0To1( ) * 0.02 - 0.01 ) ;

	/*choosing the frequency randomly*/
	meanfreq = ( maxfreq + minfreq ) / 2.0 ;
   *freq = RandGauss( ) * variance + meanfreq ;
	/*if the above value does not fall in the given range continue to guess until it does*/
   while( ( *freq < minfreq ) || ( *freq > maxfreq ) )
	   *freq = RandGauss( ) * variance + meanfreq ;

	/*working out how many whole array elements are in a half  period*/
   *numberoftimeelement = floor( 1.0 / ( 2.0 * *freq * cDT ) ) ;
   if( i + *numberoftimeelement > tlength )
      *numberoftimeelement = tlength - i ;
   }


/* Note this only generates 1 day's worth of RR data */
void RRSim1Day( double RRday, double RRnight, double *T, double *RR, double *Resp )
	{
   int Index ;
   int SleepIndex ;
   int WakeIndex ;
  /* int NumRemCycles ;*/
   int NumFalls1 ;
   int NumFalls2 ;
   int NumArouse ;
   int NumArtifacts ;
/*   int Cycle ;*/
   int Fall ;
   int Arouse ;
   int Artifact ;
/*   int RemIndex ;*/
   int FallIndex1 ;
   int FallIndex2 ;
   int ArouseIndex ;
   int i ;
   int IsSleep ;
   int IsFall1 ;
   int IsFall2 ;
   int IsArousal ;
/*   int IsREM ;*/
	double SleepTime ;
	double WakeTime ;
/*	double RemTime ;*/
   double ArouseTime ;
   double RampDuration ;
   double Duration ;
   double FallTime ;
   double FallRR ;
/*   double RemRR ;*/
   double ArouseRR ;
   double FreqWander ;
   double RandInd ;
   double amp ;
	double numberoftimeelement ;
 /*  double CoeffA[ ] = {  1.0000,   -4.0469,    6.3705,   -4.8224,    1.7212,   -0.2223 } ;
   double CoeffB[ ] = {  1.0000,   -3.4673,    4.4317,   -2.4428,    0.4608,    0.0177 } ;*/
/*	double RemDuration[ cMaxCycles ] ;
	double StableDuration[ cMaxCycles ] ;
	double RemFrom[ cMaxCycles ] ;
	double RemTo[ cMaxCycles ] ;*/
	double FallFrom1[ cMaxCycles ] ;
	double FallTo1[ cMaxCycles ] ;
	double FallFrom2[ cMaxCycles ] ;
	double FallTo2[ cMaxCycles ] ;
	double ArouseFrom[ cMaxCycles ] ;
	double ArouseTo[ cMaxCycles ] ;
	double OrderedFallFrom1[ cMaxCycles ] ;
	double OrderedFallTo1[ cMaxCycles ] ;
	double OrderedArouseFrom[ cMaxCycles ] ;
	double OrderedArouseTo[ cMaxCycles ] ;
	double OrderedFallFrom2[ cMaxCycles ] ;
	double OrderedFallTo2[ cMaxCycles ] ;
   double freqpos ;
   double freqneg ;
   double amppos;
   double ampneg ;
   static double previousamppos = 0.1 ;
   static double previousampneg = 0.1 ;

	/*
	%define time axis
	t = [0:dt:24*60*60];
	*/
	for( Index = 0 ; Index < cSamplesPerDay ; Index ++ )
   	T[ Index ] = (double) Index * cDT ;

	/*
	%t = 0 is start of recording

	%Define sleep and awake times
	%sleep 12 to 14 hrs after start time
	sleepTime = (rand*2+12)*60*60;
	*/
	SleepTime = ( Rand0To1( ) * 2.0 + 12.0 ) * 60.0 * 60.0 ;

	/*
	%wake up between 2 and 3 hrs before end of recording next day
	wakeTime = (24-(rand*1+2))*60*60;
	*/
	WakeTime = ( 24.0 - ( Rand0To1( ) * 1.0 + 2.0 ) ) * 60.0 * 60.0 ;

   /*
	RR = ones(size(t))*RRday;
	sleepIndex = find(t > (sleepTime-dt/2) & t < (sleepTime+dt/2));
	wakeIndex = find(t > (wakeTime-dt/2) & t < (wakeTime+dt/2));
	RR(sleepIndex(1):wakeIndex(1)) = RRnight;
   */
	SleepIndex = TimeToIndex( SleepTime ) ;
   WakeIndex = TimeToIndex( WakeTime ) ;
	for( Index = 0 ; Index < cSamplesPerDay ; Index ++ )
   	RR[ Index ] = RRday ;
	for( Index = SleepIndex ; Index <= WakeIndex ; Index ++ )
   	RR[ Index ] = RRnight ;



	/*
	%ramp up from day to night rr over 30 min period
	%y = mx+c
	if(1)
	rampDuration = 30*60;
	RR(sleepIndex(1):sleepIndex(1)+rampDuration/dt) = ((RRnight-RRday)/(rampDuration/dt))*[0:rampDuration/dt] + RRday;
	end
   */
/*   if( 1 )
   	{ */
/*      fprintf( stderr, "  Day-night ramp...\n" ) ;*/
      RampDuration = (Rand0To1( ) * 50.0 + 40.0 ) * 60.0 ;
		for( Index = 0 ; Index <= RampDuration / cDT ; Index ++ )
   		RR[ SleepIndex + Index ] = ( ( RRnight - RRday ) / ( RampDuration / cDT ) ) * Index + RRday ;
/*      } */
		RampDuration = (Rand0To1( ) * 15.0 + 15.0 ) * 60.0 ;
      for( Index = 0 ; Index <= RampDuration / cDT ; Index ++ )
   		RR[ WakeIndex + Index ] = ( ( RRday - RRnight ) / ( RampDuration / cDT ) ) * Index + RRnight ;

   /*
	%cycle of 50 to 70 mins of stable sleep RR (RRnight) followed by 15 to 45 mins unstable sleep RR (REM sleep)
	%reduces mean RR during REM sleep by between 5 to 20 %
	if (1)
	numCycles = fix((wakeTime-sleepTime)/(1.5*60*60));
	for (cycle = 1:numCycles)
   	remDuration(cycle) = (rand*30+15)*60;
	   stableDuration(cycle) = (rand*20+50)*60;
   	if (cycle == 1)
      	remTime = sleepTime+stableDuration(cycle);
	   else
   	   remTime = remTime+stableDuration(cycle)+remDuration(cycle-1);
	   end
	   remRR = RRnight*(1-(rand*0.15+0.05));
   	remInd(cycle).vect = find(t > remTime & t < remTime+remDuration(cycle));
	   RR(remInd(cycle).vect) = ones(size(remInd(cycle).vect))*remRR;
	end
	end
	*/
 /*	if( 1 )
   	{ */
      // fprintf( stderr, "  REM cycles...\n" ) ;
/*		NumRemCycles = (int) floor( ( WakeTime - SleepTime ) / ( 1.5 * 60.0 * 60.0 ) ) ;
		for( Cycle = 1 ; Cycle <= NumRemCycles ; Cycle ++ )
      	{
	   	RemDuration[ Cycle ] = ( Rand0To1( ) * 30.0 + 15.0 ) * 60.0 ;
		   StableDuration[ Cycle ] = ( Rand0To1( ) * 40.0 + 40.0 ) * 60.0 ;
   		if( Cycle == 1 )
      		RemTime = SleepTime + StableDuration[ Cycle ] ;
		   else
   		   RemTime = RemTime + StableDuration[ Cycle ] + RemDuration[ Cycle - 1 ] ;
	   	RemRR = RRnight * ( 1 - ( Rand0To1( ) * 0.03 +0.02 ) ) ;
         RemFrom[ Cycle ] = TimeToIndex( RemTime ) ;
         RemTo[ Cycle ] = TimeToIndex( RemTime + RemDuration[ Cycle ] ) ;
         SetArray( RR, RemFrom[ Cycle ], RemTo[ Cycle ], RemRR ) ;
         }*/
  /*    } */



   /*
	%2 to 5 sudden falls in mean RR of 10 to 20 % between start of recording and sleep
	%with durations of 3 to 30 minutes to simulate various physical activities!
	if (1)
	numFalls = fix(rand*3)+2;
	for (fall = 1:numFalls)
		duration = (rand*27+3)*60;
		fallTime = rand*sleepTime;
   	fallRR = RRday*(1-(rand*0.1+0.1));
	   fallInd1(fall).vect = find(t > fallTime & t < fallTime+duration);
   	RR(fallInd1(fall).vect) = ones(size(fallInd1(fall).vect))*fallRR;
	end
	end
   */
/*	if( 1 )
   	{ */
      // fprintf( stderr, "  Falls 1...\n" ) ;
		NumFalls1 = (int) floor( Rand0To1( ) * 3.0 ) + 2.0 ;
		for( Fall = 1 ; Fall < NumFalls1 ; Fall ++ )
      	{
			Duration = ( Rand0To1( ) * 27.0 + 3.0 ) * 60.0 ;
			FallTime = Rand0To1( ) * SleepTime ;
   		FallRR = RRday * ( 1.0 - ( Rand0To1( ) * 0.1 + 0.1 ) ) ;
         FallFrom1[ Fall ] =  TimeToIndex( FallTime ) ;
         FallTo1[ Fall ] =  TimeToIndex( FallTime + Duration ) ;
         SetArray( RR, FallFrom1[ Fall ], FallTo1[ Fall ], FallRR ) ;
         }
 /*     }
   */

	/*
	%1 to 2 sudden falls in mean RR of 10 to 20 % between awaking and end of recording
	%with durations of 3 to 30 minutes to simulate various physical activities!
	if (1)
	numFalls = fix(rand)+1;
	for (fall = 1:numFalls)
		duration = (rand*27+3)*60;
		fallTime = rand*(24*60*60-wakeTime)+wakeTime;
   	fallRR = RRday*(1-(rand*0.1+0.1));
	   fallInd2(fall).vect = find(t > fallTime & t < fallTime+duration);
   	RR(fallInd2(fall).vect) = ones(size(fallInd2(fall).vect))*fallRR;
	end
	end
   */
 /*	if( 1 )
   	{ */
      // fprintf( stderr, "  Falls 2...\n" ) ;
		NumFalls2 = (int) floor( Rand0To1( ) ) + 1.0 ;
		for( Fall = 1 ; Fall < NumFalls2 ; Fall ++ )
      	{
			Duration = ( Rand0To1( ) * 27.0 + 3.0 ) * 60.0 ;
			FallTime = Rand0To1( ) * ( 24.0 * 60.0 * 60.0 - WakeTime ) + WakeTime ;
   		FallRR = RRday * ( 1.0 - ( Rand0To1( ) * 0.1 + 0.1 ) ) ;
         FallFrom2[ Fall ] = TimeToIndex( FallTime ) ;
         FallTo2[ Fall ] = TimeToIndex( FallTime + Duration ) ;
         SetArray( RR, FallFrom2[ Fall ], FallTo2[ Fall ], FallRR ) ;
         }
  /*    }
    */

	/*
	%15 to 25 spontaneous arousels of 5 to 20 s duration during sleep
	%where mean sleep RR falls to by 10 to 20 %
	if (1)
	numArouse = fix(rand*10)+15;
	for (arouse = 1:numArouse)
		duration = rand*15+5;
		arouseTime = rand*(wakeTime-sleepTime)+sleepTime;
   	arouseRR = RRnight*(1-(rand*0.1+0.1));
	   arouseInd(fall).vect = find(t > arouseTime & t < arouseTime+duration);
   	RR(arouseInd(fall).vect) = ones(size(arouseInd(fall).vect))*arouseRR;
	end
	end
   */
  /*	if( 1 )
   	{  */
      // fprintf( stderr, "  Arousals...\n" ) ;
		NumArouse = floor( Rand0To1( ) * 10.0 ) + 15.0 ;
      for( Arouse = 1 ; Arouse <= NumArouse ; Arouse ++ )
      	{
			Duration = Rand0To1( ) * 15.0 + 5.0 ;
			ArouseTime = Rand0To1( ) * ( WakeTime - SleepTime ) + SleepTime ;
	   	ArouseRR = RRnight * ( 1.0 - ( Rand0To1( ) * 0.1 + 0.1 ) ) ;
         ArouseFrom[ Arouse ] = TimeToIndex( ArouseTime ) ;
         ArouseTo[ Arouse ] = TimeToIndex( ArouseTime + Duration ) ;
         SetArray( RR, ArouseFrom[ Arouse ], ArouseTo[ Arouse ], ArouseRR ) ;
         }
  /*    }
    */

   /*
	if (1)
	%Baseline wander -  low frequency 2 to 5 hrs
	freqWander = 1/((rand*3+2)*60*60);
	sleepIndex = fix(sleepTime/dt);
	amp = 0.1;
	RR(1:sleepIndex) = RR(1:sleepIndex)+amp*sin(2*pi*freqWander*t(1:sleepIndex));
	end
   */
 /*	if( 1 )
   	{ */
      // fprintf( stderr, "  Baseline wander...\n" ) ;
//		FreqWander = 1.0 / ( ( Rand0To1( ) * 3.0 + 2.0 ) * 24 * 60.0 * 60.0 ) ;
		FreqWander = 1.0/(13.0*60.0*60.0);
      SleepIndex = floor( SleepTime / cDT ) ;
		amp = 0.1 ;
		for( Index = 1 ; Index <= SleepIndex ; Index ++ )
   		RR[ Index ] += amp * (sin( 2.0 * cPI * FreqWander * IndexToTime( Index )+(cPI/6.5) ) +sin( 2.0 * cPI * (FreqWander-(1.0/(5.0*60.0*60.0))) * IndexToTime( Index ) ) +sin( 2.0 * cPI * (FreqWander-(1.0/(8.0*60.0*60.0))) * IndexToTime( Index ) ));
 /*     }   */


	/*
	if (1)
	%measurement artifact
	%0 or 5 artifacts during period of activity with high (1.5 to 2) rr
	for (fall = 1:size(fallInd1))
   	numArtifacts = fix(rand*5);
	   for (artifact = 1:numArtifacts)
   	   randInd = fix(rand*(fallInd1(fall).vect(end)-fallInd1(fall).vect(1))+1);
	   	ind = fallInd1(fall).vect(randInd);
	      RR(ind) = rand*0.5+1.5;
   	end
	end
	end
   */
	for( Fall = 1 ; Fall <= NumFalls1 ; Fall ++ )
   	{
      // fprintf( stderr, "  Artifacts...\n" ) ;
   	NumArtifacts = floor( Rand0To1( ) * 5.0 ) ;
	   for( Artifact = 1 ; Artifact <= NumArtifacts ; Artifact ++ )
      	{
   	   RandInd = floor( Rand0To1( ) * ( FallTo1[ Fall ] - FallFrom1[ Fall ] ) + 1.0 ) ;
	      RR[ (int)RandInd ] = Rand0To1( ) * 0.5 + 1.5 ;
         }
      }




	/* This is where Emma's respiration bit starts */
/*   if( 1)
   	{*/
//      fprintf( stderr, "  Emma's Resp bit...\n" ) ;

      /*initialise variables */
      SetArray( Resp, 1, cSamplesPerDay - 1, 0 ) ;
		i=1 ;
		FallIndex1 = 1 ;
		FallIndex2 = 1 ;
/*      RemIndex = 1 ;*/
		ArouseIndex = 1 ;

		/*ordering the random activities and arousal*/
		Sort( FallFrom1, NumFalls1, OrderedFallFrom1 ) ;
		Sort( FallTo1, NumFalls1, OrderedFallTo1 ) ;
		Sort( ArouseFrom, NumArouse, OrderedArouseFrom ) ;
		Sort( ArouseTo, NumArouse, OrderedArouseTo ) ;
 		Sort( FallFrom2, NumFalls2, OrderedFallFrom2 ) ;
		Sort( FallTo2, NumFalls2, OrderedFallTo2 ) ;




		/* while (1) */
	   while( 1 )
   		{
			/* calculating what is happening at each time instant*/
         IsSleep = ( ( T[ i ] > SleepTime ) && ( T[ i ] < WakeTime ) ) ;
         IsFall1 = ( i <= OrderedFallTo1[ FallIndex1 ] ) && ( i >= OrderedFallFrom1[ FallIndex1 ] ) ;
         IsFall2 = ( i <= OrderedFallTo2[ FallIndex2 ] ) && ( i >= OrderedFallFrom2[ FallIndex2 ] ) ;
			IsArousal = ( i <= OrderedArouseTo[ ArouseIndex ] ) && ( i >= OrderedArouseFrom[ ArouseIndex ] ) ;
/*			IsREM = ( i <= RemTo[ RemIndex ] ) && ( i >= RemFrom[ RemIndex ] ) ;*/

		   if( !IsSleep && IsFall1 )
         	{
            /*this is day (before sleep) and active
		      %so breathing rate increases so RR changes are modulated by a faster less varying frequency
				%the positive part of the cycle
			   %calculating the amp, freq and number of points the half cycle passes through*/
//				fprintf( stderr, "  Emma's Resp bit...\n" ) ;
			   choosefreqamp( cDayActiveAmpMaxPos, cDayActiveAmpMinPos, cDayActiveAmpVarPos, previousamppos, cDayActiveFrqMaxPos, cDayActiveFrqMinPos, cDayActiveFrqVarPos, cSamplesPerDay, i, &numberoftimeelement, &freqpos, &amppos ) ;

			   /*appending samples to the end of the resp signal*/
            for( Index = 0 ; Index < numberoftimeelement ; Index ++ )
            	Resp[ Index + i ] = amppos * sin( 2.0 * cPI * freqpos * T[ Index + 1 ] ) ;

			   i += numberoftimeelement ;

			   /*fall out of the loop if the array is full */
			   if( i >= cSamplesPerDay )
			      break ;

	         /*remembering previous values so changes can be related to previous values */
 			   previousamppos = amppos ;

	 	     	/*negative part of the cycle
				%calculating the amp, freq and number of points the half cycle passes through */
			   choosefreqamp( cDayActiveAmpMaxNeg, cDayActiveAmpMinNeg, cDayActiveAmpVarNeg, previousampneg, cDayActiveFrqMaxNeg, cDayActiveFrqMinNeg, cDayActiveFrqVarNeg, cSamplesPerDay, i, &numberoftimeelement, &freqneg, &ampneg ) ;
				/*%appending samples to the end of the resp signal*/
         	for( Index = 0 ; Index < numberoftimeelement ; Index ++ )
	         	Resp[ Index + i ] = -ampneg * sin( 2.0 * cPI * freqneg *T[ Index + 1 ] ) ;

				i += numberoftimeelement ;
		   	/*%fall out of the loop if the array is full*/
	   		if( i >= cSamplesPerDay )
		   	   break ;

				/*remembering previous values so changes can be related to previous values */
				previousampneg = ampneg ;

				/*got to the end of the active time look at next one*/
		      if( ( i >= OrderedFallTo1[ FallIndex1 ] ) && ( FallIndex1 != NumFalls1 ) )
		         FallIndex1 ++ ;

		      } /* end if */
		   else if( !IsSleep && IsFall2 )
         	{
            /*this is day (after wake up) and active
		      %so breathing rate increases so RR changes are modulated by a faster less varying frequency
				%the positive part of the cycle
			   %calculating the amp, freq and number of points the half cycle passes through*/
				fprintf( stderr, "  Emma's Resp bit2...\n" ) ;
            choosefreqamp( cDayActiveAmpMaxPos, cDayActiveAmpMinPos, cDayActiveAmpVarPos, previousamppos, cDayActiveFrqMaxPos, cDayActiveFrqMinPos, cDayActiveFrqVarPos, cSamplesPerDay, i, &numberoftimeelement, &freqpos, &amppos ) ;

			   /*appending samples to the end of the resp signal*/
            for( Index = 0 ; Index < numberoftimeelement ; Index ++ )
            	Resp[ Index + i ] = amppos * sin( 2.0 * cPI * freqpos * T[ Index + 1 ] ) ;

			   i += numberoftimeelement ;

			   /*fall out of the loop if the array is full */
			   if( i >= cSamplesPerDay )
			      break ;

	         /*remembering previous values so changes can be related to previous values */
 			   previousamppos = amppos ;

	 	     	/*negative part of the cycle
				%calculating the amp, freq and number of points the half cycle passes through */
			   choosefreqamp( cDayActiveAmpMaxNeg, cDayActiveAmpMinNeg, cDayActiveAmpVarNeg, previousampneg, cDayActiveFrqMaxNeg, cDayActiveFrqMinNeg, cDayActiveFrqVarNeg, cSamplesPerDay, i, &numberoftimeelement, &freqneg, &ampneg ) ;
				/*%appending samples to the end of the resp signal*/
         	for( Index = 0 ; Index < numberoftimeelement ; Index ++ )
	         	Resp[ Index + i ] = -ampneg * sin( 2.0 * cPI * freqneg *T[ Index + 1 ] ) ;

				i += numberoftimeelement ;
		   	/*%fall out of the loop if the array is full*/
	   		if( i >= cSamplesPerDay )
		   	   break ;

				/*remembering previous values so changes can be related to previous values */
				previousampneg = ampneg ;

				/*got to the end of the active time look at next one*/
		      if( ( i >= OrderedFallTo2[ FallIndex2 ] ) && ( FallIndex2 != NumFalls2 ) )
		         FallIndex2 ++ ;

		      } /* end if */

		   else if( !IsSleep )
         	{
				/*%this is day and not active
		      %large variation in breathing frequency so RR interval changes are modulated by a large number of frequencies
				%the positive part of the cycle
			   %calculating the amp, freq and number of points the half cycle passes through*/
			   choosefreqamp( cDayPassiveAmpMaxPos, cDayPassiveAmpMinPos , cDayPassiveAmpVarPos, previousamppos, cDayPassiveFrqMaxPos, cDayPassiveFrqMinPos, cDayPassiveFrqVarPos, cSamplesPerDay, i, &numberoftimeelement, &freqpos, &amppos ) ;

			   /*appending samples to the end of the resp signal*/
            for( Index = 0 ; Index < numberoftimeelement ; Index ++ )
            	Resp[ Index + i ] = amppos * sin( 2.0 * cPI * freqpos * T[ Index + 1 ] ) ;

			   i += numberoftimeelement ;

			   /*%fall out of the loop if the array is full*/
			   if( i >= cSamplesPerDay )
			      break ;

				/*%remembering previous values so changes can be related to previous values*/
            previousamppos = amppos ;

				/*negative part of the cycle
				%calculating the amp, freq and number of points the half cycle passes through*/
  			   choosefreqamp( cDayPassiveAmpMaxNeg, cDayPassiveAmpMinNeg, cDayPassiveAmpVarNeg, previousampneg, cDayPassiveFrqMaxNeg, cDayPassiveFrqMinNeg, cDayPassiveFrqVarNeg, cSamplesPerDay, i, &numberoftimeelement, &freqneg, &ampneg ) ;

				/*appending samples to the end of the resp signal*/
            for( Index = 0 ; Index < numberoftimeelement ; Index ++ )
            	Resp[ Index + i ] = -ampneg * sin( 2.0 * cPI * freqneg * T[ Index + 1 ] ) ;

			   i += numberoftimeelement ;

				/*%fall out of the loop if the array is full*/
			   if( i >= cSamplesPerDay )
			      break ;

				/*%remembering previous values so changes can be related to previous values*/
		   	previousampneg = ampneg ;

            } /* end else if */

		   else if( IsSleep && IsArousal )
         	{
            /*%this is night and aroused
		      %when aroused breathing frequency increases so RR changes are modulated by a faster frequency
				%the positive part of the cycle
            %calculating the amp, freq and number of points the half cycle passes through*/
			   choosefreqamp( cNightActiveAmpMaxPos, cNightActiveAmpMinPos, cNightActiveAmpVarPos, previousamppos, cNightActiveFrqMaxPos, cNightActiveFrqMinPos, cNightActiveFrqVarPos, cSamplesPerDay, i, &numberoftimeelement, &freqpos, &amppos ) ;

				/*%appending samples to the end of the resp signal*/
				for( Index = 0 ; Index < numberoftimeelement ; Index ++ )
            	Resp[ Index + i ] = amppos * sin( 2.0 * cPI * freqpos * T[ Index + 1 ] ) ;

			   i += numberoftimeelement ;

	   		/*%fall out of the loop if the array is full*/
            if( i >= cSamplesPerDay )
			      break ;

				/*%remembering previous values so changes can be related to previous values*/
            previousamppos = amppos ;


/*			   %appending samples to the end of the resp signal
				resp(i:i+numberoftimeelement-1)=-ampneg*sin(2*pi*freqneg*t(1:numberoftimeelement));
			   %incrementing i
				i=i+numberoftimeelement;
			   %fall out of the loop if the array is full
			   if (i >= length(t))
			      break;
			   end
				%remembering previous values so changes can be related to previous values
	   		previousampneg=ampneg;
            */

            /*negative part of the cycle
			   %calculating the amp, freq and number of points the half cycle passes through*/
			   choosefreqamp( cNightActiveAmpMaxNeg, cNightActiveAmpMinNeg, cNightActiveAmpVarNeg, previousampneg, cNightActiveFrqMaxNeg, cNightActiveFrqMinNeg, cNightActiveFrqVarNeg, cSamplesPerDay, i, &numberoftimeelement, &freqneg, &ampneg ) ;
				for( Index = 0 ; Index < numberoftimeelement ; Index ++ )
            	Resp[ Index + i ] = -ampneg * sin( 2.0 * cPI * freqneg * T[ Index + 1 ] ) ;
			   i += numberoftimeelement ;
			   if( i >= cSamplesPerDay )
			      break ;

		   	previousampneg = ampneg ;

            /*
		      %got to the end of the aroused time look at next one
		      if (t(i)>=t(orderarouseend(ArouseIndex)) & ArouseIndex~=numArouse)
      		   ArouseIndex=ArouseIndex+1;
		      end
            */
		      if( ( i >= OrderedArouseTo[ ArouseIndex ] ) && ( ArouseIndex != NumArouse ) )
      		   ArouseIndex ++ ;

				} /* end else if */


         /*
		   elseif ((t(i)>sleepTime) & (t(i)<wakeTime)) & (t(i)<=t(orderremend(RemIndex)) & t(i)>=t(orderremstart(RemIndex)))
		      %this is night and REM
      		%when aroused breathing frequency increases so RR changes are modulated by a faster frequency

				%the positive part of the cycle
	   		%calculating the amp, freq and number of points the half cycle passes through
			   [numberoftimeelement,freqpos,amppos]=choosefreqamp(0.2,0.1,0.15,previousamppos,0.2666,0.0666,1,dt,length(t),i);
	   		%appending samples to the end of the resp signal
			   resp(i:i+numberoftimeelement-1)= amppos*sin(2*pi*freqpos*t(1:numberoftimeelement));
	   		%incrementing i
			   i=i+numberoftimeelement;
	   		%fall out of the loop if the array is full
			   if (i >= length(t))
	   		   break;
			   end
				%remembering previous values so changes can be related to previous values
			   previousamppos=amppos;
         */
/*		   else if( IsSleep && IsREM )
				{
			   choosefreqamp( cNightActiveAmpMaxPos, cNightActiveAmpMinPos, cNightActiveAmpVarPos, previousamppos, cNightActiveFrqMaxPos, cNightActiveFrqMinPos, cNightActiveFrqVarPos, cSamplesPerDay, i, &numberoftimeelement, &freqpos, &amppos ) ;
				for( Index = 0 ; Index < numberoftimeelement ; Index ++ )
            	Resp[ Index + i ] = amppos * sin( 2.0 * cPI * freqpos * T[ Index + 1 ] ) ;
			   i += numberoftimeelement ;
			   if( i >= cSamplesPerDay )
			      break ;

			   previousamppos = amppos ;*/

				/*
		  		%negative part of the cycle
	   		%calculating the amp, freq and number of points the half cycle passes through
			   [numberoftimeelement,freqneg,ampneg]=choosefreqamp(0.2,0.1,0.15,previousampneg,0.2666,0.0666,1,dt,length(t),i);
	   		%appending samples to the end of the resp signal
				resp(i:i+numberoftimeelement-1)=-ampneg*sin(2*pi*freqneg*t(1:numberoftimeelement));
	   		%incrementing i
				i=i+numberoftimeelement;
	   		%fall out of the loop if the array is full
			   if (i >= length(t))
	   		   break;
			   end
				%remembering previous values so changes can be related to previous values
			   previousampneg=ampneg;
            */
/*			   choosefreqamp( cNightActiveAmpMaxNeg, cNightActiveAmpMinNeg, cNightActiveAmpVarNeg, previousamppos, cNightActiveFrqMaxNeg, cNightActiveFrqMinNeg, cNightActiveFrqVarNeg, cSamplesPerDay, i, &numberoftimeelement, &freqpos, &amppos ) ;
				for( Index = 0 ; Index < numberoftimeelement ; Index ++ )
            	Resp[ Index + i ] = -ampneg * sin( 2.0 * cPI * freqneg * T[ Index + 1 ] ) ;
			   i += numberoftimeelement ;
			   if( i >= cSamplesPerDay )
			      break ;

		   	previousampneg = ampneg ;*/

            /*
		      %got to the end of the aroused time look at next one
		      if (t(i)>=t(orderremend(RemIndex)) & RemIndex~=numRem)
      		   RemIndex=RemIndex+1;
		      end
            */
/*		      if( ( i >= RemTo[ RemIndex ] ) && ( RemIndex != NumRemCycles ) )
      		   RemIndex ++ ;

            }*/ /* end else if */

         /*
		   elseif ((t(i)>sleepTime) & (t(i)<wakeTime))
      		%this is night and not rem or aroused
		      %less variation in frequency lower breathing rates

				%the positive part of the cycle
			   %calculating the amp, freq and number of points the half cycle passes through
	   		[numberoftimeelement,freqpos,amppos]=choosefreqamp(0.3,0.1,0.15,previousamppos,0.2,0.0666,0.15,dt,length(t),i);
			   %appending samples to the end of the resp signal
	   		resp(i:i+numberoftimeelement-1)= amppos*sin(2*pi*freqpos*t(1:numberoftimeelement));
			   %incrementing i
	   		i=i+numberoftimeelement;
			   %fall out of the loop if the array is full
	   		if (i >= length(t))
			      break;
	   		end
				%remembering previous values so changes can be related to previous values
	   		previousamppos=amppos;
         */
		   else if( IsSleep )
         	{
			   choosefreqamp( cNightPassiveAmpMaxPos, cNightPassiveAmpMinPos, cNightPassiveAmpVarPos, previousamppos, cNightPassiveFrqMaxPos, cNightPassiveFrqMinPos, cNightPassiveFrqVarPos, cSamplesPerDay, i, &numberoftimeelement, &freqpos, &amppos ) ;
				for( Index = 0 ; Index < numberoftimeelement ; Index ++ )
            	Resp[ Index + i ] = amppos * sin( 2.0 * cPI * freqpos * T[ Index + 1 ] ) ;
			   i += numberoftimeelement ;
			   if( i >= cSamplesPerDay )
			      break ;

			   previousamppos = amppos ;

            /*
		  		%negative part of the cycle
	   		%calculating the amp, freq and number of points the half cycle passes through
			   [numberoftimeelement,freqneg,ampneg]=choosefreqamp(0.3,0.05,0.15,previousampneg,0.2,0.0666,0.15,dt,length(t),i);
	   		%appending samples to the end of the resp signal
				resp(i:i+numberoftimeelement-1)=-ampneg*sin(2*pi*freqneg*t(1:numberoftimeelement));
	   		%incrementing i
				i=i+numberoftimeelement;
	   		%fall out of the loop if the array is full
			   if (i >= length(t))
	   		   break;
			   end
				%remembering previous values so changes can be related to previous values
			   previousampneg=ampneg;
            */
			   choosefreqamp( cNightPassiveAmpMaxNeg, cNightPassiveAmpMinNeg, cNightPassiveAmpVarNeg, previousampneg, cNightPassiveFrqMaxNeg, cNightPassiveFrqMinNeg, cNightPassiveFrqVarNeg, cSamplesPerDay, i, &numberoftimeelement, &freqneg, &ampneg ) ;

				for( Index = 0 ; Index < numberoftimeelement ; Index ++ )
            	Resp[ Index + i ] = -ampneg * sin( 2.0 * cPI * freqneg * T[ Index + 1 ] ) ;
			   i += numberoftimeelement ;
			   if( i >= cSamplesPerDay )
			      break ;

		   	previousampneg = ampneg ;

            } /* end else if */

         } /* end while( 1 ) */

     	for( Index = 0 ; Index < cSamplesPerDay ; Index ++ )
			RR[ Index ] += Resp[ Index ] ;

	  /*	}  end if( 1 ) */
	}


/* function rrsim(dummy) */
/* Calls RRSim1Day as many times as necessary */
void RRSim( void )
	{
   int Day ;
   int StartIndex ;
   int Index ;
   double RRday ;
   double RRnight ;
   double amp ;
   double phase ;
   double deltaphase ;
   double freq1 ;
/*   double ThisFreq ;*/

	/*
	%Define mean night and day RR
	%day RR between 0.6 & 0.9 s
	%night RR between 0.8 & 1.2 s
	%difference between night & day RR greater than 0.2 s
	while (1)
		RRday = rand*0.3+0.6;
   	RRnight = rand*0.4+0.8;
	   if ((RRnight - RRday) > 0.2)
   	   break;
	   end
	end
	*/

   // fprintf( stderr, "\nSimulating background...\n" ) ;
	while( 1 )
		{
	   RRday = Rand0To1( ) * 0.3 + 0.8 ;
   	RRnight = Rand0To1( ) * 0.4 + 0.8 ;
	   if ( ( RRnight - RRday ) > 0.2 )
   	   break ;
	   }


	/*
	if (1)
	%0.1 Hz effect
	%amplitude 5 to 10 % of mean daytime RR
	%frequency 0.05 to 0.2 Hz - multiple components to smooth peak out
	freq1 = (rand*0.15)+0.05;
	amp = RRday*(rand*0.008+0.008);
	freq=[freq1-0.005:0.001:freq1+0.005];
	lfeffect = zeros(size(t));
	for (freqNo = 1:length(freq))
   	lfeffect = lfeffect+amp*sin(2*pi*freq(freqNo)*t);
	end
	RR = RR+lfeffect;
	subplot(4,1,2),plot([0:dt/(60*60):24],RR);
	set(gca,'xtick',[0:6:24]);
	xlim([0 24]);
	end
   */
 /*	if( 0 )
   	{
      fprintf( stderr, "  Phil's 0.1 Hz effect...\n" ) ;
		freq1 = ( Rand0To1( ) * 0.15 ) + 0.05 ;
      for( Index = 0 ; Index < cSamplesTotal ; Index ++ )
			nLFeffect[ Index ] = 0.0 ;
      for( ThisFreq = freq1 - 0.005 ; ThisFreq <= freq1 + 0.005 ; ThisFreq +=0.001 )
      	{
			amp = RRday * ( Rand0To1( ) * 0.008 + 0.008 ) ;
			phase = Rand0To1( ) * 2.0 * cPI ;
         for( Index = 0 ; Index < cSamplesPerDay ; Index ++ )
         	nLFeffect[ Index ] = nLFeffect[ Index ] + amp * sin( 2.0 * cPI * ThisFreq * IndexToTime( Index ) + phase );
         }

      for( Index = 0 ; Index < cSamplesTotal ; Index ++ )
			nRR[ Index ] += nLFeffect[ Index ] ;
		}

   */
	/* My effort at a 0.1 Hz effect */
/*	if( 1 )
   	{ */
      // fprintf( stderr, "  MD's 0.1 Hz effect...\n" ) ;
      amp = RRday * ( Rand0To1( ) * 0.008 + 0.008 ) ;
		freq1 = ( Rand0To1( ) * 0.15 ) + 0.05 ;
      for( Index = 0 ; Index < cSamplesTotal ; Index ++ )
			nLFeffect[ Index ] = 0.0 ;
      phase = 0.0 ;
      for( Index = 0 ; Index < cSamplesTotal ; Index ++ )
      	{
         if( Rand0To1( ) < 0.01 )
         	freq1 += ( Rand0To1( ) - 0.5 ) * 0.01 ;
         if( freq1 == 0.2 )
         	freq1 = 0.17 ;
         if( freq1 == 0.05 )
         	freq1 = 0.08 ;

	      deltaphase = 2.0 * cPI * freq1 * cDT ;
			amp += Rand0To1( ) * 0.002 - 0.001 ;
         amp = Min( amp, 0.016 ) ;
         amp = Max( amp, 0.008 ) ;
			phase += deltaphase * ( 0.95 + Rand0To1( ) * 0.1 ) ;
         nLFeffect[ Index ] = amp * sin( phase );
         }

      for( Index = 0 ; Index < cSamplesTotal ; Index ++ )
			nRR[ Index ] += nLFeffect[ Index ] ;
/*		}   */



	/*
	if (1)
	%1/f effect (pink noise)
	%amplitude 1 to 2 % of mean daytime RR
	whitenoise=randn(size(t));
	%poles and zero from Mr. Tom Bruhns of HP
	%for white to pink noise filter for audio signals
	%poles=[0.9986823 0.9914651 0.9580812 0.8090598 0.2896591]';
	%zeros=[0.9963594 0.9808756 0.9097290 0.6128445 -0.0324723]';
	%get the numerator and denominator of the filter
	%[b,a]=zp2tf(zeros,poles,1);
	a = [  1.0000   -4.0469    6.3705   -4.8224    1.7212   -0.2223];
	b = [  1.0000   -3.4673    4.4317   -2.4428    0.4608    0.0177];
	%filter to get pink noise
	pinknoise=filter(b,a,whitenoise);
	amp = RRday*(rand*0.025+0.025);
	RR = RR+amp*pinknoise;
	end
   */
  /*	if( 1 )
   	{  */
      // fprintf( stderr, "  White noise...\n" ) ;
     	for( Index = 0 ; Index < cSamplesTotal ; Index ++ )
      	nWhitenoise[ Index ] = RandGauss( ) ;

		// Filter( whitenoise, cSamplesPerDay, CoeffA, CoeffB, 5, Pinknoise ) ;
      // fprintf( stderr, "  Pink noise...\n" ) ;
		FirstOrderLowPassFilter( nWhitenoise, cSamplesTotal, 0.1, nPinknoise ) ;
		amp = RRday * ( Rand0To1( ) * 0.025 + 0.025 ) ;
     	for( Index = 0 ; Index < cSamplesPerDay ; Index ++ )
			nRR[ Index ] += amp * nPinknoise[ Index ] ;
 /*  	}
   */

   for( Day = 0 ; Day < cMaxDays ; Day ++ )
   	{
      // fprintf( stderr, "\nSimulating day %d...\n", Day + 1 ) ;
      StartIndex = Day * cSamplesPerDay ;
	   RRSim1Day( RRday, RRnight, nT + StartIndex, nRR + StartIndex, nResp + StartIndex ) ;
      }
   }



/* This function is called once, before generate() is called.  See main()
   for a description of its arguments.
*/

void initialize(long seed, long tmax)
	{
   int Index=0 ;
  /* FILE *OutputFile ; */

   /* tmax not used here, added to aviod compiler warning */
   Index = tmax < Index ? tmax : Index;
   /* Allocate memory for the arrays */
   /* Abandon hope if you can't */
   nT = malloc( cMaxAlloc * sizeof( double ) ) ;
   if( nT != NULL )
		nRR = malloc( cMaxAlloc * sizeof( double ) ) ;
   if( nRR != NULL )
		nLFeffect = malloc( cMaxAlloc * sizeof( double ) ) ;
   if( nLFeffect != NULL )
	   nWhitenoise = malloc( cMaxAlloc * sizeof( double ) ) ;
   if( nWhitenoise != NULL )
		nPinknoise = malloc( cMaxAlloc * sizeof( double ) ) ;
   if( nPinknoise != NULL )
		nResp = malloc( cMaxAlloc * sizeof( double ) ) ;

   if( nResp == NULL )
   	{
      fprintf( stderr, "Couldn't allocate enough memory to continue." ) ;
      return ;
      }

   srand((unsigned int) seed ) ;
   RRSim( ) ;

   /* Write the intermediate data to the file RR.num in the working directory */
/*   fprintf( stderr, "\nWriting to RR.num...\n" ) ;
   OutputFile = fopen( "RR.num", "wt" ) ;
   if( OutputFile == NULL )
      return ;
   for( Index = 1 ; Index < cSamplesTotal ; Index ++ )
   	fprintf( OutputFile, "%8.3lf %8.3lf %8.3lf %8.3lf %8.3lf %8.3lf\n",
      	nT[ Index ], nLFeffect[ Index ], nWhitenoise[ Index ], nPinknoise[ Index ], nResp[ Index ], nRR[ Index ] ) ;

   fclose( OutputFile ) ;*/
	}


/* This function is called once per RR interval.  It should return the length
   of the next (simulated) RR interval in seconds.

   The example code generates samples of a noisy sine wave.
*/

/*
float generate(void)
{
    float rr;
    static float t;

    rr = 0.8 + 0.05*sin(t/5.5) + ((float)rand() - RAND_MAX)/(RAND_MAX*100.0);
    t += rr;
    return (rr);
}
*/


/* Generate a series of RR intervals */
/* Start with the first element in the RR array */
/* This is the first RR interval - 'ThisRR' */
/* Following this heartbeat, time has moved on by 'ThisRR' seconds */
/* Which tells us where to find the next RR interval from */
float generate( void )
	{
   static double CurrentTime = 0.0 ;
   double ThisRR ;
   int Index ;
   Index = TimeToIndex( CurrentTime ) + 1 ;
   if( Index < cSamplesTotal )
   	ThisRR = nRR[ Index ] ;
	CurrentTime += ThisRR ;
   return ThisRR ;
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
/*#if 0*/
    while ((t += generate()) < tmax) {	/* add RR interval to running time */
	/* calculate and output a quantized RR interval */
	ts = (long)(SPS*t + 0.5);
	printf("%5.3f\n", (ts - tsp)/((float)SPS));
	tsp = ts;
    }
/*#endif   */

   free( nT ) ;
	free( nRR ) ;
	free( nLFeffect ) ;
   free( nWhitenoise ) ;
	free( nPinknoise ) ;
	free( nResp ) ;

    exit(0);
}

