/* KIrand -- wrapper of random number generator MT19937
 * $Id: KIrand.c,v 1.1 2007/09/29 20:41:49 kichiki Exp $
 */

/* 
   A C-program for MT19937, with initialization improved 2002/1/26.
   Coded by Takuji Nishimura and Makoto Matsumoto.

   Before using, initialize the state by using init_genrand(seed)  
   or init_by_array(init_key, key_length).

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)
*/

#include <stdio.h>
#include <stdlib.h> /* malloc() */
#include <math.h> /* sqrt(), log() */

#include "KIrand.h" /* KIrand */


struct KIrand *
KIrand_init (void)
{
  struct KIrand *r = NULL;

  r = (struct KIrand *) malloc (sizeof (struct KIrand));
  r->mti = MTRNG_N + 1;

  r->Gaussian_has_saved = 0; /* not saved */
  r->Gaussian_saved = 0.0;

  return (r);
}

void
KIrand_free (struct KIrand *r)
{
  free (r);
}


/* double random numbers with gaussian distribution
 * taken from gromacs-3.2.1/src/gmxlib/gmx_random.c in GROMACS;
 * gmx_rng_gaussian_real()
 */
double
KIrand_Gaussian (struct KIrand * rng)
{
  double x, y, r;

  if (rng->Gaussian_has_saved)
    {
      rng->Gaussian_has_saved = 0;
      return (rng->Gaussian_saved);
    }
  else
    {
      do
	{
	  /* use uniform random number on [0, 1) */
	  x = 2.0 * KIrand_genrand_real2 (rng) - 1.0;
	  y = 2.0 * KIrand_genrand_real2 (rng) - 1.0;
	  r = x * x + y * y;
	}
      while (r > 1.0 || r == 0.0);

      r = sqrt (- 2.0 * log (r) / r);
      rng->Gaussian_saved = y * r; /* save second random number */
      rng->Gaussian_has_saved = 1;
      return (x * r); /* return first random number */
    }
}


/* initializes mt[N] with a seed */
void KIrand_init_genrand (struct KIrand *r,
			  unsigned long s)
{
  r->mt[0]= s & 0xffffffffUL;
  for (r->mti = 1; r->mti < MTRNG_N; r->mti ++)
    {
      r->mt[r->mti] = 
	(1812433253UL * (r->mt[r->mti - 1] ^ (r->mt[r->mti - 1] >> 30))
	 + r->mti);
      /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
      /* In the previous versions, MSBs of the seed affect   */
      /* only MSBs of the array mt[].                        */
      /* 2002/01/09 modified by Makoto Matsumoto             */
      r->mt[r->mti] &= 0xffffffffUL;
      /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void KIrand_init_by_array (struct KIrand *r,
			   unsigned long init_key[], int key_length)
{
  int i, j, k;

  KIrand_init_genrand (r, 19650218UL);
  i=1;
  j=0;
  k = (MTRNG_N > key_length ? MTRNG_N : key_length);
  for (; k; k--)
    {
      r->mt[i] = (r->mt[i] ^ ((r->mt[i-1] ^ (r->mt[i-1] >> 30)) * 1664525UL))
	+ init_key[j] + j; /* non linear */
      r->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++;
      j++;
      if (i >= MTRNG_N)
	{
	  r->mt[0] = r->mt[MTRNG_N - 1];
	  i=1;
	}
      if (j >= key_length) j=0;
    }
  for (k = MTRNG_N - 1; k; k--)
    {
      r->mt[i] =
	(r->mt[i] ^ ((r->mt[i-1] ^ (r->mt[i-1] >> 30)) * 1566083941UL))
	- i; /* non linear */
      r->mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
      i++;
      if (i >= MTRNG_N)
	{
	  r->mt[0] = r->mt[MTRNG_N - 1];
	  i=1;
	}
    }

  r->mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long
KIrand_genrand_int32 (struct KIrand *r)
{
  unsigned long y;
  static unsigned long mag01[2]={0x0UL, MTRNG_MATRIX_A};
  /* mag01[x] = x * MATRIX_A  for x=0,1 */

  if (r->mti >= MTRNG_N)
    { /* generate N words at one time */
      int kk;

      if (r->mti == MTRNG_N + 1)   /* if init_genrand() has not been called, */
	KIrand_init_genrand(r, 5489UL); /* a default initial seed is used */

      for (kk = 0; kk < MTRNG_N - MTRNG_M; kk ++)
	{
	  y = (r->mt[kk]   & MTRNG_UPPER_MASK)
	    | (r->mt[kk+1] & MTRNG_LOWER_MASK);
	  r->mt[kk] = r->mt[kk + MTRNG_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
      for (; kk < MTRNG_N - 1; kk ++)
	{
	  y = (r->mt[kk]   & MTRNG_UPPER_MASK)
	    | (r->mt[kk+1] & MTRNG_LOWER_MASK);
	  r->mt[kk] = r->mt[kk + (MTRNG_M - MTRNG_N)]
	    ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
      y = (r->mt[MTRNG_N - 1] & MTRNG_UPPER_MASK)
	| (r->mt[0]           & MTRNG_LOWER_MASK);
      r->mt[MTRNG_N-1] = r->mt[MTRNG_M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

      r->mti = 0;
    }
  
  y = r->mt[r->mti++];

  /* Tempering */
  y ^= (y >> 11);
  y ^= (y << 7) & 0x9d2c5680UL;
  y ^= (y << 15) & 0xefc60000UL;
  y ^= (y >> 18);

  return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long
KIrand_genrand_int31 (struct KIrand *r)
{
  return (long)(KIrand_genrand_int32 (r) >> 1);
}

/* generates a random number on [0,1]-real-interval */
double
KIrand_genrand_real1(struct KIrand *r)
{
  return (double) KIrand_genrand_int32 (r) * (1.0 / 4294967295.0); 
  /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double
KIrand_genrand_real2 (struct KIrand *r)
{
  return (double) KIrand_genrand_int32 (r) * (1.0 / 4294967296.0); 
  /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double
KIrand_genrand_real3 (struct KIrand *r)
{
  return (((double)KIrand_genrand_int32 (r)) + 0.5) * (1.0 / 4294967296.0);
  /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double
KIrand_genrand_res53 (struct KIrand *r)
{
  unsigned long a = KIrand_genrand_int32 (r) >> 5,
    b = KIrand_genrand_int32 (r) >> 6;
  return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
}
/* These real versions are due to Isaku Wada, 2002/01/09 added */
