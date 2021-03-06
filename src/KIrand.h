/* header file for Kirand.c --
 * KIrand -- wrapper of random number generator MT19937
 * $Id: KIrand.h,v 1.2 2008/06/01 17:04:22 kichiki Exp $
 */
#ifndef	_KIRAND_H_
#define	_KIRAND_H_

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

/* Period parameters */  
#define MTRNG_N 624
#define MTRNG_M 397
#define MTRNG_MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define MTRNG_UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define MTRNG_LOWER_MASK 0x7fffffffUL /* least significant r bits */

struct KIrand {
  unsigned long mt[MTRNG_N]; /* the array for the state vector  */
  int mti; /* mti==N+1 means mt[N] is not initialized */

  /* for KIrand_Gaussian() */
  int Gaussian_has_saved;
  double Gaussian_saved;
};

struct KIrand *
KIrand_init (void);

void
KIrand_free (struct KIrand *r);


/* load the state for the random number generator
 */
void
KIrand_load_mt (struct KIrand *r,
		unsigned long *mt, int mti);

/* restore the state for the Gaussian random number generator
 */
void
KIrand_load_Gaussian (struct KIrand *rng,
		      unsigned long *mt, int mti,
		      int flag_saved,
		      double x_saved);


/* double random numbers with gaussian distribution
 * taken from gromacs-3.2.1/src/gmxlib/gmx_random.c in GROMACS;
 * gmx_rng_gaussian_real()
 */
double
KIrand_Gaussian (struct KIrand * rng);


/* initializes mt[N] with a seed */
void KIrand_init_genrand (struct KIrand *r,
			  unsigned long s);

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void KIrand_init_by_array (struct KIrand *r,
			   unsigned long init_key[], int key_length);

/* generates a random number on [0,0xffffffff]-interval */
unsigned long
KIrand_genrand_int32 (struct KIrand *r);

/* generates a random number on [0,0x7fffffff]-interval */
long
KIrand_genrand_int31 (struct KIrand *r);

/* generates a random number on [0,1]-real-interval */
double
KIrand_genrand_real1(struct KIrand *r);

/* generates a random number on [0,1)-real-interval */
double
KIrand_genrand_real2 (struct KIrand *r);

/* generates a random number on (0,1)-real-interval */
double
KIrand_genrand_real3 (struct KIrand *r);

/* generates a random number on [0,1) with 53-bit resolution*/
double
KIrand_genrand_res53 (struct KIrand *r);

/* These real versions are due to Isaku Wada, 2002/01/09 added */


#endif /* !_KIRAND_H_ */
