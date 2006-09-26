/* header file for fbsys.c --
 * structure for system parameters of stokes library.
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes.h,v 1.1 2006/09/26 01:00:32 ichiki Exp $
 */
#ifndef	_STOKES_H_
#define	_STOKES_H_

struct stokes {
  int np; /* number of all particles */
  int nm; /* number of mobile particles */
  double * pos; /* position of particles */

  /* for ewald codes */
  int pcellx, pcelly, pcellz; /* # of cell in real space */
  int kmaxx, kmaxy, kmaxz; /* # of cell in reciprocal space */

  double zeta,zeta2,zaspi,za2;
  double pi2,pivol;

  double lx, ly, lz;
  double llx[27], lly[27], llz[27]; /* for regist and lub */

  /* for zeta program */
  double cpu1, cpu2, cpu3;
};


/* all elements are zero-cleared
 */
struct stokes *
stokes_init (void);

void
stokes_free (struct stokes * sys);

void
stokes_set_np (struct stokes * sys,
	       int np, int nm);

void
stokes_set_ll (struct stokes * sys,
	       double lx, double ly, double lz);

void
stokes_set_zeta (struct stokes * sys,
		 double zeta, double cutlim);

double
zeta_by_tratio (struct stokes * sys,
		double tratio);

#endif /* !_STOKES_H_ */
