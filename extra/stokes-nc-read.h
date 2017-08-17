/* header file for stokes-nc-read.c --
 * NetCDF interface for libstokes
 * Copyright (C) 2006-2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes-nc-read.h,v 5.6 2008/06/03 02:34:14 kichiki Exp $
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 */
#ifndef	_STOKES_READ_NC_H_
#define	_STOKES_READ_NC_H_


/* open stokes_nc file in NC_NOWRITE mode
 * this is for usual analysis
 */
struct stokes_nc *
stokes_nc_open (const char * filename);

/* open stokes_nc file in NC_WRITE mode
 * this is for continuation (appending the results)
 */
struct stokes_nc *
stokes_nc_reopen (const char * filename);

/* read 1d array [vec/stt/np/npf]
 * INPUT
 *  name : either one of them, Ui0, Oi0, Ei0, Ui, Oi, Ei, a, af, l
 * OUTPUT
 *  x[]
 */
void
stokes_nc_get_array1d (const struct stokes_nc *nc,
		       const char * name,
		       double * x);
/* read constant data for particles in 2d array [np/npf][vec/stt]
 */
void
stokes_nc_get_data0 (const struct stokes_nc *nc,
		     const char * name,
		     double * x);
/* read time-dep. particle data at step in 3d array [step][np/npf][vec/stt]
 */
void
stokes_nc_get_data (const struct stokes_nc *nc,
		    const char * name,
		    int step,
		    double * x);

/* read (the whole) time vector
 * INPUT
 *  time[nc->ntime]
 * OUTPUT
 *  time[nc->ntime]
 */
void
stokes_nc_get_time (const struct stokes_nc *nc,
		    double * time);

/* read time at a step
 * INPUT
 *  step
 * OUTPUT
 *  returned value : time[step]
 */
double
stokes_nc_get_time_step (const struct stokes_nc *nc,
			 int step);

/* read rng data at time (step)
 * INPUT
 *  step
 */
void
stokes_nc_get_rng (struct stokes_nc *nc,
		   int step,
		   struct KIrand *rng);


#endif /* !_STOKES_READ_NC_H_ */
