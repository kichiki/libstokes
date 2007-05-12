/* header file for stokes-nc-read.c --
 * NetCDF interface for libstokes
 * Copyright (C) 2006-2007 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes-nc-read.h,v 5.4 2007/05/12 04:29:38 kichiki Exp $
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


struct stokes_nc *
stokes_nc_open (const char * filename);

/* read 1d array [vec/stt]
 * INPUT
 *  name : either one of them, Ui0, Oi0, Ei0, Ui, Oi, Ei
 * OUTPUT
 *  x[]
 */
void
stokes_nc_get_array1d (struct stokes_nc * nc,
		       const char * name,
		       double * x);
/* read constant data for particles in 2d array [np/npf][vec/stt]
 */
void
stokes_nc_get_data0 (struct stokes_nc * nc,
		     const char * name,
		     double * x);
/* read time-dep. particle data at step in 3d array [step][np/npf][vec/stt]
 */
void
stokes_nc_get_data (struct stokes_nc * nc,
		    const char * name,
		    int step,
		    double * x);

/* read lattice vector
 * INPUT
 *  l[nc->nvec]
 * OUTPUT
 *  l[nc->nvec]
 */
void
stokes_nc_get_l (struct stokes_nc * nc,
		 double * l);

/* read (the whole) time vector
 * INPUT
 *  time[nc->ntime]
 * OUTPUT
 *  time[nc->ntime]
 */
void
stokes_nc_get_time (struct stokes_nc * nc,
		    double * time);

/* read time at a step
 * INPUT
 *  step
 * OUTPUT
 *  returned value : time[step]
 */
double
stokes_nc_get_time_step (struct stokes_nc * nc,
			 int step);

#endif /* !_STOKES_READ_NC_H_ */
