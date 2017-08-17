/* header file for ev-dh-grid.c --
 * excluded-volume interactions by Debye-Huckel with RYUON_grid
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: ev-dh-grid.h,v 1.2 2008/11/01 05:45:40 kichiki Exp $
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
#ifndef	_EV_DH_GRID_H_
#define	_EV_DH_GRID_H_

#include "stokes.h" // struct stokes
#include "ev-dh.h"  // struct EV_DH
#include "grid.h"   // struct RYUON_grid


/*
 * for non-periodic system
 * INPUT
 *  ev_dh      : struct EV_DH
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_DH_calc_force_grid (struct EV_DH *ev_dh,
		       struct stokes *sys,
		       double *f,
		       int flag_add);

/*
 * for periodic system
 * INPUT
 *  ev_dh      : struct EV_DH
 *  sys        : struct stokes (only nm and pos are used)
 *  f [nm * 3] : force is assigned only for the mobile particles
 *  flag_add   : if 0 is given, zero-clear and set the force
 *               otherwise, add the bond force into f[]
 * OUTPUT
 */
void
EV_DH_calc_force_grid_periodic (struct EV_DH *ev_dh,
				struct stokes *sys,
				double *f,
				int flag_add);


#endif /* !_EV_DH_GRID_H_ */
