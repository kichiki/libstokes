/* header file for lub.c --
 * lubrication routines -- atimes procedure
 * Copyright (C) 1993-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: lub.h,v 5.1 2006/10/12 15:02:12 ichiki Exp $
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
#ifndef	_LUB_H_
#define	_LUB_H_


/* calculate lubrication f by uoe for all particles
 * under the periodic boundary condition
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 *   u [np * 3] : velocity, angular velocity, strain
 * OUTPUT
 *   f [np * 3] : force, torque, stresslet
 */
void
calc_lub_ewald_3f (struct stokes * sys,
		   const double * u,
		   double * f);
/* calculate lubrication ft by uoe for all particles
 * under the periodic boundary condition
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 *   uo [np * 6] : velocity, angular velocity, strain
 * OUTPUT
 *   ft [np * 6] : force, torque, stresslet
 */
void
calc_lub_ewald_3ft (struct stokes * sys,
		    const double * uo, double * ft);
/* calculate lubrication fts by uoe for all particles
 * under the periodic boundary condition
 * INPUT
 *   sys : system parameters. following entries are used;
 *         sys->pos
 *         sys->ll[xyz]
 *   uoe [np * 11] : velocity, angular velocity, strain
 * OUTPUT
 *   fts [np * 11] : force, torque, stresslet
 */
void
calc_lub_ewald_3fts (struct stokes * sys,
		     const double * uoe, double * fts);


#endif /* !_LUB_H_ */
