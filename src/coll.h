/* header file for coll.c --
 * collision handling routines
 * Copyright (C) 1995-2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: coll.h,v 1.1 2006/10/22 04:05:33 kichiki Exp $
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
#ifndef	_COLL_H_
#define	_COLL_H_


/*
 * INPUT
 *  sys : system parameters
 *  x [np * 3] : position of particles
 *  v [nm * 3] : velocity of particles before collisions
 *  en : elastic constant
 * OUTPUT
 *  v [nm * 3] : velocity of particles after collisions
 */
void
collide_particles (struct stokes * sys,
		   double * x, double * v, double en);

/*
 * INPUT
 *  sys : system parameters
 *  x [np * 3] : position of particles
 *  v [nm * 3] : velocity of particles before collisions
 *  en : elastic constant
 *  x_wall : position of the wall
 *  v_wall : 
 * OUTPUT
 *  v [nm * 3] : velocity of particles after collisions
 */
void
collide_wall_x (struct stokes * sys,
		double * x, double * v, double en,
		double x_wall, double v_wall);
void
collide_wall_y (struct stokes * sys,
		double * x, double * v, double en,
		double y_wall, double v_wall);
void
collide_wall_z (struct stokes * sys,
		double * x, double * v, double en,
		double z_wall, double v_wall);

/*
 * INPUT
 *  sys : system parameters
 *  x [np * 2] : position of particles
 *  v [nm * 2] : velocity of particles before collisions
 *  en : elastic constant
 * OUTPUT
 *  v [nm * 2] : velocity of particles after collisions
 */
void
collide_particles_2d (struct stokes * sys,
		      double * x, double * v, double en);


#endif /* !_COLL_H_ */
