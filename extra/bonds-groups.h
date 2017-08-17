/* header file for bonds-groups.c --
 * grouping by bonds among particles
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds-groups.h,v 1.1 2008/07/17 02:16:09 kichiki Exp $
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
#ifndef	_BONDS_GROUPS_H_
#define	_BONDS_GROUPS_H_


#include "bonds.h" // struct BONDS


struct BONDS_GROUP {
  int np;     // number of particles in this group
  int *ip;    // ip[np] : particle index of the group particles
  int *bonds;  // bonds[np-1] : independent bond list
};

struct BONDS_GROUPS {
  int n;                    // number of groups
  struct BONDS_GROUP **group; // grp[ng]
};



struct BONDS_GROUP *
BONDS_GROUP_init (void);

void
BONDS_GROUP_free (struct BONDS_GROUP *g);

void
BONDS_GROUP_set (struct BONDS_GROUP *g,
		 int np,
		 const int *ip,
		 const int *bonds);


struct BONDS_GROUPS *
BONDS_GROUPS_init (void);

void
BONDS_GROUPS_free (struct BONDS_GROUPS *gs);

/* add a group into struct BONDS_GROUPS
 * INPUT
 *  gs         : struct BONDS_GROUPS
 *  np         : number of particles for the group
 *  ip[np]     : particle-index list for the group, with which x of 
 *               the particle i is accessed by pos[ip[i]*3]
 *  bonds[np-1]: bond-index list for the group, with which particles 
 *               of i-th bond is accessed by ia[bonds[i]] and ib[bonds[i]]
 * OUTPUT
 *  gs
 */
void
BONDS_GROUPS_add (struct BONDS_GROUPS *gs,
		  int np, const int *ip, const int *bonds);


/**
 * goupring routines
 */

/* pickup "ig" group out of the table gid[np] and return the table ip[].
 * the particles picked up are erased as gid[i] = -1.
 * INPUT
 *  np : number of particles in the system
 *  gid[np] : group id table. 0  means the single particles
 *                            -1 means erased (unassigned) particles
 * OUTPUT
 *  *n : number of particles in the group
 *  returned value : ip[*n], the particle indices of the group
 */
int *
BONDS_GROUPS_gid_pickup (int np, int *gid, int ig,
			 int *n);

struct BONDS_GROUPS *
BONDS_GROUPS_make (struct BONDS *b, int np);


/**
 * utility routines
 */

/*
 * INPUT
 *  pos[NP*3] : positions of all particles (NP = total number of particles)
 * OUTPUT
 *  conn[np*3] : connector vectors for the bonds "b"
 *               where np is the number of the beads in the group.
 *               the last vector (for (np-1)th element) is the COM.
 */
void
BONDS_pos_to_conn_1 (struct BONDS *b,
		     struct BONDS_GROUP *g,
		     const double *pos,
		     double *conn);

/*
 * INPUT
 *  pos[NP*3] : positions of all particles (NP = total number of particles)
 * OUTPUT
 *  conn[NP*3] : connector vectors for all system
 */
void
BONDS_pos_to_conn (struct BONDS *b,
		   struct BONDS_GROUPS *gs,
		   const double *pos,
		   double *conn);

/*
 * INPUT
 *  conn[np*3] : connector vectors for the bonds "b"
 *               where np is the number of the beads in the group.
 *               the last vector (for (np-1)th element) is the COM.
 * OUTPUT
 *  pos[NP*3] : positions of all particles (NP = total number of particles)
 */
void
BONDS_conn_to_pos_1 (struct BONDS *b,
		     struct BONDS_GROUP *g,
		     const double *conn,
		     double *pos);

/*
 * INPUT
 *  conn[np*3] : connector vectors for the bonds "b"
 *               where np is the number of the beads in the group.
 *               the last vector (for (np-1)th element) is the COM.
 * OUTPUT
 *  pos[NP*3] : positions of all particles (NP = total number of particles)
 */
void
BONDS_conn_to_pos (struct BONDS *b,
		   struct BONDS_GROUPS *gs,
		   const double *conn,
		   double *pos);


#endif /* !_BONDS_GROUPS_H_ */
