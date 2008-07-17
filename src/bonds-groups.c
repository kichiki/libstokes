/* grouping by bonds among particles
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: bonds-groups.c,v 1.1 2008/07/17 02:15:34 kichiki Exp $
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
#include <stdlib.h>
#include <math.h> // fabs
#include "memory-check.h" // macro CHECK_MALLOC

#include "bonds.h" // struct BONDS

#include "bonds-groups.h"


struct BONDS_GROUP *
BONDS_GROUP_init (void)
{
  struct BONDS_GROUP *g = (struct BONDS_GROUP *)malloc (sizeof (struct BONDS_GROUP));
  CHECK_MALLOC (g, "BONDS_GROUP_init");

  g->np = 0;
  g->ip = NULL;
  g->bonds = NULL;

  return (g);
}

void
BONDS_GROUP_free (struct BONDS_GROUP *g)
{
  if (g == NULL) return;

  if (g->ip != NULL) free (g->ip);
  if (g->bonds != NULL) free (g->bonds);
  free (g);
}

void
BONDS_GROUP_set (struct BONDS_GROUP *g,
		 int np,
		 const int *ip,
		 const int *bonds)
{
  if (np == 0) return;

  g->np = np;
  g->ip = (int *)malloc (sizeof (int) * np);
  CHECK_MALLOC (g->ip, "BONDS_GROUP_set");
  int i;
  for (i = 0; i < np; i ++)
    {
      g->ip[i] = ip[i];
    }

  if (np == 1)
    {
      g->bonds = NULL;
    }
  else
    {
      g->bonds = (int *)malloc (sizeof (int) * (np-1));
      CHECK_MALLOC (g->bonds, "BONDS_GROUP_set");
    }
  for (i = 0; i < np-1; i ++)
    {
      g->bonds[i] = bonds[i];
    }
}



struct BONDS_GROUPS *
BONDS_GROUPS_init (void)
{
  struct BONDS_GROUPS *gs
    = (struct BONDS_GROUPS *)malloc (sizeof (struct BONDS_GROUPS));
  CHECK_MALLOC (gs, "BONDS_GROUPS_init");

  gs->n = 0;
  gs->group = NULL;

  return (gs);
}

void
BONDS_GROUPS_free (struct BONDS_GROUPS *gs)
{
  if (gs == NULL) return;

  if (gs->group != NULL)
    {
      int i;
      for (i = 0; i < gs->n; i ++)
	{
	  if (gs->group[i] != NULL) BONDS_GROUP_free (gs->group[i]);
	}
      free (gs->group);
    }
  free (gs);
}

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
		  int np, const int *ip, const int *bonds)
{
  gs->n ++;
  gs->group
    = (struct BONDS_GROUP **)realloc (gs->group,
				      sizeof (struct BONDS_GROUP *) * gs->n);
  CHECK_MALLOC (gs->group, "BONDS_GROUPS_add");

  int n = gs->n - 1;
  gs->group[n] = BONDS_GROUP_init ();
  CHECK_MALLOC (gs->group[n], "BONDS_GROUPS_add");

  BONDS_GROUP_set (gs->group[n], np, ip, bonds);
}


/**
 * goupring routines
 */

static void
BONDS_GROUPS_gid_flip (int np, int *gid,
		       int ig0, int ig)
{
  int i;
  for (i = 0; i < np; i ++)
    {
      if (gid[i] == ig0) gid[i] = ig;
    }
}

struct loop_bonds {
  int n;
  int *b; // b[n] : loop bond index
};

struct loop_bonds *
loop_bonds_init (void) {
  struct loop_bonds *lb
    = (struct loop_bonds *)malloc (sizeof (struct loop_bonds));
  CHECK_MALLOC (lb, "loop_bonds_init");

  lb->n = 0;
  lb->b = NULL;

  return (lb);
}

void
loop_bonds_free (struct loop_bonds *lb) {
  if (lb == NULL) return;

  if (lb->b != NULL) free (lb->b);
  free (lb);
}

void
loop_bonds_add (struct loop_bonds *lb, int ib)
{
  lb->n ++;
  lb->b = (int *)realloc (lb->b, sizeof (int) * lb->n);

  int n = lb->n - 1;
  lb->b[n] = ib;
}


/*
 */
static void
BONDS_GROUPS_gid_check (struct BONDS *b,
			int np, int *gid,
			int *nb,
			int **bond,
			int i, int ig,
			struct loop_bonds *lb)
{
  int j;
  for (j = 0; j < nb[i]; j ++)
    {
      int ib = bond[i][j];
      if (b->ib[ib] == i)
	{
	  // we are looking for a one-way map of the connection ia->ib
	  // and this is the reverse map (ib->ia) for "i", so skip
	  continue;
	}
      else if (b->ia[ib] != i)
	{
	  // hum..., it should be b->ia[ib] == i
	  fprintf (stderr, "# something is wrong...\n");
	  exit (1);
	}
      int k = b->ib[ib];
      // now k is the other particle connecting to the particle "i"

      if (gid[k] == -1)
	{
	  /*
	  fprintf (stdout, "# i=%d -> k=%d"
		   " : new member of the group %d\n",
		   i, k, ig);
	  */
	  // particle "k" is a new member for the group "ig"
	  gid[k] = ig;
	  BONDS_GROUPS_gid_check (b,
				  np, gid,
				  nb, bond,
				  k, ig, lb);
	}
      else if (gid[k] == ig)
	{
	  /*
	  fprintf (stdout, "# i=%d -> k=%d"
		   " : already the member of the group %d\n",
		   i, k, ig);
	  */
	  // this must be the loop
	  // that is, "ib" is redundant bond.
	  loop_bonds_add (lb, ib);
	}
      else
	{
	  // particle "k" belong to another group
	  /*
	  fprintf (stdout, "# i=%d -> k=%d"
		   " : another group %d of the group %d\n",
		   i, k, gid[k], ig);
	  */
	  BONDS_GROUPS_gid_flip (np, gid, gid[k], ig);
	}
    }
}

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
			 int *n)
{
  int *ip = NULL;
  (*n) = 0;

  int i;
  for (i = 0; i < np; i ++)
    {
      if (gid[i] == ig)
	{
	  (*n) ++;
	  ip = (int *)realloc (ip, sizeof (int) * (*n));
	  CHECK_MALLOC (ip, "BONDS_GROUPS_gid_pickup");
	  int in = (*n) - 1;
	  ip[in] = i;

	  // clear gid[i] for the later counting process
	  gid[i] = -1;
	}
    }
  return (ip);
}

/*
 * OUTPUT
 *  returned value : 0 == not found
 *                   1 == found
 */
static int
check_ip (int n, const int *ip, int iq)
{
  int i;
  for (i = 0; i < n; i ++)
    {
      if (ip[i] == iq) return (1);
    }
  return (0);
}


struct BONDS_GROUPS *
BONDS_GROUPS_make (struct BONDS *b, int np)
{
  // sort bonds by ia[]
  BONDS_sort_by_ia (b);
  /*
  // check BONDS after sorting
  int i;
  fprintf (stdout, "# b->n = %d =? %d\n", b->n, np);
  for (i = 0; i < b->n; i ++)
    {
      fprintf (stdout, "# b->(ia,ib)[%d] = %d, %d\n",
	       i, b->ia[i], b->ib[i]);
    }
  */

  // make nb[np] : number of bonds for each particle
  int *nb = (int *)calloc (np, sizeof (int));
  CHECK_MALLOC (nb, "BONDS_make_groups");

  int i;
  for (i = 0; i < b->n; i ++)
    {
      if (b->ia[i] < 0 || b->ia[i] >= np)
	{
	  fprintf (stderr, "# BONDS_make_groups()"
		   " : ia=%d is out of range for np=%d\n",
		   b->ia[i], np);
	  exit (1);
	}
      if (b->ib[i] < 0 || b->ib[i] >= np)
	{
	  fprintf (stderr, "# BONDS_make_groups()"
		   " : ib=%d is out of range for np=%d\n",
		   b->ib[i], np);
	  exit (1);
	}

      nb[b->ia[i]] ++;
      nb[b->ib[i]] ++;
    }
  /*
  // check
  for (i = 0; i < np; i ++)
    {
      fprintf (stdout, "# nb[%d] = %d\n", i, nb[i]);
    }
  */

  // make bond[np][nb] : bond list for each particle
  int **bond = (int **)malloc (np * sizeof (int *));
  CHECK_MALLOC (bond, "BONDS_make_groups");
  for (i = 0; i < np; i ++)
    {
      bond[i] = (int *)malloc (nb[i] * sizeof (int));
      CHECK_MALLOC (bond[i], "BONDS_make_groups");
      int j;
      for (j = 0; j < nb[i]; j ++)
	{
	  bond[i][j] = -1;
	}
    }
  for (i = 0; i < b->n; i ++)
    {
      int ia = b->ia[i];
      int j;
      for (j = 0; j < nb[ia]; j ++)
	{
	  if (bond[ia][j] == -1)
	    {
	      bond[ia][j] = i;
	      //fprintf (stdout, "# set (ia) bond[%d][%d] = %d\n", ia, j, i);
	      break;
	    }
	}
      if (bond[ia][j] == -1)
	{
	  fprintf (stderr, "# BONDS_make_groups()"
		   " : something is wrong on bond[ia]...\n");
	  exit (1);
	}

      int ib = b->ib[i];
      for (j = 0; j < nb[ib]; j ++)
	{
	  if (bond[ib][j] == -1)
	    {
	      bond[ib][j] = i;
	      //fprintf (stdout, "# set (ib) bond[%d][%d] = %d\n", ib, j, i);
	      break;
	    }
	}
      if (bond[ib][j] == -1)
	{
	  fprintf (stderr, "# BONDS_make_groups()"
		   " : something is wrong on bond[ib]...\n");
	  exit (1);
	}
    }
  // check
  for (i = 0; i < np; i ++)
    {
      //fprintf (stdout, "# bond[%d] :", i);
      int j;
      for (j = 0; j < nb[i]; j ++)
	{
	  //fprintf (stdout, " %d", bond[i][j]);
	  if (bond[i][j] == -1)
	    {
	      fprintf (stderr, "# BONDS_make_groups()"
		       " : bond[%d][%d] is not assigned.\n",
		       i, j);
	      exit (1);
	    }
	}
      //fprintf (stdout, "\n");
    }


  struct loop_bonds *lb = loop_bonds_init ();

  int *gid = (int *)malloc (np * sizeof (int));
  CHECK_MALLOC (gid, "BONDS_make_groups");
  for (i = 0; i < np; i ++)
    {
      gid[i] = -1;
    }

  int ig = 1;
  for (i = 0; i < np; i ++)
    {
      if (gid[i] != -1) continue;

      if (nb[i] == 0)
	{
	  // single particle
	  gid[i] = 0;
	  continue;
	}

      gid[i] = ig;
      BONDS_GROUPS_gid_check (b,
			      np, gid,
			      nb, bond,
			      i, ig,
			      lb);
      ig ++;
    }
  // check
  for (i = 0; i < np; i ++)
    {
      //fprintf (stdout, "# gid[%d] = %d\n", i, gid[i]);
      if (gid[i] == -1)
	{
	  fprintf (stderr, "# undefined particle %d\n", i);
	}
    }
  /*
  // check
  fprintf (stdout, "# lb->n = %d\n", lb->n);
  for (i = 0; i < lb->n; i ++)
    {
      fprintf (stdout, "# lb->b[%d] = %d\n", i, lb->b[i]);
    }
  */

  // make BONDS_GROUPS
  struct BONDS_GROUPS *gs = BONDS_GROUPS_init ();
  for (i = 0; i < np; i ++)
    {
      if (gid[i] == -1)
	{
	  continue;
	}
      else if (gid[i] == 0)
	{
	  // single particle
	  BONDS_GROUPS_add (gs, 1, &i, NULL);
	}
      else
	{
	  int nn;
	  int *ip = BONDS_GROUPS_gid_pickup (np, gid, gid[i], &nn);
	  /*
	  fprintf (stdout, "# np = %d\n", nn);
	  int ii;
	  for (ii = 0; ii < nn; ii ++)
	    {
	      fprintf (stdout, "# ip[%d] = %d\n", ii, ip[ii]);
	    }
	  exit (1);
	  */

	  // make bonds[]
	  int *bonds = (int *)malloc (sizeof (int) * (nn - 1));
	  CHECK_MALLOC (bonds, "BONDS_make_groups");

	  int k = 0;
	  int j;
	  for (j = 0; j < b->n; j ++)
	    {
	      if (check_ip (lb->n, lb->b, j) == 0 &&
		  (check_ip (nn, ip, b->ia[j]) == 1 ||
		   check_ip (nn, ip, b->ib[j]) == 1))
		{
		  if (k == (nn - 1))
		    {
		      fprintf (stderr, "too many bonds\n");
		      exit (1);
		    }
		  /*
		  fprintf (stdout, "# bond %d is the %d-th member\n",
			   j, k);
		  */
		  bonds[k] = j;
		  k ++;
		}
	    }
	  if (k != (nn - 1))
	    {
	      fprintf (stderr, "too few bonds, isn't it?"
		       " %d for %d particles\n", k, nn);
	      exit (1);
	    }

	  BONDS_GROUPS_add (gs, nn, ip, bonds);

	  free (ip);
	  free (bonds);
	}
    }


  // house keeping
  loop_bonds_free (lb);
  free (gid);
  for (i = 0; i < np; i ++)
    {
      free (bond[i]);
    }
  free (bond);
  free (nb);

  return (gs);
}


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
		     double *conn)
{
  if (g->np == 1)
    {
      int ix = g->ip[0] * 3;
      conn [0] = pos[ix  ];
      conn [1] = pos[ix+1];
      conn [2] = pos[ix+2];
    }
  else
    {
      int i;
      int ix;
      for (i = 0; i < g->np - 1; i ++)
	{
	  ix = i * 3;
	  int ibond = g->bonds[i];
	  int ia3 = b->ia[ibond] * 3;
	  int ib3 = b->ib[ibond] * 3;
	  conn[ix  ] = pos[ia3  ] - pos[ib3  ];
	  conn[ix+1] = pos[ia3+1] - pos[ib3+1];
	  conn[ix+2] = pos[ia3+2] - pos[ib3+2];
	}

      // calc COM of the group
      double cx = 0.0;
      double cy = 0.0;
      double cz = 0.0;
      for (i = 0; i < g->np; i ++)
	{
	  ix = g->ip[i] * 3;
	  cx += pos[ix  ];
	  cy += pos[ix+1];
	  cz += pos[ix+2];
	}
      cx /= (double)g->np;
      cy /= (double)g->np;
      cz /= (double)g->np;

      i = g->np - 1;
      ix = i * 3;
      conn[ix  ] = cx;
      conn[ix+1] = cy;
      conn[ix+2] = cz;
    }
}


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
		   double *conn)
{
  double *conn_i = conn;
  int i;
  for (i = 0; i < gs->n; i ++)
    {
      struct BONDS_GROUP *g = gs->group[i];
      BONDS_pos_to_conn_1 (b, g, pos, conn_i);
      conn_i += g->np * 3;
    }
}

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
		     double *pos)
{
  if (g->np == 1)
    {
      int ix = g->ip[0] * 3;
      pos[ix  ] = conn [0];
      pos[ix+1] = conn [1];
      pos[ix+2] = conn [2];
    }
  else
    {
      int ax;
      ax = b->ia[g->bonds[0]] * 3;
      pos[ax  ] = 0.0;
      pos[ax+1] = 0.0;
      pos[ax+2] = 0.0;

      double cx = 0.0;
      double cy = 0.0;
      double cz = 0.0;

      int i;
      for (i = 0; i < g->np - 1; i ++)
	{
	  int ix = i * 3;
	  int ibond = g->bonds[i];
	  ax     = b->ia[ibond] * 3;
	  int bx = b->ib[ibond] * 3;

	  pos[bx  ] = pos[ax  ] - conn[ix  ];
	  pos[bx+1] = pos[ax+1] - conn[ix+1];
	  pos[bx+2] = pos[ax+2] - conn[ix+2];

	  cx += pos[bx  ];
	  cy += pos[bx+1];
	  cz += pos[bx+2];
	}
      cx /= (double)g->np;
      cy /= (double)g->np;
      cz /= (double)g->np;

      // check
      double cx_ = 0.0;
      double cy_ = 0.0;
      double cz_ = 0.0;
      for (i = 0; i < g->np; i ++)
	{
	  ax = g->ip[i] * 3;
	  cx_ += pos[ax  ];
	  cy_ += pos[ax+1];
	  cz_ += pos[ax+2];
	}
      cx_ /= (double)g->np;
      cy_ /= (double)g->np;
      cz_ /= (double)g->np;

      double tiny = 1.0e-10;
      if (fabs (cx - cx_) > tiny ||
	  fabs (cy - cy_) > tiny ||
	  fabs (cz - cz_) > tiny)
	{
	  fprintf (stderr, "# COM differ! (%e %e %e) != (%e %e %e)\n",
		   cx, cy, cz,
		   cx_, cy_, cz_);
	  exit (1);
	}

      // adjust COM
      int ix = (g->np - 1) * 3;
      cx = conn[ix  ] - cx;
      cy = conn[ix+1] - cy;
      cz = conn[ix+2] - cz;
      for (i = 0; i < g->np; i ++)
	{
	  ax = g->ip[i] * 3;
	  pos[ax  ] += cx;
	  pos[ax+1] += cy;
	  pos[ax+2] += cz;
	}
    }
}

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
		   double *pos)
{
  const double *conn_i = conn;
  int i;
  for (i = 0; i < gs->n; i ++)
    {
      struct BONDS_GROUP *g = gs->group[i];
      BONDS_conn_to_pos_1 (b, g, conn_i, pos);
      conn_i += g->np * 3;
    }
}
