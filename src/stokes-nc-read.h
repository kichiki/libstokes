/* header file for stokes-nc-read.c --
 * NetCDF interface for libstokes
 * Copyright (C) 2006 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: stokes-nc-read.h,v 5.1 2006/10/19 02:45:38 ichiki Exp $
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


void
stokes_nc_print_actives (struct stokes_nc * nc,
			 FILE * out);
struct stokes_nc *
stokes_nc_open (const char * filename);

void
stokes_nc_get_data0 (struct stokes_nc * nc,
		     const char * name,
		     double * x);
void
stokes_nc_get_data (struct stokes_nc * nc,
		    const char * name,
		    int step,
		    double * x);

#endif /* !_STOKES_READ_NC_H_ */
