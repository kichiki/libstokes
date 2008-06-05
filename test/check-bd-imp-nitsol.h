/* header file for check-bd-im-nitsol.c --
 * test code for bd-imp-nitsol.c
 * Copyright (C) 2008 Kengo Ichiki <kichiki@users.sourceforge.net>
 * $Id: check-bd-imp-nitsol.h,v 1.1 2008/06/05 03:24:21 kichiki Exp $
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
#ifndef	_CHECK_BD_IMP_NITSOL_H_
#define	_CHECK_BD_IMP_NITSOL_H_


int
check_BD_imp_NITSOL (int version, int np,
		     int flag_lub, int flag_mat, int flag_Q, double dt,
		     int verbose, double tiny);


#endif /* !_CHECK_BD_IMP_NITSOL_H_ */
