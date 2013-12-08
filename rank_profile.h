/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2013 Martin Raum
 
******************************************************************************/

#ifndef RANK_PROFILE_H
#define RANK_PROFILE_H

#include "ulong_extras.h"
#include "flint.h"

#ifdef __cplusplus
 extern "C" {
#endif

typedef struct
{
  slong * profile;
  slong length;
}
rank_profile_struct;

typedef rank_profile_struct rank_profile_t[1];


#define rank_profile_entry(rk_prof, i) (rk_prof->profile[(i)])

void rank_profile_init(rank_profile_t rk_prof, slong len);
void rank_profile_clear(rank_profile_t rk_prof);

/* todo: implement */
void rank_profile_set(rank_profile_t rk_prof1, const rank_profile_t rk_prof2);

int rank_profile_cmp(const rank_profile_t rk_prof1, const rank_profile_t rk_prof2);

#ifdef __cplusplus
}
#endif

#endif
