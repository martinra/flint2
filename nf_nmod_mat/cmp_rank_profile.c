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

#include "nf_nmod_mat.h"

int
nf_nmod_mat_cmp_rank_profile(const nf_nmod_mat_rank_profile_t rk_prof1, const nf_nmod_mat_rank_profile_t rk_prof2, slong r)
{
  int i;

  for (i = 0; i < r; ++i)
    {
      if (rk_prof1[i] < rk_prof2[i])
	return -1;
      else if (rk_prof1[i] > rk_prof2[i])
	return 1;
    }

  return 0;
}









