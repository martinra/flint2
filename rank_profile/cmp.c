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

#include "rank_profile.h"

int
rank_profile_cmp(const rank_profile_t rk_prof1, const rank_profile_t rk_prof2)
{
  int i;

  for (i = 0; i < rk_prof1->length; ++i)
    {
      if (rk_prof1->profile[i] < rk_prof2->profile[i])
	return -1;
      else if (rk_prof1->profile[i] > rk_prof2->profile[i])
	return 1;
    }

  return 0;
}
