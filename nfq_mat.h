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

    Copyright (C) 2014 Martin Raum
 
******************************************************************************/

#ifndef NFQ_MAT_H
#define NFQ_MAT_H

#include "ulong_extras.h"
#include "flint.h"
#include "fmpq.h"
#include "rank_profile.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  fmpq * entries;
  slong r;
  slong c;
  fmpq ** rows;
  fmpq *** poly_coeffs;
} nfq_mat_struct;

typedef nfq_mat_struct nfq_mat_t[1];

#ifdef __cplusplus
}
#endif

#endif
