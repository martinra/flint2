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

#include "nfz_mat.h"

slong
nfz_mat_rank(const nfz_mat_t A, const nfz_ctx_t ctx)
{
    nfz_mat_t tmp;
    fmpz_t den;
    slong rank;

    if (nfz_mat_is_empty(A, ctx))
        return 0;

    nfz_mat_init(tmp, A->r, A->c, ctx);
    fmpz_init(den);

    rank = nfz_mat_rref_mod(tmp, den, A, ctx);

    nfz_mat_clear(tmp, ctx);
    fmpz_clear(den);

    return rank;
}
