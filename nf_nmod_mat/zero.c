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

void
nf_nmod_mat_zero(nf_nmod_mat_t mat, const nf_nmod_ctx_t ctx)
{
    slong i;

    if (mat->c < 1)
        return;

    for (i = 0; i < mat->r * ctx->deg ; ++i)
        _nmod_vec_zero(mat->rows[i], mat->c);
}
