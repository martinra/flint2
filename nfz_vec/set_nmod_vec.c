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

#include "flint.h"
#include "nf_nmod.h"
#include "nfz_vec.h"

void _nfz_vec_set_nmod_vec(nfz * res, 
			   const nf_nmod * src, slong len,
			   const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_set_nmod(res + i, src + i, ctx_nmod, ctx);
}
