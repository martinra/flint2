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

#include "nfz.h"

void
nfz_ctx_get_nmod_ctx(nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx, ulong n)
{
  slong deg = fmpz_poly_degree(ctx->modulus);

  nmod_poly_t modulus_nmod;
  nmod_t mod;

  fmpz_poly_get_nmod_poly(modulus_nmod, ctx->modulus);
  mod = modulus_nmod->mod;

  if (2 * deg - 2 < n) {
    nmod_mat_t evl_mat;
    nmod_mat_t intrpl_mat;
    mp_limb_t intrpl_den;

    nmod_mat_init(evl_mat, ctx->deg, ctx->deg, n);
    nmod_mat_init(intrpl_mat, ctx->deg, ctx->deg, n);

    fmpz_mat_get_nmod_mat(evl_mat, ctx->evl_mat);
    fmpz_mat_get_nmod_mat(intrpl_mat, ctx->intrpl_mat);
    intrpl_den = fmpz_fdiv_ui(ctx->intrpl_den, n);
    intrpl_den = nmod_inv(intrpl_den, mod);
    nmod_mat_scalar_mul(intrpl_mat, intrpl_mat, intrpl_den);

    _nf_nmod_ctx_init_with_eval(ctx_nmod, modulus_nmod, evl_mat, intrpl_mat, ctx->var);
  }
  else 
    _nf_nmod_ctx_init_with_eval(ctx_nmod, modulus_nmod, NULL, NULL, ctx->var);
}
