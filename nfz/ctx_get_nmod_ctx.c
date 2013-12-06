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
  slong deg = nmod_poly_degree(modulus);
  slong big_size = 2 * deg - 1;

  modulus_nmod modulus_nmod;
  nmod_t mod;


  fmpz_poly_get_nmod_poly(modulus_nmod, ctx->modulus);
  mod = modulus_nmod->mod;

  fmpz_mat_det(det, ctx->ev_mat);
  fmpz_mod_ui(det, det, n);

  if (!fmpz_is_zero(det))
    {
      nmod_mat_t ev_mat;
      nmod_mat_t int_mat;
      mp_limb_t int_den;

      nmod_mat_init(ev_mat, ctx->deg, ctx->deg, n);
      nmod_mat_init(int_mat, ctx->deg, ctx->deg, n);

      fmpz_mat_get_nmod_mat(ev_mat, ctx->ev_mat);
      fmpz_mat_get_nmod_mat(int_mat, ctx->int_mat);
      int_den = fmpz_get_nmod(ctx->int_den);
      int_den = nmod_inv(int_den, mod);
      nmod_mat_scalar_mul(int_mat, int_mat, int_den);

      _nf_nmod_ctx_init_with_eval(ctx_nmod, modulus_nmod, ev_mat, int_mat, ctx->var);

      nmod_mat_clear(ev_mat);
      nmod_mat_clear(int_mat);
    }
  else if (big_size <= n)
    nf_nmod_ctx_init(ctx_nmod, modulus_nmod, ctx->var);
  else 
    _nf_nmod_ctx_init_with_eval(ctx_nmod, modulus_nmod, NULL, NULL, ctx->var);
}
