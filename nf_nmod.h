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

#ifndef NF_NMOD_H
#define NF_NMOD_H

#include "ulong_extras.h"
#include "flint.h"

typedef struct
{
  nmod_poly_t modulus;
  int separable;

  slong nfp;
  slong nfq;
  mp_limb_t * fp_moduli;
  nmod_poly_t * fq_moduli;

  nmod_mat_t decomp_mat;
  nmod_mat_t reconst_mat;

  nmod_mat_t ev_mat;
  nmod_mat_t int_mat;

  char * var;
}
nf_nmod_ctx_struct;

typedef nf_nmod_ctx_struct nf_nmod_ctx_t[1];

/* Memory managment  *********************************************************/

void nf_nmod_ctx_init(nf_nmod_ctx_t ctx, const nmod_poly_t modulus, const char *var);
void _nf_nmod_ctx_init_with_eval(nf_nmod_ctx_t ctx, const nmod_poly_t modulus, const nmod_mat_t ev_mat, const nmod_mat_t int_mat, const char *var);
void nf_nmod_ctx_clear(nf_nmod_ctx_t ctx);

int
nf_nmod_ctx_is_separable(const nf_nmod_ctx_t ctx)
{
  return ctx->separable;
}

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif
