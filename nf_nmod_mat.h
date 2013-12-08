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

#ifndef NF_NMOD_MAT_H
#define NF_NMOD_MAT_H

#include "ulong_extras.h"
#include "flint.h"
#include "rank_profile.h"
#include "nmod_poly.h"
#include "fq.h"
#include "fq_nmod_mat.h"
#include "nf_nmod.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    mp_limb_t * entries;
    slong r;
    slong c;
    mp_limb_t ** rows;
    mp_limb_t *** poly_coeffs;
}
nf_nmod_mat_struct;

typedef nf_nmod_mat_struct nf_nmod_mat_t[1];

#define nf_nmod_mat_entry(mat,n,i,j) ((mat)->poly_coeffs[(n)][(i)][(j)])
#define nf_nmod_mat_nrows(mat) ((mat)->r)
#define nf_nmod_mat_ncols(mat) ((mat)->c)

void nf_nmod_mat_init(nf_nmod_mat_t mat, slong rows, slong cols, const nf_nmod_ctx_t ctx);
void nf_nmod_mat_init_set(nf_nmod_mat_t mat, const nf_nmod_mat_t src, const nf_nmod_ctx_t ctx);
void nf_nmod_mat_clear(nf_nmod_mat_t mat, const nf_nmod_ctx_t ctx);

void nf_nmod_mat_set(nf_nmod_mat_t mat1, const nf_nmod_mat_t mat2, const nf_nmod_ctx_t ctx);

void nf_nmod_mat_init_coeff_mat(nmod_mat_t B, const nf_nmod_mat_t A, slong n, const nf_nmod_ctx_t ctx);

void nf_nmod_mat_clear_coeff_mat(nmod_mat_t B, const nf_nmod_ctx_t ctx);


void nf_nmod_mat_zero(nf_nmod_mat_t A, const nf_nmod_ctx_t ctx);


void nf_nmod_mat_decompose_components(nmod_mat_t * fp_comps, fq_nmod_mat_t * fq_comps, fq_nmod_ctx_t * fq_ctxs, const nf_nmod_mat_t A, const nf_nmod_ctx_t ctx);

void nf_nmod_mat_reconstruct_components(nf_nmod_mat_t A, const nmod_mat_t * fp_comps, const fq_nmod_mat_t * fq_comps, const fq_nmod_ctx_t * fq_ctxs, const nf_nmod_ctx_t ctx);

void _nf_nmod_mat_init_components(nmod_mat_t * fp_comps, fq_nmod_mat_t * fq_comps, fq_nmod_ctx_t * fq_ctxs, const nf_nmod_mat_t A, const nf_nmod_ctx_t ctx);

void _nf_nmod_mat_clear_components(nmod_mat_t * fp_comps, fq_nmod_mat_t * fq_comps, fq_nmod_ctx_t * fq_ctxs, const nf_nmod_ctx_t ctx);


void nf_nmod_mat_rref_components(nf_nmod_mat_t B, const nf_nmod_mat_t A, const nf_nmod_ctx_t ctx);


slong nf_nmod_mat_rank_profile(rank_profile_t rk_prof, const nf_nmod_mat_t A, const nf_nmod_ctx_t ctx);

slong _nf_nmod_mat_first_non_zero_entry(mp_limb_t * entry, const nf_nmod_mat_t A, slong r, slong c, const nf_nmod_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
