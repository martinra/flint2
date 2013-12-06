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
    mp_limb_t * entries;
    slong r;
    slong c;
    mp_limb_t ** rows;
    mp_limb_t *** poly_coeffs;
}
nf_nmod_mat_struct;

typedef nf_nmod_mat_struct nf_nmod_mat_t[1];

#define nf_nmod_mat_entry(mat,e,i,j) ((mat)->poly_coeffs[(e)][(i)][(j)])
#define nf_nmod_mat_nrows(mat) ((mat)->r)
#define nf_nmod_mat_ncols(mat) ((mat)->c)

void nf_nmod_mat_init(nf_nmod_mat_t mat, slong rows, slong cols, const nf_nmod_ctx_t ctx);
void nf_nmod_mat_init_set(nf_nmod_mat_t mat, const nf_nmod_mat_t src, const nf_nmod_ctx_t ctx);
void nf_nmod_mat_set(nf_nmod_mat_t mat1, const nf_nmod_mat_t mat2, const nf_nmod_ctx_t ctx);
void nf_nmod_mat_clear(nf_nmod_mat_t mat, const nf_nmod_ctx_t ctx);

void nf_nmod_rref_components(nf_nmod_mat_t B, const nf_nmod_mat_t A, const nf_nmod_ctx_t ctx);


void _nf_nmod_reduction_mat(nmod_mat_t mat, nmod_poly_t modulus, ulong deg_bd);

void _nf_nmod_mat_init_components(nmod_mat_t * fp_comps, fq_nmod_mat_t * fq_comps, fq_ctx_t * fq_ctxs, const nf_nmod_mat_t A, const nf_nmod_ctx_t ctx);

void _nf_nmod_mat_clear_components(nmod_mat_t * fp_comps, fq_nmod_mat_t * fq_comps, fq_ctx_t * fq_ctxs, const nf_nmod_ctx_t ctx);


void nf_nmod_mat_decompose_components(nmod_mat_t * fp_comps, fq_nmod_mat_t * fq_comps, fq_ctx_t * fq_ctxs, const nf_nmod_mat_t A, const nf_nmod_ctx_t ctx);

void nf_nmod_mat_reconstruct_components(nf_nmod_mat_t A, const nmod_mat_t * fp_comps, const fq_nmod_mat_t * fq_comps, const fq_ctx_t fq_ctx, const nf_nmod_ctx_t ctx);


#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
}
#endif
