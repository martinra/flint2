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

#ifndef NFZ_MAT_H
#define NFZ_MAT_H

#include "ulong_extras.h"
#include "flint.h"
#include "rank_profile.h"
#include "fmpq_mat.h"
#include "nf_nmod_mat.h"
#include "nfz.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    fmpz * entries;
    slong r;
    slong c;
    fmpz ** rows;
    fmpz *** poly_coeffs;
}
nfz_mat_struct;

typedef nfz_mat_struct nfz_mat_t[1];


#define nfz_mat_entry(mat,n,i,j) ((mat)->poly_coeffs[(n)][(i)] + (j))
#define nfz_mat_nrows(mat) ((mat)->r)
#define nfz_mat_ncols(mat) ((mat)->c)


void nfz_mat_init(nfz_mat_t mat, slong rows, slong cols, const nfz_ctx_t ctx);
void nfz_mat_init_set(nfz_mat_t mat, const nfz_mat_t src, const nfz_ctx_t ctx);
void nfz_mat_clear(nfz_mat_t mat, const nfz_ctx_t ctx);


void nfz_mat_set(nfz_mat_t mat1, const nfz_mat_t mat2, const nfz_ctx_t ctx);


slong nfz_mat_rref_mod(nfz_mat_t B, fmpz_t den, const nfz_mat_t A, const nfz_ctx_t ctx);

slong _nfz_mat_rref_mod_prime_generator(nfz_mat_t B, fmpz_t den, const nfz_mat_t A, const nfz_ctx_t ctx, int (*next_prime)(const mp_limb_t));

slong nfz_mat_rank_profile(rank_profile_t rk_prof, const nfz_mat_t A, const nfz_ctx_t ctx);

slong _nfz_mat_first_non_zero_entry(fmpz_t * entry, const nfz_mat_t A, slong r, slong c, const nfz_ctx_t ctx);


void nfz_mat_get_nmod_mat(nf_nmod_mat_t B, const nfz_mat_t A, const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx);

void nfz_mat_CRT_nmod(nfz_mat_t out, const nfz_mat_t in1, const fmpz_t m1, const nf_nmod_mat_t in2, int sign, const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx);


void nfz_mat_coeff_bound(fmpz_t bound, const nfz_mat_t A, const nfz_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
