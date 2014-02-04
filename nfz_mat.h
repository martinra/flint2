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
    nfz * entries;
    slong r;
    slong c;
    nfz ** rows;
    nfz *** poly_coeffs;
}
nfz_mat_struct;

typedef nfz_mat_struct nfz_mat_t[1];


#define nfz_mat_entry(mat,n,i,j) ((mat)->poly_coeffs[(n)][(i)] + (j))
#define nfz_mat_nrows(mat) ((mat)->r)
#define nfz_mat_ncols(mat) ((mat)->c)


void nfz_mat_init(nfz_mat_t mat, slong rows, slong cols, const nfz_ctx_t ctx);
void nfz_mat_init_set(nfz_mat_t mat, const nfz_mat_t src, const nfz_ctx_t ctx);
void nfz_mat_swap(nfz_mat_t mat1, nfz_mat_t mat2, const nfz_ctx_t ctx);
void nfz_mat_set(nfz_mat_t mat1, const nfz_mat_t mat2, const nfz_ctx_t ctx);
void nfz_mat_clear(nfz_mat_t mat, const nfz_ctx_t ctx);

int nfz_mat_equal(const nfz_mat_t mat1, const nfz_mat_t mat2, const nfz_ctx_t ctx);
int nfz_mat_is_zero(const nfz_mat_t mat, const nfz_ctx_t ctx);

static __inline__ int
nfz_mat_is_empty(const nfz_mat_t mat, const nfz_ctx_t ctx)
{
    return (mat->r == 0) || (mat->c == 0);
}

static __inline__ int
nfz_mat_is_square(const nfz_mat_t mat, const nfz_ctx_t ctx)
{
    return (mat->r == mat->c);
}

void nfz_mat_zero(nfz_mat_t mat, const nfz_ctx_t ctx);
void nfz_mat_one(nfz_mat_t mat, const nfz_ctx_t ctx);


// int nfz_mat_fprint(FILE * file, const nfz_mat_t mat);


void nfz_mat_transpose(nfz_mat_t B, const nfz_mat_t A, const nfz_ctx_t ctx);


/* Addition and subtraction */

void nfz_mat_add(nfz_mat_t C, const nfz_mat_t A, const nfz_mat_t B,
		 const nfz_ctx_t ctx);
void nfz_mat_sub(nfz_mat_t C, const nfz_mat_t A, const nfz_mat_t B, const nfz_ctx_t ctx);
void nfz_mat_neg(nfz_mat_t B, const nfz_mat_t A, const nfz_ctx_t ctx);

/* Scalar operations */
void nfz_mat_scalar_mul_nfz(nfz_mat_t B, const nfz_mat_t A, const nfz_t c, const nfz_ctx_t ctx);
void nfz_mat_scalar_mul_fmpz(nfz_mat_t B, const nfz_mat_t A, const fmpz_t c, const nfz_ctx_t ctx);
void nfz_mat_scalar_mul_si(nfz_mat_t B, const nfz_mat_t A, slong c, const nfz_ctx_t ctx);
void nfz_mat_scalar_mul_ui(nfz_mat_t B, const nfz_mat_t A, ulong c, const nfz_ctx_t ctx);

void nfz_mat_scalar_addmul_nfz(nfz_mat_t B, const nfz_mat_t A, const nfz_t c, const nfz_ctx_t ctx);
void nfz_mat_scalar_addmul_fmpz(nfz_mat_t B, const nfz_mat_t A, const fmpz_t c, const nfz_ctx_t ctx);
void nfz_mat_scalar_addmul_si(nfz_mat_t B, const nfz_mat_t A, slong c, const nfz_ctx_t ctx);
void nfz_mat_scalar_addmul_ui(nfz_mat_t B, const nfz_mat_t A, ulong c, const nfz_ctx_t ctx);

void nfz_mat_scalar_submul_nfz(nfz_mat_t B, const nfz_mat_t A, const nfz_t c, const nfz_ctx_t ctx);
void nfz_mat_scalar_submul_fmpz(nfz_mat_t B, const nfz_mat_t A, const fmpz_t c, const nfz_ctx_t ctx);
void nfz_mat_scalar_submul_si(nfz_mat_t B, const nfz_mat_t A, slong c, const nfz_ctx_t ctx);
void nfz_mat_scalar_submul_ui(nfz_mat_t B, const nfz_mat_t A, ulong c, const nfz_ctx_t ctx);

void nfz_mat_scalar_divexact_fmpz(nfz_mat_t B, const nfz_mat_t A, const fmpz_t c, const nfz_ctx_t ctx);
void nfz_mat_scalar_divexact_si(nfz_mat_t B, const nfz_mat_t A, slong c, const nfz_ctx_t ctx);
void nfz_mat_scalar_divexact_ui(nfz_mat_t B, const nfz_mat_t A, ulong c, const nfz_ctx_t ctx);

/* Multiplication */

void nfz_mat_mul(nfz_mat_t C, const nfz_mat_t A, const nfz_mat_t B, const nfz_ctx_t ctx);

void nfz_mat_sqr(nfz_mat_t B, const nfz_mat_t A, const nfz_ctx_t ctx);

void nfz_mat_pow(nfz_mat_t B, const nfz_mat_t A, ulong exp, const nfz_ctx_t ctx);


/* Modular gaussian elimination *********************************************/

slong nfz_mat_rref_mod(nfz_mat_t B, nfz_t den, const nfz_mat_t A, const nfz_ctx_t ctx);

slong _nfz_mat_rref_mod_prime_generator(nfz_mat_t B, nfz_t den, const nfz_mat_t A, const nfz_ctx_t ctx, int (*next_prime)(const mp_limb_t));

slong _nfz_mat_first_non_zero_entry(nfz_t * entry, const nfz_mat_t A, slong r, slong c, const nfz_ctx_t ctx);

void nfz_mat_coeff_bound(nfz_t bound, const nfz_mat_t A, const nfz_ctx_t ctx);


/* Trace ********************************************************************/

void nfz_mat_trace(nfz_t trace, const nfz_mat_t mat, const nfz_ctx_t ctx);

/* Determinant **************************************************************/

void nfz_mat_det(nfz_t det, const nfz_mat_t A, const nfz_ctx_t ctx);

/* Characteristic polynomial ************************************************/

void nfz_mat_charpoly(nfz_poly_t cp, const nfz_mat_t mat, const nfz_ctx_t ctx);

/* Rank *********************************************************************/

slong nfz_mat_rank(const nfz_mat_t A, const nfz_ctx_t ctx);

slong nfz_mat_rank_profile(rank_profile_t rk_prof, const nfz_mat_t A, const nfz_ctx_t ctx);

/* Nonsingular solving ******************************************************/

int nfz_mat_solve(nfz_mat_t X, nfz_t den,
        const nfz_mat_t A, const nfz_mat_t B, const nfz_ctx_t ctx);

int nfz_mat_solve_fflu(nfz_mat_t X, nfz_t den,
        const nfz_mat_t A, const nfz_mat_t B, const nfz_ctx_t ctx);

/* Nullspace ****************************************************************/

slong nfz_mat_nullspace(nfz_mat_t res, const nfz_mat_t mat, const nfz_ctx_t ctx);

/* Inverse ******************************************************************/

int nfz_mat_inv(nfz_mat_t B, nfz_t den, const nfz_mat_t A, const nfz_ctx_t ctx);

/* Modular reduction and reconstruction *************************************/

void nfz_mat_set_nmod_mat(nfz_mat_t A, const nmod_mat_t Amod, const nfz_ctx_t ctx);

void nfz_mat_get_nmod_mat(nf_nmod_mat_t B, const nfz_mat_t A, const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx);

void nfz_mat_CRT_nmod(nfz_mat_t out, const nfz_mat_t in1, const nfz_t m1, const nf_nmod_mat_t in2, int sign, const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
