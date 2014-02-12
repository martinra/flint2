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
#include "fmpz_mat.h"
#include "fmpq_mat.h"
#include "nf_nmod_mat.h"
#include "nfz.h"

#ifdef __cplusplus
extern "C" {
#endif

/******************************************************************************/
/* Matrices over number fields ************************************************/
/******************************************************************************/

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

/* Memory managment ***********************************************************/

#define nfz_mat_entry(mat,n,i,j) ((mat)->poly_coeffs[(n)][(i)] + (j))

void nfz_mat_entry_nfz(nfz_t e, const nfz_mat_t A, slong r, slong c,
		       const nfz_ctx_t ctx);

#define nfz_mat_nrows(mat) ((mat)->r)

#define nfz_mat_ncols(mat) ((mat)->c)

void nfz_mat_init(nfz_mat_t mat, slong rows, slong cols, const nfz_ctx_t ctx);

void nfz_mat_init_set(nfz_mat_t mat, const nfz_mat_t src, const nfz_ctx_t ctx);

void nfz_mat_clear(nfz_mat_t mat, const nfz_ctx_t ctx);


/* Basic properties ***********************************************************/

static __inline__ int
nfz_mat_is_empty(const nfz_mat_t mat, const nfz_ctx_t ctx)
{
    return (mat->r == 0) || (mat->c == 0);
};

static __inline__ int
nfz_mat_is_square(const nfz_mat_t mat, const nfz_ctx_t ctx)
{
    return (mat->r == mat->c);
};

/* Evaluation and interpolation ***********************************************/

static __inline__
fmpz_mat_struct * _nfz_mat_eval_init(long r, long c, const nfz_ctx_t ctx)
{
  fmpz_mat_struct * evl =
    (fmpz_mat_struct *) flint_malloc(ctx->evl_mat->r * sizeof(fmpz_mat_struct));
  for (long n = 0; n < ctx->evl_mat->r; ++n)
    fmpz_mat_init(evl + n, r, c);
  return evl;
};

static __inline__
void _nfz_mat_eval_clear(fmpz_mat_struct * evl, const nfz_ctx_t ctx)
{
  for (long n = 0; n < ctx->evl_mat->r; ++n)
    fmpz_mat_clear(evl + n);
  flint_free(evl);
};

void _nfz_mat_eval(fmpz_mat_struct * evl, const nfz_mat_t A,
		   const nfz_ctx_t ctx);

void _nfz_mat_interpolate(nfz_mat_t A, const fmpz_mat_struct * evl,
			  const nfz_ctx_t ctx);

// todo: implement
void _nfz_mat_eval_entry(fmpz * evl, nfz_mat_t mat, slong r, slong c,
			 const nfz_ctx_t ctx);

// todo: implement
void _nfz_mat_eval_row(fmpz ** evl, nfz_mat_t mat, slong r,
		       slong start_col, slong end_col, const nfz_ctx_t ctx);

// todo: implement
// this is used indirectly in fflu
void _nfz_mat_interpolate_sub_row(nfz_mat_t mat, fmpz ** row_evl, slong r,
				  slong start_col, slong end_col,
				  const nfz_ctx_t ctx);

/* Assignment and basic manipulation ******************************************/

void nfz_mat_set(nfz_mat_t mat1, const nfz_mat_t mat2, const nfz_ctx_t ctx);

void nfz_mat_swap(nfz_mat_t mat1, nfz_mat_t mat2, const nfz_ctx_t ctx);

void nfz_mat_zero(nfz_mat_t mat, const nfz_ctx_t ctx);

void nfz_mat_one(nfz_mat_t mat, const nfz_ctx_t ctx);

void nfz_mat_neg(nfz_mat_t B, const nfz_mat_t A, const nfz_ctx_t ctx);

/* Comparison *****************************************************************/

int nfz_mat_equal(const nfz_mat_t mat1, const nfz_mat_t mat2, const nfz_ctx_t ctx);

int nfz_mat_is_zero(const nfz_mat_t mat, const nfz_ctx_t ctx);

/* Output *********************************************************************/

// int nfz_mat_fprint(FILE * file, const nfz_mat_t mat);

/* Windows ********************************************************************/

void nfz_mat_init_window(nfz_mat_t window, const nfz_mat_t mat,
			 slong r1, slong c1, slong r2, slong c2,
			 const nfz_ctx_t ctx);

void nfz_mat_clear_window(nfz_mat_t window, const nfz_ctx_t ctx);

void nfz_mat_init_window_fmpz(fmpz_mat_t window, const nfz_mat_t mat, long n,
			      const nfz_ctx_t ctx);

void nfz_mat_clear_window_fmpz(fmpz_mat_t window, const nfz_ctx_t ctx);

/* Transpose ******************************************************************/

void nfz_mat_transpose(nfz_mat_t B, const nfz_mat_t A, const nfz_ctx_t ctx);

/* Addition and subtraction ***************************************************/

void nfz_mat_add(nfz_mat_t C, const nfz_mat_t A, const nfz_mat_t B,
		 const nfz_ctx_t ctx);
void nfz_mat_sub(nfz_mat_t C, const nfz_mat_t A, const nfz_mat_t B,
		 const nfz_ctx_t ctx);

/* Scalar operations **********************************************************/

void nfz_mat_scalar_mul_nfz(nfz_mat_t B, const nfz_mat_t A, const nfz_t c,
			    const nfz_ctx_t ctx);

void nfz_mat_scalar_mul_fmpz(nfz_mat_t B, const nfz_mat_t A, const fmpz_t c,
			     const nfz_ctx_t ctx);

void nfz_mat_scalar_mul_si(nfz_mat_t B, const nfz_mat_t A, slong c,
			   const nfz_ctx_t ctx);

void nfz_mat_scalar_mul_ui(nfz_mat_t B, const nfz_mat_t A, ulong c,
			   const nfz_ctx_t ctx);

void nfz_mat_scalar_addmul_nfz(nfz_mat_t B, const nfz_mat_t A, const nfz_t c,
			       const nfz_ctx_t ctx);

void nfz_mat_scalar_addmul_fmpz(nfz_mat_t B, const nfz_mat_t A, const fmpz_t c,
				const nfz_ctx_t ctx);

void nfz_mat_scalar_addmul_si(nfz_mat_t B, const nfz_mat_t A, slong c,
			      const nfz_ctx_t ctx);

void nfz_mat_scalar_addmul_ui(nfz_mat_t B, const nfz_mat_t A, ulong c,
			      const nfz_ctx_t ctx);

void nfz_mat_scalar_submul_nfz(nfz_mat_t B, const nfz_mat_t A, const nfz_t c,
			       const nfz_ctx_t ctx);

void nfz_mat_scalar_submul_fmpz(nfz_mat_t B, const nfz_mat_t A, const fmpz_t c,
				const nfz_ctx_t ctx);

void nfz_mat_scalar_submul_si(nfz_mat_t B, const nfz_mat_t A, slong c,
			      const nfz_ctx_t ctx);

void nfz_mat_scalar_submul_ui(nfz_mat_t B, const nfz_mat_t A, ulong c,
			      const nfz_ctx_t ctx);

void nfz_mat_scalar_divexact_nfz(nfz_mat_t B, const nfz_mat_t A, const nfz_t c,
				 const nfz_ctx_t ctx);

void nfz_mat_scalar_divexact_fmpz(nfz_mat_t B, const nfz_mat_t A, const fmpz_t c,
				  const nfz_ctx_t ctx);

void nfz_mat_scalar_divexact_si(nfz_mat_t B, const nfz_mat_t A, slong c,
				const nfz_ctx_t ctx);

void nfz_mat_scalar_divexact_ui(nfz_mat_t B, const nfz_mat_t A, ulong c,
				const nfz_ctx_t ctx);

/* Multiplication *************************************************************/

void nfz_mat_mul(nfz_mat_t C, const nfz_mat_t A, const nfz_mat_t B,
		 const nfz_ctx_t ctx);

void nfz_mat_sqr(nfz_mat_t B, const nfz_mat_t A, const nfz_ctx_t ctx);

void nfz_mat_pow(nfz_mat_t B, const nfz_mat_t A, ulong exp, const nfz_ctx_t ctx);

/* Trace **********************************************************************/

void nfz_mat_trace(nfz_t trace, const nfz_mat_t A, const nfz_ctx_t ctx);

/* Determinant ****************************************************************/

void nfz_mat_det(nfz_t det, const nfz_mat_t A, const nfz_ctx_t ctx);

/* Characteristic polynomial **************************************************/

// void nfz_mat_charpoly(nfz_poly_t cp, const nfz_mat_t mat, const nfz_ctx_t ctx);

/* Rank ***********************************************************************/

slong nfz_mat_rank(const nfz_mat_t A, const nfz_ctx_t ctx);

slong nfz_mat_rank_profile(rank_profile_t rk_prof, const nfz_mat_t A,
			   const nfz_ctx_t ctx);

/* Permutations ***************************************************************/

static __inline__ void
nfz_mat_swap_rows(nfz_mat_t mat, slong * perm, slong r, slong s, const nfz_ctx_t ctx)
{
  if (r != s)
  {
    fmpz * u;
    slong t;

    if (perm)
    {
      t = perm[s];
      perm[s] = perm[r];
      perm[r] = t;
    }

    for (long n = 0; ctx->deg; ++n) {
      u = mat->poly_coeffs[n][s];
      mat->poly_coeffs[n][s] = mat->poly_coeffs[n][r];
      mat->poly_coeffs[n][r] = u; 
    }
  }
};

/* Gaussian elimination *******************************************************/

slong nfz_mat_find_pivot_any(const nfz_mat_t mat,
			     slong start_row, slong end_row, slong c,
			     const nfz_ctx_t ctx);

void _nfz_mat_row_submul_entry(nfz_mat_t mat,
			       slong r1, slong r2, slong start_col, slong end_col,
			       slong er, slong ec, const nfz_ctx_t ctx);
			 
slong nfz_mat_fflu(nfz_mat_t B, nfz_t den, slong * perm,
		   const nfz_mat_t A, int rank_check,
		   const nfz_ctx_t ctx);

// todo: implement
slong nfz_mat_rref(nfz_mat_t B, fmpz_t den, const nfz_mat_t A,
		   const nfz_ctx_t ctx);

/* Modular gaussian elimination ***********************************************/

slong nfz_mat_rref_mod(nfz_mat_t B, fmpz_t den, const nfz_mat_t A,
		       const nfz_ctx_t ctx);

slong _nfz_mat_rref_mod_prime_generator(nfz_mat_t B, fmpz_t den, const nfz_mat_t A,
					const nfz_ctx_t ctx,
					int (*next_prime)(const mp_limb_t));

slong _nfz_mat_first_non_zero_entry(nfz_t * entry, const nfz_mat_t A,
				    slong r, slong c, const nfz_ctx_t ctx);

void nfz_mat_coeff_bound(fmpz_t bound, const nfz_mat_t A, const nfz_ctx_t ctx);

/* Nonsingular solving ********************************************************/

// todo: implement
int nfz_mat_solve_fflu(nfz_mat_t X, nfz_t den,
        const nfz_mat_t A, const nfz_mat_t B, const nfz_ctx_t ctx);

static __inline__
int nfz_mat_solve(nfz_mat_t X, nfz_t den,
        const nfz_mat_t A, const nfz_mat_t B, const nfz_ctx_t ctx)
{
  return nfz_mat_solve_fflu(X, den, A, B, ctx);
};

/* Nullspace ******************************************************************/

// todo: implement
slong nfz_mat_nullspace(nfz_mat_t res, const nfz_mat_t mat, const nfz_ctx_t ctx);

/* Inverse ********************************************************************/

// todo: implement
int nfz_mat_inv(nfz_mat_t B, nfz_t den, const nfz_mat_t A, const nfz_ctx_t ctx);

/* Modular reduction and reconstruction ***************************************/

void nfz_mat_set_nmod_mat(nfz_mat_t A, const nf_nmod_mat_t Amod,
			  const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx);

void nfz_mat_get_nmod_mat(nf_nmod_mat_t B, const nfz_mat_t A,
			  const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx);

void nfz_mat_CRT_nmod(nfz_mat_t out, const nfz_mat_t in1, const fmpz_t m1,
		      const nf_nmod_mat_t in2, int sign,
		      const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
