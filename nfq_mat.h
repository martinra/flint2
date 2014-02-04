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

#ifndef NFQ_MAT_H
#define NFQ_MAT_H

#include "ulong_extras.h"
#include "flint.h"
#include "nfq.h"
#include "nfz_mat.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  nfq * entries;
  slong r;
  slong c;
  nfq ** rows;
} nfq_mat_struct;

typedef nfq_mat_struct nfq_mat_t[1];

#define nfq_mat_entry(mat,i,j) ((mat)->rows[(i)] + (j))
#define nfq_mat_entry_num(mat,i,j) ((nfz *)(&((*nfq_mat_entry(mat,i,j)).num)))
#define nfq_mat_entry_den(mat,i,j) ((nfz *)(&((*nfq_mat_entry(mat,i,j)).den)))

#define nfq_mat_nrows(mat) ((mat)->r)
#define nfq_mat_ncols(mat) ((mat)->c)

void nfq_mat_init(nfq_mat_t mat, slong rows, slong cols, const nfz_ctx_t ctx);

void nfq_mat_clear(nfq_mat_t mat, const nfz_ctx_t ctx);


// void nfq_mat_print(const nfq_mat_t mat, const nfz_ctx_t ctx);

/* Basic assignment **********************************************************/

void nfq_mat_set(nfq_mat_t dest, const nfq_mat_t src, const nfz_ctx_t ctx);

void nfq_mat_zero(nfq_mat_t mat, const nfz_ctx_t ctx);

void nfq_mat_one(nfq_mat_t mat, const nfz_ctx_t ctx);

void nfq_mat_transpose(nfq_mat_t rop, const nfq_mat_t op, const nfz_ctx_t ctx);

/* Addition, scalar multiplication  ******************************************/

void nfq_mat_add(nfq_mat_t mat, const nfq_mat_t mat1, const nfq_mat_t mat2,
		 const nfz_ctx_t ctx);

void nfq_mat_sub(nfq_mat_t mat, const nfq_mat_t mat1, const nfq_mat_t mat2,
		 const nfz_ctx_t ctx);

void nfq_mat_neg(nfq_mat_t rop, const nfq_mat_t op, const nfz_ctx_t ctx);

void nfq_mat_scalar_mul_nfz(nfq_mat_t rop, const nfq_mat_t op, const nfz_t x,
			    const nfz_ctx_t ctx);

void nfq_mat_scalar_div_nfz(nfq_mat_t rop, const nfq_mat_t op, const nfz_t x,
			    const nfz_ctx_t ctx);

/* Basic comparison and properties *******************************************/

int nfq_mat_equal(const nfq_mat_t mat1, const nfq_mat_t mat2,
		  const nfz_ctx_t ctx);

// int nfq_mat_is_integral(const nfq_mat_t mat, const nfz_ctx_t ctx);

int nfq_mat_is_zero(const nfq_mat_t mat, const nfz_ctx_t ctx);

static __inline__ int
nfq_mat_is_empty(const nfq_mat_t mat, const nfz_ctx_t ctx)
{
    return (mat->r == 0) || (mat->c == 0);
}

static __inline__ int
nfq_mat_is_square(const nfq_mat_t mat, const nfz_ctx_t ctx)
{
    return (mat->r == mat->c);
}

/* Integer matrix conversion *************************************************/

int nfq_mat_get_nfz_mat(nfz_mat_t dest, const nfq_mat_t mat, const nfz_ctx_t ctx);

void nfq_mat_get_nfz_mat_entrywise(nfz_mat_t num, nfz_mat_t den,
    const nfq_mat_t mat, const nfz_ctx_t ctx);

void nfq_mat_get_nfz_mat_matwise(nfz_mat_t num, nfz_t den,
    const nfq_mat_t mat, const nfz_ctx_t ctx);

void nfq_mat_get_nfz_mat_rowwise(nfz_mat_t num, nfz * den,
    const nfq_mat_t mat, const nfz_ctx_t ctx);

void nfq_mat_get_nfz_mat_colwise(nfz_mat_t num, nfz * den,
    const nfq_mat_t mat, const nfz_ctx_t ctx);

void nfq_mat_get_nfz_mat_rowwise_2(nfz_mat_t num, nfz_mat_t num2,
        nfz * den, const nfq_mat_t mat, const nfq_mat_t mat2, const nfz_ctx_t ctx);

void nfq_mat_get_nfz_mat_mod_nfz(nfz_mat_t dest, const nfq_mat_t mat,
    const nfz_t mod, const nfz_ctx_t ctx);

void nfq_mat_set_nfz_mat(nfq_mat_t dest, const nfz_mat_t src, const nfz_ctx_t ctx);

void nfq_mat_set_nfz_mat_div_nfz(nfq_mat_t X, const nfz_mat_t Xmod,
    const nfz_t div, const nfz_ctx_t ctx);

/* Matrix multiplication *****************************************************/

void nfq_mat_mul(nfq_mat_t C, const nfq_mat_t A, const nfq_mat_t B,
		 const nfz_ctx_t ctx);

void nfq_mat_mul_nfz_mat(nfq_mat_t C, const nfq_mat_t A,
    const nfz_mat_t B, const nfz_ctx_t ctx);

void nfq_mat_mul_r_nfz_mat(nfq_mat_t C, const nfz_mat_t A,
    const nfq_mat_t B, const nfz_ctx_t ctx);

/* Trace *********************************************************************/

void nfq_mat_trace(nfq_t trace, const nfq_mat_t mat, const nfz_ctx_t ctx);

/* Determinant ***************************************************************/

void nfq_mat_det(nfq_t det, const nfq_mat_t mat, const nfz_ctx_t ctx);

/* Nonsingular solving *******************************************************/

// int nfq_mat_solve_fraction_free(nfq_mat_t X, const nfq_mat_t A,
//    const nfq_mat_t B, const nfz_ctx_t ctx);

// int nfq_mat_solve_dixon(nfq_mat_t X, const nfq_mat_t A, const nfq_mat_t B,
//			const nfz_ctx_t ctx);

/* Inverse *******************************************************************/

int nfq_mat_inv(nfq_mat_t B, const nfq_mat_t A, const nfz_ctx_t ctx);

/* Echelon form **************************************************************/

// int nfq_mat_pivot(slong * perm, nfq_mat_t mat, slong r, slong c, const nfz_ctx_t ctx);

// slong nfq_mat_rref_classical(nfq_mat_t B, const nfq_mat_t A, const nfz_ctx_t ctx);

slong nfq_mat_rref_fraction_free(nfq_mat_t B, const nfq_mat_t A, const nfz_ctx_t ctx);

slong nfq_mat_rref(nfq_mat_t B, const nfq_mat_t A, const nfz_ctx_t ctx);


#ifdef __cplusplus
}
#endif

#endif
