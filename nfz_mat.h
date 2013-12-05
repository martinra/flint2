/* nf.h ---  */

/* Copyright (C) 2013 Martin Raum */

/* Author: Martin Raum */

/* This program is free software; you can redistribute it and/or */
/* modify it under the terms of the GNU General Public License */
/* as published by the Free Software Foundation; either version 3 */
/* of the License, or (at your option) any later version. */

/* This program is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with this program. If not, see <http://www.gnu.org/licenses/>. */

#ifndef NFZ_MAT_H
#define NFZ_MAT_H

#include "ulong_extras.h"
#include "flint.h"
#include "fmpq_mat.h"
#include "nfz.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
    mp_limb_t * entries;
    slong r;
    slong c;
    mp_limb_t ** rows;
    mp_limb_t *** poly_coefficients;
    fmpz_poly_t mod;
}
nfz_mat_struct;

typedef nfz_mat_struct nfz_mat_t[1];

/* Memory managment  *********************************************************/

#define nfz_mat_entry(mat,e,i,j) ((mat)->poly_coefficients[(e)][(i)][(j)])
#define nfz_mat_nrows(mat) ((mat)->r)
#define nfz_mat_ncols(mat) ((mat)->c)

void nfz_mat_init(fmpz_mat_t mat, slong rows, slong cols, const nfz_ctx_t ctx);
void nfz_mat_init_set(fmpz_mat_t mat, const fmpz_mat_t src, const nfz_ctx_t ctx);
void nfz_mat_set(nfz_mat_t mat1, const nfz_mat_t mat2, const nfz_ctx_t ctx);
void nfz_mat_clear(nfz_mat_t mat, const nfz_ctx_t ctx);

slong _nfz_mat_rref_mod_prime_generator(nfz_mat_t B, fmpz_t den, const nfz_mat_t A, const nfz_ctx_t ctx, int (*next_prime)(int));

#ifdef __cplusplus
}
#endif

#endif

