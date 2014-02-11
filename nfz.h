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

#ifndef NFZ_H
#define NFZ_H

#include "ulong_extras.h"
#include "flint.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"
#include "nf_nmod.h"
#include "nf_nmod_mat.h"

#ifdef __cplusplus
extern "C" {
#endif

/*****************************************************************************/
/* Number field context ******************************************************/
/*****************************************************************************/

typedef struct
{
  fmpz_poly_t modulus;
  slong deg;

  // evaluation and interpolation matrices act on columns
  fmpz_mat_t evl_mat;  
  fmpz_mat_t intrpl_mat;
  fmpz_t intrpl_den;

  char * var;
}
nfz_ctx_struct;

typedef nfz_ctx_struct nfz_ctx_t[1];

/* Memory managment ***********************************************************/

void nfz_ctx_init(nfz_ctx_t ctx, const fmpz_poly_t, const char *var);

void nfz_ctx_clear(nfz_ctx_t ctx);

/* Parameters *****************************************************************/

static __inline__
slong nfz_ctx_degree(const nfz_ctx_t ctx)
{
  return fmpz_poly_degree(ctx->modulus);
}

/* Output *********************************************************************/

static __inline__
int nfz_ctx_fprint(FILE * file, const nfz_ctx_t ctx)
{
    int r;

    r = flint_fprintf(file, "modulus =");
    if (r <= 0)
	return r;

    r = fmpz_poly_fprint_pretty(file, ctx->modulus, ctx->var);
    if (r <= 0)
	return r;

    r = flint_fprintf(file, "\n");

    return r;
}

static __inline__
void nfz_ctx_print(const nfz_ctx_t ctx)
{
    nfz_ctx_fprint(stdout, ctx);
}

/* Recuction mod number field modulus *****************************************/

void _nfz_reduction_mat(fmpz_mat_t mat, const fmpz_poly_t modulus, ulong deg_bd);

void _nfz_reduction_coeff_bound(fmpz_t bound, const fmpz_poly_t modulus, ulong deg_bd);

/* Reduction mod p ************************************************************/

void nfz_ctx_get_nmod_ctx(nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx, ulong mod);

/******************************************************************************/
/* Number field elements ******************************************************/
/******************************************************************************/

/* We implement elements of number fields (defined by monic
 * polynomials, say p) whose representation as a polynomial mod p has
 * integral coefficients. (This is, of cause, not the maximal order,
 * and depends on p).
 */

typedef fmpz_poly_struct nfz;

typedef nfz nfz_t[1];

/* Memory managment ***********************************************************/

static __inline__
void nfz_init(nfz_t f, const nfz_ctx_t ctx)
{
  fmpz_poly_init(f);
};

static __inline__
void nfz_init2(nfz_t f, const nfz_ctx_t ctx)
{
    fmpz_poly_init2(f, nfz_ctx_degree(ctx));
};

static __inline__
void nfz_clear(nfz_t f, const nfz_ctx_t ctx)
{
    fmpz_poly_clear(f);
};

/* Reduction of representatives ***********************************************/

static __inline__
void nfz_reduce(nfz_t f, const nfz_ctx_t ctx)
{
  long deg = fmpz_poly_degree(f);
  if (deg > ctx->deg) return;

  fmpz * coeffs = f->coeffs;

  for (long i = deg - 1; i >= ctx->deg; --i) {
    for (long k = ctx->deg - 1; k >= 0; --k)
      fmpz_submul(coeffs + i - (ctx->deg - k), coeffs + i, ctx->modulus->coeffs + k);
    fmpz_zero(coeffs + i);
  }
  _fmpz_poly_set_length(f, ctx->deg);
  _fmpz_poly_normalise(f);
};

/* Evaluation and interpolation ***********************************************/

void _nfz_eval(fmpz * evl, const fmpz * f, long length, const nfz_ctx_t ctx);

void _nfz_interpolate(fmpz * f, long * length, const fmpz * evl,
		      const nfz_ctx_t ctx);

void _nfz_eval_pptr(fmpz ** evl, const fmpz ** f, long length, const nfz_ctx_t ctx);

void _nfz_interpolate_pptr(fmpz ** f, long * length, const fmpz ** evl,
			   const nfz_ctx_t ctx);

/* Conversions ****************************************************************/

static __inline__
void nfz_get_nmod(nf_nmod_t f, const nfz_t g,
		  const nf_nmod_ctx_t ctx_nmod,
		  const nfz_ctx_t ctx)
{
  fmpz_poly_get_nmod_poly(f, g);
};

static __inline__
void nfz_set_nmod(nfz_t f, const nf_nmod_t g,
		  const nf_nmod_ctx_t ctx_nmod,
		  const nfz_ctx_t ctx)
{
  fmpz_poly_set_nmod_poly(f, g);
};

/* Reduction mod p ************************************************************/

static __inline__
void nfz_scalar_mod_fmpz(nfz_t f, const nfz_t g, const fmpz_t p,
			 const nfz_ctx_t ctx)
{
  fmpz_poly_scalar_mod_fmpz(f, g, p);
};

static __inline__
void nfz_scalar_smod_fmpz(nfz_t f, const nfz_t g, const fmpz_t p,
			  const nfz_ctx_t ctx)
{
  fmpz_poly_scalar_smod_fmpz(f, g, p);
};

/* Assignment and basic manipulation ******************************************/

static __inline__
void nfz_set(nfz_t f, const nfz_t g, const nfz_ctx_t ctx)
{
  fmpz_poly_set(f, g);
};

static __inline__
void nfz_swap(nfz_t f, nfz_t g, const nfz_ctx_t ctx)
{
  fmpz_poly_swap(f, g);
};

static __inline__
void nfz_zero(nfz_t f, const nfz_ctx_t ctx)
{
  fmpz_poly_zero(f);
};

static __inline__
void nfz_one(nfz_t f, const nfz_ctx_t ctx)
{
  fmpz_poly_one(f);
};

static __inline__
void nfz_neg(nfz_t f, const nfz_t g, const nfz_ctx_t ctx)
{
  fmpz_poly_neg(f, g);
};

/* Comparison *****************************************************************/

static __inline__
int nfz_equal(const nfz_t f, const nfz_t g, const nfz_ctx_t ctx)
{
  return fmpz_poly_equal(f, g);
};

static __inline__
int nfz_is_zero(const nfz_t f, const nfz_ctx_t ctx)
{
  return fmpz_poly_is_zero(f);
};

/* Addition *******************************************************************/

static __inline__
void nfz_add(nfz_t f, const nfz_t g, const nfz_t h, const nfz_ctx_t ctx)
{
  fmpz_poly_add(f, g, h);
};

// void nfz_add_fmpz(nfz_t f, const nfz_t g, const fmpz_t h, const nfz_ctx_t ctx);

// void nfz_add_si(nfz_t f, const nfz_t g, slong x, const nfz_ctx_t ctx);

// void nfz_add_ui(nfz_t f, const nfz_t g, ulong x, const nfz_ctx_t ctx);

static __inline__
void nfz_sub(nfz_t f, const nfz_t g, const nfz_t h, const nfz_ctx_t ctx)
{
  fmpz_poly_sub(f, g, h);
};

// void nfz_sub_fmpz(nfz_t f, const nfz_t g, const fmpz_t h, const nfz_ctx_t ctx);

// void nfz_sub_si(nfz_t f, const nfz_t g, slong x, const nfz_ctx_t ctx);

// void nfz_sub_ui(nfz_t f, const nfz_t g, ulong x);

/* Multiplication and division ************************************************/

void nfz_mul(nfz_t f, const nfz_t g, const nfz_t h, const nfz_ctx_t ctx);

static __inline__
void nfz_mul_fmpz(nfz_t f, const nfz_t g, fmpz_t x, const nfz_ctx_t ctx)
{
  fmpz_poly_scalar_mul_fmpz(f, g, x);
};

static __inline__
void nfz_mul_si(nfz_t f, const nfz_t g, slong x, const nfz_ctx_t ctx)
{
  fmpz_poly_scalar_mul_si(f, g, x);
};

static __inline__
void nfz_mul_ui(nfz_t f, const nfz_t g, ulong x, const nfz_ctx_t ctx)
{
  fmpz_poly_scalar_mul_ui(f, g, x);
};

// void nfz_divexact(nfz_t f, const nfz_t g, const nfz_t h, const nfz_ctx_t ctx);

static __inline__
void nfz_divexact_fmpz(nfz_t f, const nfz_t g, const fmpz_t x, const nfz_ctx_t ctx)
{
  fmpz_poly_scalar_divexact_fmpz(f, g, x);
};

static __inline__
void nfz_divexact_si(nfz_t f, const nfz_t g, ulong x, const nfz_ctx_t ctx)
{
  fmpz_poly_scalar_divexact_si(f, g, x);
};

static __inline__
void nfz_divexact_ui(nfz_t f, const nfz_t g, slong x, const nfz_ctx_t ctx)
{
  fmpz_poly_scalar_divexact_ui(f, g, x);
};

void nfz_addmul(nfz_t f, const nfz_t g, const nfz_t h, const nfz_ctx_t ctx);

static __inline__
void nfz_addmul_fmpz(nfz_t f, const nfz_t g, const fmpz_t h, const nfz_ctx_t ctx)
{
  fmpz_poly_scalar_addmul_fmpz(f, g, h);
};

// void nfz_addmul_si(nfz_t f, const nfz_t g, slong x, const nfz_ctx_t ctx);

// void nfz_addmul_ui(nfz_t f, const nfz_t g, ulong x, const nfz_ctx_t ctx);

void nfz_submul(nfz_t f, const nfz_t g, const nfz_t h, const nfz_ctx_t ctx);

static __inline__
void nfz_submul_fmpz(nfz_t f, const nfz_t g, const fmpz_t h, const nfz_ctx_t ctx)
{
  fmpz_poly_scalar_submul_fmpz(f, g, h);
};

// void nfz_submul_si(nfz_t f, const nfz_t g, slong x, const nfz_ctx_t ctx);

// void nfz_submul_ui(nfz_t f, const nfz_t g, ulong x, const nfz_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif



