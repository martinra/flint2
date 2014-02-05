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

#ifndef NFZ_VEC_H
#define NFZ_VEC_H

#include "nfz.h"
#include "nf_nmod.h"

#ifdef __cplusplus
 extern "C" {
#endif

/*  Memory management  *******************************************************/

static __inline__
nfz * _nfz_vec_init(slong len, const nfz_ctx_t ctx)
{
  nfz * res = (nfz *) flint_malloc(len * sizeof(nfz));
  for (long i = 0; i < len; ++i)
    nfz_init(res + i, ctx);
  return res;
};

static __inline__
nfz * _nfz_vec_init_set(const nfz * vec, slong len, const nfz_ctx_t ctx)
{
  nfz * res = _nfz_vec_init(len, ctx);
  for (long i = 0; i < len; ++i)
    nfz_set(res + i, vec + i, ctx);
  return res;
};

static __inline__
void _nfz_vec_clear(nfz * vec, slong len, const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; i++)
    nfz_clear(vec + i, ctx);
  flint_free(vec);
};

/*  Conversions  *************************************************************/

static __inline__
void _nfz_vec_set_nmod_vec(nfz * res, 
			   const nf_nmod * src, slong len,
			   const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_set_nmod(res + i, src + i, ctx_nmod, ctx);
};

static __inline__
void _nfz_vec_get_nmod_vec(nf_nmod * res, 
			   const nfz * src, slong len,
			   const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_get_nmod(res + i, src + i, ctx_nmod, ctx);
};

/*  Reduction mod p **********************************************************/

static __inline__
void _nfz_vec_scalar_mod_fmpz(nfz * res, const nfz * vec, slong len,
			      const fmpz_t p,
			      const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_scalar_mod_fmpz(res + i, vec + i, p, ctx);
};

static __inline__
void _nfz_vec_scalar_smod_fmpz(nfz * res, const nfz * vec, slong len,
			       const fmpz_t p,
			       const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_scalar_smod_fmpz(res + i, vec + i, p, ctx);
};

/*  Assignment and basic manipulation  ***************************************/

static __inline__
void _nfz_vec_set(nfz * vec1, const nfz * vec2, slong len, const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_set(vec1 + i, vec2 + i, ctx);
};

static __inline__
void _nfz_vec_swap(nfz * vec1, nfz * vec2, slong len, const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_swap(vec1 + i, vec2 + i, ctx);
};

static __inline__
void _nfz_vec_zero(nfz * vec, slong len, const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_zero(vec + i, ctx);
};

static __inline__
void _nfz_vec_one(nfz * vec, slong len, const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_one(vec + i, ctx);
};

static __inline__
void _nfz_vec_neg(nfz * vec1, const nfz * vec2, slong len, const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_neg(vec1 + i, vec2 + i, ctx);
};

/*  Comparison  **************************************************************/

static __inline__
int _nfz_vec_equal(const nfz * vec1, const nfz * vec2, slong len,
		   const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    if (!nfz_equal(vec1 + i, vec2 + i, ctx))
      return 0;

  return 1;
};

static __inline__
int _nfz_vec_is_zero(const nfz * vec, slong len, const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    if (!nfz_is_zero(vec + i, ctx))
      return 0;

  return 1;
};

/*  Addition  ****************************************************************/

static __inline__
void _nfz_vec_add(nfz * res, const nfz * vec1, 
		  const nfz * vec2, slong len,
		  const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_add(res + i, vec1 + i, vec2 + i, ctx);
};

// void _nfz_vec_add_fmpz(nfz * res, const nfz * g, const fmpz * h, slong len, const nfz_ctx_t ctx);

static __inline__
void _nfz_vec_sub(nfz * res, const nfz * vec1, 
		  const nfz * vec2, slong len,
		  const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_sub(res + i, vec1 + i, vec2 + i, ctx);
};

// void _nfz_vec_sub_fmpz(nfz * res, const nfz * g, const fmpz * h, slong len, const nfz_ctx_t ctx);

/*  Scalar multiplication and division  **************************************/

static __inline__
void _nfz_vec_scalar_mul(nfz * res, 
			 const nfz * vec, slong len, const nfz_t c,
			 const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_mul(res + i, vec + i, c, ctx);
};

static __inline__
void _nfz_vec_scalar_mul_fmpz(nfz * res, 
			      const nfz * vec, slong len, const fmpz_t c,
			      const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_mul_fmpz(res + i, vec + i, c, ctx);
};

static __inline__
void _nfz_vec_scalar_mul_si(nfz * res, 
			    const nfz * vec, slong len, slong c,
			    const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_mul_si(res + i, vec + i, c, ctx);
};

static __inline__
void _nfz_vec_scalar_mul_ui(nfz * res, 
			    const nfz * vec, slong len, ulong c,
			    const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_mul_ui(res + i, vec + i, c, ctx);
};

static __inline__
void _nfz_vec_scalar_divexact(nfz * res, const nfz * vec, 
				   slong len, const nfz_t c,
				   const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_divexact(res + i, vec + i, c, ctx);
};

static __inline__
void _nfz_vec_scalar_divexact_fmpz(nfz * res, const nfz * vec, 
				   slong len, const fmpz_t c,
				   const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_divexact_fmpz(res + i, vec + i, c, ctx);
};

static __inline__
void _nfz_vec_scalar_divexact_si(nfz * res, 
				 const nfz * vec, slong len, slong c,
				 const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_divexact_si(res + i, vec + i, c, ctx);
};

static __inline__
void _nfz_vec_scalar_divexact_ui(nfz * res, 
				 const nfz * vec, slong len, ulong c,
				 const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_divexact_ui(res + i, vec + i, c, ctx);
};

static __inline__
void _nfz_vec_scalar_addmul_nfz(nfz * res,
				const nfz * vec, slong len, const nfz_t c,
				const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_addmul(res + i, vec + i, c, ctx);
};

static __inline__
void _nfz_vec_scalar_addmul_fmpz(nfz * res,
				 const nfz * vec, slong len, const fmpz_t c,
				 const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_addmul_fmpz(res + i, vec + i, c, ctx);
};


/* static __inline__ */
/* void _nfz_vec_scalar_addmul_si(nfz * res, */
/* 			       const nfz * vec, slong len, slong c, */
/* 			       const nfz_ctx_t ctx) */
/* { */
/*   for (long i = 0; i < len; ++i) */
/*     nfz_addmul_si(res + i, vec + i, c, ctx); */
/* }; */


/* static __inline__ */
/* void _nfz_vec_scalar_addmul_ui(nfz * res, */
/* 			       const nfz * vec, slong len, ulong c, */
/* 			       const nfz_ctx_t ctx) */
/* { */
/*   for (long i = 0; i < len; ++i) */
/*     nfz_addmul_ui(res + i, vec + i, c, ctx); */
/* }; */


static __inline__
void _nfz_vec_scalar_submul_nfz(nfz * res,
				const nfz * vec, slong len, const nfz_t c,
				const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_submul(res + i, vec + i, c, ctx);
};

static __inline__
void _nfz_vec_scalar_submul_fmpz(nfz * res,
				 const nfz * vec, slong len, const fmpz_t c,
				 const nfz_ctx_t ctx)
{
  for (long i = 0; i < len; ++i)
    nfz_submul_fmpz(res + i, vec + i, c, ctx);
};


/* static __inline__ */
/* void _nfz_vec_scalar_submul_si(nfz * res, */
/* 			       const nfz * vec, slong len, slong c, */
/* 			       const nfz_ctx_t ctx) */
/* { */
/*   for (long i = 0; i < len; ++i) */
/*     nfz_submul_si(res + i, vec + i, c, ctx); */
/* }; */


/* static __inline__ */
/* void _nfz_vec_scalar_submul_ui(nfz * res, */
/* 			       const nfz * vec, slong len, ulong c, */
/* 			       const nfz_ctx_t ctx) */
/* { */
/*   for (long i = 0; i < len; ++i) */
/*     nfz_submul_ui(res + i, vec + i, c, ctx); */
/* }; */

/*  Vector sum and product  **************************************************/

static __inline__
void _nfz_vec_sum(nfz_t res, const nfz * vec, slong len,
		  const nfz_ctx_t ctx)
{
  nfz_zero(res, ctx);
  for (long i = 0; i < len; ++i)
    nfz_add(res, res, vec + i, ctx);
};

static __inline__
void _nfz_vec_prod(nfz_t res, const nfz * vec, slong len,
		   const nfz_ctx_t ctx)
{
  nfz_one(res, ctx);
  for (long i = 0; i < len; ++i)
    nfz_mul(res, res, vec + i, ctx);
};

#ifdef __cplusplus
}
#endif

#endif

