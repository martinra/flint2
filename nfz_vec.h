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

nfz * _nfz_vec_init(slong len, const nfz_ctx_t ctx);

void _nfz_vec_clear(nfz * vec, slong len, const nfz_ctx_t ctx);

/*  Conversions  *************************************************************/

void _nfz_vec_set_nmod_vec(nfz * res, 
			   const nf_nmod * src, slong len,
			   const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx);

void _nfz_vec_get_nmod_vec(nf_nmod * res, 
			   const nfz * poly, slong len,
			   const nf_nmod_ctx_t ctx_nmod, const nfz_ctx_t ctx);

/*  Reduction mod p **********************************************************/

void _nfz_vec_scalar_mod_fmpz(nfz * res, const nfz * vec, slong len,
			      const fmpz_t p,
			      const nfz_ctx_t ctx);

void _nfz_vec_scalar_smod_fmpz(nfz * res, const nfz * vec, slong len,
			       const fmpz_t p,
			       const nfz_ctx_t ctx);

/*  Assignment and basic manipulation  ***************************************/

void _nfz_vec_set(nfz * vec1, const nfz * vec2, slong len2, const nfz_ctx_t ctx);

void _nfz_vec_swap(nfz * vec1, nfz * vec2, slong len2, const nfz_ctx_t ctx);

void _nfz_vec_zero(nfz * vec, slong len, const nfz_ctx_t ctx);

void _nfz_vec_neg(nfz * vec1, const nfz * vec2, slong len2, const nfz_ctx_t ctx);

/*  Comparison  **************************************************************/

int _nfz_vec_equal(const nfz * vec1, const nfz * vec2, slong len,
		   const nfz_ctx_t ctx);

int _nfz_vec_is_zero(const nfz * vec, slong len, const nfz_ctx_t ctx);

/*  Addition  ****************************************************************/

void _nfz_vec_add(nfz * res, const nfz * vec1, 
		  const nfz * vec2, slong len2,
		  const nfz_ctx_t ctx);

void _nfz_vec_sub(nfz * res, const nfz * vec1, 
		  const nfz * vec2, slong len2,
		  const nfz_ctx_t ctx);

/*  Scalar multiplication and division  **************************************/

void _nfz_vec_scalar_mul_nfz(nfz * vec1, 
			     const nfz * vec2, slong len2, const nfz_t x,
			     const nfz_ctx_t ctx);

void _nfz_vec_scalar_mul_fmpz(nfz * vec1, 
			      const nfz * vec2, slong len2, const fmpz_t x,
			      const nfz_ctx_t ctx);

void _nfz_vec_scalar_mul_si(nfz * vec1, 
			    const nfz * vec2, slong len2, slong c,
			    const nfz_ctx_t ctx);

void _nfz_vec_scalar_mul_ui(nfz * vec1, 
			    const nfz * vec2, slong len2, ulong c,
			    const nfz_ctx_t ctx);


// void _nfz_vec_scalar_mul_2exp(fmpz * vec1, 
// 			      const fmpz * vec2, slong len2, ulong exp,
// 			      const nfz_ctx_t ctx);

void _nfz_vec_scalar_divexact_fmpz(nfz * vec1, const nfz * vec2, 
				   slong len2, const fmpz_t x,
				   const nfz_ctx_t ctx);

void _nfz_vec_scalar_divexact_si(nfz * vec1, 
				 const nfz * vec2, slong len2, slong c,
				 const nfz_ctx_t ctx);

void _nfz_vec_scalar_divexact_ui(nfz * vec1, 
				 const nfz * vec2, slong len2, ulong c,
				 const nfz_ctx_t ctx);

void _nfz_vec_scalar_addmul_nfz(nfz * poly1, const nfz * poly2, 
				slong len2, const nfz_t x,
				const nfz_ctx_t ctx);

void _nfz_vec_scalar_addmul_fmpz(nfz * poly1, const nfz * poly2, 
				 slong len2, const fmpz_t x,
				 const nfz_ctx_t ctx);

void _nfz_vec_scalar_addmul_si(nfz * vec1, 
			       const nfz * vec2, slong len2, slong c,
			       const nfz_ctx_t ctx);

void _nfz_vec_scalar_addmul_ui(nfz * vec1, 
			       const nfz * vec2, slong len2, ulong c,
			       const nfz_ctx_t ctx);

void _fmpz_vec_scalar_addmul_si_2exp(fmpz * vec1, const fmpz * vec2, 
                                               slong len2, slong c, ulong exp);


void _nfz_vec_scalar_submul_nfz(nfz * poly1, const nfz * poly2, 
				slong len2, const nfz_t x,
				const nfz_ctx_t ctx);

void _nfz_vec_scalar_submul_fmpz(nfz * poly1, const nfz * poly2, 
				 slong len2, const fmpz_t x,
				 const nfz_ctx_t ctx);

void _nfz_vec_scalar_submul_si(nfz * vec1, 
			       const nfz * vec2, slong len2, slong c,
			       const nfz_ctx_t ctx);

void _nfz_vec_scalar_submul_ui(nfz * vec1, 
			       const nfz * vec2, slong len2, ulong c,
			       const nfz_ctx_t ctx);

/*  Vector sum and product  **************************************************/

void _nfz_vec_sum(nfz_t res, const nfz * vec, slong len,
		  const nfz_ctx_t ctx);

void _nfz_vec_prod(nfz_t res, const nfz * vec, slong len,
		   const nfz_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif

