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

#ifndef NFQ_H
#define NFQ_H

#include "ulong_extras.h"
#include "flint.h"
#include "fmpq_poly.h"
#include "nfz.h"

#ifdef __cplusplus
extern "C" {}
#endif

typedef fmpq_poly_struct nfq;

typedef nfq nfq_t[1];

/* Memory managment **************************************************/

static __inline__ void nfq_init(nfq_t f, const nfz_ctx_t ctx)
{
  fmpq_poly_init(f);
};

static __inline__ void nfq_init2(nfq_t f, const nfz_ctx_t ctx)
{
  fmpq_poly_init2(f, nfz_ctx_degree(ctx));
};

static __inline__ void nfq_clear(nfq_t f, const nfz_ctx_t ctx)
{
  fmpq_poly_clear(f);
};

static __inline__ void nfq_reduce(nfq_t f, const nfz_ctx_t ctx)
{
  long deg = fmpq_poly_degree(f);
  if (deg > ctx->deg) return;

  fmpq * coeffs = f->coeffs;

  for (long i = deg - 1; i >= ctx->deg; --i) {
    for (long k = ctx->deg - 1; k >= 0; --k)
      fmpq_submul(coeffs + i - (ctx->deg - k), coeffs + i, ctx->modulus->coeffs + k);
    fmpq_zero(coeffs + i);
  }
};

static __inline__ void nfq_add(nfq_t f, const nfq_t g, const nfq_t h, const nfz_ctx_t ctx)
{
  fmpq_poly_add(f, g, h);
};

// void nfq_add_fmpq(nfq_t f, const nfq_t g, const fmpq_t h, const nfz_ctx_t ctx);

// void nfq_add_ui(nfq_t f, const nfq_t g, ulong x, const nfz_ctx_t ctx);

static __inline__ void nfq_sub(nfq_t f, const nfq_t g, const nfq_t h, const nfz_ctx_t ctx)
{
  fmpq_poly_sub(f, g, h);
};

// void nfq_sub_fmpq(nfq_t f, const nfq_t g, const fmpq_t h, const nfz_ctx_t ctx);

// void nfq_sub_ui(nfq_t f, const nfq_t g, ulong x);

void nfq_mul(nfq_t f, const nfq_t g, const nfq_t h, const nfz_ctx_t ctx);

static __inline__ void nfq_mul_ui(nfq_t f, const nfq_t g, ulong x, const nfz_ctx_t ctx)
{
  fmpq_poly_scalar_mul_ui(f, g, x);
};

static __inline__ void nfq_mul_si(nfq_t f, const nfq_t g, slong x, const nfz_ctx_t ctx)
{
  fmpq_poly_scalar_mul_si(f, g, x);
};

void nfq_addmul(nfq_t f, const nfq_t g, const nfq_t h, const nfz_ctx_t ctx);

// void nfq_addmul_fmpq(nfq_t f, const nfq_t g, ulong x, const nfz_ctx_t ctx);

// void nfq_addmul_ui(nfq_t f, const nfq_t g, ulong x, const nfz_ctx_t ctx);

void nfq_submul(nfq_t f, const nfq_t g, const nfq_t h, const nfz_ctx_t ctx);

// void nfq_submul_fmpq(nfq_t f, const nfq_t g, ulong x, const nfz_ctx_t ctx)

// void nfq_submul_ui(nfq_t f, const nfq_t g, ulong x, const nfz_ctx_t ctx);

void nfq_div(nfq_t f, const nfq_t g, const nfq_t h, const nfz_ctx_t ctx);

// void nfq_divexact_si(nfq_t f, const nfq_t g, slong x, const nfz_ctx_t ctx)

// void nfq_divexact_ui(nfq_t f, const nfq_t g, slong x, const nfz_ctx_t ctx)


#ifdef __cplusplus
}
#endif

#endif
