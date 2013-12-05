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

#ifdef __cplusplus
extern "C" {
#endif

typedef struct
{
  fmpz_poly_t modulus;

  fmpz_mat_t ev_mat;  
  fmpz_mat_t int_mat;
  fmpz_t int_den;

  char * var;
}
nfz_ctx_struct;

typedef nfz_ctx_struct nfz_ctx_t[1];

/* Memory managment  *********************************************************/

void nfz_ctx_init(nfz_ctx_t ctx, const fmpz_poly_t, const char *var);
void nfz_ctx_clear(nfz_ctx_t ctx);

statix __inline__
slong nfz_ctx_degree(const nfz_ctx_t ctx)
{
  return fmpz_poly_degree(ctx->modulus);
}

static __inline__
int nfz_ctx_fprint(FILE * file, const nfz_ctx_t ctx)
{
    int r;

    r = flint_fprintf(file, "modulus =");
    if (r <= 0)
	return r;

    r = fmpz_mod_poly_fprint_pretty(file, ctx->modulus, ctx->var);
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

#ifdef __cplusplus
}
#endif

#endif

