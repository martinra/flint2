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

#include "nf_nmod_mat.h"

void
_nf_nmod_mat_clear_components(nmod_mat_t * fp_comps, fq_nmod_mat_t * fq_comps, fq_nmod_ctx_t * fq_ctxs, const nf_nmod_ctx_t ctx)
{
  int i;

  for (i = 0; i < ctx->nfp; ++i)
    nmod_mat_clear(fp_comps[i]);

  for (i = 0; i < ctx->nfq; ++i)
    {
      fq_nmod_mat_clear(fq_comps[i], fq_ctxs[i]);
      fq_nmod_ctx_clear(fq_ctxs[i]);
    }

  flint_free(fp_comps);
  flint_free(fq_comps);
  flint_free(fq_ctxs);
}
