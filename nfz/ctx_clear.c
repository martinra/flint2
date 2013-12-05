/* ctx_clear.c ---  */

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

#include "nfz.h"

void
nfz_ctx_clear(nfz_ctx_t ctx)
{
  fmpz_mat_clear(ctx->ev_mat);
  fmpz_mat_clear(ctx->int_mat);
  fmpz_clear(ctx->int_den);
  flint_free(ctx->var);
}

