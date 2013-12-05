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

#include "ulong_extras.h"
#include "nfz_mat.h"

int next_prime(const int);

slong
nfz_mat_rref_mod(nfz_mat_t B, fmpz_t den, const nfz_mat_t A, const nfz_ctx_t ctx)
{
  return _nfz_mat_rref_mod_prime_generator(nfz_mat_t B, fmpz_t den, const nfz_mat_t A, const nfz_ctx_t ctx, &next_prime)
}

// todo: should we declare this as static?
int next_prime(const ulong n)
{
  return n_nextprime(n, 1);
}









