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
#include "fmpz.h"
#include "nfz.h"

#ifdef __cplusplus
extern "C" {}
#endif
typedef struct
{
  nfz num;
  fmpz den;
} nfq;

typedef nfq nfq_t[1];

#define nfq_numref(__x) (&(__x)->num)
#define nfq_denref(__y) (&(__y)->den)


#ifdef __cplusplus
}
#endif

#endif
