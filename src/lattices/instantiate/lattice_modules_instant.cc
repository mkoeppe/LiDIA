/*
 *
 *	This file is part of:
 *
 *		LiDIA --- a library for computational number theory
 *
 *	Copyright (c) 1994--2002 the LiDIA Group.  All rights reserved.
 *
 *	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
 *
 *	$Id$
 *
 */

/*
 * g++-2.95.2 occasionally does not instantiate inline specializations
 * of template class member functions when compiling with
 * -fno-implicit-templates (though that should not affect inline
 * definitions).  To work around this, we control the instantiations
 * ourselves, using the pragmas `interface' and `implementation',
 * which have the additional advantage of avoiding the duplication of
 * instantiations.
 */
#if __GNUC__ && __GNUC__ < 3
#pragma implementation "lattice_modules.cc"
#include "LiDIA/lattices/lattice_modules.h"
#endif

