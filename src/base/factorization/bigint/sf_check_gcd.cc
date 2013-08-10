//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//	$Id$
//
//	Author	: Emre Binisik (EB), Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/single_factor.h"
#include	"LiDIA/factorization.h"
#include	"LiDIA/bigmod.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bool single_factor< bigint >::check_gcd(factorization< bigint > & f,
					 const bigmod & v, const bigint & N)
{
	bigint d;

	d = abs(gcd(mantissa(v), N));

	if (! d.is_one() && d != N)        // proper factor found
	{
		if (info)
			std::cout << "\nfactor: " << d << std::flush;

		single_factor< bigint > a(d);
		if (is_prime(d, 8))
			a.set_prime_flag(prime);
		else
			a.set_prime_flag(not_prime);
		f.append(a);

		divide(a.rep, N, d);

		if (is_prime(a.rep, 8))
			a.set_prime_flag(prime);
		else
			a.set_prime_flag(not_prime);
		f.append(a);
		rep.assign_one();
		return true;
	}
	return false;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
