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
//	Author	: Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/gf_polynomial.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



gf_poly_argument::gf_poly_argument() :
	vec(0), len(-1)
{
	debug_handler("gf_poly_modulus", "gf_poly_modulus()");
}



gf_poly_argument::gf_poly_argument(const gf_poly_argument &x)
{
	debug_handler("gf_poly_modulus", "gf_poly_modulus(gf_poly_modulus&)");
	if (x.len < 1) {
		delete[] vec;
		len = -1;
	}
	else
		lidia_error_handler("gf_poly_modulus", "gf_poly_modulus(gf_poly_modulus&)::copy constructor not implemented");
}



gf_poly_argument::~gf_poly_argument()
{
	debug_handler("gf_poly_modulus", "~gf_poly_modulus()");
	delete[] vec;
}



void gf_poly_argument::
build(const gf_polynomial &h, const gf_poly_modulus &F, lidia_size_t m)
	//computes and stores h^0, h, h^2, ..., h^m mod f
{
	debug_handler("gf_poly_modulus", "build(...)");
	if (m <= 0)
		lidia_error_handler("gf_poly_modulus", "build::bad input");

	if (len != m+1) {
		delete[] vec;
		vec = new gf_polynomial[m+1];
		len = m+1;
	}

	vec[0].assign_one(h.get_field());
	vec[1].assign(h);
	for (lidia_size_t i = 2; i <= m; i++)
		multiply(vec[i], vec[i-1], h, F);

	debug_handler_c("gf_poly_modulus", "build", 8,
			std::cout << "gf_poly_argument: table\n";
			for (lidia_size_t j = 0; j < len; j++)
			std::cout << vec[j] << std::endl;);
}



void gf_poly_argument::
compose(gf_polynomial &x, const gf_polynomial &g,
	const gf_poly_modulus &F) const
	//x = g(h) mod F
{
	debug_handler("gf_poly_modulus", "compose(...)");
	if (g.degree() <= 0) {
		x.assign(g);
		return;
	}
	lidia_size_t m = len - 1;
	lidia_size_t l = ((g.degree()+m)/m) - 1;

	gf_polynomial t, s;

	inner_prod(t, g, l*m, l*m + m - 1);

	for (lidia_size_t i = l-1; i >= 0; i--) {
	//####### gf_poly_multiplier ???
		inner_prod(s, g, i*m, i*m + m - 1);
		multiply(t, t, vec[m], F); // t = t * h^m  
		add(t, t, s);
	}

	x.assign(t);
}



void gf_poly_argument::
inner_prod(gf_polynomial &x, const gf_polynomial &g,
	   lidia_size_t lo, lidia_size_t hi) const
{
	debug_handler("gf_poly_modulus", "inner_prod(...)");

	gf_polynomial t(vec[0].get_field()), s;

	hi = comparator< lidia_size_t >::min(hi, g.degree());
	for (lidia_size_t i = lo; i <= hi; i++) {
		multiply(s, vec[i-lo], g[i]);
		add(t, t, s);
	}
	x.assign(t);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
