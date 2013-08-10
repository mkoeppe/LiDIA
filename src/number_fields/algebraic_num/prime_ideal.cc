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
//	Author	: Stefan Neis (SN)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/prime_ideal.h"
#include	"LiDIA/debug.h"
#include	<cctype>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// Constructors & destructor
prime_ideal::prime_ideal ()
	: gen1(0),
	  gen2(),
	  e(0),
	  f(0),
	  valu(bigint(0))
{
	debug_handler_l("prime_ideal", "constructor (void)", 1);
}



prime_ideal::prime_ideal (const bigint & prim, const alg_number & a,
			  lidia_size_t ram, lidia_size_t inert)
	: gen1(prim),
	  gen2(a),
	  e(ram),
	  f(inert),
	  valu(0, a.which_base())
{
	debug_handler_l("prime_ideal", "constructor (const bigint &, "
			"const alg_number &, lidia_size_t, lidia_size_t)", 1);
	if (!is_prime(gen1, 10) || !a.denominator().is_one())
		lidia_error_handler("prime_ideal", "constructor can only be called with\n"
				    "\t a prime number  and\n"
				    "\t an integral algebraic number\n as arguments");

	if (f == 0 || e == 0) {
		module M = alg_ideal(gen1, gen2);
		if (f == 0) {
			f = a.degree() - M.coeff_matrix().get_no_of_columns();
		}
		if (e == 0) {
			M %= gen1;
			module P1;
			module P2 = M;

			do {
				P1 = P2;
				e++;
				multiply(P2, P1, M);
				P2 %= gen1;
			} while (P2.coeff_matrix().get_no_of_columns() <
				 P1.coeff_matrix().get_no_of_columns() ||
				 (P2.coeff_matrix().is_column_zero(0) &&
				  !P1.coeff_matrix().is_column_zero(0)));
		}
	}
}



prime_ideal::prime_ideal (const bigint & prim, const nf_base * O)
	: gen1(prim),
	  gen2(O),
	  e(1),
	  f(O->degree()),
	  valu(0, O)
{
	debug_handler_l("prime_ideal", "constructor (const bigint &, "
			"const nf_base &)", 1);
	if (!is_prime(gen1, 10))
		lidia_error_handler("prime_ideal", "constructor can only be called with\n"
				    "\t a prime number and\n"
				    "\t (optionally) an nf_base/order/number_field");
}



//Computing the help variable for computing valuations(Cohen, p. 199-201)
void
prime_ideal::compute_valu () const
{
	debug_handler("prime_ideal", "in member-function compute_valu()");

	// Set a to the second generator of the prime_ideal
	lidia_size_t n = gen2.degree();
	bigint fact;

	bigmod_matrix LGS(n, n, gen1);

	bigint *tmp = new bigint[n];

	for (register lidia_size_t i = 0; i < n; tmp[i++] = 0) {
		tmp[i] = 1; //tmp = e_i;
		alg_number b (tmp, 1, gen2.which_base()); //b = w_i
		multiply(b, b, gen2);
		LGS.sto_column_vector(b.coeff_vector(), n, i);
	}
	delete[] tmp;

	debug_handler_c("prime_ideal", "in member-function compute_valu()", 2,
			std::cout << "Computing kernel of " << LGS << " mod " << gen1);
	LGS.kernel(LGS, fact);
	debug_handler_c("prime_ideal", "in member-function compute_valu()", 2,
			std::cout << "kernel is " << LGS);

	if (!fact.is_one()) {
		lidia_error_handler("prime_ideal", "compute_valu()::"
				    "internal error 1 while computing kernel");
	}
	valu = alg_number(tmp = LGS.get_column(0), 1, gen2.which_base());
	delete[] tmp;
}



// swap
void
swap (prime_ideal &a, prime_ideal &b) {
	swap(a.gen1, b.gen1);
	swap(a.gen2, b.gen2);
	swap(a.valu, b.valu);

	lidia_size_t tmp;

	tmp = a.e;
	a.e = b.e;
	b.e = tmp;

	tmp = a.f;
	a.f = b.f;
	b.f = tmp;
}



// High-level functions:
long
ord (const prime_ideal & prim, const alg_ideal & MM)
{
	debug_handler_c("prime_ideal", "in ord(prime_ideal & alg_ideal &)", 3,
			std::cout << "Computing valuation of " << MM << " at " << prim);
	//  const bigmod_matrix & Mcoeff = M.coeff_matrix();
	// const bigint & Mmod = Mcoeff.get_modulus();
	long v = 0;
	if (prim.gen2.is_zero()) {
		bigint res, tmp(MM.denominator());
		div_rem(tmp, res, tmp, prim.gen1);
		while (res.is_zero()) {
			v --;
			div_rem(tmp, res, tmp, prim.gen1);
		}
		bigmod_matrix M = MM.coeff_matrix();
		if (v == 0)
			// special algorithm for inert primes
			while (divide(M, M, prim.gen1))
				v++;
		return v;
	}
	if (prim.valu.is_zero()) prim.compute_valu();

	bigint res, tmp2(MM.denominator());

	// Initialize with the valuation obtained from the denominator;
	div_rem(tmp2, res, tmp2, prim.gen1);
	while (res.is_zero()) {
		v -= prim.e;
		div_rem(tmp2, res, tmp2, prim.gen1);
	}
	debug_handler_c("prime_ideal", "in ord(prime_ideal &, alg_ideal &)", 3,
			std::cout << "now v is " << v);

	// check whether prim.gen1 divides all pseudo--diagonal elements.
	remainder(res, MM.coeff_matrix().get_modulus(), prim.gen1);
	if (!res.is_zero()) {
		debug_handler_c("prime_ideal", "in ord(prime_ideal &, alg_ideal &)", 3,
				std::cout << "returning(1) " << v << std::endl);
		return v;
	}

	// Now add valuation of the numerator:
	alg_number b;
	alg_ideal M (MM);

	debug_handler_c("prime_ideal", "in ord(prime_ideal &, alg_ideal &)", 3,
			std::cout << "valu is " << prim.valu << std::endl);

	do {
		debug_handler_c("prime_ideal", "in ord(prime_ideal &, alg_ideal &)", 3,
				std::cout << "At start of loop M is " << M << std::endl);
		multiply(M, M, prim.valu);

		debug_handler_c("prime_ideal", "in ord(prime_ideal &, alg_ideal &)", 3,
				std::cout << "Now(1) M is " << M.base << std::endl);

//     remainder(res, prim.O->base_denominator(), prim.p);
//     if (!res.is_zero()){
//       v++;
//       std::cout << "now(1) v is "<<v<<","<<std::endl;
//       std::cout << "  since O is "<<(*prim.O)<<" i.e. field is ";
//       std::cout << (*(prim.O->which_field())) << " with TRAFO";
//       std::cout << prim.O->base_numerator()<<"/"<<prim.O->base_denominator();
//       std::cout <<std::endl<<std::flush;
//       divide(A, A, prim.p);
//       std::cout << "Now(2) A is "<<A<<std::endl;
//       continue;
//     }
		if (!divide(M.base, M.base, prim.gen1)) {
			debug_handler_c("prime_ideal", "in ord(prime_ideal &, alg_ideal &)", 3,
					std::cout << "returning(3) " << v << std::endl);
			return v;
		}
		v++;
		debug_handler_c("prime_ideal", "in ord(prime_ideal &, alg_ideal &)", 3,
				std::cout << "Now(2) v is " << v << std::endl;
				std::cout << "and    M is " << M << std::endl);
	} while (true);

}



long
ord (const prime_ideal & prim, const alg_number & a)
{
	return ord(prim, a*order(prim.gen2.which_base()));
}



void
power (alg_ideal & c, const prime_ideal & a, const bigint & b)
{
	debug_handler("prime_ideal", "in function power(alg_ideal &, "
		      "const prime_ideal &, const bigint &)");

	const nf_base * O = a.gen2.which_base();
	if (b.is_negative())
		power(c, inverse(alg_ideal(alg_number(a.gen1, a.gen2.which_base()),
					   a.gen2)), -b);
	else if (b.is_zero())
		c = order(O);
	else {
		bigint gen1;
		bigint expo;
		div_rem(expo, gen1, b, a.e);
		if (!gen1.is_zero()) ++expo;
		power(gen1, a.gen1, expo);
		if (!a.gen2.is_zero()) {
			alg_number gen2;
			power(gen2, a.gen2, b);
			c = alg_ideal(gen1, gen2);
		}
		else
			c = alg_ideal(gen1, alg_number(0, O));
	}
}


// Input/Output:
std::ostream &
operator << (std::ostream & s, const prime_ideal & a)
{
	s << "<";
	s << a.gen1;
	if (!a.gen2.is_zero()) {
		s << ", ";
		s << a.gen2;
	}
	s << ">" << std::flush;
	return s;
}



std::istream &
operator >> (std::istream & s, prime_ideal & p)
{
	bigint a;
	alg_number b;
	char c;

	do {
		s.get(c);
	} while (isspace(c));
	if (c != '<')
		lidia_error_handler("prime_ideal", "operator >>:: A prime ideal must be"
				    "of form `<prime, alg_number>'");
	s >> a;

	do {
		s.get(c);
	} while (isspace(c));
	if (c != ',')
		lidia_error_handler("prime_ideal", "operator >>:: A prime ideal must be"
				    "of form `<prime, alg_number>'");

	s >> b;
	do {
		s.get(c);
	} while (isspace(c));
	if (c != '>')
		lidia_error_handler("prime_ideal", "operator >>:: A prime ideal must be"
				    "of form `<prime, alg_number>'");

	p = prime_ideal(a, b);
	return s;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
