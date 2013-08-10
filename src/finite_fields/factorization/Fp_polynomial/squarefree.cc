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
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/factorization.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// assumes f(x) = g(x^p), p = f.modulus();
static void
sq_fr_special(factorization< Fp_polynomial > & u, const Fp_polynomial& f)
{
	lidia_size_t p;
	f.modulus().sizetify(p);
	lidia_size_t d = f.degree()/p;

	Fp_polynomial t;
	t.set_modulus(f);
	t.set_max_degree(d);
	lidia_size_t i;
	for (i = 0; i <= d; i++)
		t.set_coefficient(f[i*p], i);

	square_free_decomp(u, t);
	u.power(p);
}



void square_free_decomp(factorization< Fp_polynomial > & u, const Fp_polynomial& f)
{
	// Performs square-free decomposition.
	// ASSUMPTION: deg(f) > 0, f must be monic
	// see: D.Knuth, The Art Of Computer Programming, Vol.2, 4.6.4, Ex. 34+36

	debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
			std::cout << "input: f = " << f << std::endl;);

	if (!f.is_monic()) {
		lidia_error_handler("Fp_polynomial", "square_free_decomp(...)::input"
				    "polynomial must be monic");
		return; // LC
	}

	if (f.degree() <= 1) {
		u.assign(single_factor< Fp_polynomial > (f));
		return;
	}

	Fp_polynomial t1, q;

	derivative(t1, f);
	debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
			std::cout << "derivation: f' = " << t1 << std::endl;);
	if (t1.is_zero()) {
		// f = g^p;
		sq_fr_special(u, f);
		return;
	}

	gcd(q, f, t1);
	debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
			std::cout << "gcd(f, f') = " << q << std::endl;);
	if (q.is_one()) {
		// f is square_free
		debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
				std::cout << "f is square_free" << std::endl;);
		u.assign(single_factor< Fp_polynomial > (f));
		return;
	}


	// gcd( f, f' ) != 1
	Fp_polynomial v, w, s, d;
	lidia_size_t i;

	divide(v, f, q);
	divide(w, t1, q);
	i = 0;
	lidia_size_t control = 0, f_old_degree = f.degree();
	u.kill();

	for (;;) {
		i = i + 1;
		debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
				std::cout << "  c(" << i << ") = " << v << std::endl;);

		derivative(t1, v);
		subtract(s, w, t1);
		debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
				std::cout << "  d(" << i << ") = " << s << std::endl;);

		if (s.is_zero()) {
			debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
					std::cout << "  p(" << i << ") = " << v << std::endl;);

			if (!v.is_one()) {
				u.append(v, i);
				control += i*v.degree();
			}
			break;
		}

		gcd(d, v, s);

		debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
				std::cout << "  p(" << i << ") = " << d << std::endl;);

		if (!d.is_one()) {
			divide(v, v, d);
			divide(w, s, d);
			u.append(d, i);
			control += i*d.degree();
		}
		else
			w.assign(s);
	}

	if (control == f_old_degree) {
		debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
				std::cout << "return " << u << std::endl;);
		return;
	}

	debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
			std::cout << "deg = " << control << " instead of " << f_old_degree << std::endl
			<< "factorization to be corrected :\n" << u << std::endl;);

	// now, some corrections  
	Fp_polynomial z, U, tmp;
	lidia_size_t j1, j2, J, J2;

	// t1 = prod_i u_i^(i-1)  
	power(t1, u.composite_base(0).base(), u.composite_exponent(0)-1);
	for (j1 = 1; j1 < u.no_of_composite_components(); j1++) {
		power(tmp, u.composite_base(j1).base(), u.composite_exponent(j1)-1);
		multiply(t1, t1, tmp);
	}
	debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
			std::cout << "t1 = " << t1 << std::endl;);

	divide(z, q, t1);
	factorization< Fp_polynomial > u2;
	sq_fr_special(u2, z);

	debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
			std::cout << "z = " << z << std::endl << "with squarefree decomposition :\n"
			<< u2 << std::endl;);

	factorization< Fp_polynomial > u3, empty, tmp_fact;

	for (j1 = 0; j1 < u.no_of_composite_components(); j1++) {
		U.assign(u.composite_base(j1).base());
		J = u.composite_exponent(j1);

		debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
				std::cout << " correcting exponent " << J << " : " << U << std::endl;);

		for (j2 = 0; j2 < u2.no_of_composite_components(); j2++) {
			const Fp_polynomial & U2 = u2.composite_base(j2).base();
			J2 = u2.composite_exponent(j2);

			debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
					std::cout << "  comparing with u2(" << j2 << ") = " << U2 << " ^" <<
					J2 << std::endl;);

			gcd(d, U, U2);

			debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
					std::cout << "  ggT = " << d << std::endl;);

			if (!d.is_one()) {
				u3.append(d, J + J2);

				// u2[j2] /= d  
				if (d == U2) {
					u2.replace(j2, empty);
					j2--;
				}
				else {
					divide(tmp, U2, d);
					tmp_fact.assign(single_factor< Fp_polynomial > (tmp));
					tmp_fact.power(J2);
					u2.replace(j2, tmp_fact);
				}

				// U /= d  
				if (d == U) {
					U.assign_one();
					break;
				}
				else
					divide(U, U, d);
			}
		}
		if (!U.is_one())
			u3.append(U, J);
	}

	multiply(u, u2, u3);
	debug_handler_c("Fp_polynomial", "square_free_decomp(...)", 8,
			std::cout << "return " << u << std::endl;);
}



factorization< Fp_polynomial > square_free_decomp(const Fp_polynomial &f)
{
	factorization< Fp_polynomial > F;
	square_free_decomp(F, f);
	return F;
}



factorization< Fp_polynomial >
single_factor< Fp_polynomial >::square_free_decomp() const
{
	factorization< Fp_polynomial > F;
	LiDIA::square_free_decomp(F, rep);
	return F;
}



#if 0
//obsolete

void old_square_free_decomp(factorization< Fp_polynomial > &u,
			    const Fp_polynomial &f)
{
	debug_handler_c("Fp_polynomial", "old_square_free_decomp(...)", 8,
			std::cout << "INPUT = " << f << std::endl;);
	if (f.degree() <= 1) {
		u.assign(single_factor< Fp_polynomial > (f));
		return;
	}

	Fp_polynomial d, t;
	derivative(t, f);

	if (t.is_zero()) {
		lidia_size_t i, p;
		f.modulus().sizetify(p);
		t.set_max_degree(f.degree()/p);
		for (i = 0; i <= f.degree()/p; i++)
			t[i] = f[i*p];
		debug_handler_c("Fp_polynomial", "old_square_free_decomp(...)", 8,
				std::cout << "F = " << t << " ^ P" << std::endl;);
		u.kill();
		old_square_free_decomp(u, t);
		u.power(p);
		return;
	}

	gcd(d, t, f);
	if (d.is_one()) {
		debug_handler_c("Fp_polynomial", "old_square_free_decomp(...)", 8,
				std::cout << "F IS SQUAREFREE" << std::endl;);
		u.assign(single_factor< Fp_polynomial > (f));
		return;
	}

	divide(t, f, d);
	debug_handler_c("Fp_polynomial", "old_square_free_decomp(...)", 8,
			std::cout << "APPEND t = " << t << std::endl;);
	u.kill();
	old_square_free_decomp(u, d);
	u.append(t);
	u.refine();
}
#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
