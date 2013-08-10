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
//	Author	: Victor Shoup (VS), Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/Fp_poly_multiplier.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//******************************************************************
//
//			  class Fp_poly_multiplier
//
//******************************************************************

Fp_poly_multiplier::Fp_poly_multiplier(const Fp_poly_multiplier &B)
{
	if (B.poly_mod_ptr != 0)
		build(B.multiplier(), *B.poly_mod_ptr);
}



Fp_poly_multiplier & Fp_poly_multiplier::operator = (const Fp_poly_multiplier &B)
{
	if (B.poly_mod_ptr != 0)
		build(B.multiplier(), *B.poly_mod_ptr);
	else
		lidia_error_handler("Fp_poly_multiplier",
				    "operator = (Fp_poly_multiplier&)::operand must be initialized");
	return *this;
}



void Fp_poly_multiplier::build(const Fp_polynomial& x, const Fp_poly_modulus& F)
{
	debug_handler("Fp_poly_multiplier", "build (Fp_polynomial&, Fp_poly_modulus&)");

	poly_mod_ptr = &F;

	lidia_size_t db;
	lidia_size_t n = F.n;

	remainder(b, x, F);
	db = x.degree();

	if (!F.use_FFT || db <= F.crov) {
		use_FFT = false;
		return;
	}

	use_FFT = true;

	Fp_polynomial P1;
	P1.set_max_degree(n-1);

	modular_fft_rep R1;
	R1.init(F.fd);
	R1.set_size(F.l);

	B2.init(F.fd);
	B2.set_size(F.k);
	for (lidia_size_t index = 0; index < F.fd.number_of_primes(); index++) {
		R1.to_modular_fft_rep(b, index);
		reduce(B2, R1, F.k, index); //initializing B2
		multiply(R1, F.HRep, R1, index);
		R1.from_modular_fft_rep(n-1, 2*n-3, index);
	}
	R1.get_result(P1, n-1, 2*n-3);

	B1.init(F.fd);
	B1.set_size(F.l);
	B1.to_fft_rep(P1); //initializing B1

}



// If you need to compute a * b % f for a fixed b, but for many a's
// (for example, computing powers of b modulo f), it is
// much more efficient to first build a Fp_poly_multiplier B for b,
// and then use the routine below.
// x = (a * b) % f
void multiply(Fp_polynomial& x, const Fp_polynomial& a,
	      const Fp_poly_multiplier& B, const Fp_poly_modulus& F)
{
	debug_handler("Fp_poly_multiplier", "multiply(Fp_polynomial&, Fp_polynomial&, Fp_poly_multiplier&, Fp_poly_modulus&)");

	a.comp_modulus(F.f, "Fp_poly_multiplier::multiply");

	if (B.poly_mod_ptr == 0 || B.poly_mod_ptr != &F) {
		lidia_error_handler("Fp_poly_multiplier", "multiply(Fp_polynomial&, "
				    "Fp_polynomial&, Fp_poly_multiplier&, Fp_poly_modulus&)::wrong "
				    "Fp_poly_modulus, or Fp_poly_multiplier not initialized");
		return;
	}

	lidia_size_t n = F.n;
	lidia_size_t da, index;

	da = a.degree();

	if (da >= n) {
		lidia_error_handler("Fp_poly_multiplier", "multiply(Fp_polynomial&, "
				    "Fp_polynomial&, Fp_poly_multiplier&, Fp_poly_modulus&)::degree of "
				    "Fp_polynomial must be < degree of Fp_poly_modulus");
		return;
	}

	if (da <= F.crov || !B.use_FFT) {
		Fp_polynomial P1;
		if (F.use_FFT)
			fft_mul(P1, a, B.b);
		else
			plain_mul(P1, a, B.b);
		plain_rem(x, P1, F.f);
		return;
	}

	Fp_polynomial P1;
	Fp_polynomial P2;
	P1.set_max_degree(n-1);
	P2.set_max_degree(n-1);

	lidia_size_t num_primes = F.fd.number_of_primes();
	if (num_primes == 0) {
		fft_rep Ra, Rb;
		Ra.init(F.fd);
		Rb.init(F.fd);

		Ra.set_size(F.l);
		Rb.set_size(F.l);
		Ra.to_fft_rep(a);
		multiply(Rb, B.B1, Ra);
		Rb.from_fft_rep(P1, n-1, 2*n-3);

		Rb.set_size(F.k);
		reduce(Rb, Ra, F.k);
		multiply(Rb, B.B2, Rb);

		Ra.set_size(F.k);
		Ra.to_fft_rep(P1);
		multiply(Ra, F.FRep, Ra);
		subtract(Rb, Rb, Ra);

		Rb.from_fft_rep(x, 0, n-1);
	}
	else {
		modular_fft_rep Ra, Rb; // more space efficient version
		fft_rep R;
		Ra.init(F.fd);
		Rb.init(F.fd);
		R.init(F.fd);

		Ra.set_size(F.l);
		Rb.set_size(F.l);
		R.set_size(F.k);
		for (index = 0; index < F.fd.number_of_primes(); index++) {
			Ra.to_modular_fft_rep(a, index);
			reduce(R, Ra, F.k, index);
			multiply(Rb, B.B1, Ra, index);
			Rb.from_modular_fft_rep(n-1, 2*n-3, index);
		}
		Rb.get_result(P1, n-1, 2*n-3);

		multiply(R, B.B2, R);

		Ra.set_size(F.k);
		Rb.set_size(F.k);
		for (index = 0; index < F.fd.number_of_primes(); index++) {
			Ra.to_modular_fft_rep(P1, index);
			multiply(Ra, F.FRep, Ra, index);
			subtract(Rb, R, Ra, index);

			Rb.from_modular_fft_rep(0, n-1, index);
		}
		Rb.get_result(x, 0, n-1);
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
