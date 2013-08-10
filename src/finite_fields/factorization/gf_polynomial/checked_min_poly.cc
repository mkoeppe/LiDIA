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
#include	"LiDIA/gf_polynomial.h"
#include	"LiDIA/Fp_poly_modulus.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bool checked_min_poly(gf_polynomial& h, const gf_polynomial& g, lidia_size_t r,
		      const gf_poly_modulus& F)
	// computes the monic minimal polynomial of (g mod f) iff it has degree r;
	// then true is returned
	// otherwise, h is a multiple of the minimal polynomial; false is returned
	//
	// r = a bound on the degree of the minimal polynomial, must be < deg(F).
{
	debug_handler("gf_polynomial", "checked_min_poly(...)");

	const gf_polynomial &f = F.modulus();
	lidia_size_t i, j, row, d, k, n = f.degree();

	Fp_polynomial **M = new Fp_polynomial *[n];
	if (!M) lidia_error_handler("gf_polynomial",
				    "checked_min_poly(...)::out of memory");
	for (i = 0; i < n; i++) {
		M[i] = new Fp_polynomial[r+1];
		if (!M[i]) lidia_error_handler("gf_polynomial",
					       "checked_min_poly(...)::out of memory");
	}

	const galois_field &K = f.get_field();
	Fp_polynomial ZERO;
	ZERO.set_modulus(K.characteristic());
	Fp_polynomial t1(ZERO), t2(ZERO), *ptr;

	const Fp_polynomial &f_mod = K.irred_polynomial();
	Fp_poly_modulus F2(f_mod);

	gf_polynomial hh;
	hh.assign_one(K);

	for (j = 0; j < r+1; j++) {
		d = hh.degree();
		for (i = 0; i < n; i++) {
			if (i <= d) M[i][j].assign(hh[i].polynomial_rep());
			else        M[i][j].assign(ZERO);
		}
		if (j < r) multiply(hh, hh, g, F);
	}

	for (i = 0; i < r+1; i++) {
		row = i;
		while (row < n && M[row][i].is_zero())
			row++;
		if (row < n) {
			ptr = M[row]; M[row] = M[i]; M[i] = ptr; //swap M[i], M[row]

			invert_mod(t1, M[i][i], f_mod);
			t1.negate();
			for (j = i+1; j < r+1; j++)
				multiply(M[i][j], M[i][j], t1, F2);
			M[i][i].assign_one();
			M[i][i].negate();

			for (j = i+1; j < n; j++)
			{//M[j] += M[j][i]*M[i]
				t1.assign(M[j][i]);
				for (k = i+1; k < r+1; k++) {
					multiply(t2, M[i][k], t1, F2);
					add(M[j][k], M[j][k], t2);
				}
				M[j][i].assign_zero();
			}
		}
	}

	bool ret = true;
	h.ffield = K;
	gf_polynomial::build_frame(K);
	h.set_degree(r);
	h[r].assign_one(K);
	for (i = r-1; i >= 0; i--) {
		if (M[i][i].is_zero()) {
			ret = false;
			h[i].randomize();
		}
		else {
			t1.assign_zero();
			for (j = i+1; j < r+1; j++) {
				multiply(t2, h[j].polynomial_rep(), M[i][j]);
				add(t1, t1, t2);
			}
			h[i].set_polynomial_rep(t1);
		}
	}
	gf_polynomial::delete_frame();
	h.remove_leading_zeros();

	for (i = 0; i < n; i++) delete[] M[i];
	delete[] M;
	return ret;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
