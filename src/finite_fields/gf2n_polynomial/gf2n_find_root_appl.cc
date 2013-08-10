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
//	Author	: Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/gf2n_poly_modulus.h"
#include	"LiDIA/base_vector.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	int degree;
	int j , i, deg, tests;
	gf2n_polynomial f, res, xq;
	gf2n_poly_modulus fpm;

	std::cout << "\nTest program for finding roots of polynomials over Gf2n";
	std::cout << "\n=======================================================\n\n";
	std::cout << "\nThe program computes for random polynomials all roots in the";
	std::cout << "\nbase field and checks their correctness. Then it constructs";
	std::cout << "\npolynomials from their roots and checks the root finding routine.\n";
	std::cout << "\nInput Extension Degree of Field : "; std::cin >> degree;
	std::cout << "\nInput Polynomial Degree : "; std::cin >> deg;
	std::cout << "\nInput # Tests : "; std::cin >> tests;

	gf2n_init(degree);
	gf2n b;
	gf2n_polynomial g, gg, xqi;

	for (i = 1; i <= tests; i++) {
		f.randomize(deg);
		f.make_monic();

		std::cout << "\n-------------------------------------------------";
		std::cout << "\n\nTest No. " << i << "\n";

		fpm.build(f);
		Xq(xq, fpm);
		xqi.assign(xq);

		std::cout << "\n\nf = " << f;

		for (j = 1; j <= deg && f.degree() != 1; j++) {
			add(g, xqi, gf2n_polynomial(1));
			gcd(g, g, f);
			std::cout << "\ngcd(X^(q^" << j << ") + X, f) = " << g << "\n" << std::flush;

			if (g.degree() > j) {
				gf2n_polynomial xxq;
				fpm.build(g);
				remainder(xxq, xq, fpm.modulus());
				EDF_one_factor(gg, xxq, fpm, j);
			}
			else
				if (!g.is_one())
					std::cout << "==> single factor of degree " << j;

			if (g.degree() >= 1) {
				divide(f, f, g);
				if (f.is_one())
					break;
				fpm.build(f);
				remainder(xq, xq, f);
				remainder(xqi, xqi, f);
			}
			compose(xqi, xqi, xq, fpm);
		}
	}
	std::cout << std::endl;
	return 0;
}

//  base_vector< gf2n > bv;

//  find_all_roots(bv, f);
//  if (bv.size() != f.degree())
//  {
//  	std::cerr << "\nERROR: not all roots were found, aborting.\n\n";
//  	exit(1);
//  }

//  std::cout << "\nRoots : " << bv << std::flush;

//  for (j = 0; j < bv.size(); j++)
//  	if (!f(bv[j]).is_zero())
//  {
//  	std::cerr << "\nERROR: Computed element " << bv[j] << " is not a root, ";
//  	std::cerr << " aborting.\n\n";
//  	exit(1);
//  }
//  }
//  std::cout << "\n\n";
//  std::cout << "\nNow constructiopn of polynomials with linear factors ...\n";
//  std::cout << "\nInput new Polynomial Degree : "; std::cin >> deg;

//  for (i = 1; i <= tests; i++)
//  {
//  f.assign_one();

//  for (j = 1; j <= deg; j++)
//  {
//  	b.randomize();
//  	multiply_by_linear(f, f, b);
//  }

//   std::cout << "\n-------------------------------------------------";
//   std::cout << "\n\nTest No. " << i << "\n";

//   fpm.build(f);
//   Xq(xq, fpm);
//   add(xq, xq, gf2n_polynomial(1));
//   gcd(f, xq, f);

//   base_vector< gf2n > bv;

//   find_all_roots(bv, f);
//   if (bv.size() != f.degree())
//   {
//  	 std::cerr << "\nERROR: not all roots were found, aborting.\n\n";
//  	 exit(1);
//   }


//   std::cout << "\nRoots : " << bv << std::flush;

//   for (j = 0; j < bv.size(); j++)
//  	 if (!f(bv[j]).is_zero())
//  	 {
//  		 std::cerr << "\nERROR: Computed element " << bv[j] << " is not a root, ";
//  		 std::cerr << " aborting.\n\n";
//  		 exit(1);
//  	 }
//  }
//  std::cout << "\n\n==> NO ERROR FOUND.\n\n";
//  return 0;
//  }


int main(int argc, char** argv) {

#if defined(LIDIA_EXCEPTIONS)
    try {
#endif

	main_LiDIA(argc, argv);
	
#if defined(LIDIA_EXCEPTIONS)
    }
    catch(basic_error const& ex) {
	ex.traditional_error_handler();
	return 1;
    }
    catch(std::exception const& ex) {
	std::cerr << "unexpected exception: " << ex.what() << "\n";
	return 1;
    }
#endif
    
}
