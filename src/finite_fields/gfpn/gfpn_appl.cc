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


#include	"LiDIA/gf_element.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	std::cout << "gf_appl   -   test application for finite field arithmetic\n\n";

	bigint p, q;
	lidia_size_t k;

	std::cout << "Enter the characteristic of the field (a prime number): ";
	std::cin >> p;
	std::cout << "Enter the degree of the field (over its prime field)  : ";
	std::cin >> k;

	p = next_prime(p-1);
	galois_field ff(p, k);
	std::cout << "\nThe field\nK = " << ff << "\nhas characteristic "
		  << ff.characteristic() << ", degree " << ff.degree()
		  << " and\n" << ff.number_of_elements() << " elements.\n";
	std::cout << "The factorization of the order of the multiplicative group is\n"
		  << std::flush << ff.factorization_of_mult_order() << std::endl;

	std::cout << "Arithmetic in K\n";
	gf_element a(ff), b(ff), c(ff);
	std::cout << "Generating random elements ...\n";
	a.randomize();
	b.randomize();
	if (a.is_zero())
		do a.randomize(); while (a.is_zero());
	std::cout << "a = " << a << std::endl;
	if (b.is_zero())
		do b.randomize(); while (b.is_zero());
	std::cout << "b = " << b << std::endl;
	std::cout << std::endl;

	std::cout << "a has order " << a.order() << std::endl;
	std::cout << "a is ";
	if (!a.is_primitive_element()) std::cout << "not ";
	std::cout << "a primitive element\n";
	std::cout << "a is ";
	if (!a.is_free_element()) std::cout << "not ";
	std::cout << "a free element\n";
	if (!a.is_square())
		std::cout << "a is not a square\n";
	else {
		c = sqrt(a);
		std::cout << "a^(1/2) = " << c << std::endl;
	}

	std::cout << "- a = " << -a << std::endl;
	std::cout << "a + b = " << a+b << std::endl;
	std::cout << "a - b = " << a-b << std::endl;
	std::cout << "a * b = " << a*b << std::endl;
	std::cout << "a / b = " << a/b << std::endl;
	std::cout << "a^(-1) = " << inverse(a) << std::endl;
	std::cout << "Tr(a) = " << a.trace() << std::endl;
	std::cout << "N(a) = " << a.norm() << std::endl;

	std::cout << "abs. degree(a) = " << a.absolute_degree() << std::endl;
	int rd = a.relative_degree();
	std::cout << "rel. degree(a) = " << rd << std::endl;
	if (rd == 1)
		std::cout << "lift(a) = " << a.lift_to_Z() << std::endl;

	std::cout << "hash(a) = " << hash(a) << std::endl;
	if (c.solve_quadratic(a, b))
		std::cout << c << " is a root of the polynomial\n"
			  << "Z^2 + " << a << "*Z + " << b << std::endl;
	else
		std::cout << "The polynomial \nZ^2 + " << a << "*Z + " << b
			  << "\ndoes not have a root in K\n";

        c = ff.generator();
	std::cout << "primitive element = " << c << std::endl;
        bool ok = true;
        for(int i = 0; i < 10; ++i) {
          if(ff.generator() != c) {
            ok = false;
            break;
          }
        }
        if(ok) {
          std::cout << "OK, galois_field::generator() returns always "
                    << "the same element.\n";
        }
        else {
          std::cout << "Error: galois_field::generator() is supposed to "
                    << "return always the same element!\n";
        }

        ok = false;
        for(int i = 0; i < 10; ++i) {
          c.assign_primitive_element(ff);
          if(c != ff.generator()) {
            ok = true;
            break;
          }
        }
        if(ok) {
          std::cout << "OK, gf_element.assign_primitive_element() "
                    << "computes different values on subsequent calls.\n";
        }
        else {
          std::cout << "Error: gf_element.assign_primitive_element() "
                    << "computed always the same value on 10 "
                    << "subsequent calls.\n";
        }

#if 0
	std::cout << "\nGenerating a primitive element :\n";
	std::cout << a.assign_primitive_elem();
	std::cout << "\nGenerating a free element :\n";
	std::cout << a.gen_free_elem();
	std::cout << "\nGenerating a free primitive element :\n";
	std::cout << a.gen_free_primitive_elem();
	std::cout << "\n";
#endif

	return 0;
}


int main(int argc, char** argv) {

#if defined(LIDIA_EXCEPTIONS)
    try {
#endif

	return main_LiDIA(argc, argv);
	
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
