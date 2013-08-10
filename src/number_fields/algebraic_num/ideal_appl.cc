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


#include	"LiDIA/alg_number.h"
#include	<fstream>
#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



#ifdef LIDIA_IMPLICIT_CAST_EXPLICIT
#define alg_ideal_cast(O) alg_ideal(O)
#else
#define alg_ideal_cast(O) O
#endif



std::ifstream in("alg_appl_input");



void identitytest(const alg_ideal& a, const alg_ideal& b, const alg_ideal& c)
{
	alg_ideal res1;
	module res2;
	res1 = a/b;
	res2 = (module(a))/b;
	assert(res1 == res2);
	assert((a / b) == (module(a) / b));
	assert((a * b) == (module(b) * module(a)));
	assert((a / (b / c)) == ((a * c) / b));
	assert(((a / b) * c) == ((a * c) / b));
	assert(((a * b) / b) == a);
	debug_handler_c("module_appl::identitytest", "preparing assertion 7", 3,
			std::cout << a << b << (a/b) << ((a/b)*b));
	assert(((module(a) / b) * b) == a);

	alg_ideal d;
	module d1;

	multiply(d1, b, a);
	debug_handler_c("module_appl::identitytest", "preparing assertion 8", 3,
			std::cout << a << b << (a*b) << d1);
	assert((a * b) == d1);
	multiply(d, b, a);
	debug_handler_c("module_appl::identitytest", "preparing assertion 9", 3,
			std::cout << d);
	assert(d == d1);
	divide(d1, a, b);
	debug_handler_c("module_appl::identitytest", "preparing assertion 10", 3,
			std::cout << a << b << (a/b) << d1);
	assert((a / b) == d1);
	divide(d, a, b);
	assert(d == d1);
	assert(d == d);
}



void all_test(const bigint & ii)
{
	alg_ideal a, b;

	debug_handler_l("ideal_appl", "start all_test(const bigint &)", 4);
	do {
		a.randomize(ii);
		// MM } while (a*a.denominator()== order(a.which_base()));
	} while (a*a.denominator() == alg_ideal_cast(order(a.which_base())));

	do {
		b.randomize(ii);
		// MM } while (b*b.denominator()==order(b.which_base()));
	} while (b*b.denominator() == alg_ideal_cast(order(b.which_base())));

	std::cout << "\n\n Test with \na = " << a << "\nb = " << b << std::endl;

	debug_handler_l("ideal_appl", "call identitytest from all_test", 4);
	identitytest(a, b, b);
	debug_handler_l("ideal_appl", "call identitytest again from all_test", 4);
	identitytest(a, a, b);
}



int main_LiDIA(int argc, char** argv)
{
	order O;
	std::cout << "Order:\n" << std::flush;
	in >> O;

	order ext = O.maximize();

	std::cout << O << std::endl << ext << std::endl;
	bigint i = 129099121;

	for (int j = 0; j < 10; j++)
		all_test(i);

	std::cout << "\nEnd of test\n";
	return 0;
}


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
