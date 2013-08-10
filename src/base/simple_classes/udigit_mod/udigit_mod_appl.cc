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
//	Author	: Thorsten Rottschaefer (TR)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/udigit_mod.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	udigit_mod::set_modulus(17);

	udigit_mod c, d, e;
	c = 10;
	d = 9;

	std::cout << " \n small Test of udigit_mod: \n\n";

	add(e, d, c);
	std::cout << "compute: (10+9) mod 17 result = " << e << "\n\n"; // output should be 2 (10+9%17 = 2)

	subtract(e, c, d);
	std::cout << "compute: (10-9) mod 17 result = " << e << "\n\n"; // output should be 1 (10-9%17 = 1)

	multiply(e, c, d);
	std::cout << "compute: (10*9) mod 17 result = " << e << "\n\n"; // output should be 5 (10*9%17 = 5)

	divide(e, c, d);
	std::cout << "compute: (10/9) mod 17 result = " << e << "\n\n"; // output should be 3 (inverse of 9 is 2 and 2*10%17 = 3)

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
