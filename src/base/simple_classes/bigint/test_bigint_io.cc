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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/bigint.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main ()
{
	bigint a;
	int    cpl;

	std::cout << std::endl;
	std::cout << "bigint::chars_per_line is " << bigint::get_chars_per_line() << std::endl;
	std::cout << std::endl;
	std::cout << "Please, enter a bigint a." << std::endl << std::endl;
	std::cout << "a = ";
	std::cin >> a;
	std::cout << std::endl;

	std::cout << "Default output: a = " << a << std::endl;
	std::cout << std::endl << std::endl;

	std::cout << "Please, enter a new value for chars_per_line : ";
	std::cin >> cpl;

	bigint::set_chars_per_line (cpl);
	std::cout << "Setting chars_per_line to " << cpl << std::endl;
	std::cout << std::endl;
	std::cout << "Now printing a again. a = " << a << std::endl << std::endl;

	cpl = 0;
	bigint::set_chars_per_line (cpl);
	std::cout << "Setting chars_per_line to " << cpl << std::endl;
	std::cout << std::endl;
	std::cout << "Now printing a again. a = " << a << std::endl;
	std::cout << std::endl;

	return 0;
}
