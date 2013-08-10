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
//	$Id: bit_reverse_table_appl.cc,v 2.3 2004/06/15 10:19:38 lidiaadm Exp $
//
//	Author	: Thorsten Rottschaefer (TR)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/finite_fields/bit_reverse_table.h"



const int LOG_SIZE = 3;

// Changed by G.A.: SIZE replaced with SIZE_

const int SIZE_ = 1 << LOG_SIZE; // must be a power of 2



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



void print(const udigit a[SIZE_])
{
	lidia_size_t i;
	for (i = 0; i < SIZE_; i++)
		std::cout << a[i] << " ";
	std::cout << std::endl;
}



int main_LiDIA(int argc, char** argv)
{
	lidia_size_t i;
	udigit a[SIZE_], b[SIZE_], c[SIZE_];
	for (i = 0; i < SIZE_; i++)
		a[i] = i;

	std::cout << "\noriginal table:\n";
	print(a);

	// initialize an object
	bit_reverse_table bitcopy;

	// perform a bitreverse copy of a in b
	bitcopy.copy(b, a, LOG_SIZE);

	std::cout << "\nbitreverse copy of the original table:\n";
	print(b);

	// copy backwards...
	bitcopy.copy(c, b, LOG_SIZE);

	std::cout << "\nbitreverse copy of the previous table:\n";
	print(c);

	bool ok = true;
	for (i = 0; i < SIZE_; i++)
		if (a[i] != c[i])
			ok = false;

	if (ok)
		std::cout << "\nresult is ok.\n";
	else
		std::cout << "\nERROR: result is not correct !\n";

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
