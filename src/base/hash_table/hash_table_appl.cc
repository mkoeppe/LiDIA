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
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/hash_table.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



bigint
ikey(const int & G)
{
	return bigint(G);
}



int main_LiDIA(int argc, char** argv)
{
	hash_table< int > HT, HT2;
	int i, x, *y;

	HT.initialize(20);
	HT.set_key_function(&ikey);

	std::cout << "#buckets = " << HT.no_of_buckets() << "\n";
	std::cout << "current size = " << HT.no_of_elements() << "\n\n";

	for (i = 0; i < 30; ++i)
		HT.hash(i);
	for (i = -50; i > -60; --i)
		HT.hash(i);

	std::cout << "\ncurrent size = " << HT.no_of_elements() << "\n";
	x = HT.last_entry();
	std::cout << "last entry = " << x << "\n";
	y = HT.search(x);
	if (y)
		std::cout << "search succeeded!  Last entry is in the table.\n";
	else
		std::cout << "search failed!  Please report this bug!\n";

	std::cout << "\nContents of HT:\n";
	std::cout << HT << "\n";

	std::cout << "\ncopying...\n\n";
	HT2 = HT;
	std::cout << "contents of new table:\n";
	std::cout << HT2 << "\n";
	std::cout << "#buckets = " << HT2.no_of_buckets() << "\n";
	std::cout << "current size = " << HT2.no_of_elements() << "\n\n";

	std::cout << "remove function...\n";
	for (i = 7; i < 50; i += 7) {
		std::cout << "remove " << i << "\n";
		HT2.remove(i);
	}
	std::cout << "contents of new table:\n";
	std::cout << HT2 << "\n";
	std::cout << "#buckets = " << HT2.no_of_buckets() << "\n";
	std::cout << "current size = " << HT2.no_of_elements() << "\n\n";

	std::cout << "emptying it...\n";
	std::cout << "contents of new table:\n";
	HT2.empty();
	std::cout << HT2 << "\n";
	std::cout << "#buckets = " << HT2.no_of_buckets() << "\n";
	std::cout << "current size = " << HT2.no_of_elements() << "\n\n";

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
