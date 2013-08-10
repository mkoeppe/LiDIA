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


#include	"LiDIA/indexed_hash_table.h"



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
	indexed_hash_table< int > IHT, IHT2;
	int i, x, *y;
	lidia_size_t j;

	IHT.initialize(20);
	IHT.set_key_function(&ikey);

	std::cout << "#buckets = " << IHT.no_of_buckets() << "\n";
	std::cout << "current size = " << IHT.no_of_elements() << "\n\n";

	for (i = 0; i < 30; ++i)
		IHT.hash(i);
	for (i = -50; i > -60; --i)
		IHT.hash(i);

	std::cout << "\ncurrent size = " << IHT.no_of_elements() << "\n";
	x = IHT.last_entry();
	std::cout << "last entry = " << x << "\n";
	y = IHT.search(x);
	if (y)
		std::cout << "search succeeded!  Last entry is in the table.\n";
	else
		std::cout << "search failed!  Please report this bug!\n";

	std::cout << "\nsequential access:\n";
	for (j = 0; j < 20; j += 5)
		std::cout << "element " << j << " = " << IHT.member(j) << "\n";
	for (j = 20; j <= 35; j += 5)
		std::cout << "element " << j << " = " << IHT[j] << "\n";

	std::cout << "Contents of IHT:\n";
	std::cout << IHT << "\n";

	std::cout << "\ncopying...\n\n";
	IHT2 = IHT;
	std::cout << "contents of new table:\n";
	std::cout << IHT2 << "\n";
	std::cout << "#buckets = " << IHT2.no_of_buckets() << "\n";
	std::cout << "current size = " << IHT2.no_of_elements() << "\n\n";

	std::cout << "remove_from function...\n";
	for (j = 35; j >= 0; j -= 5) {
		std::cout << "remove element " << j << " = " << IHT2.member(j) << "\n";
		IHT2.remove_from(j);
	}
	std::cout << "remove function...\n";
	for (i = 7; i < 50; i += 7) {
		std::cout << "remove " << i << "\n";
		IHT2.remove(i);
	}
	std::cout << "contents of new table:\n";
	std::cout << IHT2 << "\n";
	std::cout << "#buckets = " << IHT2.no_of_buckets() << "\n";
	std::cout << "current size = " << IHT2.no_of_elements() << "\n\n";
	std::cout << "contents of new table (sequentially):\n";
	IHT2.output_style(0);
	std::cout << IHT2 << "\n";

	std::cout << "emptying it...\n";
	std::cout << "contents of new table:\n";
	IHT2.empty();
	std::cout << IHT2 << "\n";
	std::cout << "#buckets = " << IHT2.no_of_buckets() << "\n";
	std::cout << "current size = " << IHT2.no_of_elements() << "\n\n";

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
