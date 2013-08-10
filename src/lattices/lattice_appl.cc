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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/bigint_lattice.h"
#include	"LiDIA/random_generator.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	random_generator rg;
	bigint_lattice A, B, C, G, Old, T;
	bigint values;
	lidia_size_t dim, ran1, ran2;
	lattice_info li;
	double y;
	timer t;
	std::cout << "Intervall for entries : ";
	std::cin >> values;
	std::cout << "Intervall for dimension : ";
	std::cin >> dim;
	std::cout << "y (0.5 < y <= 1) : ";
	std::cin >> y;
	std::cout << "generating lattice A ..." << std::endl;
	rg >> ran1 >> ran2;
	A.resize(dim+(ran1%dim), dim+(ran2%dim));
	A.randomize(values);
	std::cout << "A: " << A << std::endl;
	rg >> ran1;
	if (ran1%2)
	{
		A.set_red_orientation_columns();
		std::cout << "Reducing over columns ..." << std::endl;
	}
	else
	{
		A.set_red_orientation_rows();
		std::cout << "Reducing over rows ..." << std::endl;
	}
	std::cout << "Dimension of A : " << A.get_no_of_rows() << " x "
		  << A.get_no_of_columns() << std::endl;
	std::cout << "Checking A for basis : ";
	if (A.check_basis())
		std::cout << "ok" << std::endl;
	else
		std::cout << "failed" << std::endl;

	std::cout << "Checking A for gram : ";
	if (A.check_gram())
		std::cout << "ok" << std::endl;
	else
		std::cout << "failed" << std::endl;
	if (A.get_red_orientation())
	{
		std::cout << "Generating Gram matrix (A^T*A) ..." << std::endl;
#if defined (_MSC_VER)
		G.assign(bigint_matrix(A.trans())*A);
#else
		G.assign(A.trans()*A);
#endif
	}
	else
	{
		std::cout << "Generating Gram matrix (A*A^T) ..." << std::endl;
#if defined (_MSC_VER)
		G.assign(A*bigint_matrix(A.trans()));
#else
		G.assign(A*A.trans());
#endif
	}

	G.set_gram_flag();
	std::cout << "Checking G for basis : ";
	if (G.check_basis())
		std::cout << "ok" << std::endl;
	else
		std::cout << "failed" << std::endl;
	Old.assign(A);

	A.lll(T, y, li);
	G.lll(y, li);
	if (A.get_red_orientation())
		C.assign(Old*T);
	else
		C.assign(T*Old);
	if (A == C)
		std::cout << "1. test ok" << std::endl;
	else
		std::cout << "1. test failed" << std::endl;

	if (A.get_red_orientation())
#if defined (_MSC_VER)
		C.assign(bigint_matrix(A.trans())*A);
	else
		C.assign(A*bigint_matrix(A.trans()));
#else
		C.assign(A.trans()*A);
	else
		C.assign(A*A.trans());
#endif
	if (G == C)
		std::cout << "2. test ok" << std::endl;
	else
		std::cout << "2. test failed" << std::endl;
	if (G.check_gram())
		std::cout << "3. test ok" << std::endl;
	else
		std::cout << "3. test failed" << std::endl;;
	std::cout << "(chosen) : " << y << std::endl;
	if (A.get_red_orientation())
		A.set_no_of_columns(li.lll.rank);
	else
		A.set_no_of_rows(li.lll.rank);
	std::cout << "(lll_check_search) : " << A.lll_check_search() << std::endl;
	std::cout << "(lll_check) : "
		  << ((A.lll_check(y) == true) ? "OK" : "failed") << std::endl;
	std::cout << std::endl;
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
