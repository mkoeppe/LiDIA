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


#include	"LiDIA/bigint.h"
#include	"LiDIA/arith.inl"
#include	"LiDIA/lidia_vector.h"
#include	<cstdio>
#include	<cstdlib>


#ifndef elem_type
# define elem_type bigint
#endif



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



// Compile with option -DEXP_MODE to  create
// example application for expanding vectors

#ifdef EXP_MODE
# define MODE EXPAND
#else
# define MODE FIXED
#endif



//
// Building lidia_vectors of type 'elem_type'
//


int
less_than(const elem_type & a, const elem_type & b)
{
	if (a < b)
		return -1;
	if (a == b)
		return  0;
	return  1;
}



int
greater_than (const elem_type & a, const elem_type & b)
{
	if (a > b)
		return -1;
	if (a == b)
		return  0;
	return  1;
}



int main_LiDIA(int argc, char** argv)
{
	int i;

	int test, rc, pos;

	char *fname;
	FILE *fp;

	elem_type el, tmp;

	lidia_vector< elem_type > a(0, MODE);

	test = 0;



	// ***** input, output, and operator[] *****

	if (a.get_info_mode() == vector_flags::fixed)
		std::cout << "\n Please enter a fixed vector with non-zero entries : ";
	else
		std::cout << "\n Please enter an expanding vector with non-zero entries : ";
	std::cin >> a;


	lidia_vector< elem_type > b(0, EXPAND);

	std::cout << "\n\n Enter entries of another vector with non-zero entries : \n";

	for (i = 0; i < a.size(); i++) {
		std::cout << " " << i << " : ";
		std::cin >> b[i];
	}


	elem_type *x = new elem_type [ b.size() + 1 ];

	for (i = 0; i < b.size(); i++) {
		x[i] = b[i];
	}


	std::cout << "\n if this is your second input vector : [ ";

	for (i = 0; i < b.size(); i++) {
		std::cout << x[i] << " ";
	}

	std::cout << "], \n then operator[] is ok.\n\n" << std::flush;



	lidia_vector< elem_type > c(a, MODE);
	lidia_vector< elem_type > d(x, b.size(), MODE);

	std::cout << "\n This is your first input  : " << a;
	std::cout << "\n and a copy of it          : " << c;
	std::cout << "\n This is your second input : " << b;
	std::cout << "\n and a copy of it          : " << d;

	std::cout << "\n\n ------------------------- \n\n";




	// ***** changing size and capacity *****

	std::cout << " Changing the size of a vector : \n\n";
	std::cout << " original       : " << " size = " << a.size() << " : cap. = " << a.capacity() << " : " << a << "\n";

	a.set_size (a.size() / 2);

	std::cout << " halfsized      : " << " size = " << a.size() << " : cap. = " << a.capacity() << " : " << a << "\n";

	a.set_capacity (a.capacity() * 2);

	std::cout << " double space   : " << " size = " << a.size() << " : cap. = " << a.capacity() << " : " << a << "\n";

	a.set_size (a.size() * 2);

	std::cout << " size reset     : " << " size = " << a.size() << " : cap. = " << a.capacity() << " : " << a << "\n";

	a.kill ();

	std::cout << " delete vector  : " << " size = " << a.size() << " : cap. = " << a.capacity() << " : " << a << "\n";

	a = c;

	std::cout << " and recover it : " << " size = " << a.size() << " : cap. = " << a.capacity() << " : " << a << "\n";

	std::cout << " \n Please check all the results.\n";

	std::cout << "\n ------------------------- \n\n";




	// ***** == and != *****

	std::cout << " operator == and operator != ... \n\n" << std::flush;

	if (!(a == c)) {
		std::cout << " ERROR : operator == failed < ----- \n" << std::flush;
		test++;
	}


	a[0] = b[0] + static_cast<elem_type>(1);

	if (!(a != b)) {
		std::cout << " ERROR : operator != failed < ------ \n" << std::flush;
		test++;
	}




	std::cout << " arithmetical operations ... \n\n";

	lidia_vector< elem_type > z1 (0, MODE);
	lidia_vector< elem_type > z2 (0, MODE);
	lidia_vector< elem_type > z3 (0, MODE);


	// *****  ADDITION AND SUBTRACTION  *****

	a = c;

	z1 = a  +  b;
	z2 = z1 -  a;

	if (z2 != b) {
		std::cout << " ERROR : operator+  or operator- failed < ----- \n" << std::flush;
		test ++;
	}


	z1 = a  +  static_cast<elem_type>(10);
	z2 = z1 -  static_cast<elem_type>(10);

	if (z2 != a) {
		std::cout << " ERROR : operator+  or operator- failed  (scalar) < ----- \n" << std::flush;
		test ++;
	}


	z1 = -a;
	z2 = b - a;
	z3 = b + z1;

	if (!(z2 == z3)) {
		std::cout << " ERROR :  operator- (unary) failed < ----- \n" << std::flush;
		test++;
	}


	add      (z1, a , b);
	subtract (z2, z1, a);

	if (z2 != b) {
		std::cout << " ERROR : addition or subtraction failed < ----- \n" << std::flush;
		test ++;
	}



	add      (z1, a , static_cast<elem_type>(10));
	subtract (z2, z1, static_cast<elem_type>(10));

	if (z2 != a) {
		std::cout << " ERROR : addition or subtraction failed  (scalar) < ----- \n" << std::flush;
		test ++;
	}


	negate   (z1, a);
	subtract (z2, b, a);
	add      (z3, b, z1);

	if (!(z2 == z3)) {
		std::cout << " ERROR :  negate failed < ----- \n" << std::flush;
		test++;
	}



	// *****  MULTIPLICATION AND DIVISION  *****


	compwise_multiply(z1, a, b);
	compwise_divide(z2, z1, a);

	if (z2 != b) {
		std::cout << " ERROR : operator*  or operator/ failed < ----- \n" << std::flush;
		test ++;
	}


	z1 = a  *  static_cast<elem_type>(10);
	z2 = z1 /  static_cast<elem_type>(10);

	if (z2 != a) {
		std::cout << " ERROR : operator*  or operator/ failed  (scalar) < ----- \n" << std::flush;
		test ++;
	}



	compwise_multiply (z1, a , b);
	compwise_divide   (z2, z1, a);

	if (z2 != b) {
		std::cout << " ERROR : multiplication or division failed < ----- \n" << std::flush;
		test ++;
	}


	multiply (z1, a , static_cast<elem_type>(10));
	divide   (z2, z1, static_cast<elem_type>(10));

	if (z2 != a) {
		std::cout << " ERROR : multiplication or division failed (scalar) < ----- \n" << std::flush;
		test ++;
	}



	// *****  SUM OF SQUARES / INNER PRODUCT  *****

	i = 20;

	z3.set_capacity (i);
	z3[0] = static_cast<elem_type>(1);

	for (i = 1; i < z3.capacity(); i++)

		z3[i] = z3[i-1] + static_cast<elem_type>(1); // z3[i] == i+1


	tmp = (i * (i + i + 1) * (i + 1)) / 6;

	z2 = -z3;


	if (tmp != sum_of_squares(z3)) {
		std::cout << " ERROR : sum_of_squares failed < ----- \n" << std::flush;
		test++;
	}

	tmp = -tmp;


	if (tmp != z2 * z3) {
		std::cout << " ERROR : operator* failed < ----- \n" << std::flush;
		test++;
	}



	// *****  SORTING AND SEARCHING  *****


	std::cout << " sorting and searching ... \n\n";

	z1 = a;
	z2 = b;

	z1.sort (SORT_VECTOR_UP);
	z2.sort (SORT_VECTOR_DOWN, 0, z2.size() / 2);

	z3 = z2;
	z3.set_size ((z2.size() / 2) + 1);

	std::cout << " First input in ascending order                 : " << z1 << std::endl;
	std::cout << " First part of second input in descending order : " << z3 << std::endl;

	std::cout << "\n Please check whether the vectors are sorted correctly.\n\n";



	z2.sort (0, z2.size()-1);

	std::cout << " linear search ... " << std::endl;

	tmp = sum_of_squares (z1);
	el = z1[ z1.size()-1 ];
	rc = z1.linear_search (el, pos);

	if (! rc) {
		std::cout << " ERROR : linear search failed (" << el << " not found) < ----- \n" << std::endl;
		test++;
	}


	el = tmp + static_cast<elem_type>(1);
	rc = z1.linear_search (el, pos);

	if (rc) {
		std::cout << " ERROR : linear search failed (non existing element " << el << " found) < ----- \n" << std::endl;
		test++;
	}



	std::cout << " binary search ... \n" << std::endl;

	el = z1[ z1.size()-1 ];
	rc = z1.bin_search (el, pos);

	if (! rc) {
		std::cout << " ERROR : binary (up) search failed (" << el << " not found) < ----- \n" << std::endl;
		test++;
	}


	el = tmp + static_cast<elem_type>(1);
	rc = z1.bin_search (el, pos);

	if (rc) {
		std::cout << " ERROR : binary (up) search failed (non existing element ";
		std::cout << el << " found) < ----- \n" << std::endl;
		test++;
	}


	z1.sort (SORT_VECTOR_DOWN);

	el = z1[ z1.size()-1 ];
	rc = z1.bin_search (el, pos, SORT_VECTOR_DOWN);

	if (! rc) {
		std::cout << " ERROR : binary (down) search failed (" << el << " not found) < ----- \n" << std::endl;
		test++;
	}


	el = tmp + static_cast<elem_type>(1);
	rc = z1.bin_search (el, pos);

	if (rc) {
		std::cout << " ERROR : binary (down) search failed (non existing element ";
		std::cout << el << " found) < ----- \n" << std::endl;
		test++;
	}




	//  *****  SORTING/SEARCHING WITH COMPARE-FUNCTIONS  *****

	z1 = a;
	z2 = b;

	z1.sort (greater_than);
	z2.sort (less_than);

	std::cout << " using compare functions ...\n" << std::endl;

	std::cout << " " << a << " in descending  order : " << z1 << "\n" << std::endl;
	std::cout << " " << b << " in ascending   order : " << z2 << std::endl;


	std::cout << "\n binary search ... \n" << std::endl;

	el = z1[ z1.size()-1 ];
	rc = z1.bin_search (el, pos, greater_than);

	if (! rc) {
		std::cout << " ERROR : binary (gt) search failed (" << el << " not found) < ----- \n" << std::endl;
		test++;
	}


	el = tmp + static_cast<elem_type>(1);
	rc = z1.bin_search (el, pos, greater_than);

	if (rc) {
		std::cout << " ERROR : binary (gt) search failed (non existing element ";
		std::cout << el << " found) < ----- \n" << std::endl;
		test++;
	}

	z1.sort (less_than);

	el = z1[ z1.size()-1 ];
	rc = z1.bin_search (el, pos, less_than);


	if (! rc) {
		std::cout << " ERROR : binary (lt) search failed (";
		std::cout << el << " not found) < ----- \n" << std::endl;
		test++;
	}


	el = tmp + static_cast<elem_type>(1);
	rc = z1.bin_search (el, pos, less_than);

	if (rc) {
		std::cout << " ERROR : binary (lt) search failed (non existing element ";
		std::cout << el << " found) < ----- \n" << std::endl;
		test++;

		std::cout << z1 << "\n" << el << "\n" << pos << std::endl;
		exit(1);
	}



	// ***** READING FROM AND WRITING TO FILE *****

	std::cout << " reading from and writing to file ... \n" << std::endl;

	fname = tempnam (".", "LVECT");


	// write a to file fname twice, read the to copies of a in z1 and z2, and
	// compare these vectors with a  ( ASCII - format )

	if ((fp = fopen (fname, "w")) != NULL) {
		a.print_to_file (fp);
		a.print_to_file (fp);

		fclose (fp);

		if ((fp = fopen (fname, "r")) != NULL) {
			z1.scan_from_file (fp);
			z2.scan_from_file (fp);

			if (a != z1 || a != z2) {
				std::cout << "a.print_to_file (FILE*) or a.scan_from_file (FILE*) failed !\n" << std::endl;
				test++;
			}

			fclose (fp);
		}
		else {
			std::cout << "sorry, can't open " << fname << " for reading in ASCII - format ! \n";
			std::cout << "test was stopped. \n\n" << std::endl;
			exit (1);
		}
	}
	else {
		std::cout << "sorry, can't open " << fname << " for writing in ASCII - format ! \n";
		std::cout << "test was stopped. \n\n" << std::endl;
		exit (1);
	}


	// the same test for binary - format

	if ((fp = fopen (fname, "wb")) != NULL) {
		a.write_to_file (fp);
		a.write_to_file (fp);

		fclose (fp);

		if ((fp = fopen (fname, "rb")) != NULL) {
			z1.read_from_file (fp);
			z2.read_from_file (fp);

			if (a != z1 || a != z2) {
				std::cout << "a.write_to_file (FILE*) or a.read_from_file (FILE*) failed !\n" << std::endl;
				test++;
			}

			fclose (fp);
		}
		else {
			std::cout << "sorry, can't open " << fname << " for reading in binary - format ! \n";
			std::cout << "test was stopped. \n\n" << std::endl;
			exit (1);
		}
	}
	else {
		std::cout << "sorry, can't open " << fname << " for writing in binary - format ! \n";
		std::cout << "test was stopped. \n\n" << std::endl;
		exit (1);
	}


	std::remove (fname);
	std::free(fname);


	if (test == 0) {
		std::cout << "\n No errors detected. \n\n" << std::endl;
	}
	else {
		std::cout << "\n " << test << " error(s) detected. \n\n" << std::endl;
	}

	delete[] x;

	return 0;
}



// must instantiate lidia_vector:
template class lidia_vector< elem_type >;


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
