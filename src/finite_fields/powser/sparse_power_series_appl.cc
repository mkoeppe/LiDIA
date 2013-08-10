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


#include	"LiDIA/bigmod.h"
#include	"LiDIA/bigrational.h"
#include	"LiDIA/sparse_power_series.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif


typedef bigmod elem_type;

int main_LiDIA(int argc, char** argv)
{
        lidia_size_t const max_exponent = 15;
        bigmod::set_modulus (bigint(65537));

	sparse_power_series< elem_type > e;

	base_vector< elem_type > x;
	base_vector< lidia_size_t > exp_v;

	elem_type *x_pt;
	lidia_size_t *exp_v_pt;

	elem_type one = static_cast<elem_type>(1);
	elem_type zero = static_cast<elem_type>(0);
	elem_type tmp;

	long  err_num = 0;
	lidia_size_t  f, l;
	lidia_size_t  i, j, k;


	std::cout << "\n\n";
	std::cout << " This program runs several tests to verfiy\n";
	std::cout << " the correctness of the class template "
		  << "sparse_power_series< T > .\n";
	std::cout << " The  template is instantiated for LiDIA::bigmod. All\n";
	std::cout << " computations are done mod "
		  << bigmod::modulus() << "\n\n";
	std::cout << " If an error occurs, a corresponding message\n";
	std::cout << " is displayed. At the end, the program reports\n";
	std::cout << " the number of detected errors.\n\n" << std::endl;





	// *************************************************************
	// ********* plain constructor, input, output ******************
	// *************************************************************

	std::cout << " plain constructor, input, output ..." << "\n" << std::endl;

	sparse_power_series< elem_type > a;

	std::cout << " Please enter a series with at least two non-zero coefficients\n";
	std::cout << " (format: [ [ [bigint1 exp1] [bigint2 exp2] ... [bigintn expn] ] (int) last_exponent ]) :\n\n ";

	std::cin >> a;
	std::cout << "\n\n" << " a = " << a
		  << " mod " << bigmod::modulus() << "\n" << std::endl;

	std::cout << " Please, check whether the output for a is correct." << std::endl;
	std::cout << " ==================================================\n" << std::endl;

	// *************************************************************
	// ********* operator==, operator!=, constructor ***************
	// *************************************************************

	std::cout << " operator == , operator != , constructor ..." << "\n" << std::endl;

	sparse_power_series< elem_type > b (a);

	if (a == b) {
		if (a != b) {
			err_num++;
			std::cout << "operator != failed\n" << std::endl;
		}
	}
	else {
		err_num++;
		std::cout << "constructor "
			  << "sparse_power_series(sparse_power_series const&) "
			  << "or operator == failed\n" << std::endl;
	}



	// *************************************************************
	// constructor, get_first(), get_last(), operator[]
	// *************************************************************

	std::cout << " constructor, get_first(), get_last(), operator[] ..." << "\n" << std::endl;

	l = max_exponent;

	sparse_power_series< elem_type > c   (one, l);
	sparse_power_series< elem_type > One (one, l);



	// test the constructor with an other coefficient

	l = max_exponent;

	sparse_power_series< elem_type > d    (zero, l);
	sparse_power_series< elem_type > Zero (zero, l);

	if (d.get_first() == l && d.get_last() == l) {
		if (d[l] != zero) {
			err_num++;
			std::cout << "sparse_power_series(T, lidia_size_t) or operator[] failed.\n" << std::endl;
		}
		else if (d[-1] != zero || d[l-1] != zero) {
			err_num++;
			std::cout << "operator[] failed.\n" << std::endl;
		}
	}
	else {
		err_num++;
		std::cout << "sparse_power_series(T, lidia_size_t), get_first() or get_last() failed.\n" << std::endl;
	}



	// *************************************************************
	// **************** get(), get_coeff () ************************
	// *************************************************************

	std::cout << " get(), get_coeff() ... \n" << std::endl;

	f = a.get_first ();
	l = a.get_last  ();


	// the same test with x_pt, exp_v_pt (pointer instead of vector)

	x_pt = new elem_type [1];
	exp_v_pt = NULL;

	a.get(x_pt, exp_v_pt, i);
	b.clear();
	b.set(x_pt, exp_v_pt, i, a.get_last());

	delete[] x_pt;

	if (a != b) {
		err_num++;
		std::cout << "get(T*&, lidia_size_t*&, lidia_size_t&) or"
			"set(T*, lidia_size_t*, lidia_size_t, lidia_size_t) failed\n" << std::endl;
	}


	// check the get_coeff() - function

	for (i = f; i <= l; i++) {
		a.get_coeff (tmp, i);

		if (tmp != a[i]) {
			i = l;
			err_num++;
			std::cout << "get_coeff(T, lidia_size_t) failed.\n" << std::endl;
		}
	}


	// *************************************************************
	// *************** clear(), set(), set_coeff () ****************
	// *************************************************************

	std::cout << " clear(), set(), set_coeff() ...\n" << std::endl;

	f = a.get_first ();
	l = a.get_last  ();

	b.clear ();

	b.set (one, max_exponent);

	if (!(b == One)) {
		err_num++;
		std::cout << "set(T, lidia_size_t) failed.\n" << std::endl;
	}


	// test set_coeff()

	b.clear ();
	// b.set_first ( f );

	for (i = f; i <= l; i++) {
		b.set_coeff (a[i], i);
	}

	if (!(b == a)) {
		err_num++;
		std::cout << "set_coeff(T, lidia_size_t) failed.\n" << std::endl;
	}


	b.set_coeff (zero, f);
	b.set_coeff (a[f], f);

	if (!(b == a)) {
		err_num++;
		std::cout << "set_coeff(T, lidia_size_t) failed.\n" << std::endl;
	}


	b.set_coeff (one , f-5);
	b.set_coeff (zero, f-5);

	if (!(b == a)) {
		err_num++;
		std::cout << "set_coeff(T, lidia_size_t) failed.\n" << std::endl;
	}

	// *************************************************************
	// ********************** operator= ****************************
	// *************************************************************

	std::cout << " operator = ...\n" << std::endl;

	b.clear ();
	b = a;

	if (b != a) {
		err_num++;
		std::cout << "operator = (sparse_power_series) failed.\n" << std::endl;
	}



	// *************************************************************
	// ************** reduce_last(), normalize() *******************
	// *************************************************************

	std::cout << " reduce_last(), normalize() ...\n" << std::endl;

	// set b = a, set coefficient zero with exponent last+5;
	// this increases the last exponent of b; undo this setting
	// by reducing the last-exponent explicitly

	b = a;
	b.set_coeff (zero, b.get_last() + 5);
	b.reduce_last (a.get_last ());
	b.normalize ();

	if (b != a) {
		err_num++;
		std::cout << "reduce_last() or normailze() failed.\n" << std::endl;
	}



	// *************************************************************
	// **************** assign_zero(), assign_one() **********************
	// *************************************************************

	std::cout << " assign_zero(), assign_one () ...\n" << std::endl;

	l = static_cast<lidia_size_t>(10); // must be non-negative
	b.assign_zero (l);

	if (b.get_first() != l || b.get_last() != l || b[l] != zero) {
		err_num++;
		std::cout << "assign_zero(lidia_size_t) failed.\n" << std::endl;
	}

	b.assign_one (l);

	if (b.get_first() == 0 && b.get_last() == l && b[0] == one) {
		for (i = 1; i <= l; i++) {
			if (b[i] != zero) {
				i = l;
				err_num++;
				std::cout << "assign_one(lidia_size_t) failed.\n" << std::endl;
			}
		}
	}
	else {
		err_num++;
		std::cout << "assign_one(lidia_size_t) failed.\n" << std::endl;
	}


	l = static_cast<lidia_size_t>(-10); // must be negative
	b.assign_one  (l);
	e.assign_zero (l);

	if (b != e) {
		err_num++;
		std::cout << "assign_one(lidia_size_t) failed.\n" << std::endl;
	}



	// *************************************************************
	// ********** arithmetic via friend - functions  ***************
	// *************************************************************

	std::cout << " arithmetic via friend - functions ...\n" << std::endl;

	tmp = static_cast<elem_type>(2);
	add (b, a, a);
	multiply (c, tmp, a);
	multiply (d, a, tmp);

	if (b != c) {
		err_num++;
		std::cout << " add(sparse_power_series, sparse_power_series, sparse_power_series) or ";
		std::cout << " multiply(sparse_power_series, T, sparse_power_series) (1) failed.\n" << std::endl;
	}

	if (b != d) {
		err_num++;
		std::cout << " add(sparse_power_series, sparse_power_series, sparse_power_series) or ";
		std::cout << " multiply(sparse_power_series, sparse_power_series, T) (2) failed.\n" << std::endl;
	}


	add (b, b, a);
	add (b, b, a);
	add (c, c, c);

	if (b != c) {
		err_num++;
		std::cout << " add(sparse_power_series, sparse_power_series, sparse_power_series) (1) failed\n " << std::endl;
	}



	f = a.get_first();
	l = a.get_last ();

	add (b, a, a);



	c = a;
	d = a;
	c.set_coeff (one, f-5);
	c.set_coeff (one, l+6);
	add (e, c, d);
	e.set_coeff (zero, f-5);

	if (e != b) {
		err_num++;
		std::cout << " add(sparse_power_series, sparse_power_series, sparse_power_series) (2) failed\n " << std::endl;
	}

	tmp = static_cast<elem_type>(2);
	c = a;
	d = a;
	c.set_coeff (one , f-5);
	c.set_coeff (zero, l);
	add (e, c, d);
	e.set_coeff (zero, f-5);
	e.set_coeff (tmp * a[l], l);

	if (e != b) {
		err_num++;
		std::cout << " add(sparse_power_series, sparse_power_series, sparse_power_series) (3) failed\n " << std::endl;
	}

	tmp = static_cast<elem_type>(2);
	c = a;
	d = a;
	c.set_coeff (zero, f);
	c.set_coeff (zero, l);
	add (e, c, d);
	e.set_coeff (tmp * a[f], f);
	e.set_coeff (tmp * a[l], l);

	if (e != b) {
		err_num++;
		std::cout << " add(sparse_power_series, sparse_power_series, sparse_power_series) (4) failed\n " << std::endl;
	}



	// b = 0 * x^a.last()
	subtract (b, a, a);

	if (b.get_first() != a.get_last() ||
	    b.get_last () != a.get_last() ||
	    b[b.get_first()] != zero) {
		err_num++;
		std::cout << " subtract(sparse_power_series, sparse_power_series, sparse_power_series) (1) failed\n " << std::endl;
	}


	// c = 2 * a, b = 4 * a
	add (c, a, a);
	add (b, c, c);

	subtract (d, b, c);

	if (d != c) {
		err_num++;
		std::cout << " subtract(sparse_power_series, sparse_power_series, sparse_power_series) (2) failed\n " << std::endl;
	}


	// b = 4*a + x^(a.last()+1)
	b.set_coeff (one, l+1);
	subtract (d, b, c);

	if (d != c) {
		err_num++;
		std::cout << " subtract(sparse_power_series, sparse_power_series, sparse_power_series) (3) failed\n " << std::endl;
	}

	// b = 4*a
	b.reduce_last (l);

	// b = x^(a.first()-5) + 4*a
	b.set_coeff (one, f-5);

	subtract (d, b, c);
	d.set_coeff (zero, f-5);

	if (d != c) {
		err_num++;
		std::cout << " subtract(sparse_power_series, sparse_power_series, sparse_power_series) (4) failed\n " << std::endl;
	}

	// b = 4*a
	b.set_coeff (zero, f-5);

	// b = 4*a - 4*a[f]*X^f - 4*a[l]*X^l, f = a.first(), l = a.last()
	b.set_coeff (zero, f);
	b.reduce_last (l-1);

	// d = b - 2*a
	subtract (d, b, c);

	// c = 2*a - 2*a[l]*X^l
	c.reduce_last (l-1);

	negate (tmp, c[f]);
	c.set_coeff (tmp, f);


	if (d != c) {
		err_num++;
		std::cout << " subtract(sparse_power_series, sparse_power_series, sparse_power_series) (5) failed\n " << std::endl;
	}


	add (b, a, a);
	subtract (b, b, a);

	if (b != a) {
		err_num++;
		std::cout << " subtract(sparse_power_series, sparse_power_series, sparse_power_series) (6) failed\n " << std::endl;
	}


	subtract (b, b, b);

	if (b.get_first() != a.get_last() ||
	    b.get_last () != a.get_last() ||
	    b[b.get_first()] != zero) {
		err_num++;
		std::cout << " subtract(sparse_power_series, sparse_power_series, sparse_power_series) (7) failed.\n" << std::endl;
	}




	//
	// ***** multiply, square, and power *****
	//

	b = a;
	b.set_coeff (one, l+5);
	b.set_coeff (one, f-6);

	multiply (c, a, b);
	multiply (d, b, a);

	if (c != d) {
		err_num++;
		std::cout << " multiply (sparse_power_series, sparse_power_series, sparse_power_series) failed.\n" << std::endl;
	}


	c = a;
	multiply (b, a, c);
	square   (c, a);

	if (c != b) {
		err_num++;
		std::cout << " multiply (sparse_power_series, sparse_power_series, sparse_power_series) or ";
		std::cout << " square (sparse_power_series, sparse_power_series) (1) failed\n " << std::endl;
	}


	multiply (b, b, c);
	square   (c, c);

	if (c != b) {
		err_num++;
		std::cout << " multiply (sparse_power_series, sparse_power_series, sparse_power_series) or ";
		std::cout << " square (sparse_power_series, sparse_power_series) (2) failed.\n" << std::endl;
	}


	multiply (b, b, a);
	power    (c, a, 5);

	if (c != b) {
		err_num++;
		std::cout << " multiply (sparse_power_series, sparse_power_series, sparse_power_series) or ";
		std::cout << " power (sparse_power_series, sparse_power_series, lidia_size_t) (1) failed.\n" << std::endl;
	}

	c = a;
	power (c, c, 5);

	if (c != b) {
		err_num++;
		std::cout << " multiply (sparse_power_series, sparse_power_series, sparse_power_series) or ";
		std::cout << " power (sparse_power_series, sparse_power_series, lidia_size_t) (2) failed.\n" << std::endl;
	}


	power (c, a, 1);

	if (c != a) {
		err_num++;
		std::cout << " power (sparse_power_series, sparse_power_series, lidia_size_t) (^1) failed.\n" << std::endl;
	}


	power (c, a, 0);

	if (c.get_first () != 0 ||
	    c.get_last  () != (a.get_last()-a.get_first()) ||
	    c[0] != one) {
		err_num++;
		std::cout << " power (sparse_power_series, sparse_power_series, lidia_size_t) (^0) failed.\n" << std::endl;
	}




	//
	//	***** multiply and negate *****
	//

	negate   (c, a);
	negate   (tmp, one);
	multiply (b, a, tmp);

	if (b != c) {
		err_num++;
		std::cout << " multiply (sparse_power_series, sparse_power_series, T) or ";
		std::cout << " negate (sparse_power_series, sparse_power_series) failed.\n" << std::endl;
	}

	c = a;
	negate   (c, c);

	if (b != c) {
		err_num++;
		std::cout << " multiply (sparse_power_series, sparse_power_series, T) or ";
		std::cout << " negate (sparse_power_series, sparse_power_series) failed.\n" << std::endl;
	}




	//
	// ***** invert and divide *****
	//

	invert (c, a);
	multiply (d, c, a);

	if (d.get_first () != 0 ||
	    d.get_last  () != (a.get_last()-a.get_first()) ||
	    d[0] != one) {
		err_num++;
		std::cout << " invert (sparse_power_series, sparse_power_series) failed.\n" << std::endl;
	}


	c = a;
	invert (c, c);
	multiply (d, c, a);

	if (d.get_first () != 0 ||
	    d.get_last  () != (a.get_last()-a.get_first()) ||
	    d[0] != one) {
		err_num++;
		std::cout << " invert (sparse_power_series, sparse_power_series) failed.\n" << std::endl;
	}


	square (c, c);
	b = c;
	square (c, c);
	multiply (c, b, c);
	invert (b, a);
	multiply (c, b, c); // now : c = (1/a)^7

	power (d, a, -7);

	if (c != d) {
		err_num++;
		std::cout << " power (sparse_power_series, sparse_power_series, lidia_size_t) failed.\n" << std::endl;
	}



	square (b, a);
	divide (c, b, a);

	if (c != a) {
		err_num++;
		std::cout << " divide (sparse_power_series, sparse_power_series, sparse_power_series) (1) failed.\n" << std::endl;
	}


	c = b;
	divide (c, c, a);

	if (c != a) {
		err_num++;
		std::cout << " divide (sparse_power_series, sparse_power_series, sparse_power_series) (2) failed.\n" << std::endl;
	}


	c = a;
	divide (c, b, c);

	if (c != a) {
		err_num++;
		std::cout << " divide (sparse_power_series, sparse_power_series, sparse_power_series) (3) failed.\n" << std::endl;
	}


	divide (c, c, c);

	if (c.get_first () != 0 ||
	    c.get_last  () != (a.get_last()-a.get_first()) ||
	    c[0] != one) {
		err_num++;
		std::cout << " divide (sparse_power_series, sparse_power_series, sparse_power_series) (4) failed.\n" << std::endl;
	}




	b = a;

	add (b, b, zero);

	if (b != a) {
		err_num++;
		std::cout << " add (sparse_power_series, sparse_power_series, T) failed.\n" << std::endl;
	}


	add (b, zero, b);

	if (b != a) {
		err_num++;
		std::cout << " add (sparse_power_series, T, sparse_power_series) failed.\n" << std::endl;
	}


	subtract (b, b, zero);

	if (b != a) {
		err_num++;
		std::cout << " subtract (sparse_power_series, sparse_power_series, T) failed.\n" << std::endl;
	}


	subtract (b, zero, b);
	negate (b, b);

	if (b != a) {
		err_num++;
		std::cout << " subtract (sparse_power_series, T, sparse_power_series) failed.\n" << std::endl;
	}


	c = a;
	// c.set_first ( 5 );
	c.multiply_by_xn (-c.get_first() + 5);

	add (b, c, one);
	add (b, b, one);
	c = b;
	subtract (b, c, one);
	subtract (b, b, one);

	// b.set_first ( a.get_first() );
	b.multiply_by_xn (-b.get_first() + a.get_first ());

	if (b != a) {
		err_num++;
		std::cout << " add (sparse_power_series, sparse_power_series, T) or ";
		std::cout << " subtract (sparse_power_series, sparse_power_series, T) failed.\n" << std::endl;
	}


	c = a;
	// c.set_first ( 5 );
	c.multiply_by_xn (-c.get_first() + 5);

	add (b, one, c);
	add (b, one, b);
	c = b;
	subtract (b, c, one);
	subtract (b, b, one);

	// b.set_first ( a.get_first() );
	b.multiply_by_xn (-b.get_first() + a.get_first ());

	if (b != a) {
		err_num++;
		std::cout << " add (sparse_power_series, T, sparse_power_series) or ";
		std::cout << " subtract (sparse_power_series, sparse_power_series, T) failed.\n" << std::endl;
	}

	c = a;
	// c.set_first ( 5 );
	c.multiply_by_xn (-c.get_first() + 5);

	subtract (b, one, c);
	subtract (b, one, b);

	// b.set_first ( a.get_first() );
	b.multiply_by_xn (-b.get_first() + a.get_first ());

	if (b != a) {
		err_num++;
		std::cout << " subtract (sparse_power_series, T, sparse_power_series) failed.\n" << std::endl;
	}



	tmp = static_cast<elem_type>(3);

	b = a;
	multiply (b, b, zero);

	multiply (c, a, tmp);
	subtract (d, c, a);
	subtract (d, d, a);

	e = a;
	multiply (e, e, tmp);

	if (d != a || e != c ||
	    !(b.get_first() == a.get_last() &&
	      b.get_last() == a.get_last() &&
	      b[a.get_last()] == zero)
		) {
		err_num++;
		std::cout << " multiply (sparse_power_series, sparse_power_series, T) failed.\n" << std::endl;
	}


	b = a;
	multiply (b, zero, a);

	multiply (d, tmp, a);

	e = a;
	multiply (e, tmp, e);

	if (d != c || e != c ||
	    !(b.get_first() == a.get_last() &&
	      b.get_last() == a.get_last() &&
	      b[a.get_last()] == zero)
		) {
		err_num++;
		std::cout << " multiply (sparse_power_series, T, sparse_power_series) failed.\n" << std::endl;
	}



	divide   (b, a, tmp);
	multiply (b, b, tmp);

	c = a;
	divide   (c, c, tmp);
	multiply (c, c, tmp);

	d = a;
	divide   (d, d, one);

	if (b != a || c != a || d != a) {
		err_num++;
		std::cout << " divide (sparse_power_series, sparse_power_series, T) failed.\n" << std::endl;
	}


	divide   (b, tmp, a);
	invert   (c, a);
	multiply (c, c, tmp);

	d = a;
	divide   (d, tmp, d);

	if (b != c || d != c) {
		err_num++;
		std::cout << " divide (sparse_power_series, T, sparse_power_series) failed.\n" << std::endl;
	}






	// *************************************************************
	// ***************** arithmetic via operators ******************
	// *************************************************************

	std::cout << " arithmetic via operators ...\n" << std::endl;

	tmp = static_cast<elem_type>(2);

	b = a + a;
	c = tmp * a;
	d = a * tmp;

	if (b != c) {
		err_num++;
		std::cout << " operator+(sparse_power_series, sparse_power_series, sparse_power_series) or ";
		std::cout << " operator*(sparse_power_series, T, sparse_power_series) (1) failed.\n" << std::endl;
	}

	if (b != d) {
		err_num++;
		std::cout << " operator+(sparse_power_series, sparse_power_series, sparse_power_series) or ";
		std::cout << " operator*(sparse_power_series, sparse_power_series, T) (2) failed.\n" << std::endl;
	}


	b = b + a;
	b = b + a;
	c = c + c;

	if (b != c) {
		err_num++;
		std::cout << " operator+(sparse_power_series, sparse_power_series, sparse_power_series) (1) failed\n " << std::endl;
	}


	f = a.get_first();
	l = a.get_last ();

	b = a + a;

	c = a;
	d = a;
	c.set_coeff (one, f-5);
	c.set_coeff (one, l+6);
	e = c + d;
	e.set_coeff (zero, f-5);

	if (e != b) {
		err_num++;
		std::cout << " operator+(sparse_power_series, sparse_power_series, sparse_power_series) (2) failed\n " << std::endl;
	}

	tmp = static_cast<elem_type>(2);
	c = a;
	d = a;
	c.set_coeff (one , f-5);
	c.set_coeff (zero, l);
	e = c + d;
	e.set_coeff (zero, f-5);
	e.set_coeff (tmp * a[l], l);

	if (e != b) {
		err_num++;
		std::cout << " operator+(sparse_power_series, sparse_power_series, sparse_power_series) (3) failed\n " << std::endl;
	}

	tmp = static_cast<elem_type>(2);
	c = a;
	d = a;
	c.set_coeff (zero, f);
	c.set_coeff (zero, l);
	e = c + d;
	e.set_coeff (tmp * a[f], f);
	e.set_coeff (tmp * a[l], l);

	if (e != b) {
		err_num++;
		std::cout << " operator+(sparse_power_series, sparse_power_series, sparse_power_series) (4) failed\n " << std::endl;
	}



	b = a - a;

	if (b.get_first() != a.get_last() ||
	    b.get_last () != a.get_last() ||
	    b[b.get_first()] != zero) {
		err_num++;
		std::cout << " operator-(sparse_power_series, sparse_power_series, sparse_power_series) (1) failed\n " << std::endl;
	}


	c = a + a;
	b = c + c;
	d = b - c;

	if (d != c) {
		err_num++;
		std::cout << " operator-(sparse_power_series, sparse_power_series, sparse_power_series) (2) failed\n " << std::endl;
	}


	b.set_coeff (one, l+1);
	d = b - c;

	if (d != c) {
		err_num++;
		std::cout << " operator-(sparse_power_series, sparse_power_series, sparse_power_series) (3) failed\n " << std::endl;
	}

	b.reduce_last (l);
	b.set_coeff (one, f-5);
	d = b - c;
	d.set_coeff (zero, f-5);

	if (d != c) {
		err_num++;
		std::cout << " operator-(sparse_power_series, sparse_power_series, sparse_power_series) (4) failed\n " << std::endl;
	}

	b.set_coeff (zero, f-5);
	b.set_coeff (zero, f);
	b.reduce_last (l-1);
	d = b - c;
	c.reduce_last (l-1);
	negate (tmp, c[f]);
	c.set_coeff (tmp, f);

	if (d != c) {
		err_num++;
		std::cout << " operator-(sparse_power_series, sparse_power_series, sparse_power_series) (5) failed\n " << std::endl;
	}


	b = a + a;
	b = b - a;

	if (b != a) {
		err_num++;
		std::cout << " operator-(sparse_power_series, sparse_power_series, sparse_power_series) (6) failed\n " << std::endl;
	}


	b = b - b;

	if (b.get_first() != a.get_last() ||
	    b.get_last () != a.get_last() ||
	    b[b.get_first()] != zero) {
		err_num++;
		std::cout << " operator-(sparse_power_series, sparse_power_series, sparse_power_series) (7) failed.\n" << std::endl;
	}




	b = a;
	b.set_coeff (one, l+5);
	b.set_coeff (one, f-6);

	c = a * b;
	d = b * a;

	if (c != d) {
		err_num++;
		std::cout << " operator* (sparse_power_series, sparse_power_series, sparse_power_series) (1) failed.\n" << std::endl;
	}


	c = a;
	b = a * c;
	square (c, a);

	if (c != b) {
		err_num++;
		std::cout << " operator* (sparse_power_series, sparse_power_series, sparse_power_series) (2) failed.\n" << std::endl;
	}


	b = b * c;
	square (c, c);

	if (c != b) {
		err_num++;
		std::cout << " operator* (sparse_power_series, sparse_power_series, sparse_power_series) (3) failed.\n" << std::endl;
	}


	b = b * a;
	power (c, a, 5);

	if (c != b) {
		err_num++;
		std::cout << " operator* (sparse_power_series, sparse_power_series, sparse_power_series) (4) failed.\n" << std::endl;
	}

	c = a;
	power (c, c, 5);

	if (c != b) {
		err_num++;
		std::cout << " operator* (sparse_power_series, sparse_power_series, sparse_power_series) (5) failed.\n" << std::endl;
	}




	square (b, a);
	c = b / a;

	if (c != a) {
		err_num++;
		std::cout << " operator/ (sparse_power_series, sparse_power_series, sparse_power_series) (1) failed.\n" << std::endl;
	}

	c = a;
	c = c / c;

	if (c.get_first () != 0 ||
	    c.get_last  () != (a.get_last()-a.get_first()) ||
	    c[0] != one) {
		err_num++;
		std::cout << " operator/ (sparse_power_series, sparse_power_series, sparse_power_series) (2) failed.\n" << std::endl;
	}




	b = a;

	b = b + zero;

	if (b != a) {
		err_num++;
		std::cout << " operator+ (sparse_power_series, sparse_power_series, T) failed.\n" << std::endl;
	}


	b = zero + b;

	if (b != a) {
		err_num++;
		std::cout << " operator+ (sparse_power_series, T, sparse_power_series) failed.\n" << std::endl;
	}


	b = b - zero;

	if (b != a) {
		err_num++;
		std::cout << " operator- (sparse_power_series, sparse_power_series, T) failed.\n" << std::endl;
	}


	b = zero - b;
	negate (b, b);

	if (b != a) {
		err_num++;
		std::cout << " operator- (sparse_power_series, T, sparse_power_series) failed.\n" << std::endl;
	}


	c = a;
	// c.set_first ( 5 );
	c.multiply_by_xn (-c.get_first() + 5);

	b = c + one;
	b = b + one;
	c = b;
	b = c - one;
	b = b - one;

	// b.set_first ( a.get_first() );
	b.multiply_by_xn (-b.get_first() + a.get_first ());

	if (b != a) {
		err_num++;
		std::cout << " operator+ (sparse_power_series, sparse_power_series, T) or ";
		std::cout << " operator- (sparse_power_series, sparse_power_series, T) failed.\n" << std::endl;
	}


	c = a;
	// c.set_first ( 5 );
	c.multiply_by_xn (-c.get_first() + 5);

	b = one + c;
	b = one + b;
	c = b;
	b = c - one;
	b = b - one;

	//  b.set_first ( a.get_first() );
	b.multiply_by_xn (-b.get_first() + a.get_first ());

	if (b != a) {
		err_num++;
		std::cout << " operator+ (sparse_power_series, T, sparse_power_series) or ";
		std::cout << " operator- (sparse_power_series, sparse_power_series, T) failed.\n" << std::endl;

		std::cout << a << std::endl;
		std::cout << b << std::endl;
	}

	c = a;
	// c.set_first ( 5 );
	c.multiply_by_xn (-c.get_first() + 5);

	b = one - c;
	b = one - b;

	// b.set_first ( a.get_first() );
	b.multiply_by_xn (-b.get_first() + a.get_first ());

	if (b != a) {
		err_num++;
		std::cout << " operator- (sparse_power_series, T, sparse_power_series) failed.\n" << std::endl;
	}




	tmp = static_cast<elem_type>(3);

	b = a;
	b = b * zero;

	c = a * tmp;
	subtract (d, c, a);
	subtract (d, d, a);

	if (d != a ||
	    !(b.get_first() == a.get_last() &&
	      b.get_last() == a.get_last() &&
	      b[a.get_last()] == zero)
		) {
		err_num++;
		std::cout << " operator* (sparse_power_series, sparse_power_series, T) failed.\n" << std::endl;
	}


	b = a;
	b = zero * a;

	d = tmp * a;

	if (d != c ||
	    !(b.get_first() == a.get_last() &&
	      b.get_last() == a.get_last() &&
	      b[a.get_last()] == zero)
		) {
		err_num++;
		std::cout << " operator* (sparse_power_series, T, sparse_power_series) failed.\n" << std::endl;
	}



	b = a / tmp;
	b = b * tmp;

	d = a / one;

	if (b != a || d != a) {
		err_num++;
		std::cout << " operator/ (sparse_power_series, sparse_power_series, T) failed.\n" << std::endl;
	}


	b = tmp / a;
	invert   (c, a);
	c = c * tmp;

	if (b != c) {
		err_num++;
		std::cout << " operator/ (sparse_power_series, T, sparse_power_series) failed.\n" << std::endl;
	}



	// *************************************************************
	// ********************* miscellaneous *************************
	// *************************************************************


	f = a.get_first ();
	l = a.get_last  ();


	b = a;
	b.multiply_by_xn (f);

	if (b.get_first() != 2*f) {
		err_num++;
		std::cout << " multiply_by_xn(lidia_size_t) failed.\n" << std::endl;
	}



	b = a;
	b.compose (5);

	if (b.get_first() != 5*f ||
	    b.get_last () != 5*l ||
	    b[b.get_last()] != a[a.get_last()]) {
		err_num++;
		std::cout << " compose(lidia_size_t) (1) failed.\n" << std::endl;
	}
	else {
		for (i = f; i < l; i++) {
			if (a[i] != b[5*i]) {
				i = l;
				err_num++;
				std::cout << " compose(lidia_size_t) (2) failed.\n" << std::endl;
			}
			else {
				k = 5*i;

				for (j = 1; j < 5; j++)

					if (b[k+j] != zero) {
						j = 5;
						i = l;
						err_num++;
						std::cout << " compose(lidia_size_t) (3) failed.\n" << std::endl;
					}
			}
		}
	}



	c = a;
	b = a;
	b.set_coeff (one, b.get_first() - 1);
	d = b;

	swap (c, d);

	if (c != b || d != a) {
		err_num++;
		std::cout << " swap(sparse_power_series, sparse_power_series) failed.\n" << std::endl;
	}



	// *************************************************************
	// ************************** summary **************************
	// *************************************************************


	if (err_num == 0) {
		std::cout << " No error detected ! :-)\n" << std::endl;
	}
	else {
		std::cout << " " << err_num << " error(s) detected ! :-(\n" << std::endl;
	}

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
