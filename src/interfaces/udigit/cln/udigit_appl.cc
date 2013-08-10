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


#include	"LiDIA/udigit.h"
#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



void add_special_values ()
{
	udigit max_u;
	udigit a, b, c;
	udigit carry;

	max_u = max_udigit();

	a = 0; b = 0;
	carry = add (c, a, b);
	assert (c == 0);
	assert (carry == 0);

	a = max_u; b = 1;
	carry = add (c, a, b);
	assert (c == 0);
	assert (carry == 1);

	a = 1; b = max_u;
	carry = add (c, a, b);
	assert (c == 0);
	assert (carry == 1);

	if (max_u >= 2) {
		a = max_u; b = 2;
		carry = add (c, a, b);
		assert (c == 1);
		assert (carry == 1);

		a = 2; b = max_u;
		carry = add (c, a, b);
		assert (c == 1);
		assert (carry == 1);
	}
}



void subtract_special_values ()
{
	udigit max_u;
	udigit a, b, c;
	udigit carry;

	max_u = max_udigit();

	// carry = 0

	a = 0; b = 0;
	carry = subtract (c, a, b);
	assert (c == 0);
	assert (carry == 0);

	a = max_u; b = max_u;
	carry = subtract (c, a, b);
	assert (c == 0);
	assert (carry == 0);

	a = 0; b = max_u;
	carry = subtract (c, a, b);
	assert (c == 1);
	assert (carry == 1);

	if (max_u >= 2) {
		a = 1; b = max_u;
		carry = subtract (c, a, b);
		assert (c == 2);
		assert (carry == 1);
	}

	a = 0; b = 1;
	carry = subtract (c, a, b);
	assert (c == max_u);
	assert (carry == 1);

	a = 0; b = 0;
	carry = subtract (c, a, b);
	assert (c == 0);
	assert (carry == 0);


	// carry = 1

	a = max_u; b = max_u;
	carry = subtract (c, a, b, 1);
	assert (c == max_u);
	assert (carry == 1);

	a = 0; b = max_u;
	carry = subtract (c, a, b, 1);
	assert (c == 0);
	assert (carry == 1);

	a = 1; b = max_u;
	carry = subtract (c, a, b, 1);
	assert (c == 1);
	assert (carry == 1);

	a = 0; b = 1;
	carry = subtract (c, a, b, 1);
	assert (c == (max_u-1));
	assert (carry == 1);
}



void multiply_special_values ()
{
	udigit max_u;
	udigit a, b, c;
	udigit carry;

	max_u = max_udigit();

	a = 0; b = 0;
	carry = multiply (c, a, b);
	assert (c == 0);
	assert (carry == 0);

	a = 1; b = 0;
	carry = multiply (c, a, b);
	assert (c == 0);
	assert (carry == 0);

	a = 0; b = 1;
	carry = multiply (c, a, b);
	assert (c == 0);
	assert (carry == 0);

	a = 1; b = 1;
	carry = multiply (c, a, b);
	assert (c == 1);
	assert (carry == 0);

	a = max_u; b = 1;
	carry = multiply (c, a, b);
	assert (c == max_u);
	assert (carry == 0);
}



void add_mod_special_values ()
{
	udigit max_m = max_udigit_modulus();
	udigit a, b, c;

	a = max_m -1;
	b = 0;
	c = add_mod(a, b, max_m);
	assert (c == (max_m -1));

	a = max_m -1;
	b = 1;
	c = add_mod(a, b, max_m);
	assert (c == 0);

	if (max_m > 2) {
		a = max_m -1;
		b = 2;
		c = add_mod(a, b, max_m);
		assert (c == 1);

		a = max_m -1;
		b = max_m -1;
		c = add_mod(a, b, max_m);
		assert (c == (max_m-2));
	}
}



void subtract_mod_special_values ()
{
	udigit max_m = max_udigit_modulus();
	udigit a, b, c;

	a = max_m -1;
	b = 0;
	c = subtract_mod(a, b, max_m);
	assert (c == a);

	a = 0;
	b = max_m -1;
	c = subtract_mod(a, b, max_m);
	assert (c == 1);

	a = 0;
	b = 1;
	c = subtract_mod(a, b, max_m);
	assert (c == (max_m -1));

	if (max_m > 2) {
		a = 2;
		b = 1;
		c = subtract_mod(a, b, max_m);
		assert (c == 1);

		a = 1;
		b = 2;
		c = subtract_mod(a, b, max_m);
		assert (c == (max_m-1));
	}
}



void multiply_mod_special_values ()
{
	udigit max_m = max_udigit_modulus();
	udigit a, b, c;

	a = 0;
	b = 0;
	c = multiply_mod(a, b, max_m);
	assert (c == 0);

	a = 0;
	b = max_m -1;
	c = multiply_mod(a, b, max_m);
	assert (c == 0);

	a = 1;
	b = max_m -1;
	c = multiply_mod(a, b, max_m);
	assert (c == b);

	a = 1;
	b = 1;
	c = multiply_mod(a, b, max_m);
	assert (c == b);

	a = max_m -1;
	b = max_m -1;
	c = multiply_mod(a, b, max_m);
	assert (c == 1);
}



//
// taken from bigint_appl.c, identity_test()
//

void mod_identity_test_3param (udigit a, udigit b, udigit m)
{
	udigit c1, c2, c3, c4;

	//
	// -(-a) == a
	//
	c1 = negate_mod(a, m);
	c1 = negate_mod(c1, m);
	assert (c1 == a);

	//
	// (a + b) ==  (b + a)
	//
	c1 = add_mod(a, b, m);
	c2 = add_mod(b, a, m);
	assert (c1 == c2);

	//
	// (a + (-b)) == (a - b)
	//
	c1 = negate_mod(b, m);
	c1 = add_mod(a, c1, m);
	c2 = subtract_mod(a, b, m);
	assert (c1 == c2);

	//
	// (a * b) ==  (b * a)
	//
	c1 = multiply_mod(a, b, m);
	c2 = multiply_mod(b, a, m);
	assert (c1 == c2);

	//
	// (a * (-b)) == -(a * b)
	//
	c1 = negate_mod(b, m);
	c1 = multiply_mod(a, c1, m);
	c2 = multiply_mod(a, b, m);
	c2 = negate_mod(c2, m);
	assert (c1 == c2);

	//
	// (a - b) ==  -(b - a)
	//
	c1 = subtract_mod(a, b, m);
	c2 = subtract_mod(b, a, m);
	c2 = negate_mod(c2, m);
	assert (c1 == c2);

	//
	// ((a - b) + b)) ==  a
	//
	c1 = subtract_mod(a, b, m);
	c1 = add_mod(c1, b, m);
	assert (c1 == a);

	//
	// ((a + b) - b)) ==  a
	//
	c1 = add_mod(a, b, m);
	c1 = subtract_mod(c1, b, m);
	assert (c1 == a);

	//
	// (a+b)^2 == a^2 + 2ab + b^2 mod m
	//
	c1 = add_mod(a, b, m);
	c1 = multiply_mod(c1, c1, m);

	c2 = multiply_mod(a, b, m);
	c2 = add_mod(c2, c2, m);
	c3 = multiply_mod(a, a, m);
	c4 = multiply_mod(b, b, m);

	c2 = add_mod(c2, c3, m);
	c2 = add_mod(c2, c4, m);

	assert(c1 == c2);
}



void mod_identity_test_4param (udigit a, udigit b, udigit c, udigit m)
{
	udigit c1, c2, c3;

	//
	// (a + (b + c)) ==  ((a + b) + c)
	//
	c1 = add_mod(b, c, m);
	c1 = add_mod(a, c1, m);
	c2 = add_mod(a, b, m);
	c2 = add_mod(c2, c, m);
	assert (c1 == c2);

	//
	// (a * (b * c)) ==  ((a * b) * c)
	//
	c1 = multiply_mod(b, c, m);
	c1 = multiply_mod(a, c1, m);
	c2 = multiply_mod(a, b, m);
	c2 = multiply_mod(c2, c, m);
	assert (c1 == c2);

	//
	// (a * (b + c)) ==  ((a * b) + (a * c))
	//
	c1 = add_mod(b, c, m);
	c1 = multiply_mod(a, c1, m);
	c2 = multiply_mod(a, b, m);
	c3 = multiply_mod(a, c, m);
	c2 = add_mod(c2, c3, m);
	assert (c1 == c2);
}



int main_LiDIA(int argc, char** argv)
{
	udigit max_u;
	udigit max_m;
	udigit a, b, c;
	udigit delta;

	max_u = max_udigit();

	if (max_u == 0) {
		std::cout << "max_udigit() == 0->No tests started !\n";
		std::cout.flush();
	}
	else {
		add_special_values();
		subtract_special_values();
		multiply_special_values();
		add_mod_special_values();
		subtract_mod_special_values();
		multiply_mod_special_values();

		max_m = max_udigit_modulus();
		delta = max_u / 5;

		if (delta == 0)
			delta = 1;

		for (a = 0; a < max_m; a += delta)
			for (b = 0; b < max_m; b += delta) {
				mod_identity_test_3param(a, b, max_m);

				for (c = 0; c < max_m; c += delta)
					mod_identity_test_4param(a, b, c, max_m);
			}

		std::cout << "Tests passed.\n";
		std::cout.flush();
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
