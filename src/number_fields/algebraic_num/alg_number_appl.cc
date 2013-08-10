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
#include	<cassert>
#include	<cstdlib>
#include	<fstream>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



std::ifstream in("alg_appl_input");



void identitytest(const alg_number& a, const alg_number& b, const alg_number& c)
{
	assert(-(-a) == a);
	assert((a + b) == (b + a));
	assert((a + (-b)) == (a - b));
	assert((a + (-b)) == (b - a+2*a - 2*b));
	assert((a * b) == (b * a));
	assert((a * (-b)) == -(a * b));
	assert((a / (-b)) == -(a / b));
	assert((a - b) == -(b - a));
	assert((a + (b + c)) == ((a + b) + c));
	assert((a + (b + c)) == ((a + c) + b));
	assert((a * (b * c)) == ((a * b) * c));
	assert((a * (b * c)) == ((c * a) * b));
	assert((a * (b + c)) == ((a * b) + (a * c)));
	assert(((a - b) + b) == a);
	assert(((a + b) - b) == a);
	assert(((a * b) / b) == a);
	assert(((a * b) / a) == b);
	assert(((a / b) / a) == (1/b));

	alg_number d;

	negate(d, a); assert(-(d) == a);
	add(d, b, a); assert((a + b) == d);
	subtract(d, a , b); assert((a + (-b)) == d);
	negate(d, b); assert((a + d) == (b - a+2*a - 2*b));
	d.assign(a); d.multiply_by_2(); assert((a + a) == d);
	multiply(d, a, b); assert(d == (b * a));
	d.assign(b); d.negate(); multiply(d, a, d); assert(d == -(a * b));
	divide(d, a, b); assert((a / (-b)) == -d);
	add(d, b, c); multiply(d, a, d); assert(d == ((a * b) + (a * c)));
	d = b; d.invert(); assert(((a / b) / a) == d);
}



void utiltest(const alg_number& a)
{
	alg_number b, c;

	square(b, a);
	multiply(c, a, a);
	assert(b == c);

	alg_number x = bigint(1), y = bigint(1);

	for (int i = 0; i < 10; ++i) {
		power(b, a, i);
		assert(b == x);
		x *= a;
		y = a * y;
		assert(y == x);
	}
	x.assign_one();
	assert(x.is_one());
	x.assign(4);
	x.negate();
	x.negate();
	assert(x == bigint(4));
}



void accumtest(const alg_number& a, const alg_number& b, const alg_number& c)
{
	alg_number x = a;
	x *= b;
	assert(x == (b * a));
	x += c;
	assert(x == (c+(b * a)));
	x -= a;
	assert(x == (-a + c + (b * a)));
	x /= b;
	assert(x == (((b *a) -a + c) / b));

	x.assign(a);
	assert(x == a);
	multiply(x, x, b);
	assert(x == (b * a));
	add(x, x, c);
	assert(x == (c+(b * a)));
	subtract(x, x, a);
	assert(x == (-a + c + (b * a)));
	divide(x, x, b);
	assert(x == (((b *a) -a + c) / b));
}



void longidentitytest(const alg_number& a, long b, long c)
{
	assert((a + b) == (b + a));
	assert((a + (-b)) == (a - b));
	assert((a * b) == (b * a));
	assert((a * (-b)) == -(a * b));
	assert((a / (-b)) == -(a / b));
	assert((a - b) == -(b - a));
	assert((a + (b + c)) == ((a + b) + c));
	assert((a + (b + c)) == ((c + a) + b));
	assert((a * (b * c)) == ((b * a) * c));
	assert((a * (b + c)) == ((a * b) + (a * c)));
	assert(((a - b) + b) == a);
	assert(((a + b) - b) == a);
	assert(((a * b) / b) == a);
	assert((b * (a / b)) == a);

	alg_number d;

	negate(d, a); assert(-(d) == a);
	add(d, b, a); assert((a + b) == d);
	subtract(d, a , b); assert((a + (-b)) == d);
	negate(d, bigint(b)); assert((a + d) == (b - a+2*a - 2*b));
	d.assign(a); d.multiply_by_2(); assert((a + a) == d);
	multiply(d, a, b); assert(d == (b * a));
	d.assign(b); d.negate(); multiply(d, a, d); assert(d == -(a * b));
	divide(d, a, b); assert((a / (-b)) == -d);
	add(d, alg_number(bigint(b)), c); multiply(d, a, d);
	assert(d == ((a * b) + (a * c)));
	add(d, b, alg_number(bigint(c))); multiply(d, a, d);
	assert(d == ((a * b) + (a * c)));
	d = bigint(b); d.invert(); assert(((a / b) / a) == d);
}



void longaccumtest(const alg_number& a, long b, long c)
{
	alg_number x = a;
	x *= b;
	assert(x == (a * b));
	x += c;
	assert(x == ((a * b) + c));
	x -= a;
	assert(x == (((a * b) + c) - a));
	x /= b;
	assert(x == ((((a * b) + c) - a) / b));

	x.assign(a);
	assert(x == a);
	multiply(x, x, b);
	assert(x == (a * b));
	add(x, x, c);
	assert(x == ((a * b) + c));
	subtract(x, x, a);
	assert(x == (((a * b) + c) - a));
	divide(x, x, b);
	assert(x == ((((a * b) + c) - a) / b));
}



void anothertest()
{
	alg_number pow64;
	power(pow64, alg_number(bigint(2)), 64);
	std::cout << "power(pow64, 2, 64) = " << pow64 << "\n";

	bigint s64 = 1;
	s64 <<= 64;
	std::cout << "s64 = 1 << 64 = " << s64 << "\n";

	assert(s64 == pow64);
	assert(!(s64 != pow64));

	bigint s32 = s64 >> 32;
	std::cout << "s32 = (s64 >> 32) = " << s32 << "\n";
	assert(!(pow64 == s32));
	assert(pow64 != s32);

//  identitytest(s64, s32, pow64);
	accumtest(pow64, s32, pow64);
	utiltest(s32);
	longidentitytest(s64, 1000, 50);
	longaccumtest(s32, 100000, 11);
}



void iotest()
{
	alg_number result;

	std::cout << "\nenter an alg_number: ";
	in >> result;
	std::cout << "number = " << result << "\n";
}



void all_test(const bigint & ii)
{
	alg_number a, b;
	static long x = 1, y = 27;
	alg_number aa, bb;

	aa.randomize(ii);
	bb.randomize(ii);
	divide(a, alg_number(aa), bb);
	aa.randomize(ii);
	bb.randomize(ii);
	divide(b, alg_number(aa), bb);

	srand(static_cast<unsigned int>(x)); srand(static_cast<unsigned int>(y));
	x = rand() & 0x7fff; y = rand() & 0x7fff;

	std::cout << "\n\n Test with \na = " << a << "\nb = " << b;
	std::cout << "\nx = " << x << "\ny = " << y << std::flush;

	utiltest(a);
	identitytest(a, b, b);
	identitytest(a, a, b);

	accumtest(a, a, b);
	accumtest(b, a, b);

	utiltest(a);
	utiltest(b);

	longidentitytest(a, x, y);
	longidentitytest(b, x, x);

	longaccumtest(b, x, y);
	longaccumtest(a, x, x);
}



int main_LiDIA(int argc, char** argv)
{
	number_field F;
	order O;

	in >> O;
	std::cout << "Read order !" << std::endl;
	std::cout << O;
	for (long zaehl = 1; zaehl <= 2; zaehl++) {
		alg_number one = bigint(1);
		std::cout << "one = " << one << "\n";
		assert(one == bigint(1));
		std::cout << "one + 1 = " << (one + 1) << "\n";

		alg_number two = bigint(2);
		std::cout << "two = " << two << "\n";
		assert(two == bigint(2));

		std::cout << std::endl;

		alg_number n(bigint(0));
		assert(n.is_zero());
		anothertest();
		iotest();

		bigint i = 129099121;

		for (int j = 0; j < 10; j++)
			all_test(i);
		std::cout << "\nEnd of test\n";
		if (zaehl != 2) {
			in >> F;
			O.assign(order(F));
			std::cout << "Starting second iteration!!\n";
			std::cout << "===========================\n";
		}
		else std::cout << std::endl;
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
