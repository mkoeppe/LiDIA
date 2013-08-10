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
//	Author	: Volker M"uller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigrational.h"
#include	"LiDIA/random_generator.h"
#include	<cassert>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



void identitytest(bigrational& a, bigrational& b, bigrational& c)
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

	bigrational d;

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
	d = b; d.invert(); d.assign(b.numerator(), b.denominator()); assert(b == d);
}



void utiltest(bigrational& a)
{
	bigrational b, c;

	square(b, a);
	multiply(c, a, a);
	assert(b == c);

	bigrational x = 1, y = 1;

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
	assert(x.is_positive());
	x.negate();
	assert(x.is_negative());
}



void accumtest(bigrational& a, bigrational& b, bigrational& c)
{
	bigrational x = a;
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



void longidentitytest(bigrational& a, long b, long c)
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

	bigrational d;

	negate(d, a); assert(-(d) == a);
	add(d, b, a); assert((a + b) == d);
	subtract(d, a , b); assert((a + (-b)) == d);
	negate(d, b); assert((a + d) == (b - a+2*a - 2*b));
	d.assign(a); d.multiply_by_2(); assert((a + a) == d);
	multiply(d, a, b); assert(d == (b * a));
	d.assign(b); d.negate(); multiply(d, a, d); assert(d == -(a * b));
	divide(d, a, b); assert((a / (-b)) == -d);
	add(d, b, bigrational(c)); multiply(d, a, d);
	assert(d == ((a * b) + (a * c)));
	add(d, bigrational(b), c); multiply(d, a, d);
	assert(d == ((a * b) + (a * c)));
	d = b; d.invert(); assert(((a / b) / a) == d);
}



void longaccumtest(bigrational& a, long b, long c)
{
	bigrational x = a;
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
	bigrational pow64;
	power(pow64, 2, 64);
	std::cout << "power(pow64, 2, 64) = " << pow64 << "\n";

	bigrational s64 = 1;
	s64 <<= 64;
	std::cout << "s64 = 1 << 64 = " << s64 << "\n";

	assert(s64 == pow64);
	assert(s64 >= pow64);
	assert(s64 <= pow64);
	assert(!(s64 != pow64));
	assert(!(s64 > pow64));
	assert(!(s64 < pow64));

	bigrational s32 = s64 >> 32;
	std::cout << "s32 = (s64 >> 32) = " << s32 << "\n";
	assert(!(pow64 == s32));
	assert(!(pow64 < s32));
	assert(!(pow64 <= s32));
	assert(pow64 != s32);
	assert(pow64 >= s32);
	assert(pow64 > s32);

	identitytest(s64, s32, pow64);
	accumtest(pow64, s32, pow64);
	utiltest(s32);
	longidentitytest(s64, 1000, 50);
	longaccumtest(s32, 100000, 11);
}



void iotest()
{
	bigrational result;

	std::cout << "\nenter an bigrational: ";
	std::cin >> result;
	std::cout << "number = " << result << "\n";

// hier noch string-Fkt einfuegen

}



void all_test(const bigint & ii)
{
	bigrational a, b;
	random_generator rg;
	static long x = 1, y = 27;
	bigint aa, bb;

	aa = randomize(ii);
	bb = randomize(ii);
	divide(a, aa, bb);
	aa = randomize(ii);
	bb = randomize(ii);
	divide(b, aa, bb);

	random_generator::seed(y);
	rg >> x >> y;
	x &= 0x7fff;
	y &= 0x7fff;

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


void div_test() {
        int a_i(2);
        unsigned a_u(2);
        long a_l(2);
        unsigned long a_ul(2);
        bigint a_bi(2);
        bigrational a_r(2, 1);

        bigrational b(2, 3);

        bigrational result(3, 1);

        assert(a_i / b == result);
        assert(a_u / b == result);
        assert(a_l / b == result);
        assert(a_ul / b == result);
        assert(a_bi / b == result);
        assert(a_r / b == result);

        a_i = 0;
        a_u = 0;
        a_l = 0;
        a_ul = 0;
        a_bi.assign_zero();
        a_r.assign_zero();

        result.assign_zero();

        assert(a_i / b == result);
        assert(a_u / b == result);
        assert(a_l / b == result);
        assert(a_ul / b == result);
        assert(a_bi / b == result);
        assert(a_r / b == result);
}


int main_LiDIA(int argc, char** argv)
{
	bigrational one = 1;
	std::cout << "one = " << one << "\n";
	assert(one == 1);
	std::cout << "one + 1 = " << (one + 1) << "\n";

	bigrational two = 2;
	std::cout << "two = " << two << "\n";
	assert(two == 2);

	bigrational n(0);
	assert(n.is_zero());
	anothertest();
	iotest();

	bigint i = 129099121;

	for (int j = 0; j < 10; j++)
		all_test(i);

        div_test();

	std::cout << "\nEnd of test\n";
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
