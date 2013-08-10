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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/multi_bigmod.h"
#include	"LiDIA/random_generator.h"
#include	<cassert>
#include	<cstdio>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



void
identitytest (multi_bigmod & a, multi_bigmod & b, multi_bigmod & c)
{
	assert(-(-a) == a);
	assert((a + b) == (b + a));
	assert((a + (-b)) == (a - b));
	assert((a + (-b)) == (b - a + 2*a - 2*b));
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

	multi_bigmod d, d2;

	negate(d, a);
	assert(d == (-a));
	add(d, b, a);
	assert((a + b) == d);
	subtract(d, a , b);
	assert((a + (-b)) == d);
	negate(d, b);
	assert((a + d) == (b - a+2*a - 2*b));
	d.assign(a); d.multiply_by_2();
	assert((a + a) == d);
	multiply(d, a, b);
	assert(d == (b * a));
	square(d, a); multiply(d2, a , a);
	assert(d == d2);
	d.assign(b); d.negate(); multiply(d, a, d);
	assert(d == -(a * b));
	divide(d, a, b);
	assert((a / (-b)) == -d);
	add(d, b, c); multiply(d, a, d);
	assert(d == ((a * b) + (a * c)));
	d = b;
	d.invert();
	assert(((a / b) / a) == d);
}



void
utiltest (multi_bigmod & a)
{
	multi_bigmod b, c;
	multi_bigmod x, y;

	x.set_modulus (a.modulus());
	x = 1;
	y.set_modulus (a.modulus());
	y = 1;

	int i;
	for (i = 0; i < 10; ++i) {
		power(b, a, i);
		assert(b == x);
		x *= a;
		y = a * y;
		assert(y == x);
	}
	x.assign_one();
	assert(x.is_one());

	invert(c, a);

	for (i = 0; i > -10; --i) {
		power(b, a, i);
		assert(b == x);
		x *= c;
	}
}



void
accumtest (multi_bigmod & a, multi_bigmod & b, multi_bigmod & c)
{
	multi_bigmod x = a;
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



void
longidentitytest (multi_bigmod& a, long b, long c)
{
	assert((a + b) == (b + a));
	assert((a + (-b)) == (a - b));
	assert((a * b) == (b * a));
	assert((a * (-b)) == -(a * b));
	assert((a / (-b)) == -(a / b));
	assert((a - b) == -(b - a));
	assert((a + (b + c)) == ((a + b) + c));
	assert((a + (b + c)) == ((c + a) + b));
	assert((a * (b + c)) == ((a * b) + (a * c)));
	assert(((a - b) + b) == a);
	assert(((a + b) - b) == a);
	assert(((a * b) / b) == a);
	assert((b * (a / b)) == a);

	multi_bigmod d;

	negate(d, a);
	assert(-(d) == a);
	add(d, a, b);
	assert((a + b) == d);
	subtract(d, a , b);
	assert((a + (-b)) == d);
	negate(d, multi_bigmod(b, d.modulus()));
	assert((a + d) == (b - a+2*a - 2*b));
	d.assign(a); d.multiply_by_2();
	assert((a + a) == d);
	multiply(d, a, b);
	assert(d == (b * a));
	d.set_mantissa(b); d.negate(); multiply(d, a, d);
	assert(d == -(a * b));
	divide(d, a, b);
	assert((a / (-b)) == -d);
	add(d, multi_bigmod(c, a.modulus()), b); multiply(d, a, d);
	assert(d == ((a * b) + (a * c)));
	add(d, multi_bigmod(b, a.modulus()), c); multiply(d, a, d);
	assert(d == ((a * b) + (a * c)));

	bigint h;

	d = b; h = d.invert(1);
	assert((multi_bigmod(1, d.modulus()) / multi_bigmod(b, d.modulus())) == d);
	assert(((a / b) / a) == d);
}



void
longaccumtest (multi_bigmod & a, long b, long c)
{
	multi_bigmod x = a;
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



void
anothertest (bigint m)
{
	multi_bigmod pow64;
	pow64.set_modulus(m);

	multi_bigmod two(2, m);
	power(pow64, two, 64);

	multi_bigmod s64(1, m);

	power(s64, two, 65);
	divide(s64, s64, 2);

	assert(s64 == pow64);
	assert(!(s64 != pow64));

	multi_bigmod s32;
	power(s32, two, 32);
	assert(!(pow64 == s32));
	assert(pow64 != s32);

	identitytest(s64, s32, pow64);
	accumtest(pow64, s32, pow64);
	utiltest(s32);
	longidentitytest(s64, 1000, 50);
	longaccumtest(s32, 100000, 11);
}



void
iotest ()
{
	multi_bigmod result;

	std::cout << "\n please, enter an multi_bigmod (format: (mantissa, modulus)): ";
	std::cin >> result;
	std::cout << "number = ->" << result << " < -\n";

	char * string = new char[1000];

	multi_bigmod_to_string(result, string);
	std::cout << "\n multi_bigmod_to_string:->" << string << " < - \n";

	result.assign_zero();
	result.set_modulus(bigint(1));

	string_to_multi_bigmod(string, result);
	std::cout << " string_to_multi_bigmod:->" << result << " < -  \n\n";
	std::cout.flush();

	FILE *fp;

	fp = fopen("test.txt", "w");

	int i;
	for (i = 0; i < 10; i++) {
		inc(result);
		result.print_to_file(fp);
	}
	fclose (fp);

	fp = fopen("test.txt", "r");

	while (!feof(fp)) {
		result.scan_from_file(fp);
	}
	fclose(fp);
	std::remove("test.txt");
}



void
bigmod_and_multi_bigmod (bigint m)
{
	bigmod::set_modulus(m);
	bigint m1, m2;

	bigmod a, b;
	multi_bigmod c, d;

	m1 = randomize(m);
	m2 = randomize(m);

	a.assign(m1);
	b.assign(m2);

	c.set_modulus(m);
	c.set_mantissa(m1);

	d.assign(m2, m);

	assert(a == c);
	assert(b == d);
	assert(a+d == b+c);

	assert(b-a == b-c);
	assert(b*a == b*c);
	assert(b/a == b/c);

	assert(a-b == c-b);
	assert(a*b == c*b);
	assert(a/b == c/b);

	d = c = a;
	assert(a == d);

	d += a;
	c += c;
	assert(c == d);

	d = c = a;
	d -= a;
	c -= a;
	assert(c == d);

	d = c = a;
	d *= a;
	c *= c;
	assert(c == d);

	d = c = a;
	d /= a;
	c /= c;
	assert(c == d);
}



void
multi_bigmod_moduli (bigint m)
{
	multi_bigmod a(m, m);
	multi_bigmod b(m+bigint(1), m);

	multi_bigmod c(bigint(1), m+bigint(1));
	multi_bigmod d(bigint(1), m+bigint(1));

	multi_bigmod g(m, m);

	multi_bigmod e(bigint(2), m+bigint(2));
	multi_bigmod f(bigint(2), m+bigint(2));

	multi_bigmod h(m, m);
	multi_bigmod i;


	assert(a.is_zero());
	assert(g.is_zero());
	assert(h.is_zero());
	assert(b.is_one());
	assert(c.is_one());
	assert(d.is_one());

	assert(a == g);
	assert(c == d);
	assert(e == f);

	assert(a+b+multi_bigmod(bigint(2), m) == multi_bigmod(bigint(3), m));
	assert(c+d == c + multi_bigmod(bigint(1), m+bigint(1)));
	assert(e+f == f + multi_bigmod(bigint(2), m+bigint(2)));

	i.set_modulus(h);
	i.set_mantissa(h.mantissa());
	assert(h == i);
}



void
all_test (bigint m)
{
	random_generator rg;
	multi_bigmod a, b;
	long x, y;

	a.randomize(multi_bigmod(m - multi_bigmod(1, m)));
	b.randomize(a);

	rg >> x >> y;
	x %= 10000000;
	y %= 1000000;

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
	bigmod_and_multi_bigmod(m);
	multi_bigmod_moduli(m);
}



int main_LiDIA (int, char**)
{
	bigint modul;

	bigint::seed(234534);

	power(modul, bigint(10), 20);

	modul = randomize(modul);

	while (!is_prime(modul, 5)) {
		inc(modul);
	}

	std::cout << "\n Modul " << modul << " (prime number) \n";

	multi_bigmod one (1, modul);
	multi_bigmod two (2, modul);
	multi_bigmod n	 (0, modul);
	assert(n.is_zero());

	assert(one == multi_bigmod(1, modul));
	assert(two == multi_bigmod((1+1), modul));

	assert(n.is_zero());
	anothertest(modul);
	iotest();

	int j;
	for (j = 0; j < 10; j++) {
		all_test(modul);
		std::cout << "\n " << (j+1) << "th test succeeded ";
		std::cout.flush();
	}


	std::cout << "\nEnd of multi_bigmod_appl\n";
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
