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
#include	"LiDIA/bigmod.h"
#include	"LiDIA/random_generator.h"
#include	<cassert>
#include	<cstdio>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



void
identitytest (bigmod & a, bigmod & b, bigmod & c)
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

	bigmod d, d2;

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
utiltest (bigmod & a)
{
	bigmod b, c;

	bigmod x = 1, y = 1;

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
accumtest (bigmod & a, bigmod & b, bigmod & c)
{
	bigmod x = a;
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
longidentitytest (bigmod& a, long b, long c)
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

	bigmod d;

	negate(d, a);
	assert(-(d) == a);
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
	d.assign(b); d.negate(); multiply(d, a, d);
	assert(d == -(a * b));
	divide(d, a, b);
	assert((a / (-b)) == -d);
	add(d, b, bigmod(c)); multiply(d, a, d);
	assert(d == ((a * b) + (a * c)));
	add(d, bigmod(b), c); multiply(d, a, d);
	assert(d == ((a * b) + (a * c)));

	bigint h;

	d = b; h = d.invert(1);
	assert((1 / bigmod(b)) == d);
	assert(((a / b) / a) == d);
}



void
longaccumtest (bigmod & a, long b, long c)
{
	bigmod x = a;
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
anothertest ()
{
	bigmod pow64;
	power(pow64, 2, 64);

	bigmod s64 = 1;
	power(s64, 2, 65);
	divide(s64, s64, 2);

	assert(s64 == pow64);
	assert(!(s64 != pow64));

	bigmod s32;
	power(s32, 2, 32);
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
	bigmod result;

	std::cout << "\nenter an bigmod (format : mantissa): ";
	std::cin >> result;
	std::cout << "number = ->" << result << " < -\n";

	char * string = new char[1000];

	bigmod_to_string(result, string);
	std::cout << "\n bigmod_to_string:->" << string << " < - \n";

	string[0] = '7';
	string_to_bigmod(string, result);

	std::cout << "\n changed string[0] = '7'  -->";
	std::cout << " string_to_bigmod = ->" << result << " < -  \n\n";
	std::cout.flush();

	FILE *fp;

	fp = fopen("test.txt", "w");

	int i;
	for (i = 0; i < 10; i++) {
		inc(result);
		result.print_to_file(fp);
		fprintf(fp, " ");
	}
	fclose (fp);

	fp = fopen("test.txt", "r");

	while (!feof(fp)) {
		result.scan_from_file(fp);
	}
	fclose(fp);
	std::remove("test.txt");
	delete[] string;
}



void
all_test ()
{
	random_generator rg;
	bigmod a, b;
	long x, y;

	a = randomize(bigmod(bigmod::modulus() - bigmod(1)));
	b = randomize(a);

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

	bigmod::set_modulus(modul);

	bigmod one = 1;
	bigmod two = 2;
	bigmod n(0);
	assert(n.is_zero());

	assert(one == 1);
	assert(two == (1+1));

	assert(n.is_zero());
	anothertest();
	iotest();

	int j;
	for (j = 0; j < 10; j++) {
		all_test();
		std::cout << "\n " << (j+1) << "th test succeeded ";
		std::cout.flush();
	}


	std::cout << "\nEnd of bigmod_appl\n";
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
