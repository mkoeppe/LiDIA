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
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/bigint.h"
#include        "LiDIA/isstream.h"
#include        "LiDIA/osstream.h"
#include	<cassert>
#include        <cctype>
#include        <cstdlib>
#include        <string>
#include        <vector>
#include        <algorithm>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



bigint factorial(bigint n)
{
	bigint f;
	if (n < 0)
		f = 0;
	else {
		f = 1;
		while (n > 0) {
			f *= n;
			--n;
		}
	}
	return f;
}



bigint fibonacci(long n)
{
	bigint f;
	if (n <= 0)
		f = 0;
	else {
		f = 1;
		bigint prev = 0;
		bigint tmp;
		while (n > 1) {
			tmp = f;
			f += prev;
			prev = tmp;
			--n;
		}
	}
	return f;
}



void identitytest(bigint& a, bigint& b, bigint& c)
{
	assert(-(-a) == a);
	assert((a + b) == (b + a));
	assert((a + (-b)) == (a - b));
	assert((a * b) == (b * a));
	assert((a * (-b)) == -(a * b));
	assert((a / (-b)) == -(a / b));
	assert((a - b) == -(b - a));
	assert((a + (b + c)) == ((a + b) + c));
	assert((a * (b * c)) == ((a * b) * c));
	assert((a * (b + c)) == ((a * b) + (a * c)));
	assert(((a - b) + b) == a);
	assert(((a + b) - b) == a);
	assert(((a * b) / b) == a);
	assert(((a * b) % b) == 0);
	assert((b * (a / b) + (a % b)) == a);
	assert(((a + b) % c) == ((a % c) + (b % c)) % c);
}



void utiltest(bigint& a)
{
	bigint b, c;
	square(b, a);
	sqrt(b, b);
	assert(b == a);

	sqrt(b, a);
	square(b, b);
	assert(b <= a);

	bigint x = 1;
	for (int i = 0; i < 10; ++i) {
		power(b, a, i);
		assert(b == x);
		x *= a;
	}

	bigint j;
	x = 1;
	for (j = 0; j < 10; ++j) {
		power(b, a, j);
		assert(b == x);
		x *= a;
	}

	x.assign_one();
	assert(x.is_one());
	assert(x.is_odd());
	assert(!x.is_even());
	x.assign(4);
	assert(x.is_even());
	assert(!x.is_odd());
	assert(x % 4 == 0);

	b = a;
	b.divide_by_2(); assert(b == (a/bigint(2)));

	b = a; b.negate();
	b.divide_by_2(); assert(b == ((-a)/bigint(2)));
}



void bittest(bigint& a, bigint& b, bigint& c)
{
	assert((a | b) == (b | a));
	assert((a & b) == (b & a));
	assert((a ^ b) == (b ^ a));
	assert((a | (b | c)) == ((a | b) | c));
	assert((a & (b & c)) == ((a & b) & c));
	assert((a & (b | c)) == ((a & b) | (a & c)));
	assert((a | (b & c)) == ((a | b) & (a | c)));
	assert((a & (a | b)) == a);
	assert((a | (a & b)) == a);
}



void accumtest(bigint& a, bigint& b, bigint& c)
{
	bigint x = a;
	x *= b;
	assert(x == (a * b));
	x += c;
	assert(x == ((a * b) + c));
	x -= a;
	assert(x == (((a * b) + c) - a));
	x /= b;
	assert(x == ((((a * b) + c) - a) / b));
	x %= c;
	assert(x == (((((a * b) + c) - a) / b) % c));
	x &= a;
	assert(x == ((((((a * b) + c) - a) / b) % c) & a));
	x |= b;
	assert(x == (((((((a * b) + c) - a) / b) % c) & a) | b));
	x ^= c;
	assert(x == ((((((((a * b) + c) - a) / b) % c) & a) | b) ^ c));
}



void longidentitytest(bigint& a, long b, long c)
{
	assert((a + b) == (b + a));
	assert((a + (-b)) == (a - b));
	assert((a * b) == (b * a));
	assert((a * (-b)) == -(a * b));
	assert((a / (-b)) == -(a / b));
	assert((a - b) == -(b - a));
	assert((a + (b + c)) == ((a + b) + c));
	assert((a * (b * c)) == ((a * b) * c));
	assert((a * (b + c)) == ((a * b) + (a * c)));
	assert(((a - b) + b) == a);
	assert(((a + b) - b) == a);
	assert(((a * b) / b) == a);
	assert(((a * b) % b) == 0);
	assert((b * (a / b) + (a % b)) == a);
	assert(((a + b) % c) == ((a % c) + (b % c)) % c);
}



void longbittest(bigint& a, long b, long c)
{
	assert((a | b) == (b | a));
	assert((a & b) == (b & a));
	assert((a ^ b) == (b ^ a));
	assert((a | (b | c)) == ((a | b) | c));
	assert((a & (b & c)) == ((a & b) & c));
	assert((a & (b | c)) == ((a & b) | (a & c)));
	assert((a | (b & c)) == ((a | b) & (a | c)));
	assert((a & (a | b)) == a);
	assert((a | (a & b)) == a);
}



void longaccumtest(bigint& a, long b, long c)
{
	bigint x = a;
	x *= b;
	assert(x == (a * b));
	x += c;
	assert(x == ((a * b) + c));
	x -= a;
	assert(x == (((a * b) + c) - a));
	x /= b;
	assert(x == ((((a * b) + c) - a) / b));
	x %= c;
	assert(x == (((((a * b) + c) - a) / b) % c));
	x &= a;
	assert(x == ((((((a * b) + c) - a) / b) % c) & a));
	x |= b;
	assert(x == (((((((a * b) + c) - a) / b) % c) & a) | b));
	x ^= c;
	assert(x == ((((((((a * b) + c) - a) / b) % c) & a) | b) ^ c));
}



void anothertest()
{
	bigint pow64;
	power(pow64, 2, 64);
	std::cout << "power(pow64, 2, 64) = " << pow64 << "\n";
	std::cout << "pow64.bit_length() = " << pow64.bit_length() << "\n";
	assert(pow64.bit_length() == 65);

	bigint s64 = 1;
	s64 <<= 64;
	std::cout << "s64 = 1 << 64 = " << s64 << "\n";

	assert(s64 == pow64);
	assert(s64 >= pow64);
	assert(s64 <= pow64);
	assert(!(s64 != pow64));
	assert(!(s64 > pow64));
	assert(!(s64 < pow64));

	bigint s32 = s64 >> 32;
	std::cout << "s32 = s64 >> 32 = " << s32 << "\n";
	assert(s32.bit_length() == 33);
	assert(!(pow64 == s32));
	assert(!(pow64 < s32));
	assert(!(pow64 <= s32));
	assert(pow64 != s32);
	assert(pow64 >= s32);
	assert(pow64 > s32);

	bigint comps64 = ~s64;
	std::cout << "comps64 = ~s64 = " << comps64 << "\n";

	bigint result = (comps64 & s32);
	std::cout << "comps64 & s32 = " << result << "\n";
	result = (comps64 | s32);
	std::cout << "comps64 | s32 = " << result << "\n";
	result = (comps64 ^ s32);
	std::cout << "comps64 ^ s32 = " << result << "\n";

	assert((s64 / s32) == s32);

	identitytest(s64, s32, comps64);
	bittest(s32, s64, comps64);
	accumtest(comps64, s32, pow64);
	utiltest(s32);
	longidentitytest(s64, 1000, 50);
	longbittest(s64, 12345, 67890);
	longaccumtest(s32, 100000, 1);
}



void iotest()
{
	bigint result;

	std::cout << "\nenter an Integer: ";
	std::cin >> result;
	std::cout << "number = " << result << "\n";
}



void fibtest()
{
	bigint fib50 = fibonacci(50);
	std::cout << "fib50 = fibonacci(50) = " << fib50 << "\n";
	bigint fib48 = fibonacci(48);
	std::cout << "fib48 = fibonacci(48) = " << fib48 << "\n";

	bigint result = fib48 + fib50;
	std::cout << "fib48 + fib50 = " << result << "\n";
	result = fib48 - fib50;
	std::cout << "fib48 - fib50 = " << result << "\n";
	result = fib48 * fib50;
	std::cout << "fib48 * fib50 = " << result << "\n";
	result = fib48 / fib50;
	std::cout << "fib48 / fib50 = " << result << "\n";
	result = fib48 % fib50;
	std::cout << "fib48 % fib50 = " << result << "\n";
	result = gcd(fib50, fib48);
	std::cout << "gcd(fib50, fib48) = " << result << "\n";
	sqrt(result, fib50);
	std::cout << "sqrt(fib50) = " << result << "\n";

	identitytest(result, fib50, fib48);
	bittest(result, fib50, fib48);
	accumtest(result, fib50, fib48);
	utiltest(fib48);
	longidentitytest(fib50, 1000, 50);
	longaccumtest(fib48, 100000, 1);
}



void facttest(bigint& one, bigint& two)
{
	bigint fact30(factorial(30));
	std::cout << "fact30 = factorial(30) = " << fact30 << "\n";

	bigint fact28(factorial(28));
	std::cout << "fact28 = factorial(28) = " << fact28 << "\n";
	assert(fact30 == fact28 * 870);

	bigint result = fact30 + fact28;
	std::cout << "fact30 + fact28 = " << result << "\n";
	result = fact30 - fact28;
	std::cout << "fact30 - fact28 = " << result << "\n";
	result = fact30 * fact28;
	std::cout << "fact30 * fact28 = " << result << "\n";
	result = fact30 / fact28;
	std::cout << "fact30 / fact28 = " << result << "\n";
	result = fact30 % fact28;
	std::cout << "fact30 % fact28 = " << result << "\n";

	result = -fact30;
	std::cout << "-fact30 = " << result << "\n";
	assert(abs(result) == fact30);

	std::cout << "lg(fact30) = " << fact30.bit_length() << "\n";
	assert(fact30.bit_length() == 108);

	result = gcd(fact30, fact28);
	std::cout << "gcd(fact30, fact28) = " << result << "\n";
	assert(result == fact28);

	sqrt(result, fact30);
	std::cout << "sqrt(fact30) = " << result << "\n";

	bigint negfact31 = fact30 * -31;
	bigint posfact31 = abs(negfact31);
	std::cout << "negfact31 = " << negfact31 << "\n";
	result = fact30 + negfact31;
	std::cout << "fact30 + negfact31 = " << result << "\n";
	result = fact30 - negfact31;
	std::cout << "fact30 - negfact31 = " << result << "\n";
	result = fact30 * negfact31;
	std::cout << "fact30 * negfact31 = " << result << "\n";
	result = fact30 / negfact31;
	std::cout << "fact30 / negfact31 = " << result << "\n";
	result = fact30 % negfact31;
	std::cout << "fact30 % negfact31 = " << result << "\n";
	result = gcd(fact30, negfact31);
	std::cout << "gcd(fact30, negfact31) = " << result << "\n";
	assert(result == fact30);

	identitytest(one, one, one);
	identitytest(one, one, one);
	identitytest(one, two, fact30);
	identitytest(fact30, posfact31, fact28);
	identitytest(fact30, negfact31, fact28);
	identitytest(negfact31, posfact31, fact28);

	bittest(one, one, one);
	bittest(one, one, one);
	bittest(one, two, fact30);
	bittest(fact30, posfact31, fact28);

	accumtest(one, one, one);
	accumtest(one, one, one);
	accumtest(one, two, fact30);
	accumtest(fact30, posfact31, fact28);

	utiltest(one);
	utiltest(fact30);
	utiltest(posfact31);

	longidentitytest(one, 1, 1);
	longidentitytest(one, 2, 3);
	longidentitytest(fact30, 3, -20);
	longidentitytest(fact30, 4, 20000);
	longidentitytest(negfact31, -100, 20000);

	longbittest(one, 1, 1);
	longbittest(one, 2, 3);
	longbittest(fact30, 4, 20000);
	longbittest(fact28, 1000, 50);

	longaccumtest(one, 1, 1);
	longaccumtest(one, 2, 3);
	longaccumtest(fact30, 4, 20000);
	longaccumtest(fact30, 1000, 50);
	longaccumtest(fact28, 10000000, 100000000);
}



void modtest()
{
	bigint b, e, m;

	m = 1; m <<= 32;
	b = m + 1;

	power(e, 2, 32);
	power(b, 2, 32); // use b as a comparison
	std::cout << "2^32 = " << e << "\n";

	e %= (e-1); // do same op two ways...
	b = b % (b - 1);

	std::cout << "2^32 % (2^32-1) = " << e << "\n"; // e is incorrect here
	std::cout << "2^32 % (2^32-1) = " << b << "\n"; // but b is ok
}



void misctest ()
{
	bigint a, c;

	// exponent is a long

	a.assign_zero(); // a = 0
	c.assign_one(); // just to be on the safe side
	power(c, a, 5L); assert(c.is_zero()); c.assign_one();
	power(c, a, -5L); assert(c.is_zero());
	power(c, a, 0L); assert(c.is_one());

	a.assign_one(); // a = 1
	c.assign_zero();
	power(c, a, 5L); assert(c.is_one()); c.assign_zero();
	power(c, a, -5L); assert(c.is_one()); c.assign_zero();
	power(c, a, 0L); assert(c.is_one()); c.assign_zero();

	a.negate(); // a = -1
	power(c, a, 5L); assert(c == a); c.assign_zero();
	power(c, a, -5L); assert(c == a); c.assign_zero();
	power(c, a, 0L); assert(c.is_one()); c.assign_zero();
	power(c, a, 4L); assert(c.is_one()); c.assign_zero();
	power(c, a, -4L); assert(c.is_one());

	add(a, a, a); // a = -2
	power(c, a, -5L); assert(c.is_zero());


	// exponent is a bigint

	a.assign_zero(); // a = 0
	c.assign_one(); // just to be on the safe side
	power(c, a, bigint(5)); assert(c.is_zero()); c.assign_one();
	power(c, a, bigint(-5)); assert(c.is_zero());
	power(c, a, bigint(0)); assert(c.is_one());

	a.assign_one(); // a = 1
	c.assign_zero();
	power(c, a, bigint(5)); assert(c.is_one()); c.assign_zero();
	power(c, a, bigint(-5)); assert(c.is_one()); c.assign_zero();
	power(c, a, bigint(0)); assert(c.is_one()); c.assign_zero();

	a.negate(); // a = -1
	power(c, a, bigint(5)); assert(c == a); c.assign_zero();
	power(c, a, bigint(-5)); assert(c == a); c.assign_zero();
	power(c, a, bigint(0)); assert(c.is_one()); c.assign_zero();
	power(c, a, bigint(4)); assert(c.is_one()); c.assign_zero();
	power(c, a, bigint(-4)); assert(c.is_one());

	add(a, a, a); // a = -2
	power(c, a, bigint(-5)); assert(c.is_zero());
}



void operator_test (const bigint & s)
{
	bigint a, b, c;

	a = s;
	b = s;
	c = ++a;
	inc(b);

	assert (c == b);
	assert (a == b);

	a = s;
	b = s;
	c = a++;
	inc(b);

	assert (c == s);
	assert (a == b);

	a = s;
	b = s;
	c = --a;
	dec(b);

	assert (c == b);
	assert (a == b);

	a = s;
	b = s;
	c = a--;
	dec(b);

	assert (c == s);
	assert (a == b);
}



void gcd_test (const bigint & A, const bigint & B)
{
	bigint g1, g2, u, v, a, b;

	a = A;
	b = B;

	// gcd
	g1 = gcd(a, b);
	g2 = gcd(b, a);

	if (g1 != 0) {
		assert(a % g1 == 0);
		assert(b % g1 == 0);
	}
	else
		assert(a == 0 && b == 0);

	assert(g1 == g2);
	assert(g1 >= 0);

	// bgcd
	g1 = bgcd(a, b);
	g2 = bgcd(b, a);

	if (g1 != 0) {
		assert(a % g1 == 0);
		assert(b % g1 == 0);
	}
	else
		assert(a == 0 && b == 0);

	assert(g1 == g2);
	assert(g1 >= 0);

	// dgcd
	g1 = dgcd(a, b);
	g2 = dgcd(b, a);

	if (g1 != 0) {
		assert(a % g1 == 0);
		assert(b % g1 == 0);
	}
	else
		assert(a == 0 && b == 0);

	assert(g1 == g2);
	assert(g1 >= 0);

	// xgcd
	g2 = xgcd (u, v, b, a);
	assert(g2 == (u*b+v*a));

	g1 = xgcd (u, v, a, b);
	assert(g1 == (u*a+v*b));

	if (g1 != 0) {
		assert(a % g1 == 0);
		assert(b % g1 == 0);
	}
	else
		assert(a == 0 && b == 0);

	assert(g1 == g2);
	assert(g1 >= 0);

	if (a == 0) {
		assert(u == 0);
		assert(v == b.sign());
	}
	else if (b == 0) {
		assert(u == a.sign());
		assert(v == 0);
	}
	else if (abs(a) == abs(b)) {
		assert(u == 0);
		assert(v == b.sign());
	}
	else {
		assert(2*abs(u)*g1 <= abs(b));
		assert(2*abs(v)*g1 <= abs(a));
	}
}



void constants_test ()
{
	bigint x, y;

	x = -10;
	y = 5;
	x /= y;
	assert (x == -2);

	x = -11;
	y = 5;
	x /= y;
	assert(x == -2);

	string_to_bigint("-238770661957910065368000", x);
	string_to_bigint("3270830985724795416000", y);
	x /= y;
	assert (x == -73);

	string_to_bigint("-20483", x);
	x %= 9;
	assert(x == -8);
}


void io_check(std::string const& s,
	      std::ios::fmtflags imask, std::ios::fmtflags omask) {
    typedef std::vector<std::string> StringVec;

    StringVec prefix;
    prefix.push_back("");
    prefix.push_back(" ");
    prefix.push_back(" \n ");
    prefix.push_back("\t ");
    prefix.push_back(" \\\n");

    StringVec suffix;
    suffix.push_back("");
    suffix.push_back(" ");
    suffix.push_back(" 1");
    suffix.push_back(" non-digit");
    suffix.push_back("+ ");
    suffix.push_back("GH ");
    suffix.push_back("\\\n");
    suffix.push_back("\\ ");
    suffix.push_back("\n");

    for(StringVec::size_type i = 0; i < prefix.size(); ++i) {
	for(StringVec::size_type j = 0; j < suffix.size(); ++j) {
	    isstream iss(prefix[i] + s + suffix[j]);
	    iss.unsetf(std::ios::basefield | std::ios::showbase);
	    iss.setf(imask);
	    
	    bigint x = 42;
	    iss >> x;
	    
	    osstream oss;
	    oss.unsetf(std::ios::basefield | std::ios::showbase);
	    oss.setf(omask);
	    
	    oss << x;
	    
	    std::string xs = extractString(oss);
	    for(std::string::iterator iter = xs.begin();
		iter != xs.end(); ++iter) {
		*iter = std::tolower(*iter);
	    }
	    
	    std::string ss = s;
	    for(std::string::iterator iter = ss.begin();
		iter != ss.end(); ++iter) {
		*iter = std::tolower(*iter);
	    }
	    
	    if(xs != ss) {
		std::cerr << "\n\nError while testing formatted bigint I/O:\n"
			  << "  input:    \""
			  << prefix[i] << s << suffix[j] << "\"\n"
			  << "  output:   \"" << xs << "\"\n"
			  << "  expected: \"" << ss << "\"\n\n";
		exit(EXIT_FAILURE);
	    }
	}
    }
}

	
void io_test() {
    // decimal
    std::string s0 = "0";
    std::string s1 = "123423453456456756786789";
    std::string s2 = "-123423453456456756786789";

    // octal
    std::string s3 = "0123423453456456756706701";
    std::string s4 = "-0123423453456456756706701";
    std::string s5 = "123423453456456756706701";
    std::string s6 = "-123423453456456756706701";

    // hexadecimal

    std::string s7 = "0x123423453456456756786789789a89ab9abcabcdbcdecdef";
    std::string s8 = "-0x123423453456456756786789789a89ab9abcabcdbcdecdef";
    std::string s9 = "123423453456456756786789789a89ab9abcabcdbcdecdef";
    std::string s10 = "-123423453456456756786789789a89ab9abcabcdbcdecdef";
    std::string s11 = "0x123423453456456756786789789A89ABb9abcabcdBCDEcdef";

    std::cout << "iostream basefield test ...  " << std::flush;

    std::ios::fmtflags oct  = std::ios::oct;
    std::ios::fmtflags dec  = std::ios::dec;
    std::ios::fmtflags hex  = std::ios::hex;
    std::ios::fmtflags showbase = std::ios::showbase;
    std::ios::fmtflags none = oct & hex;

    io_check(s0, none, none);
    io_check(s1, none, none);
    io_check(s2, none, none);
    io_check(s0, dec, none);
    io_check(s1, dec, none);
    io_check(s2, none, dec);
    io_check(s0, none, showbase);
    io_check(s1, none, showbase);
    io_check(s2, showbase, none);
    io_check(s0, dec | showbase, dec | showbase);
    io_check(s1, dec | showbase, dec | showbase);
    io_check(s2, dec | showbase, dec | showbase);

    io_check(s3, none, oct | showbase);
    io_check(s4, showbase, oct | showbase);
    io_check(s5, oct | showbase, oct);
    io_check(s6, oct, oct);

    io_check(s7, showbase, hex | showbase);
    io_check(s8, none, hex | showbase);
    io_check(s9, hex, hex);
    io_check(s10, hex | showbase, hex);
    io_check(s11, none, hex | showbase);

    std::cout << "succeeded" << std::endl;
}


int main_LiDIA(int argc, char** argv)
{
	bigint one = 1;
	std::cout << "one = " << one << "\n";
	assert(one == 1);
	std::cout << "one + 1 = " << (one + 1) << "\n";

	bigint two = 2;
	std::cout << "two = " << two << "\n";
	assert(two == 2);

	bigint n(0);
	assert(n.is_zero());

	bigint a = randomize (bigint(12345));
	bigint b = randomize (bigint(678));
	power (a, a, b);

	facttest(one, two);
	fibtest();
	anothertest();
	iotest();
	modtest();
	misctest();
	operator_test(a);
	constants_test();

	io_test();

//  gcd_test (0,0);
//  gcd_test (1,1);
//  gcd_test (1,-1);
//  gcd_test (-1,-1);
//  gcd_test (0,b+1);
//  gcd_test (0,-b-1);
//  gcd_test (a,b);
//  gcd_test (b,2*b);
//  gcd_test (a,3*a);

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
