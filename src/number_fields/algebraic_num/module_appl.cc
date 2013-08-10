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


#include	<cassert>
#include	<cstdlib>
#include	<fstream>

#include	"LiDIA/alg_number.h"
#include	"LiDIA/timer.h"



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



std::ifstream in("alg_appl_input");



void identitytest(module a, module b, module c) // MM
{
	debug_handler_c("module_appl::identitytest", "preparing assertions 1-9", 0,
			std::cout << a << b << c);
	assert((a + b) == (b + a));
	assert((a & b) == (b & a));
	assert((a * b) == (b * a));

	debug_handler_l("module_appl::identitytest", "assertions 1 to 3 succeded", 3);

	assert((a + (b + c)) == ((a + b) + c));
	assert((a + (b + c)) == ((a + c) + b));
	assert((a * (b * c)) == ((a * b) * c));
	assert((a * (b * c)) == ((c * a) * b));
	assert((a & (b & c)) == ((a & b) & c));
	assert((a & (b & c)) == ((c & a) & b));

	debug_handler_l("module_appl::identitytest", "assertions 4 to 9 succeded", 3);

	debug_handler_c("module_appl::identitytest", "preparing assertions 10-13", 0,
			std::cout << a << b << (a+b) << (a&b));
	assert((a + b) >= a);
	assert((a + b) >= b);

	assert((a & b) <= a);
	assert((a & b) <= b);

	debug_handler_l("module_appl::identitytest", "assertions 10 to 13 succeded", 3);
	debug_handler_c("module_appl::identitytest", "preparing assertion 14", 0,
			std::cout << a << b << c << b+c << a*(b+c) << a*b << a*c << (a*b)+(a*c));
	assert((a * (b + c)) == ((a * b) + (a * c)));
	debug_handler_c("module_appl::identitytest", "preparing assertion 15", 3,
			std::cout << a << b << c << (b&c) << (a*(b&c)) << (a*b) << (a*c) << ((a*b)&(a*c)));
	assert((a * (b & c)) <= ((a * b) & (a * c)));
	debug_handler_c("module_appl::identitytest", "preparing assertion 16", 3,
			std::cout << a << b << c << (b+c) << (a&(b+c)) << (a&b) << (a&c) << ((a&b)+(a&c)));
	assert((a & (b + c)) == ((a & b) + (a & c)));

	debug_handler_l("module_appl::identitytest",
			"assertions 14 to 16 succeded", 3);

	debug_handler_c("module_appl::identitytest", "preparing assertion 17", 3,
			std::cout << a << b << (a&b) << ((a&b)+b));
	assert(((a & b) + b) == b);
	debug_handler_c("module_appl::identitytest", "preparing assertion 18", 3,
			std::cout << a << b << (a+b) << ((a+b)&b));
	assert(((a + b) & b) == b);
	debug_handler_c("module_appl::identitytest", "preparing assertion 19", 3,
			std::cout << a << b << (a*b) << ((a*b)&b) << (a&b) << ((a&b)*b));
	//  assert( ((a * b) & b) ==  ((a&b)*b));
	debug_handler_c("module_appl::identitytest", "preparing assertion 20", 3,
			std::cout << a << b << (a*b) << ((a*b)&a) << (a&b) << (a*(a&b)));
	//  assert( ((a * b) & a) ==  (a*(a&b)));
	debug_handler_c("module_appl::identitytest", "preparing assertion 21", 3,
			std::cout << a << b << (a&b) << ((a&b)&a));
	assert(((a & b) & a) == (a & b));

	debug_handler_l("module_appl::identitytest",
			"assertions 17 to 21 succeded", 3);

	if (!exponent(b).is_zero()) {
		debug_handler_c("module_appl::identitytest", "preparing assertion 22", 3,
				std::cout << a << b << (a*b) << (a*b)/b);
		assert(((a * b) / b) >= a);
		debug_handler_c("module_appl::identitytest", "preparing assertion 23", 3,
				std::cout << a << b << (a/b) << std::flush; std::cout << (a/b)*b << std::flush);
		assert(((a / b) * b) <= a);
	}
	if (!exponent(a).is_zero()) {
		debug_handler_c("module_appl::identitytest", "preparing assertion 24", 3,
				std::cout << a << b << (a*b) << std::flush; std::cout << (a*b)/a << std::flush);
		assert(((a * b) / a) >= b);
		debug_handler_c("module_appl::identitytest", "preparing assertion 25", 3,
				std::cout << a << b << (b/a) << (b/a)*a << std::flush);
		assert(((b / a) * a) <= b);
	}

	debug_handler_l("module_appl::identitytest",
			"assertions 22 to 25 succeded", 3);

	module d;

	add(d, b, a); assert((a + b) == d);
	intersect(d, a , b); assert((a & b) == d);
	multiply(d, a, b); assert(d == (b * a));
	d.assign(b); multiply(d, a, d); assert(d == (a * b));
	intersect(d, a, b); assert((a & b) == d);
	if (!exponent(b).is_zero()) {
		divide(d, a, b); assert((a / b) == d);
	}
	add(d, b, c); multiply(d, a, d); assert(d == ((a * b) + (a * c)));
	debug_handler_l("module_appl::identitytest",
			"procedural versions completed successfully", 3);
}



void utiltest(module a) // MM, removed "&" for MS VC++ 5.0
{
	module b, c;

	square(b, a);
	multiply(c, a, a);
	assert(b == c);

	module x = alg_number(1), y = alg_number(1);

	for (int i = 0; i < 10; ++i) {
		power(b, a, i);
		debug_handler_c("module_appl::utiltest", "power", 2,
				std::cout << a << " ^ " << i << " is " << b << std::endl);
		assert(b == x);
		x *= a;
		y = a * y;
		assert(y == x);
	}
	x.assign_one();
	assert(x.is_one());
}



void accumtest(module& a, module& b, module& c)
{
	module x = a;
	x *= b;
	assert(x == (b * a));
	x += c;
	assert(x == (c+(b * a)));
	x &= a;
	assert(x == ((c & a) + ((b * a) & a)));
	if (!exponent(b).is_zero()) {
		x /= b;
		debug_handler_c("module_appl::acummtest", "preparing assertion 4", 3,
				std::cout << x << (((c+b*a)&a)/b) << ((c&a)/b+((b*a)&a)/b) << std::flush);
		assert(x == (((c & a) + ((b * a)& a))/b));
	}

	x.assign(a);
	assert(x == a);

	debug_handler_l("module_appl::accumtest",
			"assertions 1 to 5 succeded", 3);

	multiply(x, x, b);
	debug_handler_c("module_appl::acummtest", "preparing assertion 6", 3,
			std::cout << a << b << x << (b*a) << std::flush);
	assert(x == (b * a));
	add(x, x, c);
	debug_handler_c("module_appl::acummtest", "preparing assertion 7", 3,
			std::cout << c << x << c+(b*a) << std::flush);
	assert(x == (c+(b * a)));
	intersect(x, x, a);
	debug_handler_c("module_appl::acummtest", "preparing assertion 8", 3,
			std::cout << a << x << ((c+(b*a))&a) << std::flush);
	assert(x == ((c + (b * a))&a));
	debug_handler_c("module_appl::acummtest", "preparing assertion 9", 3,
			std::cout << ((c&a)+((b*a)&a)) << std::flush);
	assert(x == ((c & a) + ((b * a) & a)));
	if (!exponent(b).is_zero()) {
		divide(x, x, b);
		debug_handler_c("module_appl::acummtest", "preparing assertion 10", 3,
				std::cout << x << (((c+(b*a))&a)/b) << std::flush);
		assert(x == (((c + (b * a))&a)/b));
	}
}



void anothertest()
{
	module pow64;
	power(pow64, module(bigint(2)), 64);

	bigint s64 = 1;
	s64 <<= 64;

	assert(alg_number(s64) == pow64);
	assert(!(alg_number(s64) != pow64));

	bigint s32 = s64 >> 32;
	assert(!(pow64 == alg_number(s32)));
	assert(pow64 != alg_number(s32));

	debug_handler_l("module_appl::anothertest",
			"simple tests succeded, starting identitytest", 4);
	identitytest(module(s64), module(s32), pow64);
	debug_handler_l("module_appl::anothertest",
			"identitytest succeded, starting accumtest", 4);
	//  accumtest(pow64, alg_number(s32), pow64);
	debug_handler_l("module_appl::anothertest",
			"accumtest succeded, starting utiltest", 4);
	utiltest(alg_number(s32));
}



#if 0
void iotest()
{
	module result;

	std::cout << "\nenter an module: ";
	std::cin >> result;
	std::cout << "module = " << result << "\n";
}
#endif



void all_test(const bigint & ii)
{
	module a, b;

	debug_handler_l("module_appl", "start all_test(const bigint &)", 4);
	a.randomize(ii);
	b.randomize(ii);
//   std::cin >> a >> b;

	std::cout << "\n\n Test with \na = " << a << "\nb = " << b << std::endl;

	debug_handler_l("module_appl", "call utiltest from all_test", 4);
	utiltest(a);
	debug_handler_l("module_appl", "call identitytest from all_test", 4);
	identitytest(a, b, b);
	debug_handler_l("module_appl", "call identitytest again from all_test", 4);
	identitytest(a, a, b);

	debug_handler_l("module_appl", "call accumtest from all_test", 4);
	//  accumtest(a, a, b);
	debug_handler_l("module_appl", "call accumtest again from all_test", 4);
	//  accumtest(b, a, b);

	debug_handler_l("module_appl", "call utiltest again from all_test", 4);
	utiltest(a);
	debug_handler_l("module_appl", "call utiltest (3) from all_test", 4);
	utiltest(b);
}



int main_LiDIA(int argc, char** argv)
{
	order O;
	std::cout << "Order:\n" << std::flush;
	in >> O;

	std::cout << "Computing the maximal order needs:" << std::flush;
	timer T;
	T.start_timer();
	order ext = O.maximize();
	T.stop_timer();
	std::cout << T << std::endl << std::flush;

// The following 3 lines lead to an internal compiler error,
// if gcc-2.6.2 is used:
	module one(alg_number(static_cast< bigint > (1)),
		   alg_number(static_cast< bigint > (0)));
	std::cout << "one = " << one << "\n";
	assert(one == one + alg_number(static_cast< bigint > (1)));
	assert(one.is_one());
	assert(!one.is_zero());

	module n(alg_number(static_cast< bigint > (0)),
		 alg_number(static_cast< bigint > (0)));
	assert(n.is_zero());
	assert(!n.is_one());
	debug_handler_l("module_appl",
			"simple tests succeded, starting anothertest", 5);
	anothertest();
	debug_handler_l("module_appl", "anothertest succeded, starting iotest", 5);
	// iotest();
	debug_handler_l("module_appl", "simple tests succeded, "
			"starting random tests", 5);

	bigint i = 129099121;

	for (int j = 0; j < 10; j++)
		all_test(i);

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
