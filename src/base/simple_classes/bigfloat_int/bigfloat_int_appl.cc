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
//      Note: Some of these tests may fail due to rounding errors.
//
//	$Id$
//
//	Author	: Thomas Papanikolaou (TP), Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#include	"LiDIA/bigfloat_int.h"
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



#define LIDIA_DEBUG

int return_value = 0;

#define lidia_assert(ex) if (ex == 0) \
                         { \
		           std::cerr << "\nbigfloat_int::lidia_assertion failed::"; \
		           std::cerr << "file \"" << __FILE__ << "\"::"; \
		           std::cerr << "line \"" << __LINE__ << "\"\n"; \
			   abort(); \
                         }

#define CLOSE_ENOUGH(x, y) \
  do { \
    bigfloat_int X(x); \
    bigfloat_int Y(y); \
    if(!close_enough(X, Y)) { \
	std::clog << "FAILURE: close_enough(" # x ", " # y  ") = false " \
		  << "in line " << __LINE__ << " of " << __FILE__ << ".\n"; \
	std::clog << "         " # x " == " << X << "\n"; \
	std::clog << "         " # y " == " << Y << "\n"; \
	std::clog << "         difference == " << (X - Y) << std::endl; \
        return_value = 1; \
    } \
  } while(false)

inline bool
close_enough(const bigfloat_int & a, const bigfloat_int & b)
{
	bigfloat_int r = a - b;
	bigfloat d;
	r.approx(d);
	return d.is_approx_zero();
}



void accu_test(const bigfloat_int &a,
	       const bigfloat_int &b,
	       const bigfloat_int &c)
{
	bigfloat_int x = a;
	x *= b;
	CLOSE_ENOUGH(x, (a * b));
	x += c;
	CLOSE_ENOUGH(x, ((a * b) + c));
	x -= a;
	CLOSE_ENOUGH(x, (((a * b) + c) - a));
	x /= b;
	CLOSE_ENOUGH(x, ((((a * b) + c) - a) / b));
}



void identity_test(const bigfloat_int &a,
		   const bigfloat_int &b,
		   const bigfloat_int &c)
{

	// (|R, +, *) is a field
	// (|R, +), (|R, *) are abelian groups and moreover * is
	// distributive in relation to +

	// Test (|R, +)
	// additive inverse
	CLOSE_ENOUGH(-(-a) , a);

	// additive abelian
	CLOSE_ENOUGH((a + b) , (b + a));

	// additive associative
	CLOSE_ENOUGH(((a + b) + c) , (a + (b + c)));

	// additive combinations
	CLOSE_ENOUGH((a + b) , -((-b) - a));
	CLOSE_ENOUGH((a - b) , -(b - a));
	CLOSE_ENOUGH((a + (-b)) , (a - b));
	CLOSE_ENOUGH(((a - b) + b) , a);
	CLOSE_ENOUGH(((a + b) - b) , a);

	// Test (|R, *)
	// multiplicative inverse
	bigfloat_int one = bigfloat_int(1);
	CLOSE_ENOUGH((one / (one / a)) , a);

	// multiplicative abelian
	CLOSE_ENOUGH((a * b) , (b * a));

	// multiplicative associative
	CLOSE_ENOUGH(((a * b) * c) , (a * (b * c)));
	CLOSE_ENOUGH(((a / b) / c) , (a / (b * c)));
	CLOSE_ENOUGH((a / (b / c)) , ((a / b) * c));

	// multiplicative combinations
	CLOSE_ENOUGH(((a * b) / b) , a);
	CLOSE_ENOUGH((b * (a / b)) , a);
	CLOSE_ENOUGH((a * (-b)) , -(a * b));
	CLOSE_ENOUGH((a / (-b)) , -(a / b));

	// multiplicative distributive
	CLOSE_ENOUGH((a * (b + c)) , ((a * b) + (a * c)));
	CLOSE_ENOUGH((a * (b - c)) , ((a * b) - (a * c)));
	CLOSE_ENOUGH(((a + b) / c) , ((a / c) + (b / c)));
}



void util_test(const bigfloat_int &a,
	       const bigfloat_int &b,
	       const bigfloat_int &c)
{
	bigfloat_int abs_b = abs(b);
	bigfloat_int abs_c = abs(c);

	// test sqrt(x)
	CLOSE_ENOUGH(sqrt(abs_c) * sqrt(abs_c) , abs_c);
	CLOSE_ENOUGH(abs_c / sqrt(abs_c) , sqrt(abs_c));

	// test exp(x), log(x), E()
	bigfloat_int y;

	CLOSE_ENOUGH(exp(b + c) , exp(b) * exp(c));
	CLOSE_ENOUGH(exp(-abs_c) , 1 / exp(abs_c));
	power(y, abs_c, a);
	CLOSE_ENOUGH(exp(a * log(abs_c)) , y);
	power(y, exp(abs_c), 3);
	CLOSE_ENOUGH(exp(3 * abs_c) , y);
	power(y, E(), 3);
	CLOSE_ENOUGH(exp(bigfloat_int(3)) , y);
	CLOSE_ENOUGH(log(E()) , 1);
	CLOSE_ENOUGH(log(abs_b * abs_b) , 2 * log(abs_b));
	CLOSE_ENOUGH(log(abs_b * abs_c) , log(abs_b) + log(abs_c));
	CLOSE_ENOUGH(log(1/abs_b) , -log(abs_b));
	CLOSE_ENOUGH(exp(log(E())) , E());

	// test Pi, sin(x), cos(x)
	CLOSE_ENOUGH(sin(0) , 0);
	CLOSE_ENOUGH(sin(Pi()/2) , 1);
	CLOSE_ENOUGH(sin(Pi()) , 0);
	CLOSE_ENOUGH(sin(3 * (Pi()/2)) , -1);
	CLOSE_ENOUGH(sin(2 * Pi()) , 0);

	CLOSE_ENOUGH(sin(a + 2 * Pi()) , sin(a));
	CLOSE_ENOUGH(sin(a + Pi()) , -sin(a));
	CLOSE_ENOUGH(sin(Pi()/2 - a) , cos(a));

	CLOSE_ENOUGH(cos(0) , 1);
	CLOSE_ENOUGH(cos(Pi()/2) , 0);
	CLOSE_ENOUGH(cos(Pi()) , -1);
	CLOSE_ENOUGH(cos(3 * (Pi()/2)) , 0);
	CLOSE_ENOUGH(cos(2 * Pi()) , 1);

	CLOSE_ENOUGH(cos(a + 2 * Pi()) , cos(a));
	CLOSE_ENOUGH(cos(a + Pi()) , -cos(a));
	CLOSE_ENOUGH(cos(Pi()/2 - a) , sin(a));

	// test tan(x), cot(x)
	CLOSE_ENOUGH(tan(a) , sin(a) / cos(a));
	CLOSE_ENOUGH(tan(a + b) , (tan(a)+tan(b))/(1-tan(a)*tan(b)));
	CLOSE_ENOUGH(cot(b) , cos(b) / sin(b));
	CLOSE_ENOUGH(cot(2 * b) , (1-tan(b)*tan(b))/(tan(b)+tan(b)));

	// sinh(x), cosh(x), tan(x), cot(x)
	CLOSE_ENOUGH(sinh(a + b) , cosh(a)*sinh(b)+sinh(a)*cosh(b));
	CLOSE_ENOUGH(cosh(a + b) , cosh(a)*cosh(b)+sinh(a)*sinh(b));
	CLOSE_ENOUGH(cosh(a)*cosh(a)-sinh(a)*sinh(a) , 1);
	CLOSE_ENOUGH(tanh(a) , sinh(a) / cosh(a));
	CLOSE_ENOUGH(coth(a) , cosh(a) / sinh(a));

	// test floor(x), ceil(x), truncate(x), round(x)
	CLOSE_ENOUGH(floor(Pi()) , 3);
	CLOSE_ENOUGH(floor(Pi()) , truncate(Pi()));
	CLOSE_ENOUGH(round(Pi()) , truncate(Pi()));
	CLOSE_ENOUGH(ceil(Pi()) , 4);
	CLOSE_ENOUGH(ceil(a) , -floor(-a));

	bigfloat_int x = 1;
	for (int i = 0; i < 10; i++) {
		power(y, b, i);
		CLOSE_ENOUGH(y , x);
		x *= b;
	}
}



void accu_test_int(const bigfloat_int &a,
		   const bigfloat_int &b,
		   const int &c)
{
	bigfloat_int x = a;
	x *= b;
	CLOSE_ENOUGH(x, (a * b));
	x += c;
	CLOSE_ENOUGH(x, ((a * b) + c));
	x -= a;
	CLOSE_ENOUGH(x, (((a * b) + c) - a));
	x /= b;
	CLOSE_ENOUGH(x, ((((a * b) + c) - a) / b));
}



void identity_test_int(const bigfloat_int &a,
		       const bigfloat_int &b,
		       const int &c)
{
	// (|R, +, *) is a field
	// (|R, +), (|R, *) are abelian groups and moreover * is
	// distributive in relation to +

	// Test (|R, +)
	// additive inverse
	CLOSE_ENOUGH(-(-a) , a);

  // additive abelian
	CLOSE_ENOUGH((a + b) , (b + a));

	// additive associative
	CLOSE_ENOUGH(((a + b) + c) , (a + (b + c)));

	// additive combinations
	CLOSE_ENOUGH((a + b) , -((-b) - a));
	CLOSE_ENOUGH((a - b) , -(b - a));
	CLOSE_ENOUGH((a + (-b)) , (a - b));
	CLOSE_ENOUGH(((a - b) + b) , a);
	CLOSE_ENOUGH(((a + b) - b) , a);

	// Test (|R, *)
	// multiplicative inverse
	bigfloat_int one = bigfloat_int(1);
	CLOSE_ENOUGH((one / (one / a)) , a);

	// multiplicative abelian
	CLOSE_ENOUGH((a * b) , (b * a));

	// multiplicative associative
	CLOSE_ENOUGH(((a * b) * c) , (a * (b * c)));
	CLOSE_ENOUGH(((a / b) / c) , (a / (b * c)));
	CLOSE_ENOUGH((a / (b / c)) , ((a / b) * c));

	// multiplicative combinations
	CLOSE_ENOUGH(((a * b) / b) , a);
	CLOSE_ENOUGH((b * (a / b)) , a);
	CLOSE_ENOUGH((a * (-b)) , -(a * b));
	CLOSE_ENOUGH((a / (-b)) , -(a / b));

  // multiplicative distributive
	CLOSE_ENOUGH((a * (b + c)) , ((a * b) + (a * c)));
	CLOSE_ENOUGH((a * (b - c)) , ((a * b) - (a * c)));
	CLOSE_ENOUGH(((a + b) / c) , ((a / c) + (b / c)));
}



void util_test_int(const bigfloat_int &a,
		   const bigfloat_int &b,
		   const int &c)
{
	bigfloat_int x = bigfloat(1.0), y;
	int i;
	for (i = 0; i < 10; i++) {
		power(y, b, i);
		CLOSE_ENOUGH(y, x);
		x *= b;
	}

	x = bigfloat(1.0);
	int abs_c = (c > 0) ? c : -c;
	for (i = 0; i < abs_c; i++) {
		power(y, a, i);
		CLOSE_ENOUGH(y , x);
		x *= a;
	}
}



void test(const bigfloat_int &a,
	  const bigfloat_int &b,
	  const bigfloat_int &c)
{
#ifdef LIDIA_DEBUG
	std::cout << "\ntest(" << a << ", " << b << ", " << c << ")";
#endif

#ifdef LIDIA_DEBUG
	std::cout << "\n    i) identity test    ";
#endif
	identity_test(a, b, c);
#ifdef LIDIA_DEBUG
	std::cout << "succeded.";
#endif

#ifdef LIDIA_DEBUG
	std::cout << "\n   ii) accumulator test ";
#endif
	accu_test(a, b, c);
#ifdef LIDIA_DEBUG
	std::cout << "succeded.";
#endif

#ifdef LIDIA_DEBUG
	std::cout << "\n  iii) utilities test   ";
#endif
	util_test(a, b, c);
#ifdef LIDIA_DEBUG
	std::cout << "succeded.";
#endif
}



void test_int(const bigfloat_int &a,
	      const int &b,
	      const int &c)
{
#ifdef LIDIA_DEBUG
	std::cout << "\ntest int(" << a << ", " << b << ", " << c << ")";
#endif

#ifdef LIDIA_DEBUG
	std::cout << "\n    i) identity test    ";
#endif
	identity_test_int(a, b, c);
#ifdef LIDIA_DEBUG
	std::cout << "succeded.";
#endif

#ifdef LIDIA_DEBUG
	std::cout << "\n   ii) accumulator test ";
#endif
	accu_test_int(a, b, c);
#ifdef LIDIA_DEBUG
	std::cout << "succeded.";
#endif

#ifdef LIDIA_DEBUG
	std::cout << "\n  iii) utilities test   ";
#endif
	util_test_int(a, b, c);
#ifdef LIDIA_DEBUG
	std::cout << "succeded.";
#endif
}



int main_LiDIA(int argc, char** argv)
{
	bigfloat::set_precision(40);

	bigfloat_int a = bigfloat_int(1);
	bigfloat_int b = bigfloat_int(1) / 13;
	bigfloat_int c = bigfloat_int(7);

	test(a, a, a);
	test(b, b, b);
	test(c, c, c);
	test(a, b, c);
	test(b, c, a);
	test(c, a, b);

	a = bigfloat_int(-2.625);
	b = bigfloat_int(4.0);
	c = bigfloat_int(-6.5);

	test(a, a, a);
	test(b, b, b);
	test(c, c, c);
	test(a, b, c);
	test(b, c, a);
	test(c, a, b);

	a = bigfloat_int(-0.5);
	b = bigfloat_int(-1.0);
	c = bigfloat_int(0.5);

	test(a, a, a);
	test(b, b, b);
	test(c, c, c);
	test(a, b, c);
	test(b, c, a);
	test(c, a, b);

	int f = -48;
	int g = 6;

	test_int(a, f, g);
	test_int(b, f, g);
	test_int(c, f, g);
	test_int(a, g, f);
	test_int(b, g, f);
	test_int(c, g, f);
	test_int(a, f, f);
	test_int(b, f, f);
	test_int(c, f, f);
	test_int(a, g, g);
	test_int(b, g, g);
	test_int(c, g, g);

	f /= 6;
	g /= 3;

	test(a, f, g);
	test(b, f, g);
	test(c, f, g);
	test_int(a, f, g);
	test_int(b, f, g);
	test_int(c, f, g);
	test(a, g, f);
	test(b, g, f);
	test(c, g, f);
	test_int(a, g, f);
	test_int(b, g, f);
	test_int(c, g, f);

	std::cout << "\nEnd of test\n";

	return return_value;;
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
