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


#include	"LiDIA/bigfloat.h"
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



#define LIDIA_DEBUG

#define lidia_assert(ex) if (ex == 0) \
                         { \
		           std::cerr << "\nbigfloat::lidia_assertion failed::"; \
		           std::cerr << "file \"" << __FILE__ << "\"::"; \
		           std::cerr << "line \"" << __LINE__ << "\"\n"; \
			   abort(); \
                         }
inline static void
close_enough(const bigfloat & a, const bigfloat & b)
{
	bigfloat r = a - b;
	bool x = r.is_approx_zero();
	if(!x) {
	  std::clog << "FAILURE: " << a << " - " << b << " = " << r
		    << " is not sufficiently close to zero." << std::endl; 
	}
}



void accu_test(const bigfloat &a,
	       const bigfloat &b,
	       const bigfloat &c)
{
	bigfloat x = a;
	x *= b;
	close_enough(x, (a * b));
	x += c;
	close_enough(x, ((a * b) + c));
	x -= a;
	close_enough(x, (((a * b) + c) - a));
	x /= b;
	close_enough(x, ((((a * b) + c) - a) / b));
}



void identity_test(const bigfloat &a,
		   const bigfloat &b,
		   const bigfloat &c)
{

	// (|R, +, *) is a field
	// (|R, +), (|R, *) are abelian groups and moreover * is
	// distributive in relation to +

	// Test (|R, +)
	// additive inverse
	close_enough(-(-a) , a);

	// additive abelian
	close_enough((a + b) , (b + a));

	// additive associative
	close_enough(((a + b) + c) , (a + (b + c)));

	// additive combinations
	close_enough((a + b) , -((-b) - a));
	close_enough((a - b) , -(b - a));
	close_enough((a + (-b)) , (a - b));
	close_enough(((a - b) + b) , a);
	close_enough(((a + b) - b) , a);

	// Test (|R, *)
	// multiplicative inverse
	bigfloat one = bigfloat(1.0);
	close_enough((one / (one / a)) , a);

	// multiplicative abelian
	close_enough((a * b) , (b * a));

	// multiplicative associative
	close_enough(((a * b) * c) , (a * (b * c)));
	close_enough(((a / b) / c) , (a / (b * c)));
	close_enough((a / (b / c)) , ((a / b) * c));

	// multiplicative combinations
	close_enough(((a * b) / b) , a);
	close_enough((b * (a / b)) , a);
	close_enough((a * (-b)) , -(a * b));
	close_enough((a / (-b)) , -(a / b));

	// multiplicative distributive
	close_enough((a * (b + c)) , ((a * b) + (a * c)));
	close_enough((a * (b - c)) , ((a * b) - (a * c)));
	close_enough(((a + b) / c) , ((a / c) + (b / c)));
}



void util_test(const bigfloat &a,
	       const bigfloat &b,
	       const bigfloat &c)
{
	bigfloat abs_b = abs(b);
	bigfloat abs_c = abs(c);

	// test sqrt(x)
	close_enough(sqrt(abs_c) * sqrt(abs_c) , abs_c);
	close_enough(abs_c / sqrt(abs_c) , sqrt(abs_c));

	// test exp(x), log(x), E()
	bigfloat y, z;

	exp(y, b, 8 * b.get_precision());
	exp(z, c, 8 * c.get_precision());
	z *= y;
	exp(y, b + c, 8 * b.get_precision());
	close_enough(y , z);

	exp(y, -abs_c, 8 * abs_c.get_precision());
	exp(z, abs_c, 8 * abs_c.get_precision());
	z = 1/z;
	close_enough(y , z);

	power(y, abs_c, a);
	log(z, abs_c, 8 * abs_c.get_precision());
	exp(z, a * z, 8 * z.get_precision());
	close_enough(z , y);

	power(y, exp(abs_c), 3);
	exp(z, 3 * abs_c, 8 * abs_c.precision());
	close_enough(z , y);

	power(y, E(), 3);
	exp(z, bigfloat(3.0), 8 * E().get_precision());
	close_enough(z , y);

	log(z, E(), 8 * E().get_precision());
	close_enough(z , 1);

	log(y, abs_b * abs_b, 8 * abs_b.get_precision());
	log(z, abs_b, 8 * abs_b.get_precision());
	close_enough(y , 2 * z);

	log(y, abs_b, 8 * abs_b.get_precision());
	log(z, abs_c, 8 * abs_c.get_precision());
	z += y;
	log(y, abs_b * abs_c, 8 * abs_b.get_precision());
	close_enough(y , z);

	log(y, 1/abs_b, 8 * abs_b.get_precision());
	log(z, abs_b, 8 * abs_b.get_precision());
	z = -z;
	close_enough(y , z);

	log(z, E(), 8 * E().get_precision());
	exp(y, z, 8 * z.get_precision());
	close_enough(y , E());

	// test Pi, sin(x), cos(x)
	close_enough(sin(0) , 0);
	close_enough(sin(Pi()/2) , 1);
	close_enough(sin(Pi()) , 0);
	close_enough(sin(3 * (Pi()/2)) , -1);
	close_enough(sin(2 * Pi()) , 0);

	close_enough(sin(a + 2 * Pi()) , sin(a));
	close_enough(sin(a + Pi()) , -sin(a));
	close_enough(sin(Pi()/2 - a) , cos(a));

	close_enough(cos(0) , 1);
	close_enough(cos(Pi()/2) , 0);
	close_enough(cos(Pi()) , -1);
	close_enough(cos(3 * (Pi()/2)) , 0);
	close_enough(cos(2 * Pi()) , 1);

	close_enough(cos(a + 2 * Pi()) , cos(a));
	close_enough(cos(a + Pi()) , -cos(a));
	close_enough(cos(Pi()/2 - a) , sin(a));

	// test tan(x), cot(x)
	close_enough(tan(a) , sin(a) / cos(a));
	close_enough(tan(a + b) , (tan(a)+tan(b))/(1-tan(a)*tan(b)));
	close_enough(cot(b) , cos(b) / sin(b));
	close_enough(cot(2 * b) , (1-tan(b)*tan(b))/(tan(b)+tan(b)));

	// sinh(x), cosh(x), tan(x), cot(x)
	close_enough(sinh(a + b) , cosh(a)*sinh(b)+sinh(a)*cosh(b));
	close_enough(cosh(a + b) , cosh(a)*cosh(b)+sinh(a)*sinh(b));
	close_enough(cosh(a)*cosh(a)-sinh(a)*sinh(a) , 1);
	close_enough(tanh(a) , sinh(a) / cosh(a));
	close_enough(coth(a) , cosh(a) / sinh(a));

	// test floor(x), ceil(x), truncate(x), round(x)
	close_enough(floor(Pi()) , 3);
	close_enough(floor(Pi()) , truncate(Pi()));
	close_enough(round(Pi()) , truncate(Pi()));
	close_enough(ceil(Pi()) , 4);
	close_enough(ceil(a) , -floor(-a));

	bigfloat x = 1;
	for (int i = 0; i < 10; i++) {
		power(y, b, i);
		close_enough(y , x);
		x *= b;
	}
}



void accu_test_int(const bigfloat &a,
		   const bigfloat &b,
		   const int &c)
{
	bigfloat x = a;
	x *= b;
	close_enough(x, (a * b));
	x += c;
	close_enough(x, ((a * b) + c));
	x -= a;
	close_enough(x, (((a * b) + c) - a));
	x /= b;
	close_enough(x, ((((a * b) + c) - a) / b));
}



void identity_test_int(const bigfloat &a,
		       const bigfloat &b,
		       const int &c)
{
	// (|R, +, *) is a field
	// (|R, +), (|R, *) are abelian groups and moreover * is
	// distributive in relation to +

	// Test (|R, +)
	// additive inverse
	close_enough(-(-a) , a);

  // additive abelian
	close_enough((a + b) , (b + a));

	// additive associative
	close_enough(((a + b) + c) , (a + (b + c)));

	// additive combinations
	close_enough((a + b) , -((-b) - a));
	close_enough((a - b) , -(b - a));
	close_enough((a + (-b)) , (a - b));
	close_enough(((a - b) + b) , a);
	close_enough(((a + b) - b) , a);

	// Test (|R, *)
	// multiplicative inverse
	bigfloat one = bigfloat(1.0);
	close_enough((one / (one / a)) , a);

	// multiplicative abelian
	close_enough((a * b) , (b * a));

	// multiplicative associative
	close_enough(((a * b) * c) , (a * (b * c)));
	close_enough(((a / b) / c) , (a / (b * c)));
	close_enough((a / (b / c)) , ((a / b) * c));

	// multiplicative combinations
	close_enough(((a * b) / b) , a);
	close_enough((b * (a / b)) , a);
	close_enough((a * (-b)) , -(a * b));
	close_enough((a / (-b)) , -(a / b));

  // multiplicative distributive
	close_enough((a * (b + c)) , ((a * b) + (a * c)));
	close_enough((a * (b - c)) , ((a * b) - (a * c)));
	close_enough(((a + b) / c) , ((a / c) + (b / c)));
}



void util_test_int(const bigfloat &a,
		   const bigfloat &b,
		   const int &c)
{
	bigfloat x = 1, y;
	int i;
	for (i = 0; i < 10; i++) {
		power(y, b, i);
		close_enough(y , x);
		x *= b;
	}

	x = 1;
	int abs_c = (c > 0) ? c : -c;
	for (i = 0; i < abs_c; i++) {
		power(y, a, i);
		close_enough(y , x);
		x *= a;
	}
}



void test(const bigfloat &a,
	  const bigfloat &b,
	  const bigfloat &c)
{
#ifdef LIDIA_DEBUG
	std::cout << "\ntest(" << a << ", " << b << ", " << c << ")";
#endif

#ifdef LIDIA_DEBUG
	std::cout << "\n    i) identity test    ";
#endif
	identity_test(a, b, c);
#ifdef LIDIA_DEBUG
	std::cout << "succeeded.";
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
	std::cout << "succeeded.";
#endif
}



void test_int(const bigfloat &a,
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
	std::cout << "succeeded.";
#endif

#ifdef LIDIA_DEBUG
	std::cout << "\n   ii) accumulator test ";
#endif
	accu_test_int(a, b, c);
#ifdef LIDIA_DEBUG
	std::cout << "succeeded.";
#endif

#ifdef LIDIA_DEBUG
	std::cout << "\n  iii) utilities test   ";
#endif
	util_test_int(a, b, c);
#ifdef LIDIA_DEBUG
	std::cout << "succeeded.";
#endif
}



int main_LiDIA(int argc, char** argv)
{
	bigfloat::set_precision(40);

	bigfloat a(1.0);
	bigfloat b(bigfloat(1.0) / 13.0);
	bigfloat c(7.0);

	test(a, a, a);
	test(b, b, b);
	test(c, c, c);
	test(a, b, c);
	test(b, c, a);
	test(c, a, b);

	a = bigfloat(-2.625);
	b = bigfloat(4.0);
	c = bigfloat(-6.5);

	test(a, a, a);
	test(b, b, b);
	test(c, c, c);
	test(a, b, c);
	test(b, c, a);
	test(c, a, b);

	a = bigfloat(-0.5);
	b = bigfloat(-1.0);
	c = bigfloat(0.5);

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
