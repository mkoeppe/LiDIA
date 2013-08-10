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


#include	"LiDIA/bigcomplex.h"
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



#define LIDIA_DEBUG

#define lidia_assert(ex) if (ex == 0) \
                         { \
		           std::cerr << "\nbigcomplex::lidia_assertion failed::"; \
		           std::cerr << "file \"" << __FILE__ << "\"::"; \
		           std::cerr << "line \"" << __LINE__ << "\"\n"; \
			   abort(); \
                         }
static void
close_enough(const bigcomplex & a, const bigcomplex & b)
{
	bigfloat r = real(a) - real(b);
	bigfloat i = imag(a) - imag(b);
	int x = r.is_approx_zero();
	int y = i.is_approx_zero();
	lidia_assert(x && y);
}



void accu_test(const bigcomplex &a,
	       const bigcomplex &b,
	       const bigcomplex &c)
{
	bigcomplex x = a;
	x *= b;
	close_enough(x, (a * b));
	x += c;
	close_enough(x, ((a * b) + c));
	x -= a;
	close_enough(x, (((a * b) + c) - a));
	x /= b;
	close_enough(x, ((((a * b) + c) - a) / b));
}



void identity_test(const bigcomplex &a,
		   const bigcomplex &b,
		   const bigcomplex &c)
{

	// (|C, +, *) is a field
	// (|C, +), (|C, *) are abelian groups and moreover * is
	// distributive in relation to +

	// Test (|C, +)
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

	// Test (|C, *)
	// multiplicative inverse
	bigcomplex one = bigcomplex(1, 0);
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



void util_test(const bigcomplex &a,
	       const bigcomplex &b,
	       const bigcomplex &c)
{
	bigfloat one = 1;
	bigfloat one_half = one / 2;
	bigcomplex one_o_two_i = one / bigcomplex(0, 2);

	close_enough(conj(conj(a)) , a);
	close_enough(conj(a + b) , conj(a) + conj(b));
	close_enough(conj(a * b) , conj(a) * conj(b));

	close_enough(real(a) , one_half * (a + conj(a)));
	close_enough(imag(a) , one_o_two_i * (a - conj(a)));

	close_enough(norm(a) , a * conj(a));
	close_enough(hypot(a) , sqrt(norm(a)));
	close_enough(abs(a) , hypot(a));
	close_enough(abs(a * b) , abs(a) * abs(b));

	close_enough(exp(a * b) , exp(a) * exp(b));
	close_enough(exp(conj(a)) , conj(exp(a)));

	bigcomplex x = one, y;
	for (int i = 0; i < 10; ++i) {
		y = power(c, i);
		close_enough(y , x);
		x *= c;
	}
}



void accu_test_bigfloat(const bigcomplex &a,
	                const bigfloat &b,
	                const bigfloat &c)
{
	bigcomplex x = a;
	x *= b;
	close_enough(x, (a * b));
	x += c;
	close_enough(x, ((a * b) + c));
	x -= a;
	close_enough(x, (((a * b) + c) - a));
	x /= b;
	close_enough(x, ((((a * b) + c) - a) / b));
}



void identity_test_bigfloat(const bigcomplex &a,
		            const bigfloat &b,
		            const bigfloat &c)
{
	// (|C, +, *) is a field
	// (|C, +), (|C, *) are abelian groups and moreover * is
	// distributive in relation to +

	// Test (|C, +)
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

	// Test (|C, *)
	// multiplicative inverse
	bigcomplex one = bigcomplex(1, 0);
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



void util_test_bigfloat(const bigcomplex &a,
	                const bigfloat &b,
	                const bigfloat &c)
{
	bigfloat one = 1;
	bigfloat one_half = one / 2;
	bigcomplex one_o_two_i = one / bigcomplex(0, 2);

	close_enough(conj(conj(a)) , a);
	close_enough(conj(a + b) , conj(a) + conj(b));
	close_enough(conj(a * b) , conj(a) * conj(b));

	close_enough(real(a) , one_half * (a + conj(a)));
	close_enough(imag(a) , one_o_two_i * (a - conj(a)));

	close_enough(norm(a) , a * conj(a));
	close_enough(hypot(a) , sqrt(norm(a)));
	close_enough(abs(a) , hypot(a));
	close_enough(abs(a * b) , abs(a) * abs(b));

	close_enough(exp(a * b) , exp(a) * exp(b));
	close_enough(exp(conj(a)) , conj(exp(a)));

	bigcomplex x = one, y;
	for (int i = 0; i < 10; i++) {
		y = power(c, i);
		close_enough(y , x);
		x *= c;
	}
}



void test(const bigcomplex &a,
	  const bigcomplex &b,
	  const bigcomplex &c)
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
	std::cout << "succeeded.";
#endif

#ifdef LIDIA_DEBUG
	std::cout << "\n  iii) utilities test   ";
#endif
	util_test(a, b, c);
#ifdef LIDIA_DEBUG
	std::cout << "succeeded.";
#endif
}



void test_bigfloat(const bigcomplex &a,
	           const bigfloat &b,
	           const bigfloat &c)
{
#ifdef LIDIA_DEBUG
	std::cout << "\ntest bigfloat(" << a << ", " << b << ", " << c << ")";
#endif

#ifdef LIDIA_DEBUG
	std::cout << "\n    i) identity test    ";
#endif
	identity_test_bigfloat(a, b, c);
#ifdef LIDIA_DEBUG
	std::cout << "succeeded.";
#endif

#ifdef LIDIA_DEBUG
	std::cout << "\n   ii) accumulator test ";
#endif
	accu_test_bigfloat(a, b, c);
#ifdef LIDIA_DEBUG
	std::cout << "succeeded.";
#endif

#ifdef LIDIA_DEBUG
	std::cout << "\n  iii) utilities test   ";
#endif
	util_test_bigfloat(a, b, c);
#ifdef LIDIA_DEBUG
	std::cout << "succeeded.";
#endif
}



void simple_test()
{
	bigcomplex a = bigcomplex(1, 0);
	bigcomplex b = bigcomplex(0, 1);
	bigcomplex c = bigcomplex(0, 0);

	lidia_assert(a.is_real() == true);
	lidia_assert(a.is_imaginary() == false);

	lidia_assert(b.is_real() == false);
	lidia_assert(b.is_imaginary() == true);

	lidia_assert(c.is_zero() == true);
	lidia_assert(a.is_zero() == false);
	lidia_assert(b.is_zero() == false);
}



int main_LiDIA(int argc, char** argv)
{
	bigfloat::set_precision(15);

	bigcomplex a = bigcomplex(1, 0);
	bigcomplex b = bigcomplex(0, 1);
	bigcomplex c = bigcomplex(1, 1);

	simple_test();

	test(a, a, a);
	test(b, b, b);
	test(c, c, c);
	test(a, b, c);
	test(b, c, a);
	test(c, a, b);

	a = bigcomplex(-2, 3);
	b = bigcomplex(4, -5);
	c = bigcomplex(-6, -7);

	test(a, a, a);
	test(b, b, b);
	test(c, c, c);
	test(a, b, c);
	test(b, c, a);
	test(c, a, b);

	bigfloat f = 2;
	bigfloat g = 3;

	test_bigfloat(a, f, g);
	test_bigfloat(b, f, g);
	test_bigfloat(c, f, g);
	test_bigfloat(a, g, f);
	test_bigfloat(b, g, f);
	test_bigfloat(c, g, f);
	test_bigfloat(a, f, f);
	test_bigfloat(b, f, f);
	test_bigfloat(c, f, f);
	test_bigfloat(a, g, g);
	test_bigfloat(b, g, g);
	test_bigfloat(c, g, g);

	f /= 3;
	g /= -17;

	test(a, f, g);
	test(b, f, g);
	test(c, f, g);
	test_bigfloat(a, f, g);
	test_bigfloat(b, f, g);
	test_bigfloat(c, f, g);
	test(a, g, f);
	test(b, g, f);
	test(c, g, f);
	test_bigfloat(a, g, f);
	test_bigfloat(b, g, f);
	test_bigfloat(c, g, f);

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
