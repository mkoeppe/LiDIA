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
#include	"LiDIA/quadratic_number_standard.h"
#include	"LiDIA/quadratic_order.h"
#include	"LiDIA/number_fields/qo_list.h"
#include	"LiDIA/power_functions.h"
#include	"LiDIA/base/b_value.h"
#include	"LiDIA/random_generator.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



char quadratic_number_standard::output_mode = 0;

const char quadratic_number_standard::WITHOUT_DISCRIMINANT = 0;
const char quadratic_number_standard::WITH_DISCRIMINANT = 1;


//
// <implementation>
//

void
quadratic_number_standard::check_and_normalize (char *s)
{
	debug_handler("quadratic_number_standard",
		      "check_and_normalize(char*)");

	if (d.is_zero())
		lidia_error_handler("quadratic_number_standard", s);
	else {
		if (d.is_negative()) {
			a.negate();
			b.negate();
			d.negate();
		}

		bigint g;

		g = gcd (a, b);
		g = gcd (g, d);

		if (g.is_negative())
			g.negate();

		if (!g.is_one()) {
			divide (a, a, g);
			divide (b, b, g);
			divide (d, d, g);
		}
	}
}



//
// Default constructor.
//

quadratic_number_standard::quadratic_number_standard ()
{
	debug_handler("quadratic_number_standard",
		      "quadratic_number_standard ()");

	a.assign_zero();
	b.assign_zero();
	d.assign_one();

	if (quadratic_order::qo_l().last() == &(quadratic_order::zero_QO))
		QO = NULL;
	else {
		QO = quadratic_order::qo_l().add_last_to_list(NULL);
		if (QO == NULL)
			lidia_error_handler("quadratic_number_standard::quadratic_number_standard ()",
					    "Initialization failed.");
	}
}



//
// Copy constructor.
//

quadratic_number_standard::quadratic_number_standard (const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "quadratic_number_standard (const quadratic_number_standard&)");
	a.assign(x.a);
	b.assign(x.b);
	d.assign(x.d);

	if (x.QO == NULL)
		lidia_error_handler("quadratic_number_standard::quadratic_number_standard (qn)",
				    "Quadratic order not initialized !");
	else
		QO = quadratic_order::qo_l().add_to_list (NULL, x.QO);
}



quadratic_number_standard::~quadratic_number_standard ()
{
	debug_handler("quadratic_number_standard",
		      "~quadratic_number_standard ()");

	if (QO != NULL)
		quadratic_order::qo_l().clear(QO);
}



//
//  Input / Output
//

//
// format: (a,b,d)
//

std::istream &
operator >> (std::istream & in, quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "operator >> (std::istream &, quadratic_number_standard &)");

	char c;

	// read white spaces

	in >> c;
	while (c == ' ') in >> c;

	// read '('

	if (c != '(') {
		lidia_error_handler("quadratic_number_standard::operator >>",
				    "'(' expected");
	}
	else {
		// read a
		in >> x.a;

		// read white spaces

		in >> c;
		while (c == ' ') in >> c;

		// read ','

		if (c != ',') {
			lidia_error_handler ("quadratic_number_standard::operator >>",
					     "(a read; ',' expected.");
		}
		else {
			// read b
			in >> x.b;

			// read white spaces
			in >> c;
			while (c == ' ') in >> c;

			// read ','

			if (c != ',') {
				lidia_error_handler ("quadratic_number_standard::operator >>",
						     "(a, b read; ',' expected.");
			}
			else {
				// read d
				in >> x.d;

				// read white spaces
				in >> c;
				while (c == ' ') in >> c;

				// read ')'
				if (c != ')') {
					lidia_error_handler ("quadratic_number_standard::operator >>",
							     "(a, b, d read; ')' expected.");
				}
			}
			// end if ( c != ',' )
		}
		// end if ( c != ',' )
	}
	// end if ( c != '(' )

	x.check_and_normalize("operator >>");
	return in;
}



std::ostream &
operator << (std::ostream & out, const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "operator << (std::ostream &, const quadratic_number_standard &)");

	switch (quadratic_number_standard::output_mode) {
	case quadratic_number_standard::WITHOUT_DISCRIMINANT:
		out << "(" << x.a << ", " << x.b << ", " << x.d << ")";
		break;

	case quadratic_number_standard::WITH_DISCRIMINANT:
		if (x.QO == NULL)
			out << "(" << x.a << ", " << x.b << ", " << x.d << ", NULL)";
		else {
			out << "(" << x.a << ", " << x.b << ", " << x.d << ", ";
			out << x.QO->get_qo()->discriminant() << ")";
		}
		break;

	default:
		lidia_error_handler("quadratic_number_standard::operator << (std::ostream, qn)",
				    "Invalid output mode.");
	}
	return out;
}



//
// Task:
//
//  Reads the quadratic_number_standard x from s and
//  returns the number of read characters.
//

int
string_to_quadratic_number_standard (const char* s, quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "string_to_quadratic_number_standard"
		      "(const char* s, quadratic_number_standard & x)");

	int   i, n;
	const char *t;

	t = s;
	n = 0;

	// read white spaces
	while (*t == ' ') {
		t++;
		n++;
	}

	// read '('
	if (*t != '(') {
		lidia_error_handler("quadratic_number_standard::"
				    "string_to_quadratic_number_standard",
				    "'(' expected");
	}
	else {
		t++;
		n++;

		// read a
		i = string_to_bigint (t, x.a);
		t += i;
		n += i;

		// read white spaces
		while (*t == ' ') {
			t++;
			n++;
		}

		// read ','
		if (*t != ',') {
			lidia_error_handler ("quadratic_number_standard::"
					     "string_to_quadratic_number_standard",
					     "(a read; ',' expected.");
		}
		else {
			t++;
			n++;

			// read b
			i = string_to_bigint (t, x.b);
			t += i;
			n += i;

			// read white spaces
			while (*t == ' ') {
				t++;
				n++;
			}

			// read ','
			if (*t != ',') {
				lidia_error_handler ("quadratic_number_standard::"
						     "string_to_quadratic_number_standard",
						     "(a, b read; ',' expected.");
			}
			else {
				t++;
				n++;

				// read d
				i = string_to_bigint (t, x.d);
				t += i;
				n += i;

				// read white spaces
				while (*t == ' ') {
					t++;
					n++;
				}

				// read ')'
				if (*t != ')') {
					lidia_error_handler("quadratic_number_standard::"
							    "string_to_quadratic_number_standard",
							    "(a, b, d read; ')' expected.");
				}
				else
					n++;
			}
			// end if ( *t != ',' )
		}
		// end if ( *t != ',' )
	}
	// end if ( *t != '(' )

	x.check_and_normalize("string_to_quadratic_number_standard");
	return n;
}



void
swap (quadratic_number_standard & x, quadratic_number_standard & y)
{
	debug_handler("quadratic_number_standard",
		      "swap(quadratic_number_standard & x, quadratic_number_standard & y)");

	swap (x.a, y.a);
	swap (x.b, y.b);
	swap (x.d, y.d);
	swap (x.QO, y.QO);
}



bigint
quadratic_number_standard::get_discriminant () const
{
	debug_handler("quadratic_number_standard",
		      "get_discriminant() const");

	if (QO == NULL) {
		lidia_error_handler("quadratic_number_standard::get_discriminant()",
				    "QO == NULL");
		return bigint(0);
	}

	if (!QO->get_qo()) {
		lidia_error_handler("quadratic_number_standard::get_discriminant()",
				    "Quadratic order has been deleted.");
		return bigint(0);
	}
	return QO->get_qo()->discriminant();
}



void
quadratic_number_standard::get (bigint & pa, bigint & pb, bigint & pd)
{
	debug_handler("quadratic_number_standard",
		      "get (bigint&, bigint&, bigint&)");
	pa.assign(a);
	pb.assign(b);
	pd.assign(d);
}



const bigint &
quadratic_number_standard::get_a () const
{
	debug_handler("quadratic_number_standard",
		      "get_a() const");
	return a;
}



const bigint &
quadratic_number_standard::get_b () const
{
	debug_handler("quadratic_number_standard",
		      "get_b() const");
	return b;
}



const bigint &
quadratic_number_standard::get_d () const
{
	debug_handler("quadratic_number_standard",
		      "get_d() const");
	return d;
}



//
// quadratic_number_standard::which_order()
//
// Task:
//   returns a reference to the quadratic_order to which the ideal belongs.
//

const quadratic_order &
quadratic_number_standard::which_order () const
{
	debug_handler("quadratic_number_standard", "which_order() const");

	if (QO == NULL) {
		lidia_error_handler("quadratic_number_standard::which_order()",
				    "Quadratic order not initialized.!");
		return quadratic_order::zero_QO;
	}

	if (!QO->get_qo()) {
		lidia_error_handler("quadratic_number_standard", "which_order() - the quadratic "
				    "order of this number has been deleted");
		return quadratic_order::zero_QO;
	}

	return *(QO->get_qo());
}



const quadratic_order &
quadratic_number_standard::get_order () const
{
	debug_handler("quadratic_number_standard", "get_order() const");

	if (QO == NULL) {
		lidia_error_handler("quadratic_number_standard::get_order()",
				    "Quadratic order not initialized.!");
		return quadratic_order::zero_QO;
	}

	if (!QO->get_qo()) {
		lidia_error_handler("quadratic_number_standard", "get_order() - the quadratic "
				    "order of this number has been deleted");
		return quadratic_order::zero_QO;
	}

	return *QO->get_qo();
}



//
// Task:
//
//  change the quadratic order to which the number belongs.
//

void
quadratic_number_standard::assign_order (const quadratic_order & QO2)
{
	debug_handler("quadratic_number_standard",
		      "assign_order(const quadratic_number&)");

	QO = quadratic_order::qo_l().add_to_list(QO, *((quadratic_order*)&QO2));
}



void
quadratic_number_standard::set_order (const quadratic_order & QO2)
{
	debug_handler("quadratic_number_standard",
		      "set_order(const quadratic_order&)");

	QO = quadratic_order::qo_l().add_to_list(QO, *((quadratic_order*)&QO2));
}



//
//  Assignment
//

quadratic_number_standard &
quadratic_number_standard::operator = (const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "operator = (const quadratic_number_standard &)");
	if (&x != this) {
		this->assign(x);
	}
	return *this;
}



void
quadratic_number_standard::assign_one ()
{
	debug_handler("quadratic_number_standard",
		      "assign_one()");
	a.assign_one();
	b.assign_zero();
	d.assign_one();
}



void
quadratic_number_standard::multiply_by_denominator ()
{
	debug_handler("quadratic_number_standard", "multiply_by_denominator()");
	d.assign_one();
}



void
quadratic_number_standard::assign_one (const quadratic_order & QO2)
{
	debug_handler("quadratic_number_standard", "assign_one(quadratic_order)");

	QO = quadratic_order::qo_l().add_to_list(QO, *((quadratic_order*)&QO2));
	a.assign_one();
	b.assign_zero();
	d.assign_one();
}



void
quadratic_number_standard::assign_zero ()
{
	debug_handler("quadratic_number_standard",
		      "assign_zero()");

	a.assign_zero();
	b.assign_zero();
	d.assign_one();
}



void
quadratic_number_standard::assign (const bigint & pa, const bigint & pb, const bigint & pd)
{
	debug_handler("quadratic_number_standard",
		      "assign (const bigint&, const bigint&, const bigint&)");

	a.assign(pa);
	b.assign(pb);
	d.assign(pd);

	this->check_and_normalize("assign(const bigint&, const bigint&, const bigint&)");
}



void
quadratic_number_standard::assign (const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "assign (const quadratic_number_standard & x)");

	a.assign(x.a);
	b.assign(x.b);
	d.assign(x.d);

	if (x.QO == NULL)
		lidia_error_handler("quadratic_number_standard::assign(qn)",
				    "Quadratic order not initialized!");
	else
		QO = quadratic_order::qo_l().add_to_list(QO, x.QO);
}



void
quadratic_number_standard::randomize ()
{
	debug_handler("quadratic_number_standard", "randomize()");

	if (QO == NULL)
		lidia_error_handler("quadratic_number_standard::randomize()",
				    "Quadratic order not initialized!");
	else {
		bigint Delta = QO->get_qo()->discriminant();
		a.randomize(Delta);
		b.randomize(Delta);
		d.randomize(Delta);
		if (d.is_zero())
			d.assign_one();

		this->check_and_normalize("randomize()");
	}
}



void
quadratic_number_standard::absolute_value ()
{
	debug_handler("quadratic_number_standard", "absolute_value()");

	if (QO == NULL)
		lidia_error_handler("quadratic_number_standard::absolute_value()",
				    "Quadratic order not initialized!");
	else {
		if (this->is_negative())
			this->negate();
	}
}



void
quadratic_number_standard::set (const bigint & pa, const bigint & pb, const bigint & pd)
{
	debug_handler("quadratic_number_standard",
		      "set(const bigint&, const bigint&, const bigint&");
	this->assign(pa, pb, pd);
}



//
// Comparison
//

bool
operator == (const quadratic_number_standard & x,
	     const quadratic_number_standard & y)
{
	debug_handler("quadratic_number_standard",
		      "operator == (const qns &, const qns &)");

	if (x.QO == NULL) {
		lidia_error_handler("quadratic_number_standard::which_order()",
				    "Quadratic order of x not initialized!");
		return false;
	}

	if (y.QO == NULL) {
		lidia_error_handler("quadratic_number_standard::which_order()",
				    "Quadratic order of y not initialized!");
		return false;
	}

	if ((x.a == y.a) &&
	    (x.b == y.b) &&
	    (x.d == y.d) &&
	    (x.QO == y.QO))
		return true;
	else
		return false;
}



bool
operator == (const quadratic_number_standard & x,
	     const bigint & n)
{
	debug_handler("quadratic_number_standard",
		      "operator == (qns, const bigint &)");

	if (x.QO == NULL) {
		lidia_error_handler("quadratic_number_standard::operator == (qn, bigint)",
				    "Quadratic order not initialized!");
		return false;
	}

	if (!x.d.is_one())
		return false;
	else if (!x.b.is_zero())
		return false;
	else
		return (x.a == n);
}



bool
operator == (const bigint & n,
	     const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "operator == (const bigint&, qns)");

	return (x == n);
}



bool
operator == (const quadratic_number_standard & x,
	     const bigrational & n)
{
	debug_handler("quadratic_number_standard",
		      "operator == (qns, const bigrational &)");

	if (x.QO == NULL) {
		lidia_error_handler("quadratic_number_standard::operator == (qn, bigrational)",
				    "Quadratic order not initialized!");
		return false;
	}

	if (!x.b.is_zero())
		return false;
	else
		return (n.numerator() == x.a && n.denominator() == x.d);
}



bool
operator == (const bigrational & n, const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "operator == (const bigrational &, qns)");
	return (x == n);
}



bool
operator != (const quadratic_number_standard & x,
   	     const quadratic_number_standard & y)
{
	debug_handler("quadratic_number_standard",
		      "operator != (const qns&, const qns&)");
	return !(x == y);
}



bool
operator != (const bigint & n,
	     const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "operator != (const bigint&, const qns&)");

	return !(x == n);
}



bool
operator != (const quadratic_number_standard & x,
	     const bigint & n)
{
	debug_handler("quadratic_number_standard",
		      "operator != (const qns&, const bigint&)");

	return !(x == n);
}



bool
operator != (const bigrational & n,
	     const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "operator != (const bigrational&, const qns&)");

	return !(x == n);
}



bool
operator != (const quadratic_number_standard & x,
	     const bigrational & n)
{
	debug_handler("quadratic_number_standard",
		      "operator != (const qns&, const bigrational&)");

	return !(x == n);
}



bool
quadratic_number_standard::is_one () const
{
	debug_handler("quadratic_number_standard",
		      "is_one()");

	if (QO == NULL) {
		lidia_error_handler("quadratic_number_standard::is_one()",
				    "Quadratic order not initialized!");
		return false;
	}
	else if (!d.is_one())
		return false;
	else if (!b.is_zero())
		return false;
	else
		return a.is_one();
}



bool
quadratic_number_standard::is_zero () const
{
	debug_handler("quadratic_number_standard", "is_zero()");

	if (QO == NULL) {
		lidia_error_handler("quadratic_number_standard::is_zero()",
				    "Quadratic order not initialized!");
		return false;
	}
	else
		return (a.is_zero() && b.is_zero());
}



//
// is_negative()
//
// Returns true, if the number is less than zero.
//

bool
quadratic_number_standard::is_negative () const
{
	debug_handler("quadratic_number_standard", "is_negative()");

	if (QO == NULL) {
		lidia_error_handler("quadratic_number_standard::is_negative()",
				    "Quadratic order not initialized!");
		return false;
	}
	else {
		// d always is positive.
		//

		// negative discriminants
		//
		if ((this->get_discriminant()).is_negative())
			return a.is_negative();

		// positive discriminants
		//

		// a+b*sqrt(D) = sign(a)|a|+sign(b)|b|*sqrt(D)
		//
		else if (a.is_positive() && b.is_positive())
			return false;
		else if (a.is_negative() && b.is_negative())
			return true;
		else if (b.is_zero())
			return a.is_negative();
		else {
			// a+b*sqrt(D) = sign(b) * (-|a|+|b|*sqrt(D))
			//
			bigint a_squared, b_squared_D;
			square (a_squared, a);
			square (b_squared_D, b);
			multiply (b_squared_D, b_squared_D, this->get_discriminant());

			if (b_squared_D > a_squared)
				return b.is_negative();
			else
				return b.is_positive();
		}
	}
}



//
// is_positive()
//
// Returns true, if the number is greater than or equal to zero.
//

bool
quadratic_number_standard::is_positive () const
{
	debug_handler("quadratic_number_standard", "is_positive()");
	if (is_negative())
		return false;
	else
		return true;
}



//
// arithmetic
//

int
quadratic_number_standard::get_sign () const
{
	debug_handler("quadratic_number_standard", "get_sign() const");

	if (this->is_negative())
		return -1;
	else
		return 1;
}



void
quadratic_number_standard::negate ()
{
	debug_handler("quadratic_number_standard", "negate()");

	if (QO == NULL)
		lidia_error_handler("quadratic_number_standard::negate()",
				    "Quadratic order not initialized!");
	else {
		a.negate();
		b.negate();
	}
}



void
negate (quadratic_number_standard & z,
	const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "negate(qn, qn)");

	if (x.QO == NULL)
		lidia_error_handler("negate(qn, qn)",
				    "Quadratic order not initialized!");
	negate (z.a, x.a);
	negate (z.b, x.b);
	z.d.assign(x.d);
	z.QO = quadratic_order::qo_l().add_to_list(z.QO, x.QO);
}



void
add (quadratic_number_standard & z,
     const quadratic_number_standard & x,
     const quadratic_number_standard & y)
{
	debug_handler("quadratic_number_standard",
		      "add(qn, qn, qn)");

	if (x.QO != y.QO) {
		lidia_error_handler ("quadratic_number_standard::add(qn, qn, qn)",
				     "quadratic orders are different.");
		return;
	}
	else if (x.QO == NULL) {
		lidia_error_handler("add(qn, qn, qn)",
				    "Quadratic order not initialized!");
		return;
	}

	if (x.is_zero()) {
		z.assign(y);
	}
	else if (y.is_zero()) {
		z.assign(x);
	}
	else {
		bigint g = gcd(x.d, y.d);
		bigint h1;
		if (g.is_one()) {
			//     y.d * x.a + x.d * y.a + sqrt(D) * (y.d * x.b + x.d * y.b)
			// z = ---------------------------------------------------------
			//           		 x.d * y.d

			multiply (h1, y.d, x.a);
			multiply (z.a, x.d, y.a);
			add      (z.a, z.a, h1);

			multiply (h1, y.d, x.b);
			multiply (z.b, x.d, y.b);
			add      (z.b, z.b, h1);

			multiply(z.d, x.d, y.d);
		}
		else {
			// h1 = (y.d / g) * x.a + (x.d / g) * y.a
			// h2 = (y.d / g) * x.b + (x.d / g) * y.b
			// h = gcd(h1, h2, g);
			//
			//       (h1/h) + sqrt(D) * (h2/h)
			// z =  ---------------------------
			//        (x.d / g) * (y.d / h)

			bigint h2, h, s;

			divide   (s, y.d, g);

			multiply (h1, s, x.a);
			multiply (h2, s, x.b);

			divide   (s, x.d, g);

			multiply (h, s, y.a);
			add      (h1, h1, h);

			multiply (h, s, y.b);
			add      (h2, h2, h);

			h = gcd(h1, g);
			h = gcd(h2, h);

			divide   (z.d, y.d, h);
			multiply (z.d, z.d, s);

			divide (z.a, h1, h);
			divide (z.b, h2, h);
		}
	}

	z.QO = quadratic_order::qo_l().add_to_list(z.QO, x.QO);
}



void
subtract (quadratic_number_standard & z,
	  const quadratic_number_standard & x,
	  const quadratic_number_standard & y)
{
	debug_handler("quadratic_number_standard",
		      "subtract(qn, qn, qn)");

	quadratic_number_standard my;
	my = y;
	my.negate ();
	add (z, x, my);
}



//
// Multiplication
//

void
multiply (quadratic_number_standard & z,
	  const quadratic_number_standard & x,
	  const quadratic_number_standard & y)
{
	debug_handler("quadratic_number_standard", "multiply(qn, qn, qn)");

	if (x.QO != y.QO) {
		lidia_error_handler ("quadratic_number_standard::multiply(qn, qn, qn)",
				     "Quadratic orders are different.");
		return;
	}
	else if (x.QO == NULL) {
		lidia_error_handler("multiply(qn, qn, qn)",
				    "Quadratic order not initialized!");
		return;
	}

	bigint h1, h2, h3;

	// h2 = x.b y.b Delta
	multiply (h2, x.b, y.b);
	h2 *= y.get_discriminant();

	// h1 = x.a y.a + x.b y.b Delta
	multiply (h1, x.a, y.a);
	h1 += h2;

	// h2 = x.a y.b
	multiply (h2, x.a, y.b);

	// h3 = x.b y.a
	multiply (h3, x.b, y.a);

	// h2 = (x.a y.b + x.b y.a)
	add (z.b, h2, h3);

	z.a.assign(h1);
	multiply(z.d, x.d, y.d);

	z.QO = quadratic_order::qo_l().add_to_list(z.QO, x.QO);

	z.check_and_normalize("multiply(quadratic_number_standard&, const qn&, const qn&)");
}



void
multiply (quadratic_number_standard & z,
	  const bigint & n,
	  const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard", "multiply(qn, bigint, qn)");

	if (x.QO == NULL) {
		lidia_error_handler("multiply(qn, qn, bigint)",
				    "Quadratic order not initialized!");
		return;
	}
	z.QO = quadratic_order::qo_l().add_to_list(z.QO, x.QO);

	multiply (z.a, n, x.a);
	multiply (z.b, n, x.b);
	z.d.assign (x.d);

	z.check_and_normalize("multiply(quadratic_number_standard&, const bigint&, const qn&)");
}



void
multiply (quadratic_number_standard & z,
	  const quadratic_number_standard & x,
	  const bigint & n)
{
	debug_handler("quadratic_number_standard",
		      "multiply(qn, qn, bigint)");
	if (x.QO == NULL) {
		lidia_error_handler("multiply(qn, qn, bigint)",
				    "Quadratic order not initialized!");
		return;
	}
	z.QO = quadratic_order::qo_l().add_to_list(z.QO, x.QO);

	multiply (z.a, n, x.a);
	multiply (z.b, n, x.b);
	z.d.assign (x.d);

	z.check_and_normalize("multiply(quadratic_number_standard&, const qn&, const bigint&)");
}



void
multiply (quadratic_number_standard & z,
	  const bigrational & n,
	  const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard", "multiply(qn, bigrational, qn)");

	if (x.QO == NULL) {
		lidia_error_handler("multiply(qn, bigrational, qn)",
				    "Quadratic order not initialized!");
		return;
	}
	z.QO = quadratic_order::qo_l().add_to_list(z.QO, x.QO);

	multiply (z.a, n.numerator(), x.a);
	multiply (z.b, n.numerator(), x.b);
	multiply (z.d, n.denominator(), x.d);

	z.check_and_normalize("multiply(qn, bigrational, qn)");
}



void
multiply (quadratic_number_standard & z,
	  const quadratic_number_standard & x,
	  const bigrational & n)
{
	debug_handler("quadratic_number_standard", "multiply(qn, qn, bigrational)");

	if (x.QO == NULL) {
		lidia_error_handler("multiply(qn, qn, bigrational)",
				    "Quadratic order not initialized!");
		return;
	}
	z.QO = quadratic_order::qo_l().add_to_list(z.QO, x.QO);

	multiply (z.a, n.numerator(), x.a);
	multiply (z.b, n.numerator(), x.b);
	multiply (z.d, n.denominator(), x.d);

	z.check_and_normalize("multiply(qn, bigrational, qn)");
}



//
// Division
//

void
divide (quadratic_number_standard & c,
	const quadratic_number_standard & a,
	const quadratic_number_standard & b)
{
	debug_handler("quadratic_number_standard",
		      "divide(qn, qn, qn)");

	quadratic_number_standard inv_b;

	inverse  (inv_b, b);
	multiply (c, a, inv_b);
}



void
divide (quadratic_number_standard       & x,
	const quadratic_number_standard & y,
	const bigrational      & n)
{
	debug_handler("quadratic_number_standard",
		      "divide(qn, qn, bigrational)");

	if (y.QO == NULL) {
		lidia_error_handler("divide(qn, qn, bigrational)",
				    "Quadratic order not initialized!");
		return;
	}

	multiply (x.a, y.a, n.numerator());
	multiply (x.b, y.b, n.numerator());
	multiply (x.d, y.d, n.denominator());
	x.QO = quadratic_order::qo_l().add_to_list(x.QO, y.QO);
	x.check_and_normalize("divide(qn, qn, bigrational)");
}



void
divide (quadratic_number_standard       & x,
	const bigrational      & n,
	const quadratic_number_standard & y)
{
	debug_handler("quadratic_number_standard",
		      "divide(qn, qn, bigrational)");
	divide (x, y, n);
}



void
square (quadratic_number_standard & z,
	const quadratic_number_standard & x)
{
	debug_handler ("quadratic_number_standard", "square(qn, qn)");

	if (x.QO == NULL) {
		lidia_error_handler("square(qn, qn)",
				    "Quadratic order not initialized!");
		return;
	}

	bigint h1, h2;

	// h1 = b^2
	square (h1, x.b);

	// b = 2ab
	multiply (z.b, x.a, x.b);
	shift_left (z.b, z.b, 1);

	// h1 = b^2 Delta
	multiply (h1, h1, x.get_discriminant());

	// h2 = a^2
	square (h2, x.a);

	// a = (h1 + h2)
	add(z.a, h1, h2);

	square (z.d, x.d);
	z.QO = quadratic_order::qo_l().add_to_list(z.QO, x.QO);
	z.check_and_normalize("square(quadratic_number_standard&, qn&, qn&)");
}



//
// ::power
//

void
power (quadratic_number_standard & z,
       const quadratic_number_standard & x,
       unsigned long e)
{
	debug_handler("quadratic_number_standard", "power(qn, qn, unsigned long)");

	if (x.QO == NULL) {
		lidia_error_handler("power(qn, qn, unsigned long)",
				    "Quadratic order not initialized!");
		return;
	}

	if (&x != &z)
		z.QO = quadratic_order::qo_l().add_to_list(z.QO, x.QO);

	quadratic_order* QO2 = quadratic_order::qo_l().last();
	quadratic_order::qo_l().set_last((x.QO)->get_qo());
	lidia_power (z, x, e);
	quadratic_order::qo_l().set_last(QO2);
}



void
power (quadratic_number_standard & z,
       const quadratic_number_standard & x,
       const bigint & e)
{
	debug_handler("quadratic_number_standard", "power(qn, qn, bigint)");

	if (x.QO == NULL) {
		lidia_error_handler("power(qn, qn, bigint)",
				    "Quadratic order not initialized!");
		return;
	}

	if (&x != &z)
		z.QO = quadratic_order::qo_l().add_to_list(z.QO, x.QO);

	quadratic_order* QO2 = quadratic_order::qo_l().last();
	quadratic_order::qo_l().set_last((x.QO)->get_qo());

	if (x.is_one())
		z.assign_one();

	else if (e.is_negative()) {
		lidia_power (z, x, -e);
		z.invert();
	}
	else {
		lidia_power (z, x, e);
	}
	quadratic_order::qo_l().set_last(QO2);
}



void
inverse (quadratic_number_standard & x,
	 const quadratic_number_standard & y)
{
	debug_handler("quadratic_number_standard",
		      "inverse (qn, qn)");
	bigint h1, h2;

	// h2 = (y.b)^2 Delta
	square (h2, y.b);
	h2 *= y.get_discriminant();

	// h1 = (y.a)^2
	square (h1, y.a);

	// x.b = - y.d * y.b
	multiply (x.b, y.d, y.b);
	x.b.negate ();

	// x.a = y.d * y.a
	multiply (x.a, y.d, y.a);

	// x.d = h1 - h2;
	subtract (x.d, h1, h2);

	x.QO = quadratic_order::qo_l().add_to_list(x.QO, y.QO);
	x.check_and_normalize("quadratic_number_standard::inverse(qn, qn)");
}



void
quadratic_number_standard::invert ()
{
	debug_handler("quadratic_number_standard", "invert()");

	bigint h1, h2;

	// h2 = b^2 Delta
	square (h2, b);
	h2 *= get_discriminant();

	// h1 = a^2
	square (h1, a);

	// b = - d * b
	multiply (b, d, b);
	b.negate ();

	// a = d * a
	multiply (a, d, a);

	// d = h1 - h2;
	subtract (d, h1, h2);

	check_and_normalize("quadratic_number_standard::invert()");
}



quadratic_number_standard
operator - (const quadratic_number_standard & a)
{
	debug_handler("quadratic_number_standard",
		      "operator -(const qns&)");

	quadratic_number_standard c(a);
	c.negate();
	return c;
}



quadratic_number_standard
operator + (const quadratic_number_standard & a,
	    const quadratic_number_standard & b)
{
	debug_handler("quadratic_number_standard",
		      "operator + (const qns&, const qns&)");

	quadratic_number_standard c;
	add(c, a, b);
	return c;
}



quadratic_number_standard
operator - (const quadratic_number_standard & a,
	    const quadratic_number_standard & b)
{
	debug_handler("quadratic_number_standard",
		      "operator - (const qns&, const qns&)");

	quadratic_number_standard c;
	subtract(c, a, b);
	return c;
}



quadratic_number_standard
operator * (const quadratic_number_standard & a,
	    const quadratic_number_standard & b)
{
	debug_handler("quadratic_number_standard",
		      "operator * (const qns&, const qns&)");

	quadratic_number_standard c;
	multiply(c, a, b);
	return c;
}



quadratic_number_standard
operator * (const bigint & a,
	    const quadratic_number_standard & b)
{
	debug_handler("quadratic_number_standard",
		      "operator * (const bigint&, const qns&)");

	quadratic_number_standard c;
	multiply(c, a, b);
	return c;
}



quadratic_number_standard
operator * (const quadratic_number_standard & a,
	    const bigint & b)
{
	debug_handler("quadratic_number_standard",
		      "operator - (const qns&, const bigint&)");

	quadratic_number_standard c;
	multiply(c, b, a);
	return c;
}



quadratic_number_standard
operator / (const quadratic_number_standard & a,
	    const quadratic_number_standard & b)
{
	debug_handler("quadratic_number_standard",
		      "operator / (const qns&, const qns&)");

	quadratic_number_standard c;
	divide(c, a, b);
	return c;
}



quadratic_number_standard
operator / (const bigrational & a,
	    const quadratic_number_standard & b)
{
	debug_handler("quadratic_number_standard",
		      "operator / (const bigrational&, const qns&)");

	quadratic_number_standard c;
	divide(c, a, b);
	return c;
}



quadratic_number_standard
operator / (const quadratic_number_standard & a,
	    const bigrational & b)
{
	debug_handler("quadratic_number_standard",
		      "operator / (const qns&, const bigrational&)");

	quadratic_number_standard c;
	divide(c, b, a);
	return c;
}



quadratic_number_standard &
quadratic_number_standard::operator += (const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "operator += (const qns&)");

	add(*this, *this, x);
	return *this;
}



quadratic_number_standard &
quadratic_number_standard::operator -= (const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "operator -= (const qns&)");

	subtract(*this, *this, x);
	return *this;
}



quadratic_number_standard &
quadratic_number_standard::operator *= (const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "operator *= (const qns&)");

	multiply(*this, *this, x);
	return *this;
}



quadratic_number_standard &
quadratic_number_standard::operator /= (const quadratic_number_standard & x)
{
	debug_handler("quadratic_number_standard",
		      "operator /= (const qns&)");

	divide(*this, *this, x);
	return *this;
}



//
// norm
//

void
quadratic_number_standard::norm (bigrational & n) const
{
	debug_handler("quadratic_number_standard",
		      "norm(bigrational &) const");

	bigint h1, h2;

	// h2 = b^2 * Delta
	square   (h2, b);
	h2 *= this->get_discriminant();

	// h1 = a^2 - b^2 * Delta
	square (h1, a);
	h1 -= h2;

	square (h2, d);

	n.assign(h1, h2);
}



bigrational
quadratic_number_standard::norm () const
{
	debug_handler("quadratic_number_standard",
		      "norm() const");

	bigrational n;
	this->norm(n);
	return n;
}



bigrational
norm (const quadratic_number_standard & y)
{
	debug_handler("quadratic_number_standard",
		      "norm(const qns&)");

	bigrational n;
	y.norm(n);
	return n;
}



//
// trace
//

void
quadratic_number_standard::trace (bigrational & n) const
{
	debug_handler("quadratic_number_standard",
		      "trace(bigint &) const");

	if (QO == NULL) {
		lidia_error_handler("quadratic_number_standard::trace(bigrational)",
				    "Quadratic order not initialized!");
		return;
	}

	bigint h;

	shift_left(h, a, 1);
	n.assign(h, d);
}



bigrational
quadratic_number_standard::trace () const
{
	debug_handler("quadratic_number_standard",
		      "trace () const");

	bigrational t;
	this->trace(t);
	return t;
}



bigrational
trace (const quadratic_number_standard & y)
{
	debug_handler("quadratic_number_standard",
		      "trace(const qns&)");

	bigrational t;
	y.trace(t);
	return t;
}



//
// conjugate
//

void
quadratic_number_standard::conjugate ()
{
	debug_handler("quadratic_number_standard",
		      "conjugate()");

	if (QO == NULL) {
		lidia_error_handler("conjugate()",
				    "Quadratic order not initialized!");
		return;
	}

	b.negate();
}



void
conjugate (quadratic_number_standard & x,
	   const quadratic_number_standard & y)
{
	debug_handler("quadratic_number_standard",
		      "conjugate (qns&, const qns&)");

	if (y.QO == NULL) {
		lidia_error_handler("conjugate(qn, qn)",
				    "Quadratic order not initialized!");
		return;
	}

	x.a = y.a;
	x.b = y.b;
	x.b.negate();
	x.d = y.d;

	x.QO = quadratic_order::qo_l().add_to_list(x.QO, y.QO);
}



//
// Test for rational number, integer, unit
//

bool
quadratic_number_standard::is_rational_number () const
{
	debug_handler("quadratic_number_standard",
		      "is_rational_number() const");

	return (b.is_zero());
}



bool
quadratic_number_standard::is_integer () const
{
	debug_handler("quadratic_number_standard",
		      "is_integer() const");

	bigrational n;

	this->norm(n);
	n.absolute_value();

	return (n.denominator().is_one());
}



bool
quadratic_number_standard::is_unit () const
{
	debug_handler("quadratic_number_standard",
		      "is_unit() const");

	bigrational n;

	this->norm(n);
	n.absolute_value();

	return (n.is_one());
}



//
// special purpose stuff
//

long
quadratic_number_standard::b_value () const
{
	std::cout << "quadratic_number_standard::b_value not yet implemented";
	return 0;
}



long
quadratic_number_standard::b_value (const bigint & sq_Delta) const
{
	(void)sq_Delta;
	std::cout << "quadratic_number_standard::b_value not yet implemented";
	return 0;
}



//
//
//  Approximations
//
//

void
quadratic_number_standard::get_relative_approximation (xbigfloat & r, long k) const
{
	debug_handler("quadratic_number_standard",
		      "get_relative_approximation(xbigfloat&, long)");

	xbigfloat D, y;

	if (QO == NULL) {
		lidia_error_handler("quadratic_number_standard::get_relative_approximation"
				    "(xbigfloat, long)",
				    "Quadratic order not initialized!");
		return;
	}

	if ((this->a).is_zero() &&
	    (this->b).is_zero()) {
		r.assign_zero();
		return;
	}

	if ((this->a).is_negative() ||
	    (this->b).is_negative()) {
		lidia_error_handler ("quadratic_number_standard::"
				     "get_relative_approximation",
				     "Negative coefficients.");
		return;
	}

	if (k < 1) {
		lidia_warning_handler ("quadratic_number_standard::"
				       "get_relative_approximation",
				       "Increasing precision to 1.");
		k = 1;
	}

	// truncating to k+6 yields relative k+5 approx.
	sqrt (D, this->get_discriminant(), k+6);
	truncate (y, this->b, k+6);

	// r is relative k+4 approx.
	multiply (r, y, D);

	// truncating to k+6 yields relative k+3 approx.
	r.truncate (k+6);

	// division with k+3 yields relative k+1 approx.
	add (r, r, this->a);
	divide (r, r, this->d, k+3);

	// truncating to k+3 yields relative k approx.
	r.truncate (k+3);
}



//
// obsolete. Does two log computations. Replaced by new
// get_absolute_Ln_approximation with only one log computation.
//

void
quadratic_number_standard::get_absolute_Ln_approximation_old (xbigfloat & l, long k) const
{
	int s;
	bigint num, h;
	quadratic_number_standard p;
	xbigfloat r, den;

	if (QO == NULL) {
		lidia_error_handler("quadratic_number_standard::absolute_Ln_approximation"
				    "(xbigfloat, long) const",
				    "Quadratic order not initialized!");
		return;
	}

	if (k < -1) {
		lidia_warning_handler ("quadratic_number_standard::"
				       "absolute_Ln_approximation",
				       "Increasing precision k to -1");
		k = -1;
	}

	if (a.sign() >= 0)
		s = 1;
	else
		s = -1;

	if (b.sign() < 0)
		s *= -1;

	p.assign(abs(a), abs(b), d);
	p.get_relative_approximation(r, k+4);
	log(l, r, k+3);

	// num = |x^2 - y^2 * Delta|;
	square (num, b);
	multiply (num, num, this->get_discriminant());
	square (h, a);
	subtract (num, num, h);
	num.absolute_value();

	// den = Truncate ((Truncate (z, k+8))^2, k+5);
	truncate (den, d, k+8);
	square (den, den);
	truncate (den, den, k+5);

	divide (r, num, den, k+4);
	log (r, r, k+2);

	// l = s * (l - r/2);
	r.divide_by_2();
	subtract (l, l, r);
	if (s == -1)
		l.negate ();

	// truncate l
	l.truncate (l.b_value()+k+1);
}



xbigfloat
quadratic_number_standard::get_absolute_Ln_approximation (long k) const
{
	debug_handler("quadratic_number_standard",
		      "get_absolute_Ln_approximation(long) const");
	xbigfloat l;
	this->get_absolute_Ln_approximation(l, k);
	return l;
}



void
quadratic_number_standard::get_absolute_Ln_approximation (xbigfloat & l, long k) const
{
	debug_handler("quadratic_number_standard",
		      "get_absolute_Ln_approximation(xbigfloat, long) const");

	long m, t, b1;
	xbigfloat a1, a2, a3, n, sigma_n;

	if (QO == NULL) {
		lidia_error_handler("quadratic_number_standard::"
				    "get_absolute_Ln_approximation"
				    "(xbigfloat, long) const",
				    "Quadratic order not initialized!");
		return;
	}

	// Ln is zero for rational numbers.
	//
	if (b.is_zero()) {
		l.assign_zero();
		return;
	}

	// k at least 0.
	//
	if (k < 0) {
		k = 0;
	}

	// b1 = ceil(b(Delta)/2) + b(b)
	b1 = LiDIA::b_value(this->get_discriminant());
	if (b1&1)
		b1 = (b1 >> 1) + 1;
	else
		b1 = b1 >> 1;

	b1 += LiDIA::b_value(b);

	// m = b1+max(b(a),b1)+1;
	m = LiDIA::b_value(a);
	if (m > b1)
		m += b1+1;
	else
		m = 2*b1+1;

	// abs. k+7 approximations to a+b*sqrt(D) and a-b*sqrt(D)
	//
	t = k+6;

	//std::cout << "qns::get_abs_Ln_approx:: m = " << m << std::endl;
	//std::cout << "qns::get_abs_Ln_approx:: t = " << t << std::endl;

	// relative t+m+4 approx. to b
	truncate(a1, b, t+m+5);

	// relative t+m+4 approx. to sqrt(Delta)
	sqrt(a2, this->get_discriminant(), t+m+4);

	// relative t+m+2 approx. to b*sqrt(Delta)
	multiply(a3, a1, a2);

	// relative t+m+1 approx. to b*sqrt(Delta)
	a3.truncate(t+m+4);

	// relative t+1 approx. to a+b*sqrt(Delta)
	add(n, a, a3);

	// relative t+1 approx. to a-b*sqrt(Delta)
	subtract(sigma_n, a, a3);


	// relative k+3 approx. to this/sigma(this)
	//

	// relative k+4 approx. to this/sigma(this)
	divide(n, n, sigma_n, k+6);

	// relative k+3 approx. to this/sigma(this)
	n.truncate(k+6);
	n.absolute_value();

	// absolute k+1 approx. to 1/2 ln|this/sigma(this)|
	log(l, n, k+2);
	shift_right(l, l, 1);

	// absolute k approx. to Ln(this)
	l.truncate(k+1+l.b_value());
}



//
// for debugging only
//
// compare with bigfloat computation.
// sufficient bigfloat precision required.
//

void
quadratic_number_standard::get_absolute_Ln_approximation_with_check_old (xbigfloat & l, long k) const
{
	int s;
	bigint num, h;
	quadratic_number_standard p;
	xbigfloat r, den;

	if (QO == NULL) {
		lidia_error_handler("quadratic_number_standard::absolute_Ln_approximation"
				    "_with_check(xbigfloat, long) const",
				    "Quadratic order not initialized!");
		return;
	}


	bigfloat l_big, r_big, delta_big;

	std::cout << "quadratic_number_standard::absolute_Ln_approximation_with_check (xbigfloat & l, long k)";
	std::cout << std::endl;
	std::cout << "Looking for absolute Ln precision k = " << k << std::endl;
	std::cout << "*this = " << *this << std::endl;

	if (k < -1) {
		lidia_warning_handler ("quadratic_number_standard::"
				       "absolute_Ln_approximation_with_check",
				       "Increasing precision k to -1");
		k = -1;
	}

	if (a.sign() >= 0)
		s = 1;
	else
		s = -1;

	if (b.sign() < 0)
		s *= -1;

	p.assign(abs(a), abs(b), d);

	std::cout << "s = " << s << std::endl;
	std::cout << "p = " << p << std::endl;

	// absolute k+2 approx. l to log|*this|

	p.get_relative_approximation(r, k+4);
	std::cout << "relative k+4 = " << k+4 << " approx. to p is ";
	std::cout << r << std::endl;

	log(l, r, k+3);

	// check
	sqrt(delta_big, bigfloat(p.get_discriminant()));
	multiply (r_big, p.get_b(), delta_big);
	add (r_big, r_big, p.get_a());
	divide (r_big, r_big, p.get_d());
	r_big.absolute_value();
	log (l_big, r_big);

	if (!(LiDIA::b_value(l-xbigfloat(l_big)) <= -k-2)) {
		std::cout << "ERROR: report follows:" << std::endl;
		std::cout << "p = " << p << std::endl;
		std::cout << "Required absolute precision is k+2 = " << k+2 << std::endl;
		std::cout << "log(p) as xbigfloat = " << l << std::endl;
		std::cout << "log(p) as  bigfloat = " << l_big << std::endl;
		exit(1);
	}

	// absolute k+1 approx. r to log|N(p)|

	// num = |x^2 - y^2 * Delta|;
	square (num, b);
	multiply (num, num, this->get_discriminant());
	square (h, a);
	subtract (num, num, h);
	num.absolute_value();

	// den = Truncate ((Truncate (z, k+8))^2, k+5);
	truncate (den, d, k+8);
	square (den, den);
	truncate (den, den, k+5);

	divide (r, num, den, k+4);
	log (r, r, k+2);

	// check
	square (den, den);
	divide (r_big, bigfloat(num), bigfloat(d*d));
	log (r_big, r_big);

	if (!(LiDIA::b_value(r-xbigfloat(r_big)) <= -k-1)) {
		std::cout << "ERROR: report follows:" << std::endl;
		std::cout << "p = " << p << std::endl;
		std::cout << "norm(p) = " << num << " / " << den << std::endl;
		std::cout << "Required absolute precision is k+1 = " << k+1 << std::endl;
		std::cout << "log(Np) as xbigfloat = " << r << std::endl;
		std::cout << "log(Np) as  bigfloat = " << r_big << std::endl;
		exit(1);
	}

	// l = s * (l - r/2);
	r.divide_by_2();
	subtract (l, l, r);
	if (s == -1)
		l.negate ();

	// truncate l
	l.truncate (l.b_value()+k+1);


	// check
	r_big.divide_by_2();
	subtract (l_big, l_big, r_big);
	if (s == -1)
		l_big.negate();

	if (!(LiDIA::b_value(l-xbigfloat(l_big)) <= -k)) {
		std::cout << "ERROR: report follows:" << std::endl;
		std::cout << "*this = " << *this << std::endl;
		std::cout << "Required absolute precision is k = " << k << std::endl;
		std::cout << "Ln(*this) as xbigfloat = " << l << std::endl;
		std::cout << "Ln(*this) as  bigfloat = " << l_big << std::endl;
		exit(1);
	}
}



void
quadratic_number_standard::get_absolute_Ln_approximation_with_check (xbigfloat & l, long k) const
{
	debug_handler("quadratic_number_standard",
		      "get_absolute_Ln_approximation_with_check"
		      "(xbigfloat, long) const");

	long m, t, b1;
	xbigfloat a1, a2, a3, n, sigma_n;

	if (QO == NULL) {
		lidia_error_handler("quadratic_number_standard::"
				    "get_absolute_Ln_approximation"
				    "(xbigfloat, long) const",
				    "Quadratic order not initialized!");
		return;
	}

	if (k < 0) {
		k = 0;
	}

	// b1 = ceil(b(Delta)/2) + b(b)
	b1 = LiDIA::b_value(this->get_discriminant());
	if (b1&1)
		b1 = (b1 >> 1) + 1;
	else
		b1 = b1 >> 1;

	b1 += LiDIA::b_value(b);

	// m = b1+max(b(a),b1)+1;
	m = LiDIA::b_value(a);
	if (m > b1)
		m += b1+1;
	else
		m = 2*b1+1;

	// abs. k+7 approximations to a+b*sqrt(D) and a-b*sqrt(D)
	//
	t = k+6;

	//std::cout << "qns::get_abs_Ln_approx:: m = " << m << std::endl;
	//std::cout << "qns::get_abs_Ln_approx:: t = " << t << std::endl;

	// relative t+m+4 approx. to b
	truncate(a1, b, t+m+5);

	// relative t+m+4 approx. to sqrt(Delta)
	sqrt(a2, this->get_discriminant(), t+m+4);

	// relative t+m+2 approx. to b*sqrt(Delta)
	multiply(a3, a1, a2);

	// relative t+m+1 approx. to b*sqrt(Delta)
	a3.truncate(t+m+4);

	// relative t+1 approx. to a+b*sqrt(Delta)
	add(n, a, a3);

	// relative t+1 approx. to a-b*sqrt(Delta)
	subtract(sigma_n, a, a3);


	// relative k+3 approx. to this/sigma(this)
	//

	// relative k+4 approx. to this/sigma(this)
	divide(n, n, sigma_n, k+6);

	// relative k+3 approx. to this/sigma(this)
	n.truncate(k+6);
	n.absolute_value();

	// absolute k+1 approx. to 1/2 ln|this/sigma(this)|
	log(l, n, k+2);
	shift_right(l, l, 1);

	// absolute k approx. to Ln(this)
	l.truncate(k+1+l.b_value());


	// check
	bigfloat l_big, n_big, sigma_n_big, delta_big;

	sqrt(delta_big, this->get_discriminant());
	n_big = (a + b * delta_big);
	sigma_n_big = (a - b * delta_big);

	n_big /= sigma_n_big;
	n_big.absolute_value();

	log(l_big, n_big);
	shift_right(l_big, l_big, 1);

	if (!(LiDIA::b_value(l-xbigfloat(l_big)) <= -k)) {
		std::cout << "ERROR: report follows:" << std::endl;
		std::cout << "this = " << *this << std::endl;
		std::cout << "Required absolute precision is k = " << k << std::endl;
		std::cout << "Ln as xbigfloat is l = " << l << std::endl;
		std::cout << "       l as bigfloat = "; l.print_as_bigfloat(); std::cout << std::endl;
		std::cout << "Ln as  bigfloat is l_big = " << l_big << std::endl;
		exit(1);
	}
	else {
		std::cout << "k = " << k << std::endl;
		std::cout << "Ln as xbigfloat is l = " << l << std::endl;
		std::cout << "       l as bigfloat = "; l.print_as_bigfloat(); std::cout << std::endl;
		std::cout << "Ln as  bigfloat is l_big = " << l_big << std::endl;
	}
}



//
// get_Ln_estimate(l,u)
//
// Let *this = alpha = (x+ysqrt(Delta))/z.
// This function returns l and u such that l <= |Ln alpha| <= u.
//
// Define
//             |x|/(|y|sqrt(Delta)), if |x| > |y|sqrt(Delta)
// c(alpha) =
//             |y|sqrt(Delta)/|x|,   if |x| < |y|sqrt(Delta)
// and
//
// d(alpha) = min{d in \Z | c(alpha) \geq 1 + 2^{-d}}.
//
// Then u/l \leq 8 if c(alpha) > 3
// and  u/l \leq 8+d(alpha) for c(alpha) < 3.
//

void
quadratic_number_standard::get_Ln_estimate (xbigfloat & l, xbigfloat & u) const
{
	debug_handler("quadratic_number_standard",
		      "get_Ln_estimate(xbigfloat&, xbigfloat&)");

	xbigfloat xa, xb, xD, c, L, U;
	bigint Delta, q, r, a2, b2Delta;
	bool   c_is_gt_3;
	long   k, d;

	// Ln(alpha) = 0, if x = 0 or y = 0.
	//
	if (a.is_zero() || b.is_zero()) {
		l.assign_zero();
		u.assign_zero();
		return;
	}

	// get discriminant
	Delta.assign(this->get_discriminant());

	// determine a^2 and b^2 * Delta
	//
	square(a2, a);
	square(b2Delta, b);
	multiply(b2Delta, b2Delta, Delta);

	// a^2 = q * b^2 * Delta + r
	//
	if (a2 > b2Delta) {
		div_rem (q, r, a2, b2Delta);
		d = - LiDIA::b_value(r) + LiDIA::b_value(b2Delta) + 3;
	}

	// b^2 * Delta = q * a^2 + r
	//
	else {
		div_rem (q, r, b2Delta, a2);
		d = - LiDIA::b_value(r) + LiDIA::b_value(a2) + 3;
	}

	// decide whether c > 3
	//
	c_is_gt_3 = (q >= 10) || (q == 9 && r > 0);

	// accuracy for approximating c
	//
	if (c_is_gt_3) {
		k = 2;
	}
	else {
		k = d+3;
	}

	// approximate c
	//
	if (a2 > b2Delta) {
		xa.assign(a);
		xa.absolute_value();
		xa.truncate(k+3);

		xb.assign(b);
		xb.absolute_value();
		xb.truncate(k+4);

		sqrt(xD, xbigfloat(Delta), k+4);
		multiply (xb, xb, xD);

		divide(c, xa, xb, k+3);
	}
	else {
		xa.assign(a);
		xa.absolute_value();
		xa.truncate(k+4);

		xb.assign(b);
		xb.absolute_value();
		xb.truncate(k+3);

		sqrt(xD, xbigfloat(Delta), k+3);
		multiply (xb, xb, xD);

		divide(c, xb, xa, k+3);
	}

	// approximate c+3 and c-1
	//
	add(L, c, 3);
	subtract(U, c, 1);

	// determine l and u
	//
	l.assign_one();
	shift_right(l, l, L.b_value()+1);

	u.assign_one();

	if (c_is_gt_3)
		shift_right(u, u, U.b_value()-2);
	else
		shift_right(u, u, U.b_value()-3);
}



// check correctness of Ln approximation
//

bool
quadratic_number_standard::check_Ln_correctness (lidia_size_t trials) const
{
	debug_handler("quadratic_number_standard",
		      "check_Ln_correctness(lidia_size_t) const");

	xbigfloat l1, l2;
	long k, c;
	lidia_size_t i;
	random_generator rg;

	if (this->is_rational_number()) {
		return true;
	}

	for (i = 0; i < trials; i++) {
		rg >> k;
		rg >> c;

		k %= 2000;
		if (k < 0) k = -k;

		c %= 2000;
		if (c < 0) c = -c;

		this->get_absolute_Ln_approximation(l1, k);
		this->get_absolute_Ln_approximation(l2, k+c);

		if (!xbigfloat::check_absolute_error(l1, l2, k, c)) {
			std::cout << "qns::check_Ln_correctness::Incorrect Ln approximation." << std::endl;
			return false;
		}
	}
	return true;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
