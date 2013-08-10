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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_BASE_ELLIPTIC_CURVE_REP_CC_GUARD_
#define LIDIA_BASE_ELLIPTIC_CURVE_REP_CC_GUARD_


#ifndef LIDIA_ELLIPTIC_CURVE_FLAGS_H_GUARD_
# include	"LiDIA/elliptic_curve_flags.h"
#endif
#ifndef LIDIA_POINT_OPERATIONS_H_GUARD_
# include	"LiDIA/elliptic_curves/point_operations.h"
#endif
#ifndef LIDIA_BASE_ELLIPTIC_CURVE_REP_H_GUARD_
# include	"LiDIA/elliptic_curves/base_elliptic_curve_rep.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// initialization
//
template< class T >
void
base_elliptic_curve_rep< T >::compute_invariants ()
{
	debug_handler ("base_elliptic_curve_rep< T >", "compute_invariants()");

	if (this->a6.characteristic() == 2) {
		T aa4(this->get_a4());

		square(this->b2, this->a1);
		multiply(this->b4, this->a1, this->a3);
		square(this->b6, this->a3);
		this->b8 = this->a1*(this->a1 * this->a6 + this->a3 * aa4) + this->a2 * this->a3 * this->a3 + aa4 * aa4;
		square(this->c4, this->b2);
		multiply(this->c6, this->c4, this->b2);
		this->delta = this->b2 * (this->b2 * this->b8 + this->b4 * this->b6) + this->b6 * this->b6;
	}
	else {
		this->b2 = this->a1 * this->a1 + 4  *  this->a2;
		this->b4 = this->a1 * this->a3 + 2 * this->a4;
		this->b6 = this->a3 * this->a3 + 4 * this->a6;
		this->b8 = this->a1 * (this->a1 * this->a6 - this->a3 * this->a4) + this->a2 * (4 * this->a6 +this->a3 * this->a3) - this->a4 * this->a4;
		this->c4 = this->b2 * this->b2 - 24 * this->b4;
		this->c6 = this->b2 * (36 * this->b4 - this->b2 * this->b2) - 216 * this->b6;
		this->delta = (this->c4 * this->c4 * this->c4 - this->c6 * this->c6) / 1728;
	}

	if (this->delta.is_zero())
		lidia_error_handler("base_elliptic_curve_rep< T >::compute_invariants",
				    "Discriminant is zero");
}



//
// assignment
//

template< class T >
void
base_elliptic_curve_rep< T >::assign (const base_elliptic_curve_rep< T > & e)
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "assign(const base_elliptic_curve_rep< T >");

	if (this != &e) {
		this->cp = e.cp; this->cm = e.cm; this->bcp = e.bcp;
		this->info = e.info;
		this->a1 = e.a1; this->a2 = e.a2; this->a3 = e.a3; this->a4 = e.a4; this->a6 = e.a6;
		this->b2 = e.b2; this->b4 = e.b4; this->b6 = e.b6; this->b8 = e.b8;
		this->c4 = e.c4; this->c6 = e.c6; this->delta = e.delta;
		this->output_mode = e.output_mode;
	}
}



template< class T >
void
base_elliptic_curve_rep< T >::set_coefficients (const T & x4, const T & x6,
						elliptic_curve_flags::curve_model m)
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "set_coefficients(const T& x4, const T& x6, "
		       "elliptic_curve_flags::curve_model)");

	if (x4.characteristic() == 3)
		lidia_error_handler("base_elliptic_curve< T >", "elliptic curves in"
				    " char. 3 are not yet supported");

	if (x4.characteristic() == 2) {
		// a1 = 1, a3 = 0;

		// FIXME:
		// bigint, bigrational, bigfloat, and bigcomplex don't have T::get_field()
#ifdef HAS_NO_GET_FIELD
		this->a1.assign_one();
		this->a3.assign_zero();
#else
		this->a1.assign_one(x4.get_field());
		this->a3.assign_zero(x4.get_field());
#endif

		this->a2.assign(x4);
		this->a6.assign(x6);
		this->cp = elliptic_curve_flags::GF2N_F;

		// this is a hack: for projective coordinates in GF2N we need
		// the fourth root of a6. Since I do not want to add a new variable
		// for this case, I store it in a4.

#ifndef HAS_NO_GET_FIELD
		bigint e;

		shift_left(e, bigint(1), this->a2.absolute_degree() - 2);
		power(this->a4, this->a6, e);
#endif
	}
	else {
		// a1 = a2 = a3 = 0;

		// FIXME:
		// bigint, bigrational, bigfloat, and bigcomplex don't have T::get_field()
#ifdef HAS_NO_GET_FIELD
		this->a1.assign_zero();
		this->a2.assign_zero();
		this->a3.assign_zero();
#else
		this->a1.assign_zero(x4.get_field());
		this->a2.assign_zero(x4.get_field());
		this->a3.assign_zero(x4.get_field());
#endif

		this->a4.assign(x4);
		this->a6.assign(x6);
		this->cp = elliptic_curve_flags::SHORT_W;
	}


	this->cm = m;
	this->bcp.set(this->cp, this->cm, this->a6.characteristic());

	this->output_mode = elliptic_curve_flags::SHORT;
	this->info = 0;
	compute_invariants();
}



template< class T >
void
base_elliptic_curve_rep< T >::set_coefficients (const T & x1, const T & x2, const T & x3,
						const T & x4, const T & x6,
						elliptic_curve_flags::curve_model m)
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "set_coefficients(const T&x1, const T&x2, const T&x3"
		       "const T&x4, const T&x6, "
		       "elliptic_curve_flags::curve_model)");

	if (x4.characteristic() == 3)
		lidia_error_handler("base_elliptic_curve< T >", "elliptic curves in"
				    " char. 3 are not yet supported");

	if (x4.characteristic() == 2) {
		this->a1.assign(x1);
		this->a2.assign(x2);
		this->a3.assign(x3);
		this->a4.assign(x4);
		this->a6.assign(x6);

		if (this->a1.is_one() && this->a3.is_zero() && this->a4.is_zero()) {
			this->cp = elliptic_curve_flags::GF2N_F;

#ifndef HAS_NO_GET_FIELD
			// this is a hack: for projective coordinates in GF2N we need
			// the fourth root of a6. Since I do not want to add a new variable
			// for this case, I store it in a4.

			bigint e;
			shift_left(e, bigint(1), this->a2.absolute_degree() - 2);

			power(this->a4, this->a6, e);
#endif
		}
		else {
			this->cp = elliptic_curve_flags::LONG_W;
			if (this->a1.is_zero() && this->a2.is_zero() && this->a3.is_zero())
				lidia_error_handler("base_elliptic_curve< T >", "short Weierstrass"
						    " form impossible for fields of char. 2");
		}
	}
	else {
		this->a1.assign(x1);
		this->a2.assign(x2);
		this->a3.assign(x3);
		this->a4.assign(x4);
		this->a6.assign(x6);

		if (this->a1.is_zero() && this->a2.is_zero() && this->a3.is_zero())
			this->cp = elliptic_curve_flags::SHORT_W;
		else
			this->cp = elliptic_curve_flags::LONG_W;
	}

	this->cm = m;
	this->bcp.set(this->cp, this->cm, this->a6.characteristic());
	this->output_mode = elliptic_curve_flags::SHORT;
	this->info = 0;
	compute_invariants();
}



//
// Transformation
//

template< class T >
void
base_elliptic_curve_rep< T >::transform (const T& r, const T& s, const T& t)
{
	// NB  u = 1; the more general case is not implemented here

	debug_handler ("base_elliptic_curve_rep< T >",
		       "transform(constT&, constT&, constT&)");

	if (this->a6.characteristic() == 2) {
		this->a6 += r*(this->a4 + r*(this->a2 + r)) + t*(this->a3 + r*this->a1 + t);
		this->a4 += s*this->a3 + (t + r*s)*this->a1 + r*r;
		this->a3 += r*this->a1 +t+t;
		this->a2 += s*this->a1 + r + s*s;
		this->a1 += s+s;
		if (this->a1.is_one() && this->a3.is_zero() && this->a4.is_zero())
			this->cp = elliptic_curve_flags::GF2N_F;
		else
			this->cp = elliptic_curve_flags::LONG_W;
	}
	else {
		this->a6 += r*(this->a4 + r*(this->a2 + r)) - t*(this->a3 + r*this->a1 + t);
		this->a4 += -s*this->a3 + 2*r*this->a2 - (t + r*s)*this->a1 + 3*r*r - 2*s*t;
		this->a3 += r*this->a1 + t + t;
		this->a2 += -s*this->a1 + 3*r - s*s;
		this->a1 += s+s;
		if (this->a1 == 0 && this->a2 == 0 && this->a3 == 0)
			this->cp = elliptic_curve_flags::SHORT_W;
		else
			this->cp = elliptic_curve_flags::LONG_W;
	}
	compute_invariants();
}



//
// input / output
//

template< class T >
void
base_elliptic_curve_rep< T >::set_verbose (int i)
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "set_verbose(int)");
	this->info = i;
}



template< class T >
int
base_elliptic_curve_rep< T >::verbose () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "verbose()");
	return this->info;
}



template< class T >
void
base_elliptic_curve_rep< T >::set_output_mode(long i)
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "set_output_mode(long i)");

	switch(i) {
	case 0:
		this->output_mode = elliptic_curve_flags::SHORT;
		break;
	case 1:
		this->output_mode = elliptic_curve_flags::LONG;
		break;
	case 2:
		this->output_mode = elliptic_curve_flags::PRETTY;
		break;
	case 3:
		this->output_mode = elliptic_curve_flags::TEX;
		break;
	default:
		lidia_error_handler("base_elliptic_curve_rep< T >",
				    "wrong output mode");
	}
}



template< class T >
long
base_elliptic_curve_rep< T >::get_output_mode () const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "get_output_mode()");

	return this->output_mode;
}



template< class T >
void
base_elliptic_curve_rep< T >::output_short (std::ostream & out) const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "output_short(std::ostream&)");

	if (this->cp != elliptic_curve_flags::GF2N_F)
		out << "[" << this->a1 << ", " << this->a2 << ", " << this->a3 << ", " << this->a4 << ", " << this->a6 << "]" << std::flush;
	else
		out << "[" << this->a1 << ", " << this->a2 << ", " << this->a3 << ", " << this->a3 << ", " << this->a6 << "]" << std::flush;
}



template< class T >
void
base_elliptic_curve_rep< T >::output_long (std::ostream & out) const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "output_long(std::ostream&)");

	if (this->cp != elliptic_curve_flags::GF2N_F)
		out << "[" << this->a1 << ", " << this->a2 << ", " << this->a3 << ", " << this->a4 << ", " << this->a6 << "; "
		    << this->b2 << ", " << this->b4 << ", " << this->b6 << ", " << this->b8 << "; "
		    << this->c4 << ", " << this->c6 << "; " << this->delta << "]" << std::flush;
	else
		out << "[" << this->a1 << ", " << this->a2 << ", " << this->a3 << ", " << this->a3 << ", " << this->a6 << "; "
		    << this->b2 << ", " << this->b4 << ", " << this->b6 << ", " << this->b8 << "; "
		    << this->c4 << ", " << this->c6 << "; " << this->delta << "]" << std::flush;
}



template< class T >
inline void
base_elliptic_curve_rep_ec_print (std::ostream & out, const T & a)
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "ec_print(std::ostream&, const T&)");
	out << " + " << a;
}



template <>
inline void
base_elliptic_curve_rep_ec_print< bigint > (std::ostream & out, const bigint & a)
{
	debug_handler ("base_elliptic_curve_rep< bigint >",
		       "ec_print(std::ostream&, const bigint&)");
	if (a.is_negative())
		out << " - " << abs(a);
	else
		out << " + " << a;
}



template <>
inline void
base_elliptic_curve_rep_ec_print< bigrational > (std::ostream & out, const bigrational & a)
{
	debug_handler ("base_elliptic_curve_rep< bigrational >",
		       "ec_print(std::ostream&, const bigrational&)");
	if (a.is_negative())
		out << " - " << abs(a);
	else
		out << " + " << a;
}



template< class T >
void
base_elliptic_curve_rep< T >::output_tex (std::ostream & out) const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "output_tex(std::ostream&)");

	out << "$y^2";
	if (!this->a1.is_zero()) {
		if (this->a1.is_one())
			out << " + ";
		else
			base_elliptic_curve_rep_ec_print(out, this->a1);

		out << "xy";
	}

	if (!this->a3.is_zero()) {
		if (this->a3.is_one())
			out << " + ";
		else
			base_elliptic_curve_rep_ec_print(out, this->a3);

		out << " y";
	}

	out << " = x^3";

	if (!this->a2.is_zero()) {
		if (this->a2.is_one())
			out << " + ";
		else
			base_elliptic_curve_rep_ec_print(out, this->a2);

		out << "x^2";
	}

	if (!this->a4.is_zero() && this->cp != elliptic_curve_flags::GF2N_F) {
		if (this->a4.is_one())
			out << " + ";
		else
			base_elliptic_curve_rep_ec_print(out, this->a4);

		out << "x";
	}

	if (!this->a6.is_zero()) {
		if (this->a6.is_one())
			out << " + 1";
		else
			base_elliptic_curve_rep_ec_print(out, this->a6);
	}
	out << "$" << std::flush;
}



template< class T > void
base_elliptic_curve_rep< T >::output_pretty(std::ostream & out) const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "output_pretty(std::ostream&)");

	out << "Y^2";

	if (!this->a1.is_zero()) {
		if (this->a1.is_one())
			out << " +";
		else {
			base_elliptic_curve_rep_ec_print(out, this->a1);
			out << " *";
		}
		out << " X*Y";
	}

	if (!this->a3.is_zero()) {
		if (this->a3.is_one()) out << " +";
		else {
			base_elliptic_curve_rep_ec_print(out, this->a3);
			out << " *";
		}
		out << " Y";
	}


	out << " = X^3";

	if (!this->a2.is_zero()) {
		if (this->a2.is_one())
			out << " +";
		else {
			base_elliptic_curve_rep_ec_print(out, this->a2);
			out << " *";
		}
		out << " X^2";
	}

	if (!this->a4.is_zero() && this->cp != elliptic_curve_flags::GF2N_F) {
		if (this->a4.is_one())
			out << " +";
		else {
			base_elliptic_curve_rep_ec_print(out, this->a4);
			out << " *";
		}

		out << " X";
	}

	if (!this->a6.is_zero()) {
		if (this->a6.is_one())
			out << " + 1";
		else
			base_elliptic_curve_rep_ec_print(out, this->a6);
	}

	out.flush();
}



template< class T >
void
base_elliptic_curve_rep< T >::write (std::ostream & out) const
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "operator << (std::ostream&, const base_elliptic_curve_rep< T > &)");

	switch(this->output_mode) {
	case elliptic_curve_flags::LONG:
		output_long(out);
		break;
	case elliptic_curve_flags::PRETTY:
		output_pretty(out);
		break;
	case elliptic_curve_flags::TEX:
		output_tex(out);
		break;
	default:
	case elliptic_curve_flags::SHORT:
		output_short(out);
		break;
	}
}



template< class T >
void
base_elliptic_curve_rep< T >::read (std::istream& in)
{
	debug_handler ("base_elliptic_curve_rep< T >",
		       "input(std::istream&)");

	char c; // to eat commas and detect [;
	// `[' flags input format [1, 2, 3, 4, 6]
	// seperators and terminators must then be present
	// (any nonnumeric will do after the first { or [)
	// else assumes this->a1 a2 this->a3 a4 a6 separated by whitespace
	in >> c; // This eats any whitespace first, and c is the next char input
	switch (c) {
	case '[':
		in >> this->a1 >> c;
		if (c != ',') lidia_error_handler("elliptic_curve", "input: comma expected");
		in >> this->a2 >> c;
		if (c != ',') lidia_error_handler("elliptic_curve", "input: comma expected");
		in >> this->a3 >> c;
		if (c != ',') lidia_error_handler("elliptic_curve", "input: comma expected");
		in >> this->a4 >> c;
		if (c != ',') lidia_error_handler("elliptic_curve", "input: comma expected");
		in >> this->a6 >> c;
		if (c != ']') lidia_error_handler("elliptic_curve", "input: ']' expected");
		break;
	default:
		in.putback(c);
		in >> this->a1 >> this->a2 >> this->a3 >> this->a4 >> this->a6;
	}

	this->set_coefficients(this->a1, this->a2, this->a3, this->a4, this->a6);
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_BASE_ELLIPTIC_CURVE_REP_CC_GUARD_
