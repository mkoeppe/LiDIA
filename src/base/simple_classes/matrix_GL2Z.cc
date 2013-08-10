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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/matrix_GL2Z.h"
#include	<cstdlib>
#include	<cctype>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// c'tors
//

matrix_GL2Z::matrix_GL2Z (const bigint & a, const bigint & c,
			  const bigint & b, const bigint & d)
{
	bigint d_2, sv, tu;

	s.assign(a);
	u.assign(c);
	t.assign(b);
	v.assign(d);
	multiply(sv, s, v);
	multiply(tu, t, u);
	subtract(d_2, sv, tu);
	if (d_2.intify(determinant) || (abs(determinant) != 1)) {
		lidia_error_handler("matrix_GL2Z",
				    "matrix_GL2Z(bigint, bigint, bigint, bigint) - abs(det) != 1::not a valid GL2Z matrix");
		return;
	}
}



matrix_GL2Z::matrix_GL2Z (const base_matrix< bigint > & A)
{
	bigint d_2, sv, tu;

	if ((A.get_no_of_rows() != 2) || (A.get_no_of_columns() != 2)) {
		lidia_error_handler("matrix_GL2Z",
				    "matrix_GL2Z(base_matrix< bigint > ) - the matrix must be 2x2");
		return;
	}

	s.assign(A.member(0, 0));
	u.assign(A.member(0, 1));
	t.assign(A.member(1, 0));
	v.assign(A.member(1, 1));
	multiply(sv, s, v);
	multiply(tu, t, u);
	subtract(d_2, sv, tu);
	if (d_2.intify(determinant) || (abs(determinant) != 1)) {
		lidia_error_handler("matrix_GL2Z",
				    "matrix_GL2Z(base_matrix< bigint > ) - abs(det) != 1::not a valid GL2Z matrix");
		return;
	}
}



int
matrix_GL2Z::det () const
{
	if (determinant == 0) {
		bigint d, sv, tu;
		multiply(sv, s, v);
		multiply(tu, t, u);
		subtract(d, sv, tu);
		if (!d.intify(determinant) || abs(determinant) != 1) {
			lidia_error_handler("matrix_GL2Z", "det() - abs(det) != 1::not a valid GL2Z matrix");
			return 0;
		}
	}
	return determinant;
}



//
// assigners
//

void
matrix_GL2Z::assign (const base_matrix< bigint > & A)
{
	if ((A.get_no_of_rows() != 2) || (A.get_no_of_columns() != 2)) {
		lidia_error_handler("matrix_GL2Z", "operator = the operand "
				    "must be a 2x2 matrix");
	}

	s.assign(A.member(0, 0));
	u.assign(A.member(0, 1));
	t.assign(A.member(1, 0));
	v.assign(A.member(1, 1));

	determinant = det();
}



//
// modifiers
//

void
matrix_GL2Z::invert ()
{
	bigint d, sv, tu;

	multiply(sv, s, v);
	multiply(tu, t, u);
	subtract(d, sv, tu);
	d.intify(determinant);

	if (determinant == 0) {
		lidia_error_handler("matrix_GL2Z", "invert() - det = 0::cannot invert");
		return;
	}
	else if (determinant == 1) {
		swap(s, v);
		t.negate();
		u.negate();
	}
	else if (determinant == -1) {
		swap(s, v);
		s.negate();
		v.negate();
	}
}



//
// accessors
//

const bigint &
matrix_GL2Z::at (int i, int j) const
{
	static bigint zero(0L);

	if (i == 0 && j == 0)
		return s;
	if (i == 0 && j == 1)
		return u;
	if (i == 1 && j == 0)
		return t;
	if (i == 1 && j == 1)
		return v;

	lidia_error_handler("matrix_GL2Z", "operator () - indices not in range");
	return zero;
}



//
// arithmetic procedures
//

void
multiply (matrix_GL2Z & C, const matrix_GL2Z & A, const matrix_GL2Z & B)
{
	bigint tmp1, tmp2, cs, ct, cu, cv;

	multiply(tmp1, A.s, B.s);
	multiply(tmp2, A.u, B.t);
	add(cs, tmp1, tmp2);

	multiply(tmp1, A.s, B.u);
	multiply(tmp2, A.u, B.v);
	add(cu, tmp1, tmp2);

	multiply(tmp1, A.t, B.s);
	multiply(tmp2, A.v, B.t);
	add(ct, tmp1, tmp2);

	multiply(tmp1, A.t, B.u);
	multiply(tmp2, A.v, B.v);
	add(cv, tmp1, tmp2);

	C.s.assign(cs);
	C.t.assign(ct);
	C.u.assign(cu);
	C.v.assign(cv);

	C.determinant = A.determinant * B.determinant;
}



std::ostream &
operator << (std::ostream & out, const matrix_GL2Z & a)
{
	char ss[1024], st[1024], su[1024], sv[1024];
	int ls, lt, lu, lv, max, i;

	ls = bigint_to_string(a.s, ss);
	max = ls;
	lt = bigint_to_string(a.t, st);
	max = (max > lt) ? max : lt;
	lu = bigint_to_string(a.u, su);
	max = (max > lu) ? max : lu;
	lv = bigint_to_string(a.v, sv);
	max = (max > lv) ? max : lv;

	out << "(";
	for (i = ls; i < max; i++)
		out << " ";
	out << a.s;
	out << " ";
	for (i = lu; i < max; i++)
		out << " ";
	out << a.u << ")\n";

	out << "(";
	for (i = lt; i < max; i++)
		out << " ";
	out << a.t;
	out << " ";
	for (i = lv; i < max; i++)
		out << " ";
	out << a.v << ")\n";

	return out;
}



//
// I/O
//

std::istream &
operator >> (std::istream & in, matrix_GL2Z & A)
{
	debug_handler("matrix_GL2Z", "operator >>");
	char c;

	in >> c;
	if (c != '(') {
		lidia_error_handler("matrix_GL2Z", "operator >> - (expected");
		return in;
	}
	in >> A.s;
	in >> A.u;
	do {
		in >> c;
	} while (isspace(c));
	if (c != ')') {
		lidia_error_handler("matrix_GL2Z", "operator >> -) expected");
		return in;
	}
	do {
		in >> c;
	} while (isspace(c));
	if (c != '(') {
		lidia_error_handler("matrix_GL2Z", "operator >> - (expected");
		return in;
	}
	in >> A.t;
	in >> A.v;
	do {
		in >> c;
	} while (isspace(c));
	if (c != ')') {
		lidia_error_handler("matrix_GL2Z", "operator >> -) expected");
		return in;
	}

	bigint d_2, sv, tu;
	multiply(sv, A.s, A.v);
	multiply(tu, A.t, A.u);
	subtract(d_2, sv, tu);
	if (d_2.intify(A.determinant) || (abs(A.determinant) != 1)) {
		lidia_error_handler("matrix_GL2Z",
				    "operator >> - abs(det) != 1::not a valid GL2Z matrix");
		return in;
	}

	return in;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
