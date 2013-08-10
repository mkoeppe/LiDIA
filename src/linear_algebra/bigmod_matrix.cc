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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigmod_matrix.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#ifdef MNF_DEBUG
void mnf_debug_info(const bigmod_matrix & TEMP, const bigmod_matrix & T,
                    const bigmod_matrix & IN, const bigmod_matrix & RES,
                    const bigmod_matrix &TRAFO, lidia_size_t res_j, const bigint & p)
{
	register lidia_size_t i, j;
	bigint * tmp, * TEMPtmp, *Ttmp, factor;
	std::cout << std::endl << TEMP << T << IN << RES << std::flush;
	for (i = 0; i < IN.get_no_of_columns(); i++) {
		tmp = TRAFO.value[i];
		TEMPtmp = TEMP.value[i];
		Ttmp = T.value[i];
		for (j = 0; j <= res_j; j++)
			tmp[j].assign(TEMPtmp[j]);
		for (j = res_j+1; j < IN.get_no_of_columns(); j++)
			tmp[j].assign(Ttmp[j]);
	}
	factor = bigint_matrix(TRAFO).det();
	std::cout << "Aktuelle Trafo-Matrix ist: " << TRAFO << std::flush
		  << " mit det = " << factor << " und ggT (" << p << ", " << factor
		  << ") = " << gcd(p, factor) << std::endl;
}



#endif

//
// base operations
//

void bigmod_matrix::
add_mod(bigint &c, const bigint &a, const bigint &b) const
{
	LiDIA::add(c, a, b);
	if (c.abs_compare(p) >= 0)
		LiDIA::subtract(c, c, p);
}



void bigmod_matrix::
sub_mod(bigint &c, const bigint &a, const bigint &b) const
{
	LiDIA::subtract(c, a, b);
	if (c.is_lt_zero())
		LiDIA::add(c, c, p);
}



void bigmod_matrix::
div_mod(bigint &c, bigint &d, const bigint &a, const bigint &b) const
// d returns factor of a, if we did find one, otherwise d = 1.
{
	bigint u;
	if (&d == &a) {
		bigint g = LiDIA::xgcd_left(u, b, p);
		if (!g.is_one()) {
			if (u.is_lt_zero() || !(gcd(u, g)).is_one()) {
				bigint v;
				LiDIA::divide(v, p, g);
				do {
					LiDIA::add(u, u , v);
				} while (!(gcd(u, g)).is_one());
			}
		}
		else if (u.is_lt_zero()) {
			LiDIA::add(u, u, p);
		}
//  std::cout << "While inverting "<<b<<" mod "<<p << " found factor "<< g<<std::endl;
		LiDIA::multiply(c, u, a);
		if (!p.is_zero())
			LiDIA::remainder(c, c, p);
		d = g;
	}
	else {
		d = LiDIA::xgcd_left(u, b, p);
		if (!d.is_one()) {
			if (u.is_lt_zero() || !(gcd(u, d)).is_one()) {
				bigint v;
				LiDIA::divide(v, p, d);
				do {
					LiDIA::add(u, u , v);
				} while (!(gcd(u, d)).is_one());
			}
		}
		else if (u.is_lt_zero()) {
			LiDIA::add(u, u, p);
		}
//  std::cout << "While inverting "<<b<<" mod "<<p << " found factor "<< d<<std::endl;
		LiDIA::multiply(c, u, a);
		if (!p.is_zero())
			LiDIA::remainder(c, c, p);
	}
}



void bigmod_matrix::
mult_mod(bigint &c, const bigint &a, const bigint &b) const
{
	LiDIA::multiply(c, a, b);
	if (!p.is_zero())
		LiDIA::remainder(c, c, p);
}



void bigmod_matrix::
inv_mod(bigint &c, bigint &d, const bigint &a) const
{
	d = LiDIA::xgcd_left(c, a, p);
	if (!d.is_one()) {
//    std::cout << "While inverting "<<b<<" mod "<<p << " found factor "<< d<<std::endl;
		if (c.is_lt_zero() || !(gcd(c, d)).is_one()) {
			bigint v;
			LiDIA::divide(v, p, d);
			do {
				LiDIA::add(c, c , v);
			} while (!(gcd(c, d)).is_one());
		}
	}
	else if (c.is_lt_zero()) {
		LiDIA::add(c, c, p);
	}
}



//********************************************
//** constructors ****************************
//********************************************

bigmod_matrix::bigmod_matrix():bigint_matrix(), p(0)
{
	debug_handler("bigmod_matrix", "in constructor bigmod_matrix()");
}



bigmod_matrix::bigmod_matrix(lidia_size_t r, lidia_size_t c):
	bigint_matrix(r, c), p(0)
{
	debug_handler("bigmod_matrix", "in constructor "
		      "bigmod_matrix(lidia_size_t, lidia_size_t)");
}



bigmod_matrix::bigmod_matrix(lidia_size_t r, lidia_size_t c, const bigint &mod):bigint_matrix(r, c), p(mod)
{
	debug_handler("bigmod_matrix", "in constructor "
		      "bigmod_matrix(lidia_size_t, lidia_size_t, const bigint &)");
}



bigmod_matrix::bigmod_matrix(const bigmod_matrix &M):bigint_matrix(M), p(M.p)
// copy-constructor
{
	debug_handler("bigmod_matrix", "in constructor bigmod_matrix(const bigmod_matrix &)");
}



bigmod_matrix::bigmod_matrix(const base_matrix< bigint > &M,
                             const bigint &mod): bigint_matrix(M), p(mod)
{
	debug_handler("bigmod_matrix", "in constructor "
		      "bigmod_matrix(const bigint_matrix &, const bigint &)");
	if (!p.is_zero()) {
		register lidia_size_t i, j;

		bigint * tmp;

		for (i = 0; i < rows; i++) {
			tmp = value[i];
			for (j = 0; j < columns; j++) {
				LiDIA::remainder(tmp[j], tmp[j], mod);
				if (tmp[j].is_negative())
					LiDIA::add(tmp[j], tmp[j], mod);
			}
		}
	}
}



//
// destructor
//

bigmod_matrix::~bigmod_matrix()
{
	debug_handler("bigmod_matrix", "in destructor ~bigmod_matrix()");
	// memory - deallocation is done by Destructor of base class.
}



//********************************************
//** Input / Output **************************
//********************************************

std::istream & operator >> (std::istream & in, bigmod_matrix &M)
{
	debug_handler("bigmod_matrix", "in operator >> (std::istream &, bigmod_matrix &)");
	lidia_size_t i, j;
	bigint p;
	in >> i >> j >> p;
	M.set_no_of_rows(i);
	M.set_no_of_columns(j);
	M.p.assign(p);
	bigint *Mtmp;
	if (!p.is_zero())
		for (i = 0; i < M.rows; i++) {
			Mtmp = M.value[i];
			for (j = 0; j < M.columns; j++) {
				in >> Mtmp[j];
				LiDIA::remainder(Mtmp[j], Mtmp[j], p);
				if (Mtmp[j].is_negative())
					LiDIA::add(Mtmp[j], Mtmp[j], p);
			}
		}
	else
		for (i = 0; i < M.rows; i++) {
			Mtmp = M.value[i];
			for (j = 0; j < M.columns; j++)
				in >> Mtmp[j];
		}
	return in;
}



std::ostream & operator << (std::ostream & out, const bigmod_matrix & M)
{
	debug_handler("bigmod_matrix", "in operator << (std::ostream &, const bigmod_matrix &)");
	register lidia_size_t i, j;
	bigint *Mtmp;
	for (i = 0; i < M.rows; i++) {
		Mtmp = M.value[i];
		out << "\n(" << std::flush;
		for (j = 0; j < M.columns; j++)
			out << Mtmp[j] << " " << std::flush;
		out << ")" << std::flush;
	}
	out << " mod " << M.p << "\n" << std::flush;
	return out;
}



//********************************************
//** stream handling *************************
//********************************************

void bigmod_matrix::
write_to_stream(std::ostream &out) const
{
	debug_handler("bigmod_matrix", "in member - function "
		      "write_to_stream(std::ofstream &)");
	register lidia_size_t i, j, TMP;
	bigint *tmp;

	out << rows << " " << columns << "\n" << std::flush;
	out << p << "\n" << std::flush;
	for (i = 0; i < rows; i++) {
		tmp = value[i];
		TMP = columns - 1;
		for (j = 0; j < TMP; j++)
			out << tmp[j] << " " << std::flush;
		out << tmp[j] << "\n" << std::flush;
	}
}



void bigmod_matrix::
read_from_stream(std::istream &in)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "read_from_stream(std::ifstream &)");
	lidia_size_t i, j;
	bigint *tmp;
	in >> i >> j >> p;

	set_no_of_rows(i);
	set_no_of_columns(j);
	if (!p)
		for (i = 0; i < rows; i++) {
			tmp = value[i];
			for (j = 0; j < columns; j++)
				in >> tmp[j];
		}
	else
		for (i = 0; i < rows; i++) {
			tmp = value[i];
			for (j = 0; j < columns; j++) {
				in >> tmp[j];
				LiDIA::remainder(tmp[j], tmp[j], p);
				if (tmp[j].is_negative())
					LiDIA::add(tmp[j], tmp[j], p);
			}
		}
}



//
// handling the modulus
//

void bigmod_matrix::reduce(const bigint & mod)
{
	bigint tmp;
	if (mod == p)
		return;
	if (p.is_zero()) {
		p = mod;
		bigint * tmp1;
		for (register lidia_size_t i = 0; i < rows; i++) {
			tmp1 = value[i];
			for (register lidia_size_t j = 0; j < columns; j++) {
				LiDIA::remainder(tmp1[j], tmp1[j], mod);
				if (tmp1[j].is_negative())
					LiDIA::add(tmp1[j], tmp1[j], mod);
			}
		}
		image(*this);
	}
	else {
		LiDIA::remainder(tmp, p, mod);
		if (!tmp) {
			p = mod;
			bigint * tmp1;
			for (register lidia_size_t i = 0; i < rows; i++) {
				tmp1 = value[i];
				for (register lidia_size_t j = 0; j < columns; j++)
					LiDIA::remainder(tmp1[j], tmp1[j], mod);
			}
			image(*this);
		}
		else
			lidia_error_handler("bigmod_matrix", "reduce(const bigint &)::"
					    "New modulus must divide the old one!");
	}
}



void bigmod_matrix::lift(const bigint & mod)
{
	bigint tmp;
	if (mod == p)
		return;
	LiDIA::remainder(tmp, mod, p);
	if (!tmp) {
		set_no_of_columns(rows + columns);
		for (register lidia_size_t i = rows; i; i--)
			value[rows - i][columns-i].assign(p);
		p = mod;
		image(*this);
	}
	else
		lidia_error_handler("bigmod_matrix", "lift(const bigint &)::"
				    "New modulus must be a multiple of the old one!");
}



//
// split functions
//

void bigmod_matrix::
split(bigmod_matrix &A, bigmod_matrix &B,
      bigmod_matrix &C, bigmod_matrix &D) const
{
	debug_handler("bigmod_matrix", "in member - function "
		      "split(bigmod_matrix &, bigmod_matrix &, bigmod_matrix &, bigmod_matrix &)");
	register lidia_size_t i, j, diff, diff1;
	bigint *tmp, *tmp1;
	if (A.rows > rows || A.columns > columns ||
	    B.rows > rows || B.columns > columns ||
	    C.rows > rows || C.columns > columns ||
	    D.rows > rows || D.columns > columns) {
		lidia_error_handler("bigmod_matrix",
				    "split :: Error in matrix dimension !");
		return;
	}

	if (&A != this)
		if (p.is_zero() && !A.p.is_zero())
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++) {
					LiDIA::remainder(tmp1[j], tmp[j], A.p);
					if (tmp1[j].is_negative())
						LiDIA::add(tmp1[j], tmp1[j], A.p);
				}
			}
		else if (A.p < p && !A.p.is_zero())
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++)
					LiDIA::remainder(tmp1[j], tmp[j], A.p);
			}
		else
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++)
					tmp1[j].assign(tmp[j]);
			}

	if (&B != this) {
		diff = columns - B.columns;
		if (p.is_zero() && !B.p.is_zero())
			for (i = 0; i < B.rows; i++) {
				tmp = value[i];
				tmp1 = B.value[i];
				for (j = 0; j < B.columns; j++) {
					LiDIA::remainder(tmp1[j], tmp[diff+j], B.p);
					if (tmp1[j].is_negative())
						LiDIA::add(tmp1[j], tmp1[j], B.p);
				}
			}
		else if (B.p < p && !B.p.is_zero())
			for (i = 0; i < B.rows; i++) {
				tmp = value[i];
				tmp1 = B.value[i];
				for (j = 0; j < B.columns; j++)
					LiDIA::remainder(tmp1[j], tmp[diff+j], B.p);
			}
		else
			for (i = 0; i < B.rows; i++) {
				tmp = value[i];
				tmp1 = B.value[i];
				for (j = 0; j < B.columns; j++)
					tmp1[j].assign(tmp[diff+j]);
			}
	}

	if (&C != this) {
		diff1 = rows-C.rows;
		if (p.is_zero() && !C.p.is_zero())
			for (i = 0; i < C.rows; i++) {
				tmp = value[diff1+i];
				tmp1 = C.value[i];
				for (j = 0; j < C.columns; j++) {
					LiDIA::remainder(tmp1[j], tmp[j], C.p);
					if (tmp1[j].is_negative())
						LiDIA::add(tmp1[j], tmp1[j], C.p);
				}
			}
		else if (C.p < p && !C.p.is_zero())
			for (i = 0; i < C.rows; i++) {
				tmp = value[diff1+i];
				tmp1 = C.value[i];
				for (j = 0; j < C.columns; j++)
					LiDIA::remainder(tmp1[j], tmp[j], C.p);
			}
		else
			for (i = 0; i < C.rows; i++) {
				tmp = value[diff1+i];
				tmp1 = C.value[i];
				for (j = 0; j < C.columns; j++)
					tmp1[j].assign(tmp[j]);
			}
	}
	if (&D != this) {
		diff = columns-D.columns;
		diff1 = rows-D.rows;
		if (p.is_zero() && !D.p.is_zero())
			for (i = 0; i < D.rows; i++) {
				tmp = value[diff1+i];
				tmp1 = D.value[i];
				for (j = 0; j < D.columns; j++) {
					LiDIA::remainder(tmp1[j], tmp[diff+j], D.p);
					if (tmp1[j].is_negative())
						LiDIA::add(tmp1[j], tmp1[j], D.p);
				}
			}
		else if (D.p < p && !D.p.is_zero())
			for (i = 0; i < D.rows; i++) {
				tmp = value[diff1+i];
				tmp1 = D.value[i];
				for (j = 0; j < D.columns; j++)
					LiDIA::remainder(tmp1[j], tmp[diff+j], D.p);
			}
		else
			for (i = 0; i < D.rows; i++) {
				tmp = value[diff1+i];
				tmp1 = D.value[i];
				for (j = 0; j < D.columns; j++)
					tmp1[j].assign(tmp[diff+j]);
			}
	}
}



void bigmod_matrix::
split_h(bigmod_matrix &A, bigmod_matrix &B) const
{
	debug_handler("bigmod_matrix", "in member - function "
		      "split_h(bigmod_matrix&, bigmod_matrix&)");

	register lidia_size_t i, j, diff;
	bigint *tmp, *tmp1;
	if (A.rows > rows || A.columns > columns ||
	    B.rows > rows || B.columns > columns) {
		lidia_error_handler("bigmod_matrix",
				    "split_h :: Error in matrix dimension !");
		return;
	}
	if (&A != this)
		if (p.is_zero() && !A.p.is_zero())
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++) {
					LiDIA::remainder(tmp1[j], tmp[j], A.p);
					if (tmp1[j].is_negative())
						LiDIA::add(tmp1[j], tmp1[j], A.p);
				}
			}
		else if (A.p < p && !A.p.is_zero())
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++)
					LiDIA::remainder(tmp1[j], tmp[j], A.p);
			}
		else
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++)
					tmp1[j].assign(tmp[j]);
			}

	if (&B != this) {
		diff = columns - B.columns;
		if (p.is_zero() && !B.p.is_zero())
			for (i = 0; i < B.rows; i++) {
				tmp = value[i];
				tmp1 = B.value[i];
				for (j = 0; j < B.columns; j++) {
					LiDIA::remainder(tmp1[j], tmp[diff+j], B.p);
					if (tmp1[j].is_negative())
						LiDIA::add(tmp1[j], tmp1[j], B.p);
				}
			}
		else if (B.p < p && !B.p.is_zero())
			for (i = 0; i < B.rows; i++) {
				tmp = value[i];
				tmp1 = B.value[i];
				for (j = 0; j < B.columns; j++)
					LiDIA::remainder(tmp1[j], tmp[diff+j], B.p);
			}
		else
			for (i = 0; i < B.rows; i++) {
				tmp = value[i];
				tmp1 = B.value[i];
				for (j = 0; j < B.columns; j++)
					tmp1[j].assign(tmp[diff+j]);
			}
	}
}



void bigmod_matrix::
split_v(bigmod_matrix &A, bigmod_matrix &B) const
{
	debug_handler("bigmod_matrix", "in member - function "
		      "split_v(bigmod_matrix &, bigmod_matrix &)");
	register lidia_size_t i, j, diff;
	bigint *tmp, *tmp1;
	if (A.rows > rows || A.columns > columns ||
	    B.rows > rows || B.columns > columns) {
		lidia_error_handler("bigmod_matrix",
				    "split_v :: Error in matrix dimension !");
		return;
	}
	if (&A != this)
		if (p.is_zero() && !A.p.is_zero())
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++) {
					LiDIA::remainder(tmp1[j], tmp[j], A.p);
					if (tmp1[j].is_negative())
						LiDIA::add(tmp1[j], tmp1[j], A.p);
				}
			}
		else if (A.p < p && !A.p.is_zero())
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++)
					LiDIA::remainder(tmp1[j], tmp[j], A.p);
			}
		else
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++)
					tmp1[j].assign(tmp[j]);
			}

	if (&B != this) {
		diff = rows-B.rows;
		if (p.is_zero() && !B.p.is_zero())
			for (i = 0; i < B.rows; i++) {
				tmp = value[diff+i];
				tmp1 = B.value[i];
				for (j = 0; j < B.columns; j++) {
					LiDIA::remainder(tmp1[j], tmp[j], B.p);
					if (tmp1[j].is_negative())
						LiDIA::add(tmp1[j], tmp1[j], B.p);
				}
			}
		else if (B.p < p && !B.p.is_zero())
			for (i = 0; i < B.rows; i++) {
				tmp = value[diff+i];
				tmp1 = B.value[i];
				for (j = 0; j < B.columns; j++)
					LiDIA::remainder(tmp1[j], tmp[j], B.p);
			}
		else
			for (i = 0; i < B.rows; i++) {
				tmp = value[diff+i];
				tmp1 = B.value[i];
				for (j = 0; j < B.columns; j++)
					tmp1[j].assign(tmp[j]);
			}
	}
}



void bigmod_matrix::
compose(const bigmod_matrix &A, const bigmod_matrix &B, const bigmod_matrix &C, const bigmod_matrix &D)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "compose(const bigmod_matrix&, const bigmod_matrix&, "
		      "const bigmod_matrix&, const bigmod_matrix&)");
	register lidia_size_t i, j, diff, diff1;
	bigint *tmp, *tmp1;
	if (A.rows > rows || A.columns > columns ||
	    B.rows > rows || B.columns > columns ||
	    C.rows > rows || C.columns > columns ||
	    D.rows > rows || D.columns > columns) {
		lidia_error_handler("bigmod_matrix",
				    "compose :: Error in matrix dimension !");
		return;
	}

	if (&D == this)
		return;
	if (&C == this)
		goto label1;
	if (&B == this)
		goto label2;
	if (&A != this) {
		if (A.p.is_zero() && !p.is_zero())
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++) {
					LiDIA::remainder(tmp[j], tmp1[j], p);
					if (tmp[j].is_negative())
						LiDIA::add(tmp[j], tmp[j], p);
				}
			}
		else if (p < A.p && !p.is_zero())
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++)
					LiDIA::remainder(tmp[j], tmp1[j], p);
			}
		else
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++)
					tmp[j].assign(tmp1[j]);
			}
	}

	diff = columns - B.columns;
	if (B.p.is_zero() && !p.is_zero())
		for (i = 0; i < B.rows; i++) {
			tmp = value[i];
			tmp1 = B.value[i];
			for (j = 0; j < B.columns; j++) {
				LiDIA::remainder(tmp[diff+j], tmp1[j], p);
				if (tmp[diff+j].is_negative())
					LiDIA::add(tmp[diff+j], tmp[diff+j], p);
			}
		}
	else if (p < B.p && !p.is_zero())
		for (i = 0; i < B.rows; i++) {
			tmp = value[i];
			tmp1 = B.value[i];
			for (j = 0; j < B.columns; j++)
				LiDIA::remainder(tmp[diff+j], tmp1[j], p);
		}
	else
		for (i = 0; i < B.rows; i++) {
			tmp = value[i];
			tmp1 = B.value[i];
			for (j = 0; j < B.columns; j++)
				tmp[diff+j].assign(tmp1[j]);
		}

 label2: diff1 = rows-C.rows;
	if (C.p.is_zero() && !p.is_zero())
		for (i = 0; i < C.rows; i++) {
			tmp = value[diff1+i];
			tmp1 = C.value[i];
			for (j = 0; j < C.columns; j++) {
				LiDIA::remainder(tmp[j], tmp1[j], p);
				if (tmp[j].is_negative())
					LiDIA::add(tmp[j], tmp[j], p);
			}
		}
	else if (p < C.p && !p.is_zero())
		for (i = 0; i < C.rows; i++) {
			tmp = value[diff1+i];
			tmp1 = C.value[i];
			for (j = 0; j < C.columns; j++)
				LiDIA::remainder(tmp[j], tmp1[j], p);
		}
	else
		for (i = 0; i < C.rows; i++) {
			tmp = value[diff1+i];
			tmp1 = C.value[i];
			for (j = 0; j < C.columns; j++)
				tmp[j].assign(tmp1[j]);
		}

 label1: diff = columns-D.columns;
	diff1 = rows-D.rows;
	if (D.p.is_zero() && !p.is_zero())
		for (i = 0; i < D.rows; i++) {
			tmp = value[diff1+i];
			tmp1 = D.value[i];
			for (j = 0; j < D.columns; j++) {
				LiDIA::remainder(tmp[diff+j], tmp1[j], p);
				if (tmp[diff+j].is_negative())
					LiDIA::add(tmp[diff+j], tmp[diff+j], p);
			}
		}
	else if (p < D.p && !p.is_zero())
		for (i = 0; i < D.rows; i++) {
			tmp = value[diff1+i];
			tmp1 = D.value[i];
			for (j = 0; j < B.columns; j++)
				LiDIA::remainder(tmp[diff+j], tmp1[j], p);
		}
	else
		for (i = 0; i < D.rows; i++) {
			tmp = value[diff1+i];
			tmp1 = D.value[i];
			for (j = 0; j < B.columns; j++)
				tmp[diff+j].assign(tmp1[j]);
		}

}



void bigmod_matrix::
compose_h(const bigmod_matrix &A, const bigmod_matrix &B)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "compose_h(const bigmod_matrix &, const bigmod_matrix &)");
	register lidia_size_t i, j, diff;
	bigint *tmp, *tmp1;
	if (A.rows > rows || A.columns > columns ||
	    B.rows > rows || B.columns > columns) {
		lidia_error_handler("bigmod_matrix",
				    "compose_h :: Error in matrix dimension !");
		return;
	}
	if (&B == this)
		return;
	if (&A != this) {
		if (A.p.is_zero() && !p.is_zero())
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++) {
					LiDIA::remainder(tmp[j], tmp1[j], p);
					if (tmp[j].is_negative())
						LiDIA::add(tmp[j], tmp[j], p);
				}
			}
		else if (p < A.p && !p.is_zero())
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++)
					LiDIA::remainder(tmp[j], tmp1[j], p);
			}
		else
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++)
					tmp[j].assign(tmp1[j]);
			}
	}

	diff = columns - B.columns;
	if (B.p.is_zero() && !p.is_zero())
		for (i = 0; i < B.rows; i++) {
			tmp = value[i];
			tmp1 = B.value[i];
			for (j = 0; j < B.columns; j++) {
				LiDIA::remainder(tmp[diff+j], tmp1[j], p);
				if (tmp[diff+j].is_negative())
					LiDIA::add(tmp[diff+j], tmp[diff+j], p);
			}
		}
	else if (p < B.p && !p.is_zero())
		for (i = 0; i < B.rows; i++) {
			tmp = value[i];
			tmp1 = B.value[i];
			for (j = 0; j < B.columns; j++)
				LiDIA::remainder(tmp[diff+j], tmp1[j], p);
		}
	else
		for (i = 0; i < B.rows; i++) {
			tmp = value[i];
			tmp1 = B.value[i];
			for (j = 0; j < B.columns; j++)
				tmp[diff+j].assign(tmp1[j]);
		}
}



void bigmod_matrix::
compose_v(const bigmod_matrix &A, const bigmod_matrix &B)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "compose_v(const bigmod_matrix &, const bigmod_matrix &)");
	register lidia_size_t i, j, diff;
	bigint *tmp, *tmp1;
	if (A.rows > rows || A.columns > columns ||
	    B.rows > rows || B.columns > columns) {
		lidia_error_handler("bigmod_matrix",
				    "compose_v :: Error in matrix dimension !");
		return;
	}
	if (&B == this)
		return;
	if (&A != this) {
		if (A.p.is_zero() && !p.is_zero())
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++) {
					LiDIA::remainder(tmp[j], tmp1[j], p);
					if (tmp[j].is_negative())
						LiDIA::add(tmp[j], tmp[j], p);
				}
			}
		else if (p < A.p && !p.is_zero())
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++)
					LiDIA::remainder(tmp[j], tmp1[j], p);
			}
		else
			for (i = 0; i < A.rows; i++) {
				tmp = value[i];
				tmp1 = A.value[i];
				for (j = 0; j < A.columns; j++)
					tmp[j].assign(tmp1[j]);
			}
	}

	diff = rows-B.rows;
	if (B.p.is_zero() && !p.is_zero())
		for (i = 0; i < B.rows; i++) {
			tmp = value[diff+i];
			tmp1 = B.value[i];
			for (j = 0; j < B.columns; j++) {
				LiDIA::remainder(tmp[j], tmp1[j], p);
				if (tmp[j].is_negative())
					LiDIA::add(tmp[j], tmp[j], p);
			}
		}
	else if (p < B.p && !p.is_zero())
		for (i = 0; i < B.rows; i++) {
			tmp = value[diff+i];
			tmp1 = B.value[i];
			for (j = 0; j < B.columns; j++)
				LiDIA::remainder(tmp[j], tmp1[j], p);
		}
	else
		for (i = 0; i < B.rows; i++) {
			tmp = value[diff+i];
			tmp1 = B.value[i];
			for (j = 0; j < B.columns; j++)
				tmp[j].assign(tmp1[j]);
		}
}



//********************************************
//** Component Input / Output ****************
//********************************************

void bigmod_matrix::
sto(lidia_size_t x, lidia_size_t y, const bigint &e)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "sto(lidia_size_t, lidia_size_t, const bigint &)");
	if (x< 0 || x >= rows || y< 0 || y >= columns) {
		lidia_error_handler("bigmod_matrix",
				    "sto :: Parameter out of range !");
		return;
	}
	bigint tmp;
	if (!p.is_zero()) {
		LiDIA::remainder(tmp, e, p);
		if (tmp.is_negative())
			LiDIA::add(tmp, tmp, p);
		value[x][y].assign(tmp);
	}
	else
		value[x][y].assign(e);
}



void bigmod_matrix::
sto_column(const bigint *v, lidia_size_t l, lidia_size_t j, lidia_size_t from)
{
	debug_handler("Matrix", "in member - function "
		      "sto_column(const bigint *, lidia_size_t, lidia_size_t)");
	register lidia_size_t k;
	if (j >= columns || rows < from + l || from < 0 || l < 0 || j < 0) {
		lidia_error_handler("bigmod_matrix",
				    "sto_column :: Parameter out of range !");
		return;
	}
	bigint tmp;
	if (!p.is_zero())
		for (k = 0; k < l; k++) {
			LiDIA::remainder(tmp, v[k], p);
			if (tmp.is_negative())
				LiDIA::add(tmp, tmp, p);
			value[k+from][j].assign(tmp);
		}
	else
		for (k = 0; k < l; k++)
			value[k+from][j].assign(v[k]);
}



void bigmod_matrix::
sto_column_vector(const base_vector< bigint > &v, lidia_size_t l,
		  lidia_size_t j, lidia_size_t from)
{
	debug_handler("Matrix", "in member - function "
		      "sto_column_vector(const base_vector< bigint > &, lidia_size_t, lidia_size_t, lidia_size_t)");
	register lidia_size_t k;
	if (j >= columns || rows < from + l || from < 0 || l < 0 || j < 0) {
		lidia_error_handler("bigmod_matrix",
				    "sto_column :: Parameter out of range !");
		return;
	}
	bigint tmp;
	if (!p.is_zero())
		for (k = 0; k < l; k++) {
			LiDIA::remainder(tmp, v[k], p);
			if (tmp.is_negative())
				LiDIA::add(tmp, tmp, p);
			value[k+from][j].assign(tmp);
		}
	else
		for (k = 0; k < l; k++)
			value[k+from][j].assign(v[k]);
}



void bigmod_matrix::
sto_row(const bigint *v, lidia_size_t l, lidia_size_t i, lidia_size_t from)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "sto_row(const bigint *, lidia_size_t, lidia_size_t)");
	register lidia_size_t k;
	if (i >= rows || from + l > columns || from < 0 || l < 0 || i < 0) {
		lidia_error_handler("bigmod_matrix",
				    "sto_row :: Parameter out of range !");
		return;
	}
	bigint *tmp, tmp1;
	tmp = value[i];
	if (!p.is_zero())
		for (k = 0; k < l; k++) {
			LiDIA::remainder(tmp1, v[k], p);
			if (tmp1.is_negative())
				LiDIA::add(tmp1, tmp1, p);
			tmp[k+from].assign(tmp1);
		}
	else
		for (k = 0; k < l; k++)
			tmp[k+from].assign(v[k]);
}



void bigmod_matrix::
sto_row_vector(const base_vector< bigint > &v, lidia_size_t l,
	       lidia_size_t i, lidia_size_t from)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "sto_row(const bigint *, lidia_size_t, lidia_size_t)");
	register lidia_size_t k;
	if (i >= rows || from + l > columns || from < 0 || l < 0 || i < 0) {
		lidia_error_handler("bigmod_matrix",
				    "sto_row :: Parameter out of range !");
		return;
	}
	bigint *tmp, tmp1;
	tmp = value[i];
	if (!p.is_zero())
		for (k = 0; k < l; k++) {
			LiDIA::remainder(tmp1, v[k], p);
			if (tmp1.is_negative())
				LiDIA::add(tmp1, tmp1, p);
			tmp[k+from].assign(tmp1);
		}
	else
		for (k = 0; k < l; k++)
			tmp[k+from].assign(v[k]);
}



//********************************************
//** exchange functions **********************
//********************************************

void
swap(bigmod_matrix &A, bigmod_matrix &B)
{
	//
	// DESCRIPTION: swap(A, B) exchanges matrix A with matrix B.
	// VERSION: 1.62
	//

	debug_handler("bigmod_matrix", "in function "
		      "swap(bigmod_matrix &, bigmod_matrix &)");
	lidia_size_t TMP;
	bigint **tmp;

	// swap no_of_columns
	TMP = A.columns;
	A.columns = B.columns;
	B.columns = TMP;

	// swap no_of_rows
	TMP = A.rows;
	A.rows = B.rows;
	B.rows = TMP;

	// swap values
	tmp = A.value;
	A.value = B.value;
	B.value = tmp;


	// swap moduli (TPf)
	swap(A.p, B.p);
}



//
// BEGIN: matrix arithmetic
//

//
// procedures
//

void
add(bigmod_matrix &RES, const bigmod_matrix &M, const bigmod_matrix &N)
{
	debug_handler("Matrix", "in function "
		      "add(bigmod_matrix &, const bigmod_matrix &, const bigmod_matrix &)");

	register lidia_size_t i, j;
	if (M.rows != N.rows || M.columns != N.columns) {
		lidia_error_handler("bigmod_matrix", "add :: Error in matrix dimensions");
		return;
	}
	if (M.p != N.p) {
		lidia_error_handler("bigmod_matrix", "add :: different moduli");
		return;
	}
	RES.p = M.p;
	if (RES.rows != N.rows)
		RES.set_no_of_rows(N.rows);
	if (RES.columns != N.columns)
		RES.set_no_of_columns(N.columns);
	bigint *Mtmp, *Ntmp, *REStmp;
	for (i = 0; i < N.rows; i++) {
		Ntmp = N.value[i];
		Mtmp = M.value[i];
		REStmp = RES.value[i];
		for (j = 0; j < N.columns; j++)
			RES.add_mod(REStmp[j], Ntmp[j], Mtmp[j]);
	}
}



void
add(bigmod_matrix &RES, const bigmod_matrix &M, const bigint &a)
{
	debug_handler("bigmod_matrix", "in function "
		      "add(bigmod_matrix &, const bigmod_matrix &, const bigint &)");
	if (RES.rows != M.rows)
		RES.set_no_of_rows(M.rows);
	if (RES.columns != M.columns)
		RES.set_no_of_columns(M.columns);
	RES.p = M.p;

	register lidia_size_t i, j;
	bigint *REStmp, *Mtmp;
	for (i = 0; i < RES.rows; i++) {
		REStmp = RES.value[i];
		Mtmp = M.value[i];
		for (j = 0; j < RES.columns; j++)
			RES.add_mod(REStmp[j], Mtmp[j], a);
	}
}



void
negate(bigmod_matrix & A, const bigmod_matrix & B)
{
	debug_handler("bigmod_matrix", "in function negate()");
	register lidia_size_t i, j;
	if (A.rows != B.rows)
		A.set_no_of_rows(B.rows);
	if (A.columns != B.columns)
		A.set_no_of_columns(B.columns);
	A.p.assign(B.p);

	bigint *Btmp, *Atmp;
	for (i = 0; i < B.rows; i++) {
		Btmp = B.value[i];
		Atmp = A.value[i];
		for (j = 0; j < B.columns; j++)
			if (!Btmp[j].is_zero())
				LiDIA::add(Atmp[j], -Btmp[j], B.p);
	}
}



void
subtract(bigmod_matrix &RES, const bigmod_matrix &M, const bigmod_matrix &N)
{
	debug_handler("bigmod_matrix", "in function "
		      "subtract(bigmod_matrix &, const bigmod_matrix &, const bigmod_matrix &)");
	register lidia_size_t i, j;
	if (M.rows != N.rows || M.columns != N.columns) {
		lidia_error_handler("bigmod_matrix", "subtract :: Error in matrix dimensions !");
		return;
	}
	if (M.p != N.p) {
		lidia_error_handler("bigmod_matrix", "subtract :: different moduli");
		return;
	}
	RES.p = M.p;
	if (RES.rows != N.rows)
		RES.set_no_of_rows(N.rows);
	if (RES.columns != N.columns)
		RES.set_no_of_columns(N.columns);

	bigint *REStmp, *Mtmp, *Ntmp;
	for (i = 0; i < N.rows; i++) {
		REStmp = RES.value[i];
		Mtmp = M.value[i];
		Ntmp = N.value[i];
		for (j = 0; j < N.columns; j++)
			RES.sub_mod(REStmp[j], Mtmp[j], Ntmp[j]);
	}
}



void
subtract(bigmod_matrix &RES, const bigmod_matrix &M, const bigint &a)
{
	debug_handler("bigmod_matrix", "in function "
		      "subtract(bigmod_matrix &, const bigmod_matrix &, const bigint &)");
	if (RES.rows != M.rows)
		RES.set_no_of_rows(M.rows);
	if (RES.columns != M.columns)
		RES.set_no_of_columns(M.columns);
	RES.p = M.p;

	register lidia_size_t i, j;
	bigint *REStmp, *Mtmp;
	for (i = 0; i < RES.rows; i++) {
		REStmp = RES.value[i];
		Mtmp = M.value[i];
		for (j = 0; j < RES.columns; j++)
			RES.sub_mod(REStmp[j], Mtmp[j], a);
	}
}



void
multiply(bigmod_matrix &RES, const bigmod_matrix &A, const bigmod_matrix &B)
{
	debug_handler("bigmod_matrix", "in function "
		      "multiply(bigmod_matrix &, const bigmod_matrix &, const bigmod_matrix &)");

	if (A.columns != B.rows) {
		lidia_error_handler("bigmod_matrix", "multiply :: Error in matrix dimensions ");
		return;
	}
	if (A.p != B.p) {
		lidia_error_handler("bigmod_matrix", "multiply :: different moduli");
		return;
	}
	RES.p = A.p;

	register lidia_size_t j, i, z;
	bigint tmp, *Atmp, *REStmp, tmp1;

	if (A.value != RES.value && B.value != RES.value) {
		if (RES.rows != A.rows)
			RES.set_no_of_rows(A.rows);
		if (RES.columns != B.columns)
			RES.set_no_of_columns(B.columns);
		for (j = 0; j < A.rows; j++) {
			Atmp = A.value[j];
			REStmp = RES.value[j];
			for (i = 0; i < B.columns; i++) {
				tmp.assign_zero();
				for (z = 0; z < B.rows; z++) {
					RES.mult_mod(tmp1, Atmp[z], B.value[z][i]);
					RES.add_mod(tmp, tmp, tmp1);
				}
				REStmp[i] = tmp;
			}
		}
	}
	else {
		bigmod_matrix RES1(A.rows, B.columns, A.p);
		for (j = 0; j < A.rows; j++) {
			Atmp = A.value[j];
			REStmp = RES1.value[j];
			for (i = 0; i < B.columns; i++) {
				tmp.assign_zero();
				for (z = 0; z < B.rows; z++) {
					RES1.mult_mod(tmp1, Atmp[z], B.value[z][i]);
					RES1.add_mod(tmp, tmp, tmp1);
				}
				REStmp[i].assign(tmp);
			}
		}
		if (RES.rows != A.rows)
			RES.set_no_of_rows(A.rows);
		if (RES.columns != B.columns)
			RES.set_no_of_columns(B.columns);
		RES.assign(RES1);
	}
}



void
multiply(bigmod_matrix & RES, const bigmod_matrix &A, const bigint & k)
{
	debug_handler("bigmod_matrix", "in function "
		      "multiply(bigmod_matrix &, const bigmod_matrix &, const bigint &)");
	if (k.is_negative()) {
		multiply(RES, A, -k);
		negate(RES, RES);
		return;
	}
	if (RES.rows != A.rows)
		RES.set_no_of_rows(A.rows);
	if (RES.columns != A.columns)
		RES.set_no_of_columns(A.columns);

	register lidia_size_t j, i;
	bigint *REStmp, *Atmp;
	for (j = 0; j < A.rows; j++) {
		REStmp = RES.value[j];
		Atmp = A.value[j];
		for (i = 0; i < A.columns; i++)
			RES.mult_mod(REStmp[i], Atmp[i], k);
	}
}



void
multiply(bigint *c, const bigmod_matrix &A, const bigint *v)
{
	debug_handler("bigmod_matrix", "in function "
		      "multiply(bigint *, const bigmod_matrix &, const bigint *)");
	register lidia_size_t i, j;
	bigint tmp, tmp1, *tmp2;
	for (i = 0; i < A.rows; i++) {
		tmp.assign_zero();
		tmp2 = A.value[i];
		for (j = 0; j < A.columns; j++) {
			A.mult_mod(tmp1, tmp2[j], v[j]);
			A.add_mod(tmp, tmp, tmp1);
		}
		c[i].assign(tmp);
	}
}



bool divide(bigmod_matrix &C, const bigmod_matrix &A, const bigint &n)
{
	debug_handler("bigmod_matrix", "in function "
		      "divides(const bigint &, const bigmod_matrix &)");
	if (C.rows != A.rows)
		C.set_no_of_rows(A.rows);
	if (C.columns != A.columns)
		C.set_no_of_columns(A.columns);

	register lidia_size_t i, j;
	bigint tmp, *tmp1, *tmp2;
	LiDIA::div_rem(C.p, tmp, A.p, n);
	if (!tmp.is_zero()) return false;
	for (i = 0; i < A.rows; i++) {
		tmp2 = A.value[i];
		tmp1 = C.value[i];
		for (j = 0; j < A.columns; j++) {
			if (tmp2[j].is_zero())
				tmp1[j].assign_zero();
			else {
				LiDIA::div_rem(tmp1[j], tmp, tmp2[j], n);
				if (!tmp.is_zero())
					return false;
			}
		}
	}
	return true;
}



//
// Operators
//

bigmod_matrix & bigmod_matrix::
operator += (const bigmod_matrix &M)
{
	LiDIA::add(*this, *this, M); return *this;
}



bigmod_matrix & bigmod_matrix::
operator += (const bigint &a)
{
	LiDIA::add(*this, *this, a); return *this;
}



bigmod_matrix bigmod_matrix::
operator + (const bigmod_matrix &M) const
{
	bigmod_matrix RES(rows, columns, p);
	LiDIA::add(RES, *this, M); return RES;
}



bigmod_matrix bigmod_matrix::
operator + (const bigint &a) const
{
	bigmod_matrix RES(rows, columns, p);
	LiDIA::add(RES, *this, a); return RES;
}



bigmod_matrix bigmod_matrix::
operator - () const
{
	bigmod_matrix A(rows, columns, p);
	LiDIA::negate(A, *this);
	return A;
}



bigmod_matrix & bigmod_matrix::
operator -= (const bigmod_matrix &M)
{
	LiDIA::subtract(*this, *this, M);
	return *this;
}



bigmod_matrix & bigmod_matrix::
operator -= (const bigint &a)
{
	LiDIA::subtract(*this, *this, a);
	return *this;
}



bigmod_matrix bigmod_matrix::
operator - (const bigmod_matrix &M) const
{
	bigmod_matrix RES(rows, columns, p);
	LiDIA::subtract(RES, *this, M);
	return RES;
}



bigmod_matrix bigmod_matrix::
operator - (const bigint &a) const
{
	bigmod_matrix RES(rows, columns, p);
	LiDIA::subtract(RES, *this, a);
	return RES;
}



bigmod_matrix &bigmod_matrix::
operator *= (const bigmod_matrix &A)
{
	LiDIA::multiply(*this, *this, A);
	return *this;
}



bigmod_matrix & bigmod_matrix::
operator *= (const bigint &m)
{
	LiDIA::multiply(*this, *this, m);
	return *this;
}



bigmod_matrix bigmod_matrix::
operator * (const bigmod_matrix & M) const
{
	bigmod_matrix RES(rows, M.columns, p);
	LiDIA::multiply(RES, *this, M);
	return RES;
}



bigmod_matrix bigmod_matrix::
operator * (const bigint &m) const
{
	bigmod_matrix B(rows, columns, p);
	LiDIA::multiply(B, *this, m);
	return B;
}



bigint * bigmod_matrix::
operator * (const bigint *v) const
{
	bigint *b = new bigint[rows];
	memory_handler(b, "bigmod_matrix", "operator * :: "
		       "Error in memory allocation (b)");
	LiDIA::multiply(b, *this, v);
	return b;
}



//********************************************
//** assign operator *************************
//********************************************

bigmod_matrix & bigmod_matrix::
operator = (const bigmod_matrix &M)
{
	debug_handler("bigmod_matrix", "in operator = (const bigmod_matrix &)");
	assign(M);
	return *this;
}



void bigmod_matrix::
assign(const bigmod_matrix &M)
{
	debug_handler("bigmod_matrix", "in member - function assign(const bigmod_matrix &)");
	register lidia_size_t i, j;
	if (rows != M.rows)
		set_no_of_rows(M.rows);
	if (columns != M.columns)
		set_no_of_columns(M.columns);
	bigint *tmp, *Mtmp;
	p = M.p;
	for (i = 0; i < rows; i++) {
		tmp = value[i];
		Mtmp = M.value[i];
		for (j = 0; j < columns; j++)
			tmp[j].assign(Mtmp[j]);
	}
}



void
assign(bigmod_matrix &RES, const bigmod_matrix &M)
{
	debug_handler("bigmod_matrix", "in function "
		      "assign(bigmod_matrix &, const bigmod_matrix &)");
	RES.assign(M);
}



void bigmod_matrix::
assign(const base_matrix< bigint > &M, const bigint &mod)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "assign(const base_matrix< bigint > &, const bigint &)");

	register lidia_size_t i, j;
	if (rows != M.get_no_of_rows())
		set_no_of_rows(M.get_no_of_rows());
	if (columns != M.get_no_of_columns())
		set_no_of_columns(M.get_no_of_columns());
	p = mod;

	bigint *tmp, tmp1;
	if (!p)
		for (i = 0; i < rows; i++) {
			tmp = value[i];
			for (j = 0; j < columns; j++) {
				tmp[j].assign(M.member(i, j));
			}
		}
	else
		for (i = 0; i < rows; i++) {
			tmp = value[i];
			for (j = 0; j < columns; j++) {
				LiDIA::remainder(tmp1, M.member(i, j), p);
				if (tmp1.is_negative())
					LiDIA::add(tmp1, tmp1, p);
				tmp[j].assign(tmp1);
			}
		}
}



void
assign(bigmod_matrix &RES, const base_matrix< bigint > &M, const bigint &mod)
{
	debug_handler("bigmod_matrix", "in function "
		      "assign(bigmod_matrix &, "
		      "const base_matrix< bigint > &, const bigint &)");
	RES.assign(M, mod);
}



//********************************************
//** boolean operators ***********************
//********************************************

bool bigmod_matrix::
operator == (const bigmod_matrix &N) const
{
	debug_handler("bigmod_matrix", "in operator == (const bigmod_matrix &)");
	return equal(N);
}



bool bigmod_matrix::
equal(const bigmod_matrix &N) const
{
	debug_handler("bigmod_matrix", "in member - function "
		      "equal(const bigmod_matrix &)");
	register lidia_size_t i, j;
	bigint *tmp, *Ntmp;
	if (rows != N.rows)
		return false;
	if (columns != N.columns)
		return false;
	if (p != N.p)
		return false;
	for (i = 0; i < rows; i++) {
		tmp = value[i];
		Ntmp = N.value[i];
		for (j = 0; j < columns; j++)
			if (tmp[j] != Ntmp[j])
				return false;
	}
	return true;
}



bool
equal(const bigmod_matrix &A, const bigmod_matrix &B)
{
	debug_handler("bigmod_matrix", "in function "
		      "equal(const bigmod_matrix &, const bigmod_matrix &)");
	return A.equal(B);
}



bool bigmod_matrix::
operator != (const bigmod_matrix &N) const
{
	debug_handler("bigmod_matrix", "in operator != (const bigmod_matrix &)");
	return !(equal(N));
}



bool bigmod_matrix::
unequal(const bigmod_matrix &A) const
{
	debug_handler("bigmod_matrix", "in member - function "
		      "unequal(const bigmod_matrix &)");
	return !(equal(A));
}



bool
unequal(const bigmod_matrix &A, const bigmod_matrix &B)
{
	debug_handler("bigmod_matrix", "in function "
		      "unequal(const bigmod_matrix &, const bigmod_matrix &)");
	return !(A.equal(B));
}



//
// randomize
//

void bigmod_matrix::
randomize(const bigint & S)
{
	//
	// DESCRIPTION: RES.random(S);
	// =  > 0 <= RES.value[i][j] <= S, i = 0, ..., RES.rows-1,
	//                 j = 0, ..., RES.columns-1
	// =  > p = S
	// VERSION: 1.62
	//

	debug_handler("bigmod_matrix", "in member - function "
		      "randomize(const bigint &)");
	register lidia_size_t i, j;
	bigint *tmp;
	p.assign(S);
	bigint::seed(S);
	for (i = 0; i < rows; i++) {
		tmp = value[i];
		for (j = 0; j < columns; j++)
			tmp[j].assign(LiDIA::randomize(S));
	}
}



//********************************************
//** diag function ***************************
//********************************************

void bigmod_matrix::
diag(const bigint & a, const bigint & b)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "diag(const bigint &, const bigint &)");
	register lidia_size_t i, j;
	bigint *tmp;
	bigint a1, b1;
	if (!p) {
		a1.assign(a);
		b1.assign(b);
	}
	else {
		LiDIA::remainder(a1, a, p);
		if (a1.is_negative())
			LiDIA::add(a1, a1, p);
		LiDIA::remainder(b1, b, p);
		if (b1.is_negative())
			LiDIA::add(b1, b1, p);
	}
	for (i = 0; i < rows; i++) {
		tmp = value[i];
		for (j = 0; j < columns; j++)
			if (i == j)
				tmp[j].assign(a1);
			else
				tmp[j].assign(b1);
	}
}



void
diag(bigmod_matrix &A, const bigint & a, const bigint & b)
{
	debug_handler("bigmod_matrix", "in function "
		      "diag(bigmod_matrix &, const bigint &, const bigint &)");
	A.diag(a, b);
}



//
// transpose function
//

bigmod_matrix bigmod_matrix::
trans() const
{
	debug_handler("bigmod_matrix", "in member - function "
		      "trans()");
	register lidia_size_t i, j;
	bigint *tmp;
	bigmod_matrix AT(columns, rows, p);
	for (i = 0; i < rows; i++) {
		tmp = value[i];
		for (j = 0; j < columns; j++)
			AT.value[j][i].assign(tmp[j]);
	}
	return AT;
}



bigmod_matrix
trans(const bigmod_matrix &A)
{
	debug_handler("bigmod_matrix", "in function "
		      "trans(const bigmod_matrix &)");
	return A.trans();
}



//********************************************
//** regular expansion ***********************
//********************************************

void bigmod_matrix::
regexpansion(const lidia_size_t *v)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "regexpansion(const lidia_size_t *)");
	register lidia_size_t k = v[0];
	if (columns > k) {
		register lidia_size_t i = 0, j = 0;
		lidia_size_t diff = columns-rows;
		bigmod_matrix A(*this);
		set_no_of_rows(columns);
		set_no_of_columns(columns);
		bigmod_matrix ERG(diff, columns, p);
		while (i< columns && k > 0 && j < diff) {
			if (i != v[k]) {
				ERG.value[j][i].assign_one();
				j++;
			}
			else
				k--;
			i++;
		}
		compose_v(ERG, A);
	}
}



void
regexpansion(bigmod_matrix &A, const lidia_size_t *v)
{
	debug_handler("bigmod_matrix", "in function "
		      "regexpansion(const bigmod_matrix &, lidia_size_t *)");
	A.regexpansion(v);
}



//
// rank and linearly independent rows
//

lidia_size_t *bigmod_matrix::
lininr(bigint & factor) const
{
	debug_handler("bigmod_matrix", "in member - function "
		      "lininr()");
	register lidia_size_t i, j;
	lidia_size_t *l = new lidia_size_t[columns+1];

	bigmod_matrix B(*this);
	// Step 1,2
	LiDIA::stf(B, factor);
	if (!factor.is_one()) {
//    std::cout << "While computing stf mod "<<p << " found factor "<< d<<std::endl;
		return l;
	}

	// Step 3 - 12
	for (j = 0; j < columns; j++) {
		i = rows - 1;
		while (i >= 0 && value[i][j].is_zero())
			i--;
		l[j] = i;
	}

	// Step 13 - 24
	lidia_size_t j0 = 0;
	while (j0 < columns && l[j0] == -1)
		j0++;
	lidia_size_t r = columns - j0 + 1;
	lidia_size_t *IND = new lidia_size_t[r+1];
	for (j = 0; j < r - 1; j++)
		IND[r - j-1] = l[j + j0];
	IND[0] = r - 1;
	delete[] l;
	return IND;
}



lidia_size_t *
lininr(const bigmod_matrix & B, bigint & factor)
{
	debug_handler("bigmod_matrix", "in function "
		      "lininr(const bigmod_matrix &)");
	return B.lininr(factor);
}



//
// inverse and adjoint matrix
//

void bigmod_matrix::
inv(const bigmod_matrix& A, bigint & factor)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "inv(const bigmod_matrix &)");

	if (A.columns != A.rows) {
		lidia_error_handler("bigmod_matrix", "adj :: non square matrix");
		return;
	}
	register lidia_size_t i, j, z;

	// Step 1,2
	if (columns != A.rows)
		set_no_of_columns(A.rows);
	if (rows != A.rows)
		set_no_of_rows(A.rows);
	p.assign(A.p);
	bigmod_matrix C(A);

	bigint *Btmp, tmp, tmp1, q;
	for (i = 0; i < C.rows; i++) {
		Btmp = value[i];
		for (j = 0; j < C.rows; j++)
			if (i == j)
				Btmp[j].assign_one();
			else
				Btmp[j].assign_zero();
	}

	// Step 3 - 5
	for (i = C.rows - 1; i >= 0; i--) {

		// Step 6 - 11
		j = i;
		while (j >= 0 && C.value[i][j].is_zero())
			j--;

		// Step 12 - 19
		if (j != i) {
			C.swap_columns(i, j);
			swap_columns(i, j);
		}

		// Step 20 - 32
		for (j = 0; j < C.rows; j++) {
			if (j != i) {
				C.div_mod(q, factor, C.value[i][j], C.value[i][i]);
				if (!factor.is_one()) {
					return;
				}
				for (z = 0; z < C.rows; z++) {
					C.mult_mod(tmp, C.value[z][i], q);
					C.sub_mod(C.value[z][j], C.value[z][j], tmp);
					mult_mod(tmp, value[z][i], q);
					sub_mod(value[z][j], value[z][j], tmp);
				}
			}
		}
	}

	// Step 33 - 43
	for (j = 0; j < C.rows; j++) {
		C.div_mod(q, factor, 1, C.value[j][j]);
		if (!factor.is_one()) {
			return;
		}
		for (z = 0; z < C.rows; z++) {
			C.mult_mod(C.value[z][j], C.value[z][j], q);
			mult_mod(value[z][j], value[z][j], q);
		}
	}
}



bigmod_matrix
inv(const bigmod_matrix& A, bigint & factor)
{
	debug_handler("bigmod_matrix", "in function "
		      "inv(const bigmod_matrix &)");
	bigmod_matrix B(A.rows, A.columns, A.p);
	B.inv(A, factor);
	return B;
}



void bigmod_matrix::
adj(const bigmod_matrix &A, bigint & factor)
{
	bigint DET;
	A.det(DET, factor);
	(*this).inv(A, factor);
	if (!factor.is_one()) {
//    std::cout << "While computing INV mod "<<p << " found factor "<< d<<std::endl;
		return;
	}

	(*this) *= DET;
}



bigmod_matrix
adj(const bigmod_matrix &A, bigint & factor)
{
	debug_handler("bigmod_matrix", "in function "
		      "adj(bigmod_matrix &)");
	bigmod_matrix B(A.rows, A.columns, A.p);
	B.adj(A, factor);
	return B;
}



//
// determinant
//

void bigmod_matrix::
det(bigint& ret, bigint & factor) const
{
	debug_handler("bigmod_matrix", "in member - function "
		      "det(bigmod &)");
	bigmod_matrix A(*this);
	lidia_size_t i, j, z, n = rows;
	bigint q, q1, *tmp, tmp1, *tmp2;

	// Step 1 - 4
	int ex = 1;
	ret.assign_one();

	// Step 5 - 8
	for (i = 0; i < n; i++) {

		// Step 9 - 13
		for (j = i; j < n && A.value[j][i].is_zero(); j++);

		// Step 14 - 26
		if (j == n) {
			ret.assign_zero();
			return;
		}
		if (j != i) {
			ex = -ex;
			tmp2 = A.value[j];
			A.value[j] = A.value[i];
			A.value[i] = tmp2;
		}
		tmp = A.value[i];

		// Step 27 - 29
		div_mod(q1, factor, 1, tmp[i]);
		if (!factor.is_one()) {
//    std::cout << "While computing det mod "<<p << " found factor "<< d<<std::endl;
			ret.assign_zero();
			return;
		}
		for (j = i+1; j < n; j++) {

			// Step 30 - 40
			tmp2 = A.value[j];
			mult_mod(q, tmp2[i], q1);
			for (z = i+1; z < n; z++) {
				mult_mod(tmp1, tmp[z], q);
				sub_mod(tmp2[z], tmp2[z], tmp1);
			}
		}
		mult_mod(ret, ret, tmp[i]);
	}
	if (ex < 0)
		ret = p-ret;
}



bigint bigmod_matrix::
det(bigint & factor) const
{
	debug_handler("bigmod_matrix", "in member - function det()");
	bigint DET;
	det(DET, factor);
	return DET;
}



bigint
det(const bigmod_matrix & A, bigint & factor)
{
	debug_handler("bigmod_matrix", "in function "
		      "det(bigmod_matrix &)");
	bigint DET;
	A.det(DET, factor);
	return DET;
}



//********************************************
//** characteristic polynom ******************
//********************************************

bigint * bigmod_matrix::
charpoly(bigint & factor) const
{
	debug_handler("bigmod_matrix", "in member - function "
		      "charpoly(bigint &)");
	register lidia_size_t i, j, r;
	bigint tmp;
	bigint sign;
	bigmod_matrix B(*this);

	// Step 1 - 5
	hbf(B, factor);
	if (!factor.is_one()) {
//    std::cout << "While computing hbf mod "<<p << " found factor "<< d<<std::endl;
		return new bigint[B.columns+1];
	}
	lidia_size_t n = B.columns;
	bigint *K = new bigint[n]; //size = n
	for (i = 0; i < n; i++)
		K[i].assign_one();

	// Step 6 - 8
	bigmod_matrix P(n + 1, n + 1, B.p); // default initial with zero
	P.value[0][0].assign_one();

	// Step 9 - 11
	for (r = 1; r <= n; r++) {

		// Step 12 - 16
		for (j = 1; j <= r-1; j++)
			B.mult_mod(K[j-1], K[j-1], B.value[r-1][r - 2]);

		// Step 17 - 23
		LiDIA::subtract(P.value[r][r], B.p, P.value[r-1][r-1]);
		for (i = 1; i <= r-1; i++) {
			B.mult_mod(tmp, B.value[r-1][r-1], P.value[i][r-1]);
			B.sub_mod(P.value[i][r], tmp, P.value[i - 1][r - 1]);
		}
		B.mult_mod(P.value[0][r], B.value[r-1][r-1], P.value[0][r - 1]);

		// Step 24 - 34
		sign.assign_one();
		for (j = r - 1; j >= 1; j--) {
			sign = -sign;
			for (i = 0; i <= j - 1; i++) {
				B.mult_mod(tmp, sign, P.value[i][j-1]);
				B.mult_mod(tmp, tmp, B.value[j-1][r-1]);
				B.mult_mod(tmp, tmp, K[j-1]);
				B.add_mod(P.value[i][r], P.value[i][r], tmp);
			}
		}
	}

	// Step 35 - 40
	bigint *RES = new bigint[n+1];
	for (i = 0; i < n+1; i++)
		RES[i].assign(P.value[i][n]);
	delete[] K;
	return RES;
}



bigint *
charpoly(const bigmod_matrix & A, bigint & factor)
{
	debug_handler("bigmod_matrix", "in function "
		      "charpoly(bigmod_matrix &, bigint &)");
	return A.charpoly(factor);
}



//********************************************
//** special matrixforms *********************
//********************************************

int
stf(bigmod_matrix &A, bigmod_matrix &TRANS, bigint & factor)
{
	debug_handler("bigmod_matrix", "in function "
		      "stf(bigmod_matrix &, bigmod_matrix &)");
	lidia_size_t n = A.columns;
	lidia_size_t m = A.rows;
	if (TRANS.columns != n || TRANS.rows != n)
		TRANS.resize(n, n);

	TRANS.diag(1, 0);
	factor.assign_one(); // Just in case no elimination step at all is needed!
	lidia_size_t i, j, z;
	bigint q = 1, q1;

	// Step 1 -4
	int exchange = 1;
	lidia_size_t j0 = n - 1;
	bigint *Atmp, tmp, tmp1;

	// Step 5 - 8
	for (i = m - 1; i >= 0; i--) {
		Atmp = A.value[i];
		// Step 9 - 13
		for (j = j0; j >= 0 && Atmp[j].is_zero(); j--);

		// Step 14 - 26
		if (j != -1) {
			if (j != j0) {
				exchange = -exchange;
				A.swap_columns(j0, j);
				TRANS.swap_columns(j0, j);
			}

			// Step 27 - 29
			for (j = j0 - 1; j >= 0; j--) {
				// Step 30 - 40
				if (!Atmp[j].is_zero()) {
					A.div_mod(q, factor, Atmp[j], Atmp[j0]);
					if (!factor.is_one()) {
//    std::cout << "While computing stf mod "<<p << " found factor "<< d<<std::endl;
						return 0;
					}
					for (z = 0; z <= i; z++) {
						A.mult_mod(tmp, A.value[z][j0], q);
						A.sub_mod(A.value[z][j], A.value[z][j], tmp);
					}
					for (z = 0; z < n; z++) {
						A.mult_mod(tmp, TRANS.value[z][j0], q);
						A.sub_mod(TRANS.value[z][j], TRANS.value[z][j], tmp);
					}
				}
			}
			// Step 41 - 48
			j0--;
		}
	}
	return exchange;
}



void bigmod_matrix::
stf(int exchange, bigint & factor)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "stf(int)");
	factor.assign_one(); // Just in case no elimination step at all is needed!
	lidia_size_t i, j, z, j0 = columns-1;
	bigint q = 1;
	bigint *tmp, tmp1, *tmp2;

	// Step 1 - 4
	exchange = 1;

	// Step 5 - 8
	for (i = rows - 1; i >= 0; i--) {
		tmp = value[i];

		// Step 9 - 13
		for (j = j0; j >= 0 && tmp[j].is_zero(); j--);

		// Step 14 - 26
		if (j == -1)
			return;
		if (j != j0) {
			exchange = -exchange;
			for (z = 0; z < rows; z++) {
				tmp2 = value[z];
				tmp1.assign(tmp2[j0]);
				tmp2[j0].assign(tmp2[j]);
				tmp2[j].assign(tmp1);
			}
		}

		// Step 27 - 29
		for (j = j0 - 1; j >= 0; j--) {

			// Step 30 - 40
			div_mod(q, factor, tmp[j], tmp[j0]); //q = tmp[j] / tmp[j0];
			if (!factor.is_one()) {
//    std::cout << "While computing stf mod "<<p << " found factor "<< d<<std::endl;
				return;
			}
			for (z = 0; z <= i; z++) {
				tmp2 = value[z];
				mult_mod(tmp1, tmp2[j0], q);
				sub_mod(tmp2[j], tmp2[j], tmp1);
			}
		}

		// Step 41 - 48
		j0--;
	}
}



bigint
stf(bigmod_matrix &A, bigint & factor)
{
	lidia_size_t n = A.columns;
	lidia_size_t m = A.rows;
	debug_handler("bigmod_matrix", "in function "
		      "stf(bigmod_matrix &)");
	factor.assign_one(); // Just in case no elimination step at all is needed!
	lidia_size_t i, j, z;
	bigint q = 1, tmp;

	// Step 1 - 4
	bigint exchange = 1;
	lidia_size_t j0 = n - 1;

	// Step 5 - 8
	for (i = m - 1; i >= 0; i--) {
		j = j0;

		// Step 9 - 13
		while (j >= 0 && A.value[i][j].is_zero())
			j--;

		// Step 14 - 26
		if (j != -1) {
			if (j != j0) {
				exchange = -exchange;
				A.swap_columns(j0, j);
			}

			// Step 27 - 29
			for (j = j0 - 1; j >= 0; j--) {

				// Step 30 - 40
				A.div_mod(q, factor, A.value[i][j], A.value[i][j0]);
				if (!factor.is_one()) {
//    std::cout << "While computing stf mod "<<p << " found factor "<< d<<std::endl;
					return 0;
				}
				for (z = 0; z < m; z++) {
					A.mult_mod(tmp, A.value[z][j0], q);
					A.sub_mod(A.value[z][j], A.value[z][j], tmp);
				}
			}

			// Step 41 - 48
			j0--;
		}
	}
	return exchange;
}



void
hbf(bigmod_matrix &A, bigint & factor)
{
	debug_handler("bigmod_matrix", "in function "
		      "hbf(bigmod_matrix &)");

	// Step 1,2
	lidia_size_t i, j, z;
	factor.assign_one(); // Just in case no elimination step at all is needed!
	bigint q = 1, tmp;
	lidia_size_t n = A.columns;
	lidia_size_t m = A.rows;

	// Step 3 - 11
	for (i = n - 1; i >= 1; i--) {
		j = i - 1;
		while (j >= 0 && A.value[i][j].is_zero())
			j--;
		if (j != -1) {

			// Step 12,13
			if (j != i - 1) {

				// Step 14 - 18
				A.swap_columns(i - 1, j);

				// Step 19 - 24
				A.swap_rows(i - 1, j);
			}

			// Step 25 - 41
			for (j = i - 2; j >= 0; j--) {
				A.div_mod(q, factor, A.value[i][j], A.value[i][i - 1]);
				if (!factor.is_one()) {
//    std::cout << "While computing hbf mod "<<p << " found factor "<< d<<std::endl;
					return;
				}

				for (z = 0; z < m; z++) {
					A.mult_mod(tmp, A.value[z][i-1], q);
					A.sub_mod(A.value[z][j], A.value[z][j], tmp);
				}
				for (z = 0; z < n; z++) {
					A.mult_mod(tmp, A.value[j][z], q);
					A.add_mod(A.value[i-1][z], A.value[i-1][z], tmp);
				}
			}
		}
	}
}



//
// Kernel
//

void bigmod_matrix::
kernel(const bigmod_matrix & A, bigint & factor)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "kernel(const bigmod_matrix &, bigint &)");
	if (A.p.is_zero()) {
		bigint_matrix tmp = (*this);
		tmp.kernel(A);
		assign(tmp);
		factor.assign_one();
		return;
	}

	bigmod_matrix TRANS(A.columns, A.columns, A.p);
	bigmod_matrix COPY(A);

	// Step 1
	LiDIA::stf(COPY, TRANS, factor);
	if (!factor.is_one()) {
//    std::cout << "While computing kernel mod "<<p << " found factor "<< d<<std::endl;
		return;
	}


	// Step 2
	register lidia_size_t a = 0;
	while (a < COPY.columns && COPY.is_column_zero(a))
		a++;
	if (a == 0) {
		bigmod_matrix KERN(COPY.columns, 1, COPY.p);
		assign(KERN);
		return;
	}

	// Step 3
	set_no_of_rows(COPY.columns);
	set_no_of_columns(a);
	p = COPY.p;
	bigint *KERNtmp, *TRANStmp;
	register lidia_size_t i, z;
	for (i = 0; i < COPY.columns; i++) {
		KERNtmp = value[i];
		TRANStmp = TRANS.value[i];
		for (z = 0; z < a; z++)
			KERNtmp[z].assign(TRANStmp[z]);
	}
}



bigmod_matrix
kernel(const bigmod_matrix & A, bigint & factor)
{
	debug_handler("bigmod_matrix",
		      "in function kernel(const bigmod_matrix &, bigint &)");
	bigmod_matrix RET;
	RET.kernel(A, factor);
	return RET;
}



//
// Image
//

void bigmod_matrix::
image(const bigmod_matrix &A, bigint & factor)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "image(const bigmod_matrix &, bigint &)");
	if (A.p.is_zero()) {
		factor.assign_one();
		if (A.is_matrix_zero()) {
			if (rows != A.rows)
				set_no_of_rows(A.rows);
			if (!p.is_zero())
				p.assign_zero();
			if (columns != 1)
				set_no_of_columns(1);
			for (lidia_size_t i = 0; i < rows; i++)
				value[i][0].assign_zero();
			return;
		}
		bigint_matrix tmp = (*this);
		tmp.image(A);
		assign(tmp);
		return;
	}
	bigmod_matrix TRANS(A.columns, A.columns, A.p);

	assign(A);
	// Step 1
	LiDIA::stf(*this, TRANS, factor);
	if (!factor.is_one()) {
//    std::cout << "While computing image mod "<<p << " found factor "<< d<<std::endl;
		return;
	}

	register lidia_size_t a = 0;
	// Step 2
	while (a < columns && is_column_zero(a))
		a++;

	// Step 3
	bigint *BILDtmp, *tmp;
	register lidia_size_t i, z;
	for (i = 0; i < rows; i++) {
		BILDtmp = value[i];
		tmp = value[i];
		for (z = 0; z < columns - a; z++)
			BILDtmp[z].assign(tmp[a+z]);
	}
	set_no_of_columns(A.columns - a);
}



bigmod_matrix
image(const bigmod_matrix & A, bigint & factor)
{
	debug_handler("bigmod_matrix", "in function "
		      "image(const bigmod_matrix &, bigint &)");
	bigmod_matrix BILD;
	BILD.image(A, factor);
	return BILD;
}



//
// solve
//

void bigmod_matrix::
solve(const bigmod_matrix &B, const bigint *b, bigint & factor)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "solve(const bigmod_matrix &, const bigint *, bigint &)");
	if (B.p.is_zero()) {
		bigint_matrix tmp = (*this);
		tmp.solve(B, b);
		assign(tmp);
		factor.assign_one();
		return;
	}
	register lidia_size_t i, k, l1 = 0;
	bigint tmprem;

	// Step 1
	bigmod_matrix A = B;
	A.set_no_of_columns(B.columns + 1);
	for (i = 0; i < B.rows; i++) {
		LiDIA::remainder(tmprem, b[i], A.p);
		if (tmprem.is_negative())
			LiDIA::add(tmprem, tmprem, A.p);
		A.value[i][B.columns].assign(tmprem);
	}
	bigmod_matrix TRANS(A.columns, A.columns, A.p);
	LiDIA::stf(A, TRANS, factor);
	if (!factor.is_one()) {
//    std::cout << "While solving mod "<<A.p << " found factor "<< d<<std::endl;
		return;
	}

	// Step 2,3
	while (A.is_column_zero(l1) && l1 < A.columns)
		l1++;
	if (l1 == 0) {
		bigmod_matrix C(B.columns, 1, A.p);
		assign(C);
		return;
	}

	// Step 4
	i = 0;
	while (TRANS.value[B.columns][i].is_zero() && i < l1)
		i++;
	if (i == l1) {
		bigmod_matrix C(B.columns, 1, A.p);
		assign(C);
		return;
	}

	// Step 5
	bigint *x = new bigint[B.columns];
	bigint divisor = A.p - TRANS.value[B.columns][i];
	for (k = 0; k < B.columns; k++) {
		A.div_mod(x[k], factor, TRANS.value[k][i], divisor);
		if (!factor.is_one()) {
//    std::cout << "While solving mod "<<A.p << " found factor "<< d<<std::endl;
			delete[] x;
			return;
		}
	}
	kernel(B, factor);
	if (!factor.is_one()) {
//    std::cout << "While solving mod "<<A.p << " found factor "<< d<<std::endl;
		delete[] x;
		return;
	}
	set_no_of_columns(columns + 1);
	for (i = 0; i < rows; i++)
		value[i][columns-1].assign(x[i]);
	delete[] x;
}



bigmod_matrix
solve(const bigmod_matrix & A, const bigint * b, bigint & factor)
{
	debug_handler("bigmod_matrix", "in function "
		      "solve(const bigmod_matrix &, const bigint *)");
	bigmod_matrix B;
	B.solve(A, b, factor);
	return B;
}



void bigmod_matrix::
solve(const bigmod_matrix &B, const base_vector< bigint > &b, bigint & factor)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "solve(const bigmod_matrix &, const bigint *, bigint &)");
	if (B.p.is_zero()) {
		bigint_matrix tmp = (*this);
		tmp.solve(B, b);
		assign(tmp);
		return;
	}
	register lidia_size_t i, k, l1 = 0;
	bigint tmprem;

	// Step 1
	bigmod_matrix A = B;
	A.set_no_of_columns(B.columns + 1);
	for (i = 0; i < B.rows; i++) {
		LiDIA::remainder(tmprem, b[i], A.p);
		if (tmprem.is_negative())
			LiDIA::add(tmprem, tmprem, A.p);
		A.value[i][B.columns].assign(tmprem);
	}
	bigmod_matrix TRANS(A.columns, A.columns, A.p);
	LiDIA::stf(A, TRANS, factor);
	if (!factor.is_one()) {
//    std::cout << "While solving mod "<<A.p << " found factor "<< d<<std::endl;
		return;
	}

	// Step 2,3
	while (A.is_column_zero(l1) && l1 < A.columns)
		l1++;
	if (l1 == 0) {
		bigmod_matrix C(B.columns, 1, A.p);
		assign(C);
		return;
	}

	// Step 4
	i = 0;
	while (TRANS.value[B.columns][i].is_zero() && i < l1)
		i++;
	if (i == l1) {
		bigmod_matrix C(B.columns, 1, A.p);
		assign(C);
		return;
	}

	// Step 5
	bigint *x = new bigint[B.columns];
	bigint divisor = A.p - TRANS.value[B.columns][i];
	for (k = 0; k < B.columns; k++) {
		A.div_mod(x[k], factor, TRANS.value[k][i], divisor);
		if (!factor.is_one()) {
//    std::cout << "While solving mod "<<A.p << " found factor "<< d<<std::endl;
			delete[] x;
			return;
		}
	}
	kernel(B, factor);
	if (!factor.is_one()) {
//    std::cout << "While solving mod "<<A.p << " found factor "<< d<<std::endl;
		delete[] x;
		return;
	}
	set_no_of_columns(columns + 1);
	for (i = 0; i < rows; i++)
		value[i][columns-1].assign(x[i]);
	delete[] x;
}



bigmod_matrix
solve(const bigmod_matrix & A, const base_vector< bigint > & b, bigint & factor)
{
	debug_handler("bigmod_matrix", "in function "
		      "solve(const bigmod_matrix &, const bigint *)");
	bigmod_matrix B;
	B.solve(A, b, factor);
	return B;
}



// LINEAR ALGEBRA - METHOD II:
// We are not satisfied with finding a factor.

//
// Kernel
//

void bigmod_matrix::kernel(const bigmod_matrix & A)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "kernel(const bigmod_matrix &)");
	if (A.p.is_zero()) {
		bigint_matrix tmp;
		tmp.kernel(A);
		assign(tmp);
		return;
	}
	register lidia_size_t z = 0, i;
	bigmod_matrix COPY(A);

	p.assign(A.p);
	if (columns != A.columns)
		set_no_of_columns(A.columns);
	if (rows != A.columns)
		set_no_of_rows(A.columns);
	diag(1, 0);

	register lidia_size_t j, act_row, index;
	register lidia_size_t res_j = A.columns - 1;
	register char SW;
	bigint factor, q, tmp1;
	bigint *tmp, *tmp2, *TEMPtmp;

	for (act_row = COPY.rows - 1; act_row >= 0; act_row--) {
		tmp = COPY.value[act_row];
#ifdef MNF_DEBUG
		std::cout << "working on row " << act_row << std::endl;
#endif
		for (j = res_j; j >= 0 && tmp[j].is_zero(); j--);
		if (j == -1) {
			continue; // next iteration for act_row
		}
		index = j; // "working column"
		// Attention : We don't swap columns yet!!!
		do {
			SW = 0;
			tmp1 = gcd(tmp[index], p);
			for (j = res_j; j >= 0 && !tmp1.is_one(); j--) {
				if (!tmp[j].is_zero()) {
					q = gcd(tmp[j], p);
					if (tmp1 > q) {
						index = j;
						tmp1 = q;
					}
				}
			}
#ifdef MNF_DEBUG
			std::cout << "Found working column " << index << std::endl;
#endif
			// eventuelly different working column.
			// STILL NO COLUMN SWAPPING!

			// Simplify current column
			inv_mod(q, factor, tmp[index]); //q = 1 / tmp[index];
			if (!(gcd(q, p)).is_one()) {
#ifdef MNF_DEBUG
				std::cout << "Inverting " << tmp[index] << " mod " << p;
				std::cout << "resulted in " << q << std::endl;
#endif
				lidia_error_handler("bigmod_matrix", "kernel::internal division error");
				return;
			}

			for (i = 0; i <= act_row; i++) {
				tmp2 = COPY.value[i];
				mult_mod(tmp2[index], tmp2[index], q);
			}
			for (i = 0; i < A.columns; i++) {
				TEMPtmp = value[i];
				mult_mod(TEMPtmp[index], TEMPtmp[index], q);
			}
#ifdef MNF_DEBUG
			std::cout << "Simplified working column" << std::endl;
			std::cout << COPY << (*this);
#endif

			if (tmp[index].is_one()) {
				// Step 27 - 29: Eliminate to the left
				for (j = res_j; j >= 0; j--)
					// Step 30 - 40
					if ((j != index) && !tmp[j].is_zero()) {
						q = tmp[j];
						for (i = 0; i <= act_row; i++) {
							tmp2 = COPY.value[i];
							mult_mod(tmp1, q, tmp2[index]);
							sub_mod(tmp2[j], tmp2[j], tmp1);
						}
						for (i = 0; i < columns; i++) {
							TEMPtmp = value[i];
							mult_mod(tmp1, q, TEMPtmp[index]);
							sub_mod(TEMPtmp[j], TEMPtmp[j], tmp1);
						}
					}
#ifdef MNF_DEBUG
				std::cout << "Eliminated with working column" << index << std::endl;
				std::cout << COPY << (*this);
#endif
			}
			else {
				// i.e. we have a zero divisor, so everything is difficult
				for (j = res_j; j >= 0; j--)
					if ((j != index) && !tmp[j].is_zero()) {
						SW = 1;
						div_rem(q, tmp1, tmp[j], tmp[index]);
						for (i = 0; i <= act_row; i++) {
							tmp2 = COPY.value[i];
							mult_mod(tmp1, q, tmp2[index]);
							sub_mod(tmp2[j], tmp2[j], tmp1);
						}
						for (i = 0; i < columns; i++) {
							TEMPtmp = value[i];
							mult_mod(tmp1, q, TEMPtmp[index]);
							sub_mod(TEMPtmp[j], TEMPtmp[j], tmp1);
						}
					}
#ifdef MNF_DEBUG
				std::cout << "Eliminated with working column" << std::endl;
				std::cout << COPY << (*this);
#endif
			}
		}while (SW == 1);

#ifdef MNF_DEBUG
		std::cout << "Before swapping" << std::endl;
		std::cout << COPY << (*this);
#endif
		if (index != res_j) {
			COPY.swap_columns(res_j, index);
			swap_columns(res_j, index);
		}
#ifdef MNF_DEBUG
		std::cout << "Swapped columns (if necessary)" << index << " and "
			  << res_j << std::endl;
		std::cout << COPY << (*this);
#endif
		if (tmp[res_j].is_one()) {
			res_j--;
		}
		else {
			// Now get a zero at COPY[act_row][res_j], i.e.
			// multiply the column by
			// p/COPY.value[act_row-1][res_j])
			LiDIA::remainder(tmp1, p, tmp[res_j]);
			if (!(tmp1.is_zero())) {
				lidia_error_handler("bigmod_matrix", "mnf::internal error");
				return;
			}
			LiDIA::divide(tmp1, p, tmp[res_j]);
#ifdef MNF_DEBUG
			std::cout << "\nCalling div_mod with" << tmp1 << std::endl << std::flush;
#endif
			for (i = 0; i <= act_row; i++)
				mult_mod(COPY.value[i][res_j], COPY.value[i][res_j],
					 tmp1);
			for (i = 0; i < columns; i++)
				mult_mod(value[i][res_j], value[i][res_j], tmp1);
#ifdef MNF_DEBUG
			std::cout << COPY << (*this);
#endif
		}
	}

	// Step 4: Store result
	z = 0;
	for (j = 0; j+z <= res_j; j++) {
		while (j + z <= res_j && is_column_zero(j + z)) z++;
		if (z > 0 && j + z <= res_j)
			for (i = 0; i < rows; i++)
				value[i][j].assign(value[i][j+z]);
	}

	res_j -= z-1;
	if (res_j)
		set_no_of_columns(res_j);
	else {
		set_no_of_columns(1);
		for (i = 0; i < rows; i++) value[i][0].assign_zero();
	}
}



bigmod_matrix kernel(const bigmod_matrix & A)
{
	bigmod_matrix RET;
	RET.kernel(A);
	return RET;
}



//
// Image
//

void bigmod_matrix::image(const bigmod_matrix & A)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "image(const bigmod_matrix &)");
	if (A.p.is_zero()) {
		if (A.is_matrix_zero()) {
			if (rows != A.rows)
				set_no_of_rows(A.rows);
			if (!p.is_zero())
				p.assign_zero();
			if (columns != 1)
				set_no_of_columns(1);
			for (lidia_size_t i = 0; i < rows; i++)
				value[i][0].assign_zero();
			return;
		}
		bigint_matrix tmp;
		tmp.image(A);
		assign(tmp);
		return;
	}

	if (this != &A) assign(A);

//  std::cout <<"Computing image of "<<(*this)<<std::endl;
	// Step 1
	bigmod_matrix RES(rows, columns, p);
	register lidia_size_t i, j, act_row, act_column = columns-1;
	register lidia_size_t res_j = act_column;
	register lidia_size_t index;
	register char SW, SAVE = 1;
	bigint factor, q, tmp1;
	bigint *tmp, *tmp2;

	for (act_row = rows - 1; act_row >= 0; act_row--) {
//     std::cout <<"computing on row " << act_row<<std::endl;
		tmp = value[act_row];
		for (j = act_column; j >= 0 && tmp[j].is_zero(); j--);
		if (j == -1) {
			continue; // next iteration for act_row
		}
		index = j; // "working column"
		// Attention : We don't swap columns yet!!!
		do {
			SW = 0;
			tmp1 = gcd(tmp[index], p);
			for (j = act_column; j >= 0 && !tmp1.is_one(); j--) {
				if (!tmp[j].is_zero()) {
					q = gcd(tmp[j], p);
					if (tmp1 > q) {
						index = j;
						tmp1 = q;
					}
				}
			}
			if ((SAVE == 0) && (index != act_column)) {
				SAVE = 1; // This column will be needed
				act_column--; //  in the generating system,
				// the actual column will be useless!!
			}

//       std::cout <<" Found working column "<<index<< " SAVE ="<<(int)SAVE<<std::endl;
			// eventuelly different working column.
			// STILL NO COLUMN SWAPPING!

			// Simplify current column
			inv_mod(q, factor, tmp[index]); //q = 1 / tmp[index];
			if (!(gcd(q, p)).is_one()) {
				std::cout << "Inverting " << tmp[index] << " mod " << p;
				std::cout << "resulted in " << q << std::endl;
				lidia_error_handler("bigmod_matrix", "mnf::internal division error");
				return;
			}

			for (i = 0; i <= act_row; i++) {
				tmp2 = value[i];
				mult_mod(tmp2[index], tmp2[index], q);
			}

			if (tmp[index].is_one()) {
				// Step 27 - 29: Eliminate to the left
				for (j = act_column; j >= 0; j--)
					// Step 30 - 40
					if ((j != index) && !tmp[j].is_zero()) {
// 	    std::cout <<" Eliminating column "<< j <<" of" <<(*this)<<" with column"<<index<<std::endl;
						q = tmp[j];
						for (i = 0; i <= act_row; i++) {
							tmp2 = value[i];
							mult_mod(tmp1, q, tmp2[index]);
							sub_mod(tmp2[j], tmp2[j], tmp1);
						}
//  	    std::cout <<" Result is "<<(*this)<<std::endl;
					}
			}
			else {
				// i.e. we have a zero divisor, so everything is difficult
				for (j = act_column; j >= 0; j--)
					if ((j != index) && !tmp[j].is_zero()) {
						SW = 1;
//  	    std::cout <<" Eliminating column "<< j <<" of" <<(*this)<<" with column"<<index<<std::endl;
						div_rem(q, tmp1, tmp[j], tmp[index]);
						for (i = 0; i <= act_row; i++) {
							tmp2 = value[i];
							mult_mod(tmp1, q, tmp2[index]);
							sub_mod(tmp2[j], tmp2[j], tmp1);
						}
//  	    std::cout <<" Result is "<<(*this)<<std::endl;
					}
			}
		} while (SW == 1);

		if (index != act_column)
			swap_columns(act_column, index);

		// Reducing with a column only is possible, if we don't
		// throw it away afterwards. However, we may reduce columns
		// res_j+2, ... by a multiple of column res_j+1 !!!

		for (j = res_j+1+(!SAVE); j < columns; j++) {
//      std::cout <<" Reducing column "<< j <<" of" <<(RES)<<" with column"<<act_column <<" of"<<(*this)<<std::endl;
			if (RES.value[act_row][j] >= tmp[act_column]) {
				LiDIA::div_rem(q, RES.value[act_row][j],
						  RES.value[act_row][j], tmp[act_column]);
				for (i = 0; i < act_row; i++) {
					mult_mod(tmp1, q, value[i][act_column]);
					sub_mod(RES.value[i][j], RES.value[i][j], tmp1);
				}
			}
		}
//     std::cout <<" Result is "<<(RES)<<std::endl;
		// Store this column in the result:
		if (SAVE) {
//       std::cout <<"Store column"<<act_column<<std::endl;
			for (i = 0; i <= act_row; i++)
				RES.value[i][res_j].assign(value[i][act_column]);
			res_j--;
//    	   std::cout <<"RES is"<<RES<<std::endl;
		}
		if (tmp[act_column].is_one()) {
			// Work on next column:
			act_column--;
		}
		else if (act_row != 0) {
			// Now we try, to eliminate one of the previous
			// columns, therefore multiply the column by
			// p/tmp[act_column]/value[act_row-1][act_column])
			LiDIA::remainder(tmp1, p, tmp[act_column]);
			if (!(tmp1.is_zero())) {
				lidia_error_handler("bigmod_matrix", "mnf::internal error");
				return;
			}
			LiDIA::divide(tmp1, p, tmp[act_column]);
			div_mod(tmp1, factor, tmp1, value[act_row-1][act_column]);
			for (i = 0; i <= act_row; i++)
				mult_mod(value[i][act_column], value[i][act_column], tmp1);
			if (is_column_zero(act_column)) {
				SAVE = 1;
				act_column--;
			}
			else SAVE = 0;
		}
	}
//   std::cout <<"temporary image is "<<RES<<std::endl;
	// Count zero columns.
	register lidia_size_t a = 0;
	while (a < RES.columns && RES.is_column_zero(a))
		a++;

	// Remove zero columns.
	if (a == A.columns)
		set_no_of_columns(1);
	else
		set_no_of_columns(A.columns - a);
	RES.split_h(RES, *this);
}



bigmod_matrix image(const bigmod_matrix & A)
{
	bigmod_matrix BILD;
	BILD.image(A);
	return BILD;
}



void bigmod_matrix::unique_image(const bigmod_matrix & A)
{
	debug_handler("bigmod_matrix", "in member - function "
		      "image(const bigmod_matrix &)");
	if (A.p.is_zero()) {
		bigint_matrix RES(A);
		RES.hnf_havas_cont();
		lidia_size_t a = 0;
		while (a < RES.get_no_of_columns() && RES.is_column_zero(a))
			a++;

		// Step 3
		if (rows != A.rows)
			set_no_of_rows(A.rows);
		if (a == RES.get_no_of_columns())
			set_no_of_columns(1);
		else
			set_no_of_columns(RES.get_no_of_columns() - a);
		RES.split_h(RES, *this);
		return;
	}

	if (this != &A) assign(A);

	// Step 1
	bigmod_matrix RES(rows, rows, p);
	register lidia_size_t i, j, act_row, act_column = columns - 1;
	register lidia_size_t res_j = rows;
	register lidia_size_t index;
	register char SW;
	bigint factor, q, tmp1;
	bigint *tmp, *tmp2;

	for (act_row = rows - 1; act_row >= 0; act_row--) {
		tmp = value[act_row];
		for (j = act_column; j >= 0 && tmp[j].is_zero(); j--);
		if (j == -1) {
			continue; // next iteration for act_row
		}
		index = j; // "working column"
		// Attention : We don't swap columns yet!!!
		do {
			SW = 0;
			tmp1 = gcd(tmp[index], p);
			for (j = act_column; j >= 0 && !tmp1.is_one(); j--) {
				if (!tmp[j].is_zero()) {
					q = gcd(tmp[j], p);
					if (tmp1 > q) {
						index = j;
						tmp1 = q;
					}
				}
			}

			// eventuelly different working column.
			// STILL NO COLUMN SWAPPING!

			// Simplify current column
			inv_mod(q, factor, tmp[index]); //q = 1 / tmp[index];
			if (!(gcd(q, p)).is_one()) {
				std::cout << "Inverting " << tmp[index] << " mod " << p;
				std::cout << "resulted in " << q << std::endl;
				lidia_error_handler("bigmod_matrix", "mnf::internal division error");
				return;
			}

			for (i = 0; i <= act_row; i++) {
				tmp2 = value[i];
				mult_mod(tmp2[index], tmp2[index], q);
			}

			if (tmp[index].is_one()) {
				// Step 27 - 29: Eliminate to the left
				for (j = act_column; j >= 0; j--)
					// Step 30 - 40
					if ((j != index) && !tmp[j].is_zero()) {
						q = tmp[j];
						for (i = 0; i <= act_row; i++) {
							tmp2 = value[i];
							mult_mod(tmp1, q, tmp2[index]);
							sub_mod(tmp2[j], tmp2[j], tmp1);
						}
					}
			}
			else {
				// i.e. we have a zero divisor, so everything is difficult
				for (j = act_column; j >= 0; j--)
					if ((j != index) && !tmp[j].is_zero()) {
						SW = 1;
						div_rem(q, tmp1, tmp[j], tmp[index]);
						for (i = 0; i <= act_row; i++) {
							tmp2 = value[i];
							mult_mod(tmp1, q, tmp2[index]);
							sub_mod(tmp2[j], tmp2[j], tmp1);
						}
					}
			}
		} while (SW == 1);

		if (index != act_column)
			swap_columns(act_column, index);


		// Reduce to the right
		for (j = res_j; j < RES.columns; j++) {
			if (RES.value[act_row][j] >= tmp[act_column]) {
				LiDIA::div_rem(q, RES.value[act_row][j],
						  RES.value[act_row][j], tmp[act_column]);
				for (i = 0; i < act_row; i++) {
					mult_mod(tmp1, q, value[i][act_column]);
					sub_mod(RES.value[i][j], RES.value[i][j], tmp1);
				}
			}
		}
		// Store this column in the result:
		res_j--;
		for (i = 0; i <= act_row; i++)
			RES.value[i][res_j].assign(value[i][act_column]);

		if (tmp[act_column].is_one()) {
			// Work on next column:
			act_column--;
		}
		else if (act_row != 0) {
			// Now we try, to eliminate one of the previous
			// columns, therefore multiply the column by
			// p/tmp[act_column]/value[act_row-1][act_column])
			LiDIA::remainder(tmp1, p, tmp[act_column]);
			if (!(tmp1.is_zero())) {
				lidia_error_handler("bigmod_matrix", "mnf::internal error");
				return;
			}
			LiDIA::divide(tmp1, p, tmp[act_column]);
			div_mod(tmp1, factor, tmp1, value[act_row-1][act_column]);
			for (i = 0; i <= act_row; i++)
				mult_mod(value[i][act_column], value[i][act_column], tmp1);
			if (is_column_zero(act_column))
				act_column--;
		}
	}

	// Step 2
	if (res_j == RES.columns)
		set_no_of_columns(1);
	else
		set_no_of_columns(RES.columns - res_j);
	RES.split_h(RES, *this);
}



bigmod_matrix unique_image(const bigmod_matrix & A)
{
	bigmod_matrix BILD;
	BILD.unique_image(A);
	return BILD;
}



//
// exponent
//

//
//  bigint_matrix res(n+1, n), Id(n, n);
//
//  Id.diag(1, 0);
//  res.reginvimage(a.base, Id);
//
//  bigint result = res(n, 0);
//  for (lidia_size_t i = 1; i < n; i++)
//      result = lcm(result, res(n, i));
//  return bigrational(result, a.den);



void bigmod_matrix::
exponent(bigint& ret) const
{
	debug_handler("bigmod_matrix", "in member - function "
		      "exponent(bigint &)");
	bigmod_matrix A(*this);
	A.lift(0);
	//  bigmod_matrix A;
	//  A.unique_image(*this);
	if (A.columns != A.rows) {
		ret.assign(p);
		return;
	}
	// So we know that every diagonal element is set!!
//	if (A.value[0][0] == p) {
//		ret.assign(p);
//		return;
//	}
	// Maybe exponent is a true divisor of p;
	bigint multiplier, tmp;
	bigint * tmp1;
	ret.assign_one();
	if (!A.p.is_zero()) {
#if 0
		for (lidia_size_t i = 0; i < A.rows; i++){
			tmp1 = A.value[i];
			for (lidia_size_t j = i + 1; (j<A.columns) && (!tmp1[i].is_zero()); j++){
				if (!tmp1[j].is_zero()){
					tmp = gcd(tmp1[i], tmp1[j]);
					LiDIA::divide(multiplier, tmp1[j], tmp);
					multiplier = gcd(multiplier, p);
					while (!(tmp = gcd(tmp1[i], multiplier)).is_one())
						LiDIA::divide(multiplier, multiplier, tmp);
					mult_mod(tmp1[i], tmp1[i], multiplier);
					if (tmp1[i].is_zero()) {
						ret.assign(p);
						return;
					}
					tmp = gcd(tmp1[j], tmp1[i].is_zero()?p:tmp1[i]);
					LiDIA::divide(multiplier, tmp1[i], tmp);
					for (lidia_size_t k = i + 1; k <= j; k++)
						mult_mod(A.value[k][j], A.value[k][j], multiplier);
					A.value[j][j] = gcd (A.value[j][j], p);
					if (A.value[j][j].is_zero()) {
						ret.assign(p);
						return;
					}
				}
			}
			if ((ret = lcm(ret, tmp1[i])) == p) return;
		}
	}
	else
		for (lidia_size_t i = 0; i < A.rows; i++){
			tmp1 = A.value[i];
			for (lidia_size_t j = i + 1; j < A.columns; j++){
				if (!tmp1[j].is_zero()){
					tmp = gcd(tmp1[i], tmp1[j]);
					LiDIA::divide(multiplier, tmp1[j], tmp);
					while (!(tmp = gcd(tmp1[i], multiplier)).is_one())
						LiDIA::divide(multiplier, multiplier, tmp);
					LiDIA::multiply(tmp1[i], tmp1[i], multiplier);

					tmp = gcd(tmp1[i], tmp1[j]);
					LiDIA::divide(multiplier, tmp1[i], tmp);
					for (lidia_size_t k = i + 1; k <= j; k++)
						LiDIA::multiply(A.value[k][j], A.value[k][j], multiplier);
				}
			}
			ret = lcm(ret, tmp1[i]);
		}
#endif
		for (lidia_size_t i = A.rows; i;) {
			i--;
			for (lidia_size_t j = i - 1; (j >= 0) && (!A.value[i][i].is_zero()); j--)
				if (!A.value[j][i].is_zero()) {
					tmp = gcd(A.value[i][i], A.value[j][i]);
					LiDIA::divide(multiplier, A.value[j][i], tmp);
					multiplier = gcd(multiplier, p);
					while (!(tmp = gcd(A.value[i][i], multiplier)).is_one()) {
						LiDIA::divide(multiplier, multiplier, tmp);
					}
					multiplier = gcd(multiplier, p);
					mult_mod(A.value[i][i], A.value[i][i], multiplier);
					if (A.value[i][i].is_zero()) {
						ret.assign(p);
						return;
					}
					tmp = gcd(A.value[j][i], A.value[i][i]);
					LiDIA::divide(multiplier, A.value[i][i], tmp);
					tmp1 = A.value[j];
					for (lidia_size_t k = i - 1; k >= j; k--)
						mult_mod(tmp1[k], tmp1[k], multiplier);
					tmp1[j] = gcd (tmp1[j], p);
					if (tmp1[j].is_zero()) {
						ret.assign(p);
						return;
					}
				}
			if ((ret = lcm(ret, A.value[i][i])) == p)
				return;
		}
	}
	else {
		bigint_matrix B(A.rows, A.columns);
		B.diag(1, 0);
		bigint_matrix T(A.rows + 1, A.columns);
		T.reginvimage(A, B);
		ret.assign(T.member(A.rows, 0));
		for (lidia_size_t i = 1; i < A.columns; i++)
			ret.assign(lcm(ret, T.member(A.rows, i)));
	}
}



bigint bigmod_matrix::
exponent() const
{
	debug_handler("bigmod_matrix", "in member - function exponent()");
	bigint EXP;
	exponent(EXP);
	return EXP;
}



bigint
exponent(const bigmod_matrix & A)
{
	debug_handler("bigmod_matrix", "in function "
		      "exponent(bigmod_matrix &)");
	bigint EXP;
	A.exponent(EXP);
	return EXP;

}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
