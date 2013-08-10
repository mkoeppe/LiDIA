// -*- C++ -*-
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
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_PARTIAL_RELATION_H_GUARD_
#define LIDIA_PARTIAL_RELATION_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_MATH_VECTOR_H_GUARD_
# include	"LiDIA/math_vector.h"
#endif
#ifndef LIDIA_SORT_VECTOR_H_GUARD_
# include	"LiDIA/sort_vector.h"
#endif
#ifndef LIDIA_MATRIX_H_GUARD_
# include	"LiDIA/matrix.h"
#endif
#ifndef LIDIA_QUADRATIC_NUMBER_STANDARD_H_GUARD_
# include	"LiDIA/quadratic_number_standard.h"
#endif
#include	<cstring>
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// Class: partial_relation
//
// This class represents a partial relation found from a sieving algorithm.
//      It consists of an integer representing the value of the large prime,
//      an index which distinguishes an instance of this class from other
//      instances, and a string representing the entire relation.
//

#define REL_STR_LEN 	1000
#define LINE_LEN	1000

class partial_relation
{
protected:

	int lp; // value of the large prime
	int ex; // exponent of the large prime (1 = negative)
	bigint H; // index
	char str[REL_STR_LEN]; // string representing the relation
	quadratic_number_standard g; // the principal ideal generator


public:

	partial_relation()
	{
		lp = ex = 0;
		H.assign_zero();
		memset(str, '\0', REL_STR_LEN);
		g.assign_one();
	}

	~partial_relation()
	{}

	int get_large() const
	{
		return lp;
	}

	int get_exp() const
	{
		return ex;
	}

	bigint get_index() const
	{
		return H;
	}

	quadratic_number_standard get_min() const
	{
		return g;
	}

	void strcopy(char *st)
	{
		strcpy(st, str);
	}

	void assign_zero()
	{
		lp = ex = 0;
		H.assign_zero();
		memset(str, '\0', REL_STR_LEN);
		g.assign_one();
	}

	void assign(bigint & n_H, char *n_str, quadratic_number_standard & q)
	{
		lp = 1;
		ex = 0;
		H = n_H;
		strcpy(str, n_str);
		g.assign(q);
	}

	void assign(int n_lp, int n_ex, bigint & n_H, char *n_str, quadratic_number_standard & q)
	{
		lp = n_lp;
		ex = n_ex;
		H = n_H;
		strcpy(str, n_str);
		g.assign(q);
	}

	partial_relation operator = (const partial_relation & B)
	{
		lp = B.lp;
		ex = B.ex;
		H.assign(B.H);
		strcpy(str, B.str);
		g.assign(B.g);
		return *this;
	}



	friend int compare (const partial_relation & P1, const partial_relation & P2);

	bool is_zero()
	{
		return (lp == 0);
	}

	bool is_full_relation() const
	{
		return (lp == 1);
	}



	void add_to_matrix(matrix< int > & M,
			   base_vector< quadratic_number_standard > &vec_qn,
			   lidia_size_t j);



	friend void swap(partial_relation & A, partial_relation & B);

	friend std::istream & operator >> (std::istream & in, partial_relation & A);

	friend std::ostream & operator << (std::ostream & out, const partial_relation & A);

	friend bool make_full(partial_relation & C,
			      partial_relation & A,
			      partial_relation & B,
			      bigint & N,
			      int size_exp);

};



inline void
partial_relation::add_to_matrix (matrix< int > & M,
				 base_vector< quadratic_number_standard > &vec_qn,
				 lidia_size_t j)
{

	lidia_size_t i;
	int e;
	char helpstr[LINE_LEN], *pp, *p;

	// write the relation to the last column of M
	strcpy(helpstr, str);
	pp = helpstr;
	p = strtok(pp, " \n");
	while (p != NULL) {
		e = atoi(p);
		if (!e)
			break;
		p = strtok(NULL, " \n");
		i = static_cast<lidia_size_t>(std::atoi(p));
		M.sto(i, j, e);
		p = strtok(NULL, " \n");
	}

	// write minimum to jth entry of vec
	vec_qn[j].assign(g);
}



inline std::istream &
operator >> (std::istream & in, partial_relation & A)
{
	char line[LINE_LEN], *p;


	in.getline(line, LINE_LEN);
	if (!in.eof()) {
		p = strchr(line, '@');
		if (p) {
				// reading a parital relation
			p = strtok(line, "@");
			A.lp = atoi(p);
			p = strtok(NULL, "@");
			A.ex = atoi(p);
			p = strtok(NULL, ":");
		}
		else {
				// reading a full relation
			A.lp = 1;
			A.ex = 0;
			p = strtok(line, ":");
		}
		string_to_bigint(p, A.H);
		p = strtok(NULL, ":");
		++p;
		strcpy(A.str, p);
		p = strtok(NULL, ":");
		++p;
		string_to_quadratic_number_standard(p, A.g);
	}
	else
		A.assign_zero();

	return in;
}



inline int compare (const partial_relation & P1, const partial_relation & P2)
{
	if (&P1 == &P2) {
		return 0;
	}
	if (P1.lp == P2.lp) {
		return 0;
	}
	return ((P1.lp < P2.lp) ? -1 : 1);
}



inline bool operator == (const partial_relation & P1, const partial_relation & P2)
{
	return (compare(P1, P2) == 0);
}



inline bool operator < (const partial_relation & P1, const partial_relation & P2)
{
	return (compare(P1, P2) < 0);
}



inline bool operator <= (const partial_relation & P1, const partial_relation & P2)
{
	return (compare(P1, P2) <= 0);
}



inline bool operator > (const partial_relation & P1, const partial_relation & P2)
{
	return (compare(P1, P2) > 0);
}



inline bool operator >= (const partial_relation & P1, const partial_relation & P2)
{
	return (compare(P1, P2) >= 0);
}



inline std::ostream &
operator << (std::ostream & out, const partial_relation & A)
{
	if (A.is_full_relation())
		out << A.H << " :" << A.str << ": " << A.g << "\n";
	else
		out << A.lp << " @ " << A.ex << " @ " << A.H << " : " << A.str << ": " << A.g << "\n";
	return out;
}



inline void
swap(partial_relation & A, partial_relation & B)
{
	partial_relation C;


	C = A;
	A = B;
	B = C;
}



inline bool
make_full(partial_relation & C, partial_relation & A, partial_relation & B, bigint & N, int size_exp)
{
	if (A.lp == B.lp) {
		char helpstr[REL_STR_LEN], *p, *pp;
		bigint h_inv, temp;
		int e;
		lidia_size_t i;
		bool neg, zero;
		short *exp = NULL;

		exp = new short[size_exp];
		memset(exp, '\0', size_exp*sizeof(short));

		neg = (A.ex == B.ex);

		// do exponent arithmetic
		strcpy(helpstr, A.str);
		pp = helpstr;
		p = strtok(pp, " \n");
		while (p != NULL) {
			e = atoi(p);
			if (!e)
				break;
			p = strtok(NULL, " \n");
			exp[atoi(p)] += e;
			p = strtok(NULL, " \n");
		}

		zero = true;
		strcpy(helpstr, B.str);
		pp = helpstr;
		p = strtok(pp, " \n");
		while (p != NULL) {
			e = atoi(p);
			if (!e)
				break;
			p = strtok(NULL, " \n");
			if (neg)
				exp[atoi(p)] -= e;
			else
				exp[atoi(p)] += e;
			if (exp[atoi(p)])
				zero = false;
			p = strtok(NULL, " \n");
		}

		// exit if relation is zero
		if (zero) {
			C.assign_zero();
			delete [] exp;
			return true;
		}

		temp = xgcd_left(h_inv, B.H, N);
		C.H = (h_inv * A.H) % N;

		C.ex = 0;
		C.lp = 1;
		memset(C.str, '\0', REL_STR_LEN);
		for (i = 0; i < size_exp; i++) {
			if (exp[i])
				sprintf(C.str, "%s %d %d", C.str, exp[i], i);
		}
		sprintf(C.str, "%s 0", C.str);

		if (N.is_gt_zero()) {
			if (neg)
				divide(C.g, A.g, B.g);
			else
				multiply(C.g, A.g, B.g);
		}
		else
			C.g.assign_one();

		delete [] exp;

		return true;
	}
	else
		return false;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_PARTIAL_RELATION_H_GUARD_
