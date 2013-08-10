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


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/number_fields/qo_util.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// MyTime
//
// Task:
//      outputs a time in my own special format.
//

void
MyTime(long t)
{
	long d, h, m, s, ms;

	ms = t % 100;
	t /= 100;
	s = t % 60;
	t /= 60;
	m = t % 60;
	t /= 60;
	h = t % 24;
	d = t / 24;

	if (d > 0) {
		std::cout << d << " day, ";
		std::cout << h << " hour, ";
		std::cout << m << " min, ";
		std::cout << s << " sec, ";
		std::cout << ms << " csec";
	}
	else if (h > 0) {
		std::cout << h << " hour, ";
		std::cout << m << " min, ";
		std::cout << s << " sec, ";
		std::cout << ms << " csec";
	}
	else if (m > 0) {
		std::cout << m << " min, ";
		std::cout << s << " sec, ";
		std::cout << ms << " csec";
	}
	else if (s > 0) {
		std::cout << s << " sec, ";
		std::cout << ms << " csec";
	}
	else
		std::cout << ms << " csec";

#if 0
	if (d > 0) {
		std::cout << setw(3) << d << " day, ";
		std::cout << setw(3) << h << " hour, ";
		std::cout << setw(3) << m << " min, ";
		std::cout << setw(3) << s << " sec, ";
		std::cout << setw(3) << ms << " hsec";
	}
	else if (h > 0) {
		std::cout << setw(3) << h << " hour, ";
		std::cout << setw(3) << m << " min, ";
		std::cout << setw(3) << s << " sec, ";
		std::cout << setw(3) << ms << " hsec";
	}
	else if (m > 0) {
		std::cout << setw(3) << m << " min, ";
		std::cout << setw(3) << s << " sec, ";
		std::cout << setw(3) << ms << " hsec";
	}
	else if (s > 0) {
		std::cout << setw(3) << s << " sec, ";
		std::cout << setw(3) << ms << " hsec";
	}
	else
		std::cout << setw(3) << ms << " hsec";
#endif

	return;
};



//
// decode_vector()
//
// Task:
//      given indices r and q from the sets R and Q, decodes the respective
//      exponent vectors.
//

void
decode_vector(matrix< bigint > & Bmat, const bigint & Bjj, long r,
              long q, base_vector< long > & Rvec, base_vector< long > & Qvec,
              long nR, long nQ)
{
	debug_handler("qi_class", "decode_vector");

	lidia_size_t i, rank;
	bigint temp;
	base_vector< bigint > Bj;

	// compute powers in index vector
	rank = Bmat.get_no_of_columns();

	if (rank == 0) {
		Bmat.resize(1, 1);
		Bmat.sto(0, 0, Bjj);
	}
	else {
		Bj.set_mode(EXPAND);
		Bj.reset();

		Bj[rank] = Bjj;
		for (i = rank-1; i >= 0; --i) {
			nR /= Rvec[i];
			nQ /= Qvec[i];
			temp.assign(q/nQ);
			multiply(temp, temp, Rvec[i]);
			add(Bj[i], temp, (r/nR));
			r %= nR;
			q %= nQ;
		}

		// new row and column
		++rank;
		Bmat.set_no_of_rows(rank);
		Bmat.set_no_of_columns(rank);
		Bmat.sto_column_vector(Bj, rank, rank-1);
	}
}



bigint
int_key(const int & G)
{
	return bigint(G);
}



bigint
long_key(const long & G)
{
	return bigint(G);
}



bigint
bigint_key(const bigint & G)
{
	return G;
}



bigint
pair_bigint_key(const pair< bigint, bigint > & G)
{
	return G.left();
}



bigint
bigrational_key(const bigrational & G)
{
	return G.denominator() + G.numerator();
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
