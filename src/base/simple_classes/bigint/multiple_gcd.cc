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
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/bigint.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bigint *
mgcd(const bigint * a, const int n)
{
	debug_handler("matrix_bigint", "in function "
		      "mgcd(const bigint *, const int)");

	// Step 1 - 8
	int i;
	bigint TMP, TMP1, TMP2;

	bigint *tmp = new bigint[n];
	memory_handler(tmp, "matrix_bigint", "mgcd :: "
		       "Error in memory allocation (tmp)");
	bigint *y = new bigint[n];
	memory_handler(y, "matrix_bigint", "mgcd :: "
		       "Error in memory allocation (y)");
	bigint *z = new bigint[n];
	memory_handler(z, "matrix_bigint", "mgcd :: "
		       "Error in memory allocation (z)");

	tmp[0].assign(a[0]);
	for (i = 1; i < n; i++)
		tmp[i] = xgcd(y[i], z[i], tmp[i - 1], a[i]);

	// Step 9 - 19
	bigint *x = new bigint[n];
	memory_handler(x, "matrix_bigint", "mgcd :: "
		       "Error in memory allocation (x)");
	bigint *y1 = new bigint[n];
	memory_handler(y1, "matrix_bigint", "mgcd :: "
		       "Error in memory allocation (y1)");

	x[n - 1].assign(z[n - 1]);
	y1[n - 1].assign(y[n - 1]);

	bigint *G = new bigint[n];
	memory_handler(G, "matrix_bigint", "mgcd :: "
		       "Error in memory allocation (G)");
	bigint *S = new bigint[n];
	memory_handler(S, "matrix_bigint", "mgcd :: "
		       "Error in memory allocation (S)");
	for (i = n - 2; i >= 1; i--) {
		div_rem(G[i], TMP, tmp[i - 1], tmp[i]);
		nearest(S[i], y1[i + 1] * z[i], G[i]);

		// y1[i] = y1[i+1]*y[i] + S[i]*(a[i]/tmp[i]);
		div_rem(TMP1, TMP, a[i], tmp[i]);
		multiply(TMP1, S[i], TMP1);
		multiply(TMP2, y1[i + 1], y[i]);
		add(y1[i], TMP2, TMP1);

		// x[i] = y1[i+1]*z[i] - S[i]*G[i];
		multiply(TMP1, S[i], G[i]);
		multiply(TMP2, y1[i + 1], z[i]);
		subtract(x[i], TMP2, TMP1);
	}
	x[0].assign(y1[1]);

	// Step 20,21
	bigint *RES = new bigint[n + 1];
	memory_handler(RES, "matrix_bigint", "mgcd :: "
		       "Error in memory allocation (RES)");
	RES[0].assign(tmp[n - 1]);
	for (i = 1; i <= n; i++)
		RES[i].assign(x[i - 1]);

	delete[] tmp;
	delete[] y;
	delete[] z;
	delete[] x;
	delete[] y1;
	delete[] G;
	delete[] S;

	return RES;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
