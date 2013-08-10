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
//	Author	: Franz-Dieter Berger (FDB), Patric Kirsch (PK),
//                Volker Mueller(VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/finite_fields/info_gf2n.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// ==============================================================
// Function:  partial_reduce1                                    
// ==============================================================
// 
// reduce all array elements with indices >= anzBI to zero      
// Note: for complete reduction, the coefficients for x^k with  
// k > extension degree have to be reduced, too.                
// --> partial_reduce2                                          
// 
// special versions for trinomials, pentinomials                
// ==============================================================

void info_gf2n::tri_partial_reduce1(udigit *f) const
{
	register unsigned int i, k, l, w;
	register udigit h;

	i = 2*anzBI - 1;

	while (i >= anzBI) {
		h = f[i];
		if (!h) {
			i--;
			continue;
		}

		k = BITS_PER_LONG*i - degree;
		f[i--] = static_cast<udigit>(0); // x^degree 

		w = k / (BITS_PER_LONG); // x^0     
		l = k % (BITS_PER_LONG);

		if (l != 0) {
			f[w] ^= (h << l);
			f[w+1] ^= (h >> (BITS_PER_LONG-l));
		}
		else
			f[w] ^= h;

		w = (k + exp1) / (BITS_PER_LONG);
		l = (k + exp1) % (BITS_PER_LONG);

		if (l != 0) {
			f[w] ^= (h << l);
			f[w+1] ^= (h >> (BITS_PER_LONG-l));
		}
		else
			f[w] ^= h;
	}
}



void info_gf2n::pent_partial_reduce1(udigit *f) const
{
	register unsigned int i, k, l, w;
	register udigit h;

	i = 2*anzBI - 1;

	while (i >= anzBI) {
		h = f[i];
		if (!h) {
			i--;
			continue;
		}

		k = BITS_PER_LONG*i - degree;
		f[i--] = static_cast<udigit>(0); // x^degree 

		w = k / (BITS_PER_LONG); // x^0     
		l = k % (BITS_PER_LONG);

		if (l != 0) {
			f[w] ^= (h << l);
			f[w+1] ^= (h >> (BITS_PER_LONG-l));
		}
		else
			f[w] ^= h;

		w = (k + exp1) / (BITS_PER_LONG);
		l = (k + exp1) % (BITS_PER_LONG);

		if (l != 0) {
			f[w] ^= (h << l);
			f[w+1] ^= (h >> (BITS_PER_LONG-l));
		}
		else
			f[w] ^= h;

		w = (k + exp2) / (BITS_PER_LONG);
		l = (k + exp2) % (BITS_PER_LONG);

		if (l != 0) {
			f[w] ^= (h << l);
			f[w+1] ^= (h >> (BITS_PER_LONG-l));
		}
		else
			f[w] ^= h;

		w = (k + exp3) / (BITS_PER_LONG);
		l = (k + exp3) % (BITS_PER_LONG);

		if (l != 0) {
			f[w] ^= (h << l);
			f[w+1] ^= (h >> (BITS_PER_LONG-l));
		}
		else
			f[w] ^= h;
	}
}



void info_gf2n::general_partial_reduce1(udigit *f) const
{
	register unsigned int i, j, k, l, w, s;
	register udigit h;

	i = 2*anzBI - 1;

	while (i >= anzBI) {
		h = f[i];
		if (!h) {
			i--;
			continue;
		}

		k = BITS_PER_LONG*i - degree;
		f[i] = static_cast<udigit>(0);

		w = k / (BITS_PER_LONG);
		l = k % (BITS_PER_LONG);

		if (l != 0) {
			f[w] ^= (h << l);
			f[w+1] ^= (h >> (BITS_PER_LONG-l));
		}
		else
			f[w] ^= h;

		for (j = 0; j < anz_exponents; j++) {
			s = k + exponents[j];
			w = s / (BITS_PER_LONG);
			l = s % (BITS_PER_LONG);

			if (l != 0) {
				f[w] ^= (h << l);
				f[w+1] ^= (h >> (BITS_PER_LONG-l));
			}
			else
				f[w] ^= h;
		}
	}
}



// ==============================================================
// Function:  partial_reduce2                                    
// ==============================================================
// we assume that partial_reduce1 has been used  and only the    
// word with index (anzBI-1) may be unreduced                    
// 
// special versions for trinomials, pentinomials                
// ==============================================================


void info_gf2n::tri_partial_reduce2(udigit *f) const
{
	register unsigned int l, w, deg, anz;
	register udigit h;

	anz = anzBI-1;
	deg = degree % (BITS_PER_LONG);
	if (deg != 0)
		h = f[anz] >> deg;
	else
		h = f[anz];

	if (h == static_cast<udigit>(0))   // deg(f) < deg(modulus) 
		return;

	f[0] ^= h;
	f[anz] ^= (h << deg);

	w = exp1 / (BITS_PER_LONG);
	l = exp1 % (BITS_PER_LONG);

	if (l != 0) {
		f[w] ^= (h << l);
		f[w+1] ^= (h >> (BITS_PER_LONG-l));
	}
	else
		f[w] ^= h;
}



void info_gf2n::pent_partial_reduce2(udigit *f) const
{
	register unsigned int l, w, anz, deg;
	register udigit h;

	anz = anzBI-1;
	deg = degree % (BITS_PER_LONG);
	if (deg != 0)
		h = f[anz] >> deg;
	else
		h = f[anz];

	if (h == static_cast<udigit>(0))   // deg(f) < deg(modulus) 
		return;

	f[0] ^= h;
	f[anz] ^= (h << deg);

	w = exp1 / (BITS_PER_LONG);
	l = exp1 % (BITS_PER_LONG);

	if (l != 0) {
		f[w] ^= (h << l);
		f[w+1] ^= (h >> (BITS_PER_LONG-l));
	}
	else
		f[w] ^= h;

	w = exp2 / (BITS_PER_LONG);
	l = exp2 % (BITS_PER_LONG);

	if (l != 0) {
		f[w] ^= (h << l);
		f[w+1] ^= (h >> (BITS_PER_LONG-l));
	}
	else
		f[w] ^= h;

	w = exp3 / (BITS_PER_LONG);
	l = exp3 % (BITS_PER_LONG);

	if (l != 0) {
		f[w] ^= (h << l);
		f[w+1] ^= (h >> (BITS_PER_LONG-l));
	}
	else
		f[w] ^= h;
}



void info_gf2n::general_partial_reduce2(udigit *f) const
{
	register unsigned int j, s, l, w, anz, deg;
	register udigit h;

	anz = anzBI-1;
	deg = degree % (BITS_PER_LONG);
	if (deg != 0)
		h = f[anz] >> deg;
	else
		h = f[anz];


	if (h == static_cast<udigit>(0))   // deg(f) < deg(modulus) 
		return;

	while (h != 0) {
		f[0] ^= h;
		f[anz] ^= (h << deg);

		for (j = 0; j < anz_exponents; j++) {
			s = exponents[j];
			w = s / (BITS_PER_LONG);
			l = s % (BITS_PER_LONG);

			if (l != 0) {
				f[w] ^= (h << l);

				if (w < anz)
				  f[w+1] ^= (h >> (BITS_PER_LONG-l));
			}
			else
				f[w] ^= h;
		}
		h = f[anz] >> deg;
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
