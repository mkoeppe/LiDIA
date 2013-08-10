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
#include	"LiDIA/gf2n.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



unsigned int gf2n::exp1 = 0, gf2n::exp2 = 0, gf2n::exp3 = 0; // set in gf2n_init



// ==============================================================
// Function:  partial_reduce1                                    
// ==============================================================
// 
// reduce all array elements with indices >= gf2n::anzBI to zero
// Note: for complete reduction, the coefficients for x^k with  
// k > extension degree have to be reduced, too.                
// --> partial_reduce2                                          
// 
// special versions for trinomials, pentinomials                
// ==============================================================

void tri_partial_reduce1(gf2n_word *f)
{
	register unsigned int i, k, l, w;
	register unsigned int anzBI = gf2n::anzBI, degree = gf2n::degree;
	register gf2n_word h;

	i = 2*anzBI - 1;

	while (i >= anzBI) {
		h = f[i];
		if (!h) {
			i--;
			continue;
		}

		k = CHAR_BIT*sizeof(gf2n_word)*i - degree;
		f[i--] = static_cast<gf2n_word>(0); // x^degree 

		w = k / (CHAR_BIT*sizeof(gf2n_word)); // x^0     
		l = k % (CHAR_BIT*sizeof(gf2n_word));

		if (l != 0) {
			f[w] ^= (h << l);
			f[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
		}
		else
			f[w] ^= h;

		w = (k + gf2n::exp1) / (CHAR_BIT*sizeof(gf2n_word));
		l = (k + gf2n::exp1) % (CHAR_BIT*sizeof(gf2n_word));

		if (l != 0) {
			f[w] ^= (h << l);
			f[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
		}
		else
			f[w] ^= h;
	}
}



void pent_partial_reduce1(gf2n_word *f)
{
	register unsigned int i, k, l, w;
	register unsigned int anzBI = gf2n::anzBI, degree = gf2n::degree;
	register gf2n_word h;

	i = 2*anzBI - 1;

	while (i >= anzBI) {
		h = f[i];
		if (!h) {
			i--;
			continue;
		}

		k = CHAR_BIT*sizeof(gf2n_word)*i - degree;
		f[i--] = static_cast<gf2n_word>(0); // x^degree 

		w = k / (CHAR_BIT*sizeof(gf2n_word)); // x^0     
		l = k % (CHAR_BIT*sizeof(gf2n_word));

		if (l != 0) {
			f[w] ^= (h << l);
			f[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
		}
		else
			f[w] ^= h;

		w = (k + gf2n::exp1) / (CHAR_BIT*sizeof(gf2n_word));
		l = (k + gf2n::exp1) % (CHAR_BIT*sizeof(gf2n_word));

		if (l != 0) {
			f[w] ^= (h << l);
			f[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
		}
		else
			f[w] ^= h;

		w = (k + gf2n::exp2) / (CHAR_BIT*sizeof(gf2n_word));
		l = (k + gf2n::exp2) % (CHAR_BIT*sizeof(gf2n_word));

		if (l != 0) {
			f[w] ^= (h << l);
			f[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
		}
		else
			f[w] ^= h;

		w = (k + gf2n::exp3) / (CHAR_BIT*sizeof(gf2n_word));
		l = (k + gf2n::exp3) % (CHAR_BIT*sizeof(gf2n_word));

		if (l != 0) {
			f[w] ^= (h << l);
			f[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
		}
		else
			f[w] ^= h;
	}
}



void general_partial_reduce1(gf2n_word *f)
{
	register unsigned int i, j, k, l, w, s;
	register unsigned int anzBI = gf2n::anzBI, degree = gf2n::degree;
	register unsigned int anz_expo = gf2n::anz_exponents;
	register gf2n_word h;

	i = 2*anzBI - 1;

	while (i >= anzBI) {
		h = f[i];
		if (!h) {
			i--;
			continue;
		}

		k = CHAR_BIT*sizeof(gf2n_word)*i - degree;
		f[i] = static_cast<gf2n_word>(0);

		w = k / (CHAR_BIT*sizeof(gf2n_word));
		l = k % (CHAR_BIT*sizeof(gf2n_word));

		if (l != 0) {
			f[w] ^= (h << l);
			f[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
		}
		else
			f[w] ^= h;

		for (j = 0; j < anz_expo; j++) {
			s = k + gf2n::exponents[j];
			w = s / (CHAR_BIT*sizeof(gf2n_word));
			l = s % (CHAR_BIT*sizeof(gf2n_word));

			if (l != 0) {
				f[w] ^= (h << l);
				f[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
			}
			else
				f[w] ^= h;
		}
	}
}



//******* partial_reduce1 - selector *********************************

void (*partial_reduce1[]) (gf2n_word*) =
{
	tri_partial_reduce1, pent_partial_reduce1, general_partial_reduce1
};


// ==============================================================
// Function:  partial_reduce2                                    
// ==============================================================
// we assume that partial_reduce1 has been used  and only the    
// word with index (gf2n::anzBI-1) may be unreduced              
// 
// special versions for trinomials, pentinomials                
// ==============================================================


void tri_partial_reduce2(gf2n_word *f)
{
	register unsigned int l, w, deg, anz;
	register gf2n_word h;

	anz = gf2n::anzBI-1;
	deg = gf2n::degree % (CHAR_BIT*sizeof(gf2n_word));
	if (deg != 0)
		h = f[anz] >> deg;
	else
		h = f[anz];

	if (h == static_cast<gf2n_word>(0))   // deg(f) < deg(modulus) 
		return;

	f[0] ^= h;
	f[anz] ^= (h << deg);

	w = gf2n::exp1 / (CHAR_BIT*sizeof(gf2n_word));
	l = gf2n::exp1 % (CHAR_BIT*sizeof(gf2n_word));

	if (l != 0) {
		f[w] ^= (h << l);
		// patch to avoid illegal memory access -rpw
		if (w < anz)
			f[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
	}
	else
		f[w] ^= h;
}



void pent_partial_reduce2(gf2n_word *f)
{
	register unsigned int l, w, anz, deg;
	register gf2n_word h;

	anz = gf2n::anzBI-1;
	deg = gf2n::degree % (CHAR_BIT*sizeof(gf2n_word));
	if (deg != 0)
		h = f[anz] >> deg;
	else
		h = f[anz];

	if (h == static_cast<gf2n_word>(0))   // deg(f) < deg(modulus) 
		return;

	f[0] ^= h;
	f[anz] ^= (h << deg);

	w = gf2n::exp1 / (CHAR_BIT*sizeof(gf2n_word));
	l = gf2n::exp1 % (CHAR_BIT*sizeof(gf2n_word));

	if (l != 0) {
		f[w] ^= (h << l);
		// patch -- rpw
		if (w < anz)
			f[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
	}
	else
		f[w] ^= h;

	w = gf2n::exp2 / (CHAR_BIT*sizeof(gf2n_word));
	l = gf2n::exp2 % (CHAR_BIT*sizeof(gf2n_word));

	if (l != 0) {
		f[w] ^= (h << l);
		// patch -- rpw
		if (w < anz)
			f[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
	}
	else
		f[w] ^= h;

	w = gf2n::exp3 / (CHAR_BIT*sizeof(gf2n_word));
	l = gf2n::exp3 % (CHAR_BIT*sizeof(gf2n_word));

	if (l != 0) {
		f[w] ^= (h << l);
		// patch -- rpw
		if (w < anz)
			f[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
	}
	else
		f[w] ^= h;
}



void general_partial_reduce2(gf2n_word *f)
{
	register unsigned int j, s, l, w, anz, deg;
	register gf2n_word h;

	anz = gf2n::anzBI-1;
	deg = gf2n::degree % (CHAR_BIT*sizeof(gf2n_word));
	if (deg != 0)
		h = f[anz] >> deg;
	else
		h = f[anz];


	if (h == static_cast<gf2n_word>(0))   // deg(f) < deg(modulus) 
		return;

	while (h != 0) {
		f[0] ^= h;
		f[anz] ^= (h << deg);

		for (j = 0; j < gf2n::anz_exponents; j++) {
			s = gf2n::exponents[j];
			w = s / (CHAR_BIT*sizeof(gf2n_word));
			l = s % (CHAR_BIT*sizeof(gf2n_word));

			if (l != 0) {
				f[w] ^= (h << l);
				// patch -- rpw
				if (w < anz)
					f[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
			}
			else
				f[w] ^= h;
		}
		h = f[anz] >> deg;
	}
}



//******* partial_reduce2 - selector *********************************

void (*partial_reduce2[]) (gf2n_word*) =
{
	tri_partial_reduce2, pent_partial_reduce2, general_partial_reduce2
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
