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



void tri_partial_reduce2(gf2n_word *);
void pent_partial_reduce2(gf2n_word *);
extern void (*partial_reduce2[]) (gf2n_word *);

// simple macro for calculating 2^n as a gf2n_word.
// this looks much cleaner that having the shift operations all over
// the code and fixes the missing typecast problems as well.
// -rpw


#define EXP2(x)		((static_cast<gf2n_word>(1)) << x)



inline gf2n_word max(gf2n_word a, gf2n_word b)
{
	if (a >= b)
		return a;
	else return b;
}



//************************************************************************
// shift F k bits to the left, adjust Flast
//************************************************************************

inline void shift_left(gf2n_word *F, unsigned int k, unsigned int & Flast)
{
	register int i, j, s;
	register gf2n_word *fp1, *fp2;

	s = k % (CHAR_BIT*sizeof(gf2n_word));
	j = k / (CHAR_BIT*sizeof(gf2n_word));

	if (s == 0) {
		for (i = Flast, fp1 = F+Flast+j, fp2 = F+Flast; i >= 0; i--, fp1--, fp2--)
			*fp1 = *fp2;
		Flast += j;
	}
	else {

		for (i = Flast, fp1 = F+Flast+j+1, fp2 = F+Flast; i >= 0; i--, fp2--) {
			(*fp1) |= ((*fp2) >> (CHAR_BIT*sizeof(gf2n_word)-s));
			fp1--;
			(*fp1) = ((*fp2) << s);
		}
		Flast += j+1;
	}

	i = 0; fp1 = F;
	while (i < j) {
		*fp1 = static_cast<gf2n_word>(0);
		i++; fp1++;
	}

	fp1 = F + Flast;
	while ((*fp1) == static_cast<gf2n_word>(0) && Flast > 0) {
		Flast --;
		fp1--;
	}
}



//************************************************************************
// shift F k bits to the right, adjust Flast
//************************************************************************

inline void shift_right(gf2n_word *F, unsigned int k, unsigned int & Flast)
{
	register unsigned int j, i, s;
	register gf2n_word *fp1, *fp2;

	s = k % (CHAR_BIT*sizeof(gf2n_word));
	j = k / (CHAR_BIT*sizeof(gf2n_word));

	if (s == 0) {
		for (i = j, fp1 = F, fp2 = F+j; i <= Flast; i++, fp1++, fp2++)
			*fp1 = *fp2;
	}
	else {
		// patch to avoid illegal memory access -rpw
		for (i = j, fp1 = F, fp2 = F+j; i < Flast; i++, fp1++, fp2++)
			*fp1 = ((*fp2) >> s) | ((*(fp2+1)) << (CHAR_BIT*sizeof(gf2n_word)-s));
		*fp1 = ((*fp2) >> s);
	}
	i = Flast-j+1;
	fp1 = F+i;
	while (i <= Flast) {
		*fp1 = static_cast<gf2n_word>(0);
		fp1++;
		i++;
	}

	Flast -= j;
	fp1 = F+Flast;
	while (*fp1 == static_cast<gf2n_word>(0) && Flast > 0) {
		Flast --; fp1--;
	}
}



//************************************************************************
// compute F = (F + G) / x^k, adjust Flast
//************************************************************************

inline void add_shift_right(gf2n_word *F, gf2n_word*G, unsigned int k,
			    unsigned int & Flast, unsigned int Glast)

{
	register unsigned int j, i, s, l;
	register gf2n_word h;
	register gf2n_word *fp1, *fp2, *gp;

	s = k % (CHAR_BIT*sizeof(gf2n_word));
	j = k / (CHAR_BIT*sizeof(gf2n_word));

	l = static_cast<unsigned int>(max(Flast, Glast));

	if (s == 0) {
		for (i = j, fp1 = F, fp2 = F+j, gp = G+j; i <= l;
		     i++, fp1++, fp2++, gp++)
			*fp1 = (*fp2) ^ (*gp);
	}
	else {
		fp2 = F+j; gp = G+j;
		h = (*fp2) ^ (*gp);
		for (i = j, fp1 = F; i <= l; i++, fp1++) {
			*fp1 = (h >> s);
			fp2++; gp++;
			h = (*fp2) ^ (*gp);
			*fp1 |= (h << (CHAR_BIT*sizeof(gf2n_word)-s));
		}
	}
	i = l-j+1;
	fp1 = F+i;
	while (i <= Flast) {
		*fp1 = static_cast<gf2n_word>(0);
		i++; fp1++;
	}


	Flast = l - j;
	fp1 = F+Flast;
	while (*fp1 == static_cast<gf2n_word>(0) && Flast > 0) {
		Flast --;
		fp1--;
	}
}



//*********************************************************************
// invert function, res and a are assumed to have length gf2n::anzBI
//*********************************************************************

void tri_invert(gf2n_word *res, gf2n_word *a)
{
	register gf2n_word h;
	register gf2n_word *bp, *cp, *fp, *gp, *ap;
	unsigned int i, j, w, l;
	unsigned int xdegree = 0;
	unsigned int anzBI = gf2n::anzBI;
	unsigned int Clast, Blast, Flast, Glast;

	tri_partial_reduce2(a); // make sure that a is totally reduced 

	bp = gf2n::B; cp = gf2n::C; fp = gf2n::F; gp = gf2n::G; ap = a;

	for (i = 0; i < anzBI; i++, fp++, gp++, cp++, bp++, ap++) {
                // initialize B=1, C = 0, F = a, G = modulus, res = 0
		*fp = *ap;
		*gp = *cp = *bp = static_cast<gf2n_word>(0);
	}


	gf2n::B[0] = static_cast<gf2n_word>(1);

	for (i = anzBI, bp = gf2n::B+anzBI, cp = gf2n::C+anzBI, fp = gf2n::F+anzBI,
		     gp = gf2n::G+anzBI; i < 2*anzBI; i++, bp++, cp++, fp++, gp++) {
		*fp = *gp = *cp = *bp = static_cast<gf2n_word>(0);
	}


	gf2n::G[0] = static_cast<gf2n_word>(1);
	i = gf2n::degree;
	gf2n::G[i/(CHAR_BIT*sizeof(gf2n_word))] ^= EXP2(i % (CHAR_BIT*sizeof(gf2n_word)));
	gf2n::G[gf2n::exp1/(CHAR_BIT*sizeof(gf2n_word))] ^= EXP2(gf2n::exp1 % (CHAR_BIT*sizeof(gf2n_word)));


	Clast = Blast = 0; // initialize last non zero indices
	Flast = Glast = anzBI-1;

	fp = gf2n::F+Flast;
	while (*fp == static_cast<gf2n_word>(0) && Flast > 0) {
		Flast --; fp--;
	}

	//*************************************************************
	// start with gcd-computation, use the invariant:
	//
	// B * F + ?? * m = F
	// C * F + ?? * m = G
	//
	//*************************************************************


	// make F odd, if necessary
	h = gf2n::F[0];
	fp = gf2n::F;

	if (!(h & 1)) {
		while (h == static_cast<gf2n_word>(0)) {
			fp++;
			xdegree += CHAR_BIT*sizeof(gf2n_word);
			h = *fp;
		}

		while (!(h & 1)) {
			xdegree++;
			h >>= 1;
		}
		shift_right(gf2n::F, xdegree, Flast);

		if (Flast == 0 && gf2n::F[0] == 1)
			goto final_step;
	}

	while (true) {
		// degree(G) > degree(F)
		do {
			j = static_cast<unsigned int>(max(Clast, Blast));
			for (i = 0, cp = gf2n::C, bp = gf2n::B; i <= j; i++, cp++, bp++)   // C = C+B
				(*cp) ^= (*bp);

			Clast = j;
			cp--;
			while (*cp == static_cast<gf2n_word>(0) && Clast > 0) {
				Clast --;
				cp--;
			}

			h = gf2n::G[0] ^ gf2n::F[0];
			i = 0;

			while (h == static_cast<gf2n_word>(0)) {
				xdegree += CHAR_BIT*sizeof(gf2n_word);
				i += CHAR_BIT*sizeof(gf2n_word);
				h = gf2n::G[i/(CHAR_BIT*sizeof(gf2n_word))] ^ gf2n::F[i/(CHAR_BIT*sizeof(gf2n_word))];
			}

			while (!(h & 1)) {
				xdegree++;
				i++;
				h >>= 1;
			}
			add_shift_right(gf2n::G, gf2n::F, i, Glast, Flast);
			shift_left(gf2n::B, i, Blast);

			if (Glast == 0 && gf2n::G[0] == static_cast<gf2n_word>(1))          // G == 1  ??
				goto final_step;
		} while (Glast > Flast || ((gf2n::G[Glast] >= gf2n::F[Flast]) && (Glast == Flast)));

		do {
			// deg F > deg G 
			j = static_cast<unsigned int>(max(Clast, Blast));
			for (i = 0, bp = gf2n::B, cp = gf2n::C; i <= j; i++, bp++, cp++)   // B = B + C
				(*bp) ^= (*cp);

			Blast = j;
			bp--;
			while (*bp == 0 && Blast > 0) {
				Blast --;
				bp--;
			}

			i = 0; // make F odd
			h = gf2n::F[0] ^ gf2n::G[0];

			while (h == static_cast<gf2n_word>(0)) {
				xdegree += CHAR_BIT*sizeof(gf2n_word);
				i += CHAR_BIT*sizeof(gf2n_word);
				h = gf2n::F[i/(CHAR_BIT*sizeof(gf2n_word))] ^ gf2n::G[i/(CHAR_BIT*sizeof(gf2n_word))];
			}

			while (!(h & 1)) {
				xdegree++;
				i ++;
				h >>= 1;
			}
			add_shift_right(gf2n::F, gf2n::G, i, Flast, Glast);
			shift_left(gf2n::C, i, Clast);

			if (Flast == 0 && gf2n::F[0] == static_cast<gf2n_word>(1))   // F == 1
				goto final_step;
		} while (Flast > Glast || ((gf2n::F[Flast] >= gf2n::G[Glast]) && (Flast == Glast)));
	}


 final_step:       // check whether F == 1 mod modulus or G == 1 mod modulus

	fp = gf2n::B;

	if (gf2n::G[0] == static_cast<gf2n_word>(1) && Glast == 0) {
		gf2n::B = gf2n::C;
		Blast = Clast;
	}

	// now B == x^(xdegree) mod modulus and B has last index Blast

	while (xdegree >= CHAR_BIT*sizeof(gf2n_word)) {
		xdegree -= CHAR_BIT*sizeof(gf2n_word);
		h = gf2n::B[0];

		gf2n::B[0] = 0;
		w = gf2n::degree / (CHAR_BIT*sizeof(gf2n_word));
		l = gf2n::degree % (CHAR_BIT*sizeof(gf2n_word));
		Blast = w;

		if (l != 0) {
			gf2n::B[w] ^= (h << l);
			gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
			Blast++;
		}
		else
			gf2n::B[w] ^= h;

		w = gf2n::exp1 / (CHAR_BIT*sizeof(gf2n_word));
		l = gf2n::exp1 % (CHAR_BIT*sizeof(gf2n_word));

		if (l != 0) {
			gf2n::B[w] ^= (h << l);
			gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
		}
		else
			gf2n::B[w] ^= h;

		shift_right(gf2n::B, CHAR_BIT*sizeof(gf2n_word), Blast);
	}

	h = gf2n::B[0] & (EXP2(xdegree)-1);
	gf2n::B[0] ^= h;
	w = gf2n::degree / (CHAR_BIT*sizeof(gf2n_word));
	l = gf2n::degree % (CHAR_BIT*sizeof(gf2n_word));
	Blast = w;
	if (l != 0) {
		gf2n::B[w] ^= (h << l);
		gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
		Blast++;
	}
	else gf2n::B[w] ^= h;

	w = gf2n::exp1 / (CHAR_BIT*sizeof(gf2n_word));
	l = gf2n::exp1 % (CHAR_BIT*sizeof(gf2n_word));

	if (l != 0) {
		gf2n::B[w] ^= (h << l);
		gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
	}
	else gf2n::B[w] ^= h;

	shift_right(gf2n::B, xdegree, Blast);

	for (i = 0, bp = gf2n::B, cp = res; i <= Blast; i++, bp++, cp++)  // copy B into res
		*cp = *bp;

	if (Blast < anzBI-1)
		for (i = Blast+1; i < anzBI; i++, cp++)
			*cp = static_cast<gf2n_word>(0);

	gf2n::B = fp;
}



//******************************************************************

void pent_invert(gf2n_word *res, gf2n_word *a)
{
	register gf2n_word h;
	register gf2n_word *bp, *cp, *fp, *gp, *ap;
	unsigned int i, j, w, l;
	unsigned int xdegree = 0;
	unsigned int anzBI = gf2n::anzBI;
	unsigned int Clast, Blast, Flast, Glast;

	pent_partial_reduce2(a); // make sure that a is totally reduced 

	bp = gf2n::B; cp = gf2n::C; fp = gf2n::F; gp = gf2n::G; ap = a;

	for (i = 0, ap = a; i < anzBI; i++, fp++, gp++, cp++, bp++, ap++) {
                // initialize B=1, C = 0, F = a, G = modulus, res = 0
		*fp = *ap;
		*gp = *bp = *cp = static_cast<gf2n_word>(0);
	}

	gf2n::B[0] = static_cast<gf2n_word>(1);

	for (i = anzBI, bp = gf2n::B+anzBI, cp = gf2n::C+anzBI, fp = gf2n::F+anzBI, gp = gf2n::G+anzBI;
	     i < 2*anzBI; i++, bp++, cp++, fp++, gp++) {
		*fp = *gp = *cp = *bp = static_cast<gf2n_word>(0);
	}

	gf2n::G[0] = static_cast<gf2n_word>(1);
	i = gf2n::degree;
	gf2n::G[i/(CHAR_BIT*sizeof(gf2n_word))] ^= EXP2(i % (CHAR_BIT*sizeof(gf2n_word)));
	gf2n::G[gf2n::exp1/(CHAR_BIT*sizeof(gf2n_word))] ^= EXP2(gf2n::exp1 % (CHAR_BIT*sizeof(gf2n_word)));
	gf2n::G[gf2n::exp2/(CHAR_BIT*sizeof(gf2n_word))] ^= EXP2(gf2n::exp2 % (CHAR_BIT*sizeof(gf2n_word)));
	gf2n::G[gf2n::exp3/(CHAR_BIT*sizeof(gf2n_word))] ^= EXP2(gf2n::exp3 % (CHAR_BIT*sizeof(gf2n_word)));

	Clast = Blast = 0; // initialize last non zero indices
	Flast = Glast = anzBI-1;

	fp = gf2n::F+Flast;
	while (*fp == static_cast<gf2n_word>(0) && Flast > 0) {
		Flast --; fp--;
	}

	//*************************************************************
	// start with gcd-computation, use the invariant:
	//
	// B * F + ?? * m = F
	// C * F + ?? * m = G
	//
	//*************************************************************


	// make F odd, if necessary
	h = gf2n::F[0];
	fp = gf2n::F;

	if (!(h & 1)) {
		while (h == static_cast<gf2n_word>(0)) {
			fp++;
			xdegree += CHAR_BIT*sizeof(gf2n_word);
			h = *fp;
		}

		while (!(h & 1)) {
			xdegree++;
			h >>= 1;
		}
		shift_right(gf2n::F, xdegree, Flast);

		if (Flast == 0 && gf2n::F[0] == 1)
			goto final_step;
	}



	while (true) {
		// degree(G) > degree(F)
		do {
			j = static_cast<unsigned int>(max(Clast, Blast));
			for (i = 0, cp = gf2n::C, bp = gf2n::B; i <= j; i++, cp++, bp++)   // C = C+B
				(*cp) ^= (*bp);

			Clast = j;
			cp--;
			while (*cp == static_cast<gf2n_word>(0) && Clast > 0) {
				Clast --;
				cp--;
			}

			h = gf2n::G[0] ^ gf2n::F[0];
			i = 0;

			while (h == static_cast<gf2n_word>(0)) {
				xdegree += CHAR_BIT*sizeof(gf2n_word);
				i += CHAR_BIT*sizeof(gf2n_word);
				h = gf2n::G[i/(CHAR_BIT*sizeof(gf2n_word))] ^ gf2n::F[i/(CHAR_BIT*sizeof(gf2n_word))];
			}

			while (!(h & 1)) {
				xdegree++;
				i++;
				h >>= 1;
			}
			add_shift_right(gf2n::G, gf2n::F, i, Glast, Flast);
			shift_left(gf2n::B, i, Blast);

			if (Glast == 0 && gf2n::G[0] == static_cast<gf2n_word>(1))          // G == 1  ??
				goto final_step;
		} while (Glast > Flast || ((gf2n::G[Glast] >= gf2n::F[Flast]) && (Glast == Flast)));

		do {
			// deg F > deg G 
			j = static_cast<unsigned int>(max(Clast, Blast));
			for (i = 0, bp = gf2n::B, cp = gf2n::C; i <= j; i++, bp++, cp++)   // B = B + C
				(*bp) ^= (*cp);

			Blast = j;
			bp--;
			while (*bp == 0 && Blast > 0) {
				Blast --;
				bp--;
			}

			i = 0; // make F odd
			h = gf2n::F[0] ^ gf2n::G[0];

			while (h == static_cast<gf2n_word>(0)) {
				xdegree += CHAR_BIT*sizeof(gf2n_word);
				i += CHAR_BIT*sizeof(gf2n_word);
				h = gf2n::F[i/(CHAR_BIT*sizeof(gf2n_word))] ^ gf2n::G[i/(CHAR_BIT*sizeof(gf2n_word))];
			}

			while (!(h & 1)) {
				xdegree++;
				i ++;
				h >>= 1;
			}
			add_shift_right(gf2n::F, gf2n::G, i, Flast, Glast);
			shift_left(gf2n::C, i, Clast);

			if (Flast == 0 && gf2n::F[0] == static_cast<gf2n_word>(1))   // F == 1
				goto final_step;
		}
		while (Flast > Glast || ((gf2n::F[Flast] >= gf2n::G[Glast]) && (Flast == Glast)));
	}


 final_step:       // check whether F == 1 mod modulus or G == 1 mod modulus

	fp = gf2n::B;

	if (gf2n::G[0] == static_cast<gf2n_word>(1) && Glast == 0) {
		gf2n::B = gf2n::C;
		Blast = Clast;
	}
	// now B == x^(xdegree) mod modulus and B has last index Blast

	while (xdegree >= CHAR_BIT*sizeof(gf2n_word)) {
		xdegree -= CHAR_BIT*sizeof(gf2n_word);
		h = gf2n::B[0];

		gf2n::B[0] = 0;
		w = gf2n::degree / (CHAR_BIT*sizeof(gf2n_word));
		l = gf2n::degree % (CHAR_BIT*sizeof(gf2n_word));
		Blast = w;

		if (l != 0) {
			gf2n::B[w] ^= (h << l);
			gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
			Blast++;
		}
		else
			gf2n::B[w] ^= h;

		w = gf2n::exp1 / (CHAR_BIT*sizeof(gf2n_word));
		l = gf2n::exp1 % (CHAR_BIT*sizeof(gf2n_word));

		if (l != 0) {
			gf2n::B[w] ^= (h << l);
			gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
		}
		else
			gf2n::B[w] ^= h;

		w = gf2n::exp2 / (CHAR_BIT*sizeof(gf2n_word));
		l = gf2n::exp2 % (CHAR_BIT*sizeof(gf2n_word));

		if (l != 0) {
			gf2n::B[w] ^= (h << l);
			gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
		}
		else
			gf2n::B[w] ^= h;

		w = gf2n::exp3 / (CHAR_BIT*sizeof(gf2n_word));
		l = gf2n::exp3 % (CHAR_BIT*sizeof(gf2n_word));

		if (l != 0) {
			gf2n::B[w] ^= (h << l);
			gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
		}
		else
			gf2n::B[w] ^= h;

		shift_right(gf2n::B, CHAR_BIT*sizeof(gf2n_word), Blast);
	}

	h = gf2n::B[0] & (EXP2(xdegree)-1);
	gf2n::B[0] ^= h;
	w = gf2n::degree / (CHAR_BIT*sizeof(gf2n_word));
	l = gf2n::degree % (CHAR_BIT*sizeof(gf2n_word));
	Blast = w;
	if (l != 0) {
		gf2n::B[w] ^= (h << l);
		gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
		Blast++;
	}
	else gf2n::B[w] ^= h;

	w = gf2n::exp1 / (CHAR_BIT*sizeof(gf2n_word));
	l = gf2n::exp1 % (CHAR_BIT*sizeof(gf2n_word));

	if (l != 0) {
		gf2n::B[w] ^= (h << l);
		gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
	}
	else gf2n::B[w] ^= h;

	w = gf2n::exp2 / (CHAR_BIT*sizeof(gf2n_word));
	l = gf2n::exp2 % (CHAR_BIT*sizeof(gf2n_word));

	if (l != 0) {
		gf2n::B[w] ^= (h << l);
		gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
	}
	else gf2n::B[w] ^= h;

	w = gf2n::exp3 / (CHAR_BIT*sizeof(gf2n_word));
	l = gf2n::exp3 % (CHAR_BIT*sizeof(gf2n_word));

	if (l != 0) {
		gf2n::B[w] ^= (h << l);
		gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
	}
	else gf2n::B[w] ^= h;

	shift_right(gf2n::B, xdegree, Blast);

	for (i = 0, bp = gf2n::B, cp = res; i <= Blast; i++, bp++, cp++)  // copy B into res
		*cp = *bp;

	if (Blast < anzBI-1)
		for (i = Blast+1; i < anzBI; i++, cp++)
			*cp = static_cast<gf2n_word>(0);

	gf2n::B = fp;
}



//**********************************************************************

void general_invert(gf2n_word *res, gf2n_word *a)
{
	register gf2n_word h;
	register gf2n_word *bp, *cp, *fp, *gp, *ap;
	unsigned int i, j, w, l, s;
	unsigned int xdegree = 0;
	unsigned int anzBI = gf2n::anzBI;
	unsigned int Clast, Blast, Flast, Glast;

	partial_reduce2[gf2n::invsel](a); // make sure that a is totally reduced 

	bp = gf2n::B; cp = gf2n::C; fp = gf2n::F; gp = gf2n::G; ap = a;


	for (i = 0, ap = a; i < anzBI; i++, fp++, gp++, cp++, bp++, ap++) {
                // initialize B=1, C = 0, F = a, G = modulus, res = 0
		*fp = *ap;
		*gp = *bp = *cp = static_cast<gf2n_word>(0);
	}

	gf2n::B[0] = static_cast<gf2n_word>(1);

	for (i = anzBI, bp = gf2n::B+anzBI, cp = gf2n::C+anzBI, fp = gf2n::F+anzBI, gp = gf2n::G+anzBI;
	     i < 2*anzBI; i++, bp++, cp++, fp++, gp++) {
		*fp = *gp = *cp = *bp = static_cast<gf2n_word>(0);
	}

	gf2n::G[0] = static_cast<gf2n_word>(1);
	i = gf2n::degree;
	gf2n::G[i/(CHAR_BIT*sizeof(gf2n_word))] ^= EXP2(i % (CHAR_BIT*sizeof(gf2n_word)));

	for (j = 0; j < gf2n::anz_exponents; j++) {
		i = gf2n::exponents[j];
		gf2n::G[i/(CHAR_BIT*sizeof(gf2n_word))] ^= EXP2(i % (CHAR_BIT*sizeof(gf2n_word)));
	}

	Clast = Blast = 0; // initialize last non zero indices
	Flast = Glast = anzBI-1;

	fp = gf2n::F+Flast;
	while (*fp == static_cast<gf2n_word>(0) && Flast > 0) {
		Flast --; fp--;
	}

	//*************************************************************
	// start with gcd-computation, use the invariant:
	//
	// B * F + ?? * m = F
	// C * F + ?? * m = G
	//
	//*************************************************************


	// make F odd, if necessary
	h = gf2n::F[0];
	fp = gf2n::F;

	if (!(h & 1)) {
		while (h == static_cast<gf2n_word>(0)) {
			fp++;
			xdegree += CHAR_BIT*sizeof(gf2n_word);
			h = *fp;
		}

		while (!(h & 1)) {
			xdegree++;
			h >>= 1;
		}
		shift_right(gf2n::F, xdegree, Flast);

		if (Flast == 0 && gf2n::F[0] == 1)
			goto final_step;
	}

	while (true) {
		// degree(G) > degree(F)
		do {
			j = static_cast<unsigned int>(max(Clast, Blast));
			for (i = 0, cp = gf2n::C, bp = gf2n::B; i <= j; i++, cp++, bp++)   // C = C+B
				(*cp) ^= (*bp);

			Clast = j;
			cp--;
			while (*cp == static_cast<gf2n_word>(0) && Clast > 0) {
				Clast --;
				cp--;
			}

			h = gf2n::G[0] ^ gf2n::F[0];
			i = 0;

			while (h == static_cast<gf2n_word>(0)) {
				xdegree += CHAR_BIT*sizeof(gf2n_word);
				i += CHAR_BIT*sizeof(gf2n_word);
				h = gf2n::G[i/(CHAR_BIT*sizeof(gf2n_word))] ^ gf2n::F[i/(CHAR_BIT*sizeof(gf2n_word))];
			}

			while (!(h & 1)) {
				xdegree++;
				i++;
				h >>= 1;
			}
			add_shift_right(gf2n::G, gf2n::F, i, Glast, Flast);
			shift_left(gf2n::B, i, Blast);

			if (Glast == 0 && gf2n::G[0] == static_cast<gf2n_word>(1))          // G == 1  ??
				goto final_step;
		} while (Glast > Flast || ((gf2n::G[Glast] >= gf2n::F[Flast]) && (Glast == Flast)));

		do {
			// deg F > deg G 
			j = static_cast<unsigned int>(max(Clast, Blast));
			for (i = 0, bp = gf2n::B, cp = gf2n::C; i <= j; i++, bp++, cp++)   // B = B + C
				(*bp) ^= (*cp);

			Blast = j;
			bp--;
			while (*bp == 0 && Blast > 0) {
				Blast --;
				bp--;
			}

			i = 0; // make F odd
			h = gf2n::F[0] ^ gf2n::G[0];

			while (h == static_cast<gf2n_word>(0)) {
				xdegree += CHAR_BIT*sizeof(gf2n_word);
				i += CHAR_BIT*sizeof(gf2n_word);
				h = gf2n::F[i/(CHAR_BIT*sizeof(gf2n_word))] ^ gf2n::G[i/(CHAR_BIT*sizeof(gf2n_word))];
			}

			while (!(h & 1)) {
				xdegree++;
				i ++;
				h >>= 1;
			}
			add_shift_right(gf2n::F, gf2n::G, i, Flast, Glast);
			shift_left(gf2n::C, i, Clast);

			if (Flast == 0 && gf2n::F[0] == static_cast<gf2n_word>(1))   // F == 1
				goto final_step;
		} while (Flast > Glast || ((gf2n::F[Flast] >= gf2n::G[Glast]) && (Flast == Glast)));
	}


 final_step:       // check whether F == 1 mod modulus or G == 1 mod modulus

	fp = gf2n::B;

	if (gf2n::G[0] == static_cast<gf2n_word>(1) && Glast == 0) {
		gf2n::B = gf2n::C;
		Blast = Clast;
	}

	// now B == x^(xdegree) mod modulus and B has last index Blast

	while (xdegree >= CHAR_BIT*sizeof(gf2n_word)) {
		xdegree -= CHAR_BIT*sizeof(gf2n_word);
		h = gf2n::B[0];

		while (h != 0) {
			gf2n::B[0] = 0;
			w = gf2n::degree / (CHAR_BIT*sizeof(gf2n_word));
			l = gf2n::degree % (CHAR_BIT*sizeof(gf2n_word));
			Blast = w;

			if (l != 0) {
				gf2n::B[w] ^= (h << l);
				gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
				Blast++;
			}
			else
				gf2n::B[w] ^= h;

			for (j = 0; j < gf2n::anz_exponents; j++) {
				s = gf2n::exponents[j];
				w = s / (CHAR_BIT*sizeof(gf2n_word));
				l = s % (CHAR_BIT*sizeof(gf2n_word));

				if (l != 0) {
					gf2n::B[w] ^= (h << l);
					gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
				}
				else
					gf2n::B[w] ^= h;
			}
			h = gf2n::B[0];
			while (gf2n::B[Blast] == 0 && Blast > 0)
				Blast --;
		}
		shift_right(gf2n::B, CHAR_BIT*sizeof(gf2n_word), Blast);
	}

	h = gf2n::B[0] & (EXP2(xdegree)-1);

	while (h != 0)
        {
		gf2n::B[0] ^= h;
		w = gf2n::degree / (CHAR_BIT*sizeof(gf2n_word));
		l = gf2n::degree % (CHAR_BIT*sizeof(gf2n_word));
		Blast = w;
		if (l != 0) {
			gf2n::B[w] ^= (h << l);
			gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
			Blast++;
		}
		else gf2n::B[w] ^= h;


		for (j = 0; j < gf2n::anz_exponents; j++) {
			s = gf2n::exponents[j];
			w = s / (CHAR_BIT*sizeof(gf2n_word));
			l = s % (CHAR_BIT*sizeof(gf2n_word));

			if (l != 0) {
				gf2n::B[w] ^= (h << l);
				gf2n::B[w+1] ^= (h >> (CHAR_BIT*sizeof(gf2n_word)-l));
			}
			else gf2n::B[w] ^= h;
		}
		h = gf2n::B[0] & (EXP2(xdegree)-1);
		while (gf2n::B[Blast] == 0 && Blast > 0)
			Blast --;
        }
	shift_right(gf2n::B, xdegree, Blast);

	for (i = 0, bp = gf2n::B, cp = res; i <= Blast; i++, bp++, cp++)  // copy B into res
		*cp = *bp;

	if (Blast < anzBI-1)
		for (i = Blast+1; i < anzBI; i++, cp++)
			*cp = static_cast<gf2n_word>(0);

	gf2n::B = fp;
}



//******* invert - selector *********************************

void (*uinvert[]) (gf2n_word*, gf2n_word*) =
{
	tri_invert, pent_invert, general_invert
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
