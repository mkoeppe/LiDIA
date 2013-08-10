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
//                Volker Mueller (VM), Thomas Pfahler (TPf)
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



void tri_partial_reduce2(udigit *);
void pent_partial_reduce2(udigit *);
extern void (*partial_reduce2[]) (udigit *);



static
inline udigit max(udigit a, udigit b)
{
	if (a >= b)
		return a;
	else return b;
}



//************************************************************************
// shift F k bits to the left, adjust Flast
//************************************************************************

static
inline void shift_left(udigit *F, unsigned int k, unsigned int & Flast)
{
	register int i, j, s;
	register udigit *fp1, *fp2;

	s = k % (BITS_PER_LONG);
	j = k / (BITS_PER_LONG);

	if (s == 0) {
		for (i = Flast, fp1 = F+Flast+j, fp2 = F+Flast; i >= 0; i--, fp1--, fp2--)
			*fp1 = *fp2;
		Flast += j;
	}
	else {

		for (i = Flast, fp1 = F+Flast+j+1, fp2 = F+Flast; i >= 0; i--, fp2--) {
			(*fp1) |= ((*fp2) >> (BITS_PER_LONG-s));
			fp1--;
			(*fp1) = ((*fp2) << s);
		}
		Flast += j+1;
	}

	i = 0; fp1 = F;
	while (i < j) {
		*fp1 = static_cast<udigit>(0);
		i++; fp1++;
	}

	fp1 = F + Flast;
	while ((*fp1) == static_cast<udigit>(0) && Flast > 0) {
		Flast --;
		fp1--;
	}
}



//************************************************************************
// shift F k bits to the right, adjust Flast
//************************************************************************

static
inline void shift_right(udigit *F, unsigned int k, unsigned int & Flast)
{
	register unsigned int j, i, s;
	register udigit *fp1, *fp2;

	s = k % (BITS_PER_LONG);
	j = k / (BITS_PER_LONG);

	if (s == 0) {
		for (i = j, fp1 = F, fp2 = F+j; i <= Flast; i++, fp1++, fp2++)
			*fp1 = *fp2;
	}
	else {
		for (i = j, fp1 = F, fp2 = F+j; i < Flast; i++, fp1++, fp2++)
			*fp1 = ((*fp2) >> s) | ((*(fp2+1)) << (BITS_PER_LONG-s));
		*fp1 = ((*fp2) >> s);
	}
	i = Flast-j+1;
	fp1 = F+i;
	while (i <= Flast) {
		*fp1 = static_cast<udigit>(0);
		fp1++;
		i++;
	}

	Flast -= j;
	fp1 = F+Flast;
	while (*fp1 == static_cast<udigit>(0) && Flast > 0) {
		Flast --; fp1--;
	}
}



//************************************************************************
// compute F = (F + G) / x^k, adjust Flast
//************************************************************************

static
inline void add_shift_right(udigit *F, udigit*G, unsigned int k,
			    unsigned int & Flast, unsigned int Glast)

{
	register unsigned int j, i, s, l;
	register udigit h;
	register udigit *fp1, *fp2, *gp;

	s = k % (BITS_PER_LONG);
	j = k / (BITS_PER_LONG);

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
			*fp1 |= (h << (BITS_PER_LONG-s));
		}
	}
	i = l-j+1;
	fp1 = F+i;
	while (i <= Flast) {
		*fp1 = static_cast<udigit>(0);
		i++; fp1++;
	}


	Flast = l - j;
	fp1 = F+Flast;
	while (*fp1 == static_cast<udigit>(0) && Flast > 0) {
		Flast --;
		fp1--;
	}
}



//*********************************************************************
// invert function, res and a are assumed to have length anzBI
//*********************************************************************

void info_gf2n::tri_invert(udigit *res, udigit *a) const
{
	register udigit h;
	register udigit *bp, *cp, *fp, *gp, *ap;
	unsigned int i, j, w, l;
	unsigned int xdegree = 0;
	unsigned int Clast, Blast, Flast, Glast;

	tri_partial_reduce2(a); // make sure that a is totally reduced 

	bp = info_gf2n::B; cp = info_gf2n::C; fp = info_gf2n::F; gp = info_gf2n::G; ap = a;

	for (i = 0; i < anzBI; i++, fp++, gp++, cp++, bp++, ap++) {
                // initialize B=1, C = 0, F = a, G = modulus, res = 0
		*fp = *ap;
		*gp = *cp = *bp = static_cast<udigit>(0);
	}


	info_gf2n::B[0] = static_cast<udigit>(1);

	for (i = anzBI, bp = info_gf2n::B+anzBI, cp = info_gf2n::C+anzBI, fp = info_gf2n::F+anzBI,
		     gp = info_gf2n::G+anzBI; i < 2*anzBI; i++, bp++, cp++, fp++, gp++) {
		*fp = *gp = *cp = *bp = static_cast<udigit>(0);
	}


	info_gf2n::G[0] = static_cast<udigit>(1);
	i = degree;
	info_gf2n::G[i/(BITS_PER_LONG)] ^= udigit (1 << (i % (BITS_PER_LONG)));
	info_gf2n::G[exp1/(BITS_PER_LONG)] ^= udigit (1 << (exp1 % (BITS_PER_LONG)));


	Clast = Blast = 0; // initialize last non zero indices
	Flast = Glast = anzBI-1;

	fp = info_gf2n::F+Flast;
	while (*fp == static_cast<udigit>(0) && Flast > 0) {
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
	h = info_gf2n::F[0];
	fp = info_gf2n::F;

	if (!(h & 1)) {
		while (h == static_cast<udigit>(0)) {
			fp++;
			xdegree += BITS_PER_LONG;
			h = *fp;
		}

		while (!(h & 1)) {
			xdegree++;
			h >>= 1;
		}
		shift_right(info_gf2n::F, xdegree, Flast);

		if (Flast == 0 && info_gf2n::F[0] == 1)
			goto final_step;
	}

	while (true) {
		// degree(G) > degree(F)
		do {
			j = static_cast<unsigned int>(max(Clast, Blast));
			for (i = 0, cp = info_gf2n::C, bp = info_gf2n::B; i <= j; i++, cp++, bp++)   // C = C+B
				(*cp) ^= (*bp);

			Clast = j;
			cp--;
			while (*cp == static_cast<udigit>(0) && Clast > 0) {
				Clast --;
				cp--;
			}

			h = info_gf2n::G[0] ^ info_gf2n::F[0];
			i = 0;

			while (h == static_cast<udigit>(0)) {
				xdegree += BITS_PER_LONG;
				i += BITS_PER_LONG;
				h = info_gf2n::G[i/(BITS_PER_LONG)] ^ info_gf2n::F[i/(BITS_PER_LONG)];
			}

			while (!(h & 1)) {
				xdegree++;
				i++;
				h >>= 1;
			}
			add_shift_right(info_gf2n::G, info_gf2n::F, i, Glast, Flast);
			shift_left(info_gf2n::B, i, Blast);

			if (Glast == 0 && info_gf2n::G[0] == static_cast<udigit>(1))          // G == 1  ??
				goto final_step;
		}
		while (Glast > Flast || ((info_gf2n::G[Glast] >= info_gf2n::F[Flast]) && (Glast == Flast)));

		do                                     // deg F > deg G 
		{
			j = static_cast<unsigned int>(max(Clast, Blast));
			for (i = 0, bp = info_gf2n::B, cp = info_gf2n::C; i <= j; i++, bp++, cp++)   // B = B + C
				(*bp) ^= (*cp);

			Blast = j;
			bp--;
			while (*bp == 0 && Blast > 0) {
				Blast --;
				bp--;
			}

			i = 0; // make F odd
			h = info_gf2n::F[0] ^ info_gf2n::G[0];

			while (h == static_cast<udigit>(0)) {
				xdegree += BITS_PER_LONG;
				i += BITS_PER_LONG;
				h = info_gf2n::F[i/(BITS_PER_LONG)] ^ info_gf2n::G[i/(BITS_PER_LONG)];
			}

			while (!(h & 1)) {
				xdegree++;
				i ++;
				h >>= 1;
			}
			add_shift_right(info_gf2n::F, info_gf2n::G, i, Flast, Glast);
			shift_left(info_gf2n::C, i, Clast);

			if (Flast == 0 && info_gf2n::F[0] == static_cast<udigit>(1))   // F == 1
				goto final_step;
		}
		while (Flast > Glast || ((info_gf2n::F[Flast] >= info_gf2n::G[Glast]) && (Flast == Glast)));
	}


 final_step:       // check whether F == 1 mod modulus or G == 1 mod modulus

	fp = info_gf2n::B;

	if (info_gf2n::G[0] == static_cast<udigit>(1) && Glast == 0) {
		info_gf2n::B = info_gf2n::C;
		Blast = Clast;
	}

	// now B == x^(xdegree) mod modulus and B has last index Blast

	while (xdegree >= BITS_PER_LONG) {
		xdegree -= BITS_PER_LONG;
		h = info_gf2n::B[0];

		info_gf2n::B[0] = 0;
		w = degree / (BITS_PER_LONG);
		l = degree % (BITS_PER_LONG);
		Blast = w;

		if (l != 0) {
			info_gf2n::B[w] ^= (h << l);
			info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
			Blast++;
		}
		else
			info_gf2n::B[w] ^= h;

		w = exp1 / (BITS_PER_LONG);
		l = exp1 % (BITS_PER_LONG);

		if (l != 0) {
			info_gf2n::B[w] ^= (h << l);
			info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
		}
		else
			info_gf2n::B[w] ^= h;

		shift_right(info_gf2n::B, BITS_PER_LONG, Blast);
	}

	h = info_gf2n::B[0] & ((1 << xdegree)-1);
	info_gf2n::B[0] ^= h;
	w = degree / (BITS_PER_LONG);
	l = degree % (BITS_PER_LONG);
	Blast = w;
	if (l != 0) {
		info_gf2n::B[w] ^= (h << l);
		info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
		Blast++;
	}
	else info_gf2n::B[w] ^= h;

	w = exp1 / (BITS_PER_LONG);
	l = exp1 % (BITS_PER_LONG);

	if (l != 0) {
		info_gf2n::B[w] ^= (h << l);
		info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
	}
	else info_gf2n::B[w] ^= h;

	shift_right(info_gf2n::B, xdegree, Blast);

	for (i = 0, bp = info_gf2n::B, cp = res; i <= Blast; i++, bp++, cp++)  // copy B into res
		*cp = *bp;

	if (Blast < anzBI-1)
		for (i = Blast+1; i < anzBI; i++, cp++)
			*cp = static_cast<udigit>(0);

	info_gf2n::B = fp;
}



//******************************************************************

void info_gf2n::pent_invert(udigit *res, udigit *a) const
{
	register udigit h;
	register udigit *bp, *cp, *fp, *gp, *ap;
	unsigned int i, j, w, l;
	unsigned int xdegree = 0;
	unsigned int anzBI = info_gf2n::anzBI;
	unsigned int Clast, Blast, Flast, Glast;

	pent_partial_reduce2(a); // make sure that a is totally reduced 

	bp = info_gf2n::B; cp = info_gf2n::C; fp = info_gf2n::F; gp = info_gf2n::G; ap = a;

	for (i = 0, ap = a; i < anzBI; i++, fp++, gp++, cp++, bp++, ap++) {
                // initialize B=1, C = 0, F = a, G = modulus, res = 0
		*fp = *ap;
		*gp = *bp = *cp = static_cast<udigit>(0);
	}

	info_gf2n::B[0] = static_cast<udigit>(1);

	for (i = anzBI, bp = info_gf2n::B+anzBI, cp = info_gf2n::C+anzBI, fp = info_gf2n::F+anzBI, gp = info_gf2n::G+anzBI;
	     i < 2*anzBI; i++, bp++, cp++, fp++, gp++) {
		*fp = *gp = *cp = *bp = static_cast<udigit>(0);
	}

	info_gf2n::G[0] = static_cast<udigit>(1);
	i = degree;
	info_gf2n::G[i/(BITS_PER_LONG)] ^= udigit (1 << (i % (BITS_PER_LONG)));
	info_gf2n::G[exp1/(BITS_PER_LONG)] ^= udigit (1 << (exp1 % (BITS_PER_LONG)));
	info_gf2n::G[exp2/(BITS_PER_LONG)] ^= udigit (1 << (exp2 % (BITS_PER_LONG)));
	info_gf2n::G[exp3/(BITS_PER_LONG)] ^= udigit (1 << (exp3 % (BITS_PER_LONG)));

	Clast = Blast = 0; // initialize last non zero indices
	Flast = Glast = anzBI-1;

	fp = info_gf2n::F+Flast;
	while (*fp == static_cast<udigit>(0) && Flast > 0) {
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
	h = info_gf2n::F[0];
	fp = info_gf2n::F;

	if (!(h & 1)) {
		while (h == static_cast<udigit>(0)) {
			fp++;
			xdegree += BITS_PER_LONG;
			h = *fp;
		}

		while (!(h & 1)) {
			xdegree++;
			h >>= 1;
		}
		shift_right(info_gf2n::F, xdegree, Flast);

		if (Flast == 0 && info_gf2n::F[0] == 1)
			goto final_step;
	}



	while (true) {
		/* degree(G) > degree(F)*/
		do {
			j = static_cast<unsigned int>(max(Clast, Blast));
			for (i = 0, cp = info_gf2n::C, bp = info_gf2n::B; i <= j; i++, cp++, bp++)   // C = C+B
				(*cp) ^= (*bp);

			Clast = j;
			cp--;
			while (*cp == static_cast<udigit>(0) && Clast > 0) {
				Clast --;
				cp--;
			}

			h = info_gf2n::G[0] ^ info_gf2n::F[0];
			i = 0;

			while (h == static_cast<udigit>(0)) {
				xdegree += BITS_PER_LONG;
				i += BITS_PER_LONG;
				h = info_gf2n::G[i/(BITS_PER_LONG)] ^ info_gf2n::F[i/(BITS_PER_LONG)];
			}

			while (!(h & 1)) {
				xdegree++;
				i++;
				h >>= 1;
			}
			add_shift_right(info_gf2n::G, info_gf2n::F, i, Glast, Flast);
			shift_left(info_gf2n::B, i, Blast);

			if (Glast == 0 && info_gf2n::G[0] == static_cast<udigit>(1))          // G == 1  ??
				goto final_step;
		} while (Glast > Flast || ((info_gf2n::G[Glast] >= info_gf2n::F[Flast]) && (Glast == Flast)));

		do {
			// deg F > deg G 
			j = static_cast<unsigned int>(max(Clast, Blast));
			for (i = 0, bp = info_gf2n::B, cp = info_gf2n::C; i <= j; i++, bp++, cp++)   // B = B + C
				(*bp) ^= (*cp);

			Blast = j;
			bp--;
			while (*bp == 0 && Blast > 0) {
				Blast --;
				bp--;
			}

			i = 0; // make F odd
			h = info_gf2n::F[0] ^ info_gf2n::G[0];

			while (h == static_cast<udigit>(0)) {
				xdegree += BITS_PER_LONG;
				i += BITS_PER_LONG;
				h = info_gf2n::F[i/(BITS_PER_LONG)] ^ info_gf2n::G[i/(BITS_PER_LONG)];
			}

			while (!(h & 1)) {
				xdegree++;
				i ++;
				h >>= 1;
			}
			add_shift_right(info_gf2n::F, info_gf2n::G, i, Flast, Glast);
			shift_left(info_gf2n::C, i, Clast);

			if (Flast == 0 && info_gf2n::F[0] == static_cast<udigit>(1))   // F == 1
				goto final_step;
		} while (Flast > Glast || ((info_gf2n::F[Flast] >= info_gf2n::G[Glast]) && (Flast == Glast)));
	}


 final_step:       // check whether F == 1 mod modulus or G == 1 mod modulus

	fp = info_gf2n::B;

	if (info_gf2n::G[0] == static_cast<udigit>(1) && Glast == 0) {
		info_gf2n::B = info_gf2n::C;
		Blast = Clast;
	}
	// now B == x^(xdegree) mod modulus and B has last index Blast

	while (xdegree >= BITS_PER_LONG) {
		xdegree -= BITS_PER_LONG;
		h = info_gf2n::B[0];

		info_gf2n::B[0] = 0;
		w = degree / (BITS_PER_LONG);
		l = degree % (BITS_PER_LONG);
		Blast = w;

		if (l != 0) {
			info_gf2n::B[w] ^= (h << l);
			info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
			Blast++;
		}
		else
			info_gf2n::B[w] ^= h;

		w = exp1 / (BITS_PER_LONG);
		l = exp1 % (BITS_PER_LONG);

		if (l != 0) {
			info_gf2n::B[w] ^= (h << l);
			info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
		}
		else
			info_gf2n::B[w] ^= h;

		w = exp2 / (BITS_PER_LONG);
		l = exp2 % (BITS_PER_LONG);

		if (l != 0) {
			info_gf2n::B[w] ^= (h << l);
			info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
		}
		else
			info_gf2n::B[w] ^= h;

		w = exp3 / (BITS_PER_LONG);
		l = exp3 % (BITS_PER_LONG);

		if (l != 0) {
			info_gf2n::B[w] ^= (h << l);
			info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
		}
		else
			info_gf2n::B[w] ^= h;

		shift_right(info_gf2n::B, BITS_PER_LONG, Blast);
	}

	h = info_gf2n::B[0] & ((1 << xdegree)-1);
	info_gf2n::B[0] ^= h;
	w = degree / (BITS_PER_LONG);
	l = degree % (BITS_PER_LONG);
	Blast = w;
	if (l != 0) {
		info_gf2n::B[w] ^= (h << l);
		info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
		Blast++;
	}
	else info_gf2n::B[w] ^= h;

	w = exp1 / (BITS_PER_LONG);
	l = exp1 % (BITS_PER_LONG);

	if (l != 0) {
		info_gf2n::B[w] ^= (h << l);
		info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
	}
	else info_gf2n::B[w] ^= h;

	w = exp2 / (BITS_PER_LONG);
	l = exp2 % (BITS_PER_LONG);

	if (l != 0) {
		info_gf2n::B[w] ^= (h << l);
		info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
	}
	else info_gf2n::B[w] ^= h;

	w = exp3 / (BITS_PER_LONG);
	l = exp3 % (BITS_PER_LONG);

	if (l != 0) {
		info_gf2n::B[w] ^= (h << l);
		info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
	}
	else info_gf2n::B[w] ^= h;

	shift_right(info_gf2n::B, xdegree, Blast);

	for (i = 0, bp = info_gf2n::B, cp = res; i <= Blast; i++, bp++, cp++)  // copy B into res
		*cp = *bp;

	if (Blast < anzBI-1)
		for (i = Blast+1; i < anzBI; i++, cp++)
			*cp = static_cast<udigit>(0);

	info_gf2n::B = fp;
}



//**********************************************************************

void info_gf2n::general_invert(udigit *res, udigit *a) const
{
	register udigit h;
	register udigit *bp, *cp, *fp, *gp, *ap;
	unsigned int i, j, w, l, s;
	unsigned int xdegree = 0;
	unsigned int anzBI = info_gf2n::anzBI;
	unsigned int Clast, Blast, Flast, Glast;

	general_partial_reduce2(a); // make sure that a is totally reduced 

	bp = info_gf2n::B; cp = info_gf2n::C; fp = info_gf2n::F; gp = info_gf2n::G; ap = a;


	for (i = 0, ap = a; i < anzBI; i++, fp++, gp++, cp++, bp++, ap++) {
                // initialize B=1, C = 0, F = a, G = modulus, res = 0
		*fp = *ap;
		*gp = *bp = *cp = static_cast<udigit>(0);
	}

	info_gf2n::B[0] = static_cast<udigit>(1);

	for (i = anzBI, bp = info_gf2n::B+anzBI, cp = info_gf2n::C+anzBI, fp = info_gf2n::F+anzBI, gp = info_gf2n::G+anzBI;
	     i < 2*anzBI; i++, bp++, cp++, fp++, gp++) {
		*fp = *gp = *cp = *bp = static_cast<udigit>(0);
	}

	info_gf2n::G[0] = static_cast<udigit>(1);
	i = degree;
	info_gf2n::G[i/(BITS_PER_LONG)] ^= udigit (1 << (i % (BITS_PER_LONG)));

	for (j = 0; j < anz_exponents; j++) {
		i = exponents[j];
		info_gf2n::G[i/(BITS_PER_LONG)] ^= udigit (1 << (i % (BITS_PER_LONG)));
	}

	Clast = Blast = 0; // initialize last non zero indices
	Flast = Glast = anzBI-1;

	fp = info_gf2n::F+Flast;
	while (*fp == static_cast<udigit>(0) && Flast > 0) {
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
	h = info_gf2n::F[0];
	fp = info_gf2n::F;

	if (!(h & 1)) {
		while (h == static_cast<udigit>(0)) {
			fp++;
			xdegree += BITS_PER_LONG;
			h = *fp;
		}

		while (!(h & 1)) {
			xdegree++;
			h >>= 1;
		}
		shift_right(info_gf2n::F, xdegree, Flast);

		if (Flast == 0 && info_gf2n::F[0] == 1)
			goto final_step;
	}

	while (true) {
		// degree(G) > degree(F)
		do {
			j = static_cast<unsigned int>(max(Clast, Blast));
			for (i = 0, cp = info_gf2n::C, bp = info_gf2n::B; i <= j; i++, cp++, bp++)   // C = C+B
				(*cp) ^= (*bp);

			Clast = j;
			cp--;
			while (*cp == static_cast<udigit>(0) && Clast > 0) {
				Clast --;
				cp--;
			}

			h = info_gf2n::G[0] ^ info_gf2n::F[0];
			i = 0;

			while (h == static_cast<udigit>(0)) {
				xdegree += BITS_PER_LONG;
				i += BITS_PER_LONG;
				h = info_gf2n::G[i/(BITS_PER_LONG)] ^ info_gf2n::F[i/(BITS_PER_LONG)];
			}

			while (!(h & 1)) {
				xdegree++;
				i++;
				h >>= 1;
			}
			add_shift_right(info_gf2n::G, info_gf2n::F, i, Glast, Flast);
			shift_left(info_gf2n::B, i, Blast);

			if (Glast == 0 && info_gf2n::G[0] == static_cast<udigit>(1))          // G == 1  ??
				goto final_step;
		} while (Glast > Flast || ((info_gf2n::G[Glast] >= info_gf2n::F[Flast]) && (Glast == Flast)));

		do {
			// deg F > deg G 
			j = static_cast<unsigned int>(max(Clast, Blast));
			for (i = 0, bp = info_gf2n::B, cp = info_gf2n::C; i <= j; i++, bp++, cp++)   // B = B + C
				(*bp) ^= (*cp);

			Blast = j;
			bp--;
			while (*bp == 0 && Blast > 0) {
				Blast --;
				bp--;
			}

			i = 0; // make F odd
			h = info_gf2n::F[0] ^ info_gf2n::G[0];

			while (h == static_cast<udigit>(0)) {
				xdegree += BITS_PER_LONG;
				i += BITS_PER_LONG;
				h = info_gf2n::F[i/(BITS_PER_LONG)] ^ info_gf2n::G[i/(BITS_PER_LONG)];
			}

			while (!(h & 1)) {
				xdegree++;
				i ++;
				h >>= 1;
			}
			add_shift_right(info_gf2n::F, info_gf2n::G, i, Flast, Glast);
			shift_left(info_gf2n::C, i, Clast);

			if (Flast == 0 && info_gf2n::F[0] == static_cast<udigit>(1))   // F == 1
				goto final_step;
		} while (Flast > Glast || ((info_gf2n::F[Flast] >= info_gf2n::G[Glast]) && (Flast == Glast)));
	}


 final_step:       // check whether F == 1 mod modulus or G == 1 mod modulus

	fp = info_gf2n::B;

	if (info_gf2n::G[0] == static_cast<udigit>(1) && Glast == 0) {
		info_gf2n::B = info_gf2n::C;
		Blast = Clast;
	}

	// now B == x^(xdegree) mod modulus and B has last index Blast

	while (xdegree >= BITS_PER_LONG) {
		xdegree -= BITS_PER_LONG;
		h = info_gf2n::B[0];

		while (h != 0) {
			info_gf2n::B[0] = 0;
			w = degree / (BITS_PER_LONG);
			l = degree % (BITS_PER_LONG);
			Blast = w;

			if (l != 0) {
				info_gf2n::B[w] ^= (h << l);
				info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
				Blast++;
			}
			else
				info_gf2n::B[w] ^= h;

			for (j = 0; j < anz_exponents; j++) {
				s = exponents[j];
				w = s / (BITS_PER_LONG);
				l = s % (BITS_PER_LONG);

				if (l != 0) {
					info_gf2n::B[w] ^= (h << l);
					info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
				}
				else
					info_gf2n::B[w] ^= h;
			}
			h = info_gf2n::B[0];
			while (info_gf2n::B[Blast] == 0 && Blast > 0)
				Blast --;
		}
		shift_right(info_gf2n::B, BITS_PER_LONG, Blast);
	}

	h = info_gf2n::B[0] & ((1 << xdegree)-1);

	while (h != 0) {
		info_gf2n::B[0] ^= h;
		w = degree / (BITS_PER_LONG);
		l = degree % (BITS_PER_LONG);
		Blast = w;
		if (l != 0) {
			info_gf2n::B[w] ^= (h << l);
			info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
			Blast++;
		}
		else info_gf2n::B[w] ^= h;


		for (j = 0; j < anz_exponents; j++) {
			s = exponents[j];
			w = s / (BITS_PER_LONG);
			l = s % (BITS_PER_LONG);

			if (l != 0) {
				info_gf2n::B[w] ^= (h << l);
				info_gf2n::B[w+1] ^= (h >> (BITS_PER_LONG-l));
			}
			else info_gf2n::B[w] ^= h;
		}
		h = info_gf2n::B[0] & ((1 << xdegree)-1);
		while (info_gf2n::B[Blast] == 0 && Blast > 0)
			Blast --;
        }
	shift_right(info_gf2n::B, xdegree, Blast);

	for (i = 0, bp = info_gf2n::B, cp = res; i <= Blast; i++, bp++, cp++)  // copy B into res
		*cp = *bp;

	if (Blast < anzBI-1)
		for (i = Blast+1; i < anzBI; i++, cp++)
			*cp = static_cast<udigit>(0);

	info_gf2n::B = fp;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
