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
#include	"LiDIA/gf2n_polynomial.h"
#include        "LiDIA/arith.inl"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



gf2n * gf2n_polynomial::gf2n_p_tmp = NULL;
int gf2n_polynomial::gf2n_p_tmpsize = 0;
int gf2n_polynomial::gf2n_p_top = 0; // first free position

//--------------------------------------------------------------------
// increase the stack

void gf2n_polynomial::resize_stack_for_degree(int deg)
{
	int i = gf2n_polynomial::gf2n_p_top +
		(2*deg * comparator< int >::max((static_cast<int>(LiDIA::log2(static_cast<double>(deg))) + 1), 12));

	gf2n_polynomial::resize_stack(i);
}



void gf2n_polynomial::resize_stack(int s)
{
	if (s < 0)
		lidia_error_handler("gf2n_polynomial", "negative stack size");

	if (gf2n_polynomial::gf2n_p_tmpsize == 0) {
		if (s == 0) {
			gf2n_polynomial::gf2n_p_tmp = new gf2n[512];
			gf2n_polynomial::gf2n_p_tmpsize = 512;
		}
		else {
			gf2n_polynomial::gf2n_p_tmp = new gf2n[s];
			gf2n_polynomial::gf2n_p_tmpsize = static_cast<int>(s);
		}
	}
	else {
		gf2n * n;

		if (s == 0)
			s = gf2n_polynomial::gf2n_p_tmpsize << 2;

		n = new gf2n[s];
		gf2n_polynomial::gf2n_p_tmpsize = s;

		if (s > gf2n_polynomial::gf2n_p_top)
			for (int i = 0; i < gf2n_polynomial::gf2n_p_top; i++)
				n[i].assign(gf2n_polynomial::gf2n_p_tmp[i]);
		else {
			lidia_error_handler("gf2n_polynomial", "gf2n_p_karatzuba: can"
					    " not resize stack");
			delete [] n;
			return;
		}

		delete[] gf2n_polynomial::gf2n_p_tmp;

		gf2n_polynomial::gf2n_p_tmp = n;
	}
}



void gf2n_p_mult(gf2n *r, const gf2n * a, const gf2n * b, int ad, int bd);

//=================================================================
// now we start with the multiplication functions of fixed size
// the array r MUST be big enough to store the result
// First we use a define for Karatzuba of fixed size ad, bd, later
// the compiler can optimize. This is unreadable but should be fast !!

#define KARA(ad, bd) { int old_top = gf2n_polynomial::gf2n_p_top; \
  int ad1, ad2, bd2 = 0, i; \
  gf2n* t1, *t2, *t3, *t4; \
  ad1 = ad / 2; if (ad & 1)  ad1 ++; \
  ad2 = ad - ad1; bd2 = bd - ad1; \
  gf2n_polynomial::gf2n_p_top += 6*ad1  + ad2 + bd2 + 10; \
  t1 = gf2n_polynomial::gf2n_p_tmp + old_top; \
  if (ad1 >= bd) {  gf2n_p_mult(t1, a, b, ad1, bd); \
  if (ad2 >= bd) gf2n_p_mult(t1 + ad1 + bd, a+ad1, b, ad2, bd); \
  else gf2n_p_mult(t1 + ad1 + bd, b, a+ad1, bd, ad2); \
  for (i = 0; i < ad1; i++) r[i].assign(t1[i]); \
  for (i = ad1; i < ad1 + bd - 1; i++) add(r[i], t1[i], t1[bd + i]); \
  for (i = ad1 + bd - 1; i < ad + bd - 1; i++) r[i].assign(t1[i + bd]); \
  gf2n_polynomial::gf2n_p_top = old_top; return; } \
  else {  t2 = t1 + 2*ad1; t3 = t2 + (ad2 + bd2); t4 = t3 + 2*ad1; \
  gf2n_p_mult(t1, a, b, ad1, ad1); gf2n_p_mult(t2, a+ad1, b + ad1, ad2, bd2); \
  for (i = 0; i < ad1; i++) if (i < ad2) add(t3[i], a[i], a[i+ad1]); \
  else t3[i].assign(a[i]); \
  for (i = 0; i < ad1; i++) if (i < bd2) add(t3[i+ad1], b[i], b[i+ad1]); \
  else t3[i+ad1].assign(b[i]); \
  gf2n_p_mult(t4, t3+ad1, t3, ad1, ad1); \
  for (i = 0; i < 2*ad1 - 1; i++)  r[i].assign(t1[i]); \
  for (i = 2*ad1 - 1; i < ad + bd - 1; i++) r[i].assign_zero(); \
  for (i = 0; i < ad2 + bd2 - 1; i++) { add(r[i+2*ad1], r[i+2*ad1], t2[i]); \
  add(r[i+ad1], r[i+ad1], t2[i] + t1[i] + t4[i]); } \
  for (i = ad2 + bd2 - 1; i < 2 * ad1 - 1; i++)  \
  add(r[i+ad1], r[i+ad1], t4[i] + t1[i]); \
  gf2n_polynomial::gf2n_p_top = old_top; \
  } }




//=================================================================
// now the preinstalled versions, ad >= bd > 1 always guaranteed

inline void gf2n_p_mult21(gf2n * r, const gf2n * a, const gf2n * b)
{
	multiply(r[0], a[0], b[0]);
	multiply(r[1], a[1], b[0]);
}



void gf2n_p_mult22(gf2n * r, const gf2n * a, const gf2n * b)
{
	multiply(r[0], a[0], b[0]);
	multiply(r[2], a[1], b[1]);
	add(r[1], (a[0] + a[1])*(b[1] + b[0]), r[0] + r[2]);
}



void gf2n_p_mult32(gf2n * r, const gf2n * a, const gf2n * b)
{
	gf2n_p_mult22(r, a, b);
	gf2n_p_mult21(gf2n_polynomial::gf2n_p_tmp+gf2n_polynomial::gf2n_p_top, b, a+2);
	add(r[2], r[2], gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top]);
	r[3].assign(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+1]);
}



void gf2n_p_mult33(gf2n * r, const gf2n * a, const gf2n * b)
{
	gf2n_p_mult22(gf2n_polynomial::gf2n_p_tmp+gf2n_polynomial::gf2n_p_top, a, b);
	multiply(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+3], a[2], b[2]);

	add(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+6], a[0], a[2]);
	gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+7].assign(a[1]);
	add(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+8], b[0], b[2]);
	gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+9].assign(b[1]);
	gf2n_p_mult22(r+2, gf2n_polynomial::gf2n_p_tmp + gf2n_polynomial::gf2n_p_top + 6,
		      gf2n_polynomial::gf2n_p_tmp + gf2n_polynomial::gf2n_p_top + 8);

	r[0] = gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top];
	r[1] = gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+1];
	r[2] += gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+2] + gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top]
		+ gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+3];
	r[3] += gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+1];
	r[4] += gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+3] + gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+2];
}



void gf2n_p_mult42(gf2n * r, const gf2n * a, const gf2n * b)
{
	gf2n_p_mult22(r, a, b);
	gf2n_p_mult22(gf2n_polynomial::gf2n_p_tmp+gf2n_polynomial::gf2n_p_top, a+2, b);
	r[2] += gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top];
	r[3] = gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+1];
	r[4] = gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+2];
}



void gf2n_p_mult43(gf2n * r, const gf2n * a, const gf2n * b)
{
	gf2n_p_mult22(gf2n_polynomial::gf2n_p_tmp+gf2n_polynomial::gf2n_p_top, a, b);
	gf2n_p_mult21(gf2n_polynomial::gf2n_p_tmp+gf2n_polynomial::gf2n_p_top+3, a+2, b+2);
	add(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+6], a[0], a[2]);
	add(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+7], a[1], a[3]);
	add(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+8], b[0], b[2]);
	gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+9].assign(b[1]);
	gf2n_p_mult22(r+2, gf2n_polynomial::gf2n_p_tmp+gf2n_polynomial::gf2n_p_top+6,
		      gf2n_polynomial::gf2n_p_tmp+gf2n_polynomial::gf2n_p_top + 8);
	r[0].assign(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top]);
	r[1].assign(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+1]);
	r[2] += gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+2] + gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top]
		+ gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+3];
	r[3] += gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+1] + gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+4];
	r[4] += gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+3] + gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+2];
	r[5].assign(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+4]);
}



void gf2n_p_mult44(gf2n * r, const gf2n * a, const gf2n * b)
{
	gf2n_p_mult22(gf2n_polynomial::gf2n_p_tmp+gf2n_polynomial::gf2n_p_top, a, b);
	gf2n_p_mult22(gf2n_polynomial::gf2n_p_tmp+gf2n_polynomial::gf2n_p_top+3, a+2, b+2);

	add(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+6], a[0], a[2]);
	add(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+7], a[1], a[3]);
	add(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+8], b[0], b[2]);
	add(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+9], b[1], b[3]);
	gf2n_p_mult22(r+2, gf2n_polynomial::gf2n_p_tmp + gf2n_polynomial::gf2n_p_top + 6,
		      gf2n_polynomial::gf2n_p_tmp + gf2n_polynomial::gf2n_p_top + 8);

	r[0].assign(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top]);
	r[1].assign(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+1]);
	r[2] += gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+2] + gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top]
		+ gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+3];
	r[3] += gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+1] + gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+4];
	r[4] += gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+2] + gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+5] +
		gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+3];
	r[5].assign(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+4]);
	r[6].assign(gf2n_polynomial::gf2n_p_tmp[gf2n_polynomial::gf2n_p_top+5]);
}



void gf2n_p_mult52(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(5, 2)

	void gf2n_p_mult53(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(5, 3)

	void gf2n_p_mult54(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(5, 4)

	void gf2n_p_mult55(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(5, 5)

	void gf2n_p_mult62(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(6, 2)

	void gf2n_p_mult63(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(6, 3)

	void gf2n_p_mult64(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(6, 4)

	void gf2n_p_mult65(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(6, 5)

	void gf2n_p_mult66(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(6, 6)

	void gf2n_p_mult72(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(7, 2)

	void gf2n_p_mult73(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(7, 3)

	void gf2n_p_mult74(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(7, 4)

	void gf2n_p_mult75(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(7, 5)

	void gf2n_p_mult76(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(7, 6)

	void gf2n_p_mult77(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(7, 7)

	void gf2n_p_mult82(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(8, 2)

	void gf2n_p_mult83(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(8, 3)

	void gf2n_p_mult84(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(8, 4)

	void gf2n_p_mult85(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(8, 5)

	void gf2n_p_mult86(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(8, 6)

	void gf2n_p_mult87(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(8, 7)

	void gf2n_p_mult88(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(8, 8)

	void gf2n_p_mult92(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(9, 2)

	void gf2n_p_mult93(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(9, 3)

	void gf2n_p_mult94(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(9, 4)

	void gf2n_p_mult95(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(9, 5)

	void gf2n_p_mult96(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(9, 6)

	void gf2n_p_mult97(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(9, 7)

	void gf2n_p_mult98(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(9, 8)

	void gf2n_p_mult99(gf2n * r, const gf2n * a, const gf2n * b)
	KARA(9, 9)


	void (*karamul[]) (gf2n *, const gf2n*, const gf2n*) =
{
	gf2n_p_mult22, gf2n_p_mult32, gf2n_p_mult33, gf2n_p_mult42,
	gf2n_p_mult43, gf2n_p_mult44, gf2n_p_mult52, gf2n_p_mult53,
	gf2n_p_mult54, gf2n_p_mult55, gf2n_p_mult62, gf2n_p_mult63,
	gf2n_p_mult64, gf2n_p_mult65, gf2n_p_mult66, gf2n_p_mult72,
	gf2n_p_mult73, gf2n_p_mult74, gf2n_p_mult75, gf2n_p_mult76,
	gf2n_p_mult77, gf2n_p_mult82, gf2n_p_mult83, gf2n_p_mult84,
	gf2n_p_mult85, gf2n_p_mult86, gf2n_p_mult87, gf2n_p_mult88,
	gf2n_p_mult92, gf2n_p_mult93, gf2n_p_mult94, gf2n_p_mult95,
	gf2n_p_mult96, gf2n_p_mult97, gf2n_p_mult98, gf2n_p_mult99
};


//================================================================
// the real function starts now. Small functions are selected via
// jump table.
// the result must be big enough, ad >= bd >= 1.
//


void gf2n_p_mult(gf2n *r, const gf2n * a, const gf2n * b, int ad, int bd)
{
	int i;
	if (bd == 1) {
		for (i = 0; i < ad; i++)
			multiply(r[i], a[i], b[0]);
		return;
	}

	if (ad <= 9) {
		karamul[(ad-2)*(ad-1)/2 + bd -2](r, a, b);
		return;
	}

	// now the recursive case

	int old_top = gf2n_polynomial::gf2n_p_top; // to reset the stack after computation
	int ad1, ad2, bd2 = 0;
	gf2n* t1, *t2, *t3, *t4;

	// split the polynomials in part of sizes [ad1, ad2], [ad1, bd2], resp.

	ad1 = ad / 2;
	if (ad & 1)
		ad1 ++;

	ad2 = ad - ad1; // ad1 >= ad2, bd1 = ad1 >= bd2
	bd2 = bd - ad1;

	gf2n_polynomial::gf2n_p_top += 6*ad1  + ad2 + bd2 + 10;
	t1 = gf2n_polynomial::gf2n_p_tmp + old_top;

	if (ad1 >= bd) {
		gf2n_p_mult(t1, a, b, ad1, bd);

		if (ad2 >= bd)
			gf2n_p_mult(t1 + ad1 + bd, a+ad1, b, ad2, bd);
		else
			gf2n_p_mult(t1 + ad1 + bd, b, a+ad1, bd, ad2);

		for (i = 0; i < ad1; i++)
			r[i].assign(t1[i]);

		for (i = ad1; i < ad1 + bd - 1; i++)
			add(r[i], t1[i], t1[bd + i]);

		for (i = ad1 + bd - 1; i < ad + bd - 1; i++)
			r[i].assign(t1[i + bd]);

		gf2n_polynomial::gf2n_p_top = old_top;
		return;
	}

	t2 = t1 + 2*ad1;
	t3 = t2 + (ad2 + bd2);
	t4 = t3 + 2*ad1;

	gf2n_p_mult(t1, a, b, ad1, ad1);
	gf2n_p_mult(t2, a+ad1, b + ad1, ad2, bd2);

	for (i = 0; i < ad1; i++)
		if (i < ad2)
			add(t3[i], a[i], a[i+ad1]);
		else
			t3[i].assign(a[i]);

	for (i = 0; i < ad1; i++)
		if (i < bd2)
			add(t3[i+ad1], b[i], b[i+ad1]);
		else
			t3[i+ad1].assign(b[i]);

	gf2n_p_mult(t4, t3+ad1, t3, ad1, ad1);

	for (i = 0; i < 2*ad1 - 1; i++)
		r[i].assign(t1[i]);
	for (i = 2*ad1 - 1; i < ad + bd - 1; i++)
		r[i].assign_zero();

	for (i = 0; i < ad2 + bd2 - 1; i++) {
		add(r[i+2*ad1], r[i+2*ad1], t2[i]);
		add(r[i+ad1], r[i+ad1], t2[i] + t1[i] + t4[i]);
	}

	for (i = ad2 + bd2 - 1; i < 2 * ad1 - 1; i++)
		add(r[i+ad1], r[i+ad1], t4[i] + t1[i]);

	gf2n_polynomial::gf2n_p_top = old_top;
}



//-----------------------------------------------------------------
// multiply function based on karatzuba


void multiply(gf2n_polynomial & c, const gf2n_polynomial & a,
	      const gf2n_polynomial & b)
{
	register int deg_a = a.deg, deg_b = b.deg;
	register int deg_ab = deg_a + deg_b;
	int s = c.size;
	gf2n * cpp;

	if (deg_a < 0 || deg_b < 0) {
		c.deg = -1;
		return;
	}

	if (&a == &b) {
		square(c, a);
		return;
	}

	if (deg_a == 0 || deg_b == 0) {
		plain_mul(c, a, b);
		return;
	}


	if (deg_ab >= s || c.coeff == a.coeff || c.coeff == b.coeff) {
		if (deg_ab >= static_cast<int>(gf2n_polynomial::default_size))
			c.size = deg_ab + 1;
		else
			c.size = gf2n_polynomial::default_size;

		cpp = new gf2n[c.size];
	}
	else
		cpp = c.coeff;

	int i = gf2n_polynomial::gf2n_p_top +
		(2*deg_ab * comparator< int >::max((static_cast<int>(LiDIA::log2(static_cast<double>(deg_ab))) + 1), 12));


	if (i > gf2n_polynomial::gf2n_p_tmpsize)
		gf2n_polynomial::resize_stack(i);

	if (deg_a >= deg_b)
		gf2n_p_mult(cpp, a.coeff, b.coeff, deg_a+1, deg_b+1);
	else
		gf2n_p_mult(cpp, b.coeff, a.coeff, deg_b+1, deg_a+1);

	gf2n_polynomial::gf2n_p_top = 0; // reset the stack for safety reasons !!

	if (cpp != c.coeff) {
		if (s > -1)
			delete[] c.coeff;
		c.coeff = cpp;
	}
	c.deg = deg_ab;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
