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



const unsigned int dim_mul = 256; // dimension for multiplication table 
const unsigned int dim_square = 65536; // dimension for squaring table       

static gf2n_bit16 tabmul[dim_mul][dim_mul]; // table for multiplication
static gf2n_word tabsquare[dim_square]; // table for squaring


//***************************************************************
// multiply two 8-bit numbers without carrys
//***************************************************************

gf2n_bit16 mul8bit (gf2n_bit16 a, gf2n_bit16 b)
{
	gf2n_bit16 h = 0;

	for (register int i = 1; i <= 8; i++) {
		if (!a)
			return (h);
		if (a&1) h ^= b;
		a >>= 1;
		b <<= 1;
	}
	return (h);

}



//**************************************************************
// multiply two 16 bit numbers without carrys
//**************************************************************

inline
gf2n_word mul16bit(gf2n_bit16 a, gf2n_bit16 b)
{
	gf2n_bit8 a0, a1, b0, b1;

	a0 = gf2n_bit8(a);
	a1 = gf2n_bit8(a >> 8);
	b0 = gf2n_bit8(b);
	b1 = gf2n_bit8(b >> 8);

	return(tabmul[a0][b0] ^
	       gf2n_word((tabmul[a0][b1] ^ tabmul[a1][b0]) << 8) ^
	       gf2n_word(tabmul[a1][b1] << 16));
}



//***************************************************************

void gf2n::generate_mul_table()
{
	register unsigned int i, j;

	for (i = 0; i < dim_mul; i++)
		tabmul[i][i] = mul8bit(i, i);

	for (i = 0; i < dim_mul-1; i++)
		for (j = i+1; j < dim_mul; j++) {
			tabmul[i][j] = mul8bit(i, j);
			tabmul[j][i] = tabmul[i][j];
		}
}



//**************************************************************

void gf2n::generate_square_table()
{
	register unsigned int i;

	for (i = 0; i < dim_square; i++)
		tabsquare[i] = mul16bit(i, i);
}




//**************************************************************
// multiplication of arrays of gf2n_words
//
//   mulj(...)->input is array [0, .., j-1] of gf2n_word
//   imulj(..)  is corresponding function as inline
//*************************************************************

#if GF2N_WORDSIZE == 64
#define gf2n_bit32(a)	((a) & 0xffffffff)

// a and b are considered 32 bit values. function returns a 64 bit value
inline
gf2n_word mul32bit (gf2n_word a, gf2n_word b)
{
	gf2n_word cm;
	gf2n_bit16 a0, a1, b0, b1;
	gf2n_word c0, c1;

	a0 = gf2n_bit16(a);
	b0 = gf2n_bit16(b);
	a1 = gf2n_bit16(a >> 16);
	b1 = gf2n_bit16(b >> 16);

	c0 = mul16bit(a0, b0);
	c1 = mul16bit(a1, b1);

	cm = mul16bit(a0^a1, b0^b1) ^ c0 ^ c1;

	return c0 ^ gf2n_bit32(cm << 16) ^ ((c1 ^ gf2n_bit32(cm >> 16)) << 32);
}



inline
void imul1 (gf2n_word *c, gf2n_word *a, gf2n_word *b)
{
	gf2n_word cm, a0, a1, b0, b1;

	a0 = gf2n_bit32(a[0]);
	b0 = gf2n_bit32(b[0]);
	a1 = gf2n_bit32(a[0] >> 32);
	b1 = gf2n_bit32(b[0] >> 32);

	c[0] = mul32bit(a0, b0);
	c[1] = mul32bit(a1, b1);

	cm = mul32bit(a0^a1, b0^b1) ^ c[0] ^ c[1];

	c[0] ^= (cm << 32);
	c[1] ^= (cm >> 32);
}
#else // we are running on a 32 bit architecture

inline
void imul1 (gf2n_word *c, gf2n_word *a, gf2n_word *b)
{
	gf2n_word cm;
	gf2n_bit16 a0, a1, b0, b1;

	a0 = gf2n_bit16(a[0]);
	b0 = gf2n_bit16(b[0]);
	a1 = gf2n_bit16(a[0] >> 16);
	b1 = gf2n_bit16(b[0] >> 16);

	c[0] = mul16bit(a0, b0);
	c[1] = mul16bit(a1, b1);

	cm = mul16bit(a0^a1, b0^b1) ^ c[0] ^ c[1];

	c[0] ^= (cm << 16);
	c[1] ^= (cm >> 16);
}
#endif



void mul1(gf2n_word *c, gf2n_word *a, gf2n_word *b)
{
	imul1(c, a, b);
}



inline
void imul2 (gf2n_word *c, gf2n_word *a, gf2n_word *b)
{
	gf2n_word hs0, hs1;
	gf2n_word hl2[2]; // MM

	hs0 = a[0] ^ a[1]; // a0 + a1
	hs1 = b[0] ^ b[1]; // b0 + b1

	mul1(c, a, b);
	mul1(c+2, a+1, b+1);
	mul1(hl2, & hs0, & hs1);

	hl2[0] = hl2[0] ^ c[0] ^ c[2];
	hl2[1] = hl2[1] ^ c[1] ^ c[3];

	c[1] ^= hl2[0];
	c[2] ^= hl2[1];
}



void mul2 (gf2n_word *c, gf2n_word *a, gf2n_word *b)
{
	imul2(c, a, b);
}



inline
void imul3 (gf2n_word *c, gf2n_word *a, gf2n_word *b)
{
	gf2n_word hs0[2], hs1[2];
	gf2n_word hl2[4];

	hs0[0] = a[0] ^ a[2]; // a0 + a1
	hs0[1] = a[1];
	hs1[0] = b[0] ^ b[2]; // b0 + b1
	hs1[1] = b[1];

	imul2(c, a, b); // a0 * b0
	imul1(c+4, a+2, b+2); // a1 * b1
	imul2(hl2, hs0, hs1); // () * ()

	hl2[0] = hl2[0] ^ c[0] ^ c[4]; // ()() + a0b0 + a1b1
	hl2[1] = hl2[1] ^ c[1] ^ c[5];
	hl2[2] = hl2[2] ^ c[2];
	hl2[3] = hl2[3] ^ c[3];

	c[2] ^= hl2[0];
	c[3] ^= hl2[1];
	c[4] ^= hl2[2];
	c[5] ^= hl2[3];
}



void mul3 (gf2n_word *c, gf2n_word *a, gf2n_word *b)
{
	imul3(c, a, b);
}



inline
void imul4 (gf2n_word *c, gf2n_word *a, gf2n_word *b)
{
	gf2n_word hs0[2], hs1[2];
	gf2n_word hl2[4];

	hs0[0] = a[0] ^ a[2]; // a0 + a1
	hs0[1] = a[1] ^ a[3];
	hs1[0] = b[0] ^ b[2]; // b0 + b1
	hs1[1] = b[1] ^ b[3];

	imul2(c, a, b); // a0 * b0
	imul2(c+4, a+2, b+2); // a1 * b1
	imul2(hl2, hs0, hs1); // () * ()

	hl2[0] = hl2[0] ^ c[0] ^ c[4]; // ()() + a0b0 + a1b1
	hl2[1] = hl2[1] ^ c[1] ^ c[5];
	hl2[2] = hl2[2] ^ c[2] ^ c[6];
	hl2[3] = hl2[3] ^ c[3] ^ c[7];

	c[2] ^= hl2[0];
	c[3] ^= hl2[1];
	c[4] ^= hl2[2];
	c[5] ^= hl2[3];
}



void mul4 (gf2n_word *c, gf2n_word *a, gf2n_word *b)
{
	imul4(c, a, b);
}



void mul5 (gf2n_word *c, gf2n_word *a, gf2n_word *b)
{
	gf2n_word hs0[3], hs1[3];
	gf2n_word hl2[6];

	hs0[0] = a[0] ^ a[3]; // a0 + a1
	hs0[1] = a[1] ^ a[4];
	hs0[2] = a[2];
	hs1[0] = b[0] ^ b[3]; // b0 + b1
	hs1[1] = b[1] ^ b[4];
	hs1[2] = b[2];

	imul3(c, a, b); // a0 * b0
	imul2(c+6, a+3, b+3); // a1 * b1
	imul3(hl2, hs0, hs1); // () * ()

	hl2[0] = hl2[0] ^ c[0] ^ c[6]; // ()() - a0b0 -a1b1
	hl2[1] = hl2[1] ^ c[1] ^ c[7];
	hl2[2] = hl2[2] ^ c[2] ^ c[8];
	hl2[3] = hl2[3] ^ c[3] ^ c[9];
	hl2[4] = hl2[4] ^ c[4];
	hl2[5] = hl2[5] ^ c[5];


	c[3] ^= hl2[0];
	c[4] ^= hl2[1];
	c[5] ^= hl2[2];
	c[6] ^= hl2[3];
	c[7] ^= hl2[4];
	c[8] ^= hl2[5];
}



void mul6 (gf2n_word *c, gf2n_word *a, gf2n_word *b)
{
	gf2n_word hs0[3], hs1[3];
	gf2n_word hl2[6];

	hs0[0] = a[0] ^ a[3]; // a0 + a1
	hs0[1] = a[1] ^ a[4];
	hs0[2] = a[2] ^ a[5];
	hs1[0] = b[0] ^ b[3]; // b0 + b1
	hs1[1] = b[1] ^ b[4];
	hs1[2] = b[2] ^ b[5];

	imul3(c, a, b); // a0 * b0
	imul3(c+6, a+3, b+3); // a1 * b1
	imul3(hl2, hs0, hs1); // () * ()

	hl2[0] = hl2[0] ^ c[0] ^ c[6]; // ()() + a0b0 + a1b1
	hl2[1] = hl2[1] ^ c[1] ^ c[7];
	hl2[2] = hl2[2] ^ c[2] ^ c[8];
	hl2[3] = hl2[3] ^ c[3] ^ c[9];
	hl2[4] = hl2[4] ^ c[4] ^ c[10];
	hl2[5] = hl2[5] ^ c[5] ^ c[11];

	c[3] ^= hl2[0];
	c[4] ^= hl2[1];
	c[5] ^= hl2[2];
	c[6] ^= hl2[3];
	c[7] ^= hl2[4];
	c[8] ^= hl2[5];
}



void mul7 (gf2n_word *c, gf2n_word *a, gf2n_word *b)
{
	gf2n_word hs0[4], hs1[4];
	gf2n_word hl2[8];

	hs0[0] = a[0] ^ a[4]; // a0 + a1
	hs0[1] = a[1] ^ a[5];
	hs0[2] = a[2] ^ a[6];
	hs0[3] = a[3];
	hs1[0] = b[0] ^ b[4]; // b0 + b1
	hs1[1] = b[1] ^ b[5];
	hs1[2] = b[2] ^ b[6];
	hs1[3] = b[3];

	imul4(c, a, b); // a0 * b0
	imul3(c+8, a+4, b+4); // a1 * b1
	imul4(hl2, hs0, hs1); // () * ()

	hl2[0] = hl2[0] ^ c[0] ^ c[8]; // ()() + a0b0 + a1b1
	hl2[1] = hl2[1] ^ c[1] ^ c[9];
	hl2[2] = hl2[2] ^ c[2] ^ c[10];
	hl2[3] = hl2[3] ^ c[3] ^ c[11];
	hl2[4] = hl2[4] ^ c[4] ^ c[12];
	hl2[5] = hl2[5] ^ c[5] ^ c[13];
	hl2[6] = hl2[6] ^ c[6];
	hl2[7] = hl2[7] ^ c[7];

	c[4] ^= hl2[0];
	c[5] ^= hl2[1];
	c[6] ^= hl2[2];
	c[7] ^= hl2[3];
	c[8] ^= hl2[4];
	c[9] ^= hl2[5];
	c[10] ^= hl2[6];
	c[11] ^= hl2[7];
}



void mul8 (gf2n_word *c, gf2n_word *a, gf2n_word *b)
{
	gf2n_word hs0[4], hs1[4];
	gf2n_word hl2[8];

	hs0[0] = a[0] ^ a[4]; // a0 + a1
	hs0[1] = a[1] ^ a[5];
	hs0[2] = a[2] ^ a[6];
	hs0[3] = a[3] ^ a[7];
	hs1[0] = b[0] ^ b[4]; // b0 + b1
	hs1[1] = b[1] ^ b[5];
	hs1[2] = b[2] ^ b[6];
	hs1[3] = b[3] ^ b[7];

	imul4(c, a, b); // a0 * b0
	imul4(c+8, a+4, b+4); // a1 * b1
	imul4(hl2, hs0, hs1); // () * ()

	hl2[0] = hl2[0] ^ c[0] ^ c[8]; // ()() + a0b0 + a1b1
	hl2[1] = hl2[1] ^ c[1] ^ c[9];
	hl2[2] = hl2[2] ^ c[2] ^ c[10];
	hl2[3] = hl2[3] ^ c[3] ^ c[11];
	hl2[4] = hl2[4] ^ c[4] ^ c[12];
	hl2[5] = hl2[5] ^ c[5] ^ c[13];
	hl2[6] = hl2[6] ^ c[6] ^ c[14];
	hl2[7] = hl2[7] ^ c[7] ^ c[15];

	c[4] ^= hl2[0];
	c[5] ^= hl2[1];
	c[6] ^= hl2[2];
	c[7] ^= hl2[3];
	c[8] ^= hl2[4];
	c[9] ^= hl2[5];
	c[10] ^= hl2[6];
	c[11] ^= hl2[7];
}



void kara (gf2n_word *c, gf2n_word *a, gf2n_word *b, unsigned int length)
{
	unsigned int ll, lh, i, ll2, lh2;
	gf2n_word *a0, *a1, *b0, *b1, *a01, *b01, *h;

	if (length < FIXMUL)
		gf2nmul[length](c, a, b);
	else {
		if (length & 1)       // if length odd
		{
			lh = length >> 1;
			ll = lh + 1;
		}
		else
			lh = ll = length >> 1;

		ll2 = ll << 1;
		lh2 = lh << 1;

		a01 = new gf2n_word[ll+1];
		b01 = new gf2n_word[ll+1];
		h = new gf2n_word[ll2+1];

		a0 = a;
		a1 = a+ll;
		b0 = b;
		b1 = b+ll;

		kara(c , a0, b0, ll); // c_low = a0*b0
		kara(c+ll2, a1, b1, lh); // c_high = a1*b1

		for (i = 0; i < lh; i++) {
			a01[i] = a[i] ^ a[i+ll]; // a01 = a0+a1
			b01[i] = b[i] ^ b[i+ll]; // b01 = b0+b1
		}

		if (lh < ll) {
			a01[lh] = a[lh]; // a01[ll-1] = a[ll-1], since ll-1 == lh
			b01[lh] = b[lh];
		}

		kara(h, a01, b01, ll); // h = (a0+a1)(b0+b1)

		for (i = 0; i < ll2; i++)       // h = h + a0b0
			h[i] ^= c[i];

		for (i = 0; i < lh2; i++)       // h = h + a1b1
			h[i] ^= c[i+ll2];

		for (i = 0; i < ll2; i++)       // c = c + h*2^ll
			c[i+ll] ^= h[i];

		delete [] h;
		delete [] b01;
		delete [] a01;
	}
}



void karatsuba_mul (gf2n_word *c, gf2n_word *a, gf2n_word *b)
{
	kara(c, a, b, gf2n::anzBI);
}



// ======================================================================
// MULTIPLICATION SELECTOR                                            
// ======================================================================

void (*gf2nmul[]) (gf2n_word*, gf2n_word*, gf2n_word*) =
{
	mul1, mul1, mul2, mul3, mul4, mul5, mul6, mul7, mul8, karatsuba_mul
};



// ======================================================================
// SQUARING                                
// ======================================================================

#if GF2N_WORDSIZE == 64
void square (gf2n_word *c, gf2n_word *a)
{
	for (register int i = gf2n::anzBI-1, j = 2*gf2n::anzBI-1; i >= 0; i--) {
		c[j--] = tabsquare[ gf2n_bit16(a[i] >> 48) ] << 32 |
			tabsquare[ gf2n_bit16(a[i] >> 32) ];
		c[j--] = tabsquare[ gf2n_bit16(a[i] >> 16) ] << 32 |
			tabsquare[ gf2n_bit16(a[i]) ];
	}
}
#else // we are running on a 32bit architecture
void square (gf2n_word *c, gf2n_word *a)
{
	for (register int i = gf2n::anzBI-1, j = 2*gf2n::anzBI-1; i >= 0; i--) {
		c[j--] = tabsquare[ gf2n_bit16(a[i] >> 16) ];
		c[j--] = tabsquare[ gf2n_bit16(a[i]) ];
	}
}
#endif



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
