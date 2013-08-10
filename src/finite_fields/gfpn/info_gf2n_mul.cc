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
// 		  Volker Mueller(VM), Thomas Pfahler (TPf)
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



// FIXME: not portable
typedef unsigned short gf2n_bit16; // oh well...
typedef unsigned char  gf2n_bit8;



static const int dim_mul = 256; // dimension for multiplication table
static const int dim_square = 65536; // dimension for squaring table

static gf2n_bit16 tabmul[dim_mul][dim_mul]; // table for multiplication
static udigit tabsquare[dim_square]; // table for squaring



//------------------------------------------------------------------------

//*************************************************************
// multiply two 8-bit numbers without carrys
//*************************************************************
static gf2n_bit16 mul8bit (gf2n_bit16 a, gf2n_bit16 b)
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
udigit mul16bit(gf2n_bit16 a, gf2n_bit16 b)
{
	gf2n_bit8 a0, a1, b0, b1;

	a0 = gf2n_bit8(a);
	a1 = gf2n_bit8(a >> 8);
	b0 = gf2n_bit8(b);
	b1 = gf2n_bit8(b >> 8);

	return(tabmul[a0][b0] ^
	       udigit((tabmul[a0][b1] ^ tabmul[a1][b0]) << 8) ^
	       udigit(tabmul[a1][b1] << 16));
}



//***************************************************************

void info_gf2n::gen_tables()
{
	int i, j;
	// was: generate_mul_table()
	for (i = 0; i < dim_mul; i++)
		tabmul[i][i] = mul8bit(i, i);
	for (i = 0; i < dim_mul-1; i++)
		for (j = i+1; j < dim_mul; j++)
			tabmul[j][i] = tabmul[i][j] = mul8bit(i, j);

	// was: generate_square_table()
	for (i = 0; i < dim_square; i++)
		tabsquare[i] = mul16bit(i, i);

}



//**************************************************************
// multiplication of arrays of udigits
//
//     mulj(...)->input is array [0, .., j-1] of udigit
//     imulj(..)  is corresponding function as inline
//*************************************************************

inline
void imul1 (udigit *c, udigit *a, udigit *b)
{
	udigit cm;
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



static void mul1_noinl(udigit *c, udigit *a, udigit *b)
{
	imul1(c, a, b);
}



void info_gf2n::mul1(udigit *c, udigit *a, udigit *b) const
{
	imul1(c, a, b);
}



inline
void imul2 (udigit *c, udigit *a, udigit *b)
{
	udigit hs0, hs1;
	udigit hl2[2];

	hs0 = a[0] ^ a[1]; // a0 + a1
	hs1 = b[0] ^ b[1]; // b0 + b1

	mul1_noinl(c, a, b);
	mul1_noinl(c+2, a+1, b+1);
	mul1_noinl(hl2, & hs0, & hs1);

	hl2[0] = hl2[0] ^ c[0] ^ c[2];
	hl2[1] = hl2[1] ^ c[1] ^ c[3];

	c[1] ^= hl2[0];
	c[2] ^= hl2[1];
}



void info_gf2n::mul2 (udigit *c, udigit *a, udigit *b) const
{
	imul2(c, a, b);
}



inline
void imul3 (udigit *c, udigit *a, udigit *b)
{
	udigit hs0[2], hs1[2];
	udigit hl2[4];

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



void info_gf2n::mul3 (udigit *c, udigit *a, udigit *b) const
{
	imul3(c, a, b);
}



inline
void imul4 (udigit *c, udigit *a, udigit *b)
{
	udigit hs0[2], hs1[2];
	udigit hl2[4];

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



void info_gf2n::mul4 (udigit *c, udigit *a, udigit *b) const
{
	imul4(c, a, b);
}



void info_gf2n::mul5 (udigit *c, udigit *a, udigit *b) const
{
	udigit hs0[3], hs1[3];
	udigit hl2[6];

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



void info_gf2n::mul6 (udigit *c, udigit *a, udigit *b) const
{
	udigit hs0[3], hs1[3];
	udigit hl2[6];

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



void info_gf2n::mul7 (udigit *c, udigit *a, udigit *b) const
{
	udigit hs0[4], hs1[4];
	udigit hl2[8];

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



void info_gf2n::mul8 (udigit *c, udigit *a, udigit *b) const
{
	udigit hs0[4], hs1[4];
	udigit hl2[8];

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



struct kara_struct
{
	const info_gf2n &info;

	kara_struct(const info_gf2n &a) : info(a) {}
	~kara_struct() {}

	void kara(udigit *, udigit *, udigit *, unsigned int) const;

	// inhibit:
private:
	kara_struct();
	kara_struct & operator = (const kara_struct &);
};



void
kara_struct::kara (udigit *c, udigit *a, udigit *b, unsigned int length) const
{
	unsigned int ll, lh, i, ll2, lh2;
	udigit *a0, *a1, *b0, *b1, *a01, *b01, *h;

	if (length < info.FIXMUL)
		(info.*info_gf2n::gf2nmul[length])(c, a, b);
	else {
		if (length & 1) {
			// if length odd
			lh = length >> 1;
			ll = lh + 1;
		}
		else
			lh = ll = length >> 1;

		ll2 = ll << 1;
		lh2 = lh << 1;

		a01 = new udigit[ll+1];
		b01 = new udigit[ll+1];
		h = new udigit[ll2+1];

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



void info_gf2n::karatsuba_mul (udigit *c, udigit *a, udigit *b) const
{
	kara_struct k(*this);
	k.kara(c, a, b, anzBI);
}



// ======================================================================
// MULTIPLICATION SELECTOR                                            
// ======================================================================

void (info_gf2n::*(info_gf2n::gf2nmul)[]) (udigit*, udigit*, udigit*) const =
{
	&info_gf2n::mul1, &info_gf2n::mul1, &info_gf2n::mul2, &info_gf2n::mul3,
	&info_gf2n::mul4, &info_gf2n::mul5, &info_gf2n::mul6, &info_gf2n::mul7,
	&info_gf2n::mul8, &info_gf2n::karatsuba_mul
};



// ======================================================================
// SQUARING                                
// ======================================================================

void info_gf2n::square (udigit *c, udigit *a) const
{
	for (register int i = anzBI-1, j = 2*anzBI-1; i >= 0; i--) {
		c[j--] = tabsquare[ gf2n_bit16(a[i] >> 16) ];
		c[j--] = tabsquare[ gf2n_bit16(a[i]) ];
	}
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
