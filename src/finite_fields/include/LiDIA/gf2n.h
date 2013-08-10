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
//	Author	: Volker Mueller (VM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_GF2N_H_GUARD_
#define LIDIA_GF2N_H_GUARD_


#ifndef LIDIA_GF2NIO_H_GUARD_
# include	"LiDIA/finite_fields/gf2nIO.h"
#endif
#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_RATIONAL_FACTORIZATION_H_GUARD_
# include	"LiDIA/rational_factorization.h"
#endif
#ifndef LIDIA_POWER_FUNCTIONS_H_GUARD_
# include	"LiDIA/power_functions.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class gf2n_polynomial;
class gf2nIO;

const unsigned int FIXMUL = 9; // number of multiplication routines

// we now use gf2n_word instead of gf2n_word -rpw.
typedef unsigned long gf2n_word;

// we only allow for 32 bit or 64 bit longs. -rpw
#if SIZEOF_LONG == 4
# define GF2N_WORDSIZE	32
#else
# if SIZEOF_LONG == 8
#  define GF2N_WORDSIZE	64
# else
#  error "word size not supported"
# endif
#endif

typedef unsigned short gf2n_bit16;
typedef unsigned char  gf2n_bit8;


class gf2n
{
private:

	static unsigned int anzBI; // number of gf2n_words for each gf2n  
	static unsigned int degree; // degree of extension over F_2
	static rational_factorization ord; // factorization of group order or 1
	static bool ord_is_fact; // true, if ord is set, false otherwise
	static unsigned int *exponents; // array of exponents with 1 coefficients
	// of generating polynomial, excluded
	// degree and 0
	static unsigned int exp1, exp2, exp3; // for storing the exponents of   
	// 3/5-nomials
	static unsigned int anz_exponents; // number of exponents
	static unsigned int mulsel; // Multiplication-selector
	static unsigned int invsel; // Inversion and Reduction selector
	static gf2n_word *B, *C, *F, *G; // used in inversion

	static gf2n* table_solve_quadratic; // table for fast solve_quadratic
	static bool table_solve_quadratic_init; // true iff table initialized
	static bigint * table_solve_quadratic_masks;
	static int table_solve_quadratic_number_masks;

	gf2n_word *element; // array representing element


	//=============================================================
	// PRIVATE FRIEND & STATIC FUNCTIONS, USED IN ARITHMETIC     
	//=============================================================

	bool is_reduced() const
	{
		return ((element[anzBI-1] >> (degree % (CHAR_BIT * sizeof(gf2n_word)))) == (gf2n_word) 0);
	}

	friend void karatsuba_mul (gf2n_word *, gf2n_word *, gf2n_word*);
	friend void square(gf2n_word*, gf2n_word*);
	friend void tri_partial_reduce1(gf2n_word*);
	friend void pent_partial_reduce1(gf2n_word*);
	friend void general_partial_reduce1(gf2n_word*);
	friend void tri_partial_reduce2(gf2n_word*);
	friend void pent_partial_reduce2(gf2n_word*);
	friend void general_partial_reduce2(gf2n_word*);
	friend void tri_invert(gf2n_word*, gf2n_word*);
	friend void pent_invert(gf2n_word*, gf2n_word*);
	friend void general_invert(gf2n_word*, gf2n_word*);

	static void generate_mul_table();
	static void generate_square_table();

	static char* read_modulus_from_database(unsigned int);


public:

	//=====================================================
	// FRIENDS                                 
	//=====================================================

	friend class gf2nIO;

	//=====================================================
	// INITIALISATION                          
	//=====================================================


	void re_initialize();

	friend void gf2n_init(char *, gf2nIO::base);

	friend void gf2n_init(unsigned int, gf2nIO::base);

	static unsigned int extension_degree()
	{
		return degree;
	}
	static unsigned int get_absolute_degree()
	{
		return degree;
	}

	static void initialize_table_for_solve_quadratic();

	static void delete_table_for_solve_quadratic();

	// the following function fucntion uses a table to compute a root c
	// of x^2 + a1*x + a0. If certainly_factors is true, then it is not
	// checked whether that polynomial factors, but that is assumed. It is
	// a little bit faster, but the functions returns wrong results if indeed
	// the polynomial does not factor. If certainly == false, then it is checked
	// whether the polynomial has a root, if so, true is returned and a root is
	// assigned to c, otherwise false is returned.

	friend bool solve_quadratic_with_table(gf2n & c, const gf2n & a1,
					       const gf2n & a0,
					       bool certainly_factors);




	//=====================================================
	// CONSTRUCTORS, DESTRUCTOR
	//=====================================================



	gf2n();

	gf2n(const gf2n &);

	gf2n(unsigned long);

	gf2n(const bigint &);

	~gf2n();

	//=====================================================
	// ASSIGNMENTS                          
	//=====================================================


	gf2n & operator = (const gf2n &);
	gf2n & operator = (unsigned long);
	gf2n & operator = (const bigint &);
	void assign_zero();
	void assign_one();
	void assign(const gf2n &);
	friend void swap (gf2n &, gf2n &);



	//=====================================================
	// COMPARISONS
	//=====================================================

	friend bool operator == (const gf2n &, const gf2n &);

	friend bool operator != (const gf2n &, const gf2n &);

	bool is_zero() const;

	bool is_one() const;

	bool is_square() const    // for consistency with gf_p_element
	{
		return true;
	}


	//=====================================================
	// IO
	//=====================================================


	friend std::istream & operator >> (std::istream &, gf2n &);

	friend std::ostream & operator << (std::ostream &, const gf2n &);


	//=====================================================
	// ARITHMETIC
	//=====================================================

        friend gf2n operator - (const gf2n & a);

	friend gf2n operator * (const gf2n &, const gf2n &);

	friend gf2n operator + (const gf2n &, const gf2n &);

	friend gf2n operator - (const gf2n &, const gf2n &);

	friend gf2n operator / (const gf2n &, const gf2n &);


	gf2n & operator += (const gf2n &);

	gf2n & operator -= (const gf2n &);

	gf2n & operator *= (const gf2n &);

	gf2n & operator /= (const gf2n &);

	friend void add(gf2n &, const gf2n &, const gf2n &);

	friend void subtract(gf2n &, const gf2n &, const gf2n &);

	friend void multiply(gf2n &, const gf2n &, const gf2n &);

	friend void square(gf2n &, const gf2n &);

	friend void divide(gf2n &, const gf2n &, const gf2n &);


	friend void power(gf2n & c, const gf2n & a, const long x);
	friend void power(gf2n & c, const gf2n & a, const bigint & x);

	friend void sqrt(gf2n &, const gf2n &);

	friend gf2n sqrt(const gf2n &);

	void randomize(unsigned int d = gf2n::degree);
	friend gf2n randomize(const gf2n & x, unsigned int d);

	void invert();

	friend void invert (gf2n &, const gf2n &);

	friend gf2n inverse(const gf2n &);

	friend bigint compute_order(const gf2n &);

	friend gf2n get_generator(unsigned int d);

	friend gf2n_word hash(const gf2n &);

	unsigned int relative_degree() const;

	int trace() const;

};

// declaration of gf2n's friend functions
void gf2n_init(char *, gf2nIO::base = gf2nIO::Dec);
void gf2n_init(unsigned int, gf2nIO::base = gf2nIO::Dec);
bool solve_quadratic_with_table(gf2n & c, const gf2n & a1,
				const gf2n & a0,
				bool certainly_factors = false);
void swap (gf2n &, gf2n &);
bool operator == (const gf2n &, const gf2n &);
bool operator != (const gf2n &, const gf2n &);
std::istream & operator >> (std::istream &, gf2n &);
std::ostream & operator << (std::ostream &, const gf2n &);

    //=====================================================
    // ARITHMETIC
    //=====================================================

inline
gf2n operator - (const gf2n & a)
{
    return a;
}

gf2n operator * (const gf2n &, const gf2n &);
gf2n operator + (const gf2n &, const gf2n &);
gf2n operator - (const gf2n &, const gf2n &);
gf2n operator / (const gf2n &, const gf2n &);

void add(gf2n &, const gf2n &, const gf2n &);
void subtract(gf2n &, const gf2n &, const gf2n &);
void multiply(gf2n &, const gf2n &, const gf2n &);
void square(gf2n &, const gf2n &);
void divide(gf2n &, const gf2n &, const gf2n &);

void power(gf2n & c, const gf2n & a, const long x);
void power(gf2n & c, const gf2n & a, const bigint & x);
void sqrt(gf2n &, const gf2n &);
gf2n sqrt(const gf2n &);

gf2n randomize(const gf2n & x, unsigned int d = gf2n::extension_degree());
void invert (gf2n &, const gf2n &);
gf2n inverse(const gf2n &);

bigint compute_order(const gf2n &);
gf2n get_generator(unsigned int d = gf2n::extension_degree());

gf2n_word hash(const gf2n &);


//======================================================================
// INLINE PROCEDURAL VERSIONS
//======================================================================


extern void (*gf2nmul[]) (gf2n_word*, gf2n_word*, gf2n_word*);
extern void (*partial_reduce1[]) (gf2n_word *);
extern void (*partial_reduce2[]) (gf2n_word *);
extern void (*uinvert[]) (gf2n_word *, gf2n_word *);



//
// MULTIPLICATION
// ----------------------------------------------------------------------

inline void
multiply(gf2n & c, const gf2n & a, const gf2n & b)
{
	void square(gf2n_word*, gf2n_word*);

	register unsigned int i;

	for (i = 0; i < 2*gf2n::anzBI; i++)
		gf2n::B[i] = (gf2n_word) 0;

	if (a.element == b.element) {
		square(gf2n::B, a.element);
		partial_reduce1[gf2n::invsel](gf2n::B);
	}
	else {
		gf2nmul[gf2n::mulsel](gf2n::B, a.element, b.element);
		partial_reduce1[gf2n::invsel](gf2n::B);
	}
	for (i = 0; i < gf2n::anzBI; i++)
		c.element[i] = gf2n::B[i];
}

//
// SQUARING
// ----------------------------------------------------------------------

inline void
square(gf2n & c, const gf2n & a)
{
	void square(gf2n_word*, gf2n_word*);

	square(gf2n::B, a.element);
	partial_reduce1[gf2n::invsel](gf2n::B);
	for (register unsigned i = 0; i < gf2n::anzBI; i++)
		c.element[i] = gf2n::B[i];
}


//
//   INVERSION & DIVISION
// ----------------------------------------------------------------------

inline void
invert (gf2n & b, const gf2n & a)
{
	if (a.is_one () == true) {
		b.assign_one ();
	}
	else {
		if (a.is_zero () == true)
			lidia_error_handler("gf2n", "invert(gf2n&, const gf2n&)::zero not invertible\n");

		uinvert[gf2n::invsel](b.element, a.element);
	}
}




inline void
divide(gf2n & c, const gf2n & a, const gf2n & b)
{
	gf2n d;
	LiDIA::invert(d, b);
	multiply(c, a, d);
}



//
// ADDITION
// ----------------------------------------------------------------------

inline void
add(gf2n & c, const gf2n & a, const gf2n & b)
{
	register unsigned int i;
	register gf2n_word *cp, *ap, *bp;

	for (i = 0, ap = a.element, bp = b.element, cp = c.element;
	     i < gf2n::anzBI; i++, cp++, ap++, bp++)
		*cp = (*ap) ^ (*bp);
}

//
// SUBTRACTION
// ----------------------------------------------------------------------

inline void
subtract(gf2n & c, const gf2n & a, const gf2n & b)
{
	add(c, a, b);
}



//
//   gf2n_lib
//


//**** relative degree of element over prime field ***************
unsigned int relative_degree(const gf2n & a);

//***** absolute degree over prime field *************************
unsigned int get_absolute_degree(const gf2n & x);

//****** solve_quadratic *****************************************
bool solve_quadratic(gf2n & root, const gf2n & a1, const gf2n & a0);

//****** characteristic ******************************************
bigint characteristic(const gf2n & a);

//***** number_of_elements ***************************************
bigint number_of_elements(const gf2n & a);



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_GF2N_H_GUARD_
