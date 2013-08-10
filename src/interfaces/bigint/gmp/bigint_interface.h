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
//	Author	:Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


//
// c'tors and d'tor
//

inline
bigint::bigint ()
{
	mpz_init(&I);
}



inline
bigint::bigint (int i)
{
	mpz_init_set_si(&I, i);
}



inline
bigint::bigint (long l)
{
	mpz_init_set_si(&I, l);
}



inline
bigint::bigint (unsigned long ul)
{
	mpz_init_set_ui(&I, ul);
}



inline
bigint::bigint (double d)
{
	mpz_init_set_d(&I, d);
}



inline
bigint::bigint (const bigint_rep_t & a)
{
	mpz_init_set(&I, &a);
}



inline
bigint::bigint (const bigint & a)
{
	mpz_init_set(&I, &a.I);
}



inline
bigint::~bigint()
{
	mpz_clear(&I);
}



//
// accessors
//

inline int
bigint::bit (unsigned long idx) const
{
	return mpz_tstbit(&I, idx);
}



inline unsigned long
bigint::most_significant_digit () const
{
	size_t	i = mpz_size(&I);

	return ((i == 0) ? 0 : mpz_getlimbn(&I, i - 1));
}



inline unsigned long
bigint::least_significant_digit () const
{
	return ((mpz_size(&I) == 0) ? 0 : mpz_getlimbn(&I, 0));
}



//
// assigners
//

inline void
bigint::assign_zero ()
{
	mpz_set_ui(&I, 0UL);
}



inline void
bigint::assign_one ()
{
	mpz_set_ui(&I, 1UL);
}



inline void
bigint::assign (int i)
{
	mpz_set_si(&I, i);
}



inline void
bigint::assign (long i)
{
	mpz_set_si(&I, i);
}



inline void
bigint::assign (unsigned long ui)
{
	mpz_set_ui(&I, ui);
}



inline void
bigint::assign (double d)
{
	mpz_set_d(&I, d);
}



inline void
bigint::assign (const bigint_rep_t & a)
{
	if (&a != &I) {
		mpz_set(&I, &a);
	}
}



inline void
bigint::assign (const bigint & a)
{
	if (&a != this) {
		mpz_set(&I, &a.I);
	}
}



//
// comparators
//

inline int
bigint::abs_compare (const bigint & a) const
{
	if (&a == this) {
		return 0;
	}

	int cmp = mpz_cmpabs(&I, &a.I);

	return ((cmp < 0) ? -1 : (cmp > 0));
}



inline int
bigint::abs_compare (unsigned long a) const
{
	int cmp = mpz_cmpabs_ui(&I, a);

	return ((cmp < 0) ? -1 : (cmp > 0));
}



inline int
bigint::compare (const bigint & a) const
{
	if (&a == this) {
		return 0;
	}

	int cmp = mpz_cmp(&I, &a.I);

	return ((cmp < 0) ? -1 : (cmp > 0));
}



inline int
bigint::compare (unsigned long a) const
{
	int cmp = mpz_cmp_ui(&I, a);

	return ((cmp < 0) ? -1 : (cmp > 0));
}



inline int
bigint::compare (long a) const
{
	int cmp = mpz_cmp_si(&I, a);

	return ((cmp < 0) ? -1 : (cmp > 0));
}



inline int
bigint::sign () const
{
	return mpz_sgn(&I);
}



//
// properties
//

inline lidia_size_t
bigint::length () const
{
	return mpz_size(&I);
}



inline lidia_size_t
bigint::bit_length () const
{
	return mpz_sizeinbase(&I, 2);
}



//
// predicates
//

inline bool
bigint::is_even () const
{
	return mpz_even_p(&I);
}



inline bool
bigint::is_odd () const
{
	return mpz_odd_p(&I);
}




//
// converters
//

inline bool
bigint::intify (int & i) const
{
	if (!mpz_fits_sint_p(&I)) {
		return true;
	}
	i = static_cast<int>(mpz_get_si(&I));
	return false;
}



inline bool
bigint::longify (long & l) const
{
	if (!mpz_fits_slong_p(&I)) {
		return true;
	}
	l = mpz_get_si(&I);
	return false;
}



inline double
bigint::dbl () const
{
	return mpz_get_d(&I);
}



//
// misc. functions
//

inline double
bigint::radix ()
{
	return std::ldexp(1.0, mp_bits_per_limb);
}



inline int
bigint::bits_per_digit ()
{
	return mp_bits_per_limb;
}



//
// modifiers
//

inline void
bigint::absolute_value ()
{
	mpz_abs(&I, &I);
}



inline void
bigint::abs ()
{
	mpz_abs(&I, &I);
}



inline void
bigint::negate ()
{
	mpz_neg(&I, &I);
}



inline void
bigint::inc ()
{
	mpz_add_ui(&I, &I, 1);
}



inline void
bigint::dec ()
{
	mpz_sub_ui(&I, &I, 1);
}



inline void
bigint::multiply_by_2 ()
{
	mpz_mul_2exp(&I, &I, 1);
}



inline void
bigint::divide_by_2 ()
{
	mpz_tdiv_q_2exp(&I, &I, 1);
}



inline void
bigint::swap (bigint & b)
{
	mpz_swap(&I, &b.I);
}



//
// type predicates
//

inline bool
bigint::is_char () const
{
	return ((mpz_sgn(&I) >= 0) ? (mpz_cmp_si(&I, SCHAR_MAX) <= 0) : (mpz_cmp_si(&I, SCHAR_MIN) >= 0));
}



inline bool
bigint::is_uchar () const
{
	return ((mpz_sgn(&I) >= 0) && (mpz_cmp_ui(&I, UCHAR_MAX) <= 0));
}



inline bool
bigint::is_short () const
{
	return mpz_fits_sshort_p(&I);
}



inline bool
bigint::is_ushort () const
{
	return mpz_fits_ushort_p(&I);
}



inline bool
bigint::is_int () const
{
	return mpz_fits_sint_p(&I);
}



inline bool
bigint::is_uint () const
{
	return mpz_fits_uint_p(&I);
}



inline bool
bigint::is_long () const
{
	return mpz_fits_slong_p(&I);
}



inline bool
bigint::is_ulong () const
{
	return mpz_fits_ulong_p(&I);
}



//
// arithmetic procedures
//

inline void
negate (bigint & a, const bigint & b)
{
	mpz_neg(&a.I, &b.I);
}



inline void
add (bigint & c, const bigint & a, const bigint & b)
{
	mpz_add(&c.I, &a.I, &b.I);
}



inline void
add (bigint & c, const bigint & a, unsigned long b)
{
	mpz_add_ui(&c.I, &a.I, b);
}



inline void
add (bigint & c, const bigint & a, long b)
{
	if (b >= 0) {
		mpz_add_ui(&c.I, &a.I, static_cast<unsigned long>(b));
	}
	else {
		mpz_sub_ui(&c.I, &a.I, static_cast<unsigned long>(-b));
	}
}



inline void
subtract (bigint & c, const bigint & a, const bigint & b)
{
	mpz_sub(&c.I, &a.I, &b.I);
}



inline void
subtract (bigint & c, const bigint & a, unsigned long b)
{
	mpz_sub_ui(&c.I, &a.I, b);
}



inline void
subtract (bigint & c, const bigint & a, long b)
{
	if (b >= 0) {
		mpz_sub_ui(&c.I, &a.I, b);
	}
	else {
		mpz_add_ui(&c.I, &a.I, static_cast<unsigned long>(-b));
	}
}



inline void
multiply (bigint & c, const bigint & a, const bigint & b)
{
	mpz_mul(&c.I, &a.I, &b.I);
}



inline void
multiply (bigint & c, const bigint & a, unsigned long b)
{
	mpz_mul_ui(&c.I, &a.I, b);
}



inline void
multiply (bigint & c, const bigint & a, long b)
{
	mpz_mul_si(&c.I, &a.I, b);
}



inline void
square (bigint & a, const bigint & b)
{
	mpz_mul(&a.I, &b.I, &b.I);
}



inline void
divide (bigint & c, const bigint & a, const bigint & b)
{
	mpz_tdiv_q(&c.I, &a.I, &b.I);
}



inline void
divide (bigint & c, const bigint & a, unsigned long b)
{
	mpz_tdiv_q_ui(&c.I, &a.I, b);
}



inline void
divide (bigint & c, const bigint & a, long b)
{
	if (b >= 0) {
		mpz_tdiv_q_ui(&c.I, &a.I, static_cast<unsigned long>(b));
	}
	else {
		mpz_tdiv_q_ui(&c.I, &a.I, static_cast<unsigned long>(-b));
		mpz_neg(&c.I, &c.I);
	}
}



inline void
remainder (bigint & c, const bigint & a, const bigint & b)
{
	mpz_tdiv_r(&c.I, &a.I, &b.I);
}




inline void
remainder (bigint & c, const bigint & a, unsigned long b)
{
	mpz_tdiv_r_ui(&c.I, &a.I, b);
}




inline void
remainder (bigint & c, const bigint & a, long b)
{
	mpz_tdiv_r_ui(&c.I, &a.I, static_cast<unsigned long>(std::labs(b)));
}




inline void
remainder (long & r, const bigint & a, long b)
{
	if (b > 0) {
		r = static_cast<long>(mpz_tdiv_ui(&a.I, static_cast<unsigned long>(b)));
	}
	else {
		r = static_cast<long>(mpz_tdiv_ui(&a.I, static_cast<unsigned long>(-b)));
	}

	if (r != 0 && mpz_sgn(&a.I) < 0) {
		// r must be negative
		r = -r;
	}
}



inline void
remainder (unsigned long & r, const bigint & a, unsigned long b)
{
	r = mpz_tdiv_ui(&a.I, b);

	if (r != 0 && mpz_sgn(&a.I) < 0) {
		// r must be negative
		r = -r;
	}
}



inline long
remainder (const bigint & a, long b)
{
	long r;

	if (b > 0) {
		r = static_cast<long>(mpz_tdiv_ui(&a.I, static_cast<unsigned long>(b)));
	}
	else {
		r = static_cast<long>(mpz_tdiv_ui(&a.I, static_cast<unsigned long>(-b)));
	}

	if (r != 0 && mpz_sgn(&a.I) < 0) {
		// r must be negative
		r = -r;
	}

	return r;
}



inline void
div_rem (bigint & q, bigint & r, const bigint & a, const bigint & b)
{
	mpz_tdiv_qr(&q.I, &r.I, &a.I, &b.I);
}



inline void
div_rem (bigint & q, bigint & r, const bigint & a, unsigned long b)
{
	mpz_tdiv_qr_ui(&q.I, &r.I, &a.I, b);
}



inline void
div_rem (bigint & q, bigint & r, const bigint & a, long b)
{
	if (b >= 0) {
		mpz_tdiv_qr_ui(&q.I, &r.I, &a.I, static_cast<unsigned long>(b));
	}
	else {
		mpz_tdiv_qr_ui(&q.I, &r.I, &a.I, static_cast<unsigned long>(-b));
		mpz_neg(&q.I, &q.I);
	}
}



inline void
div_rem (bigint & q, long & r, const bigint & a, long b)
{
	if (b >= 0) {
		r = static_cast<long>(mpz_tdiv_q_ui(&q.I, &a.I, static_cast<unsigned long>(b)));
	}
	else {
		r = static_cast<long>(mpz_tdiv_q_ui(&q.I, &a.I, static_cast<unsigned long>(-b)));
		mpz_neg(&q.I, &q.I);
	}

	if (r != 0 && mpz_sgn(&a.I) < 0) {
		// r must be negative
		r = -r;
	}
}



inline void
invert (bigint & a, const bigint & b)
{
	if (mpz_cmpabs_ui(&b.I, 1) == 0) {
		mpz_set(&a.I, &b.I);
	}
	else {
		lidia_error_handler("bigint", "invert::inverting of a non-unit.");
	}
}



inline void
sqrt (bigint & a, const bigint & b)
{
	if (mpz_sgn(&b.I) < 0) {
		lidia_error_handler("bigint", "sqrt(bigint&a, const bigint&b):: b < 0");
	}
	else {
		mpz_sqrt(&a.I, &b.I);
	}
}



inline void
shift_left (bigint & c, const bigint & a, long ui)
{
	if (ui < 0)
		lidia_error_handler("bigint", "shift_left()::index is negative.");
	mpz_mul_2exp(&c.I, &a.I, ui);
}



inline void
shift_right (bigint & c, const bigint & a, long ui)
{
	if (ui < 0)
		lidia_error_handler("bigint", "shift_right()::index is negative.");
	mpz_tdiv_q_2exp(&c.I, &a.I, ui);
}



inline void
bitwise_and (bigint & c, const bigint & a, const bigint & b)
{
	mpz_and(&c.I, &a.I, &b.I);
}



inline void
bitwise_or (bigint & c, const bigint & a, const bigint & b)
{
	mpz_ior(&c.I, &a.I, &b.I);
}



inline void
bitwise_xor (bigint & c, const bigint & a, const bigint & b)
{
	mpz_xor(&c.I, &a.I, &b.I);
}



inline void
bitwise_not (bigint & b, const bigint & a)
{
	mpz_com(&b.I, &a.I);
}



//
// gcd's
//

inline bigint
gcd (const bigint & a, const bigint & b)
{
	bigint c;

	mpz_gcd(&c.I, &a.I, &b.I);
	return c;
}



inline bigint
bgcd (const bigint & a, const bigint & b)
{
	bigint c;

	mpz_gcd(&c.I, &a.I, &b.I);
	return c;
}



inline bigint
dgcd (const bigint & a, const bigint & b)
{
	bigint c;

	mpz_gcd(&c.I, &a.I, &b.I);
	return c;
}



inline bigint
xgcd (bigint & u, bigint & v, const bigint & a, const bigint & b)
{
	bigint c;

	mpz_gcdext(&c.I, &u.I, &v.I, &a.I, &b.I);
	return c;
}



inline bigint
xgcd_left (bigint & u, const bigint & a, const bigint & b)
{
	bigint c;

	mpz_gcdext(&c.I, &u.I, 0, &a.I, &b.I);
	return c;
}



inline bigint
xgcd_right (bigint & v, const bigint & a, const bigint & b)
{
	bigint c;

	mpz_gcdext(&c.I, &v.I, 0, &b.I, &a.I);
	return c;
}



inline bigint
lcm (const bigint & a, const bigint & b)
{
	bigint c;

	mpz_lcm(&c.I, &a.I, &b.I);
	return c;
}



//
// conversion
//


inline int
bigint_to_string (const bigint & a, char *s)
{
	mpz_get_str(s, 10, &a.I);
	return std::strlen(s);
}



#ifdef C_STDIO

//
// using fread/fwrite
//

inline void
bigint::read_from_file (FILE * fp)
{
	mpz_inp_raw(&I, fp);
}



inline void
bigint::write_to_file (FILE * fp)
{
	mpz_out_raw(fp, &I);
}



//
// using fscanf/fprintf
//

inline void
bigint::scan_from_file (FILE * fp)
{
	mpz_inp_str(&I, fp, 10);
}



inline void
bigint::print_to_file (FILE * fp)
{
	mpz_out_str(fp, 10, &I);
}



#endif	// C_STDIO
