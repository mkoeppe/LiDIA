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
//	Author	: Thomas Pfahler (TPf), Thorsten Rottschaefer (TR),
//		  Victor Shoup (VS)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/fft_prime.h"
#include	"LiDIA/random_generator.h"
#include	"LiDIA/finite_fields/fft_mul_mod.inl"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// The file "fft_mul_mod.inl" either defines or undefines LIDIA_MUL_MOD_SPECIAL.
// The decision is based on a small test programm that is run beforehand
// to determine the faster version.
//
// If LIDIA_MUL_MOD_SPECIAL is defined, two more functions are
// defined (otherwise, multiply_mod is used):
//
// long mul_mod_sp(long a, long b, long q, double qinv)
//    returns a*b mod q
//    qinv must be  1/double(q)
//
// long mul_mod_sp2(long a, long b, long q, double bqinv)
//    returns a*b mod q
//    bqinv must be  double(b)/double(q)
//
// Since their implementation is based on double-arithmetic,
// long-arguments must not be greater than 2^26 !!


#ifdef LIDIA_MUL_MOD_SPECIAL

#define MUL_MOD_MAX ((1 << 26) - 1)

static inline
long mul_mod_sp(long a, long b, long p, double pinv)
{
	double ab = static_cast<double>(a) * static_cast<double>(b);
	register long q = static_cast<long>(ab * pinv);
	register long res = static_cast<long>(ab - (static_cast<double>(q) * static_cast<double>(p)));
//    res += (res >> ((SIZEOF_LONG*8)-1)) & p;
//    res -= p;
//    res += (res >> ((SIZEOF_LONG*8)-1)) & p;

	if (res >= p)
		res -= p;
	else if (res < 0)
		res += p;
	return
		res;
}



static inline
long mul_mod_sp2(long a, long b, long p, double bpinv)
{
	double ab = static_cast<double>(a) * static_cast<double>(b);
	register long q = static_cast<long>(a* bpinv);
	register long res = static_cast<long>(ab - (static_cast<double>(q) * static_cast<double>(p)));
//    res += (res >> ((SIZEOF_LONG*8)-1)) & p;
//    res -= p;
//    res += (res >> ((SIZEOF_LONG*8)-1)) & p;

	if (res >= p)
		res -= p;
	else if (res < 0)
		res += p;
	return
		res;
}



// We could still save another 5-10% for an FFT multiplication if we could omit
// the "if (p < MUL_MOD_MAX)" statements and the code overhead/branches/etc...
// involved.
#endif // LIDIA_MUL_MOD_SPECIAL



//--------------------------------------------------------------------------




bit_reverse_table fft_prime::rev_table;

fft_prime_t* fft_prime::A = 0;
fft_prime_t* fft_prime::B = 0;
fft_prime_t* fft_prime::C = 0;

lidia_size_t fft_prime::size_A = 0;
lidia_size_t fft_prime::size_B = 0;
lidia_size_t fft_prime::size_C = 0;



fft_prime::fft_prime()
{
	debug_handler("fft_prime", "constructor");
	p = 0;
	current_degree = -1;
	RootInvTable = 0;
	TwoInvTable = 0;
	RootTable = 0;
}



fft_prime::fft_prime(const fft_prime &a)
{
	debug_handler("fft_prime", "copy constructor");
	p = 0;
	current_degree = -1;
	RootInvTable = 0;
	TwoInvTable = 0;
	RootTable = 0;
	assign(a);
}



fft_prime::~fft_prime ()      // deletes the fields
{
	debug_handler("fft_prime", "destructor");
	delete [] RootInvTable;
	delete [] TwoInvTable;
	delete [] RootTable;
}



fft_prime & fft_prime::operator = (const fft_prime &a)
{
	debug_handler("fft_prime", "operator = (fft_prime&)");
	assign(a);
	return *this;
}



void fft_prime::assign(const fft_prime &a)
{
	debug_handler("fft_prime", "assign(fft_prime&)");

	delete[] RootTable;
	delete[] RootInvTable;
	delete[] TwoInvTable;
	p = a.p;
	current_degree = a.current_degree;
	max_degree = a.max_degree;
	if (current_degree >= 0) {
		RootTable = new fft_prime_t[current_degree + 1];
		RootInvTable = new fft_prime_t[current_degree + 1];
		TwoInvTable = new fft_prime_t[current_degree + 1];
		memory_handler(RootTable, "fft_prime", "operator = (fft_prime)::"
			       "Error in memory allocation (RootTable)");
		memory_handler(RootInvTable, "fft_prime", "operator = (fft_prime)::"
			       "Error in memory allocation (RootInvTable)");
		memory_handler(TwoInvTable, "fft_prime", "operator = (fft_prime)::"
			       "Error in memory allocation (TwoInvTable)");
		lidia_size_t i;
		for (i = 0; i <= current_degree; i++) {
			RootTable[i] = a.RootTable[i];
			RootInvTable[i] = a.RootInvTable[i];
			TwoInvTable[i] = a.TwoInvTable[i];
		}
	}
	else {
		RootTable = 0;
		RootInvTable = 0;
		TwoInvTable = 0;
	}
}



bool fft_prime::build_tables (lidia_size_t l)
{
	debug_handler("fft_prime", "build_tables (lidia_size_t)");
	bool result = true;
	fft_prime_t w, t;


	// p == 0 ?
	if (p == 0) {
		result = false;
		lidia_error_handler("fft_prime",
				    "build_tables(lidia_size_t)::no prime number set.");
	}

	// If l <= current_degree, do nothing but return true.
	else if (l <= current_degree) {
		result = true;
	}
	// check if a primitive root exists
	// it exists, if and only if 2^l | p-1
	else if (l > max_degree) {
		result = false; // stored prime is not high enough for this degree!
	}
	// tables could be computed...
	else {
		// free pointers and reinitialize them.
		delete[] RootTable;
		delete[] RootInvTable;
		delete[] TwoInvTable;

		RootTable = new fft_prime_t[l+1];
		RootInvTable = new fft_prime_t[l+1];
		TwoInvTable = new fft_prime_t[l+1];

		memory_handler(RootTable, "fft_prime", "build_tables(lidia_size_t)::"
			       "Error in memory allocation (RootTable)");
		memory_handler(RootInvTable, "fft_prime", "build_tables(lidia_size_t)::"
			       "Error in memory allocation (RootInvTable)");
		memory_handler(TwoInvTable, "fft_prime", "build_tables(lidia_size_t)::"
			       "Error in memory allocation (TwoInvTable)");

		random_generator rg;

		// Set current_degree = l.
		current_degree = l;

		//std::cout << "look for primitive root of unity:";
		do {
			//look for primitive root of unity
			rg >> w;
		} while (power_mod(w, p >> 1, p) != static_cast<udigit>(p)-1); // is equal to: if ((w^((p-1)/2))%p != p-1)

		//std::cout << "\nprimitive root of unity:\n" << w <<"\n";

		w = power_mod(w, p >> current_degree, p); // is equal to: w^((p-1)/2^l) % p
		//w now is a primitive 2^current_degree root of unity
		//std::cout << "\nprimitive " << (1<<current_degree) <<"-th root of unity:\n" << w << "\n";

		t = invert_mod(2, p);

		//compute the roottables...
		RootTable[current_degree] = w;
		RootInvTable[current_degree] = invert_mod(w, p);
		TwoInvTable[0] = 1;

		lidia_size_t j;
		t = invert_mod(static_cast<udigit>(2), p);

		//compute the roottables...

		for (j = current_degree-1; j >= 0; j--)
			RootTable[j] = multiply_mod(RootTable[j+1], RootTable[j+1], p);
		for (j = current_degree-1; j >= 0; j--)
			RootInvTable[j] = multiply_mod(RootInvTable[j+1], RootInvTable[j+1], p);
		for (j = 1; j <= current_degree; j++)
			TwoInvTable[j] = multiply_mod(TwoInvTable[j-1], t, p);

	}
	return result;
}



#if 0
//
//    OLD VERSION of FFT
//
//    -- BEGIN --
//
void fft_prime::FFT(fft_prime_t* AA, const fft_prime_t* a, lidia_size_t k, fft_prime_t q, const fft_prime_t* root)
	// Performs a 2^k-point convolution mod q.
	// Assumes, that a has size 2^k
	// and root[i] is a primitve 2^{i}-th root of unity, 1 <= i <= k.
{
	debug_handler("fft_table", "FFT(fft_prime_t*, fft_prime_t*, lidia_size_t, fft_prime_t, fft_prime_t*)");

	if (k == 0) {
		AA[0] = a[0];
		return;
	}

	if (k == 1) {
		AA[0] = add_mod(a[0], a[1], q);
		AA[1] = subtract_mod(a[0], a[1], q);
		return;
	}

	// assume k > 1

	lidia_size_t n = 1 << k;
	lidia_size_t s, m, m2, j;
	fft_prime_t t, u, v, w, z, tt;
	fft_prime_t *p1, *p2, *ub, *ub1;

	// create object of bit_reverse table and copy a in bitreverse-order into AA
	rev_table.copy(AA, a, k);

	ub = AA+n;

	p2 = AA;
	// perform lowest level of evaluating the "bitreverse"-vector
	while (p2 < ub) {
		// s = 1
		u = *p2;
		v = *(p2+1);
		*p2 = add_mod(u, v, q);
		*(p2+1) = subtract_mod(u, v, q);
		p2 += 2;
    	}

	// perform higher levels
	for (s = 2; s < k; s++) {
		m = 1 << s;
		m2 = m >> 1;

		p2 = AA;
		p1 = p2 + m2;
		while (p2 < ub) {
			u = *p2;
			v = *p1;
			*p2 = add_mod(u, v, q);
			*p1 = subtract_mod(u, v, q);
			p1 += m;
			p2 += m;
		}

		z = root[s];
		w = z;
		for (j = 1; j < m2; j++) {
			p2 = AA + j;
			p1 = p2 + m2;
			ub1 = ub-m;
			u = *p2;
			t = multiply_mod(*p1, w, q);

			while (p2 < ub1) {
				tt = multiply_mod(*(p1+m), w, q);
				(*p2) = add_mod(u, t, q);
				(*p1) = subtract_mod(u, t, q);
				p1 += m;
				p2 += m;
				u = *p2;
				t = tt;
			}
			(*p2) = add_mod(u, t, q);
			(*p1) = subtract_mod(u, t, q);
			w = multiply_mod(z, w, q);
		}
	}

	m2 = n >> 1;
	z = root[k];
	w = 1;
	p2 = AA;
	p1 = AA + m2;
	m2--;
	u = *p2;
	t = *p1;
	while (m2) {
		w = multiply_mod(w, z, q);
		tt = multiply_mod(*(p1+1), w, q);
		(*p2) = add_mod(u, t, q);
		(*p1) = subtract_mod(u, t, q);
		p2++;
		p1++;
		u = *p2;
		t = tt;
		m2--;
	}
	(*p2) = add_mod(u, t, q);
	(*p1) = subtract_mod(u, t, q);
}



//
//    -- END --
//
//    OLD VERSION of FFT
//
#endif



//
//    NEW VERSION of FFT
//    if LIDIA_MUL_MOD_SPECIAL is defined, we use a faster FFT implementation
//    for primes less than MUL_MOD_MAX
//
void fft_prime::FFT(fft_prime_t* AA, const fft_prime_t* a, lidia_size_t k, fft_prime_t q, const fft_prime_t* root)
	// Performs a 2^k-point convolution mod q.
	// Assumes, that a has size 2^k
	// and root[i] is a primitve 2^{i}-th root of unity, 1 <= i <= k.
{
	debug_handler("fft_table", "FFT(fft_prime_t*, fft_prime_t*, lidia_size_t, fft_prime_t, fft_prime_t*)");

	if (k == 0) {
		AA[0] = a[0];
		return;
	}

	if (k == 1) {
		AA[0] = add_mod(a[0], a[1], q);
		AA[1] = subtract_mod(a[0], a[1], q);
		return;
	}

	// assume k > 1

	lidia_size_t n = 1 << k;
	lidia_size_t s, m, m2, j;
	fft_prime_t t, u, v, w, z, tt;
	fft_prime_t *p1, *p2, *ub, *ub1;

	// create object of bit_reverse table and copy a in bitreverse-order into AA
	rev_table.copy(AA, a, k);

	ub = AA+n;

	p2 = AA;
	// perform lowest level of evaluating the "bitreverse"-vector
	while (p2 < ub) {
		// s = 1
		u = *p2;
		v = *(p2+1);
		*p2 = add_mod(u, v, q);
		*(p2+1) = subtract_mod(u, v, q);
		p2 += 2;
	}

#ifdef LIDIA_MUL_MOD_SPECIAL
	if (q < MUL_MOD_MAX) {
		double qinv = 1.0 / static_cast<double>(q);
		double wqinv, zqinv;

		// perform higher levels
		for (s = 2; s < k; s++) {
			m = 1 << s;
			m2 = m >> 1;

			p2 = AA;
			p1 = p2 + m2;
			while (p2 < ub) {
				u = *p2;
				v = *p1;
				*p2 = add_mod(u, v, q);
				*p1 = subtract_mod(u, v, q);
				p1 += m;
				p2 += m;
			}

			z = root[s];
			w = z;
			for (j = 1; j < m2; j++) {
				wqinv = (static_cast<double>(w))*qinv;
				p2 = AA + j;
				p1 = p2 + m2;
				ub1 = ub-m;
				u = *p2;
				t = mul_mod_sp2(*p1, w, q, wqinv);

				while (p2 < ub1) {
					tt = mul_mod_sp2(*(p1+m), w, q, wqinv);
					(*p2) = add_mod(u, t, q);
					(*p1) = subtract_mod(u, t, q);
					p1 += m;
					p2 += m;
					u = *p2;
					t = tt;
				}
				(*p2) = add_mod(u, t, q);
				(*p1) = subtract_mod(u, t, q);
				w = mul_mod_sp2(z, w, q, wqinv);
			}
		}

		m2 = n >> 1;
		z = root[k];
		zqinv = (static_cast<double>(z))*qinv;
		w = 1;
		p2 = AA;
		p1 = AA + m2;
		m2--;
		u = *p2;
		t = *p1;
		while (m2) {
			w = mul_mod_sp2(w, z, q, zqinv);
			tt = mul_mod_sp(*(p1+1), w, q, qinv);
			(*p2) = add_mod(u, t, q);
			(*p1) = subtract_mod(u, t, q);
			p2++;
			p1++;
			u = *p2;
			t = tt;
			m2--;
		}
		(*p2) = add_mod(u, t, q);
		(*p1) = subtract_mod(u, t, q);
	}
	else
#endif
	{
		for (s = 2; s < k; s++) {
			m = 1 << s;
			m2 = m >> 1;

			p2 = AA;
			p1 = p2 + m2;
			while (p2 < ub) {
				u = *p2;
				v = *p1;
				*p2 = add_mod(u, v, q);
				*p1 = subtract_mod(u, v, q);
				p1 += m;
				p2 += m;
			}

			z = root[s];
			w = z;
			for (j = 1; j < m2; j++) {
				p2 = AA + j;
				p1 = p2 + m2;
				ub1 = ub-m;
				u = *p2;
				t = multiply_mod(*p1, w, q);

				while (p2 < ub1) {
					tt = multiply_mod(*(p1+m), w, q);
					(*p2) = add_mod(u, t, q);
					(*p1) = subtract_mod(u, t, q);
					p1 += m;
					p2 += m;
					u = *p2;
					t = tt;
				}
				(*p2) = add_mod(u, t, q);
				(*p1) = subtract_mod(u, t, q);
				w = multiply_mod(z, w, q);
			}
		}

		m2 = n >> 1;
		z = root[k];
		w = 1;
		p2 = AA;
		p1 = AA + m2;
		m2--;
		u = *p2;
		t = *p1;
		while (m2) {
			w = multiply_mod(w, z, q);
			tt = multiply_mod(*(p1+1), w, q);
			(*p2) = add_mod(u, t, q);
			(*p1) = subtract_mod(u, t, q);
			p2++;
			p1++;
			u = *p2;
			t = tt;
			m2--;
		}
		(*p2) = add_mod(u, t, q);
		(*p1) = subtract_mod(u, t, q);
	}
}



bool fft_prime::evaluate(fft_prime_t* o, int type, const void* a, lidia_size_t d, lidia_size_t k)
{
        // Performs a 2^k (point,value) representation mod p.
        // Assumes, that "a" represents a degree d polynomial modulo p.
        //
        // if p  = 0, lidia_error_handler is called
	// if k <= current_degree ok = true; else ok = build_tables (k).
	// if (ok) fills "a" with zeroes to get size 2^k, if necessary.
	// if (ok) perform the convolution by calling FFT
	// return ok.

	debug_handler("fft_prime" ,
		      "evaluate(fft_prime_t*, const fft_prime_t*, lidia_size_t, lidia_size_t)");

	bool ok;
	lidia_size_t i;

	// if p  = 0, lidia_error_handler is called
	if (p == 0)
		lidia_error_handler("fft_prime",
				    "build_tables(lidia_size_t)::fft_prime equal zero!");

	if (k <= current_degree)
		ok = true; // no need to build new tables..
	else
		ok = build_tables(k);

	if (ok) {
		if (!type && (d+1) == (1 << k)) {
			FFT(o, static_cast<const fft_prime_t*>(a), k, p, RootTable);
			return ok;
		}

		// check, wether B is large enough
		if (size_B < (1 << k)) {
			delete [] B;
			size_B = 1 << k;
			B = new fft_prime_t[size_B];
			memory_handler(B, "fft_prime",
				       "evaluate(fft_prime_t* , const fft_prime_t* , lidia_size_t , lidia_size_t)::"
				       "Error in memory allocation (B)");
		}

		// if (ok) fills "empty" places in "a" with zeroes
		// to get size 2^k, if necessary.
		if (type) {
			for (i = 0; i <= d; i++)
				B[i] = (static_cast<const udigit_mod*>(a))[i].get_mantissa(); // copys a into In

			for (i = d+1; i < size_B; i++)
				B[i] = 0; // fills rest of In with 0's
		}
		else {
			for (i = 0; i <= d; i++)
				B[i] = (static_cast<const fft_prime_t*>(a))[i]; // copys a into A

			for (i = d+1; i < size_B; i++)
				B[i] = 0; // fills rest of In with 0's
		}
		//perform the convolution into point-value form by calling FFT
		FFT(o, B, k, p, RootTable);

	}
	return ok;
}



bool fft_prime::interpolate(fft_prime_t* o, const fft_prime_t* a, lidia_size_t k)
{
        // Performs a 2^k coefficient representation mod p.
	// Assumes, that "a" represents a 2^k (point,value) representation
	// modulo p.
	//
	// if p  = 0, lidia_error_handler is called
	// if k <= current_degree ok = true; else ok = build_tables (k).
	// if (ok) fills "a" with zeroes to get size 2^k, if necessary.
	// if (ok) perform the convolution by calling FFT
	// return ok.
	//
        // This function assumes, that the user exactly knows what he is doing;
        // i.e., he must not change the tables of the fft_prime by evaluating
        // another polynomial with a larger current_degree (other root of unity),
        // before calling the interpolate routine for the first polynomial.

	debug_handler("fft_prime" , "interpolate(fft_prime_t*, const fft_prime_t*, lidia_size_t, lidia_size_t)");
	lidia_size_t j;
	bool ok;

	if (k <= max_degree)
		ok = true;
	else
		ok = build_tables (k);

	if (ok) {

		//perform the convolution back to polynom by calling FFT with the "inverse"-teble
		FFT(o, a, k, p, RootInvTable);

		fft_prime_t two_inv = TwoInvTable[k];
#ifdef LIDIA_MUL_MOD_SPECIAL
		if (p < MUL_MOD_MAX) {
			double pinv = 1.0 / static_cast<double>(p);
			fft_prime_t *vp = o;
			for (j = (1 << k); j != 0; j--, vp++)
				*vp = mul_mod_sp(*vp, two_inv, p, pinv);
		}
		else
#endif
			for (j = 0; j < (1 << k); j++)
				o[j] = multiply_mod(o[j], two_inv, p);
	}

	return ok;
}



bool fft_prime::interpolate2(fft_prime_t* o, const fft_prime_t* a, lidia_size_t k, lidia_size_t lo, lidia_size_t length)
{
	debug_handler("fft_prime" , "interpolate2(fft_prime_t*, const fft_prime_t*, lidia_size_t, lidia_size_t)");
	lidia_size_t j;

	//perform the convolution back to polynom by calling FFT with the "inverse"-teble
	FFT(o, a, k, p, RootInvTable);

	fft_prime_t two_inv = TwoInvTable[k];
	fft_prime_t *vp = o+lo;
#ifdef LIDIA_MUL_MOD_SPECIAL
	if (p < MUL_MOD_MAX) {
		double pinv = 1.0 / static_cast<double>(p);
		for (j = length; j != 0; j--, vp++)
			*vp = mul_mod_sp(*vp, two_inv, p, pinv);
	}
	else
#endif
		for (j = length; j != 0; j--, vp++)
			*vp = multiply_mod(*vp, two_inv, p);

	return true;
}



void fft_prime::pointwise_multiply(fft_prime_t* x, const fft_prime_t* a, const fft_prime_t* b, lidia_size_t k) const
// x[i] = a[i]*b[i] mod p; i = 0 .. (2^k-1)i
{
	debug_handler("fft_prime", "pointwise_multiply(fft_prime_t*, "
		      "const fft_prime_t*, const fft_prime_t*, lidia_size_t)");

	lidia_size_t K = 1 << k;
	lidia_size_t i;
	fft_prime_t *xp = x;
	const fft_prime_t *ap = a, *bp = b;

#ifdef LIDIA_MUL_MOD_SPECIAL
	if (p < MUL_MOD_MAX) {
		double pinv = 1.0 / static_cast<double>(p);
		for (i = K; i != 0; i--, xp++, ap++, bp++)
			*xp = mul_mod_sp(*ap, *bp, p, pinv);
	}
	else
#endif
		for (i = K; i != 0; i--, xp++, ap++, bp++)
			*xp = multiply_mod(*ap, *bp, p);
}



void fft_prime::pointwise_add(fft_prime_t* x, const fft_prime_t* a,
			      const fft_prime_t* b, lidia_size_t k) const
// x[i] = a[i]+b[i] mod p; i = 0 .. (2^k-1)
{
	debug_handler("fft_prime", "pointwise_add(fft_prime_t*, "
		      "const fft_prime_t*, const fft_prime_t*, lidia_size_t)");

	lidia_size_t K = 1 << k;
	lidia_size_t i;
	fft_prime_t *xp = x;
	const fft_prime_t *ap = a, *bp = b;

	for (i = K; i != 0; i--, xp++, ap++, bp++)
		*xp = add_mod(*ap, *bp, p);
}



void fft_prime::pointwise_subtract(fft_prime_t* x, const fft_prime_t* a,
				   const fft_prime_t* b, lidia_size_t k) const
// x[i] = a[i]-b[i] mod p; i = 0 .. (2^k-1)
{
	debug_handler("fft_prime", "pointwise_subtract(fft_prime_t*, "
		      "const fft_prime_t*, const fft_prime_t*, lidia_size_t)");

	lidia_size_t K = 1 << k;
	lidia_size_t i;
	fft_prime_t *xp = x;
	const fft_prime_t *ap = a, *bp = b;

	for (i = K; i != 0; i--, xp++, ap++, bp++)
		*xp = subtract_mod(*ap, *bp, p);
}



void fft_prime::set_prime (fft_prime_t q)
{
	debug_handler("fft_prime", "set_prime(fft_prime_t)");

	if (p != q) {
		current_degree = -1;

		delete[] RootTable;
		delete[] RootInvTable;
		delete[] TwoInvTable;

		RootInvTable = 0;
		TwoInvTable = 0;
		RootTable = 0;

		p = q; // sets actual modulus p = q

		max_degree = 0;
		q--; // comutes the maximum possible degree for the given
		// prime q and stores it in max_degree
		while ((q & 1) == 0) {
			q >>= 1;
			max_degree++;
		}
	}
}  // assigns p = q and reinitializes the tables
// if p = q do nothing, because q is the old p!!


lidia_size_t fft_prime::get_max_degree ()
{
	debug_handler("fft_prime", "get_max_degree()");
	if (p == 0)
		lidia_error_handler("fft_prime::get_max_degree", "stored prime is 0 =  > there's no actual max_degree");
	return max_degree;
}



bool multiply_fft (void* x, int type,
		   const void* const a, lidia_size_t da,
		   const void* const b, lidia_size_t db,
		   fft_prime & fftprime)

{
	debug_handler("fft_prime",
		      "multiply_fft(fft_prime_t*, udigit_mod*, lidia_size_t, udigit_mod*, lidia_size_t, &fftprime)");

	bool a_is_equal_b = true;
	bool result = true;
	lidia_size_t i, d, k;

	if (da < 0 || db < 0) {
		lidia_error_handler ("fft_prime::multiply_fft",
				     "Invalid negative degree.");
	}
	else {
		// check, if a=b (assuming equality)
		if (da == db && a == b)
			a_is_equal_b = true;
		else
			a_is_equal_b = false;

		// evaluate degree of a*b
		d = da + db + 1;

		// evaluate for which k 2^k >= d
		// (means the next power of 2 bigger than d)
		i = 1; k = 0;
		while (i < d) {
			i = i << 1;
			k++;
		}
		//std::cout << i;

		if (a_is_equal_b) {
			//call fkt. evaluate only once because a = b
			// check, wether A is large enough
			if (fft_prime::size_A < (1 << k)) {
				delete [] fft_prime::A;
				fft_prime::size_A = 1 << k;
				fft_prime::A = new fft_prime_t[fft_prime::size_A];
				memory_handler(A, "fft_prime", "multiply_fft" "Error in memory allocation (A)");
			}
			// check, wether C is large enough
			if (fft_prime::size_C < (1 << k)) {
				delete [] fft_prime::C;
				fft_prime::size_C = 1 << k;
				fft_prime::C = new fft_prime_t[fft_prime::size_C];
				memory_handler(C, "fft_prime", "multiply_fft" "Error in memory allocation (C)");
			}

			if (!fftprime.evaluate(fft_prime::A, type, a, da, k)) {
				lidia_error_handler("fft_prime::evaluate", "tables could not be build");
				result = false;
			}
			else {
				fftprime.pointwise_multiply(fft_prime::A, fft_prime::A, fft_prime::A, k);
				fftprime.interpolate(fft_prime::C, fft_prime::A, k);
			}
		}
		else {
			// check, wether A is large enough
			if (fft_prime::size_A < (1 << k)) {
				delete [] fft_prime::A;
				fft_prime::size_A = 1 << k;
				fft_prime::A = new fft_prime_t[fft_prime::size_A];
				memory_handler(A, "fft_prime", "multiply_fft" "Error in memory allocation (A)");
			}

			// check, wether C is large enough
			if (fft_prime::size_C < (1 << k)) {
				delete [] fft_prime::C;
				fft_prime::size_C = 1 << k;
				fft_prime::C = new fft_prime_t[fft_prime::size_C];
				memory_handler(C, "fft_prime", "multiply_fft" "Error in memory allocation (C)");
			}

			if (!fftprime.evaluate(fft_prime::A, type, a, da, k) || !fftprime.evaluate(fft_prime::C, type, b, db, k)) {
				lidia_error_handler("fft_prime::evaluate", "tables could not be build");
				result = false;
			}
			else {
				fftprime.pointwise_multiply(fft_prime::A, fft_prime::A, fft_prime::C, k);
				fftprime.interpolate(fft_prime::C, fft_prime::A, k);
			}
		}

		if (result) // copy contents of C back to vector of type udigit_mog_p*
			if (type)
				for (i = 0; i < (1 << k); i++)	   (static_cast<udigit_mod*>(x))[i] = fft_prime::C[i];
			else
				for (i = 0; i < (1 << k); i++)	   (static_cast<fft_prime_t*>(x))[i] = fft_prime::C[i];


	}
	// end else if (da < 0 || db < 0)


	return result; // return true, if fft was succesful else false
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
