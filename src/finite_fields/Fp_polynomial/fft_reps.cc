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
//	Author	: Thomas Pfahler (TPf)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/Fp_polynomial.h"
#include	"LiDIA/finite_fields/Fp_polynomial_fft.h"
#include	"LiDIA/crt.h"
#include	"LiDIA/fft_prime.h"
#include	"LiDIA/udigit.h"
#include	"LiDIA/finite_fields/fft_mul_mod.inl"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



// The file "fft_mul_mod.inl" either defines or undefines LIDIA_MUL_MOD_SPECIAL.
// The decision is based on a small test programm that is run beforehand
// to determine the faster version (see simple_classes/fft_prime).
//
// If LIDIA_MUL_MOD_SPECIAL is defined, we use fft_primes < 2^26
// (see fft_prime.c, LIDIA_MUL_MOD_SPECIAL for reasons). This is faster
// than using larger primes (despite the fact that we  might need more of these
// primes).

#ifdef LIDIA_MUL_MOD_SPECIAL
#define MUL_MOD_MAX ((1 << 26) - 1)
#else
#define MUL_MOD_MAX (comparator< udigit >::min(max_udigit(), (1L << (SIZEOF_LONG*CHAR_BIT-1)) - 1))
//		 IMPORTANT: long(MUL_MOD_MAX) must be positive !
#endif



// ##############################################################
// TODO:
// - if num_primes==1, skip combine/get_result (done for fft_rep, but not for
//	modular_fft_rep)
// - do profiling: count number of calls for the functions, inline them
// 	accordingly
// - memory management: new/delete[] of fft_prime's



//
// This code heavily relies on the fact that class crt_table
// DOES NOT MODIFY THE ORDER of the primes when it is given
// a vector of primes for initialization !
//
// This is because we first compute a vector of suitable FFT-primes
// (let's say of lenght n) and then initialize a crt_table with
// the first m < n FFT-primes.
// Let ct be of type crt_table and fd be of type fft_data.
// When multiplying two polynomials, we first
// 	reduce the coefficients modulo ct.get_prime(i), then
// 	apply the (modular) FFT using fd.node->prime[i] and eventually
// 	combine the result -- assuming that has been computed modulo
// 						ct.get_prime(i)
// where (i = 0, ..., m-1).
//
// As a consequence, if crt_table (or class crt) ever "rearranges" the
// primes, something goes wrong ...
//



static
fft_prime_t next_prime(fft_prime_t old, lidia_size_t max_degree)
	// returns q where q is less than (old % 2^l) + 1 and prime,
	// q = 1 mod 2^max_degree
	// error if no more primes are found.
{
	debug_handler("fft_data", "next_prime(...)");
	fft_prime_t cand;
	fft_prime_t N = 1 << max_degree;

	if (old == fft_prime_t(~0)) {
		cand = MUL_MOD_MAX;
		if (max_degree != 0)		// make cand = 1 mod 2^max_degree
		{
			cand >>= max_degree;
			cand <<= max_degree;
			cand++; // (*)
		}
	}
	else {
		cand = old;
		if (cand % N != 1) {
			lidia_error_handler("fft_data", "next_prime(...) :: "
					    "old != 1 mod 2^max_degree");
			return 0; // LC
		}
	}
	for (;;) {
		if (cand <= N) {
			lidia_error_handler("fft_data", "next_prime(...) :: "
					    "no more primes");
			return 0; // LC
		}
		cand -= N;
		if (!is_prime(bigint(cand), 8))
			continue;
		return  cand;
	}
}



//*************************************************************************
//				class fft_data
//*************************************************************************

fft_data::item* fft_data::head = 0;

fft_data::item::~item()
{
	delete[] prime;
	fft_data::crt_item * cptr = crt_list;
	fft_data::crt_item * cptr2;
	while (cptr != 0) {
		cptr2 = cptr->next; //remove all entries from crt_list
		delete cptr->ct;
		delete cptr;
		cptr = cptr2;
	}
}



fft_data::fft_data()
	: node(0),
	  crt_node(0),
	  k(-1)
{ }



fft_data::fft_data(const bigint &p, lidia_size_t l)
	: node(0),
	  crt_node(0),
	  k(-1)
{
	init(p, l);
}



fft_data::~fft_data()
{
	release();
}



fft_data::fft_data(const fft_data& F)
	: node(F.node),
	  crt_node(F.crt_node),
	  k(F.k)
{
	F.add_ref();
}



void fft_data::release()
{
	if (node) {
		if (--(node->ref_counter) <= 0 && //if no references to node
		    node != fft_data::head)		//and node != head
		{
			fft_data::item * ptr = fft_data::head;
			while (ptr->next != node)
				ptr = ptr->next;
			if (ptr->next == 0) {
				lidia_error_handler("fft_data", "~fft_data::internal error");
				return;
			}

			ptr->next = node->next; //remove node from list

			delete node; //finally, delete node
		}
		node = 0;
	}
	crt_node = 0;
	k = -1;
}




fft_data & fft_data::operator = (const fft_data & F)
{
	release();
	node = F.node;
	crt_node = F.crt_node;
	k = F.k;
	F.add_ref();
	return *this;
}



lidia_size_t fft_data::number_of_primes() const
{
	if (crt_node == 0) {
		lidia_error_handler("fft_data", "number_of_primes() :: "
				    "fft_data not initialized");
		return 0; // LC
	}
	return crt_node->ct->number_of_primes();
}



void fft_data::init(const bigint& p, lidia_size_t l)
{
	if (node != 0) {
		lidia_error_handler("fft_data", "init(...) :: "
				    "re-initialization not allowed");
		return;
//alternatively:	release(); (anything else?)
	}


	//
	//section 1: find a suitable node
	//
	if (fft_data::head != 0)			//if the list is not empty
	{
		if (fft_data::head->max_degree >= l)	//if *head is suff. 'large'
		{
			node = fft_data::head;
			if (node->next == 0) {
				if (node->ref_counter <= 0 && node->max_degree > l+3)
					node->max_degree = l; //we can safely decrease
				//head->max_degree (bec. noone
				//currently uses head)
			}
			else {
				if ((node->next->max_degree >= l)
				    && (node->next->num_primes > node->num_primes))
					node = node->next; //if head->next is still large
				//enough and has more primes,
				//use it
			}
		}
	}

	if (node == 0) {
		//we need a new item

		fft_data::item *ptr = new fft_data::item;
		memory_handler(ptr, "fft_data", "init(...) :: "
			       "Error in memory allocation (ptr)");
		ptr->prime = 0;
		ptr->num_primes = 0;
		ptr->max_degree = l; //set max. degree
		ptr->crt_list = 0;
		ptr->ref_counter = 0;

		if (fft_data::head != 0)
			if (fft_data::head->ref_counter <= 0) {
				delete fft_data::head; //we can delete the old one
				fft_data::head = 0; //since we'll create a new head
			}
		ptr->next = fft_data::head; //and insert *ptr as the first
		fft_data::head = ptr; //element of the list
		node = ptr;

	}

	k = l;
	add_ref();


	//
	// section 2: find a suitable crt_node
	//

	crt_node = 0;
	crt_item * cptr = node->crt_list;
	while (cptr != 0) {
		if (cptr->p == p)			//if there is already a table
			crt_node = cptr; //for p
		cptr = cptr->next;
	}

	if (crt_node == 0) {
		//we need a new crt_item
		cptr = new crt_item;
		memory_handler(cptr, "fft_data", "init(...) :: "
			       "Error in memory allocation (cptr)");

		//the product of all primes must be > 2^max_degree * modulus^2
		//we do not multiply the primes explicitely, but compute logarithms;
		//thus we can simply add bit_lengths instead of multiplying large numbers
		lidia_size_t b_log = node->max_degree + 2*(p.bit_length()-1) + 2;
		// = log_2(2^max_degree * modulus^2) + safety

		lidia_size_t new_num = 0, sum = 0, i, j;
		for (i = 0; i < node->num_primes; i++) {
			sum += integer_log(node->prime[i].get_prime());
			if (sum >= b_log) {
				//no new primes needed
				new_num = i+1;
				break;
			}
		}

		if (new_num == 0) {
			//we need more fft_primes
			fft_prime_t *tmp = 0, *tmp2 = 0;
			long q;
			if (node->num_primes == 0)
				q = fft_prime_t(~0);
			else
				q = node->prime[node->num_primes - 1].get_prime();

			for (; b_log >= sum; new_num++) {
				if ((new_num % 64) == 0)	//need space for another
				{				//64 primes
					tmp2 = new fft_prime_t[new_num + 64];
					memory_handler(tmp2, "fft_data", "init(...) :: "
						       "Error in memory allocation (tmp2)");
					for (j = 0; j < new_num; j++)
						tmp2[j] = tmp[j];
					delete[] tmp;
					tmp = tmp2;
				}
				q = next_prime(q, node->max_degree); //q is prime and < p
				//and q = 1 mod 2^max_degree
				tmp[new_num] = q;
				sum += integer_log(q);
			}

			lidia_size_t old_num = node->num_primes;
			fft_prime *fptr = new fft_prime[old_num + new_num];
			memory_handler(fptr, "fft_data", "init(...) :: "
				       "Error in memory allocation (node->prime)");
			for (j = 0; j < old_num; j++) fptr[j] = node->prime[j];
			for (j = 0; j < new_num; j++) fptr[j + old_num].set_prime(tmp[j]);
			delete[] node->prime;
			delete[] tmp;

			new_num += old_num;
			node->prime = fptr;
			node->num_primes = new_num;
			//std::cout<<"delete "<<old_num<<" primes, insert "<<new_num<<" primes."<<std::endl;

			lidia_size_t new_max = fptr[old_num].get_max_degree();
			for (j = old_num+1; j < new_num; j++)
				if (new_max > fptr[j].get_max_degree())
					new_max = fptr[j].get_max_degree();
			//std::cout<<"new_max = "<<new_max<<"     l = "<<l<<"   max_degree = "<<node->max_degree<<std::endl;
			//std::cout << "ref_counter = " << node->ref_counter << std::endl;

			if (old_num == 0)
				node->max_degree = new_max;
			else
				new_max = node->max_degree;

			for (j = old_num; j < new_num; j++)
				fptr[j].build_tables(new_max);
#if 0
			for (j = 0; j < old_num; j++)
				std::cout << "..prime[" << j << "] = " << fptr[j].get_prime()
					  << ", max_degree = " << fptr[j].get_max_degree() << std::endl;
			for (j = old_num; j < new_num; j++)
				std::cout << "  prime[" << j << "] = " << fptr[j].get_prime()
					  << ", max_degree = " << fptr[j].get_max_degree() << std::endl;
			std::cout << "list of fft_data::item:" << std::endl;
			fft_data::item *ptr2 = fft_data::head;
			while (ptr2 != 0)
			{
				std::cout << "  [ np = " << ptr2->num_primes << ", md = "
					  << ptr2->max_degree << ", rc = " << ptr2->ref_counter
					  << " ]" << std::endl;
				ptr2 = ptr2->next;
			}
#endif
		}

		base_vector< sdigit > vec(new_num, new_num);
		for (i = 0; i < new_num; i++)
			vec[i] = node->prime[i].get_prime();

		cptr->ct = new crt_table(vec);
		memory_handler(cptr->ct, "fft_data", "init(...) :: "
			       "Error in memory allocation (cptr->c)");

		cptr->p = p; //set modulus
		cptr->next = node->crt_list; //insert *cptr as the first
		node->crt_list = cptr; //element in the crt_list

		crt_node = cptr;
	}
}



//*************************************************************************
//			class base_fft_rep
//************************************************************************/

fft_prime_t * base_fft_rep::stat_vec = 0;
lidia_size_t  base_fft_rep::stat_k = -1;



base_fft_rep::base_fft_rep() :
	k(-1),
	c(0)
{}



base_fft_rep::~base_fft_rep()
{
	delete c;
}



void base_fft_rep::init_(const fft_data &f)
{
	//raise an error if f is not initialized, or if fd has already been set
	if (f.node == 0 || fd.node != 0 || c != 0) {
//	lidia_error_handler("base_fft_rep", "init(...) :: internal error");
//	return;
		delete c;
		c = 0;
		k = -1; //TOM
	}
	fd = f;

	c = new crt(*f.crt_node->ct);
	memory_handler(c, "base_fft_rep", "init(fft_data&) :: "
		       "Error in memory allocation (c)");
}



void base_fft_rep::to_mod_rep(fft_prime_t *x, const Fp_polynomial &a,
			      lidia_size_t lo, lidia_size_t hi, lidia_size_t index)
	// computes an n = 2^k point convolution.
	// if deg(a) >= 2^k, then a is first reduced modulo X^n-1.
{

	if (k < 0 || fd.node == 0) {
		lidia_error_handler("base_fft_rep", "to_mod_rep(...)::not initialized");
		return; // LC
	}
	lidia_size_t K = 1 << k;
	lidia_size_t j, m, j1;

	hi = comparator< lidia_size_t >::min(hi, a.degree());
	m = comparator< lidia_size_t >::max(hi-lo + 1, 0);

	const bigint *aptr = &a.coeff[lo];
	bool ok;

//	if (m < K) {
//		c->reduce(stat_vec, aptr, m, index);
	fft_prime_t q = c->get_prime(index);
	if (m < K) {
		for (j = 0; j < m; j++)
			remainder(stat_vec[j], aptr[j], q);
//	ok = fd.node->prime[index].evaluate(x, 0, stat_vec, m-1, k);
// if we fill x with trailing zeros, we don't need copying in
// fft_prime::evaluate(...):
		for (j = m; j < K; j++)
			stat_vec[j] = 0;
		ok = fd.node->prime[index].evaluate(x, 0, stat_vec, K-1, k);
	}
	else {
		if (m < (K << 1)) {
			bigint accum;
			lidia_size_t m2 = m - K;
			const bigint *ap2 = &a.coeff[lo+K];
			for (j = 0; j < m2; j++, aptr++, ap2++) {
				add(accum, *aptr, *ap2);
				remainder(stat_vec[j], accum, q);
			}
			for (j = m2; j < K; j++)
				remainder(stat_vec[j], a.coeff[lo+j], q);
		}
		else {
			bigint accum;
			for (j = 0; j < K; j++) {
				accum.assign(aptr[j]);
				for (j1 = j + K; j1 < m; j1 += K)
					add(accum, accum, aptr[j1]);
#if 0
				c->reduce(stat_vec[j], accum, index);
#else
				remainder(stat_vec[j], accum, q);
#endif
			}
		}
		ok = fd.node->prime[index].evaluate(x, 0, stat_vec, K-1, k);
	}
	if (!ok) {
		std::cout << "prime[" << index << "] = "
			  << fd.node->prime[index].get_prime() << std::endl;
		std::cout << "k = " << k << ", get_max_degree = "
			  << fd.node->prime[index].get_max_degree() << std::endl;
		lidia_error_handler("base_fft_rep", "to_mod_rep(...)::re-init of "
				    "FFT prime");
	}
}



void base_fft_rep::combine(fft_prime_t *x, lidia_size_t l, lidia_size_t index)
{
	//if (c->number_of_primes() > 1)
	c->combine(x, l, index);
}



void base_fft_rep::get_res(bigint *x, lidia_size_t l)
	// assumes that enough space for x has been allocated
{
	bigint *x2 = x;
	c->get_result(x, l);
	if (x2 != x) {
		lidia_error_handler("base_fft_rep", "get_res(...)::"
				    "unnecessary reallocation in crt::get_result");
		return; // LC
	}

	const bigint &p = fd.crt_node->p;
	for (lidia_size_t i = l-1; i >= 0; i--)
		Remainder(x[i], x[i], p);
}



//*************************************************************************
//			class modular_fft_rep
//*************************************************************************

modular_fft_rep::modular_fft_rep() :
	max_k(-1),
	vec(0)
{}



modular_fft_rep::~modular_fft_rep()
{
	delete[] vec;
}



void modular_fft_rep::set_size(lidia_size_t l)
{
	if (fd.node ? (fd.node->max_degree < l) : 1) {
		lidia_error_handler("modular_fft_rep", "set_size(...)::internal error");
		return; // LC
	}
	if (l > max_k && l >= 0) {

//if (k != -1)	//###########
//    std::cout<<"   #in modular_fft_rep::set_size; enlarge vec: "<<k<<" ->"<<l<<std::endl;

		delete[] vec;
		vec = new fft_prime_t[1 << l];
		memory_handler(vec, "modular_fft_rep", "init(...) :: "
			       "Error in memory allocation (vec)");
		max_k = l;

		if (base_fft_rep::stat_k < l) {
			delete[] base_fft_rep::stat_vec;
			base_fft_rep::stat_vec = new fft_prime_t[1 << l];
			memory_handler(base_fft_rep::stat_vec, "modular_fft_rep",
				       "init(...) :: Error in memory allocation (vec)");
			base_fft_rep::stat_k = l;
		}
	}
	k = l;
	c->reset();
}



void modular_fft_rep::from_modular_fft_rep(lidia_size_t lo, lidia_size_t hi,
					   lidia_size_t index)
{
//    std::cout <<"@" << std::flush;
//std::cout<<"from_modular_fft_rep lo="<<lo<<"   hi="<<hi<<"   index="<<index<<std::endl;
	lidia_size_t K = 1 << k;
	hi = comparator< lidia_size_t >::min(hi, K-1);
	lidia_size_t l = comparator< lidia_size_t >::max(hi-lo+1, 0);

	bool ok;

// old code:
//    ok = fd.node->prime[index].interpolate(stat_vec, vec, k);
	ok = fd.node->prime[index].interpolate2(stat_vec, vec, k, lo, l);

	if (!ok)
		lidia_error_handler("modular_fft_rep", "from_modular_fft_rep(...)::"
				    "re-init of FFT prime");

	combine(&stat_vec[lo], l, index);
}



void modular_fft_rep::get_result_ptr(bigint *a, lidia_size_t lo,
				     lidia_size_t hi)
//used only in build_from_roots
//enough space must be allocated for result !!!!!
{
	debug_handler("modular_fft_rep", "get_result_ptr(bigint*, lidia_size_t, lidia_size_t)");

	lidia_size_t l = comparator< lidia_size_t >::max(hi-lo+1, 0);
	get_res(a, l);
}



void modular_fft_rep::get_result(Fp_polynomial &a, lidia_size_t lo,
				 lidia_size_t hi)
{
	debug_handler("modular_fft_rep", "get_result(Fp_polynomial&, lidia_size_t, lidia_size_t)");

	lidia_size_t K = 1 << k;
	hi = comparator< lidia_size_t >::min(hi, K-1);
	lidia_size_t l = comparator< lidia_size_t >::max(hi-lo+1, 0);

	a.set_modulus(fd.crt_node->p);
	a.set_degree(l-1);

	get_res(a.coeff, l);

	a.remove_leading_zeros();
}



void multiply(modular_fft_rep &x, const modular_fft_rep &a,
	      const modular_fft_rep &b, lidia_size_t index)
{
	debug_handler("modular_fft_rep", "multiply(...)");

	lidia_size_t k = a.k;
	if (a.fd != b.fd || a.fd != x.fd || k != b.k || k != x.k) {
		lidia_error_handler("modular_fft_rep",
				    "multiply(...)::FFT rep mismatch");
		return;
	}

	a.fd.node->prime[index].pointwise_multiply(x.vec, a.vec, b.vec, k);
}



void add(modular_fft_rep &x, const modular_fft_rep &a,
	 const modular_fft_rep &b, lidia_size_t index)
{
	debug_handler("modular_fft_rep", "add(...)");

	lidia_size_t k = a.k;
	if (a.fd != b.fd || a.fd != x.fd || k != b.k || k != x.k) {
		lidia_error_handler("modular_fft_rep", "add(...)::FFT rep mismatch");
		return;
	}

	a.fd.node->prime[index].pointwise_add(x.vec, a.vec, b.vec, k);
}



void subtract(modular_fft_rep &x, const modular_fft_rep &a,
	      const modular_fft_rep &b, lidia_size_t index)
{
	debug_handler("modular_fft_rep", "subtract(...)");

	lidia_size_t k = a.k;
	if (a.fd != b.fd || a.fd != x.fd || k != b.k || k != x.k) {
		lidia_error_handler("modular_fft_rep",
				    "subtract(...)::FFT rep mismatch");
		return;
	}

	a.fd.node->prime[index].pointwise_subtract(x.vec, a.vec, b.vec, k);
}



void add_mul(modular_fft_rep &x,
	     const modular_fft_rep &a, const modular_fft_rep &b,
	     const modular_fft_rep &c, const modular_fft_rep &d, lidia_size_t index)
{
	debug_handler("modular_fft_rep", "add_mul(modular_fft_rep&, modular_fft_rep&, modular_fft_rep&, modular_fft_rep&, modular_fft_rep&, lidia_size_t)");

	lidia_size_t k = a.k;
	if (a.fd != b.fd || a.fd != c.fd || a.fd != d.fd || a.fd != x.fd
	    || k != b.k || k != c.k || k != d.k || k != x.k) {
		lidia_error_handler("modular_fft_rep",
				    "add_mul(...)::FFT rep mismatch");
		return;
	}

//    a.fd.node->prime[index].pointwise_add_mul(x.vec, a.vec, b.vec, c.vec, d.vec, k);

	const fft_prime *P = &a.fd.node->prime[index];
	P->pointwise_multiply(base_fft_rep::stat_vec, a.vec, b.vec, k);
	P->pointwise_multiply(x.vec, c.vec, d.vec, k);
	P->pointwise_add(x.vec, x.vec, base_fft_rep::stat_vec, k);
}



//*************************************************************************
//			    class fft_rep
//*************************************************************************

fft_rep::fft_rep() :
	max_k(-1),
	tbl(0),
	num_primes(0)
{ }



fft_rep::~fft_rep()
{
	if (tbl != 0) {
		for (lidia_size_t i = 0; i < num_primes; i++)
			delete[] tbl[i];
		delete[] tbl;
	}
}



void fft_rep::init(const fft_data &f)
{
	if (tbl != 0) {
		for (lidia_size_t i = 0; i < num_primes; i++)
			delete[] tbl[i];
		delete[] tbl;
		tbl = 0;
	}//TOM

	init_(f);

	num_primes = fd.number_of_primes();
	tbl = new fft_prime_t*[num_primes];
	memory_handler(tbl, "fft_rep", "init(fft_rep&) :: "
		       "Error in memory allocation (tbl)");
	for (lidia_size_t i = 0; i < num_primes; i++)
		tbl[i] = 0;
	k = -1;
	max_k = -1;
}



void fft_rep::set_size(lidia_size_t l)
{
	if (fd.node ? (fd.node->max_degree < l) : 1) {
		lidia_error_handler("fft_rep", "set_size(...)::internal error");
		return;
	}
	if (l > max_k && l >= 0) {

//if (k != -1)	//#########
//    std::cout<<"   #in fft_rep::set_size; enlarge vec: "<<k<<" ->"<<l<<std::endl;

		for (lidia_size_t i = 0; i < num_primes; i++) {
			if (max_k >= 0)
				delete[] tbl[i];
			tbl[i] = new fft_prime_t[1 << l];
			memory_handler(tbl, "fft_rep", "init(...) :: "
				       "Error in memory allocation (tbl[i])");
		}
		max_k = l;

		if (base_fft_rep::stat_k < l) {
			delete[] base_fft_rep::stat_vec;
			base_fft_rep::stat_vec = new fft_prime_t[1 << l];
			memory_handler(base_fft_rep::stat_vec, "fft_rep", "init(...) :: "
				       "Error in memory allocation (vec)");
			base_fft_rep::stat_k = l;
		}
	}
	k = l;
	c->reset();
}



void fft_rep::to_fft_rep(const Fp_polynomial &f,
			 lidia_size_t lo, lidia_size_t hi)
{
	lidia_size_t i;
	for (i = 0; i < num_primes; i++)
		base_fft_rep::to_mod_rep(tbl[i], f, lo, hi, i);
}



void fft_rep::from_fft_rep(Fp_polynomial& x, lidia_size_t lo, lidia_size_t hi)
	// converts from FFT-representation to coefficient representation
	// only the coefficients lo..hi are computed
	// NOTE: this version does not destroy the data in (*this)
{
	debug_handler("fft_rep", "from_fft_rep(Fp_polynomial, lidia_size_t, lidia_size_t)");

	lidia_size_t K = (1 << k);
	hi = comparator< lidia_size_t >::min(hi, K-1);
	lidia_size_t l = comparator< lidia_size_t >::max(hi-lo+1, 0);

	fft_prime_t *sp = &stat_vec[lo];

	x.set_modulus(fd.crt_node->p);
	x.set_degree(l-1);
	bool ok = true;

	if (num_primes > 1) {
		fft_prime_t *yp;
		for (lidia_size_t index = num_primes-1; index >= 0; index--) {
			yp = tbl[index];

// old code:
//	    ok = ok && fd.node->prime[index].interpolate(stat_vec, yp, k);
			ok = ok && fd.node->prime[index].interpolate2(stat_vec, yp, k, lo, l);

			combine(sp, l, index);
		}

		get_res(x.coeff, l);
	}
	else {
// old code:
//    	ok = fd.node->prime[0].interpolate(stat_vec, tbl[0], k);
		ok = fd.node->prime[0].interpolate2(stat_vec, tbl[0], k, lo, l);

				// skip combine/get_res: we only used _one_ prime
				// => (stat_vec[lo..hi] % p) is almost already the desired coefficient
				//    vector
				// we only have to convert sp[i] from the representation as a
				// nonnegative integer in [0..fftprime] to the representation as
				// integer in [-z..z] with z=floor(fftprime/2)

		fft_prime_t prime = fd.node->prime[0].get_prime();
		fft_prime_t z = prime >> 1; // = floor(prime/2)
		long t;

		long p;
		if (fd.crt_node->p.longify(p))
			lidia_error_handler("fft_rep", "from_fft_rep(...)::internal error");

		for (lidia_size_t i = l-1; i >= 0; i--) {
			t = sp[i];
			if (t > long(z)) {
				t -= prime; // =  > -z <= t < 0
				t %= p; // =  > -p < t <= 0
				if (t != 0)
					t += p; // =  > 0 <= t < p
			}
			else
				t %= p;
			x.coeff[i].assign(t);
		}

	}
	if (!ok)
		lidia_error_handler("fft_rep",
				    "from_fft_rep(...)::re-init of FFT prime");

	x.remove_leading_zeros();
}



void multiply(fft_rep& z, const fft_rep& x, const fft_rep& y)
{
	lidia_size_t k = x.k;
	if (x.fd != y.fd || x.fd != z.fd || k != z.k || k != y.k) {
		lidia_error_handler("fft_rep", "multiply(...)::FFT rep mismatch");
		return;
	}

	lidia_size_t index;

	for (index = x.num_primes-1; index >= 0; index--) {
		fft_prime_t *zp = &z.tbl[index][0];
		const fft_prime_t *xp = &x.tbl[index][0];
		const fft_prime_t *yp = &y.tbl[index][0];

		x.fd.node->prime[index].pointwise_multiply(zp, xp, yp, k);
	}
}



void subtract(fft_rep &x, const fft_rep &a, const fft_rep &b)
{
	debug_handler("modular_fft_rep", "subtract(modular_fft_rep&, fft_rep&, modular_fft_rep&, lidia_size_t)");

	lidia_size_t k = a.k;
	if (a.fd != b.fd || a.fd != x.fd || k != b.k || k != x.k) {
		lidia_error_handler("fft_rep", "subtract(...)::FFT rep mismatch");
		return;
	}

	lidia_size_t index;
	for (index = a.num_primes-1; index >= 0; index--) {
		a.fd.node->prime[index].pointwise_subtract(x.tbl[index],
							   a.tbl[index], b.tbl[index], k);
	}
}



void reduce(fft_rep& x, const fft_rep& a, lidia_size_t k)
	// reduces a 2^l point FFT-rep to a 2^k point FFT-rep
	// input may alias output
{
	debug_handler("fft_rep", "reduce(fft_rep&, fft_rep&, lidia_size_t)");

	if (&x == &a)
		lidia_error_handler("fft_rep",
				    "reduce(...)::input may not alias output");

	if (x.fd != a.fd || k != x.k) {
		lidia_error_handler("fft_rep", "reduce(...)::FFT rep mismatch");
		return;
	}

	lidia_size_t index, j;
	fft_prime_t* xp;
	const fft_prime_t* ap;

	lidia_size_t n = 1 << k;
	lidia_size_t diff = a.k-k;
	if (diff < 0) {
		lidia_error_handler("fft_rep", "reduce(...)::bad operands");
		return;
	}

	for (index = a.num_primes-1; index >= 0; index--) {
		ap = &a.tbl[index][0];
		xp = &x.tbl[index][0];
		for (j = 0; j < n; j++, xp++)
			*xp = ap[j << diff];
	}
}



//*************************************************************************
//    three special purpose functions, used only in
//  void update_map(base_vector< bigint > & x, const base_vector< bigint > & a,
//		    const Fp_poly_multiplier& B, const Fp_poly_modulus& F)
//*************************************************************************

void fft_rep::rev_to_fft_rep(const base_vector< bigint > & x, lidia_size_t lo,
				     lidia_size_t hi, lidia_size_t offset)
	// computes an n = 2^k point convolution of X^offset*x[lo..hi]
	// using "inverted" evaluation points.
	// if deg(x) >= 2^k, then x is first reduced modulo X^n-1.
{
	debug_handler("fft_rep", "rev_to_fft_rep(base_vector< bigint > & x, lidia_size_t, lidia_size_t, lidia_size_t)");

	if (k < 0 || fd.node == 0) {
		lidia_error_handler("base_fft_rep", "to_mod_rep(...)::not initialized");
		return; // LC
	}

	if (lo < 0) {
		lidia_error_handler("fft_rep", "rev_to_fft_rep(base_vector< bigint > & "
				    "x, lidia_size_t, lidia_size_t, lidia_size_t)::bad arg");
		return; // LC
	}

	lidia_size_t n, i, j, m, j1;
	bigint accum;

	hi = comparator< lidia_size_t >::min(hi, x.size()-1);

	n = 1 << k;
	m = comparator< lidia_size_t >::max(hi-lo + 1, 0);

	const bigint *xx;
	if (x.size() == 0)
		xx = 0;
	else
		xx = &x[0];

	bool ok = true;
	if (m + offset <= n) {
		for (i = 0; i < num_primes; i++) {
			for (j = 0; j < offset; j++)
				stat_vec[j] = 0;
#if 0
			c->reduce(&stat_vec[offset], &xx[lo], m, i);
#else
			fft_prime_t q = c->get_prime(i);
			for (j = 0; j < m; j++)
				remainder(stat_vec[offset+j], xx[lo+j], q);
#endif
			for (j = m+offset; j < n; j++)
				stat_vec[j] = 0;

			ok = ok && fd.node->prime[i].interpolate(tbl[i], stat_vec, k);
		}
	}
	else {
		for (j = 0; j < n; j++) {
			if (j >= m) {
				for (i = 0; i < num_primes; i++)
					tbl[i][j] = 0;
			}
			else {
				accum = xx[j+lo];
				for (j1 = j + n; j1 < m; j1 += n)
					add(accum, accum, xx[j1+lo-offset]);
#if 0
				for (i = 0; i < num_primes; i++)
					c->reduce(tbl[i][(offset + j) % n], accum, i);
#else
				for (i = 0; i < num_primes; i++)
					remainder(tbl[i][(offset + j) % n], accum, c->get_prime(i));
#endif
			}
		}

		for (i = 0; i < num_primes; i++) {
			ok = ok && fd.node->prime[i].interpolate(stat_vec, tbl[i], k);
			for (j = 0; j < n; j++)
				tbl[i][j] = stat_vec[j];
		}
	}
	if (!ok)
		lidia_error_handler("fft_rep",
				    "rev_to_fft_rep(...)::re-init of FFT prime");
}



void fft_rep::rev_from_fft_rep(base_vector< bigint > & x, lidia_size_t lo, lidia_size_t hi)
	// converts from FFT-representation to coefficient representation
	// using "inverted" evaluation points.
	// only the coefficients lo..hi are computed
{
	debug_handler("fft_rep", "rev_from_fft_rep(lidia_size_t, lidia_size_t)")
		;
	lidia_size_t n, i, l;
	n = (1 << k);

	hi = comparator< lidia_size_t >::min(hi, n-1);
	l = comparator< lidia_size_t >::max(hi-lo+1, 0);

	if (x.capacity() < l)
		x.set_capacity(l);
	x.set_size(l);
	bool ok = true;

	if (l != 0) {
		bigint* xx = &x[0]; //x.get_data_address();

		if (num_primes > 1) {
			for (i = 0; i < num_primes; i++) {
				ok = ok && fd.node->prime[i].evaluate(stat_vec, 0, tbl[i], n-1, k);

				combine(&stat_vec[lo], l, i);
			}
			get_res(xx, l);
		}
		else {
			ok = fd.node->prime[0].evaluate(stat_vec, 0, tbl[0], n-1, k);

			// skip combine/get_res: we only used _one_ prime
			// => (stat_vec[lo..hi] % p) already is the desired coeff. vector

			const bigint &p = fd.crt_node->p;
			for (i = l-1; i >= 0; i--)
				Remainder(xx[i], (stat_vec+lo)[i], p);
		}
	}
	if (!ok)
		lidia_error_handler("fft_rep",
				    "rev_from_fft_rep(...)::re-init of FFT prime");
}



void fft_rep::add_expand(const fft_rep& a)
	//  x = x + (an "expanded" version of a)
{
	debug_handler("fft_rep", "add_expand(fft_rep&)");
	lidia_size_t j, j1, index;
	lidia_size_t a_k = a.k;
	lidia_size_t n = 1 << a_k;

	if (fd != a.fd) {
		lidia_error_handler("fft_rep", "add_expand(fft_rep&)::"
				    "Reps do not match");
		return;
	}
	if (k < a_k) {
		lidia_error_handler("fft_rep", "add_expand(fft_rep&)::bad args");
		return;
	}

	for (index = 0; index < num_primes; index++) {
		fft_prime_t q = c->get_prime(index);
		const fft_prime_t *ap = a.tbl[index];
		fft_prime_t *xp = tbl[index];

		for (j = 0; j < n; j++) {
			j1 = j << (k-a_k);
			xp[j1] = add_mod(xp[j1], ap[j], q);
		}
	}
}



//*************************************************************************
//	    arithmetic with modular_fft_rep and fft_rep
//*************************************************************************

void multiply(modular_fft_rep &x, const fft_rep &a, const modular_fft_rep &b,
	      lidia_size_t index)
{
	debug_handler("modular_fft_rep", "multiply(modular_fft_rep&, fft_rep&, modular_fft_rep&, lidia_size_t)");

	lidia_size_t k = a.k;
	if (a.fd != b.fd || a.fd != x.fd || k != b.k || k != x.k) {
		lidia_error_handler("modular_fft_rep",
				    "multiply(...)::FFT rep mismatch");
		return;
	}

	a.fd.node->prime[index].pointwise_multiply(x.vec, a.tbl[index], b.vec, k);
}



void subtract(modular_fft_rep &x, const fft_rep &a, const modular_fft_rep &b,
	      lidia_size_t index)
{
	debug_handler("modular_fft_rep", "subtract(modular_fft_rep&, fft_rep&, modular_fft_rep&, lidia_size_t)");

	lidia_size_t k = a.k;
	if (a.fd != b.fd || a.fd != x.fd || k != b.k || k != x.k) {
		lidia_error_handler("modular_fft_rep",
				    "subtract(...)::FFT rep mismatch");
		return;
	}

	a.fd.node->prime[index].pointwise_subtract(x.vec, a.tbl[index], b.vec, k);
}



void reduce(fft_rep &x, const modular_fft_rep &a, lidia_size_t k,
	    lidia_size_t index)
{
	debug_handler("fft_rep", "reduce(fft_rep&, modular_fft_rep&, lidia_size_t)");

	if (x.fd != a.fd || k != x.k) {
		lidia_error_handler("fft_rep", "reduce(...)::FFT rep mismatch");
		return;
	}

	fft_prime_t* xp;
	const fft_prime_t* ap;

	lidia_size_t n = 1 << k;
	lidia_size_t diff = a.k-k;
	if (diff < 0) {
		lidia_error_handler("fft_rep", "reduce(...)::bad operands");
		return;
	}

	ap = a.vec;
	xp = x.tbl[index];
	for (lidia_size_t j = 0; j < n; j++, xp++)
		*xp = ap[j << diff];
}



void reduce(modular_fft_rep &x, const fft_rep &a, lidia_size_t k,
	    lidia_size_t index)
{
	debug_handler("modular_fft_rep", "reduce(modular_fft_rep&, fft_rep&, lidia_size_t)");

	if (x.fd != a.fd || k != x.k) {
		lidia_error_handler("modular_fft_rep", "reduce(...)::FFT rep mismatch");
		return;
	}

	fft_prime_t* xp;
	const fft_prime_t* ap;

	lidia_size_t n = 1 << k;
	lidia_size_t diff = a.k-k;
	if (diff < 0) {
		lidia_error_handler("fft_rep", "reduce(...)::bad operands");
		return;
	}

	ap = a.tbl[index];
	xp = x.vec;
	for (lidia_size_t j = 0; j < n; j++, xp++)
		*xp = ap[j << diff];
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
