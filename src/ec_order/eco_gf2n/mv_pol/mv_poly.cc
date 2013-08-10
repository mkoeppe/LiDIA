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
//	Author	: Volker Mueller
//	Changes	: See CVS log
//
//==============================================================================================



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/mv_poly.h"
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



#define QUICKSORT 1            // as standard, quicksort is used


mv_list_entry *mv_poly::free_list = NULL;
int mv_poly::size_free_list = 0;
int mv_poly::refc = 0;
static mv_list_entry ** listcp;
static int listcp_size = 0;

const int DELETE_FREE = 200000; // the # of terms when free_list is deleted
const int mult_blocking_factor = 50; // used in multiply



//
// functions for handling the free_list
//

//--------------------------------------------------------------------
// add all the list elements of p to free_list and set p's pointer
// to NULL

void mv_poly::put_in_free(mv_poly & p)
{
	if (p.first_term != NULL) {
		mv_poly::size_free_list += p.length();
		p.last_term->next = mv_poly::free_list;
		mv_poly::free_list = p.first_term;
		p.first_term = NULL;
		p.last_term = NULL;
	}

	if (mv_poly::size_free_list > DELETE_FREE)
		delete_free_list(DELETE_FREE / 2);
}



//--------------------------------------------------------------------
// add one list element to the free list.

void mv_poly::put_in_free(mv_list_entry & t)
{
	mv_poly::size_free_list ++;

	t.next = mv_poly::free_list;
	mv_poly::free_list = & t;

	if (mv_poly::size_free_list > DELETE_FREE)
		delete_free_list(DELETE_FREE/2);

}



//--------------------------------------------------------------------
// return a pointer to a free list element, remove that element
// from the free_list.

mv_list_entry* mv_poly::get_new(unsigned int anz)
{
	mv_list_entry* temp = NULL;

	if (anz == 1) {
		if (mv_poly::free_list == NULL)
			temp = new mv_list_entry;
		else {
			mv_poly::size_free_list --;
			temp = mv_poly::free_list;
			mv_poly::free_list = temp->next;
		}
		temp->next = NULL;
	}
	else {
		mv_list_entry * p;

		if (anz == 0)
			lidia_error_handler("mv_poly", "get_new::size is zero");

		if (mv_poly::size_free_list >= static_cast<int>(anz)) {
			mv_poly::size_free_list -= anz;
			temp = p = mv_poly::free_list;
			anz --;

			while (anz >= 1) {
				p = p->next;
				anz --;
			}
			mv_poly::free_list = p->next;
			p->next = NULL;
		}
		else {
			p = mv_poly::free_list;
			anz -= mv_poly::size_free_list;
			mv_poly::size_free_list = 0;
			mv_poly::free_list = NULL;

			while (anz > 0) {
				anz --;
				temp = new mv_list_entry;
				temp->next = p;
				p = temp;
			}
		}

		p = temp;
		int i = 0;
		while (p != NULL) {
			p = p->next; i++;
		}
	}
	return temp;
}



//----------------------------------------------------------
// delete all the elements in the free_list.

void mv_poly::delete_free_list(unsigned int anz)
{
	mv_list_entry *p;

	while (mv_poly::free_list != NULL && anz > 0) {
		p = mv_poly::free_list;
		mv_poly::free_list = p->next;
		delete p;
		anz --;
		mv_poly::size_free_list --;
	}

	if (listcp_size > 0)
		delete[] listcp;
	listcp_size = 0;

}



//--------------------------------------------------------------
// assignments for mv_poly
//

void swap(mv_poly & a, mv_poly & b)
{
	mv_list_entry * first, *last;
	int l;

	first = a.first_term;
	last = a.last_term;
	l = a.len;

	a.first_term = b.first_term;
	a.last_term = b.last_term;
	a.len = b.len;

	b.first_term = first;
	b.last_term = last;
	b.len = l;
}



//-----------------------------------------------------------------

#if 0
int mv_poly::length() const
{
	mv_list_entry *p;
	int i = 0;

	p = first_term;

	while (p != NULL) {
		i++;
		p = p->next;
	}

	if (i != len) {
		std::cout<<"\n\nWARNING: computed length = "<<i<<" len = "<<len<<std::flush;
		lidia_error_handler("", "");
	}

	return i;
}
#endif



//----------------------------------------------------------------------

void mv_poly::assign_one()              // list with one list element
{
	mv_poly::put_in_free(*this);
	first_term = last_term = mv_poly::get_new(1);
	last_term->term.assign_one();
	last_term->next = NULL;
	len = 1;
}



//---------------------------------------------------------------------

mv_poly & mv_poly::operator = (const mv_poly & p)
{
	if (&p == this)
		return *this;

	if (! (this->is_zero()))
		mv_poly::put_in_free(*this);

	if (p.is_zero()) {
		assign_zero();
		return *this;
	}

	const mv_list_entry* cursor; // cursor runs through p
	mv_list_entry* temp; // temp runs through *this

	first_term = mv_poly::get_new(p.length());
	cursor = p.first_term;
	first_term->term.assign(cursor->term);
	temp = first_term;
	len = 1;

	while (cursor->next != NULL)   // copy the complete list
	{
		cursor = cursor->next;
		temp = temp->next;
		temp->term.assign(cursor->term);
		len ++;
	}
	last_term = temp;
	temp->next = NULL;
	return *this;
}



//--------------------------------------------------------------------

mv_poly & mv_poly::operator = (const mv_term & t)
{
	mv_poly::put_in_free(*this);
	first_term = last_term = mv_poly::get_new(1);
	first_term->term.assign(t);
	len = 1;
	return *this;
}



//--------------------------------------------------------------
// comparisons of mv_polynomials, equal and trivial tests like
// is_one, is_zero, is_const
//

int mv_poly::max_index_variable() const
{
	bigint h(0);
	mv_list_entry *p;

	if (is_zero())
		return 0;

	p = first_term;

	while (p != NULL) {
		bitwise_or(h, h, p->term.get_var());
		p = p->next;
	}
	return h.bit_length() - 1;
}



bool operator == (const mv_poly & a, const mv_poly & b)
{
	// if (a.len != b.len)
	// return false;

	if ((a.is_zero() && b.is_zero()) || (&a == &b))
		return true;

	if ((a.is_zero() && !b.is_zero()) || (!a.is_zero() && b.is_zero()))
		return false;

	const mv_list_entry* cursor1;
	const mv_list_entry* cursor2;

	cursor1 = a.first_term;
	cursor2 = b.first_term;

	while (cursor1 != NULL) {                             // note: not both can be NULL here !!
		if (cursor2 == NULL)
			return false;

		if (! cursor1->term.all_equal(cursor2->term))
			return false;

		cursor2 = cursor2->next;
		cursor1 = cursor1->next;
	}

	if (cursor2 != NULL)
		return false;
	else
		return true;
}



//------------------------------------------------------------------

bool mv_poly::is_one() const
{
	if (is_zero())
		return false;

	if (first_term->term.is_one() && first_term->term.get_var().is_zero()
	    && first_term->next == NULL)
		return true;
	else
		return false;
}



//-----------------------------------------------------------------
//  returns true iff mv_poly is non-zero (!!) constant

bool mv_poly::is_const() const
{
	if (is_zero())
		return false;

	if (first_term->term.is_const() && first_term->next == NULL)
		return true;
	else
		return false;
}



//-----------------------------------------------------------------
// arithmetic functions: addition just uses +=, multiplication *=
//

void add(mv_poly & c, const mv_poly & a, const mv_poly & b)
{
	if (&a == &b) {
		c.assign_zero();
		return;
	}

	if (&a == &c) {
		c += b;
		return;
	}

	if (&b == &c) {
		c += a;
		return;
	}

	if (a.is_zero()) {
		c = b;
		return;
	}

	if (b.is_zero())
		c.assign(a);
	else {
		c.assign(b);
		c += a;
	}
}



//---------------------------------------------------------------------

void multiply(mv_poly & c, const mv_poly & a, const mv_poly & b)
{
	if (a.is_zero() || b.is_zero()) {
		c.assign_zero();
		return;
	}

	if (&a == &c) {
		c *= b;
		return;
	}

	if (&b == &c)
		c *= a;
	else {
		c.assign(b);
		c *= a;
	}
}



//-----------------------------------------------------------------

void multiply(mv_poly &c, const gf2n & a, const mv_poly & b)
{
	if (a.is_zero() || b.is_zero()) {
		c.assign_zero();
		return;
	}
	else
		c.assign(b);

	mv_list_entry* cursor;

	for (cursor = c.first_term; cursor != NULL; cursor = cursor->next)
		multiply(cursor->term, cursor->term, a);

}



//-----------------------------------------------------------------

void square(mv_poly & c, const mv_poly & a)
{
	if (a.is_zero()) {
		c.assign_zero();
		return;
	}
	else
		c.assign(a);

	mv_list_entry* cursor;

	for (cursor = c.first_term; cursor != NULL; cursor = cursor->next)
		square(cursor->term, cursor->term);
}



//-----------------------------------------------------------------------

void sqrt(mv_poly & c, const mv_poly & a)
{
	if (a.is_zero()) {
		c.assign_zero();
		return;
	}
	else
		c.assign(a);

	mv_list_entry* cursor;

	for (cursor = c.first_term; cursor != NULL; cursor = cursor->next)
		sqrt(cursor->term, cursor->term);
}



//-------------------------------------------------------------
// the first difficult operator +=
//
//  simultaneous run over both lists and selective insertion at
//  the right places --> order is maintained without resorting.
//

mv_poly & mv_poly::operator += (const mv_poly & a)
{
	if (a.is_zero())
		return *this;

	if (is_zero()) {
		this->assign(a);
		return *this;
	}

	if (this == &a) {
		assign_zero();
		return *this;
	}

	const mv_list_entry* cursor1; // runs through 'a'
	mv_list_entry* cursor2, *last; // run through '*this'

	cursor1 = a.first_term;
	last = cursor2 = first_term;

	do {
		// Find the right position to insert the term to which points
		// cursor1

		while ((cursor1->term > cursor2->term) && (cursor2->next != NULL)) {
			last = cursor2;
			cursor2 = cursor2->next;
		}

		// Three cases, we have to deal with

		if (cursor1->term > cursor2->term) {
			// Insert at the end of *this

			cursor2->next = mv_poly::get_new(1);
			len++;

			last = cursor2;
			cursor2 = cursor2->next;
			cursor2->term = cursor1->term;
			cursor2->next = NULL;
			last_term = cursor2;
		}
		else {
			if (cursor1->term == cursor2->term)  // add the terms, this might
			{                                    // result in zero -> deletion
				add(cursor2->term, cursor2->term, cursor1->term);

				if (((cursor2->term).get_coeff()).is_zero()) {
					len --;
					// Remove the zero-term pointed to by cursor2 from *this

					if (cursor2 == last) {
						// remove at the beginning of the list
						last = first_term = cursor2->next;
						mv_poly::put_in_free(*cursor2);
						cursor2 = last;

						if (cursor2 == NULL && cursor1->next != NULL) {
							first_term = mv_poly::get_new(1);
							cursor1 = cursor1->next;
							first_term->term.assign(cursor1->term);
							cursor2 = first_term;
							len ++;

							while (cursor1->next != NULL) {
								// copy the complete list
								len ++;
								cursor1 = cursor1->next;
								cursor2->next = mv_poly::get_new(1);
								cursor2 = cursor2->next;
								cursor2->term.assign(cursor1->term);
							}
							last_term = cursor2;
							return *this;
						}
					}
					else {
						last->next = cursor2->next;
						mv_poly::put_in_free(*cursor2);
						cursor2 = last;
					}
					if (first_term == NULL)
						last_term = NULL;
					else
						if (last->next == NULL)
							last_term = last;
				}
			}
			else {
				// Insert a copy of *cursor1 between *last and *cursor2

				len++;

				if (cursor1->term < last->term) {
					first_term = mv_poly::get_new();
					first_term->term = cursor1->term;
					first_term->next = last;
					last = first_term;
					cursor2 = last->next;
				}
				else {
					last->next = mv_poly::get_new();
					last->next->term = cursor1->term;
					last->next->next = cursor2;
					cursor2 = last->next;
				}
			}
		}
		cursor1 = cursor1->next;
	} while (cursor1 != NULL);

	return *this;
}



//--------------------------------------------------------------
// operator *= for a mv_term only
//
// We can not be sure that the ordering is destroyed and have to
// resort.

mv_poly& mv_poly::operator *= (const mv_term & t)
{
	if (is_zero() || t.is_zero()) {
		assign_zero();
		return *this;
	}

	mv_list_entry* cursor;
	lidia_size_t counter = 0; // used in the sorting function

	for (cursor = first_term; cursor != NULL; cursor = cursor->next) {
		counter ++;
		cursor->term *= t;
	}

	sort(counter); // sort assigns len

	return *this;
}



mv_poly & mv_poly::operator *= (const mv_poly & a)
{
	if (is_zero() || a.is_zero()) {
		assign_zero();
		return *this;
	}

	if (this == &a) {
		square(*this, a);
		return *this;
	}

	mv_poly this_copy, th2, th3;
	mv_term t;
	const mv_list_entry* cursor;
	mv_list_entry *cursor2;
	int mb = mult_blocking_factor, block_factor;

	this_copy.assign(*this);
	assign_zero();

	//if (this_copy.length() > 100)
	//   mb = 1;

	for (cursor = a.first_term; cursor != NULL; ) {
		block_factor = 1;
		th2.assign_zero();

		while (block_factor <= mb && cursor != NULL) {
			block_factor ++;
			t = cursor->term;
			cursor = cursor->next;
			if (t.is_zero())
				continue;

			th3.assign(this_copy);

			for (cursor2 = th3.first_term; cursor2 != NULL; cursor2 = cursor2->next) {
				multiply(cursor2->term, cursor2->term, t);
				th2.len ++;
			}

			if(th2.is_zero()) {
				th2.first_term = th3.first_term;
				th2.last_term = th3.last_term;
			}
			else {
				th2.last_term->next = th3.first_term;
				th2.last_term = th3.last_term;
			}
			th3.first_term = th3.last_term = NULL;
			th3.len = 0;
		}
		if (th2.len > 0) {
			th2.sort(th2.len);
			*this += th2;
		}
	}

	return *this;
}



//---------------------------------------------------------------------
// the sorting routine that uses the qsort functions from the C library.
// Note that we have to copy elements first in an array of pointers.
// Double monomials are combined. Counter has to hold the number of
// list elements of *this.
//
// There are two available in functions for sorting:
//  #define QUICKSORT --> use the stdlib qsort function
//  #define HEAPSORT --> use a self written heapsort

#ifdef QUICKSORT
extern "C" int cmp_mv_term(const void * aa, const void * bb)
{
	const mv_term *at, *bt;

	at = &((*(static_cast<const mv_list_entry * const *>(aa)))->term);
	bt = &((*(static_cast<const mv_list_entry * const *>(bb)))->term);

	if (*at < *bt)
		return -1;

	if (*at == *bt)
		return 0;
	else
		return 1;
}



#else




inline
void heapify(mv_list_entry ** A, int i, int size_A)
{
	register int l, r, largest;
	mv_list_entry *p;
	bool end;

	do {
		l = 2*i+1;
		r = 2*i+2;
		end = true;

		if (l < size_A && A[l]->term > A[i]->term)
			largest = l;
		else
			largest = i;

		if (r < size_A &&  A[r]->term > A[largest]->term)
			largest = r;

		if (largest != i) {
			p = A[largest];
			A[largest] = A[i];
			A[i] = p;
			i = largest;
			end = false;
		}
	} while(!end);
}
#endif



//---------------------------------------------------------

void mv_poly::sort(lidia_size_t number_of_terms)
{
	if (number_of_terms <= 0) {
		len = 0;
		return;
	}

	mv_list_entry* cursor;
	lidia_size_t i;

	if (listcp_size < number_of_terms) {
		if (listcp_size > 0)
			delete[] listcp;

		listcp = new mv_list_entry* [number_of_terms];
		listcp_size = number_of_terms;
	}

	cursor = first_term; // first copy the polynomial in an array
	i = 0;

	do {
		listcp[i] = cursor;
		cursor = cursor->next;
		i ++;
	} while (cursor != NULL);

	// set *this to zero without destroying the list elements

	first_term = cursor = NULL;
	len = 0;

	//-------------------- now call the sorting method

#ifdef QUICKSORT
	std::qsort (&(listcp[0]), number_of_terms, sizeof(mv_list_entry*),
		    cmp_mv_term);

#else           // Heapsort

	int size_A = number_of_terms;
	mv_list_entry *p;

	for (int i = (number_of_terms - 1)/2; i >= 0; i--)
		heapify(listcp, i, size_A);

	for (int i = number_of_terms - 1; i >= 1; i--) {
		p = listcp[0];
		listcp[0] = listcp[i];
		listcp[i] = p;
		size_A --;
		heapify(listcp, 0, size_A);
	}
#endif



	// copy the array back into the list
	// in parallel we combine terms with the same variables and delete
	// zero terms which might show up because combination
	// the algorithm compares *act with the next unread term, if they are
	// differentand *act is not zero, then *act is inserted, otherwise
	// *act is deleted or the two terms are combined


	mv_list_entry *act = listcp[0];
	len = 0;

	for(i = 1; i < number_of_terms; i ++) {
		if (listcp[i]->term.get_var() ==
		    act->term.get_var()) {
			// combine the two terms
			act->term += listcp[i]->term;
			mv_poly::put_in_free(*listcp[i]);
			listcp[i] = NULL;
			continue;
		}

		if (! act->term.is_zero()) {
			len ++;
			if (first_term == NULL) {
				// insert at the beginning
				first_term = cursor = act;
				first_term->next = NULL;
			}
			else {
				cursor->next = act;
				cursor = cursor->next;
				cursor->next = NULL;
			}
		}
		else
			mv_poly::put_in_free(*act);

		act = listcp[i];
		listcp[i] = NULL;
	}

	if (!act->term.is_zero()) {
		// finally the last term, still stored in *act
		len ++;
		if (first_term == NULL) {
			first_term = cursor = act;
			first_term->next = NULL;
		}
		else {
			cursor->next = act;
			cursor = cursor->next;
			cursor->next = NULL;
		}
	}
	else
		mv_poly::put_in_free(*act);

	last_term = cursor;
}



//--------------------------------------------------------------------
// output
//

std::ostream& operator << (std::ostream & out, const mv_poly & a)
{
	const mv_list_entry* cursor;

	if (a.is_zero())
		out << "0";
	else {
		cursor = a.first_term;
		do {
			out << cursor->term;
			cursor = cursor->next;
			if (cursor != NULL)
				out<<" + ";
		} while (cursor != NULL);
	}
	return out;
}



//-------------------------------------------------------------------
// solve X_k = q if there is some term in q with a single
// variable X_k
//
// q is set to the terms that do not contain variable X_k, divided by
// the coefficient of the X_k term. The function also determines and
// returns some suitable k.
//
// Only in the case that there is some linear term in *this,
// we solve for the corresponding variable. The function then returns true,
// otherwise it returns false.

bool solve_x_k(mv_poly & q, lidia_size_t & k)
{
	mv_list_entry* before_term;
	mv_list_entry* check_term;
	mv_list_entry* temp;
	gf2n coeff_lin, coeff_term;

	// before is set to the last term before the term which has some single
	// variable

	k = -1;

	if ((before_term = q.has_one_linear_term(k)) == NULL)
		return false;

	if (before_term == q.first_term) {
		// check whether the first or the second term are the correct one
		if (before_term->term.is_linear() &&
		    before_term->term.get_var().bit_length() - 1 == k) {
			temp = before_term;
			q.first_term = temp->next;
			if (before_term == q.last_term) {
				// only one (linear) term in list
				q.assign_zero();
				mv_poly::put_in_free(*temp);
				return true;
			}
		}
		else {
			temp = before_term->next;
			before_term->next = temp->next;
		}
	}
	else {
		// linear term is not the first one
		temp = before_term->next;
		before_term->next = temp->next;
	}

	// now term 'temp' contains single variable X_k

	coeff_lin = (temp->term).get_coeff(); // eliminate that term from list

	if (temp == q.last_term)
		q.last_term = before_term;

	mv_poly::put_in_free(*temp);

	// divide all remaining terms by the coeff. of the lin. term
	if (!coeff_lin.is_one()) {
		invert(coeff_lin, coeff_lin);
		for (check_term = q.first_term; check_term != NULL;
		     check_term = check_term->next) {
			multiply(coeff_term, check_term->term.get_coeff(), coeff_lin);
			check_term->term.assign_coeff(coeff_term);
		}
	}
	q.len --;

	if (q.first_term == NULL)
		q.last_term = NULL;

	return true;
}



//------------------------------------------------------------------
//
// as above, but k is chosen by the user
// on entry, q has a term with single variable k, on exit q is changed
// and does no longer have that variable

bool solve_x_k_fixed(mv_poly & q, lidia_size_t k)
{
	mv_list_entry* before_term;
	mv_list_entry* check_term;
	mv_list_entry* temp;
	gf2n coeff_lin, coeff_term;

	if ((before_term = q.has_one_linear_term(k)) == NULL)
		return false;

	if (before_term->term.is_linear() && before_term->term.has_var_k(k)) {
		// first term of list
		temp = before_term;
		q.first_term = temp->next;
		if (before_term == q.last_term) {
			// only one (linear) term in list
			k = temp->term.get_var().bit_length() - 1;
			q.assign_zero();
			mv_poly::put_in_free(*temp);
			return true;
		}
	}
	else {
		// linear term is not the first one
		temp = before_term->next;
		before_term->next = temp->next;
	}

	// now term 'temp' contains single variable X_k

	coeff_lin = (temp->term).get_coeff(); // eliminate that term from list

	if (temp == q.last_term)
		q.last_term = before_term;

	mv_poly::put_in_free(*temp);

	// divide all remaining terms by the coeff. of the lin. term
	if (!coeff_lin.is_one()) {
		invert(coeff_lin, coeff_lin);
		for (q.len = 0, check_term = q.first_term; check_term != NULL;
		     check_term = check_term->next, q.len ++) {
			multiply(coeff_term, check_term->term.get_coeff(), coeff_lin);
			check_term->term.assign_coeff(coeff_term);
		}
	}

	if (q.first_term == NULL)
		q.last_term = NULL;
	return true;
}



//----------------------------------------------------------------------
// substitute(q, p, k)
//
// all occurances of variable X_k in q are substituted with the
// mv_poly p (which does not have terms with variable X_k),
// the result is returned in q.
//

void substitute(mv_poly & q, const mv_poly & p, lidia_size_t k)
{
	const mv_list_entry* act_term2;
	mv_poly temp;
	mv_poly q_k;

	if (q.is_zero())
		return;

	//  terms of p are not allowed to have var. x_k

	for (act_term2 = p.first_term; act_term2 != NULL;
	     act_term2 = act_term2->next)
		if ((act_term2->term).has_var_k(k)) {
			lidia_error_handler("mv_poly", "void substitute(mv_poly & q, "
					    "const mv_poly & p, lidia_size_t k)"
					    "::p includes var k.");
			return;
		}

	q.split_x_k(q_k, k); // splits q into all term without X_k (stored in
	q_k *= p; // new q) and terms having X_k (in q_k, with
	q += q_k; // deleting X_k from all terms)
}



//---------------------------------------------------------------------------
// evaluate:
//
// all terms in t will be evaluated by the values for the variables X_1 ...
// X_max that are in stored in c, i.e. if the i-th bit of c is 1, 0, then
// X_i is set to 1, 0, respectively. The result is returned in a.

void evaluate(gf2n & a, const mv_poly & t, const bigint & c)
{
	if (t.is_zero()) {
		a.assign_zero();
		return;
	}

	if (c.is_zero()) {
		// t is sorted !!
		if (t.first_term->term.is_const())
			a.assign(t.first_term->term.get_coeff());
		else
			a.assign_zero();
		return;
	}

	bigint h;
	const mv_list_entry* act_term;

	a.assign_zero();

	for (act_term = t.first_term; act_term != NULL;
	     act_term = act_term->next) {
		h =  act_term->term.get_var();
		// all the bits which are different in act_var, c
		// no 1-bit in h
		if (((h ^ c) & h).is_zero())
			add(a, a, act_term->term.get_coeff());
	}
}



//---------------------------------------------------------------------------
// evaluate:
//
// all terms in t will be evaluated by the values for the variables X_1 ...
// X_max that are in stored in c, i.e. if the i-th bit of c is 1, 0 and the
// i-th bit in bits is also 1, , then
// X_i is set to 1, 0, respectively. The result is returned in a.

void evaluate(mv_poly & a, const bigint & cc, const bigint & bits)
{
	if (bits.is_zero() || a.is_zero())
		return;

	bigint h, c;
	mv_list_entry* act_term, *p, *before_act_term;
	bool beginning = true;

	bitwise_and(c, cc, bits);

	before_act_term = act_term = a.first_term;

	do {
		substitute(act_term->term, c, bits);
		if (act_term->term.is_zero()) {
			act_term = act_term->next;
			mv_poly::put_in_free(*before_act_term);
			before_act_term = act_term;
		}
		else
			beginning = false;
	} while (beginning && act_term != NULL);

	a.first_term = act_term;
	if (a.first_term != NULL)
		a.len = 1;
	else {
		a.last_term = NULL; a.len = 0; return;
	}
	act_term = act_term->next;

	while(act_term != NULL) {
		substitute(act_term->term, c, bits);
		if (act_term->term.is_zero()) {
			p = act_term;
			act_term = act_term->next;
			mv_poly::put_in_free(*p);
			before_act_term->next = act_term;
		}
		else {
			a.len ++;
			before_act_term = act_term;
			act_term = act_term->next;
		}
	}
	a.last_term = before_act_term;
	a.clean();
}



//------------------------------------------------------------------
//  checks if there is exactly a linear term in X_k and X_k does not
//  occur in any other list element. If k = -1, then a suitable value
//  for k is computed (if it exists).
//
//  If there exists a term like described above, the function returns
//  a pointer to the term before that linear term (if a predecessor
//  exists) or the first list element (if that is a suitable term),
//  if there exists no suitable term, it returns NULL.

mv_list_entry* mv_poly::has_one_linear_term(lidia_size_t & k)
{
	if (is_zero())
		return NULL;

	mv_list_entry* act_term;
	mv_list_entry* ret_term;
	mv_list_entry* check_term;
	bool ret = true;
	lidia_size_t k_act;

	act_term = first_term;

	do {
		if (act_term->term.is_linear()) {
			if (k == -1 ||
			    act_term->term.get_var().bit_length() - 1 == k) {
				k_act = act_term->term.get_var().bit_length() - 1;

				// all other terms are not allowed to have var. X_{k_act}

				ret_term = first_term;

				for (check_term = first_term; check_term != NULL;
				     check_term = check_term->next) {
					if (check_term->next == act_term)
						ret_term = check_term;

					if (check_term != act_term &&
					    check_term->term.has_var_k(k_act)) {
						if (k != -1)     // we found a second term with X_k
							return NULL;
						ret = false;
						break;
					}
				}
				if (ret) {
					k = k_act;
					return ret_term;
				}
			}
		}

		act_term = act_term->next;
		ret = true;
	} while (act_term != NULL);
	return NULL;
}



//-----------------------------------------------------------------------
// this function splits the polynomial in *this in two parts:
// the part with the terms not containing X_k remains in *this
// terms containing X_k are collected and returned in q_k
// and the variable X_k is eliminated in all these terms.

void mv_poly::split_x_k (mv_poly & q_k, lidia_size_t k)
{
	mv_list_entry* cursor; // used for *this list
	mv_list_entry* before_cursor; // used for *this list
	mv_list_entry* last; // used for q_k list
	bigint mask, h;

	shift_left(mask, bigint(1), k);

	q_k.assign_zero();
	last = q_k.first_term;
	before_cursor = cursor = first_term;
	len = 0;
	q_k.len = 0;

	do {
		if (cursor->term.has_var_k (k)) {
			q_k.len ++;
			if (last != NULL)
				// add to q_k
				last->next = cursor;
			else
				q_k.first_term = cursor;
			last = cursor;

			if (before_cursor == cursor) {
				// delete from *this
				first_term = before_cursor = cursor->next;
			}
			else
				before_cursor->next = cursor->next;

			cursor = cursor->next;
			last->next = NULL;

			h.assign(last->term.get_var()); // eliminate variable X_k
			bitwise_or (h, h, !h);
			bitwise_xor (h, h, mask);
			last->term.assign_var(h & last->term.get_var());
		}
		else {
			len ++;

			if (before_cursor != cursor)
				before_cursor = before_cursor->next;
			cursor = cursor->next;
		}
	} while (cursor != NULL);

	if (q_k.first_term == NULL)
		q_k.last_term = NULL;
	else
		q_k.last_term = last;

	if (before_cursor != NULL)
		before_cursor->next = NULL;

	last_term = before_cursor;
}



//---------------------------------------------------------------------
// clean:
// assumes that the multivariate polynomial is sorted according the
// variable-indices, combines terms with the same variables.

void mv_poly::clean()
{
	mv_list_entry* cursor, *act;
	mv_list_entry* last, *p;

	if (is_zero()) {
		len = 0;
		return;
	}

	act = first_term;
	cursor = first_term->next;
	len = 0;
	first_term = last = NULL;

	while (cursor != NULL) {
		if(cursor->term.get_var() == act->term.get_var()) {
			act->term += cursor->term;
			p = cursor->next;
			mv_poly::put_in_free(*cursor);
			cursor = p;
		}
		else {
			if (!act->term.is_zero()) {
				len ++;
				if (first_term == NULL) {
					first_term = act;
					last = act;
				}
				else {
					last->next = act;
					last = act;
				}
				last->next = NULL;
			}
			act = cursor;
			cursor = cursor->next;
		}
	}

	if (!act->term.is_zero()) {
		len ++;
		if (first_term == NULL) {
			first_term = act;
			last = act;
		}
		else {
			last->next = act;
			last = act;
		}
		last->next = NULL;
	}

	last_term = last;
}



//-----------------------------------------------------------------------
//
// this function assumes, that the polynomial is sorted according
// variable indices and is clean (does not contain several terms in one
// variable). It returns k if the polynomial has the form  'coeff_1 +
// coeff_2 * X_k' for some k (coeff_1 arbitrary, coeff_2 != 0), otherwise
// it returns -1.


lidia_size_t mv_poly::linear_poly_in_only_one_var() const
{
	if (is_zero() || is_one())
		return -1;

	mv_list_entry *act_term;

	act_term = first_term;

	if (act_term->term.is_linear() && act_term->next == NULL)
		return (act_term->term.get_var().bit_length() - 1);
	else
		if (act_term->term.is_const() && act_term->next != NULL
		    && act_term->next->term.is_linear()
		    && act_term->next->next == NULL)
			return (act_term->next->term.get_var().bit_length() - 1);
		else
			return -1;
}



//----------------------------------------------------------------------
// has_var_k returns true iff there is some term which has variable X_k
//

bool mv_poly::has_var_k(lidia_size_t k) const
{
	if (is_zero())
		return false;

	const mv_list_entry * cursor;

	for(cursor = first_term; cursor != NULL; cursor = cursor->next)
		if(cursor->term.has_var_k (k))
			return true;
	return false;
}



//--------------------------------------------------------------------
//
// this function returns a polynomial T, whose coefficients are just
// the trace of the corresponding terms in *this

void mv_poly::trace_computation(mv_poly & Tr)
{
	mv_list_entry *term;
	gf2n Cc;
	int t;

	Tr.assign_zero();

	for (term = first_term; term != NULL; term = term->next) {
		Cc = term->term.get_coeff();
		t = Cc.trace();
		if (t == 1)
			Tr += mv_term(1, term->term.get_var());
	}
}



//--------------------------------------------------------------------
// substitute the k-th variable in q by n/d, multiply by denominator
// and return result in q

void substitute(mv_poly & q, const mv_poly & n, const mv_poly & d,
		lidia_size_t k)
{
	const mv_list_entry* act_term2;
	mv_poly temp;
	mv_poly q_k;

	if (q.is_zero())
		return;

	//  terms of p are not allowed to have var. x_k

	for (act_term2 = n.first_term; act_term2 != NULL;
	     act_term2 = act_term2->next)
		if ((act_term2->term).has_var_k(k)) {
			lidia_error_handler("mv_poly", "void substitute(mv_poly & q, "
					    "const mv_poly & n, const mv_poly &n, "
					    "lidia_size_t k)"
					    "::n includes var k.");
			return;
		}

	for (act_term2 = d.first_term; act_term2 != NULL;
	     act_term2 = act_term2->next)
		if ((act_term2->term).has_var_k(k)) {
			lidia_error_handler("mv_poly", "void substitute(mv_poly & q, "
					    "const mv_poly & n, const mv_poly &n, "
					    "lidia_size_t k)"
					    "::d includes var k.");
			return;
		}

	q.split_x_k(q_k, k); // splits q into all term without X_k (stored in
	q_k *= n; // new q) and terms having X_k (in q_k, with
	q *= d; // deleting X_k from all terms)
	q += q_k; // ==> q_k * X_k + q = q_k * (n/d) + q = 0
	// ==> q_k * n + d * q = 0

}



void substitute(mv_poly & qn, mv_poly & qd,
		const mv_poly & n, const mv_poly & d,
		lidia_size_t k)
{
	const mv_list_entry* act_term2;
	mv_poly q_kn, q_kd;

	if (qn.is_zero())
		return;

	//  terms of p are not allowed to have var. x_k

	for (act_term2 = n.first_term; act_term2 != NULL;
	     act_term2 = act_term2->next)
		if ((act_term2->term).has_var_k(k)) {
			lidia_error_handler("mv_poly", "void substitute(mv_poly & q, "
					    "const mv_poly & n, const mv_poly &n, "
					    "lidia_size_t k)"
					    "::n includes var k.");
			return;
		}

	for (act_term2 = d.first_term; act_term2 != NULL;
	     act_term2 = act_term2->next)
		if ((act_term2->term).has_var_k(k)) {
			lidia_error_handler("mv_poly", "void substitute(mv_poly & q, "
					    "const mv_poly & n, const mv_poly &n, "
					    "lidia_size_t k)"
					    "::d includes var k.");
			return;
		}

	qd.split_x_k(q_kd, k); // splits q into all term without X_k (stored in
	qn.split_x_k(q_kn, k);
	q_kn *= n;
	q_kd *= n;
	qd *= d;
	qn *= d;
	qn += q_kn;
	qd += q_kd;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
