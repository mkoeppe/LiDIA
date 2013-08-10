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
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_HASH_TABLE_CC_GUARD_
#define LIDIA_HASH_TABLE_CC_GUARD_



#ifndef LIDIA_HASH_TABLE_H_GUARD_
# include	"LiDIA/hash_table.h"
#endif
#include	<cstring>



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// hentry destructor:
//

template< class T >
hentry< T >::~hentry ()
{
	if (this->item)
		delete this->item;
}



//
// constructor:
//	- set number of buckets and number of elements to 0
//	- set all pointers to NULL
//

template< class T >
hash_table< T >::hash_table ()
{
	debug_handler("hash_table", "hash_table");

	this->size = 0;
	this->curr_size = 0;
	this->buckets = NULL;
	this->last_one = NULL;
	this->key = NULL;
}



//
// destructor:
//	if the table has been allocated, delete each bucket with the function
//	empty(), and then delete the list of buckets
//

template< class T >
hash_table< T >::~hash_table ()
{
	debug_handler("hash_table", "~hash_table");

	if (this->buckets) {
		this->empty();
		delete [] this->buckets;
		this->buckets = NULL;
		this->last_one = NULL;
		this->key = NULL;
	}
}



//
// hash_table < T >::assign()
//
// Task:
//	make a copy of an existing hash table
//

template< class T >
void
hash_table< T >::assign (const hash_table< T > & old_table)
{
	debug_handler("hash_table", "assign");

	lidia_size_t i;
	hentry< T > *ptr, *newone;
	T *tptr;

	// initialize new table
	initialize(old_table.size);
	this->key = old_table.key;

	// perform linked list traversal, and hash each element into the new table
	newone = NULL;
	for (i = 0; i < old_table.size; ++i) {
		ptr = old_table.buckets[i];
		while (ptr) {
			tptr = ptr->item;
			hash(*tptr);
			ptr = ptr->next;
		}
	}
}



//
// operator =
//
// Task:
//	make a copy of an existing hash table
//

template< class T >
hash_table< T > &
hash_table< T >::operator = (const hash_table< T > & old_table)
{
	debug_handler("hash_table", "operator = ");

	this->assign(old_table);
	return *this;
}



//
// hash_table < T >::initialize()
//
// Task:
//	set the number of buckets and allocate storage for the array of	buckets //

template< class T >
void
hash_table< T >::initialize (const lidia_size_t table_size)
{
	debug_handler("hash_table", "initialize");

	bigint temp;
	long lsize;

	if (this->buckets) {
		this->empty();
		delete [] this->buckets;
	}

	// number of buckets should be a prime
	temp = next_prime(bigint(table_size-1));
	temp.longify(lsize);
	this->size = static_cast<lidia_size_t>(lsize);

	this->buckets = new hentry< T > *[lsize];
	memory_handler(this->buckets, "hash_table::initialize", "allocating buckets");

	memset(this->buckets, '\0', static_cast<int>(lsize) * sizeof(this->buckets[0]));
}



//
// hash_table < T >::set_key_function()
//
// Task:
//	set the key function for the hash table elements
//

template< class T >
void
hash_table< T >::set_key_function (bigint (*get_key) (const T &))
{
	debug_handler("hash_table", "set_key_function");

	this->key = get_key;
}



//
// hash_table < T >::no_of_buckets()
//
// Task:
//	returns the number of buckets in the table
//

template< class T >
lidia_size_t
hash_table< T >::no_of_buckets () const
{
	debug_handler("hash_table", "no_of_buckets");

	return this->size;
}



//
// hash_table < T >::no_of_elements()
//
// Task:
//	returns the number of elements currently in the table
//

template< class T >
lidia_size_t
hash_table< T >::no_of_elements () const
{
	debug_handler("hash_table", "no_of_elements");

	return this->curr_size;
}



//
// hash_table < T >::remove()
//
// Task:
//	removes G from the table, if it is there.
//

template< class T >
void
hash_table< T >::remove (const T & G)
{
	debug_handler("hash_table", "remove");

	T *target, *tptr;
	hentry< T > *ptr, *pptr, *nptr;
	lidia_size_t i, j;

	// check whether G is in the table (linked list search)
	target = NULL;
	i = static_cast<lidia_size_t>(remainder(this->key(G), this->size));
	if (i < 0)
		i += this->size;
	ptr = this->buckets[i];
	while ((ptr) && (!target)) {
		tptr = ptr->item;
		if (G == *tptr)
			target = tptr;
		else
			ptr = ptr->next;
	}

	// if so, remove it (linked list delete)
	if (target) {
		pptr = ptr->prev;
		nptr = ptr->next;

		if (pptr) {
			pptr->next = nptr;
			if (nptr)
				nptr->prev = pptr;
		}
		else {
			// G was the first element in the list - handle specially
			tptr = ptr->item;
			j = static_cast<lidia_size_t>(remainder(this->key(*tptr), this->size));
			if (j < 0)
				j += this->size;
			this->buckets[j] = nptr;
			if (nptr)
				nptr->prev = NULL;
		}

		delete ptr;
		--curr_size;
	}
}



//
// hash_table < T >::empty()
//
// Task:
//	delete all the elements in the table.  At the end, the number of
//	buckets is the same, but they will all be empty.
//

template< class T >
void
hash_table< T >::empty ()
{
	debug_handler("hash_table", "empty");

	lidia_size_t i;
	hentry< T > *ptr, *nptr;

	// execute a linked list deleta on each bucket
	for (i = 0; i < this->size; ++i) {
		ptr = this->buckets[i];
		while (ptr) {
			nptr = ptr->next;
			delete ptr;
			ptr = nptr;
		}
		this->buckets[i] = NULL;
	}

	this->curr_size = 0;
	this->last_one = NULL;
}



//
// hash_table < T >::last_entry()
//
// Task:
//	returns a constant reference to the most recent element inserted in
//	the hash table
//
// Conditions:
//	there must be at least one element in the list
//

template< class T >
const T
hash_table< T >::last_entry () const
{
	debug_handler("hash_table", "last_entry");

	if (!last_one) {
		lidia_error_handler("hash_table", "last_entry - table is empty");
		return T();
	}

	return *last_one;
}


//
// hash_table < T >::hash()
//
// Task:
//	insert G into the hash table
//

template< class T >
void
hash_table< T >::hash (const T & G)
{
	debug_handler("hash_table", "hash");

	lidia_size_t i;
	hentry< T > *ptr, *nptr, *newone;
	T *newT;

	// allocate new list element and a copy of G
	newone = new hentry< T >;
	memory_handler(newone, "hash_table::hash", "allocating new list element");

	newT = new T;
	memory_handler(newT, "hash_table::hash", "allocating copy of item");

	(*newT) = G;
	newone->item = newT;

	++curr_size;
	this->last_one = newT;

	// compute correct bucket and append G to the end of the linked list
	i = static_cast<lidia_size_t>(remainder(this->key(G), this->size));
	if (i < 0)
		i += this->size;
	ptr = this->buckets[i];
	if (ptr) {
		while (ptr) {
			nptr = ptr;
			ptr = nptr->next;
		}
		nptr->next = newone;
		newone->prev = nptr;
	}
	else
		this->buckets[i] = newone;
}



//
// hash_table < T >::search()
//
// Task:
//	returns a pointer to the first occurance of G in the hash table.  If
//	G is not in the hash table, the NULL pointer is returned.
//

template< class T >
T *
hash_table< T >::search (const T & G) const
{
	debug_handler("hash_table", "search");

	lidia_size_t i;
	T *target, *tptr;
	hentry< T > *ptr;

	// compute correct bucket and perform linked list search
	i = static_cast<lidia_size_t>(remainder(this->key(G), this->size));
	if (i < 0)
		i += this->size;
	target = NULL;
	ptr = this->buckets[i];
	while ((ptr) && (!target)) {
		tptr = ptr->item;
		if (G == (*tptr))
			target = tptr;
		else
			ptr = ptr->next;
	}

	return target;
}



//
// hash_table < T >::get_bucket()
//
// Task:
//      returns a pointer to the bucket corresponding to the element G.
//      This pointer is essentially the head of a simply linked list.
//      Note that this function does not test if the bucket is full or
//      empty - this is the user's responsibility.
//

template< class T >
hentry< T > *
hash_table< T >::get_bucket (const T & G) const
{
	debug_handler("hash_table", "get_bucket");

	lidia_size_t i;

	// compute correct bucket index
	i = static_cast<lidia_size_t>(remainder(this->key(G), this->size));
	if (i < 0)
		i += this->size;

	return this->buckets[i];
}



//
// hash_table < T >::read()
//
// Task:
//	read in a hash table from the std::istream in.
//
// Conditions:
//	input must be the number of buckets on one line, followed by the
//	number of elements to insert into the list on a new line, and finally
//	n instances of type T.
//

template< class T >
void
hash_table< T >::read (std::istream & in)
{
	debug_handler("hash_table", "read");

	lidia_size_t new_size, num, i;
	T new_item;

	in >> new_size;
	initialize(static_cast<long>(new_size));

	in >> num;
	for (i = 0; i < num; ++i) {
		in >> new_item;
		hash(new_item);
	}
}



//
// hash_table < T >::print()
//
// Task:
//	output a hash table to the std::ostream out.
//

template< class T >
void
hash_table< T >::print (std::ostream & out) const
{
	debug_handler("hash_table", "print");

	lidia_size_t i;
	hentry< T > *ptr;
	T *tptr;

	for (i = 0; i < this->size; ++i) {
		if (this->buckets[i]) {
			out << "[" << i << "]";
			ptr = this->buckets[i];
			while (ptr) {
				tptr = ptr->item;
				out << " : " << (*tptr);
				ptr = ptr->next;
			}
			out << std::endl;
		}
	}
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_HASH_TABLE_CC_GUARD_
