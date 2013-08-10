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


#ifndef LIDIA_INDEXED_HASH_TABLE_CC_GUARD_
#define LIDIA_INDEXED_HASH_TABLE_CC_GUARD_



#ifndef LIDIA_INDEXED_HASH_TABLE_H_GUARD_
# include	"LiDIA/indexed_hash_table.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// constructor:
//	- execute hash_table constructor
//	- initialize list of pointers (expand mode, 0 elements)
//	- set default output style (hash table mode)
//

template< class T >
indexed_hash_table< T >::indexed_hash_table ()
	: hash_table< T > ()
{
	debug_handler("indexed_hash_table", "indexed_hash_table");

	this->IDX = new hentry< T > *[1];
	this->allocated = 1;
	this->outstyle = 0;
}



//
// destructor:
//	reset the list of pointers
//

template< class T >
indexed_hash_table< T >::~indexed_hash_table ()
{
	debug_handler("indexed_hash_table", "~indexed_hash_table");

	if (this->IDX)
		delete [] this->IDX;
}



//
// indexed_hash_table < T >::assign()
//
// Task:
//    make a copy of an existing hash table
//

template< class T >
void
indexed_hash_table< T >::assign (const indexed_hash_table< T > & old_table)
{
	debug_handler("indexed_hash_table", "assign");

	register lidia_size_t i;

	// initialize new table
	this->initialize(old_table.size);
	this->key = old_table.key;

	// insert every element of the old table into the new table
	for (i = 0; i < old_table.curr_size; ++i)
		hash(old_table[i]);
}



//
// operator =
//
// Task:
//      make a copy of an existing hash table
//

template< class T >
indexed_hash_table< T > &
indexed_hash_table< T >::operator = (const indexed_hash_table< T > & old_table)
{
	debug_handler("indexed_hash_table", "operator = ");

	this->assign(old_table);
	return *this;
}



//
// operator []
//
// Task:
//	return the ith element (sequentially)
//
// Conditions:
//	the index i must satisfy 0 <= i < curr_size
//

template< class T >
const T
indexed_hash_table< T >:: operator [] (lidia_size_t i) const
{
	debug_handler("indexed_hash_table", "operator []");

	hentry< T > *ptr;
	T *tptr;

	if ((i< 0) || (i >= this->curr_size)) {
		lidia_error_handler("indexed_hash_table", "operator [] - i >= current size or < 0");
		return T();
	}

	ptr = this->IDX[i];
	tptr = ptr->item;

	return *tptr;
}



//
// indexed_hash_table < T >::member()
//
// Task:
//	return the ith element (sequentially)
//
// Conditions:
//	the index i must satisfy 0 <= i < curr_size
//

template< class T >
const T
indexed_hash_table< T >::member (lidia_size_t i) const
{
	debug_handler("indexed_hash_table", "member");

	hentry< T > *ptr;
	T *tptr;

	if ((i< 0) || (i >= this->curr_size)) {
		lidia_error_handler("indexed_hash_table", "member - i >= current size or < 0");
		return T();
	}

	ptr = this->IDX[i];
	tptr = ptr->item;

	return *tptr;
}



//
// indexed_hash_table < T >::remove()
//
// Task:
//      removes G from the table, if it is there.
//

template< class T >
void
indexed_hash_table< T >::remove (const T & G)
{
	debug_handler("indexed_hash_table", "remove");

	register lidia_size_t j;
	hentry< T > *ptr;
	T *target, *tptr;

	// check whether G is in the table
	target = search(G);

	// if so, delete it
	if (target) {
		// compute j, its sequential index
		for (j = 0; j < this->curr_size; ++j) {
			ptr = this->IDX[j];
			tptr = ptr->item;
			if (tptr == target)
				break;
		}

		// remove element j
		remove_from(j);
	}
}



//
// indexed_hash_table < T >::remove_from()
//
// Task:
//	remove the item with sequential index i from the table
//
// Conditions:
//      the index i must satisfy 0 <= i < curr_size
//

template< class T >
void
indexed_hash_table< T >::remove_from (const lidia_size_t i)
{
	debug_handler("indexed_hash_table", "remove_from");

	hentry< T > *ptr, *pptr, *nptr;
	T *tptr;
	register lidia_size_t j;

	if ((i< 0) || (i >= this->curr_size)) {
		lidia_error_handler("indexed_hash_table", "remove_from - i >= current size or < 0");
		return;
	}

	ptr = this->IDX[i];
	pptr = ptr->prev;
	nptr = ptr->next;

	if (pptr) {
		pptr->next = nptr;
		if (nptr)
			nptr->prev = pptr;
		delete ptr;
	}
	else {
		// G was the first element in the list - handle specially
		tptr = ptr->item;
		j = static_cast<lidia_size_t>(remainder(this->key(*tptr), this->size));
		if (j < 0)
			j += this->size;
		delete this->buckets[j];
		this->buckets[j] = nptr;
		if (nptr)
			nptr->prev = NULL;
	}

	--this->curr_size;
	for (j = i; j < this->curr_size; ++j)
		this->IDX[j] = this->IDX[j+1];
	this->IDX[j] = (hentry< T > *) NULL;
}



//
// indexed_hash_table < T >::empty()
//
// Task:
//      delete all the elements in the table.  At the end, the number of
//      buckets is the same, but they will all be empty.
//

template< class T >
void
indexed_hash_table< T >::empty ()
{
	debug_handler("indexed_hash_table", "empty");

	register lidia_size_t i, end;

	end = this->curr_size;
	for (i = end-1; i >= 0; --i)
		remove_from(i);
}



//
// indexed_hash_table < T >::hash()
//
// Task:
//      insert G into the hash table
//

template< class T >
void
indexed_hash_table< T >::hash (const T & G)
{
	debug_handler("indexed_hash_table", "hash");

	register lidia_size_t i;
	hentry< T > *ptr, *nptr, *newone;
	T *newT;

	// allocate new list element and a copy of G
	newone = new hentry< T >;
	memory_handler(newone, "indexed_hash_table::hash", "allocating new list element");

	newT = new T;
	memory_handler(newone, "indexed_hash_table::hash", "allocating new list element");

	(*newT) = G;
	newone->item = newT;
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

	// append the pointer to the list of pointers
	if (this->curr_size == this->allocated) {
		// allocated more storage for the sequential list

		hentry< T > **tmp = new hentry< T > *[this->allocated << 1];
		memory_handler(tmp, "indexed_hash_table::hash", "allocating more sequential list elemets");

		for (i = 0; i < this->allocated; ++i)
			tmp[i] = this->IDX[i];

		if (this->IDX)
			delete [] this->IDX;
		this->IDX = tmp;
		this->allocated <<= 1;
	}

	this->IDX[this->curr_size] = newone;
	++this->curr_size;
}



//
// indexed_hash_table < T >::output_style()
//
// Task:
//	set the style of output (0 = list, 1 = hash table)
//

template< class T >
void
indexed_hash_table< T >::output_style (int style)
{
	debug_handler("indexed_hash_table", "output_style");

	if (style == 1)
		this->outstyle = 1;
	else
		this->outstyle = 0;
}



//
// indexed_hash_table < T >::read()
//
// Task:
//      read in a hash table from the std::istream in.
//
// Conditions:
//      input must be the number of buckets on one line, followed by the
//      number of elements to insert into the list on a new line, and finally
//      n instances of type T.
//

template< class T >
void
indexed_hash_table< T >::read (std::istream & in)
{
	debug_handler("indexed_hash_table", "read");

	register lidia_size_t i;
	lidia_size_t new_size, num;
	T new_item;

	in >> new_size;
	this->initialize(static_cast<long>(new_size));

	in >> num;
	for (i = 0; i < num; ++i) {
		in >> new_item;
		hash(new_item);
	}
}



//
// indexed_hash_table < T >::print()
//
// Task:
//      output a hash table to the std::ostream out.
//

template< class T >
void
indexed_hash_table< T >::print (std::ostream & out) const
{
	debug_handler("indexed_hash_table", "print");

	register lidia_size_t i;
	hentry< T > *ptr;
	T *tptr;

	if (this->outstyle == 0) {
		// list style output
		for (i = 0; i < this->curr_size; ++i) {
			ptr = this->IDX[i];
			tptr = ptr->item;
			out << i << ": " << (*tptr) << std::endl;
		}
	}
	else {
		// hash table style output
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
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_INDEXED_HASH_TABLE_CC_GUARD_
