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
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_HASH_TABLE_H_GUARD_
#define LIDIA_HASH_TABLE_H_GUARD_



#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// Class: hentry < T >
//
// This class represents one element in the hash table.  It is essentially
//	an element in a linked list.
//

template< class T >
class hentry
{
protected:

public:

	T *item; // pointer to the data item
	hentry< T > *next, *prev; // pointers to next and previous list elements

	hentry()
	{
		item = NULL;
		next = NULL;
		prev = NULL;
	};
	~hentry();
};



//
// Class: hash_table < T >
//
// This class represents a hash table, whose elements are of some type T.
//	The collision resolution scheme is bucketing, and each bucket is a
//	linked list.
// The hash function is simply the key value of the data item modulo the
//	number of buckets in the table.  The user must define a key function
//	corresponding to the type T, and initialize it with the function
//	set_key_function().
//

template< class T >
class hash_table
{
protected:

	lidia_size_t size; // number of buckets in the hash table
	lidia_size_t curr_size; // current number of elemets
	hentry< T > **buckets; // array of data buckets
	T *last_one; // pointer to most recent entry

	bigint (*key) (const T & G); // key function (must be set by user)


public:

	hash_table();
	~hash_table();

	void assign(const hash_table< T > & old_table);
	hash_table< T > & operator = (const hash_table< T > & old_table);

	void initialize(const lidia_size_t table_size);
	void set_key_function(bigint (*get_key) (const T &));

	lidia_size_t no_of_buckets() const;
	lidia_size_t no_of_elements() const;

	void remove(const T & G);
	void empty();
	const T last_entry() const;
	void hash(const T & G);
	T * search(const T & G) const;

#ifndef HEADBANGER
	hentry< T > * get_bucket(const T & G) const;
#endif

	void read  (std::istream & in);
	void print (std::ostream & out) const;

};



template< class T >
inline std::istream &
operator >> (std::istream & in, hash_table< T > & HT)
{
	HT.read(in);
	return in;
}



template< class T >
inline std::ostream &
operator << (std::ostream & out, const hash_table< T > & HT)
{
	HT.print(out);
	return out;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/hash_table.cc"
#endif



#endif	// LIDIA_HASH_TABLE_H_GUARD_
