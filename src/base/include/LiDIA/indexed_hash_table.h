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


#ifndef LIDIA_INDEXED_HASH_TABLE_H_GUARD_
#define LIDIA_INDEXED_HASH_TABLE_H_GUARD_



#ifndef LIDIA_HASH_TABLE_H_GUARD_
# include	"LiDIA/hash_table.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// Class: indexed_hash_table < T >
//
// This class represents a hash table, whose elemets are of some type T.  The
//	elements can be accessed either as in a regular hash table or
//	sequentially, as in a list.
//

template< class T >
class indexed_hash_table : public hash_table< T >
{
private:

	lidia_size_t allocated; // size of array
	hentry< T > **IDX; // array of pointers to list elements
	int outstyle; // output style (0 = list, 1 = hash table)


public:

	indexed_hash_table();
	~indexed_hash_table();

	void assign(const indexed_hash_table< T > & old_table);
	indexed_hash_table< T > & operator = (const indexed_hash_table< T > & old_table);

	const T operator[] (lidia_size_t i) const;
	const T member(lidia_size_t i) const;

	void remove(const T & G);
	void remove_from(const lidia_size_t i);
	void empty();
	void hash(const T & G);

	void output_style(int style);
	void read  (std::istream & in);
	void print (std::ostream & out) const;

};



template< class T >
inline std::istream &
operator >> (std::istream & in, indexed_hash_table< T > & HT)
{
	HT.read(in);
	return (in);
}



template< class T >
inline std::ostream &
operator << (std::ostream & out, const indexed_hash_table< T > & HT)
{
	HT.print(out);
	return (out);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/indexed_hash_table.cc"
#endif



#endif	// LIDIA_INDEXED_HASH_TABLE_H_GUARD_
