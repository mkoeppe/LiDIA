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
//	Author	: Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_FILE_VECTOR_H_GUARD_
#define LIDIA_FILE_VECTOR_H_GUARD_



#ifndef LIDIA_BASE_VECTOR_H_GUARD_
# include	"LiDIA/base_vector.h"
#endif
#ifndef LIDIA_FILE_IO_CLASS_H_GUARD_
# include	"LiDIA/base/file_io_class.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T >
class file_vector : public virtual base_vector< T >
{

	//
	// constructors
	//

public:

	file_vector():
		base_vector< T > () {}

	file_vector(vector_flags & md):
		base_vector< T > (md) {}
	file_vector(lidia_size_t all):
		base_vector< T > (all) {}

	file_vector(lidia_size_t all, vector_flags & md):
		base_vector< T > (all, md) {}
	file_vector(lidia_size_t all, lidia_size_t len):
		base_vector< T > (all, len) {}
	file_vector(lidia_size_t all, lidia_size_t len, vector_flags & md):
		base_vector< T > (all, len, md) {}
	file_vector(const base_vector< T > &v):
		base_vector< T > (v) {}
	file_vector(const base_vector< T > &v, vector_flags & md):
		base_vector< T > (v, md) {}
	file_vector(const T *v, lidia_size_t len):
		base_vector< T > (v, len) {}
	file_vector(const T *v, lidia_size_t len, vector_flags & md):
		base_vector< T > (v, len, md) {}

	//
	// destructor
	//

public:

	~file_vector() {}

	//
	// reading from / writing to file
	//

public:

#ifdef C_STDIO
	void write_to_file(FILE *fp) const; // writes the vector to file in binary-format
	void read_from_file(FILE *fp); // reads a vector from file in binary-format

	void append_b(FILE *fp) const; // appends the vector to file in binary format
	void append_b(FILE *fp, lidia_size_t n) const; // appends the (n+1)-st element of the vector to file

	void print_to_file(FILE *fp) const; // writes the vector to file in ASCII-format
	void scan_from_file(FILE *fp); // reads a vector from file in ASCII-format

	void append_a(FILE *fp) const; // appends the vector to file in ASCII-format
	void append_a(FILE *fp, lidia_size_t n) const; // appends the (n+1)-st element of the vector to file
#endif	// C_STDIO

};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/file_vector.cc"
#endif



#endif	// LIDIA_FILE_VECTOR_H_GUARD_
