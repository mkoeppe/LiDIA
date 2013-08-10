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
//	Author	:
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_LANCZOS_H_GUARD_
#define LIDIA_LANCZOS_H_GUARD_

#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif

#include <memory>
#include <vector>
#include <cstring>    // declares memcpy()

#ifdef DEBUG
#include <cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// globale lanczos definitions
//

#define WordSize BITS_PER_LONG
void bin_out(const unsigned long);

//
// Bitmask for operations
//
#define One_mask (~(0UL))
#define Bit_mask(k) (1UL<<(WordSize-1-(k)))

class index_list
{
public:
  typedef size_t size_type;
  typedef long value_type;

private:
  size_type length;
  value_type *list;

public:
  index_list(const size_type number)
  {
    if (number == 0)
      lidia_error_handler("lanczos","index_list(size_t)::length "
			  "is non-positive"); 
    length = number;
    list = new value_type [number];
    memset(list, 0, length * sizeof(value_type));
  }
  
  ~index_list()
  {  
    delete[] list;   
  }
  
  value_type get(const size_type pos) const
  {
    if (pos >= length)
      lidia_error_handler("index_list","pos >= length");
    return list[pos];
  }
  
  void put(const size_type pos, const value_type value)
  {
    if (pos >= length)
      lidia_error_handler("index_list","pos >= length");
    list[pos] = value;
  }
  
  size_type number() const
  {
    return length;
  }
  
  void print() const
  {
    printf("\nIndex_List:\t");
    for (size_type i = 0; i < length; i++)
      printf(" [%ld] %ld\n", i, list[i]);
  }
  
  void clear()
  {
    memset(list, 0, length * sizeof(value_type));
  }
};



//************************************************************************
//
// Class Lanczos_Sparse_Vector
// helping class for a sparse matrix
//
//************************************************************************

class lanczos_sparse_vector
{
friend class lanczos_sparse_matrix;

public:
  typedef size_t size_type;
  typedef unsigned long value_type;

private:
  size_type length;
  size_type number_of_entries; // number of non zero entries
  value_type *entries; // array of entries

public:
  lanczos_sparse_vector();
  lanczos_sparse_vector(const size_type len,
			const size_type number);
  lanczos_sparse_vector(const lanczos_sparse_vector & vector);
  
  ~lanczos_sparse_vector();

  const lanczos_sparse_vector &operator = (const lanczos_sparse_vector&v)
  {
    if (&v != this) 
      {
        delete[] entries;

	length = v.length;
	number_of_entries = v.number_of_entries;
	entries = new value_type [number_of_entries];

	memcpy(entries, v.entries, number_of_entries*sizeof(value_type));
      }
    return *this;
  }

  value_type *get_entries() const
  {
    return entries;
  }

  void set_entries(value_type *values, const size_type number,
		   const size_type len)
  {
    if( entries != values ) {
      delete[] entries;
      entries = values;
      number_of_entries = number;
      length = len;
    }
  }

  size_type get_length() const
  {
    return length;
  }

  size_type get_number_of_entries() const
  {
    return number_of_entries;
  }

  void put_number_of_entries(const size_type number)
  {
    number_of_entries = number;
  }
  
  value_type get_entry(const size_type i) const
  {
    if(i >= number_of_entries)
      lidia_error_handler("lanczos_sparse_vector", "out of range");
    return entries[i];
  }

  void put_entry(const size_type i, const value_type value)
  {
    if(i >= number_of_entries || value > length)
      lidia_error_handler("lanczos_sparse_vector", "out of range");
    entries[i] = value;
  }

  void clear()
  {
    memset(entries, 0, number_of_entries * sizeof(value_type));
  }

  void print() const;
  lanczos_sparse_vector &fill_random (long part, int ratio);
};



//***********************************************************************
//
// small dense matrix over GF(2) (Wordsize x Wordsize)
//
//***********************************************************************

class lanczos_small_matrix
{
public:
  typedef size_t size_type;
  typedef unsigned long value_type;

private:
  value_type rows[WordSize];
  // matrix is wordsize x wordsize over GF(2)
  // ulongs are rows

public:
  lanczos_small_matrix() {memset(rows,0,sizeof(rows));};
  ~lanczos_small_matrix() {};

  const lanczos_small_matrix& operator= (const lanczos_small_matrix& matrix)
  {
    memcpy(rows, matrix.rows, WordSize * sizeof(value_type));
  }

  void put_row(const size_type pos, const value_type row)
  {
    if (pos >= WordSize)
      lidia_error_handler("lanczos_small_matrix","put_row::index out"
			  " of range");
    rows[pos] = row;
  }

  value_type get_row(const size_type pos) const
  {
    if (pos >= WordSize)
      lidia_error_handler("lanczos_small_matrix","get_row::index out"
			  " of range");
    return rows[pos];
  }

  void clear()
  {
    memset(rows, 0, WordSize * sizeof(value_type));
  }
  lanczos_small_matrix *mult_right(const lanczos_small_matrix& m) const;
  lanczos_small_matrix *mult_right_transpose(const lanczos_small_matrix& m) const;

  void mult_right_to(const lanczos_small_matrix& m,
		     lanczos_small_matrix& result) const;
  void mult_right_transpose_to(const lanczos_small_matrix& m,
			       lanczos_small_matrix& result) const;
  
  bool is_zero() const
  {
    for (size_type i = 0; i < WordSize; i++)
      if (rows[i] != 0)
	return false;
    return true;
  }

  void eliminate (const value_type entry)
  {
    for (size_type i = 0; i < WordSize; i++)
      put_row(i, get_row(i) & entry);
  }

  void add(const lanczos_small_matrix& matrix)
  {
    for (size_type i = 0; i < WordSize; i++)
      put_row(i, matrix.get_row(i) ^ get_row(i));
  }

  void print() const;
};



//************************************************************************
//
// vector block over GF(2)
//
//************************************************************************

class lanczos_vector_block
{
public:
  typedef size_t size_type;
  typedef unsigned long value_type;

  typedef std::vector< std::vector< size_type> > result_vector_type;

private:
  size_type length;
  value_type *rows;

public:

  lanczos_vector_block(const size_type len)
  {
    if (len == 0)
      lidia_error_handler("lanczos_vector_block","ct::length is zero");
    rows = new value_type [len];
    length = len;
    memset(rows, 0, length * sizeof(value_type));
  }

  lanczos_vector_block(const lanczos_vector_block &v)
  {
    length = v.get_length();
    rows = new value_type [length];
    memcpy(rows, v.rows, length * sizeof(value_type));
  }

  ~lanczos_vector_block()
  {
    delete[] rows;
  }

  size_type get_length() const
  {
    return length;
  }

  void put_row(const size_type pos, const value_type row)
  {
    if(pos >= length)
      lidia_error_handler("lanczos_vector_block", "put_row::out of"
			  " range");
    rows[pos] = row;
  }

  value_type get_row(const size_type pos) const
  {
    if(pos >= length)
      lidia_error_handler("lanczos_vector_block", "get_row::out of"
			  " range");
    return rows[pos];
  }

  void clear()
  {
    memset(rows, 0, length * sizeof(value_type));
  }

  bool is_zero() const
  {
    for (size_type i = 0; i < length; i++)
      if (rows[i] != 0)
	return false;
    return true;
  }

  void print() const;

  void read(result_vector_type const& res);
  result_vector_type result() const;
  lanczos_small_matrix *mult(const lanczos_vector_block& v) const;
  lanczos_vector_block *mult_small (const lanczos_small_matrix& m) const;
  void mult_to(const lanczos_vector_block& v, 
               lanczos_small_matrix& result) const;
  void mult_small_to (const lanczos_small_matrix& m, 
                      lanczos_vector_block& result) const;
  
  void add (const lanczos_vector_block& vector)
  {
    for (size_type i = 0; i < length; i++)
      bitwise_xor(i, vector.get_row(i));
  }

  void bitwise_xor (const size_type row, const value_type value)
  {
    if(row >= length)
      lidia_error_handler("lanczos_vector_block", "bitwise_xor::out of"
			  " range");   
    rows[row] ^= value;
  }

  void bitwise_and (const size_type row, const value_type value)
  {
    if(row >= length)
      lidia_error_handler("lanczos_vector_block", "bitwise_and::out of"
			  " range");   
    rows[row] &= value;
  }

  void bitwise_or (const size_type row, const value_type value)
  {
    if(row >= length)
      lidia_error_handler("lanczos_vector_block", "bitwise_or::out of"
			  " range");   
    rows[row] |= value;
  }

  void eliminate(const value_type entry)
  {
    for (size_type i = 0; i < length; i++)
      bitwise_and(i, entry);
  }
};



//************************************************************************
//
// big sparse matrix over GF(2) for lanczos
//
//************************************************************************

class lanczos_sparse_matrix
{
friend struct preprocess;

public:
  typedef size_t size_type;
  typedef lanczos_sparse_vector::value_type value_type;

private:
  size_type columns; // number of columns
  size_type rows; // number of rows
  lanczos_sparse_vector *entries; // Array of sparse_vectors
  index_list *row_weight;

public:
  lanczos_sparse_matrix(const size_type row, const size_type col);
  lanczos_sparse_matrix(const char *filename);
  ~lanczos_sparse_matrix();
  
  
  size_type number_of_columns() const
  {
    return columns;
  }
  
  size_type number_of_rows() const
  {
    return rows;
  }

  void put_columns(const size_type col)
  {
    columns = col;
  }

  void put_vector(const size_type pos, const lanczos_sparse_vector &v)
  { // check apa, columns or rows??
    entries[pos] = v;
  }
	
  lanczos_sparse_vector& get_vector(const size_type pos)
  {
    return entries[pos];
  }

  const lanczos_sparse_vector& get_vector(const size_type pos) const
  {
    return entries[pos];
  }

  void delete_vector(const size_type pos);
  lanczos_vector_block *mult_vectorblock(const lanczos_vector_block& vector_block) const;
  void mult_vectorblock_to(const lanczos_vector_block& vector_block,
			   lanczos_vector_block& result) const;
  void print() const;
  void write(const char *filename) const;
  void fill_random(const size_type number, const long part, const int ratio);
  long calculate_weight() const;
  void delete_rows(const index_list& row_list);
};



//***************************************************************************
//
// class lanczos
//
//***************************************************************************

class lanczos
{
public:
  typedef size_t size_type;
  typedef unsigned long value_type;

private:
	const lanczos_sparse_matrix& A;
	// vektors
	lanczos_vector_block *V_next;
	lanczos_vector_block *V;
	lanczos_vector_block *V_prev;
	lanczos_vector_block *V_pprev;
	lanczos_vector_block *X;

	lanczos_vector_block *dummy1_vec;
	lanczos_vector_block *dummy2_vec;

	size_type result_rank; // rank of the resultvector


	size_type M, N; //rows and columns

	size_type compute_tau2(const lanczos_small_matrix& mt,
		  	       value_type& S_entry,
			       lanczos_small_matrix& mw_inv) const;

	void compute_result();
	void post_process();
	void compute_D(lanczos_small_matrix& VAAV,
		       const lanczos_small_matrix& VAV,
		       const lanczos_small_matrix& W_inv,
		       const value_type S,
		       lanczos_small_matrix& eta,
		       lanczos_small_matrix& result) const;

	void compute_E(const lanczos_small_matrix& VAV,
		       const lanczos_small_matrix& W_inv_prev,
		       const value_type S,
		       lanczos_small_matrix& result) const;

	void compute_F(const lanczos_small_matrix& VAV_prev,
		       const lanczos_small_matrix& eta_prev,
		       const lanczos_small_matrix& W_inv_prev,
		       const lanczos_small_matrix& W_inv_pprev,
		       const value_type S,
		       lanczos_small_matrix& result) const;

public:
	explicit lanczos (const lanczos_sparse_matrix& matrix);
	~lanczos ();
	void solve ();
	size_type get_result_rank() const
	{
	  return result_rank;
	}

	const lanczos_vector_block& get_result() const
	{
	  return *X;
	}

//G.A. - removed 'const' keyword
	lanczos_vector_block& get_result()
	{
	  return *X;
	}
};



struct preprocess {
        typedef lanczos_sparse_matrix::size_type size_type;
	std::auto_ptr<index_list>
        process(lanczos_sparse_matrix& matrix) const;
};



struct postprocess {
        typedef lanczos_vector_block::size_type size_type;
	std::auto_ptr<lanczos_vector_block>
        process(const lanczos_vector_block& vector,
                const index_list& list) const;
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_LANCZOS_H_GUARD_
