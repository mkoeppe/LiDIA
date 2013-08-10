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
//	Author	: Werner Backes (WB), Thorsten Lauer (TL)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_P_VECTOR_H_GUARD_
#define LIDIA_P_VECTOR_H_GUARD_


#ifndef LIDIA_BIGINT_H_GUARD_
# include	"LiDIA/bigint.h"
#endif
#ifndef LIDIA_BIGFLOAT_H_GUARD_
# include	"LiDIA/bigfloat.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// class p_vector
// ~~~~~~~~~~~~~~
//
// template class for vector operations
// using pointer to the data type
// size of vector is stored in vectsize
//
// requirements for the data - type :
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// arith. Operators : +, -, *
// other Operators  : =
// and integer - zero assignment
// specialized for bigint and bigfloat
//

template< class T > class p_vector;



template<>
class p_vector< bigint >
{
protected :
	bigint*  tempPT;
	bigint   tempT;

public :
	lidia_size_t  vectsize;

	p_vector()
	{
		this->vectsize = 0;
	}
	p_vector(lidia_size_t size)
	{
		this->vectsize = size;
	}
	~p_vector()
	{ }

	void assign(bigint*, bigint*);
	void assign_zero(bigint*);
	void add(bigint*, bigint*, bigint*);
	void subtract(bigint*, bigint*, bigint*);
	void scalmul(bigint*, const bigint, bigint*);
	void scalsub(bigint*, bigint*, const bigint, bigint*);
	void scalprod(bigint&, bigint*, bigint*);
	void swap(bigint*&, bigint*&);
};



//
// Specialization of p_vector for the class bigfloat
//

template<>
class p_vector< bigfloat >
{
protected :
	bigfloat*  tempPT;
	bigfloat   tempT;

public :
	lidia_size_t  vectsize;

	p_vector()
	{
		this->vectsize = 0;
	}
	p_vector(lidia_size_t size)
	{
		this->vectsize = size;
	}
	~p_vector()
	{ }

	void assign(bigfloat*, bigfloat*);
	void assign_zero(bigfloat*);
	void add(bigfloat*, bigfloat*, bigfloat*);
	void subtract(bigfloat*, bigfloat*, bigfloat*);
	void scalmul(bigfloat*, const bigfloat, bigfloat*);
	void scalsub(bigfloat*, bigfloat*, const bigfloat, bigfloat*);
	void scalprod(bigfloat&, bigfloat*, bigfloat*);
	void swap(bigfloat*&, bigfloat*&);
};



//
// Template Code
//
template< class T >
class p_vector
{
protected :
	T*  tempPT;
	T   tempT;

public :
	lidia_size_t  vectsize;

	p_vector()
	{
		this->vectsize = 0;
	}
	p_vector(lidia_size_t size)
	{
		this->vectsize = size;
	}
	~p_vector() { }

	void assign(T*, T*);
	void assign_zero(T*);
	void add(T*, T*, T*);
	void subtract(T*, T*, T*);
	void scalmul(T*, const T, T*);
	void scalsub(T*, T*, const T, T*);
	void scalprod(T&, T*, T*);
	void swap(T*&, T*&);
};



#include	"LiDIA/lattices/p_vector.cc"


//
// class p_vector_SP
// ~~~~~~~~~~~~~~~~~
//
// template class derived from class p_vector
// with the same functionality but allows to
// define a scalarproduct by giving a pointer
// to a function :
// void (*scal_prod_T)(T&, T*, T*, lidia_size_t)
//
template< class T >
class p_vector_SP : public p_vector< T >
{
protected :
	void (*scal_prod_T)(T&, T*, T*, lidia_size_t);
public :
	p_vector_SP():p_vector< T > ()
	{
		this->scal_prod_T = NULL;
	}
	p_vector_SP(lidia_size_t size):p_vector< T > (size)
	{
		this->scal_prod_T = NULL;
	}
	p_vector_SP(void (*spT)(T&, T*, T*, lidia_size_t),
		    lidia_size_t size = 0):p_vector< T > (size)
	{
		this->scal_prod_T = spT;
	}
	~p_vector_SP() { }


	void set_pointer(void (*spT)(T&, T*, T*, lidia_size_t))
	{
		this->scal_prod_T = spT;
	}

	void scalprod(T& res, T* a, T * b)
	{
		(*this->scal_prod_T)(res, a, b, this->vectsize);
	}
};



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_P_VECTOR_H_GUARD_
