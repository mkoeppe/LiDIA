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
//	Changes	: see CVS log
//
//==============================================================================================


#ifndef LIDIA_P_VECTOR_CC_GUARD_
#define LIDIA_P_VECTOR_CC_GUARD_



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



inline void
p_vector< bigint >::assign (bigint* c, bigint* a)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		c[counter].assign(a[counter]);
}



inline void
p_vector< bigint >::assign_zero (bigint* a)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		a[counter].assign_zero();
}



inline void
p_vector< bigint >::add (bigint* c, bigint* a, bigint* b)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		LiDIA::add(c[counter], a[counter], b[counter]);
}



inline void
p_vector< bigint >::subtract (bigint* c, bigint* a, bigint* b)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		LiDIA::subtract(c[counter], a[counter], b[counter]);
}



inline void
p_vector< bigint >::scalmul (bigint* c, const bigint d, bigint* a)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		LiDIA::multiply(c[counter], d, a[counter]);
}



inline void
p_vector< bigint >::scalsub (bigint* c, bigint* a,
			     const bigint d, bigint* b)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--) {
		LiDIA::multiply(this->tempT, d, b[counter]);
		LiDIA::subtract(c[counter], a[counter], this->tempT);
	}
}



inline void
p_vector< bigint >::scalprod (bigint& res, bigint* a, bigint * b)
{
	res.assign_zero();
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--) {
		LiDIA::multiply(this->tempT, a[counter], b[counter]);
		LiDIA::add(res, res, this->tempT);
	}
}



inline void
p_vector< bigint >::swap (bigint*& a, bigint*& b)
{
	this->tempPT = a;
	a = b;
	b = this->tempPT;
}



inline void
p_vector< bigfloat >::assign (bigfloat* c, bigfloat* a)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		c[counter].assign(a[counter]);
}



inline void
p_vector< bigfloat >::assign_zero (bigfloat* a)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		a[counter].assign_zero();
}



inline void
p_vector< bigfloat >::add (bigfloat* c, bigfloat* a, bigfloat* b)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		LiDIA::add(c[counter], a[counter], b[counter]);
}



inline void
p_vector< bigfloat >::subtract (bigfloat* c, bigfloat* a, bigfloat* b)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		LiDIA::subtract(c[counter], a[counter], b[counter]);
}



inline void
p_vector< bigfloat >::scalmul (bigfloat* c, const bigfloat d,
			       bigfloat* a)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		LiDIA::multiply(c[counter], d, a[counter]);
}



inline void
p_vector< bigfloat >::scalsub (bigfloat* c, bigfloat* a,
			       const bigfloat d, bigfloat* b)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--) {
		LiDIA::multiply(this->tempT, d, b[counter]);
		LiDIA::subtract(c[counter], a[counter], this->tempT);
	}
}



inline void
p_vector< bigfloat >::scalprod (bigfloat& res, bigfloat* a,
				bigfloat * b)
{
	res.assign_zero();
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--) {
		LiDIA::multiply(this->tempT, a[counter], b[counter]);
		LiDIA::add(res, res, this->tempT);
	}
}



inline void
p_vector< bigfloat >::swap (bigfloat*& a, bigfloat*& b)
{
	this->tempPT = a;
	a = b;
	b = this->tempPT;
}



template< class T >
inline void
p_vector< T >::assign (T* c, T* a)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		c[counter] = a[counter];
}



template< class T >
inline void
p_vector< T >::assign_zero (T* a)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		a[counter] = 0;
}



template< class T >
inline void
p_vector< T >::add (T* c, T* a, T* b)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		c[counter] = a[counter]+b[counter];
}



template< class T >
inline void
p_vector< T >::subtract (T* c, T* a, T* b)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		c[counter] = a[counter]-b[counter];
}



template< class T >
inline void
p_vector< T >::scalmul (T* c, const T d, T* a)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		c[counter] = d*a[counter];
}



template< class T >
inline void
p_vector< T >::scalsub (T* c, T* a, const T d, T* b)
{
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		c[counter] = a[counter]-d*b[counter];
}



template< class T >
inline void
p_vector< T >::scalprod (T& res, T* a, T * b)
{
	res = 0;
	for (lidia_size_t counter = this->vectsize-1; counter >= 0; counter--)
		res += a[counter]*b[counter];
}



template< class T >
inline void
p_vector< T >::swap (T*& a, T*& b)
{
	this->tempPT = a;
	a = b;
	b = this->tempPT;
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_P_VECTOR_CC_GUARD_
