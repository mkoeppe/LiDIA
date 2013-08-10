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
//	Author	: Thomas Papanikolaou (TP)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_PAIR_H_GUARD_
#define LIDIA_PAIR_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include	"LiDIA/LiDIA.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



template< class T1, class T2 > class pair
{
private:

	T1 l;
	T2 r;


public:

	pair();
	pair(T1 &pl, T2 &pr);
	pair(const pair< T1, T2 > &p);
	~pair();


	void assign(const T1 & pl, const T2 & pr);
	void assign (const pair< T1, T2 > &p);
	pair< T1, T2 > & operator = (const pair< T1, T2 > &p);


	const T1 & left() const;
	const T2 & right() const;

	T1 & left();
	T2 & right();


	void read(std::istream & in = std::cin);
	void write(std::ostream & out = std::cout) const;


	void swap(pair< T1, T2 > & p);
};



//
// c'tors and d'tor
//

template< class T1, class T2 >
inline
pair< T1, T2 >::pair ()
{
}



template< class T1, class T2 >
inline
pair< T1, T2 >::pair (T1 & pl, T2 & pr)
	: l(pl),
	  r(pr)
{
}



template< class T1, class T2 >
inline
pair< T1, T2 >::pair (const pair< T1, T2 > & p)
	: l(p.l),
	  r(p.r)
{
}



template< class T1, class T2 >
inline
pair< T1, T2 >::~pair ()
{
}



//
// assigners
//

template< class T1, class T2 >
inline void
pair< T1, T2 >::assign (const T1 & pl, const T2 & pr)
{
	l = pl;
	r = pr;
}



template< class T1, class T2 >
inline void
pair< T1, T2 >::assign (const pair< T1, T2 > & p)
{
	if (&p != this) {
		l = p.l;
		r = p.r;
	}
}



template< class T1, class T2 >
inline pair< T1, T2 > &
pair< T1, T2 >::operator = (const pair< T1, T2 > & p)
{
	assign(p);
	return *this;
}



//
// accessors
//

template< class T1, class T2 >
inline const T1 &
pair< T1, T2 >::left () const
{
	return l;
}



template< class T1, class T2 >
inline const T2 &
pair< T1, T2 >::right () const
{
	return r;
}



template< class T1, class T2 >
inline T1 &
pair< T1, T2 >::left ()
{
	return l;
}



template< class T1, class T2 >
inline T2 &
pair< T1, T2 >::right ()
{
	return r;
}



template< class T1, class T2 >
inline std::istream &
operator >> (std::istream & in, pair< T1, T2 > & p)
{
	p.read(in);
	return in;
}



template< class T1, class T2 >
inline std::ostream &
operator << (std::ostream & out, const pair< T1, T2 > & p)
{
	p.write(out);
	return out;
}



//
// comparators
//

template< class T1, class T2 >
inline bool
operator == (const pair< T1, T2 > & p, const pair< T1, T2 > & q)
{
	return ((&p == &q) || ((q.left() == p.left()) && (q.right() == p.right())));
}



template< class T1, class T2 >
inline bool
operator != (const pair< T1, T2 > & p, const pair< T1, T2 > & q)
{
	return ((&p != &q) && ((q.left() != p.left()) || (q.right() != p.right())));
}



template< class T1, class T2 >
inline bool
operator > (const pair< T1, T2 > & p, const pair< T1, T2 > & q)
{
	return (p.left() > q.left()) || ((p.left() == q.left()) && (p.right() > q.right()));
}



template< class T1, class T2 >
inline bool
operator >= (const pair< T1, T2 > & p, const pair< T1, T2 > & q)
{
	return (p.left() >= q.left());
}



template< class T1, class T2 >
inline bool
operator< (const pair < T1, T2 > & p, const pair< T1, T2 > & q)
{
	return (p.left() < q.left()) || ((p.left() == q.left()) && (p.right() < q.right()));
}



template< class T1, class T2 >
inline bool
operator <= (const pair< T1, T2 > & p, const pair< T1, T2 > & q)
{
	return (p.left() <= q.left());
}



template< class T1, class T2 >
inline void
pair< T1, T2 >::swap (pair< T1, T2 > & p)
{
	LiDIA::swap(l, p.l);
	LiDIA::swap(r, p.r);
}



template< class T1, class T2 >
inline void
swap(pair< T1, T2 > & p, pair< T1, T2 > & q)
{
	p.swap(q);
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include	"LiDIA/pair.cc"
#endif



#endif	// LIDIA_PAIR_H_GUARD_
