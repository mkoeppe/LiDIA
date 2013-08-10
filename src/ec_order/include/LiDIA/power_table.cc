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


#ifndef LIDIA_POWER_TABLE_CC_GUARD_
#define LIDIA_POWER_TABLE_CC_GUARD_



#ifndef LIDIA_POWER_TABLE_H_GUARD_
# include	"LiDIA/power_table.h"
#endif



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// constructors and destructor
//

template <class T>
power_tab< T >::power_tab()
{
	initialized = false;
	N = 0;
}



template <class T>
power_tab< T >::~power_tab()
{
	if (initialized)
		delete[] f_powers;
}



//
// member functions
//

template <class T>
void power_tab< T >::initialize(const T & f, lidia_size_t n)
{
	int j = 2, i;

	f_powers = new T[n+1];

	f_powers[0].assign_one();
	f_powers[1] = f;

	while (j <= n) {
		multiply(f_powers[j], f_powers[1], f_powers[j-1]);

		for(i = (j << 1); i <= n; i <<= 1) {
			square(f_powers[i], f_powers[i>>1]);
		}
		if (j & 1)
			j = j + 2;
		else
			j++;
	}
	N = n;
	initialized = true;
}



template <class T>
void power_tab< T >::get_power(T &f, lidia_size_t i) const
{
	if (initialized && i >= 0 && i <= N)
		f = f_powers[i];
	else
		lidia_error_handler ("power_tab::get_power(T &f, lidia_size_t &i)",
				     "Either not initialized or i out of range.");
}



template <class T>
T power_tab< T >::get_power(lidia_size_t i) const
{
	if (initialized && i >= 0 && i <= N)
		return f_powers[i];
	else {
		lidia_error_handler ("power_tab::get_power(lidia_size_t &i)",
				     "Either not initialized or i out of range.");
		return T();
	}
}



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}	// end of namespace LiDIA
# endif
#endif



#endif	// LIDIA_POWER_TABLE_CC_GUARD_
