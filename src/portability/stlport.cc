/*
 *
 *	This file is part of:
 *
 *		LiDIA --- a library for computational number theory
 *
 *	Copyright (c) 1994--2002 the LiDIA Group.  All rights reserved.
 *
 *	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
 *
 *	$Id$
 *
 */

/*
 * Since LiDIA and its apps are compiled with -fno-implicit-templates,
 * some templates from the standard template library (STL) will not be
 * instantiated automatically.  To compensate for that, this file contains
 * some explicit instantiations.  Note that these depend on the particular
 * STL implementation.
 */

#include <iostream>
#ifdef _STLP_IOSTREAM		// got STLport's <iostream>

// Template instantiations supplementing STLport 4.5.3 / gcc 2.95.2

_STLP_BEGIN_NAMESPACE

// Level 0 requisites (by LiDIA or its applications)

template
basic_istream< char, char_traits< char > > & _STLP_CALL
ws< char, char_traits< char > >(
	basic_istream< char, char_traits< char > > &
);

template
basic_ostream< char, char_traits< char > > & _STLP_CALL
operator << < char, char_traits< char >, allocator< char > >(
	basic_ostream< char, char_traits< char > > &,
	basic_string< char, char_traits< char >, allocator< char > > const &
);

// Level 1 requisites (by the above)

template
bool _STLP_CALL
__stlp_string_fill< char, char_traits< char > >(
	basic_ostream< char, char_traits< char > > &,
	basic_streambuf< char, char_traits< char > > *,
	size_t
);

_STLP_END_NAMESPACE
#endif

