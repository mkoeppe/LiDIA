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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/random_generator.h"
#ifdef HAVE_RANDOM
# include	<cstdlib>
#else
# include	"LiDIA/random.h"
#endif
#include	"LiDIA/timer.h"
#include	<unistd.h>
#ifdef HAVE_POSIX_TIMES
# include	<sys/times.h>
#endif
#ifdef HAVE_POSIX_TIME
# include	<time.h>
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



bool random_generator::initialized = false;



void
random_generator::seed ()
{
	unsigned int	s;
#ifdef HAVE_POSIX_TIMES
	struct tms	t_info;
#endif


#ifdef HAVE_POSIX_TIMES
	s ^= static_cast<unsigned int>(times(&t_info));
#endif
	s ^= static_cast<unsigned int>(getpid());
	s ^= static_cast<unsigned int>(getppid());
#ifdef HAVE_POSIX_TIME
	s ^= static_cast<unsigned int>(time(NULL));
#endif
	
	srandom(s);
	initialized = true;
}



void
random_generator::seed (unsigned int s)
{
	srandom(s);
	initialized = true;
}



//
// random number production via operator>>
//

random_generator &
operator >> (random_generator & rg, int & i)
{
	i = static_cast<int>(random());
	return rg;
}



random_generator &
operator >> (random_generator & rg, unsigned int & i)
{
	i = static_cast<unsigned int>(random());
	return rg;
}



random_generator &
operator >> (random_generator & rg, long & l)
{
	l = static_cast<long>(random());
	return rg;
}



random_generator &
operator >> (random_generator & rg, unsigned long & l)
{
	l = static_cast<unsigned long>(random());
	return rg;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
