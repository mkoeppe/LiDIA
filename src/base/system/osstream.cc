// -*-C++-*-

//=============================================================================
//
//     This file is part of LiDIA --- a library for computational number theory
//
//     Copyright (c) 1994--2002 the LiDIA Group.  All rights reserved.
//
//     See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//-----------------------------------------------------------------------------
//
//     $Id$
//
//     Author  : Christoph Ludwig
//     Changes : See CVS log
//
//=============================================================================

#include "LiDIA/osstream.h"

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

std::string extractString(osstream& oss) {
#if !defined(OSSTREAM_USE_STRSTREAM)
  return oss.str();
#else
  oss << std::ends;                       // add '\0' (implies str() != 0)
  oss.seekp(-1, std::ios::cur);           // restore old write position
  std::string result = oss.str();         // implicitly freezes buffer
  oss.freeze(false);                      // unfreeze buffer
  return result;
#endif
}
 
#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
   
    
