//=============================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2002 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//-----------------------------------------------------------------------------
//
//	$Id$
//
//	Author	: Christoph Ludwig (CL)
//	Changes	: See CVS log
//
//=============================================================================

// LiDIA's headers
#include "LiDIA/bigint.h"

#include "LiDIA/param_desc.cc"

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif

  namespace details {

  template class param_desc_impl<bool>;
  template class param_desc_impl<short>;
  template class param_desc_impl<unsigned short>;
  template class param_desc_impl<int>;
  template class param_desc_impl<unsigned int>;
  template class param_desc_impl<long>;
  template class param_desc_impl<unsigned long>;
  template class param_desc_impl<char>;
  template class param_desc_impl<signed char>;
  template class param_desc_impl<unsigned char>;
  template class param_desc_impl<const char*>;
  template class param_desc_impl<const signed char*>;
  template class param_desc_impl<const unsigned char*>;
  template class param_desc_impl<const void*>;
  template class param_desc_impl<float>;
  template class param_desc_impl<double>;
  template class param_desc_impl<bigint>;

  }

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
