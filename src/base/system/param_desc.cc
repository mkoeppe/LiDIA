// -*- C++ -*-
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


#include "LiDIA/param_desc.h"


#include "LiDIA/osstream.h"


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif

/**
 * free functions
 */

    std::ostream& operator<<(std::ostream& os,
							 const details::param_desc& pd) {
      return os << pd.nameStr() << " = " << pd.valueStr()
				<< " [Condition: " << pd.condStr() << "]";
    }
    
    void swap(details::param_desc*& ptr1,
			  details::param_desc*& ptr2) {
      details::param_desc* tmp = ptr1;
      ptr1 = ptr2;
      ptr2 = tmp;
    }
    
    std::istream& operator>>(std::istream& is,
							 details::param_desc*& ptr) {
      // dummy implementation
      // needed only in order to instantiate basic_error
      ptr = 0;
      return is;
    }
    
  namespace details {

    void swap(param_container& lhs, param_container& rhs) {
      swap(lhs.container_, rhs.container_);
    }
    
/**
 * class param_desc
 */
    
    param_desc::param_desc(const std::string& paramName,
				       const std::string& paramCond)
      : name_(paramName), cond_(paramName) {
      // nothing to do
    }
    
    param_desc::~param_desc() {
      // nothing to do
    }
    
    const std::string& param_desc::nameStr() const {
      return name_;
    }
    
    const std::string& param_desc::condStr() const {
      return cond_;
    }
    
/**
 * class param_container
 */
    
    param_container::param_container(size_type n) 
      : container_(n, 0, EXPAND) {
      // nothing to do
    }
    
    param_container::~param_container() {
      size_type noElems = container_.size();
      for(size_type index = 0; index < noElems; ++index) {
	delete container_[index];
      }
    }
    
    void param_container::push_back(param_container::elem_type pdp) {
      container_.insert_at(pdp, container_.size());
    }
    
    const param_container::elem_type param_container::operator[](size_type n) const {
      return (n < container_.size()) ? container_[n] : 0;
    }
    
    param_container::size_type param_container::size() const {
      return container_.size();
    }
    
  } // namespace details
  
#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif

//  Local Variables:    --- for emacs
//  mode: c++           --- for emacs
//  tab-width: 4        --- for emacs
//  End:                --- for emacs
