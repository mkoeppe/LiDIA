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


#include "LiDIA/bigint.h"
#include "LiDIA/precondition_error.cc"

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif

  template void precondition_error::addParam(const std::string& paramName,
					    const bool& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const short& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const unsigned short& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const int& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const unsigned int& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const long& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const unsigned long& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const char& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const signed char& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const unsigned char& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const char * const& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const signed char * const& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const unsigned char * const& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const void * const& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const float& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const double& value,
					    const std::string& paramCond);
  template void precondition_error::addParam(const std::string& paramName,
					    const bigint& value,
					    const std::string& paramCond);

  template const bool& precondition_error::
  paramValue<bool>(precondition_error::size_type n) const;
  template const short& precondition_error::
  paramValue<short>(precondition_error::size_type n) const;
  template const unsigned short& precondition_error::
  paramValue<unsigned short>(precondition_error::size_type n) const;
  template const int& precondition_error::
  paramValue<int>(precondition_error::size_type n) const;
  template const unsigned int& precondition_error::
  paramValue<unsigned int>(precondition_error::size_type n) const;
  template const long& precondition_error::
  paramValue<long>(precondition_error::size_type n) const;
  template const unsigned long& precondition_error::
  paramValue<unsigned long>(precondition_error::size_type n) const;
  template const char& precondition_error::
  paramValue<char>(precondition_error::size_type n) const;
  template const signed char& precondition_error::
  paramValue<signed char>(precondition_error::size_type n) const;
  template const unsigned char& precondition_error::
  paramValue<unsigned char>(precondition_error::size_type n) const;
  template const char* const& precondition_error::
  paramValue<const char*>(precondition_error::size_type n) const;
  template const signed char* const& precondition_error::
  paramValue<const signed char*>(precondition_error::size_type n) const;
  template const unsigned char* const& precondition_error::
  paramValue<const unsigned char*>(precondition_error::size_type n) const;
  template const void* const& precondition_error::
  paramValue<const void*>(precondition_error::size_type n) const;
  template const float& precondition_error::
  paramValue<float>(precondition_error::size_type n) const;
  template const double& precondition_error::
  paramValue<double>(precondition_error::size_type n) const;
  template const bigint& precondition_error::
  paramValue<bigint>(precondition_error::size_type n) const;


  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const unsigned long & para1, const char *name_1,
			  const char *cond1,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const char & para1, const char *name_1,
			  const char *cond1,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const char * const & para1, const char *name_1,
			  const char *cond1,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const bigint & para1, const char *name_1,
			  const char *cond1,
			  const char *proto, const char *file,
			  const char *msg);
  
  template
  void precondition_error_handler(const int & para1, const char *name_1,
				  const char *cond1,
				  const int & para2, const char *name_2,
				  const char *cond2,
				  const char *proto, const char *file,
				  const char *msg);
  template
  void precondition_error_handler(void const* const& para1, const char *name_1,
				  const char *cond1,
				  const int & para2, const char *name_2,
				  const char *cond2,
				  const char *proto, const char *file,
				  const char *msg);
  template
  void precondition_error_handler(const unsigned int& para1,
				  const char *name_1,
				  const char *cond1,
				  const unsigned int & para2,
				  const char *name_2,
				  const char *cond2,
				  const char *proto, const char *file,
				  const char *msg);
  template 
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const unsigned long & para2, const char *name_2,
			  const char *cond2,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const unsigned long & para1, const char *name_1,
			  const char *cond1,
			  const unsigned long & para2, const char *name_2,
			  const char *cond2,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const char & para2, const char *name_2,
			  const char *cond2,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const char * const & para2, const char *name_2,
			  const char *cond2,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const bigint & para1, const char *name_1,
			  const char *cond1,
			  const bigint & para2, const char *name_2,
			  const char *cond2,
			  const char *proto, const char *file,
			  const char *msg);

  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const int & para3, const char *name_3,
			  const char *cond3,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const bigint & para3, const char *name_3,
			  const char *cond3,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const char * const & para3, const char *name_3,
			  const char *cond3,
			  const char *proto, const char *file,
			  const char *msg);

  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const int & para3, const char *name_3,
			  const char *cond3,
			  const int & para4, const char *name_4,
			  const char *cond4,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const bigint & para1, const char *name_1,
			  const char *cond1,
			  const bigint & para2, const char *name_2,
			  const char *cond2,
			  const bigint & para3, const char *name_3,
			  const char *cond3,
			  const bigint & para4, const char *name_4,
			  const char *cond4,
			  const char *proto, const char *file,
			  const char *msg);

  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const int & para3, const char *name_3,
			  const char *cond3,
			  const int & para4, const char *name_4,
			  const char *cond4,
			  const int & para5, const char *name_5,
			  const char *cond5,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const int & para3, const char *name_3,
			  const char *cond3,
			  const int & para4, const char *name_4,
			  const char *cond4,
			  const char *const & para5, const char *name_5,
			  const char *cond5,
			  const char *proto, const char *file,
			  const char *msg);

  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const int & para3, const char *name_3,
			  const char *cond3,
			  const int & para4, const char *name_4,
			  const char *cond4,
			  const int & para5, const char *name_5,
			  const char *cond5,
			  const int & para6, const char *name_6,
			  const char *cond6,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const int & para3, const char *name_3,
			  const char *cond3,
			  const int & para4, const char *name_4,
			  const char *cond4,
			  const int & para5, const char *name_5,
			  const char *cond5,
			  const char *const & para6, const char *name_6,
			  const char *cond6,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const int & para3, const char *name_3,
			  const char *cond3,
			  const int & para4, const char *name_4,
			  const char *cond4,
			  const long & para5, const char *name_5,
			  const char *cond5,
			  const bigint & para6, const char *name_6,
			  const char *cond6,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const int & para3, const char *name_3,
			  const char *cond3,
			  const int & para4, const char *name_4,
			  const char *cond4,
			  const bigint & para5, const char *name_5,
			  const char *cond5,
			  const bigint & para6, const char *name_6,
			  const char *cond6,
			  const char *proto, const char *file,
			  const char *msg);

  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const int & para3, const char *name_3,
			  const char *cond3,
			  const int & para4, const char *name_4,
			  const char *cond4,
			  const int & para5, const char *name_5,
			  const char *cond5,
			  const int & para6, const char *name_6,
			  const char *cond6,
			  const bigint & para7, const char *name_7,
			  const char *cond7,
			  const bigint & para8, const char *name_8,
			  const char *cond8,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const int & para3, const char *name_3,
			  const char *cond3,
			  const int & para4, const char *name_4,
			  const char *cond4,
			  const int & para5, const char *name_5,
			  const char *cond5,
			  const char* const & para6, const char *name_6,
			  const char *cond6,
			  const bigint & para7, const char *name_7,
			  const char *cond7,
			  const bigint & para8, const char *name_8,
			  const char *cond8,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const int & para3, const char *name_3,
			  const char *cond3,
			  const int & para4, const char *name_4,
			  const char *cond4,
			  const int & para5, const char *name_5,
			  const char *cond5,
			  const char* const & para6, const char *name_6,
			  const char *cond6,
			  const long & para7, const char *name_7,
			  const char *cond7,
			  const bigint & para8, const char *name_8,
			  const char *cond8,
			  const char *proto, const char *file,
			  const char *msg);
  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const int & para3, const char *name_3,
			  const char *cond3,
			  const int & para4, const char *name_4,
			  const char *cond4,
			  const int & para5, const char *name_5,
			  const char *cond5,
			  const int & para6, const char *name_6,
			  const char *cond6,
			  const long & para7, const char *name_7,
			  const char *cond7,
			  const bigint & para8, const char *name_8,
			  const char *cond8,
			  const char *proto, const char *file,
			  const char *msg);

  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const int & para3, const char *name_3,
			  const char *cond3,
			  const int & para4, const char *name_4,
			  const char *cond4,
			  const int & para5, const char *name_5,
			  const char *cond5,
			  const int & para6, const char *name_6,
			  const char *cond6,
			  const int & para7, const char *name_7,
			  const char *cond7,
			  const int & para8, const char *name_8,
			  const char *cond8,
			  const int & para9, const char *name_9,
			  const char *cond9,
			  const int & para10, const char *name_10,
			  const char *cond10,
			  const char *proto, const char *file,
			  const char *msg);

  template
  void precondition_error_handler(const int & para1, const char *name_1,
			  const char *cond1,
			  const int & para2, const char *name_2,
			  const char *cond2,
			  const int & para3, const char *name_3,
			  const char *cond3,
			  const int & para4, const char *name_4,
			  const char *cond4,
			  const int & para5, const char *name_5,
			  const char *cond5,
			  const int & para6, const char *name_6,
			  const char *cond6,
			  const int & para7, const char *name_7,
			  const char *cond7,
			  const int & para8, const char *name_8,
			  const char *cond8,
			  const int & para9, const char *name_9,
			  const char *cond9,
			  const int & para10, const char *name_10,
			  const char *cond10,
			  const bigint & para11, const char *name_11,
			  const char *cond11,
			  const bigint & para12, const char *name_12,
			  const char *cond12,
			  const char *proto, const char *file,
			  const char *msg);

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
