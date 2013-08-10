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
//	Author	: Thomas Papanikolaou (TP), Patrick Theobald (PT)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/LiDIA.h"
#include	<cstdlib>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



void default_error_handler(const char *f, const char *m)
{
	std::cerr << "\n error_handler";
	std::cerr << "::" << f;
	std::cerr << "::" << m;
	std::cerr << "\n";
	abort();
}



error_handler_ptr lidia_error_handler_internal = default_error_handler;



error_handler_ptr set_error_handler(error_handler_ptr new_handler)
{
	error_handler_ptr old_handler = lidia_error_handler_internal;
	lidia_error_handler_internal = new_handler;
	return old_handler;
}


/**
 * class basic_error
 */

basic_error::basic_error(const std::string& what_msg)
  : std::exception(), msg_(what_msg) {
  // Nothing to do beyond member initialization
}

basic_error::~basic_error() throw() {
  // Nothing to do
}

const char* basic_error::what() const throw() {
  return msg_.c_str();
}

void basic_error::traditional_error_handler() const {
  std::string classname;
  std::string msg;
  traditional_error_handler_impl(classname, msg);
  lidia_error_handler_internal(classname.c_str(), msg.c_str());

  // We should never get here!
  std::cerr << "\n\nOuch! lidia_error_handler() did return!\n"
	    << "Bailing out...\n\n";
  abort();
}

void basic_error::traditional_error_handler_impl(std::string& classname,
						std::string& msg) const {
  classname = offendingClass();
  msg = what();
}


/**
 * class generic_error
 */

generic_error::generic_error(const std::string& what_msg)
  : basic_error(what_msg), offendingClassName_("<generic_error>") {
  // Nothing to do beyond member initialization
};

generic_error::generic_error(const std::string& classname,
			   const std::string& what_msg)
  : basic_error(what_msg), offendingClassName_(classname) {
  // Nothing to do beyond member initialization
};

generic_error::~generic_error() throw() {
  // Nothing to do
}

const std::string& generic_error::offendingClass() const {
  return offendingClassName_;
}


/**
 * class index_out_of_bounds_error
 */

index_out_of_bounds_error::
index_out_of_bounds_error(const std::string& classname, long n) 
  : basic_error("index out of bounds"),
    std::out_of_range("LiDIA: index out of bounds"),
	 offendingClassName_(classname), index_(n) {
  // nothing to do
}

index_out_of_bounds_error::~index_out_of_bounds_error() throw() {
  // Nothing to do
}

const char* index_out_of_bounds_error::what() const throw() {
  return basic_error::what();
}

const std::string& index_out_of_bounds_error::offendingClass() const {
  return offendingClassName_;
}

long index_out_of_bounds_error::offendingIndex() const{
  return index_;
}


/**
 * cast_error
 */

cast_error::cast_error(const std::string& classname,
		     const std::string& what_msg)
  : basic_error(what_msg), std::bad_cast(),
	 offendingClassName_(classname) {
  // Nothing to do beyond member initialization
};

cast_error::~cast_error() throw() {
  // Nothing to do
}

const char* cast_error::what() const throw() {
  return basic_error::what();
}

const std::string& cast_error::offendingClass() const {
  return offendingClassName_;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif






