// -*- C++ -*-
//==============================================================================================
//
//      This file is part of LiDIA --- a library for computational number theory
//
//      Copyright (c) 1994--2003 the LiDIA Group.  All rights reserved.
//
//      See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//
//  File    : certificate.cc
//  Author  : Jochen Hechler (JH)
//      Changes : See CVS log
//
//==============================================================================================

#include        <iostream>
#include        <cstdio>
#include        <fstream>
#include        <cstdlib>

#include        "LiDIA/error.h"
#include        "LiDIA/bigint.h"
#include        "LiDIA/LiDIA.h"
#include        "LiDIA/base_vector.h"
#include	"LiDIA/certificate.h"

#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif


// This class holds an certificate for a primeproof.
// it is only a base vector of bigints
// the different tests are divided by bigint with a value of -2
// -2 is used, because no test can have a parameter of -2
// the tests are different in length 

    certificate::certificate()
    {
    }

    certificate::~certificate()
    {
    }

    void certificate::set_certificate(base_vector <bigint> vector)
    {
	r = vector;
	r.set_mode(EXPAND);
	
	
    }
 
    void certificate::write_certificate(std::string const& file)
    {
	std::ofstream out(file.c_str());
	if(!out.good()) {
	    std::string msg("can't open file ");
	    msg += file;
	    lidia_error_handler("certificate::write_certificate",
				msg.c_str());
	}
	
	out << std::endl;
	for(int i = 0;i<r.size();i++)
	{
	    out<<r[i]<< std::endl;
	}
	out.close();
	if(!out.good()) {
	    std::string msg("error while writing certificate to file ");
	    msg += file;
	    lidia_error_handler("certificate::write_certificate",
				msg.c_str());
	}
    }	

    void certificate::delete_all()
    {
	r.set_capacity(0);
    }
  
    void certificate::add(certificate cert)
    {
	base_vector <bigint> bv = cert.get_certification_vector();
	base_vector <bigint> space;
	space.set_capacity(1);
	space[0]=-2;
	r.concat(r,space);
	r.concat(r,bv);
    }
	
    base_vector <bigint> certificate::get_certification_vector()
    {
	return r;
    }
        
    // the following functions returns the certification vector of the 
    // $test_number$-th test. the first test is test_number 0;

    base_vector <bigint> certificate::get_cert_vector(int test_number)
    {
	int length = r.get_size();
	int tmp = 0;
	int i = 0;
	while(i < length && tmp<test_number )
	{
	    if(r[i++]==-2)tmp++;
	}
	base_vector <bigint> ret;
	int ii=0;
	ret.set_mode(EXPAND);
	ret.set_capacity(0);
	while(i < length && r[i]!=-2 )
	{
	    ret[ii]=r[i];	
	    ii++;
	    i++;
	}
	return ret;
    }

    int certificate::get_number_of_tests()
    {
	int ret=1;
	int length = r.get_size();
	int i = 0;
	while(i<length)
	{
	    if(r[i++]==-2)ret++;
	}
	return ret;
    }

    void certificate::read_certificate(std::string const& file)
    {
	r.set_capacity(0);
	r.set_mode(EXPAND);	

	std::ifstream in;
	in.open(file.c_str());
	if(!in.good()) {
	    std::string msg("can't open file ");
	    msg += file;
	    lidia_error_handler("certificate::read_certificate",
				msg.c_str());
	}

	bigint next;
	int i=0;
	do {
	    in >> next;
	    if(!(in.fail())) {
		r[i++]=next;
	    }
	    if(in.bad() || (in.fail() && !in.eof())) {
		// there was an input error not caused by reading over the EOF
		std::string msg("error while reading certificate from file ");
		msg += file;
		lidia_error_handler("certificate::read_certificate",
				    msg.c_str());
	    }
	} while(in.good());

	in.clear();
	in.close();
    }

#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
