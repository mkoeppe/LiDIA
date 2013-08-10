// -*- C++ -*-
//==============================================================================================
//      This file is part of LiDIA --- a library for computational number theory
//
//      Copyright (c) 1994--2002 the LiDIA Group.  All rights reserved.
//
//      See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//
//  File    : prime_proof_appl.cc
//  Author  : Jochen Hechler (JH) 
//      Changes : See CVS log
//
//==============================================================================================

// dies soll nur eine test klasse sein fuer die klasse ECPP.cc die noch prgrammiert wird
// erst mal ein paar compiler tests :-)


#include <iostream>       /* This is the stream definition file */
#include <stdio.h>
#include "LiDIA/bigint.h"
#include "LiDIA/prime_proof.h"
#include "LiDIA/certificate.h"

#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif

int main_LiDIA(int, char**)
{

    bigint p=1;
    int i=0;

    bigint exponent;
    bigint basis;
    std::cout << "-------------------------------------"<<"\n";
    std::cout << " Begin prime_proof test"<<  "\n";
    std::cout << "-------------------------------------"<<"\n";
    std::cout << " For a tests the following procedure can be used: bigint.next_prime(base^exponent) \n";
    std::cout << " Please enter base of prime (enter 0 for explicit prime input): "; std::cin >>basis;
    std::cout << "\n";
    if(!(basis==0)){
	std::cout << "Please enter exponent of prime (>1): "; 
	std::cin >>exponent;
	std::cout << "\n";
	while(i++ < exponent)p = p * basis;
	std::cout << "done["<<p<<"]"<<std::endl;
	std::cout << "is prime: "<<p.is_prime()<<"\n";
	p = p-1;
	p = p.next_prime();
    }
    if(basis==0){
	std::cout << "Please enter prime : "; 
	std::cin >>p;
	std::cout << "\n";
    }

    prime_proof pp;
    certificate c;
    bool suc = false;
    std::cout << "-------------------------------------------------------------------------------"<<"\n";
    std::cout << "| Using "<<p<<" as possible prime, length: "<<p.decimal_length()<<  "\n";
    std::cout << "-------------------------------------------------------------------------------"<<std::endl;
    if(p<17){
	std::cout << "The possible prime "<<p<<" is too small for prime_proof \n";
	std::cout << "Please use greater possible primes \n";
	return 1; 
    }
    pp.set_verbose(true);
    pp.set_ecpp_mode(2);
    pp.set_prime(p);
    suc = pp.prove_prime();
    if(suc)std::cout<<"SUCCESS"<<"\n"<<std::endl;
    c = pp.get_certificate();
    std::cout << "-------------------------------------"<<"\n";
    std::cout<<"Certificate retrieved"<<"\n";
    std::cout << "-------------------------------------"<<"\n"<<std::endl;
    c.write_certificate("cert.f");
    std::cout<<"Verify certificate"<<"\n";
    std::cout << "-------------------------------------"<<std::endl;
    if(pp.verify_certificate(c))std::cout<<"verify success"<<"\n";
    else std::cout<<"the certificate is not correct"<<"\n";
    std::cout << "-------------------------------------"<<"\n"<<"\n";


    std::cout << "-------------------------------------"<<"\n";
    std::cout << "End of test"<<  "\n";
    std::cout << "-------------------------------------"<<std::endl;

    return 0;
}

int main(int argc, char** argv) {

#if defined(LIDIA_EXCEPTIONS)
    try {
#endif

	main_LiDIA(argc, argv);
	
#if defined(LIDIA_EXCEPTIONS)
    }
    catch(basic_error const& ex) {
	ex.traditional_error_handler();
	return 1;
    }
    catch(std::exception const& ex) {
	std::cerr << "unexpected exception: " << ex.what() << "\n";
	return 1;
    }
#endif
    
}
