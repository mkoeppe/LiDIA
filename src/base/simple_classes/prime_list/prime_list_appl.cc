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


#include	"LiDIA/bigint.h"
#include	"LiDIA/prime_list.h"

#include        <cstdlib>


#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



const char *tmp_file = "prime_list_appl.lst";

bool to_PLN(PRIME_LIST_NUMBER& target, bigint const& source) {
    bool result = false;

    if(!source.is_negative()) {
	long l;
	
	if(!source.longify(l)) {
	    target = static_cast<PRIME_LIST_NUMBER>(l);
	    result = true;
	}
	else {
	    bigint tmp(source);
	    tmp.divide_by_2();
	    if(!tmp.longify(l)) {
		target = 2 * static_cast<PRIME_LIST_NUMBER>(l) +
		    static_cast<PRIME_LIST_NUMBER>(source.is_odd());
		result = (source == target);
	    }
	    else {
		result = false;
	    }
	}
    }

    return result;
}

	    
int main_LiDIA(int argc, char** argv)
{
	PRIME_LIST_NUMBER lb = 2, ub = 0, m = 0, n = 0, fp = 0, lp = 0;
	bigint tmp;
	std::cout << "Please enter lower bound: " << std::flush;
	std::cin >> tmp;
	if(!to_PLN(lb, tmp)) {
	    std::cerr << "Error: input is negative or too big.\n";
	    return 1;
	}
	std::cout << "Please enter upper bound: " << std::flush;
	std::cin >> tmp;
	if(!to_PLN(ub, tmp)) {
	    std::cerr << "Error: input is negative or too big.\n";
	    return 1;
	}

	std::cout << "*****************************************************" << std::endl;
	std::cout << "**           Testprogram for prime_list            **" << std::endl;
	std::cout << "*****************************************************" << std::endl;

	std::cout << "check constructors ........................" << std::flush;
	prime_list p1(lb, ub, 'E');
	prime_list p2(lb, ub, 'K');
	prime_list p3(lb, ub, 'B');
	prime_list p4(lb, ub, '6');
	prime_list p5(lb, ub, 'I');
	std::cout << " completed" << std::endl;

	std::cout << "check get_number_of_primes ................" << std::flush;
	if (p1.get_number_of_primes() != p2.get_number_of_primes() ||
	    p1.get_number_of_primes() != p3.get_number_of_primes() ||
	    p1.get_number_of_primes() != p4.get_number_of_primes() ||
	    p1.get_number_of_primes() != p5.get_number_of_primes()) {
		std::cout << p1.get_number_of_primes() << std::endl;
		std::cout << p2.get_number_of_primes() << std::endl;
		std::cout << p3.get_number_of_primes() << std::endl;
		std::cout << p4.get_number_of_primes() << std::endl;
		std::cout << p5.get_number_of_primes() << std::endl;
		p1.get_first_prime();
		std::cout << p1.get_current_prime() << " ";
		while (p1.get_next_prime())
			std::cout << p1.get_current_prime() << " ";
		std::cout << std::endl;
		p2.get_first_prime();
		std::cout << p2.get_current_prime() << " ";
		while (p2.get_next_prime())
			std::cout << p2.get_current_prime() << " ";
		std::cout << std::endl;
		p3.get_first_prime();
		std::cout << p3.get_current_prime() << " ";
		while (p3.get_next_prime())
			std::cout << p3.get_current_prime() << " ";
		std::cout << std::endl;
		p4.get_first_prime();
		std::cout << p4.get_current_prime() << " ";
		while (p4.get_next_prime())
			std::cout << p4.get_current_prime() << " ";
		std::cout << std::endl;
		p5.get_first_prime();
		std::cout << p5.get_current_prime() << " ";
		while (p5.get_next_prime())
			std::cout << p5.get_current_prime() << " ";
		std::cout << std::endl;
		std::exit(1);
	}
	std::cout << " completed" << std::endl;

	std::cout << "check primes .............................." << std::flush;
	register lidia_size_t i;
	for (i = 0; i < p1.get_number_of_primes(); i++) {
		if (p1.get_prime(i) != p2.get_prime(i) ||
		    p1.get_prime(i) != p3.get_prime(i) ||
		    p1.get_prime(i) != p4.get_prime(i) ||
		    p1.get_prime(i) != p5.get_prime(i) ||
		    !is_prime(p1.get_prime(i)))
			std::exit(1);
	}
	std::cout << " completed" << std::endl;

	std::cout << "check get_first_prime ....................." << std::flush;
	if (p1.get_first_prime() != p2.get_first_prime() ||
	    p1.get_first_prime() != p3.get_first_prime() ||
	    p1.get_first_prime() != p4.get_first_prime() ||
	    p1.get_first_prime() != p5.get_first_prime()) {
		std::cout << p1.get_first_prime() << std::endl;
		std::cout << p2.get_first_prime() << std::endl;
		std::cout << p3.get_first_prime() << std::endl;
		std::cout << p4.get_first_prime() << std::endl;
		std::cout << p5.get_first_prime() << std::endl;
		std::exit(1);
	}
	int j;
	fp = p1.get_first_prime();
	lp = p1.get_last_prime();
	for (i = 0; i < p2.get_number_of_primes(); i++) {
		for (j = 0; j < 3; j++) {
			switch (j) {
			case 0:
				m = p2[i] - 1; break;
			case 1:
				m = p2[i]; break;
			case 2:
				m = p2[i] + 1; break;
			}
			n = p1.get_first_prime(m);
			if (m > lp) {
				if (n != 0) {
					std::cout << "get_first_prime(" << m << ") = " << n << std::endl;
					std::exit(1);
				}
			}
			else
				if ((n< m) || ((p1.get_prev_prime() > n) && (m > fp))) {
					std::cout << "get_first_prime(" << m << ") = " << n << std::endl;
					std::exit(1);
				}
		}
	}
	std::cout << " completed" << std::endl;

	std::cout << "check get_next_prime / get_prime .........." << std::flush;
	n = p1.get_first_prime();
	i = 0;
	while (n) {
		if (p2[i] != n)
			std::exit(1);
		n = p1.get_next_prime();
		i++;
		if (n && (p3[i] != n))
			std::exit(1);
		n = p1.get_next_prime();
		i++;
		if (n && (p4[i] != n))
			std::exit(1);
		n = p1.get_next_prime();
		i++;
		if (n && (p5[i] != n))
			std::exit(1);
		n = p1.get_next_prime();
		i++;
	}
	std::cout << " completed" << std::endl;

	std::cout << "check get_last_prime ......................" << std::flush;
	if (p1.get_last_prime() != p2.get_last_prime() ||
	    p1.get_last_prime() != p3.get_last_prime() ||
	    p1.get_last_prime() != p4.get_last_prime() ||
	    p1.get_last_prime() != p5.get_last_prime()) {
		std::cout << p1.get_last_prime() << std::endl;
		std::cout << p2.get_last_prime() << std::endl;
		std::cout << p3.get_last_prime() << std::endl;
		std::cout << p4.get_last_prime() << std::endl;
		std::cout << p5.get_last_prime() << std::endl;
		std::exit(1);
	}
	for (i = 0; i < p2.get_number_of_primes(); i++) {
		for (j = 0; j < 3; j++) {
			switch (j) {
			case 0:
				m = p2[i] - 1; break;
			case 1:
				m = p2[i]; break;
			case 2:
				m = p2[i] + 1; break;
			}
			n = p1.get_last_prime(m);
			if (m < fp) {
				if (n != 0) {
					std::cout << "get_last_prime(" << m << ") = " << n << std::endl;
					std::exit(1);
				}
			}
			else
				if ((n > m) || ((p1.get_next_prime() < n) && (m < lp))) {
					std::cout << "get_last_prime(" << m << ") = " << n << std::endl;
					std::exit(1);
				}
		}
	}
	std::cout << " completed" << std::endl;

	std::cout << "check get_prev_prime / get_prime .........." << std::flush;
	n = p1.get_last_prime();
	i = p1.get_number_of_primes() - 1;
	while (n) {
		if (p2[i] != n)
			std::exit(1);
		n = p1.get_prev_prime();
		i--;
		if (n && (p3[i] != n))
			std::exit(1);
		n = p1.get_prev_prime();
		i--;
		if (n && (p4[i] != n))
			std::exit(1);
		n = p1.get_prev_prime();
		i--;
		if (n && (p5[i] != n))
			std::exit(1);
		n = p1.get_prev_prime();
		i--;
	}
	std::cout << " completed" << std::endl;

	std::cout << "check get_index ..........................." << std::flush;
	unsigned long TMP1;
	lidia_size_t index;
	for (i = 0; i < p1.get_number_of_primes(); i++) {
		TMP1 = p1.get_prime(i);
		index = p1.get_index(TMP1);
		if (index != i)
			std::exit(1);
	}
	std::cout << " completed" << std::endl;

	std::cout << "check set_lower_bound ....................." << std::flush;
	if (p1.get_number_of_primes() > 0) {
		for (i = 10; i >= 2; i--) {
			m = fp + (lp - fp) / i;
			n = p2.get_first_prime(m);
			p1.set_lower_bound(m, "EKB6I"[i%5]);
			if (p1.get_first_prime() != n)
				std::exit(1);
			while (n) {
				n = p1.get_next_prime();
				if (n != p2.get_next_prime())
					std::exit(1);
			}
		}
		for (i = 2; i <= 11; i++) {
			m = fp + (lp - fp) / i;
			if (i == 11) m = lb;
			n = p2.get_first_prime(m);
			p1.set_lower_bound(m, "EKB6I"[i%5]);
			if (p1.get_first_prime() != n)
				std::exit(1);
			while (n) {
				n = p1.get_next_prime();
				if (n != p2.get_next_prime())
					std::exit(1);
			}
		}
	}
	std::cout << " completed" << std::endl;

	std::cout << "check set_upper_bound ....................." << std::flush;
	if (p1.get_number_of_primes() > 0) {
		for (i = 10; i >= 2; i--) {
			m = lp - (lp - fp) / i;
			n = p2.get_last_prime(m);
			p1.set_upper_bound(m, "EKB6I"[i%5]);
			if (p1.get_last_prime() != n)
				std::exit(1);
			while (n) {
				n = p1.get_prev_prime();
				if (n != p2.get_prev_prime())
					std::exit(1);
			}
		}
		for (i = 2; i <= 11; i++) {
			m = lp - (lp - fp) / i;
			if (i == 11) m = ub;
			n = p2.get_last_prime(m);
			p1.set_upper_bound(m, "EKB6I"[i%5]);
			if (p1.get_last_prime() != n)
				std::exit(1);
			while (n) {
				n = p1.get_prev_prime();
				if (n != p2.get_prev_prime())
					std::exit(1);
			}
		}
	}
	std::cout << " completed" << std::endl;

	std::cout << "check save_to_file / load_from_file ......." << std::flush;
	p1.save_to_file(tmp_file);
	if (p2.get_number_of_primes() > 100) {
		p1.load_from_file(tmp_file, 100);
		if (p1.get_number_of_primes() != 100)
			std::exit(1);
	}
	p1.load_from_file(tmp_file);
	n = p2.get_first_prime();
	if (p1.get_first_prime() != n)
		std::exit(1);
	while (n) {
		n = p1.get_next_prime();
		if (n != p2.get_next_prime())
			std::exit(1);
	}
	p1.save_to_file(tmp_file, true);
	p1.load_from_file(tmp_file);
	n = p2.get_first_prime();
	if (p1.get_first_prime() != n)
		std::exit(1);
	while (n) {
		n = p1.get_next_prime();
		if (n != p2.get_next_prime())
			std::exit(1);
	}
	std::cout << " completed" << std::endl;

	std::cout << "check assign .............................." << std::flush;
	p1.resize(0, 0);
	n = p2.get_prime(p2.get_number_of_primes() / 2);
	p1 = p2;
	if (p2.get_current_prime() != n)
		std::exit(1);
	n = p2.get_first_prime();
	if (p1.get_first_prime() != n)
		std::exit(1);
	while (n) {
		n = p1.get_next_prime();
		if (n != p2.get_next_prime())
			std::exit(1);
	}
	std::cout << " completed" << std::endl;

	std::cout.flush();

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
    catch(std::bad_alloc const& ex) {
	std::cerr << "\nout of memory: std::bad_alloc exception caught.\n";
	return 1;
    }
    catch(std::exception const& ex) {
	std::cerr << "\nunexpected exception: " << ex.what() << "\n";
	return 1;
    }
#endif
    
}
