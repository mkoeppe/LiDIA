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


#include	"LiDIA/dense_power_series.h"
#include	"LiDIA/Fp_polynomial.h"

#include	<cassert>
#include        <cstdlib>



#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif



int main_LiDIA(int argc, char** argv)
{
	Fp_polynomial    		     pa, pb, pc;
	unsigned long		     da, db, dc;
	unsigned long		     power_of_2, expo;
	dense_power_series< udigit_mod > sa, sb, sc;
	udigit			     u_modulus;
	udigit			     maximal_number;
	lidia_size_t			     i;
	bigint 			     b_tmp;
	long				     u_tmp;
	unsigned long j;
	//   FILE *input = (FILE*) fopen("output_fft.txt","w");

	for (j = 1000; j <= 15000; j += 500) {

		std::cout << "degree: " << j << std::endl;

		da = db = j;
		dc = da + db + 1;

		//
		// power_of_2 = 2^expo >= dc
		//
		power_of_2 = 1;
		expo = 0;

		while (power_of_2 < dc) {
			power_of_2 <<= 1;
			expo++;
		}

		//
		// Find suitable characteristic of finite prime field
		//
		maximal_number = max_udigit_modulus () >> 1;

		u_modulus = maximal_number;
		u_modulus >>= expo;
		u_modulus <<= expo;

		// now: u_modulus = k * 2^n, k maximal

		if (u_modulus == maximal_number)
			if (u_modulus > power_of_2)
				u_modulus -= power_of_2;
			else
				std::cout << "overflow: power_of_2 equals maximal possible modulus" << std::endl;


		// now: u_modulus = k* 2^n, k maximal, with
		//      u_modulus+1 <= maximal_number;


		while (u_modulus != 0 &&
		       !is_prime(u_modulus+1)) {
			u_modulus -= power_of_2;
		}

		if (u_modulus != 0) {
			u_modulus++;
			std::cout << "found modulus " << u_modulus << std::endl;
		}
		else {
			std::cout << "No suitable fft-prime found." << std::endl;
			std::exit (1);
		}


		//
		// initialize
		//
		udigit_mod::set_modulus (u_modulus);

		pa.set_modulus (bigint(u_modulus));
		pb.set_modulus (bigint(u_modulus));

		pa.randomize (da);
		pb.randomize (db);

		for (i = da; i >= 0; i--) {
			b_tmp = pa[i];
			if (b_tmp.longify(u_tmp))
				std::cout << "error pa: " << b_tmp << " too big." << std::endl;
			else
				sa(i) = static_cast<udigit>(u_tmp);
		}

		for (i = db; i >= 0; i--) {
			b_tmp = pb[i];
			if (b_tmp.longify(u_tmp))
				std::cout << "error pb: " << b_tmp << " too big." << std::endl;
			else
				sb(i) = static_cast<udigit>(u_tmp);
		}

		//
		// compare the results of the multiplications
		//
		std::cout << "Multiplying polynomials ... " << std::endl;
		multiply     (pc, pa, pb);
		std::cout << "Multiplying power series ... " << std::endl;
		multiply_fft (sc, sa, sb);
		std::cout << "Done." << std::endl;

		for (i = sc.get_last()-1; i >= 0; i--)
			assert (bigint(sc[i].get_mantissa()) == pc[i]);

		//
		// print to file
		//
		//   fprintf (input,"size: %d\n",j);
		//std::cout << sc << "\n";
		//for (i = sc.get_last()-1; i >= 0; i--)
		//  fprintf (input,"%d ",sc[i].get_mantissa());

		// fprintf (input,"\n");
	};

	//fclose (input);
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
