#include "LiDIA/bigint.h"
#include "LiDIA/gec_complex_multiplication.h"
#include "LiDIA/gec_point_counting_mod_p.h"
#include "LiDIA/gec_point_counting_mod_2n.h"

#include "LiDIA/timer.h"

#ifdef LIDIA_NAMESPACE
using namespace LiDIA;
#endif

int main_LiDIA(int, char**)
{
   lidia_size_t bitlength_r;
   std::cout << "Input lower bound of bitlength of r >> ";
   std::cin >> bitlength_r;

   bool verbosity = true;
   bool efficiency = true;

   bigint bound_k;
   std::cout << "Input upper_bound_k (k is the cofactor) >> ";
   std::cin >> bound_k;

   bool bsi = false;
   bigint field;
   std::cout << "Input field GF(q) (q ~ r * k) >> ";
   std::cin >> field;

   char method;
   std::cout << "Input method (c = CM), (p = PC mod p), (2 = PC mod 2^n) >> "; 
   std::cin >> method;

   if( method == 'c' )
   {
	   std::cout << "We use complex multiplication." << std::endl;

	   gec_complex_multiplication I;
	   I.set_lower_bound_bitlength_r( bitlength_r );
	   I.set_upper_bound_k( bound_k );
	   I.set_according_to_BSI( bsi );

	   bigint delta( -21311 );
	   std::cout << "Input delta >> "; std::cin >> delta;
	   if( delta != 0 )
		   I.set_delta(delta);
	   
	   if( field != 0 )
		   I.set_field( field );

	   I.set_verbose_level( verbosity );
//	   I.set_generation_mode( 1 );
//	   I.set_efficient_curve_parameters( efficiency );

	   timer all_time;

	   all_time.start_timer();
	   I.generate();
	   all_time.stop_timer();

//	   std::cout << I.get_class_polynomial() << std::endl;

	   std::cout << I.get_a4() << std::endl;
	   std::cout << I.get_a4().return_as_bigint() << std::endl;

	   std::cout << std::endl << "all_time = " << all_time << std::endl;

   }
   else if( method == 'p' )
   {
	   std::cout << "We use random approach." << std::endl;

	   gec_point_counting_mod_p I;
	   I.set_lower_bound_bitlength_r( bitlength_r );
	   I.set_upper_bound_k( bound_k );
	   I.set_according_to_BSI( bsi );
   
	   if( field != 0 )
		   I.set_field( field );
	   
	   I.set_verbose_level( verbosity );
//	   I.set_efficient_curve_parameters( efficiency );

	   I.generate();
   }
   else if( method == '2' )
   {
	   gec_point_counting_mod_2n I;
	   I.set_lower_bound_bitlength_r( bitlength_r );
	   I.set_upper_bound_k( bound_k );
	   I.set_according_to_BSI( bsi );

	   unsigned int degree;
	   std::cout << "Input degree >> "; std::cin >> degree;
   
	   if( degree != 0 )
		   I.set_degree( degree );
	   else if( field != 0 )
		   I.set_field( field );

	   I.set_verbose_level( verbosity );
	   I.set_efficient_curve_parameters( efficiency );

	   I.generate();
   }

   return( 0 );
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
