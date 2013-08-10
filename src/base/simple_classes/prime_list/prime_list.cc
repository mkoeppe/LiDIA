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
//	Author	: Anja Jantke (AJ), Patrick Theobald (PT),
//                Dirk Schramm (DS)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/prime_list.h"
#include	"LiDIA/error.h"
#include	"LiDIA/xerror.h"
#include        "LiDIA/precondition_error.h"
#include	<fstream>
#include	<cmath>



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//
// DEBUG DEFINES / ERROR DEFINES
//

  static const char *DM_CP = "prime_list"; // Debug message / Error message

//
// HELP FUNCTIONS
//

  static PRIME_LIST_COUNTER
  invert(register PRIME_LIST_COUNTER a,
	 register PRIME_LIST_COUNTER b)
  {
    //
    // Invert a mod b
    //

    register PRIME_LIST_COUNTER r, q, x0 = 1, x1 = 0;
    register bool sign = false;
    PRIME_LIST_COUNTER sb = b;

    while (b != 0) {
      r = a % b;
      q = a / b;
      a = b;
      b = r;
      r = x1;
      x1 = q * x1 + x0;
      x0 = r;
      sign = !sign;
    }
    if (sign)
      x0 = sb - x0;
    return x0;
  }



  inline static void
  create_bit_sieve_prime_mask (register PRIME_LIST_COUNTER prime,
			       register PRIME_LIST_COUNTER base,
			       register PRIME_LIST_BIT_SIEVE *bit_mask,
			       register PRIME_LIST_BIT_SIEVE *prime_mask)
  {
    //
    // Create prime-mask with every prime-th bit set to 1,
    // starting with base-th bit
    //

    register PRIME_LIST_COUNTER index;

    for (index = 0; index < prime; index++)
      prime_mask[index] = 0;
    for (index = base; index < prime * prime_list_bit_sieve_size; index += prime)
      prime_mask[index / prime_list_bit_sieve_size] |=
	bit_mask[index % prime_list_bit_sieve_size];
  }



  inline static void
  apply_bit_sieve_prime_mask (register PRIME_LIST_BIT_SIEVE *prime_mask,
			      register PRIME_LIST_COUNTER prime,
			      register PRIME_LIST_BIT_SIEVE *sieve,
			      register PRIME_LIST_COUNTER base_index,
			      register PRIME_LIST_COUNTER sieve_size)
  {
    //
    // Apply prime-mask to bit-sieve
    //

    register PRIME_LIST_COUNTER index, index2;

    index2 = base_index % prime;
    for (index = 0; index < sieve_size; index++) {
      sieve[index] |= prime_mask[index2];
      index2++;
      if (index2 >= prime) index2 = 0;
    }
  }



  inline static unsigned int
  compute_diff (unsigned int value)
  {
    //
    // Convert stored diff-value to actual prime difference
    //

    return (value == 0) ? 1 : value << 1;
  }



//
// INTERNAL LIST MANAGEMENT
//

  void prime_list::init_prime_list ()
  {
    //
    // Initialize data members
    //

    number_of_blocks = 0;
    first_block = NULL;
    last_block = NULL;
    number_of_primes = 0;
    first_prime = 0;
    last_prime = 0;
    first_diff_index = 0;
    last_diff_index = prime_list_block_size - 1;
    current_index = 0;
    current_prime = 0;
    current_block = NULL;
    current_diff_index = 0;
    block_list = NULL;
  }



  void prime_list::create_prime_list (PRIME_LIST_NUMBER prime)
  {
    //
    // Create prime-list and set first prime
    //

    release_prime_list();
    first_block = new prime_block;
    memory_handler(first_block, DM_CP, "prime_list::create_prime_list: "
		   "Error in memory allocation (prime_block)");
    first_block->next_block = NULL;
    first_block->prev_block = NULL;
    last_block = first_block;
    number_of_blocks = 1;
    first_prime = prime;
    last_prime = prime;
    number_of_primes = 1;
    first_diff_index = 0;
    last_diff_index = prime_list_block_size - 1;
  }



  inline void
  prime_list::add_next_prime (PRIME_LIST_NUMBER prime)
  {
    //
    // Add primes in ascending order
    //

    if (last_diff_index < prime_list_block_size - 1) {
      last_diff_index++;
    }
    else {
      if (number_of_primes == 0) {
	create_prime_list(prime);
	return;
      }
      if (number_of_primes > 1) {
	last_block->next_block = new prime_block;
	memory_handler(last_block->next_block, DM_CP,
		       "prime_list::add_next_prime: "
		       "Error in memory allocation (prime_block)");
	last_block->next_block->prev_block = last_block;
	last_block->next_block->next_block = NULL;
	last_block = last_block->next_block;
	last_block->first_prime = last_prime;
	number_of_blocks++;
      }
      last_diff_index = 0;
    }
    number_of_primes++;
    last_block->diff[last_diff_index] =
      static_cast<PRIME_LIST_DIFF>((prime - last_prime) >> 1);
    last_prime = prime;
  }



  void
  prime_list::check_and_add_next_prime (PRIME_LIST_NUMBER prime)
  {
    //
    // Check validity of primes and add them in ascending order
    //

    if (number_of_primes > 0) {
      if (prime <= last_prime) {
	lidia_error_handler_para(
	  prime, "prime", "prime <= last_prime",
	  "prime_list::check_and_add_next_prime",
	  DM_CP, "Prime is not greater than last_prime!");
      }
      if ((prime - last_prime) > static_cast<PRIME_LIST_NUMBER>(prime_list_max_diff)) {
	lidia_error_handler_para(
	  prime, "prime", "prime - last_prime > prime_list_max_diff",
	  "prime_list::check_and_add_next_prime",
	  DM_CP, "Difference between prime and last_prime is too large!");
      }
      if (((prime - last_prime) & 1) != 0) // is difference odd?
	if ((prime - last_prime) != 1) {
	  lidia_error_handler_para(
	    prime, "prime", "(last_prime - prime) is odd and != 1",
	    "prime_list::check_and_add_next_prime",
	    DM_CP, "Difference between prime and last_prime is invalid!");
	}
    }
    add_next_prime(prime);
  }



  inline void
  prime_list::add_prev_prime (PRIME_LIST_NUMBER prime)
  {
    //
    // Add primes in descending order
    //

    if (first_diff_index > 0) {
      first_diff_index--;
    }
    else {
      if (number_of_primes == 0) {
	create_prime_list(prime);
	return;
      }
      if (number_of_primes > 1) {
	first_block->first_prime = first_prime;
	first_block->prev_block = new prime_block;
	memory_handler(first_block->prev_block, DM_CP,
		       "prime_list::add_prev_prime: "
		       "Error in memory allocation (prime_block)");
	first_block->prev_block->next_block = first_block;
	first_block->prev_block->prev_block = NULL;
	first_block = first_block->prev_block;
	number_of_blocks++;
      }
      first_diff_index = prime_list_block_size - 1;
    }
    number_of_primes++;
    first_block->diff[first_diff_index] =
      static_cast<PRIME_LIST_DIFF>((first_prime - prime) >> 1);
    first_prime = prime;
  }



  void
  prime_list::check_and_add_prev_prime (PRIME_LIST_NUMBER prime)
  {
    //
    // Check validity of primes and add them in descending order
    //

    if (number_of_primes > 0) {
      if (prime >= first_prime) {
	lidia_error_handler_para(
	  prime, "prime", "prime >= first_prime",
	  "prime_list::check_and_add_prev_prime",
	  DM_CP, "Prime is not smaller than first_prime!");
      }
      if ((first_prime - prime) > static_cast<PRIME_LIST_NUMBER>(prime_list_max_diff)) {
	lidia_error_handler_para(
	  prime, "prime", "first_prime - prime > prime_list_max_diff",
	  "prime_list::check_and_add_prev_prime",
	  DM_CP, "Difference between first_prime and prime is too large!");
      }
      if (((first_prime - prime) & 1) != 0) // is difference odd?
	if ((first_prime - prime) != 1) {
	  lidia_error_handler_para(
	    prime, "prime", "(prime - first_prime) is odd and != 1",
	    "prime_list::check_and_add_prev_prime",
	    DM_CP, "Difference between first_prime and prime is invalid!");
	}
    }
    add_prev_prime(prime);
  }



  void
  prime_list::create_block_list ()
  {
    //
    // Create array of prime-blocks
    //

    lidia_size_t i;
    prime_block *block;

    if (block_list) {
      delete[] block_list;
      block_list = NULL;
    }
    if (number_of_blocks > 0) {
      block_list = new prime_block*[number_of_blocks];
      memory_handler(block_list, DM_CP,
		     "prime_list::create_block_list: "
		     "Error in memory allocation (block_list)");
      block = first_block;
      for (i = 0; i < number_of_blocks; i++) {
	block_list[i] = block;
	block = block->next_block;
      }
    }
  }



  lidia_size_t
  prime_list::find_prime (PRIME_LIST_NUMBER prime,
			  bool set_current_prime) const
  {
    //
    // Find "prime" and return index or -1 if not in list
    // Optionally set current prime to first prime that is
    // equal to or greater than "prime"
    //

    register prime_block *block;
    lidia_size_t block_index, first_block_index, last_block_index;
    register PRIME_LIST_NUMBER tmp_prime;
    lidia_size_t index;
    register lidia_size_t diff_index;

    if ((prime< first_prime) || (prime > last_prime) || (number_of_primes == 0))
      return -1;
    if (number_of_primes == 1) return 0;

    first_block_index = 0;
    last_block_index = number_of_blocks - 1;
    while (first_block_index < last_block_index) {
      block_index = first_block_index + (last_block_index - first_block_index) / 2;
      block = block_list[block_index];
      if (prime < block_first_prime(block)) {
	last_block_index = block_index - 1;
      }
      else {
	if (prime >= block_last_prime(block)) {
	  first_block_index = block_index + 1;
	}
	else {
	  first_block_index = block_index;
	  last_block_index = block_index;
	}
      }
    }
    block_index = first_block_index;
    if (block_index > number_of_blocks - 1)
      block_index = number_of_blocks - 1;
    block = block_list[block_index];

    if ((prime - block_first_prime(block)) <=
	(block_last_prime(block) - block_first_prime(block)) / 2) {
      diff_index = block_first_diff_index(block);
      tmp_prime = block_first_prime(block);
      while (tmp_prime < prime) {
	tmp_prime += compute_diff(block->diff[diff_index]);
	diff_index++;
      }
    }
    else {
      diff_index = block_last_diff_index(block) + 1;
      tmp_prime = block_last_prime(block);
      while (tmp_prime > prime) {
	diff_index--;
	tmp_prime -= compute_diff(block->diff[diff_index]);
      }
      if (tmp_prime < prime) {
	tmp_prime += compute_diff(block->diff[diff_index]);
	diff_index++;
      }
    }
    index = block_index * prime_list_block_size - first_diff_index + diff_index;
    if (set_current_prime) {
      if (diff_index >= prime_list_block_size) {
	block = block->next_block;
	diff_index = 0;
      }
      current_index = index;
      current_prime = tmp_prime;
      current_block = block;
      current_diff_index = diff_index;
    }
    if (tmp_prime != prime)
      return -1;
    return index;
  }



  void
  prime_list::release_prime_list ()
  {
    //
    // Delete prime-list and free memory
    //

    prime_block *block, *next_block;

    if (block_list) {
      delete[] block_list;
      block_list = NULL;
    }
    block = first_block;
    while (block) {
      next_block = block->next_block;
      delete block;
      block = next_block;
    }
    init_prime_list();
  }



//
// CONSTRUCTORS
//

  prime_list::prime_list ()
  {
    //
    // Create an empty prime-list
    //

    init_prime_list();
    lower_bound = 0;
    upper_bound = 0;
  }



  prime_list::prime_list (PRIME_LIST_NUMBER init_upper_bound, char mode)
  {
    //
    // Create prime-list from 2 to upper_bound
    //

    init_prime_list();
    lower_bound = 2;
    upper_bound = init_upper_bound;
    sieve(lower_bound, upper_bound, false, mode);
  }



  prime_list::prime_list (PRIME_LIST_NUMBER init_lower_bound, PRIME_LIST_NUMBER init_upper_bound, char mode)
  {
    //
    // Create prime-list from lower_bound to upper_bound
    //

    init_prime_list();
    lower_bound = init_lower_bound;
    upper_bound = init_upper_bound;
    sieve(lower_bound, upper_bound, false, mode);
  }



  prime_list::prime_list (const char *filename, lidia_size_t max_number_of_primes)
  {
    //
    // Load primes from file
    //

    init_prime_list();
    lower_bound = 0;
    upper_bound = 0;
    load_from_file(filename, max_number_of_primes);
  }



  prime_list::prime_list (const prime_list& A)
  {
    //
    // Copy constructor
    //

    init_prime_list();
    lower_bound = 0;
    upper_bound = 0;
    assign(A);
  }



//
// DESTRUCTOR
//

  prime_list::~prime_list ()
  {
    release_prime_list();
  }



//
// GENERATION (sieve-algorithms)
//

  void
  prime_list::sieve (PRIME_LIST_NUMBER lower_bound, PRIME_LIST_NUMBER upper_bound, bool descending, char mode)
  {
    //
    // Create prime-list from lower_bound to upper_bound
    //

    if (lower_bound < 2)
      lower_bound = 2;
    if (upper_bound >= lower_bound) {
      switch (mode) {
	  case 'E':
	    sieve_e(lower_bound, upper_bound, descending);
	    break;
	  case 'K':
	    sieve_6k(lower_bound, upper_bound, descending);
	    break;
	  case 'B':
	    sieve_ebit(lower_bound, upper_bound, descending);
	    break;
	  case '6':
	    sieve_6kbit(lower_bound, upper_bound, descending);
	    break;
	  case 'I':
	    sieve_int(lower_bound, upper_bound, descending);
	    break;
	  default:
	    lidia_error_handler_para(mode, "mode", "mode in {E, K, B, 6, I}",
				     "prime_list::"
				     "sieve(lower_bound, upper_bound, descending, mode)",
				     DM_CP, "Error!! Mode not supported!");
      }
      create_block_list();
      get_first_prime();
    }
  }



  void
  prime_list::sieve_e (PRIME_LIST_NUMBER lower_bound, PRIME_LIST_NUMBER upper_bound, bool descending)
  {
    //
    // Calculate primes from lower_bound to upper_bound
    //
    // Algorithm: Sieve of Erathostenes
    //
    // Memory requirements (bytes):
    //   ((upper_bound - lower_bound) + min(lower_bound, sqrt(upper_bound))) * sizeof(PRIME_LIST_SIEVE)
    //

    register PRIME_LIST_SIEVE *base_sieve, *interval_sieve;
    bool sieve_overlap;
    register PRIME_LIST_COUNTER base_sieve_size, interval_sieve_size;
    register PRIME_LIST_COUNTER upper_bound_sqrt;
    register PRIME_LIST_COUNTER i, j;

    interval_sieve_size = upper_bound - lower_bound + 1;
    upper_bound_sqrt = static_cast<PRIME_LIST_COUNTER>(std::sqrt(static_cast<PRIME_LIST_FLOAT_NUMBER>(upper_bound))) + 1;
    interval_sieve = NULL;
    if (lower_bound <= static_cast<PRIME_LIST_NUMBER>(upper_bound_sqrt)) {
      base_sieve_size = upper_bound + 1;
      sieve_overlap = true;
    }
    else {
      base_sieve_size = upper_bound_sqrt + 1;
      sieve_overlap = false;
      interval_sieve = new PRIME_LIST_SIEVE[interval_sieve_size];
      memory_handler(interval_sieve, DM_CP, "prime_list::sieve_e: "
		     "Error in memory allocation (interval_sieve)");
    }
    base_sieve = new PRIME_LIST_SIEVE[base_sieve_size];
    memory_handler(base_sieve, DM_CP, "prime_list::sieve_e: "
		   "Error in memory allocation (base_sieve)");
    if (sieve_overlap)
      interval_sieve = &(base_sieve[lower_bound]);

    // Init
    for (i = 0; i < base_sieve_size; i++)
      base_sieve[i] = 0;
    if (!sieve_overlap)
      for (i = 0; i < interval_sieve_size; i++)
	interval_sieve[i] = 0;

    // base_sieve[i]=0 => i is prime
    // interval_sieve[i]=0 => lower_bound+i is prime

    for (i = 2; i <= upper_bound_sqrt; i++)
      if (base_sieve[i] == 0) {
	// i is prime
	for (j = 2 * i; j < base_sieve_size; j += i)
	  base_sieve[j] = 1;
	if (!sieve_overlap) {
	  j = lower_bound % i;
	  if (j != 0) j = i - j;
	  for (; j < interval_sieve_size; j += i)
	    interval_sieve[j] = 1;
	}
      }

    if (descending) {
      // save primes in descending order
      for (i = interval_sieve_size; (i--) > 0;)
	if (interval_sieve[i] == 0)
	  add_prev_prime(lower_bound + i);
    }
    else {
      // save primes in ascending order
      for (i = 0; i < interval_sieve_size; i++)
	if (interval_sieve[i] == 0)
	  add_next_prime(lower_bound + i);
    }

    delete[] base_sieve;
    if (!sieve_overlap) delete[] interval_sieve;
  }



  void
  prime_list::sieve_6k (PRIME_LIST_NUMBER lower_bound, PRIME_LIST_NUMBER upper_bound, bool descending)
  {
    //
    // Calculate primes from lower_bound to upper_bound
    //
    // Algorithm: 6k+-1 Sieve
    //
    // Memory requirements (bytes):
    //   ((upper_bound - lower_bound) + min(lower_bound, sqrt(upper_bound))) / 3 * sizeof(PRIME_LIST_SIEVE)
    //

    PRIME_LIST_NUMBER lower_bound_k, upper_bound_k;
    register PRIME_LIST_SIEVE *base_sieve1, *base_sieve2;
    register PRIME_LIST_SIEVE *interval_sieve1, *interval_sieve2;
    bool sieve_overlap;
    register PRIME_LIST_COUNTER base_sieve_size, interval_sieve_size;
    register PRIME_LIST_COUNTER upper_bound_sqrt_k;
    register PRIME_LIST_COUNTER i, k;
    PRIME_LIST_COUNTER j, tmp_prime;
    PRIME_LIST_NUMBER n;

    lower_bound_k = static_cast<PRIME_LIST_NUMBER>((lower_bound + 1) / 6);
    if (lower_bound_k < 1) lower_bound_k = 1;
    upper_bound_k = static_cast<PRIME_LIST_NUMBER>((upper_bound + 1) / 6);
    if (upper_bound_k < 1) upper_bound_k = 1;
    interval_sieve_size = upper_bound_k - lower_bound_k + 1;
    upper_bound_sqrt_k = static_cast<PRIME_LIST_COUNTER>((std::sqrt(static_cast<PRIME_LIST_FLOAT_NUMBER>(upper_bound)) + 1) / 6) + 1;
    interval_sieve1 = NULL;
    interval_sieve2 = NULL;
    if (lower_bound_k <= static_cast<PRIME_LIST_NUMBER>(upper_bound_sqrt_k)) {
      base_sieve_size = upper_bound_k + 1;
      sieve_overlap = true;
    }
    else {
      base_sieve_size = upper_bound_sqrt_k + 1;
      sieve_overlap = false;
      interval_sieve1 = new PRIME_LIST_SIEVE[interval_sieve_size];
      memory_handler(interval_sieve1, DM_CP, "prime_list::sieve_6k: "
		     "Error in memory allocation (interval_sieve1)");
      interval_sieve2 = new PRIME_LIST_SIEVE[interval_sieve_size];
      memory_handler(interval_sieve2, DM_CP, "prime_list::sieve_6k: "
		     "Error in memory allocation (interval_sieve2)");
    }
    // 0/1-Sieve-Array for 6k-1
    base_sieve1 = new PRIME_LIST_SIEVE[base_sieve_size];
    memory_handler(base_sieve1, DM_CP, "prime_list::sieve_6k: "
		   "Error in memory allocation (base_sieve1)");
    // 0/1-Sieve-Array for 6k+1
    base_sieve2 = new PRIME_LIST_SIEVE[base_sieve_size];
    memory_handler(base_sieve2, DM_CP, "prime_list::sieve_6k: "
		   "Error in memory allocation (base_sieve2)");
    if (sieve_overlap) {
      interval_sieve1 = &(base_sieve1[lower_bound_k]);
      interval_sieve2 = &(base_sieve2[lower_bound_k]);
    }

    // Init
    for (i = 0; i < base_sieve_size; i++) {
      base_sieve1[i] = 0;
      base_sieve2[i] = 0;
    }
    if (!sieve_overlap)
      for (i = 0; i < interval_sieve_size; i++) {
	interval_sieve1[i] = 0;
	interval_sieve2[i] = 0;
      }

    // base_sieve1[k]=0 => 6*k-1 is prime
    // base_sieve2[k]=0 => 6*k+1 is prime
    // interval_sieve1[k]=0 => 6*(lower_bound_k+k)-1 is prime
    // interval_sieve2[k]=0 => 6*(lower_bound_k+k)+1 is prime

    for (k = 1; k <= upper_bound_sqrt_k; k++) {
      if (base_sieve1[k] == 0) {
	tmp_prime = 6 * k - 1; // 6*k-1 is prime
	for (i = k + tmp_prime; i < base_sieve_size; i += tmp_prime)
	  base_sieve1[i] = 1;
	j = (tmp_prime != 5) ? invert(tmp_prime - 6, tmp_prime) : 4;
	for (i = j; i < base_sieve_size; i += tmp_prime)
	  base_sieve2[i] = 1;
	if (!sieve_overlap) {
	  i = (lower_bound_k - k) % tmp_prime;
	  if (i != 0) i = tmp_prime - i;
	  for (; i < interval_sieve_size; i += tmp_prime)
	    interval_sieve1[i] = 1;
	  if (static_cast<PRIME_LIST_NUMBER>(j) <= lower_bound_k) {
	    i = (lower_bound_k - j) % tmp_prime;
	    if (i != 0) i = tmp_prime - i;
	  }
	  else
	    i = j - lower_bound_k;
	  for (; i < interval_sieve_size; i += tmp_prime)
	    interval_sieve2[i] = 1;
	}
      }
      if (base_sieve2[k] == 0) {
	tmp_prime = 6 * k + 1; // 6*k+1 is prime
	for (i = k + tmp_prime; i < base_sieve_size; i += tmp_prime)
	  base_sieve2[i] = 1;
	j = invert(6, tmp_prime);
	for (i = j; i < base_sieve_size; i += tmp_prime)
	  base_sieve1[i] = 1;
	if (!sieve_overlap) {
	  i = (lower_bound_k - k) % tmp_prime;
	  if (i != 0) i = tmp_prime - i;
	  for (; i < interval_sieve_size; i += tmp_prime)
	    interval_sieve2[i] = 1;
	  if (static_cast<PRIME_LIST_NUMBER>(j) <= lower_bound_k) {
	    i = (lower_bound_k - j) % tmp_prime;
	    if (i != 0) i = tmp_prime - i;
	  }
	  else
	    i = j - lower_bound_k;
	  for (; i < interval_sieve_size; i += tmp_prime)
	    interval_sieve1[i] = 1;
	}
      }
    }

    if (descending) {
      // save primes in descending order
      k = interval_sieve_size - 1;
      n = (lower_bound_k + k) * 6;
      if (k > 0) {
	if ((interval_sieve2[k] == 0) && (n + 1 <= upper_bound))
	  add_prev_prime(n + 1);
	if ((interval_sieve1[k] == 0) && (n - 1 <= upper_bound))
	  add_prev_prime(n - 1);
	for (k--, n -= 6; k > 0; k--, n -= 6) {
	  if (interval_sieve2[k] == 0)
	    add_prev_prime(n + 1);
	  if (interval_sieve1[k] == 0)
	    add_prev_prime(n - 1);
	}
      }
      if ((interval_sieve2[k] == 0) && (n + 1 >= lower_bound) && (n + 1 <= upper_bound))
	add_prev_prime(n + 1);
      if ((interval_sieve1[k] == 0) && (n - 1 >= lower_bound) && (n - 1 <= upper_bound))
	add_prev_prime(n - 1);
      if ((3 >= lower_bound) && (3 <= upper_bound)) add_prev_prime(3);
      if ((2 >= lower_bound) && (2 <= upper_bound)) add_prev_prime(2);
    }
    else {
      // save primes in ascending order
      if ((2 >= lower_bound) && (2 <= upper_bound)) add_next_prime(2);
      if ((3 >= lower_bound) && (3 <= upper_bound)) add_next_prime(3);
      k = 0;
      n = lower_bound_k * 6;
      if ((interval_sieve1[k] == 0) && (n - 1 >= lower_bound) && (n - 1 <= upper_bound))
	add_next_prime(n - 1);
      if ((interval_sieve2[k] == 0) && (n + 1 >= lower_bound) && (n + 1 <= upper_bound))
	add_next_prime(n + 1);
      for (k++, n += 6; k < interval_sieve_size - 1; k++, n += 6) {
	if (interval_sieve1[k] == 0)
	  add_next_prime(n - 1);
	if (interval_sieve2[k] == 0)
	  add_next_prime(n + 1);
      }
      if (k < interval_sieve_size) {
	if ((interval_sieve1[k] == 0) && (n - 1 <= upper_bound))
	  add_next_prime(n - 1);
	if ((interval_sieve2[k] == 0) && (n + 1 <= upper_bound))
	  add_next_prime(n + 1);
      }
    }

    delete[] base_sieve1;
    delete[] base_sieve2;
    if (!sieve_overlap) {
      delete[] interval_sieve1;
      delete[] interval_sieve2;
    }
  }



  void
  prime_list::sieve_ebit (PRIME_LIST_NUMBER lower_bound, PRIME_LIST_NUMBER upper_bound, bool descending)
  {
    //
    // Calculate primes from lower_bound to upper_bound
    //
    // Algorithm: Sieve of Erathostenes (Bit-level)
    //
    // Memory requirements (bytes):
    //   ((upper_bound - lower_bound) + min(lower_bound, sqrt(upper_bound))) / 8
    //

    PRIME_LIST_NUMBER lower_bound_index, upper_bound_index;
    PRIME_LIST_NUMBER interval_sieve_lower_bound, interval_sieve_upper_bound;
    register PRIME_LIST_BIT_SIEVE *base_sieve, *interval_sieve;
    bool sieve_overlap;
    register PRIME_LIST_COUNTER base_sieve_size, interval_sieve_size;
    register PRIME_LIST_COUNTER base_sieve_bit_size, interval_sieve_bit_size;
    register PRIME_LIST_COUNTER upper_bound_sqrt_index;
    register PRIME_LIST_COUNTER index;
    register int bit;
    register PRIME_LIST_COUNTER i, j;
    register PRIME_LIST_NUMBER n;
    PRIME_LIST_BIT_SIEVE tmp_mask;
    PRIME_LIST_BIT_SIEVE bit_mask[prime_list_bit_sieve_size];
    PRIME_LIST_BIT_SIEVE prime_mask[prime_list_bit_sieve_size];

    lower_bound_index = static_cast<PRIME_LIST_NUMBER>(lower_bound / prime_list_bit_sieve_size);
    upper_bound_index = static_cast<PRIME_LIST_NUMBER>(upper_bound / prime_list_bit_sieve_size);
    interval_sieve_size = upper_bound_index - lower_bound_index + 1;
    interval_sieve_lower_bound = lower_bound_index * prime_list_bit_sieve_size;
    interval_sieve_upper_bound = (upper_bound_index + 1) * prime_list_bit_sieve_size - 1;
    upper_bound_sqrt_index = static_cast<PRIME_LIST_COUNTER>(std::sqrt(static_cast<PRIME_LIST_FLOAT_NUMBER>(upper_bound)) / prime_list_bit_sieve_size) + 1;
    interval_sieve = NULL;
    if (lower_bound_index <= static_cast<PRIME_LIST_NUMBER>(upper_bound_sqrt_index)) {
      base_sieve_size = upper_bound_index + 1;
      sieve_overlap = true;
    }
    else {
      base_sieve_size = upper_bound_sqrt_index + 1;
      sieve_overlap = false;
      interval_sieve = new PRIME_LIST_BIT_SIEVE[interval_sieve_size];
      memory_handler(interval_sieve, DM_CP, "prime_list::sieve_ebit: "
		     "Error in memory allocation (interval_sieve)");
    }
    base_sieve = new PRIME_LIST_BIT_SIEVE[base_sieve_size];
    memory_handler(base_sieve, DM_CP, "prime_list::sieve_ebit: "
		   "Error in memory allocation (base_sieve)");
    if (sieve_overlap)
      interval_sieve = &(base_sieve[lower_bound_index]);

    base_sieve_bit_size = base_sieve_size * prime_list_bit_sieve_size;
    interval_sieve_bit_size = interval_sieve_size * prime_list_bit_sieve_size;

    // Init
    for (index = 0; index < base_sieve_size; index++)
      base_sieve[index] = 0;
    if (!sieve_overlap)
      for (index = 0; index < interval_sieve_size; index++)
	interval_sieve[index] = 0;

    // bit_mask[i] = 2^i
    tmp_mask = 1;
    for (bit = 0; bit < prime_list_bit_sieve_size; bit++) {
      bit_mask[bit] = tmp_mask;
      tmp_mask = tmp_mask << 1;
    }

    // use prime-masks for small primes
    for (bit = 2; bit < prime_list_bit_sieve_size; bit++)
      if ((base_sieve[0] & bit_mask[bit]) == 0) {
	// bit is prime
	create_bit_sieve_prime_mask(bit, 0, bit_mask, prime_mask);
	apply_bit_sieve_prime_mask(prime_mask, bit, base_sieve, 0, base_sieve_size);
	base_sieve[0] &= ~bit_mask[bit];
	if (!sieve_overlap)
	  apply_bit_sieve_prime_mask(prime_mask, bit, interval_sieve,
				     lower_bound_index, interval_sieve_size);
      }

    for (index = 1; index <= upper_bound_sqrt_index; index++)
      for (bit = 0; bit < prime_list_bit_sieve_size; bit++)
	if ((base_sieve[index] & bit_mask[bit]) == 0) {
	  i = static_cast<PRIME_LIST_COUNTER>(index) * prime_list_bit_sieve_size + bit;
				// i is prime
	  n = static_cast<PRIME_LIST_NUMBER>(i) * i;
				// start sieving at i^2
	  if (n < static_cast<PRIME_LIST_NUMBER>(base_sieve_bit_size))
	    for (j = n; j < base_sieve_bit_size; j += i)
	      base_sieve[j / prime_list_bit_sieve_size] |=
		bit_mask[j % prime_list_bit_sieve_size];
	  if ((!sieve_overlap) && (n <= interval_sieve_upper_bound)) {
	    j = interval_sieve_lower_bound % i;
	    if (j != 0) j = i - j;
	    if (n > interval_sieve_lower_bound) {
	      n -= interval_sieve_lower_bound;
	      if (n > static_cast<PRIME_LIST_NUMBER>(j)) j = n;
	    }
	    for (; j < interval_sieve_bit_size; j += i)
	      interval_sieve[j / prime_list_bit_sieve_size] |=
		bit_mask[j % prime_list_bit_sieve_size];
	  }
	}

    if (descending) {
      // save primes in descending order
      index = interval_sieve_size - 1;
      n = interval_sieve_lower_bound + static_cast<PRIME_LIST_NUMBER>(index) * prime_list_bit_sieve_size;
      if (index > 0) {
	for (bit = prime_list_bit_sieve_size - 1; bit >= 0; bit--)
	  if (((interval_sieve[index] & bit_mask[bit]) == 0) && (n + bit <= upper_bound))
	    add_prev_prime(n + bit);
	for (index--, n -= prime_list_bit_sieve_size; index > 0;
	     index--, n -= prime_list_bit_sieve_size)
	  for (bit = prime_list_bit_sieve_size - 1; bit >= 0; bit--)
	    if ((interval_sieve[index] & bit_mask[bit]) == 0)
	      add_prev_prime(n + bit);
      }
      for (bit = prime_list_bit_sieve_size - 1; bit >= 0; bit--)
	if (((interval_sieve[index] & bit_mask[bit]) == 0) &&
	    (n + bit >= lower_bound) && (n + bit <= upper_bound))
	  add_prev_prime(n + bit);
    }
    else {
      // save primes in ascending order
      index = 0;
      n = interval_sieve_lower_bound;
      for (bit = 0; bit < prime_list_bit_sieve_size; bit++)
	if (((interval_sieve[index] & bit_mask[bit]) == 0) &&
	    (n + bit >= lower_bound) && (n + bit <= upper_bound))
	  add_next_prime(n + bit);
      for (index++, n += prime_list_bit_sieve_size;
	   index < interval_sieve_size - 1;
	   index++, n += prime_list_bit_sieve_size)
	for (bit = 0; bit < prime_list_bit_sieve_size; bit++)
	  if ((interval_sieve[index] & bit_mask[bit]) == 0)
	    add_next_prime(n + bit);
      if (index < interval_sieve_size)
	for (bit = 0; bit < prime_list_bit_sieve_size; bit++)
	  if (((interval_sieve[index] & bit_mask[bit]) == 0) && (n + bit <= upper_bound))
	    add_next_prime(n + bit);
    }

    delete[] base_sieve;
    if (!sieve_overlap) delete[] interval_sieve;
  }



  void
  prime_list::sieve_6kbit (PRIME_LIST_NUMBER lower_bound, PRIME_LIST_NUMBER upper_bound, bool descending)
  {
    //
    // Calculate primes from lower_bound to upper_bound
    //
    // Algorithm: 6k+-1 Bit sieve
    //
    // Memory requirements (bytes):
    //   ((upper_bound - lower_bound) + min(lower_bound, sqrt(upper_bound))) / 3 / 8
    //

    PRIME_LIST_NUMBER lower_bound_k, upper_bound_k;
    PRIME_LIST_NUMBER lower_bound_index, upper_bound_index;
    PRIME_LIST_NUMBER interval_sieve_lower_bound_k, interval_sieve_upper_bound_k;
    register PRIME_LIST_BIT_SIEVE *base_sieve1, *base_sieve2;
    register PRIME_LIST_BIT_SIEVE *interval_sieve1, *interval_sieve2;
    bool sieve_overlap;
    register PRIME_LIST_COUNTER base_sieve_size, interval_sieve_size;
    register PRIME_LIST_COUNTER base_sieve_bit_size, interval_sieve_bit_size;
    register PRIME_LIST_COUNTER upper_bound_sqrt_index;
    register PRIME_LIST_COUNTER index;
    register int bit;
    register PRIME_LIST_COUNTER i, k;
    PRIME_LIST_COUNTER j, tmp_prime;
    register PRIME_LIST_NUMBER n;
    PRIME_LIST_BIT_SIEVE tmp_mask;
    PRIME_LIST_BIT_SIEVE bit_mask[prime_list_bit_sieve_size];
    PRIME_LIST_BIT_SIEVE prime_mask[prime_list_bit_sieve_size];

    lower_bound_k = static_cast<PRIME_LIST_NUMBER>((lower_bound + 1) / 6);
    if (lower_bound_k < 1) lower_bound_k = 1;
    upper_bound_k = static_cast<PRIME_LIST_NUMBER>((upper_bound + 1) / 6);
    if (upper_bound_k < 1) upper_bound_k = 1;
    lower_bound_index = static_cast<PRIME_LIST_NUMBER>(lower_bound_k / prime_list_bit_sieve_size);
    upper_bound_index = static_cast<PRIME_LIST_NUMBER>(upper_bound_k / prime_list_bit_sieve_size);
    interval_sieve_size = upper_bound_index - lower_bound_index + 1;
    interval_sieve_lower_bound_k = lower_bound_index * prime_list_bit_sieve_size;
    interval_sieve_upper_bound_k = (upper_bound_index + 1) * prime_list_bit_sieve_size - 1;
    upper_bound_sqrt_index = static_cast<PRIME_LIST_COUNTER>((std::sqrt(static_cast<PRIME_LIST_FLOAT_NUMBER>(upper_bound)) + 1) / 6 / prime_list_bit_sieve_size) + 1;
    interval_sieve1 = NULL;
    interval_sieve2 = NULL;
    if (lower_bound_index <= static_cast<PRIME_LIST_NUMBER>(upper_bound_sqrt_index)) {
      base_sieve_size = upper_bound_index + 1;
      sieve_overlap = true;
    }
    else {
      base_sieve_size = upper_bound_sqrt_index + 1;
      sieve_overlap = false;
      interval_sieve1 = new PRIME_LIST_BIT_SIEVE[interval_sieve_size];
      memory_handler(interval_sieve1, DM_CP, "prime_list::sieve_6kbit: "
		     "Error in memory allocation (interval_sieve1)");
      interval_sieve2 = new PRIME_LIST_BIT_SIEVE[interval_sieve_size];
      memory_handler(interval_sieve2, DM_CP, "prime_list::sieve_6kbit: "
		     "Error in memory allocation (interval_sieve2)");
    }
    base_sieve1 = new PRIME_LIST_BIT_SIEVE[base_sieve_size];
    memory_handler(base_sieve1, DM_CP, "prime_list::sieve_6kbit: "
		   "Error in memory allocation (base_sieve1)");
    base_sieve2 = new PRIME_LIST_BIT_SIEVE[base_sieve_size];
    memory_handler(base_sieve2, DM_CP, "prime_list::sieve_6kbit: "
		   "Error in memory allocation (base_sieve2)");
    if (sieve_overlap) {
      interval_sieve1 = &(base_sieve1[lower_bound_index]);
      interval_sieve2 = &(base_sieve2[lower_bound_index]);
    }

    base_sieve_bit_size = base_sieve_size * prime_list_bit_sieve_size;
    interval_sieve_bit_size = interval_sieve_size * prime_list_bit_sieve_size;

    // Init
    for (index = 0; index < base_sieve_size; index++) {
      base_sieve1[index] = 0;
      base_sieve2[index] = 0;
    }
    if (!sieve_overlap)
      for (index = 0; index < interval_sieve_size; index++) {
	interval_sieve1[index] = 0;
	interval_sieve2[index] = 0;
      }

    // bit_mask[i] = 2^i
    tmp_mask = 1;
    for (bit = 0; bit < prime_list_bit_sieve_size; bit++) {
      bit_mask[bit] = tmp_mask;
      tmp_mask = tmp_mask << 1;
    }

    // use prime-masks for small primes
    for (bit = 1; 6 * bit + 1 < prime_list_bit_sieve_size; bit++) {
      if ((base_sieve1[0] & bit_mask[bit]) == 0) {
	tmp_prime = 6 * bit - 1; // 6*bit-1 is prime
	create_bit_sieve_prime_mask(tmp_prime, bit, bit_mask, prime_mask);
	apply_bit_sieve_prime_mask(prime_mask, tmp_prime, base_sieve1, 0, base_sieve_size);
	base_sieve1[0] &= ~bit_mask[bit];
	if (!sieve_overlap)
	  apply_bit_sieve_prime_mask(prime_mask, tmp_prime, interval_sieve1,
				     lower_bound_index, interval_sieve_size);
	j = (tmp_prime != 5) ? invert(tmp_prime - 6, tmp_prime) : 4;
	create_bit_sieve_prime_mask(tmp_prime, j, bit_mask, prime_mask);
	apply_bit_sieve_prime_mask(prime_mask, tmp_prime, base_sieve2, 0, base_sieve_size);
	if (!sieve_overlap)
	  apply_bit_sieve_prime_mask(prime_mask, tmp_prime, interval_sieve2,
				     lower_bound_index, interval_sieve_size);
      }
      if ((base_sieve2[0] & bit_mask[bit]) == 0) {
	tmp_prime = 6 * bit + 1; // 6*bit+1 is prime
	create_bit_sieve_prime_mask(tmp_prime, bit, bit_mask, prime_mask);
	apply_bit_sieve_prime_mask(prime_mask, tmp_prime, base_sieve2, 0, base_sieve_size);
	base_sieve2[0] &= ~bit_mask[bit];
	if (!sieve_overlap)
	  apply_bit_sieve_prime_mask(prime_mask, tmp_prime, interval_sieve2,
				     lower_bound_index, interval_sieve_size);
	j = invert(6, tmp_prime);
	create_bit_sieve_prime_mask(tmp_prime, j, bit_mask, prime_mask);
	apply_bit_sieve_prime_mask(prime_mask, tmp_prime, base_sieve1, 0, base_sieve_size);
	if (!sieve_overlap)
	  apply_bit_sieve_prime_mask(prime_mask, tmp_prime, interval_sieve1,
				     lower_bound_index, interval_sieve_size);
      }
    }

    for (index = 0; index <= upper_bound_sqrt_index; index++) {
      for (; bit < prime_list_bit_sieve_size; bit++) {
	if ((base_sieve1[index] & bit_mask[bit]) == 0) {
	  k = static_cast<PRIME_LIST_COUNTER>(index) * prime_list_bit_sieve_size + bit;
	  tmp_prime = k * 6 - 1; // 6*k-1 is prime
	  for (i = k + tmp_prime; i < base_sieve_bit_size; i += tmp_prime)
	    base_sieve1[i / prime_list_bit_sieve_size] |=
	      bit_mask[i % prime_list_bit_sieve_size];
	  j = (tmp_prime != 5) ? invert(tmp_prime - 6, tmp_prime) : 4;
	  for (i = j; i < base_sieve_bit_size; i += tmp_prime)
	    base_sieve2[i / prime_list_bit_sieve_size] |=
	      bit_mask[i % prime_list_bit_sieve_size];
	  if (!sieve_overlap) {
	    i = (interval_sieve_lower_bound_k - k) % tmp_prime;
	    if (i != 0) i = tmp_prime - i;
	    for (; i < interval_sieve_bit_size; i += tmp_prime)
	      interval_sieve1[i / prime_list_bit_sieve_size] |=
		bit_mask[i % prime_list_bit_sieve_size];
	    if (static_cast<PRIME_LIST_NUMBER>(j) <= interval_sieve_lower_bound_k) {
	      i = (interval_sieve_lower_bound_k - j) % tmp_prime;
	      if (i != 0) i = tmp_prime - i;
	    }
	    else
	      i = j - interval_sieve_lower_bound_k;
	    for (; i < interval_sieve_bit_size; i += tmp_prime)
	      interval_sieve2[i / prime_list_bit_sieve_size] |=
		bit_mask[i % prime_list_bit_sieve_size];
	  }
	}
	if ((base_sieve2[index] & bit_mask[bit]) == 0) {
	  k = static_cast<PRIME_LIST_COUNTER>(index) * prime_list_bit_sieve_size + bit;
	  tmp_prime = k * 6 + 1; // 6*k+1 is prime
	  for (i = k + tmp_prime; i < base_sieve_bit_size; i += tmp_prime)
	    base_sieve2[i / prime_list_bit_sieve_size] |=
	      bit_mask[i % prime_list_bit_sieve_size];
	  j = invert(6, tmp_prime);
	  for (i = j; i < base_sieve_bit_size; i += tmp_prime)
	    base_sieve1[i / prime_list_bit_sieve_size] |=
	      bit_mask[i % prime_list_bit_sieve_size];
	  if (!sieve_overlap) {
	    i = (interval_sieve_lower_bound_k - k) % tmp_prime;
	    if (i != 0) i = tmp_prime - i;
	    for (; i < interval_sieve_bit_size; i += tmp_prime)
	      interval_sieve2[i / prime_list_bit_sieve_size] |=
		bit_mask[i % prime_list_bit_sieve_size];
	    if (static_cast<PRIME_LIST_NUMBER>(j) <= interval_sieve_lower_bound_k) {
	      i = (interval_sieve_lower_bound_k - j) % tmp_prime;
	      if (i != 0) i = tmp_prime - i;
	    }
	    else
	      i = j - interval_sieve_lower_bound_k;
	    for (; i < interval_sieve_bit_size; i += tmp_prime)
	      interval_sieve1[i / prime_list_bit_sieve_size] |=
		bit_mask[i % prime_list_bit_sieve_size];
	  }
	}
      }
      bit = 0;
    }

    if (descending) {
      // save primes in descending order
      index = interval_sieve_size - 1;
      k = static_cast<PRIME_LIST_COUNTER>(index + 1) * prime_list_bit_sieve_size - 1;
      n = (interval_sieve_lower_bound_k + k) * 6;
      if (index > 0) {
	for (bit = prime_list_bit_sieve_size - 1; bit >= 0; bit--, n -= 6) {
	  if (((interval_sieve2[index] & bit_mask[bit]) == 0) && (n + 1 <= upper_bound))
	    add_prev_prime(n + 1);
	  if (((interval_sieve1[index] & bit_mask[bit]) == 0) && (n - 1 <= upper_bound))
	    add_prev_prime(n - 1);
	}
	for (index--; index > 0; index--)
	  for (bit = prime_list_bit_sieve_size - 1; bit >= 0; bit--, n -= 6) {
	    if ((interval_sieve2[index] & bit_mask[bit]) == 0)
	      add_prev_prime(n + 1);
	    if ((interval_sieve1[index] & bit_mask[bit]) == 0)
	      add_prev_prime(n - 1);
	  }
      }
      for (bit = prime_list_bit_sieve_size - 1; bit >= 0; bit--, n -= 6) {
	if (((interval_sieve2[index] & bit_mask[bit]) == 0) &&
	    (n + 1 >= lower_bound) && (n + 1 <= upper_bound))
	  add_prev_prime(n + 1);
	if (((interval_sieve1[index] & bit_mask[bit]) == 0) &&
	    (n - 1 >= lower_bound) && (n - 1 <= upper_bound))
	  add_prev_prime(n - 1);
      }
      if ((3 >= lower_bound) && (3 <= upper_bound)) add_prev_prime(3);
      if ((2 >= lower_bound) && (2 <= upper_bound)) add_prev_prime(2);
    }
    else {
      // save primes in ascending order
      if ((2 >= lower_bound) && (2 <= upper_bound)) add_next_prime(2);
      if ((3 >= lower_bound) && (3 <= upper_bound)) add_next_prime(3);
      index = 0;
      k = 0;
      n = interval_sieve_lower_bound_k * 6;
      for (bit = 0; bit < prime_list_bit_sieve_size; bit++, n += 6) {
	if (((interval_sieve1[index] & bit_mask[bit]) == 0) &&
	    (n - 1 >= lower_bound) && (n - 1 <= upper_bound))
	  add_next_prime(n - 1);
	if (((interval_sieve2[index] & bit_mask[bit]) == 0) &&
	    (n + 1 >= lower_bound) && (n + 1 <= upper_bound))
	  add_next_prime(n + 1);
      }
      for (index++; index < interval_sieve_size - 1; index++)
	for (bit = 0; bit < prime_list_bit_sieve_size; bit++, n += 6) {
	  if ((interval_sieve1[index] & bit_mask[bit]) == 0)
	    add_next_prime(n - 1);
	  if ((interval_sieve2[index] & bit_mask[bit]) == 0)
	    add_next_prime(n + 1);
	}
      if (index < interval_sieve_size)
	for (bit = 0; bit < prime_list_bit_sieve_size; bit++, n += 6) {
	  if (((interval_sieve1[index] & bit_mask[bit]) == 0) && (n - 1 <= upper_bound))
	    add_next_prime(n - 1);
	  if (((interval_sieve2[index] & bit_mask[bit]) == 0) && (n + 1 <= upper_bound))
	    add_next_prime(n + 1);
	}
    }

    delete[] base_sieve1;
    delete[] base_sieve2;
    if (!sieve_overlap) {
      delete[] interval_sieve1;
      delete[] interval_sieve2;
    }
  }



  void
  prime_list::sieve_int (PRIME_LIST_NUMBER lower_bound, PRIME_LIST_NUMBER upper_bound, bool descending)
  {
    //
    // Calculate primes from lower_bound to upper_bound
    //
    // Algorithm: interval sieve
    //
    // Approximate memory requirements (bytes):
    //   min(max_sieve_size, max((upper_bound - lower_bound), sqrt(upper_bound))) * sizeof(PRIME_LIST_SIEVE)
    //   + sqrt(upper_bound) / log(sqrt(upper_bound)) * sizeof(PRIME_LIST_COUNTER)
    //

    const PRIME_LIST_COUNTER max_sieve_size = 1000000;

    PRIME_LIST_COUNTER upper_bound_sqrt;
    register PRIME_LIST_SIEVE *sieve;
    PRIME_LIST_COUNTER sieve_size;
    register PRIME_LIST_COUNTER current_sieve_size;
    register PRIME_LIST_NUMBER sieve_lower_bound;
    PRIME_LIST_NUMBER lower_bound2;
    int interval;
    register PRIME_LIST_COUNTER sieve_upper_bound_sqrt;
    register PRIME_LIST_COUNTER i, j, k;
    register PRIME_LIST_COUNTER *prime_factor_list;
    PRIME_LIST_COUNTER prime_factor_list_size;
    PRIME_LIST_COUNTER prime_factor_count;

    upper_bound_sqrt = static_cast<PRIME_LIST_COUNTER>(std::sqrt(static_cast<PRIME_LIST_FLOAT_NUMBER>(upper_bound))) + 1;
    sieve_size = upper_bound_sqrt + 1;
    if (upper_bound - lower_bound + 1 > static_cast<PRIME_LIST_NUMBER>(sieve_size))
      sieve_size = upper_bound - lower_bound + 1;
    if (sieve_size > max_sieve_size)
      sieve_size = max_sieve_size;

    sieve = new PRIME_LIST_SIEVE[sieve_size];
    memory_handler(base_sieve, DM_CP, "prime_list::sieve_int: "
		   "Error in memory allocation (sieve)");

    prime_factor_list_size = (upper_bound_sqrt < 20) ? 10 :
      static_cast<PRIME_LIST_COUNTER>(upper_bound_sqrt /
				      (std::log(double(upper_bound_sqrt))
				       - 1.5));
    prime_factor_list = new PRIME_LIST_COUNTER[prime_factor_list_size];
    memory_handler(base_sieve, DM_CP, "prime_list::sieve_int: "
		   "Error in memory allocation (prime_factor_list)");
    prime_factor_count = 0;

    for (sieve_lower_bound = 0; sieve_lower_bound < static_cast<PRIME_LIST_NUMBER>(upper_bound_sqrt); sieve_lower_bound += sieve_size) {
      current_sieve_size = sieve_size;
      if (sieve_lower_bound + sieve_size - 1 < static_cast<PRIME_LIST_NUMBER>(upper_bound_sqrt))
	current_sieve_size = sieve_size;
      else
	current_sieve_size = upper_bound_sqrt - sieve_lower_bound + 1;
      sieve_upper_bound_sqrt = static_cast<PRIME_LIST_COUNTER>
	(std::sqrt(static_cast<PRIME_LIST_FLOAT_NUMBER>(sieve_lower_bound) + current_sieve_size - 1)) + 1;
      for (i = 0; i < current_sieve_size; i++)
	sieve[i] = 0;
      if (sieve_lower_bound == 0) {
	// first interval
	sieve[0] = 1; // 0 and 1 aren't prime
	sieve[1] = 1;
	for (i = 2; i <= sieve_upper_bound_sqrt; i++) {
	  if (sieve[i] == 0)
	    // i is prime
	    for (j = 2 * i; j < current_sieve_size; j += i)
	      sieve[j] = 1;
	}
      }
      else {
	// use calculated prime factors to sieve interval
	i = 0;
	do
	{
	  k = prime_factor_list[i];
	  j = sieve_lower_bound % k;
	  if (j != 0) j = k - j;
	  for (; j < current_sieve_size; j += k)
	    sieve[j] = 1;
	  i++;
	} while ((k < sieve_upper_bound_sqrt) && (i < prime_factor_count));
      }
      // save prime factors
      for (i = 0; i < current_sieve_size; i++)
	if (sieve[i] == 0) {
	  prime_factor_list[prime_factor_count] = sieve_lower_bound + i;
	  prime_factor_count++;
	}
    }
    if (descending) {
      // sieve remaining intervals and save primes in descending order
      lower_bound2 = lower_bound;
      if (static_cast<PRIME_LIST_NUMBER>(upper_bound_sqrt) >= lower_bound)
	lower_bound2 = upper_bound_sqrt + 1;
      if (lower_bound2 <= upper_bound) {
	interval = (upper_bound - lower_bound2) / sieve_size;
	for (; interval >= 0; interval--) {
	  sieve_lower_bound = lower_bound2 + sieve_size * interval;
	  current_sieve_size = sieve_size;
	  if (sieve_lower_bound + sieve_size - 1 >= upper_bound)
	    current_sieve_size = upper_bound - sieve_lower_bound + 1;
	  sieve_upper_bound_sqrt =
	    static_cast<PRIME_LIST_COUNTER>(
	      std::sqrt(
		static_cast<PRIME_LIST_FLOAT_NUMBER>(sieve_lower_bound) +
		current_sieve_size - 1)) + 1;
	  for (i = 0; i < current_sieve_size; i++)
	    sieve[i] = 0;
	  i = 0;
	  do
	  {
	    k = prime_factor_list[i];
	    j = sieve_lower_bound % k;
	    if (j != 0) j = k - j;
	    for (; j < current_sieve_size; j += k)
	      sieve[j] = 1;
	    i++;
	  } while ((k < sieve_upper_bound_sqrt) && (i < prime_factor_count));
	  for (i = current_sieve_size; (i--) > 0;)
	    if (sieve[i] == 0)
	      add_prev_prime(sieve_lower_bound + i);
	}
      }
      for (i = prime_factor_count; (i--) > 0;)
	if (static_cast<PRIME_LIST_NUMBER>(prime_factor_list[i]) >= lower_bound)
	  add_prev_prime(prime_factor_list[i]);
    }
    else {
      // sieve remaining intervals and save primes in ascending order
      for (i = 0; i < prime_factor_count; i++)
	if (static_cast<PRIME_LIST_NUMBER>(prime_factor_list[i]) >= lower_bound)
	  add_next_prime(prime_factor_list[i]);
      sieve_lower_bound = lower_bound;
      if (static_cast<PRIME_LIST_NUMBER>(upper_bound_sqrt) >= lower_bound)
	sieve_lower_bound = upper_bound_sqrt + 1;
      for (; sieve_lower_bound <= upper_bound; sieve_lower_bound += sieve_size) {
	current_sieve_size = sieve_size;
	if (sieve_lower_bound + sieve_size - 1 >= upper_bound)
	  current_sieve_size = upper_bound - sieve_lower_bound + 1;
	sieve_upper_bound_sqrt = static_cast<PRIME_LIST_COUNTER>
	  (std::sqrt(static_cast<PRIME_LIST_FLOAT_NUMBER>(sieve_lower_bound) + current_sieve_size - 1)) + 1;
	for (i = 0; i < current_sieve_size; i++)
	  sieve[i] = 0;
	i = 0;
	do
	{
	  k = prime_factor_list[i];
	  j = sieve_lower_bound % k;
	  if (j != 0) j = k - j;
	  for (; j < current_sieve_size; j += k)
	    sieve[j] = 1;
	  i++;
	} while ((k < sieve_upper_bound_sqrt) && (i < prime_factor_count));
	for (i = 0; i < current_sieve_size; i++)
	  if (sieve[i] == 0)
	    add_next_prime(sieve_lower_bound + i);
      }
    }
    delete[] sieve;
    delete[] prime_factor_list;
  }



//
// BASIC FUNCTIONS
//

  void
  prime_list::set_lower_bound (PRIME_LIST_NUMBER new_lower_bound, char mode)
  {
    //
    // Set new lower bound for prime list,
    // delete primes / calculate new primes (with algorithm "mode")
    // as necessary
    //

    prime_block *block, *next_block;

    if (new_lower_bound > first_prime) {
      if (new_lower_bound > last_prime) {
	release_prime_list();
      }
      else {
	if (get_first_prime(new_lower_bound) == last_prime) {
	  create_prime_list(last_prime);
	}
	else {
	  block = first_block;
	  while (block != current_block) {
	    next_block = block->next_block;
	    delete block;
	    number_of_blocks--;
	    block = next_block;
	  }
	  current_block->prev_block = NULL;
	  first_diff_index = current_diff_index;
	  first_prime = current_prime;
	  first_block = current_block;
	  number_of_primes -= current_index;
	  create_block_list();
	}
      }
    }
    else {
      if (new_lower_bound < lower_bound) {
	sieve(new_lower_bound, lower_bound - 1, true, mode);
      }
    }
    get_first_prime();
    lower_bound = new_lower_bound;
  }



  void
  prime_list::set_upper_bound (PRIME_LIST_NUMBER new_upper_bound, char mode)
  {
    //
    // Set new upper bound for prime list,
    // delete primes / calculate new primes (with algorithm "mode")
    // as necessary
    //

    prime_block *block, *prev_block;
    PRIME_LIST_NUMBER new_last_prime;

    if (new_upper_bound < last_prime) {
      if (new_upper_bound < first_prime) {
	release_prime_list();
      }
      else {
	new_last_prime = get_last_prime(new_upper_bound);
	if (new_last_prime == first_prime) {
	  create_prime_list(first_prime);
	}
	else {
	  get_prev_prime();
	  block = last_block;
	  while (block != current_block) {
	    prev_block = block->prev_block;
	    delete block;
	    number_of_blocks--;
	    block = prev_block;
	  }
	  current_block->next_block = NULL;
	  last_diff_index = current_diff_index;
	  last_prime = new_last_prime;
	  last_block = current_block;
	  number_of_primes = current_index + 2;
	  create_block_list();
	}
      }
    }
    else {
      if (new_upper_bound > upper_bound) {
	sieve(upper_bound + 1, new_upper_bound, false, mode);
      }
    }
    get_first_prime();
    upper_bound = new_upper_bound;
  }



  void
  prime_list::resize (PRIME_LIST_NUMBER new_lower_bound, PRIME_LIST_NUMBER new_upper_bound, char mode)
  {
    //
    // Set new lower bound and upper bound for prime list,
    // delete primes / calculate new primes (with algorithm "mode")
    // as necessary
    //
    if ((new_lower_bound > last_prime) || (new_upper_bound < first_prime)) {
      release_prime_list();
      sieve(new_lower_bound, new_upper_bound, false, mode);
    }
    else {
      set_lower_bound(new_lower_bound, mode);
      set_upper_bound(new_upper_bound, mode);
    }
    lower_bound = new_lower_bound;
    upper_bound = new_upper_bound;
  }



  bool
  prime_list::is_element (PRIME_LIST_NUMBER prime) const
  {
    //
    // Check if "prime" is in list
    //

    return find_prime(prime, false) >= 0;
  }



  lidia_size_t
  prime_list::get_index (PRIME_LIST_NUMBER prime) const
  {
    //
    // Return index of "prime" in list
    // Return 0 if "prime" is not in list
    //

    return find_prime(prime, false);
  }



//
// ACCESS FUNCTIONS
//

  PRIME_LIST_NUMBER
  prime_list::get_first_prime () const
  {
    //
    // Return first prime in list and update current prime
    // Return 0 if list is empty
    //

    if (number_of_primes == 0) return 0;
    current_index = 0;
    current_prime = first_prime;
    current_block = first_block;
    current_diff_index = first_diff_index;
    return current_prime;
  }



  PRIME_LIST_NUMBER
  prime_list::get_first_prime (PRIME_LIST_NUMBER prime) const
  {
    //
    // Return first prime in list equal to or greater than "prime"
    // and update current prime
    // Return 0 if "prime" > last prime in list
    //

    if (prime <= first_prime) return get_first_prime();
    if (prime > last_prime) return 0;
    find_prime(prime, true);
    return current_prime;
  }



  PRIME_LIST_NUMBER
  prime_list::get_next_prime () const
  {
    //
    // Return prime that follows current prime, update current prime
    // Return 0 if current prime = last prime
    //

    if (current_index >= number_of_primes - 1) return 0;
    current_index++;
    current_prime += compute_diff(current_block->diff[current_diff_index]);
    current_diff_index++;
    if (current_diff_index >= prime_list_block_size) {
      if (current_block->next_block != NULL) {
	current_block = current_block->next_block;
	current_diff_index = 0;
      }
    }
    return current_prime;
  }



  PRIME_LIST_NUMBER
  prime_list::get_last_prime () const
  {
    //
    // Return last prime in list and update current prime
    // Return 0 if list is empty
    //

    if (number_of_primes == 0) return 0;
    current_index = number_of_primes - 1;
    current_prime = last_prime;
    current_block = last_block;
    current_diff_index = last_diff_index + 1;
    return current_prime;
  }



  PRIME_LIST_NUMBER
  prime_list::get_last_prime (PRIME_LIST_NUMBER prime) const
  {
    //
    // Return last prime in list smaller than or equal to "prime"
    // and update current prime
    // Return 0 if "prime" < first prime
    //

    if (prime < first_prime) return 0;
    if (prime >= last_prime) return get_last_prime();
    if (find_prime(prime, true) >= 0)
      return current_prime;
    else
      return get_prev_prime();
  }



  PRIME_LIST_NUMBER
  prime_list::get_prev_prime () const
  {
    //
    // Return prime before current prime, update current prime
    // Return 0 if current prime is first prime
    //

    if (current_index <= 0) return 0;
    current_index--;
    if (current_diff_index <= 0) {
      current_block = current_block->prev_block;
      current_diff_index = prime_list_block_size;
    }
    current_diff_index--;
    current_prime -= compute_diff(current_block->diff[current_diff_index]);
    return current_prime;
  }



  PRIME_LIST_NUMBER
  prime_list::get_prime (lidia_size_t index) const
  {
    //
    // Return prime at position "index" (starting with 0)
    // Return 0 if index< 0 or index >= number of primes
    //

    lidia_size_t block_index;
    register prime_block *block;
    register lidia_size_t new_diff_index, diff_index;
    register PRIME_LIST_NUMBER prime;

    if ((index< 0) || (index >= number_of_primes)) return 0;
    if (index == current_index) return current_prime;
    if (index == current_index + 1) return get_next_prime();
    if (index == current_index - 1) return get_prev_prime();

    block_index = static_cast<lidia_size_t>((index + first_diff_index) / prime_list_block_size);
    current_block = block_list[block_index];
    current_index = (block_index == 0) ? 0 :
      block_index * prime_list_block_size - first_diff_index;
    block = current_block;
    new_diff_index = index - current_index + block_first_diff_index(block);
    if ((index - current_index) <= (block_size(block) / 2)) {
      diff_index = block_first_diff_index(block);
      prime = block_first_prime(block);
      while (diff_index < new_diff_index) {
	prime += compute_diff(block->diff[diff_index]);
	diff_index++;
      }
    }
    else {
      diff_index = block_last_diff_index(block) + 1;
      prime = block_last_prime(block);
      while (diff_index > new_diff_index) {
	diff_index--;
	prime -= compute_diff(block->diff[diff_index]);
      }

    }
    current_index = index;
    current_diff_index = new_diff_index;
    current_prime = prime;
    return current_prime;
  }



//
// INPUT/OUTPUT
//

  void
  prime_list::load_from_file (const char *filename, lidia_size_t max_number_of_primes)
  {
    // enforce filename != 0
    if(filename == 0) {
      precondition_error_handler(static_cast<void const*>(filename),
				 "filename", "filename != 0",
				 max_number_of_primes, "max_number_of_primes",
				 "max_number_of_primes >= 0",
				 "void prime_list::"
				 "load_from_file(char const* filename, "
				 "lidia_size_t max_number_of_primes)",
				 "class prime_list",
				 "filename must point to a C-string!");
    }

    //
    // Load primes from file
    // Number of primes being read can be limited by setting
    // "max_number_of_primes" > 0
    //

    PRIME_LIST_NUMBER prime, first_prime;
    std::ifstream in(filename, std::ios::in);
    if (in.fail()) { // Recall: fail() checks both failbit and badbit
	std::string error_msg;
	error_msg += "Can't open file \"";
	error_msg += filename;
	error_msg += "\"!";
	lidia_error_handler("prime_list", error_msg.c_str());
    }

    release_prime_list();
    if (!in.eof()) {
      in >> first_prime;
      create_prime_list(first_prime);
      max_number_of_primes--;
      if (in.good() && (max_number_of_primes != 0)) {
	in >> prime;
	if (in.good()) {
	  max_number_of_primes--;
	  if (prime > first_prime) {
	    check_and_add_next_prime(prime);
	    while (in.good() && (max_number_of_primes != 0)) {
	      in >> prime;
	      if (in.good()) {
		check_and_add_next_prime(prime);
		max_number_of_primes--;
	      }
	    }
	  }
	  else {
	    check_and_add_prev_prime(prime);
	    while (in.good() && (max_number_of_primes != 0)) {
	      in >> prime;
	      if (in.good()) {
		check_and_add_prev_prime(prime);
		max_number_of_primes--;
	      }
	    }
	  }
	}
      }
    }
    if (in.bad() || (in.fail() && !in.eof())) {
	std::string error_msg;
	error_msg += "An error occured while reading from file \"";
	error_msg += filename;
	error_msg += "\"!";
	lidia_error_handler("prime_list", error_msg.c_str());
    }
    in.close();
    create_block_list();
    get_first_prime();
    lower_bound = first_prime;
    upper_bound = last_prime;
  }



  void
  prime_list::save_to_file (const char *filename, bool descending) const
  {
    //
    // Save primes to file in ascending/descending order
    //

    lidia_size_t old_index;
    PRIME_LIST_NUMBER prime;
    std::ofstream out(filename);

    if (out.fail())
      lidia_error_handler_para(filename, "filename", "",
			       "prime_list::save_to_file",
			       DM_CP, "Can't open file!");
    old_index = current_index; // save current state
    if (descending) {
      prime = get_last_prime();
      while (prime && out.good()) {
	out << prime << "\n";
	prime = get_prev_prime();
      }
    }
    else {
      prime = get_first_prime();
      while (prime && out.good()) {
	out << prime << "\n";
	prime = get_next_prime();
      }
    }
    if (out.fail())
      lidia_error_handler_para(filename, "filename", "",
			       "prime_list::save_to_file",
			       DM_CP, "An error occured, while writing to file!");
    out.close();
    get_prime(old_index); // restore initial state
  }



//
// ASSIGNMENT
//

  void
  prime_list::assign (const prime_list& A)
  {
    //
    // Copy primes from A
    //

    PRIME_LIST_NUMBER prime;
    lidia_size_t old_index;

    release_prime_list();
    old_index = A.get_current_index(); // save A's current state
    prime = A.get_first_prime();
    while (prime != 0) {
      add_next_prime(prime);
      prime = A.get_next_prime();
    }
    lower_bound = A.get_lower_bound();
    upper_bound = A.get_upper_bound();
    A.get_prime(old_index); // restore A's initial state
  }



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
