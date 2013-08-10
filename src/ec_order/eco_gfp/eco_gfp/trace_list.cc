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
//	Author	: Volker Mueller (VM), Markus Maurer (MM), Denis Pochuev
//	Changes	: See CVS log
//
//==============================================================================================


// A trace_list is a list of trace_mod's



#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/trace_list.h"
#include	"LiDIA/crt_table.h"
#include	"LiDIA/crt.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



int trace_list::info; // info mode
int trace_list::MAX_NOF_TRACES = 9000000;


void trace_list::set_info_mode(int i)
{
	trace_list::info = i;
}



void trace_list::set_max_nof_traces(int m)
{
	if (m > 0)
		trace_list::MAX_NOF_TRACES = m;
}



trace_list::trace_list()
{
	last_index = -1;
	M1.assign_one();
	M2.assign_one();
	M3.assign_one();
	C3.assign_zero();
}



trace_list::~trace_list()
{
}



std::ostream & operator << (std::ostream & out, const trace_list & t)
{
	out << t.E << " " << t.C3 << " " << t.M3 << " " << t.l;
	return out;
}



std::istream & operator >> (std::istream & in, trace_list & t)
{
	elliptic_curve< gf_element > e;
	t.clear();

	in >> e;
	t.set_curve(e);

	in >> t.C3;
	in >> t.M3;
	in >> t.l;
	return in;
}



void trace_list::set_curve(const elliptic_curve< gf_element > & e)
{
	clear();
	E = e;

        power (four_sqrt_q, E.get_a6 ().characteristic (),
	       E.degree_of_definition ());
	sqrt (four_sqrt_q, four_sqrt_q);

	inc(four_sqrt_q);
	shift_left(four_sqrt_q, four_sqrt_q, 2);
}



const bigint & trace_list::get_M1() const
{
	return (M1);
}



const bigint & trace_list::get_M2() const
{
	return (M2);
}



const bigint & trace_list::get_M3() const
{
	return (M3);
}



const bigint& trace_list::get_C3() const
{
	return (C3);
}



bigint trace_list::get_absolute_smallest_C3() const
{
	// return c = c3 mod M3 with |c| <= M3/2.
	bigint h;

	shift_right(h, M3, 1);

	if (C3 <= h)
		return C3;
	else
		return (C3 - M3);
}



const sort_vector< trace_mod >& trace_list::get_list() const
{
	return (l);
}



void trace_list::tl_sort (lidia_size_t left, lidia_size_t right)
{
	l.sort(compare, left, right);
}



void trace_list::clear()
{
	last_index = -1;
	l.kill();
	M1.assign_one();
	M2.assign_one();
	M3.assign_one();
	C3.assign_zero();
}



//----------------------------------------------------------------
// Return the number of combinations which are necessary until
// M = prod(i=0, last_index) m_i > 4*sqrt(q). If M > 4 sqrt(q),
// set last_index.
//
// If M3 is larger than 4 * sqrt(q) then return 0. If the product
// of all moduli is smaller than 4 sqrt(q), then return -1.
//

bigint trace_list::number_of_combinations()
{
	bigint num, M;
	lidia_size_t i, list_len;

	// Sort the relations. Increasing density.
	//
	tl_sort(0, l.get_size()-1);

	// Check for M3 >= 4 sqrt(q).
	//
	if (M3 >= four_sqrt_q) {
		last_index = -1;
		return 0;
	}

	// Multiply moduli until product is larger than 4 sqrt(p)
	// or end of list is reached.
	//
	M.assign(M3);
	num.assign_one();
	list_len = l.get_size();

	for (i = 0; i < list_len; i++) {
		multiply(M, M, l[i].get_modulus());
		multiply(num, num, l[i].get_size_of_trace_mod());
		if (M >= four_sqrt_q)
			break;
	}

	// If product is too small, return -1.
	//
	if (M < four_sqrt_q)
		return -1;

	// If product is large enough, return the number
	// of combinations for the moduli of the product.
	// Also set last_index to the index of last modul used
	// for the product.
	//
	else {
		last_index = i;
		return num;
	}
}



//----------------------------------------------------------------
// append a new trace_mod to trace_list. An Elkies prime is added to
// (C3, M3), an Atkin prime to the list of candidates.
// Then the number of candidates is checked, iff BG is possible, return
// true.

bool trace_list::append(const trace_mod & tm)
{
	if (tm.get_size_of_trace_mod() == 1) // Elkies prime
	{
		C3 = chinese_remainder(C3, M3, bigint(tm.get_first_element()),
				       bigint(tm.get_modulus()));
		multiply(M3, M3, tm.get_modulus());
	}
	else {
		l.set_mode(EXPAND);
		l[l.get_size()] = tm;
		l.set_mode(FIXED);
	}
	bigint nc = number_of_combinations();

	if (nc < trace_list::MAX_NOF_TRACES && !nc.is_negative())
		return true;
	else
		return false;
}



//------------------------------------------------------------
// split the trace_list into baby- and giantsteps, compute
// these two lists of values and corresponding moduli.
// kills the list of trace_mod's to safe memory.
//

bool trace_list::split_baby_giant(sort_vector< bigint > & baby,
		  	          sort_vector< bigint > & giant)
{
	bigint nof_comb, nof_baby;
	bigint h, up, M;
	udigit msmall;
	lidia_size_t i, s, k, e, j;
	lidia_size_t size_giant, index_baby, index_giant;
	bigint Mnew, a;
	long t;
	udigit u, hh;


	// Determine the number of combinations
	//
	nof_comb = number_of_combinations();

	// Positive result. Determine baby- and giantsteps.
	//

	if (nof_comb.is_positive()) {
		if (trace_list::info)
			std::cout << "\n# BG Test Candidates = " << nof_comb << std::endl;
	}
	// Zero. Trace uniquely determined by Elkies primes.
	//
	else if (nof_comb.is_zero()) {
		if (trace_list::info)
			std::cout << "\nNo BG algorithm necessary, group order uniquely determined" << std::endl;
		baby.kill();
		giant.kill();
		return true;
	}
	// Otherwise product of primes too small.
	//
	else {
		lidia_error_handler ("trace_list::split_baby_giant",
				     "Number of primes too small.");
		return false;
	}

	// Sort the list of relations. Increasing density.
	//
	tl_sort(0, l.size()-1);

	// last_index was set in number_of_combinations to the
	// index of the last modul in the product M of primes
	// such that M > 4 sqrt(p)
	//
	if (last_index < l.size()-1)
		l.remove_from(last_index+1, l.size()-last_index-1);

	// No giantsteps if there is only one relation.
	//
	if (last_index == 0) {
		index_baby = last_index+1;
		nof_baby.assign(nof_comb);
		M2.assign_one();
		giant.set_size(0);
	}

	// Otherwise determine number of baby/giantsteps, i.e.,
	// split moduli from 0 to last_index into babysteps and
	// giantsteps.
	//
	else {
		// Divide list l of trace_mods into babystep part and giantstep part,
		// i.e.,
		// determine index index_baby, such that the number of
		// combinations for the moduli with index 0 <= i < index_baby
		// is approx. sqrt(number of combinations).
		//
		// Set nof_baby = # babysteps
		//
		// We combine l[0, ..., index_baby-1] to babysteps,
		// l[index_baby, ..., index_giant-1] to giantsteps.
		//
		sqrt(up, nof_comb);
		h = l[0].get_size_of_trace_mod();

		if (h >= up) {
			index_baby = 1; // we know that there are at least 2 relations.
			nof_baby = h;
		}
		else {
			index_baby = 0;
			while (h < up) {
				index_baby++;
				multiply(h, h, l[index_baby].get_size_of_trace_mod());
			}
			divide(nof_baby, h, l[index_baby].get_size_of_trace_mod());
		}

		index_giant = l.get_size();


		// Chinese remaindering for giantsteps.
		//
		divide(h, nof_comb, nof_baby);
		h.longify(t);
		giant.set_capacity(t);

		// Initialize array giant of giansteps.
		// giant[0, ..., size_giant-1] = candidates for trace of l[index_baby]
		//
		// size_giant denotes the size of the giant array during chinese remaindering.
		//
		size_giant = l[index_baby].get_size_of_trace_mod();
		M2 = l[index_baby].get_modulus();

		for (i = 0; i < size_giant; i++)
			giant[i].assign(l[index_baby].get_element(i));

		// For all trace_mods index_baby <= j < index_giant belonging we apply chinese
		// remaindering, i.e.,
		// let msmall = l[j].get_modulus().
		// We combine candidates for trace mod msmall stored in l[j]
		// with candidates for trace mod M2 stored in giant array.
		//
		for (j = index_baby+1; j < index_giant; j++) {
			msmall = static_cast<udigit>(l[j].get_modulus());
			i = size_giant;
			multiply(Mnew, M2, msmall);
			u = invert_mod(remainder(M2, msmall), msmall);

			for (s = 1; s < l[j].get_size_of_trace_mod(); s++) {
				hh = static_cast<udigit>(l[j].get_element(s));

				for (k = 0; k < size_giant; k++) {
					a = giant[k];
					t = remainder((hh-a)*u, static_cast<long>(msmall));
					remainder(giant[i++], a + t*M2, Mnew);
				}
			}
			hh = static_cast<udigit>(l[j].get_element(0));

			for (k = 0; k < size_giant; k++) {
				a = giant[k];
				t = remainder((hh-a)*u, msmall);
				remainder(giant[k], a + t*M2, Mnew);
			}

			size_giant = i;
			M2.assign(Mnew);
		}

		//------ now giant list should be done completely
		//
		l.remove_from(index_baby, index_giant-index_baby); // clear giantstep part of TL
	}

	//---------------------------------------------------------------
	// Now it's babystep time.
	//
	// Same chinese remaindering as for giantsteps above.
	//
	nof_baby.longify(t);
	baby.set_capacity(t);
	e = l[0].get_size_of_trace_mod();
	M1.assign(l[0].get_modulus());

	for (i = 0; i < e; i++)
		baby[i].assign(l[0].get_element(i));

	if (index_baby > 1)
		for (j = 1; j < index_baby; j++) {
			msmall = static_cast<udigit>(l[j].get_modulus());
			i = e;
			multiply(Mnew, M1, msmall);
			u = invert_mod(remainder(M1, msmall), msmall);

			for (s = 1; s < l[j].get_size_of_trace_mod(); s++) {
				hh = static_cast<udigit>(l[j].get_element(s));

				for (k = 0; k < e; k++) {
					a = baby[k];
					t = remainder((hh-a)*u, static_cast<long>(msmall));
					remainder(baby[i++], a + t*M1, Mnew);
				}
			}
			hh = static_cast<udigit>(l[j].get_element(0));

			for (k = 0; k < e; k++) {
				a = baby[k];
				t = remainder((hh-a)*u, msmall);
				remainder(baby[k], a + t*M1, Mnew);
			}
			e = i;
			M1.assign(Mnew);
		}

	l.kill();

	transform_lists(baby, giant);
	return false;
}



void trace_list::transform_lists(sort_vector< bigint > & baby,
				 sort_vector< bigint > & giant)
{
	lidia_size_t i, e;
	bigint h1, h2, h3;


	// -----------------------------------------------------------------
	// transform babylist and giantlist (Volker Mueller, PhD Thesis page 145)
	//

	// babysteps
	//
	// h1 = 1 / (M2 M3) mod M1
	// h2 = floor(M1/2)
	// h3 = floor(-M1/2)
	//
	multiply (h1, M2, M3);
	xgcd_left(h1, h1, M1);

	shift_right(h2, M1, 1);
	negate(h3, h2);
	if (M1.is_odd())
		dec(h3);

	// Absolute smallest residues R1, -M1/2 < R1 <= M1/2
	// with R1 = (c - C3)/(M2 M3) mod M1.
	//
	e = baby.size();

	for(i = 0; i < e; i++) {
		remainder(baby[i], (baby[i]-C3)*h1, M1);
		if (baby[i] > h2)
			subtract(baby[i], baby[i], M1);
		else
			if (baby[i] <= h3)
				add(baby[i], baby[i], M1);
	}


	// giantsteps
	//
	// h1 = 1 / (M1 M3) mod M2
	//
	multiply(h1, M1, M3);
	xgcd_left(h1, h1, M2);


	// For each non-negative residue g, 0 <= g < M2, we have to check
	// all values -M2 <= g + k M2 <= M2 during giant step phase.
	// (see PhD Thesis of Volker Mueller, page 145).
	//
	//
	// If g != 0, then the values to check are g and g - M2.
	// If g == 0, then the values to check are -M2, 0, M2.
	//
	// We check g == 0 during babysteps. So it remains -M2 and M2
	// during giantsteps. For g = 0 we add the entry M2 to the giant list.
	// During giantsteps we treat this as a special case and also test for -M2.
	//
	// Note that in case of no giant steps we have M2 = 1.
	// Hence there is one g = 0 implicitely that must be checked.
	//
	e = giant.size();

	if (e == 0) {
		giant.set_capacity(1);
		giant[0] = M2;
	}
	else {
		// Smallest non-negative residue g = (c-C3)/(M1 M3) mod M2,
		// 0 <= g < M2; append additional M2 for g = 0;
		//
		for(i = 0; i < e; i++) {
			remainder(giant[i], (giant[i]-C3)*h1, M2);
			if (giant[i].is_negative())
				add(giant[i], giant[i], M2);

			if (giant[i].is_zero())
				giant[i] = M2;
		}
	}

	baby.sort();
	giant.sort();
}



//
//
//  Babystep / Giantstep algorithm
//
//

// first we define ec_tables which are used in BG algorithm.

class ec_table_entry
{
private:
	udigit h_value; // hash_value of x-coord of point
	int index; // index, -1 denotes empty
public:

	ec_table_entry()
	{
		index = -1;
	}

	~ec_table_entry() { }

	ec_table_entry & operator = (const ec_table_entry & e)
	{
		h_value = e.h_value;
		index = e.index;
		return *this;
	}

	void set(udigit h, int i)
	{
		h_value = h;
		index = i;
	}

	int get_index()
        {
		return index;
	}
	udigit get_hash()
        {
		return h_value;
	}
};



bigint trace_list::bg_search_for_order()
{
  if (number_of_combinations() < 100)
    return simple_search_for_order();

  sort_vector< bigint > baby, giant;
  bigint pn;

  power(pn, E.get_a6().characteristic(), E.degree_of_definition());
 
  if (this->split_baby_giant(baby, giant)) {
    // trace is the absolute smallest residue of c3 mod M3.
    return (pn + 1  - this->get_absolute_smallest_C3());
  }
  
  point< gf_element > P(E);
  P = E.random_point(E.degree_of_definition());
  
  bigint h;
  long number_baby = baby.size();
  long number_giant = giant.size();
  
  point< gf_element > H(P.get_curve()), H2(P.get_curve()), 
    H3(P.get_curve()), pn1mc3P(P.get_curve()), M2M3P(P.get_curve()), 
    H6(P.get_curve()), H_test(P.get_curve());
  long i, j, l, s, size_table, max_table = 0;
  udigit P_hash;
  bool use_internal_table;
  
  h = next_prime(2*number_baby);
  h.longify(size_table);
  
  ec_table_entry *HT;
  point< gf_element > *T_point = NULL; // for internal table of points used to
                                       // speed up computations

  HT = new ec_table_entry[size_table];

  
  // Initialize points.

  multiply(pn1mc3P, pn+1-C3, P);
  multiply(M2M3P, M2*M3, P);
  
  // For the babystep part, we have:

  H2 = M2M3P;

  // Compute table to allow fast multiplication n * H2
  // by table lookup for 0 < n 2s^2. See diploma thesis
  // of Markus Maurer, page 54.

  if (number_baby > 50) {
    bigint diff(0);
    use_internal_table = true;
    // Determine maximal gap for babysteps.
    // If g is a gap between two babysteps then only
    // g * H2 is computed in below and then added to
    // the previous babystep. So it is sufficient to build
    // the table only for gaps.
    //
    for (i = 1; i < number_baby; i++) {
      subtract(h, baby[i], baby[i-1]);
      if (h > diff)
	diff.assign(h);
    }
    sqrt(diff, diff);
    
    if (diff > 50)   // determine size of internal table, default value 2*50
      s = 50;
    else {
      diff.longify(s);
      s++;
    }
    
    // Build table:
    //
    // T[i]   = i * H2, 0 <= i <= s
    // T[s+i] = i * 2s * H2, 1 <= i <= s
    //
    T_point = new point< gf_element > [2*s+1];
    T_point[0] = H2;
    T_point[0].assign_zero(E);
    for (i = 1; i <= s; i++)
      add(T_point[i], T_point[i-1], H2);
    multiply_by_2(T_point[s+1], T_point[s]);
    
    for (i = s+2; i <= 2*s; i++)
      add(T_point[i], T_point[i-1], T_point[s+1]);
    max_table = (s+1)*s;
  }
  else
    use_internal_table = false;
  
  
  //****** Babysteps, store hash(x(+/- i*P))  *********************
  
  if (trace_list::info)
    std::cout << "\n# Babysteps = " << number_baby << std::endl;
  
  multiply(H3, baby[0], H2);
  subtract(H, pn1mc3P, H3);
  
  for (i = 0; i < number_baby; i++) {
    if (H.is_zero()) {
      h = pn + 1 - C3 - baby[i]*M2*M3;
      multiply(H_test, h, P);
      
      if (H_test.is_zero()) {
	delete[] HT;
	if (use_internal_table)
	  delete[] T_point;
	return h;
      }
    }
    
    P_hash = hash(H.get_x()); // insert into hash table
    j = P_hash % size_table;
    while (HT[j].get_index() != -1) {
      // HT[j] is not empty
      j++;
      if (j == size_table)
	j = 0;
    }
    HT[j].set(P_hash, i);
    
    if (i == number_baby-1)
      break;
    
    // Compute next babystep H = (p+1-c3) * H2 - baby[i+1] * H2
    //                         = H - (baby[i+1] - baby[i]) * H2
    //
    subtract(h, baby[i+1], baby[i]);
    
    // H3 = h * H2
    //
    // without table
    //
    if (!use_internal_table || (use_internal_table && h > max_table))
      multiply(H3, h, H2);
    
    // with table
    //
    else {
      h.longify(j);
      l = j % (2*s);
      j = j / (2*s);
      if (l > s)
	subtract(H3, T_point[s+j+1], T_point[2*s-l]);
      else
	if (j != 0)
	  add(H3, T_point[s+j], T_point[l]);
	else
	  H3 = T_point[l];
    }
    
    // H = H - H3
    //
    subtract(H, H, H3);
  }
  
  if (use_internal_table)
    delete[] T_point;
  
  
  //******************** Giantsteps **********************************
  
  // For each giant step 0 <= g < M2, we have to check all -M2 <= g + k M2 <= M2.
  // (see PhD Thesis of Volker Mueller, page 145).
  //
  // If g != 0, then the values to check are g and g - M2.
  // If g == 0, then the values to check are -M2, 0, M2.
  //  0 has already be tested during babysteps. So it remains -M2 and M2.
  //
  // Note that in case of no giant steps we have M2 = 1.
  // Hence there is one g = 0 implicitely that must be checked.
  //
  
  if (number_giant == 0 || giant.get_size() == 0) {
    number_giant = 2;
    giant.set_capacity(number_giant);
    giant[0] =  0;
    giant[1] = M2;
  }
  
  if (trace_list::info)
    std::cout << "# Giantsteps = " << number_giant << std::endl;
  
  // For the giant step part, we have:
  //
  // H2 = M1 M3 P
  //
  multiply(H2, M1*M3, P);
  multiply(H3, M2, H2);
  
  
  // Compute table to allow fast multiplication n * H2
  // by table lookup for 0 < n 2s^2. See diploma thesis
  // of Markus Maurer, page 54.
  //
  if (number_giant > 50) {
    bigint diff(0);
    use_internal_table = true;
    
    // Determine maximal gap for giantsteps.
    // If g is a gap between two giantsteps then only
    // g * H2 is computed in below and then added to
    // the previous point. So it is sufficient to build
    // the table only for gaps.
    //
    for (i = 1; i < number_giant; i++) {
      subtract(h, giant[i], giant[i-1]);
      if (h > diff)
	diff.assign(h);
    }
    sqrt(diff, diff);
    
    if (diff > 50)   // determine size of internal table, default value 2*50
      s = 50;
    else {
      diff.longify(s);
      s++;
    }
    
    // Build table:
    //
    // T[i]   = i * H2, 0 <= i <= s
    // T[s+i] = i * 2s * H2, 1 <= i <= s
    //
    T_point = new point< gf_element > [2*s+1];
    T_point[0] = H2;
    T_point[0].assign_zero(E);
    for (i = 1; i <= s; i++)
      add(T_point[i], T_point[i-1], H2);
    multiply_by_2(T_point[s+1], T_point[s]);
    
    for (i = s+2; i <= 2*s; i++)
      add(T_point[i], T_point[i-1], T_point[s+1]);
    max_table = (s+1)*s;
  }
  else
    use_internal_table = false;
  
  
  // For all giant steps g with 0 <= g <= M2, test all values
  // -M2 <= g + k * M2 <= M2.
  //
  // This is done as follows. For g != 0, M2 test
  //
  // g * H2      for match with babysteps during first loop.
  // (g-M2) * H2 for match with babysteps during second loop.
  //
  // We ignore g == 0, because it has already been tested during babysteps.
  // For g == M2, we test M2 and -M2.
  //
  j = 0;
  while (giant[j].is_zero())
    j++;
  multiply(H, giant[j], H2);
  
  for (; j < number_giant; j++) {
    // H = giant[j] * H2.
    //
    // If H is zero, group order probably is giant[j] M1 M3.
    //
    if (H.is_zero()) {
      h = giant[j]*M1*M3;
      multiply(H_test, h, P);
      
      if (H_test.is_zero()) {
	delete[] HT;
	if (use_internal_table)
	  delete[] T_point;
	
#ifdef DEBUG
	std::cout << "trace_list::bg_search_for_order::giant[" << j << "] entry yields result" << std::endl;
	std::cout << "giant[" << j << "] = " << giant[j]
		  << std::endl;
#endif
	return h;
      }
    }
    
    // Test for giant[j].
    //
    // Search for x-coordinate of H in the hash table.
    // Maybe a match for babysteps.
    //
    P_hash = hash(H.get_x());
    l = P_hash % size_table;
    
    while (HT[l].get_index() != -1) {
      // Hash values are the same. Reconstruct babystep in H3
      // and check for match.
      //
      if (HT[l].get_hash() == P_hash) {
	multiply(H3, baby[HT[l].get_index()], M2M3P);
	subtract(H3, pn1mc3P, H3);
	
                                // Match. Done.
                                //
	if (H3 == H) {
	  h = pn+1-C3-M2*M3*baby[HT[l].get_index()
	  ] - giant[j]*M1*M3;
	  
	  multiply(H_test, h, P);
	  if (H_test.is_zero()) {
	    
#ifdef DEBUG
	    std::cout << "trace_list::bg_search_for_order:: -giant[" << j << "] entry yields result" << std::endl;
	    std::cout << "giant[" << j << "] = " << giant[j] << std::endl;
	    std::cout << "baby[" << HT[l].get_index() << "] = " << baby[HT[l].get_index()] << std::endl;
#endif
	    
	    if (use_internal_table)
	      delete[] T_point;
	    delete [] HT;
	    
	    return h;
	  }
	}
	
      // Same x-coordinates, but different y-coordinates. Inverse points.
                                // Match for -H3.
                                //
	if (H3.get_x() == H.get_x()) {
	  h = pn+1-C3-M2*M3*baby[HT[l].get_index()
	  ] + giant[j]*M1*M3;
	  
	  multiply(H_test, h, P);
	  if (H_test.is_zero()) {
#ifdef DEBUG
	    std::cout << "trace_list::bg_search_for_order:: +giant[" << j << "] entry yields result" << std::endl;
	    std::cout << "giant[" << j << "] = " << giant[j] << std::endl;
	    std::cout << "baby[" << HT[l].get_index() << "] = " << baby[HT[l].get_index()] << std::endl;
#endif
	    
	    if (use_internal_table)
	      delete[] T_point;
	    delete [] HT;
	    
	    return h;
	  }
	}
      }
      
      // Next hash value.
      //
      l++;
      
      if (l == size_table)
	l = 0;
    }
    
    // If giant[j] == M2, test for -M2.
    // Otherwise, test for giant[j] - M2.
    //
    if (giant[j] == M2) {
      // H6 = -M2 * H2
      //
      negate(H6, H);
      giant[j].assign_zero();
    }
    else {
      // H6 = (giant[j] - M2) * H2
      //
      subtract(H6, H, H3);
      if (H6.is_zero()) {
	h = (giant[j]-M2)*M1*M3;
	
	multiply(H_test, h, P);
	if (H_test.is_zero()) {
#ifdef DEBUG
	  std::cout << "trace_list::bg_search_for_order::giant[" << j << "] entry yields result" << std::endl;
	  std::cout << "giant[" << j << "] = " << 
	    giant[j] << std::endl;
#endif
	  
	  if (use_internal_table)
	    delete[] T_point;
	  delete[] HT;
	  
	  return h;
	}
      }
    }
    
    P_hash = hash(H6.get_x());
    l = P_hash % size_table;
    
    while (HT[l].get_index() != -1) {
      // hash values are the same
      if (HT[l].get_hash() == P_hash) {
	multiply(H3, baby[HT[l].get_index()], M2M3P);
	subtract(H3, pn1mc3P, H3);
	
	if (H3 == H6) {
	  h = pn+1-C3-M2*M3*baby[HT[l].get_index()
	  ] - (giant[j]-M2)*M1*M3;
	  multiply(H_test, h, P);
	  if (H_test.is_zero()) {
#ifdef DEBUG
	    std::cout << "trace_list::bg_search_for_order::giant[" << j
		      << "] entry yields result" << std::endl;
	    std::cout << "giant[" << j << "] = " << giant[j] << std::endl;
#endif
	    
	    if (use_internal_table)
	      delete[] T_point;
	    delete [] HT;
	    
	    return h;
	  }
	}
	
	if (H3.get_x() == H6.get_x()) {
	  h = pn+1-C3-M2*M3*baby[HT[l].get_index()
	  ] + (giant[j]-M2)*M1*M3;
	  multiply(H_test, h, P);
	  if (H_test.is_zero()) {
#ifdef DEBUG
	    std::cout << "trace_list::bg_search_for_order::giant[" << j
		      << "] entry yields result" << std::endl;
	    std::cout << "giant[" << j << "] = " << giant[j] << std::endl;
#endif
	    
	    if (use_internal_table)
	      delete[] T_point;
	    delete [] HT;
	    
	    return h;
	  }
	}
      }
      l++;
      if (l == size_table)
	l = 0;
    }
    
    if (giant[j].is_zero())
      giant[j] = M2;
    
    // No more giantsteps and no match. Hmm :-(
    //
    //
    if (j+1 == number_giant) {
      std::cout << "trace_list::bg_search_for_order::P = " << 
	P << std::endl;
      lidia_error_handler ("trace_list",
			   "bg_search_for_order::Last giant step passed. No match.");
    }
    
    // Compute next giantstep, H = giant[j+1] * H2
    //                         = giant[j] * H2 + (giant[j+1] - giant[j]) * H2
    //                           = H + h * H2
    //
    subtract(h, giant[j+1], giant[j]);
    
    // H6 = h * H2
    //
    // without table
    //
    if (!use_internal_table || (use_internal_table && h > max_table))
      multiply(H6, h, H2);
    
    // with table
    //
    else {
      h.longify(i);
      l = i % (2*s);
      i = i / (2*s);
      
      if (l > s)
	subtract(H6, T_point[s+i+1], T_point[2*s-l]);
      else
	if (i != 0)
	  add(H6, T_point[s+i], T_point[l]);
	else
	  H6 = T_point[l];
    }
    
    // H = H + H6
    //
    add(H, H, H6);
  }
  
  if (use_internal_table)
    delete[] T_point;
  
  delete[] HT;
  return (-1);
}




//
//  Name: simple_search_for_order
//
//  Computes all possible traces t, generates random point P, and
//  tests for (p+1-t) * P = O. If no match is found, the function
//  returns -1. Otherwise it returns (p+1-t0), where t0 is the first
//  trace candidate for which (p+1-t0) * P = O was found.
//
//  Only used for debugging purposes, because bg_search_for_order is
//  faster.
//

bigint trace_list ::simple_search_for_order()
{
  // Compute all possible traces.
  //
  
  bigint pn;

  power (pn, E.get_a6 ().characteristic (), E.degree_of_definition ());

	// Determine the number of combinations
	//
	bigint nof_comb = number_of_combinations();

	// Positive result. Determine traces.
	//
	if (nof_comb.is_positive()) {
		if (trace_list::info)
			std::cout << "\n# Test Candidates for trace = " << nof_comb << std::endl;
	}

	// Zero. Trace uniquely determined by Elkies primes.
	//
	else if (nof_comb.is_zero()) {
		if (trace_list::info)
			std::cout << "\nNo search necessary, group order uniquely determined";

		return (pn + 1  - this->get_absolute_smallest_C3());
	}
	// Otherwise product of primes too small.
	//
	else {
		lidia_error_handler ("trace_list::simple_search",
				     "Number of primes too small.");
		return -2;
	}


	// last_index was set in number_of_combinations to the
	// index of the last modul in the product M of primes
	// such that M > 4 sqrt(p)
	//
	base_vector< lidia_size_t > index;
	lidia_size_t nof_rows;
	lidia_size_t i;

	// combine traces t[index[i]] from l[i], i = 0, ..., l.get_size-1
	//
	nof_rows = last_index + 1;
	index.set_capacity(nof_rows);

	for (i = 0; i < nof_rows; i++)
		index[i] = 0;

	// Initialize crt.
	//
	base_vector< sdigit > small_moduli;
	bigint N;

	small_moduli.set_capacity(nof_rows);
	N.assign_one();
	for (i = 0; i < nof_rows; i++) {
		small_moduli[i] = static_cast<sdigit>(l[i].get_modulus());
		multiply (N, N, l[i].get_modulus());
	}

	crt_table ct;
	ct.init (small_moduli);

	crt C;
	C.init(ct);


	bigint M, h2, h3;

	multiply (M, M3, N);

	shift_right(h2, M, 1);
	negate(h3, h2);
	if (M.is_odd())
		h3--;

	// Determine point P and precompute (p+1) * P.
	//
	point< gf_element > P(E), Q(E), R(E);
	bigint pp1;

	P = E.random_point(E.degree_of_definition ());
	pp1 = pn+1;
	multiply (Q, pp1, P);

#ifdef DEBUG
	std::cout << "Random point P = " << P << std::endl;
#endif

	// Use crt to compute trace candidate t and test (p+1-t) * P = O.
	//
	lidia_size_t completed_rows;
	lidia_size_t step;
	bigint t;

	completed_rows = 0;
	step = 1;

	while (completed_rows < nof_rows) {
		if (info >= 1000) {
			if (step % 1000 == 0)
				std::cout << " " << step << std::flush;
			step++;
		}

		// combine traces t[index[i]] from l[i], i = 0, ..., l.get_size-1,
		// and C3
		//
		for (i = 0; i < nof_rows; i++)
			C.combine(static_cast<sdigit>(l[i].get_element(index[i])), i);

		C.get_result(t);
		t = chinese_remainder (C3, M3, t, N);

		// absolute smallest residues t, -M/2 < t <= M/2
		//
		if (t > h2)
			subtract(t, t, M);
		else if (t <= h3)
			add(t, t, M);

#ifdef DEBUG
		std::cout << "Testing trace candidate t = " << t << std::endl;
#endif

		// compute point multiple R = Q - t * P
		//
		multiply (R, t, P);
		subtract (R, Q, R);

		if (R.is_zero())
			return pp1-t;

		// next combination
		//
		i = 0;
		completed_rows = 0;

		while (i < nof_rows) {
			// index[i] += 1 modulo l[i].get_size().
			// If at the end, handle next row.
			//
			if (index[i] == l[i].get_size_of_trace_mod()-1) {
				index[i] = 0;
				completed_rows++;
				i++;
			}
			else {
				index[i]++;
				i = nof_rows;
			}
		}
	}

	lidia_error_handler ("trace_list::simple_search",
			     "No match found.");
	return -1;
}



//
//  Name: baby_giant_lists_correct
//
//  ec_order is the group order of the curve E.
//
//  The functions computes the babystep and giantstep lists. Then it
//  derives the trace of the E from its group order and checks whether
//  the correct values that should yield a match during bg algorithm are
//  contained in the babystep list and in the giantstep list.
//
//  If those values can be found in the babystep and giantstep lists, the
//  function returns true. Otherwise it returns false.
//

bool trace_list ::baby_giant_lists_correct (const bigint & ec_order)
{
	sort_vector< bigint > baby, giant;

	//
	// Compute babystep and giantstep lists.
	//
	if (this->split_baby_giant(baby, giant)) {
		// trace is the absolute smallest residue of c3 mod M3.
		return true;
	}

	//
	// Determine trace.
	//
	bigint pn, t;

	power (pn, E.get_a6 ().characteristic (), E.degree_of_definition ());
	t  = pn + 1 - ec_order;

	bigint h1, h2, h3;

	//
	// Determine element of babystep list.
	//
	// h1 = 1 / (M2 M3) mod M1
	// h2 = floor(M1/2)
	// h3 = floor(-M1/2)
	//
	multiply (h1, M2, M3);
	xgcd_left(h1, h1, M1);

	shift_right(h2, M1, 1);
	negate(h3, h2);
	if (M1.is_odd())
		dec(h3);

	// Absolute smallest residues R1, -M1/2 < R1 <= M1/2
	// with R1 = (c - C3)/(M2 M3) mod M1.
	//
	bigint B;

	remainder(B, (t-C3)*h1, M1);
	if (B > h2)
		subtract(B, B, M1);
	else if (B <= h3)
		add(B, B, M1);


	//
	// Search for B in babysteps.
	//
	lidia_size_t i, n;
	bool found;

	n     = baby.get_size();
	found = false;

	i = 0;
	while ((i < n) && (!found))
		if (baby[i] == B)
			found = true;
		else
			i++;

	if (info >= 1000) {
		std::cout << "trace_list::baby_giant_lists_correct::";
		if (found == true) {
			std::cout << "baby[" << i << "] is the correct entry." << std::endl;
			std::cout << "baby[" << i << "] = " << baby[i] << std::endl;
		}
		else
			std::cout << "No correct babystep entry found." << std::endl;
	}

	if (found == false)
		return false;


	//
	// Determine element of gianstep list.
	//
	// h1 = 1 / (M1 M3) mod M2
	//
	multiply(h1, M1, M3);
	xgcd_left(h1, h1, M2);

	bigint G;

	remainder(G, (t-C3)*h1, M2);
	if (G.is_negative())
		add(G, G, M2);

	if (G.is_zero())
		G = M2;

	//
	// Search for G in giantsteps.
	//
	n     = giant.get_size();
	found = false;

	i = 0;
	while ((i < n) && (!found))
		if (giant[i] == B)
			found = true;
		else
			i++;

	if (info >= 1000) {
		std::cout << "trace_list::baby_giant_lists_correct::";
		if (found == true) {
			std::cout << "giant[" << i << "] is the correct entry." << std::endl;
			std::cout << "giant[" << i << "] = " << giant[i] << std::endl;
		}
		else
			std::cout << "No correct giantstep entry found." << std::endl;
	}

	return found;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
