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
//	Author	: Markus Maurer (MM)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/quadratic_ideal_power_product.h"
#include	"LiDIA/quadratic_number_standard.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//#define QIPP_DEBUG_REFINEMENT

//
// constructor / destructor
//
quadratic_ideal_power_product
::quadratic_ideal_power_product ()
	: base_power_product< quadratic_ideal, bigint > ()
{
	debug_handler("quadratic_ideal_power_product",
		      "quadratic_ideal_power_product");
}



quadratic_ideal_power_product
::~quadratic_ideal_power_product ()
{
	debug_handler("quadratic_ideal_power_product",
		      "~quadratic_ideal_power_product");
}



//
// assignment
//
quadratic_ideal_power_product &
quadratic_ideal_power_product::operator = (const quadratic_ideal_power_product & x)
{
	debug_handler("quadratic_ideal_power_product",
		      "operator = (const quadratic_ideal_power_product&)");

	base_power_product< quadratic_ideal, bigint >::assign(x);
	return *this;
}



void quadratic_ideal_power_product
::assign (const base_vector< quadratic_number_standard > & b,
	  const base_vector< bigint > & e)
{
	debug_handler ("quadratic_ideal_power_product",
		       "assign(const base_vector< quadratic_number_standard > &, "
		       "const base_vector< bigint > &)");

	if (b.get_size() != e.get_size()) {
		lidia_error_handler ("quadratic_ideal_power_product::assign"
				     "(const base_vector< quadratic_number_standard > &, "
				     "const base_vector< bigint > &",
				     "Different sizes.");
	}
	else if (b.get_size() == 0) {
		component.set_capacity(0);
	}
	else {
		lidia_size_t i, j, n;
		quadratic_ideal O;

		O.assign_one(b[0].get_order());
		n = b.get_size();
		component.set_capacity(n);

		for (j = 0, i = 0; i < n; i++) {
			if (!e[i].is_zero()) {
				multiply (component[j].left(), O, b[i]);
				component[j].right() = e[i];
				j++;
			}
		}
	}
}



//
// factor refinement
//
// Returns 1, if a ascent to an over order occurs and 0 otherwise.
//
int quadratic_ideal_power_product
::factor_refinement()
{
	debug_handler ("quadratic_ideal_power_product",
		       "factor_refinement()");

	lidia_size_t    i, j, k, m, n, new_i;
	quadratic_ideal gcd, inv_gcd;
	bigint f;

	int overorder = 0;

#ifdef QIPP_DEBUG_REFINEMENT
	std::cout << "qipp::factor_refinement::this = " << *this << std::endl;
#endif

	if (component.get_size() == 0)
		return overorder;

	// Separate numerator and denominator.
	//
	const quadratic_order & O1 = component[0].left().get_order();
	n = component.get_size();

	for (j = 0, i = n; j < n; j++) {
		// get denominator f of ideal
		//
		f = component[j].left().denominator();

		if (!f.is_one()) {
			// multiply ideal by denominator -> integral
			//
			component[j].left().multiply_by_denominator();

			// append (f*O)^{-e} to power product
			//
			component[i].left().assign_principal(f, bigint(0), O1);
			component[i].right().assign(-component[j].right());
			i++;
		}
	}

#ifdef QIPP_DEBUG_REFINEMENT
	std::cout << "After separating num and den: this = " << *this << std::endl;
#endif

	//
	// Remove ones
	//
	n = component.get_size();

	for (j = 0, k = 0; j < n; j++) {
		// Is the component one ?
		//
		if (!component[j].left().is_one() &&
		    !component[j].right().is_zero()) {
			// If there have been ones in between,
			// overwrite them.
			//
			if (j != k)
				component[k] = component[j];

			// index for next non-one element.
			k++;
		}
	}
	component.set_size(k);

#ifdef QIPP_DEBUG_REFINEMENT
	std::cout << "After removing ones, this = " << *this << std::endl;
#endif

	//
	// Start refinement
	//

	// Refine with elements (i,j).
	//
	for (i = 0; i < component.get_size()-1; i = new_i) {
#ifdef QIPP_DEBUG_REFINEMENT
		std::cout << "Traversing factorization with component[" << i << "] = " << component[i] << std::endl;
#endif

		for (j = i+1; j < component.get_size(); j++) {

#ifdef QIPP_DEBUG_REFINEMENT
			std::cout << "Refinement with component[" << j << "] = " << component[j] << std::endl;
#endif

			// Check for non-trivial gcd
			//
			add(gcd, component[i].left(), component[j].left());

#ifdef QIPP_DEBUG_REFINEMENT
			std::cout << "gcd = " << gcd << std::endl;
#endif

			if (!gcd.is_one()) {
				// If gcd is not invertible, lift ideals.
				//
				f = gcd.conductor();

#ifdef QIPP_DEBUG_REFINEMENT
				std::cout << "conductor = " << f << std::endl;
#endif
				if (!f.is_one()) {
					overorder = 1;

					// Determine overorder
					//
					lidia_error_handler ("quadratic_ideal_power_product"
							     "::factor_refinement()",
							     "Cannot handle non-maximal orders yet.");
				}

				inverse(inv_gcd, gcd);

#ifdef QIPP_DEBUG_REFINEMENT
				std::cout << "Inverse of gcd = " << inv_gcd << std::endl;
#endif
				// Now divide elements i and j by gcd
				// and append gcd.
				//
				component[i].left() *= inv_gcd;
				component[j].left() *= inv_gcd;
				m = component.get_size();
				component[m].left().assign(gcd);
				component[m].right().assign(component[i].right()+
							    component[j].right());
#ifdef QIPP_DEBUG_REFINEMENT
				std::cout << "New component[" << i << "] = " << component[i] << std::endl;
				std::cout << "New component[" << j << "] = " << component[j] << std::endl;
				std::cout << "New component[" << m << "] = " << component[m] << std::endl;
#endif

				// Stop, if element i is one.
				//
				if (component[i].left().is_one())
					j = n;
			}
		}

#ifdef QIPP_DEBUG_REFINEMENT
		std::cout << "Stopped refinement with component " << i << std::endl;
		std::cout << "Factorization = " << *this << std::endl;
#endif
		if (component[i].left().is_one())
			new_i = i;
		else
			new_i = i+1;


		// Remove ones.
		//
		n = component.get_size();

		for (j = i, k = i; j < n; j++) {
			// Is the component one ?
			//
			if (!component[j].left().is_one() &&
			    !component[j].right().is_zero()) {
				// If there have been ones in between,
				// overwrite them.
				//
				if (j != k)
					component[k] = component[j];

				// index for next non-one element.
				k++;
			}
		}
		component.set_size(k);

#ifdef QIPP_DEBUG_REFINEMENT
		std::cout << "Removed ones, factorization = " << *this << std::endl;
#endif
	}

	return overorder;
}



//
// is_gcd_free
//
bool quadratic_ideal_power_product
::is_gcd_free() const
{
	debug_handler ("quadratic_ideal_power_product",
		       "is_gcd_free()");

	lidia_size_t i, j , n;

	n = component.get_size();

	// special cases n == 0,1
	//
	if (n == 0 || n == 1)
		return true;

	// general case
	//
	quadratic_ideal gcd;

	for (i = 0; i < n-1; i++) {
		for (j = i+1; j < n; j++) {
			add(gcd, component[i].left(), component[j].left());

			if (!gcd.is_one()) {
				std::cout << "component[" << i << "] = " << component[i] << std::endl;
				std::cout << "component[" << j << "] = " << component[j] << std::endl;
				std::cout << "gcd = " << gcd << std::endl;
				return false;
			}
		}
	}

	return true;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
