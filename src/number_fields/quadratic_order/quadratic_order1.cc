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
//	Author	: Michael Jacobson (MJ)
//	Changes	: See CVS log
//
//==============================================================================================


#ifdef HAVE_CONFIG_H
# include	"config.h"
#endif
#include	"LiDIA/quadratic_order.h"
#include	"LiDIA/qi_class.h"
#include	"LiDIA/qi_class_real.h"
#include	"LiDIA/number_fields/qo_list.h"



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif



//<TPf>
qo_node *qi_class::current_order = NULL;
bigint qi_class::Delta = 0;
bigint qi_class::rootD = 0;
xbigfloat qi_class::rd = 0.0;
bigint qi_class::PEA_L = 0;
bigint qi_class::ordmult = 0;
rational_factorization qi_class::omfact;
int qi_class::info = 0;
int qi_class::do_verify = 0;
long qi_class::xprec = 0;
//</TPf>


quadratic_order quadratic_order::zero_QO;

qo_list& quadratic_order::qo_l() {
    static qo_list quadratic_order_list;

    return quadratic_order_list;
}

int quadratic_order::info = 0;
int quadratic_order::do_verify = 0;

// FOR TABLE GENERATION
bool quadratic_order::special2 = false;


const long OQvals[] = {1663, 4177, 7523, 11657, 16477, 21991, 28151, 34913, 42293,
		       50261, 58787, 67883, 77527, 87719, 98419, 109661, 121403, 133649, 146407, 159667};

const long OQvals_cnum[] = {2269, 5741, 10427, 16183, 22901, 30631,
			    39209, 48731, 59063, 70237, 82223, 95009, 108571, 122921, 137983, 153817, 170341,
			    187631, 205589, 224261};

const long OQvals_cfunc[] = {630000, 825000, 995000, 1145000,
			     1290000, 1425000, 1555000, 1675000, 1795000, 1910000, 2020000, 2125000, 2230000,
			     2335000, 2435000, 2530000, 2625000, 2720000, 2810000, 2900000};



//
// constructor
//     - sets resets all the class members which are vectors
//

quadratic_order::quadratic_order()
{
	debug_handler("quadratic_order", "entering quadratic_order()");

	SPL.resize(2, 70000);
	CL.set_mode(EXPAND);
	gens.set_mode(EXPAND);
	fact_base.set_mode(EXPAND);
	contributors.set_mode(EXPAND);
	minima_qn.set_mode(EXPAND);
	CL.reset();
	gens.reset();
	fact_base.reset();
	contributors.reset();
	minima_qn.reset();
	subused = false;
	U.resize(1, 1);
	UINV.resize(1, 1);

	R_xbig.assign_zero(); // MM

	if (!qi_class::get_current_order_ptr())
		qi_class::set_current_order(zero_QO);

	debug_handler("quadratic_order", "leaving quadratic_order()");
}



//
// constructor
//     - sets resets all the class members which are vectors and initializes
//       with discriminant D.  If D is not a quadratic discriminant, 0 will
//       be used.
//

quadratic_order::quadratic_order(const long D)
{
	debug_handler("quadratic_order", "entering quadratic_order(long)");

	long oprec = bigfloat::get_precision();
	bigfloat temp;

	// initialization
	SPL.resize(2, 70000);
	CL.set_mode(EXPAND);
	gens.set_mode(EXPAND);
	fact_base.set_mode(EXPAND);
	contributors.set_mode(EXPAND);
	minima_qn.set_mode(EXPAND);
	CL.reset();
	gens.reset();
	fact_base.reset();
	contributors.reset();
	minima_qn.reset();
	DET = 1;
	AMAT.kill();
	AMAT.set_zero_element(0);
	AMAT.set_representation(matrix_flags::sparse_representation);
	AMAT.set_orientation(matrix_flags::row_oriented);
	subused = false;

	RHNF.kill();
	RHNF.set_representation(matrix_flags::sparse_representation);
	U.kill();
	UINV.kill();
	U.resize(1, 1);
	UINV.resize(1, 1);

	if (is_quadratic_discriminant(D))
		Delta.assign(D);
	else
		Delta.assign_zero();

	prec = static_cast<long>(decimal_length(Delta) >> 1) + 10;
	xprec = static_cast<long>(std::ceil(static_cast<double>(prec)*std::log(10.0L) / std::log(2.0L)));

	bigfloat::set_precision(prec);

	R.assign_zero();
	R_xbig.assign_zero(); // MM
	h.assign_zero();
	L.assign_zero();
	Cfunc.assign_zero();
	FI.assign_zero();
	Q = 0;

	if (Delta.is_le_zero())
		newton_root(nucomp_bound, abs(Delta >> 2), 4);
	else {
		nucomp_bound.assign_zero();
		min_log = log((sqrt(bigfloat(Delta-4)) + sqrt(bigfloat(Delta))) /
			      bigfloat(2.0)) / bigfloat(2.0);
	}

	bigfloat::set_precision(oprec);

	if (!qi_class::get_current_order_ptr())
		qi_class::set_current_order(zero_QO);

	qo_l().set_last(this);

	debug_handler("quadratic_order", "leaving quadratic_order(long)");
}



//
// constructor
//     - sets resets all the class members which are vectors and initializes
//       with discriminant D.  If D is not a quadratic discriminant, 0 will
//       be used.
//

quadratic_order::quadratic_order(const bigint & D)
{
	debug_handler("quadratic_order", "entering quadratic_order(bigint)");

	long oprec = bigfloat::get_precision();
	bigfloat temp;

	// initialization
	SPL.resize(2, 70000);
	CL.set_mode(EXPAND);
	gens.set_mode(EXPAND);
	fact_base.set_mode(EXPAND);
	contributors.set_mode(EXPAND);
	minima_qn.set_mode(EXPAND);

	CL.reset();
	gens.reset();
	fact_base.reset();
	contributors.reset();
	minima_qn.reset();
	DET = 1;
	AMAT.kill();
	AMAT.set_zero_element(0);
	AMAT.set_representation(matrix_flags::sparse_representation);
	AMAT.set_orientation(matrix_flags::row_oriented);
	subused = false;

	RHNF.kill();
	RHNF.set_representation(matrix_flags::sparse_representation);
	U.kill();
	UINV.kill();
	U.resize(1, 1);
	UINV.resize(1, 1);

	if (is_quadratic_discriminant(D))
		Delta.assign(D);
	else
		Delta.assign_zero();

	prec = static_cast<long>(decimal_length(Delta) >> 1) + 10;
	xprec = static_cast<long>(std::ceil(static_cast<double>(prec)*std::log(10.0L) / std::log(2.0L)));

	bigfloat::set_precision(prec);

	R.assign_zero();
	R_xbig.assign_zero(); // MM
	h.assign_zero();
	L.assign_zero();
	Cfunc.assign_zero();
	FI.assign_zero();
	Q = 0;

	if (Delta.is_le_zero())
		newton_root(nucomp_bound, abs(Delta >> 2), 4);
	else {
		nucomp_bound.assign_zero();
		min_log = log((sqrt(bigfloat(Delta-4)) + sqrt(bigfloat(Delta))) /
			      bigfloat(2.0)) / bigfloat(2.0);
	}

	bigfloat::set_precision(oprec);

	if (!qi_class::get_current_order_ptr())
		qi_class::set_current_order(zero_QO);

	qo_l().set_last(this);

	debug_handler("quadratic_order", "leaving quadratic_order(bigint)");
}



//
// constructor
//     - sets resets all the class members which are vectors and initializes
//       with a copy of QO.
//

quadratic_order::quadratic_order(const quadratic_order & QO)
{
	debug_handler("quadratic_order", "quadratic_order(quadratic_order)");

	long oprec = bigfloat::get_precision();

	SPL.resize(2, 70000);
	CL.set_mode(EXPAND);
	gens.set_mode(EXPAND);
	fact_base.set_mode(EXPAND);
	contributors.set_mode(EXPAND);
	minima_qn.set_mode(EXPAND);

	Delta.assign(QO.Delta);
	prec = QO.prec;
	xprec = QO.xprec;
	bigfloat::set_precision(prec);

	R.assign(QO.R);
	R_xbig.assign(QO.R_xbig); //MM
	R_xbig_prec = QO.R_xbig_prec;
	h.assign(QO.h);
	L.assign(QO.L);
	Cfunc.assign(QO.Cfunc);
	CL = QO.CL;
	gens = QO.gens;
	hfact.assign(QO.hfact);
	disc_fact.assign(QO.disc_fact);
	fact_base = QO.fact_base;
	contributors = QO.contributors;
	minima_qn = QO.minima_qn;
	fund_unit = QO.fund_unit;
	DET = QO.DET;
	AMAT = QO.AMAT;

	subused = QO.subused;
	prin_list = QO.prin_list;
	FI = QO.FI;
	Q = QO.Q;
	nucomp_bound = QO.nucomp_bound;
	RHNF = QO.RHNF;
	U = QO.U;
	UINV = QO.UINV;

	bigfloat::set_precision(oprec);

	qo_l().set_last(this);
}



//
// destuctor
//     - deletes this order from the list of current orders.
//

quadratic_order::~quadratic_order()
{
	debug_handler("quadratic_order", "~quadratic_order");

	TR.kill();
	qo_l().nullify(this);
}



//
// quadratic_order::verbose()
//
// Task:
//      sets the verbosity of commands.  Currently, the following levels are
//      supported:
//         0 - nothing
//         1 - run times and some run-time data (data only for subexp algs.)
//         >1 - debug information for subexponential algorithms
//

void
quadratic_order::verbose(int state)
{
	debug_handler("quadratic_order", "verbose");

	if (state <= 0) {
		info = 0;
	}
	else {
		info = state;
	}

	//MM, seems to be unused, qo_info = info;
	qo_info = info;
}



//
// quadratic_order::verification()
//
// Task:
//      sets the level of verifications.  Currently, the following levels are
//      supported:
//         0 - nothing
//         1-3 - varies levels of verification
//

void
quadratic_order::verification(int level)
{
	debug_handler("quadratic_order", "verification");

	if (level <= 0)
		do_verify = 0;
	else
		do_verify = level;
}



//
// quadratic_order::assign(long)
//
// Task:
//    set to the quadratic order of discriminant D.  If D is not a quadratic
//    discriminant, 0 the quadratic order will be set to the empty order
//    (discriminant 0) and false will be returned.
//

bool
quadratic_order::assign(const long D)
{
	debug_handler("quadratic_order", "assign(long)");

	return assign(bigint(D));
}



//
// quadratic_order::assign(bigint)
//
// Task:
//    set to the quadratic order of discriminant D.  If D is not a quadratic
//    discriminant, 0 the quadratic order will be set to the empty order
//    (discriminant 0) and false will be returned.
//

bool
quadratic_order::assign(const bigint & D)
{
	debug_handler("quadratic_order", "assign(bigint)");

	long oprec = bigfloat::get_precision();
	bigfloat temp;
	bool is_disc;

	is_disc = is_quadratic_discriminant(D);

	if (is_disc)
		Delta.assign(D);
	else
		Delta.assign_zero();

	prec = static_cast<long>(decimal_length(Delta) >> 1) + 10;
	xprec = static_cast<long>(std::ceil(static_cast<double>(prec)*std::log(10.0L) / std::log(2.0L)));

	bigfloat::set_precision(prec);

	CL.set_mode(EXPAND);
	gens.set_mode(EXPAND);
	fact_base.set_mode(EXPAND);
	contributors.set_mode(EXPAND);
	minima_qn.set_mode(EXPAND);

	R.assign_zero();
	R_xbig.assign_zero(); // MM
	h.assign_zero();
	L.assign_zero();
	Cfunc.assign_zero();
	CL.reset();
	gens.reset();
	fact_base.reset();
	contributors.reset();
	TR.kill();
	minima_qn.reset();
	fund_unit.reset();
	DET = 1;
	AMAT.kill();
	AMAT.set_zero_element(0);
	AMAT.set_representation(matrix_flags::sparse_representation);
	AMAT.set_orientation(matrix_flags::row_oriented);

	subused = false;
	hfact = rational_factorization();
	disc_fact = rational_factorization();
	RHNF.kill();
	RHNF.set_representation(matrix_flags::sparse_representation);
	U.kill();
	UINV.kill();
	U.resize(1, 1);
	UINV.resize(1, 1);
	prin_list.empty();
	FI.assign_zero();
	Q = 0;

	if (Delta.is_le_zero())
		newton_root(nucomp_bound, abs(Delta >> 2), 4);
	else {
		nucomp_bound.assign_zero();
		min_log = log((sqrt(bigfloat(Delta-4)) + sqrt(bigfloat(Delta))) /
			      bigfloat(2.0)) / bigfloat(2.0);
	}

	sieve.reset();

	bigfloat::set_precision(oprec);

	// set last-used quadratic_order
	qo_l().set_last(this);

	return is_disc;
}



//
// quadratic_order::assign(quadratic_order)
//
// Task:
//    set to a copy of QO.
//

void
quadratic_order::assign(const quadratic_order & QO)
{
	debug_handler("quadratic_order", "assign(quadratic_order)");

	long oprec = bigfloat::get_precision();

	Delta.assign(QO.Delta);
	prec = QO.prec;
	xprec = QO.xprec;
	bigfloat::set_precision(prec);

	R.assign(QO.R);
	R_xbig.assign(QO.R_xbig); //MM
	R_xbig_prec = QO.R_xbig_prec;
	h.assign(QO.h);
	L.assign(QO.L);
	Cfunc.assign(QO.Cfunc);
	CL = QO.CL;
	gens = QO.gens;
	fact_base = QO.fact_base;
	contributors = QO.contributors;
	minima_qn = QO.minima_qn;
	fund_unit = QO.fund_unit;
	DET = QO.DET;
	AMAT = QO.AMAT;

	subused = QO.subused;
	hfact.assign(QO.hfact);
	disc_fact.assign(QO.disc_fact);
	prin_list = QO.prin_list;
	FI = QO.FI;
	Q = QO.Q;
	nucomp_bound = QO.nucomp_bound;
	RHNF = QO.RHNF;
	U = QO.U;
	UINV = QO.UINV;
	min_log = QO.min_log;

	bigfloat::set_precision(oprec);

	// set last-used quadratic_order
	qo_l().set_last(this);
}



//
// operator =
//
// Task:
//    set to a copy of QO.
//

quadratic_order & quadratic_order::operator = (const quadratic_order & QO)
{
	debug_handler("quadratic_order", "operator = ");

	this->assign(QO);
	return *this;
}



//
// quadratic_order::discriminant()
//
// Task:
//      returns the discriminant
//

bigint
quadratic_order::discriminant() const
{
	return Delta;
}



//
// quadratic_order::nu_bound()
//
// Task:
//      returns the bound used for NUCOMP and NUDUPL
//

bigint
quadratic_order::nu_bound() const
{
	return nucomp_bound;
}



//
// quadratic_order::is_zero()
//
// Task:
//      tests if the quadratic order is the empty order
//

bool
quadratic_order::is_zero() const
{
	debug_handler("quadratic_order", "is_zero");

	return Delta.is_zero();
}



//
// quadratic_order::is_equal()
//
// Task:
//      tests if the quadratic orders are equal
//

bool
quadratic_order::is_equal(const quadratic_order & QO) const
{
	debug_handler("quadratic_order", "is_equal");

	return !Delta.compare(QO.Delta);
}



//
// quadratic_order::is_subset()
//
// Task:
//      tests if the quadratic order is a subset of QO2.
//

bool
quadratic_order::is_subset(quadratic_order & QO2)
{
	debug_handler("quadratic_order", "is_subset");

	return ((this->is_proper_subset(QO2)) || (this->is_equal(QO2)));
}



//
// quadratic_order::is_proper_subset()
//
// Task:
//      tests if the quadratic order is a proper subset of QO2.
//

bool
quadratic_order::is_proper_subset(quadratic_order & QO2)
{
	debug_handler("quadratic_order", "is_proper_subset");

	bigint c1, c2, m1, m2, r;

	// set last-used quadratic_order
	qo_l().set_last(this);

	c1.assign(conductor());
	c2.assign(QO2.conductor());
	divide(m1, Delta, c1*c1);
	divide(m2, QO2.Delta, c2*c2);
	remainder(r, c1, c2);

	return ((m1 == m2) && (r.is_zero()) && (c1 != c2));
}



//
// operator ==
//
// Task:
//      tests if QO1 and QO2 are equal (same discriminant)
//

bool operator == (const quadratic_order & QO1, const quadratic_order & QO2)
{
	debug_handler("quadratic_order", "operator == ");

	return !QO1.Delta.compare(QO2.Delta);
}



//
// operator !=
//
// Task:
//      tests if QO1 and QO2 are not equal
//

bool operator != (const quadratic_order & QO1, const quadratic_order & QO2)
{
	debug_handler("quadratic_order", "operator != ");

	return QO1.Delta.compare(QO2.Delta);
}



//
// operator <=
//
// Task:
//      tests if QO1 is a subset of QO2
//

bool operator <= (quadratic_order & QO1, quadratic_order & QO2)
{
	debug_handler("quadratic_order", "operator <= ");

	return (QO1.is_subset(QO2));
}



//
// operator <
//
// Task:
//      tests if QO1 is a proper subset of QO2
//

bool operator < (quadratic_order & QO1, quadratic_order & QO2)
{
	debug_handler("quadratic_order", "operator < ");

	return (QO1.is_proper_subset(QO2));
}



//
// operator >=
//
// Task:
//      tests if QO2 is a subset of QO1
//

bool operator >= (quadratic_order & QO1, quadratic_order & QO2)
{
	debug_handler("quadratic_order", "operator >= ");

	return (QO2.is_subset(QO1));
}



//
// operator >
//
// Task:
//      tests if QO2 is a proper subset of QO1
//

bool operator > (quadratic_order & QO1, quadratic_order & QO2)
{
	debug_handler("quadratic_order", "operator >");

	return (QO2.is_proper_subset(QO1));
}



//
// operator !
//
// Task:
//      tests if QO is the empty order
//

bool operator ! (const quadratic_order & QO)
{
	debug_handler("quadratic_order", "operator !");

	return (QO.is_zero());
}



//
// is_quadratic_discriminant(long)
//
// Task:
//     returns true if D is a quadratic discriminant (congruent to 0 or 1 mod
//     4 and not a perfect square).
//

bool
is_quadratic_discriminant(const long D)
{
	debug_handler("quadratic_order", "is_quadratic_discriminant(long)");

	return is_quadratic_discriminant(bigint(D));
}



//
// is_quadratic_discriminant(bigint)
//
// Task:
//     returns true if D is a quadratic discriminant (congruent to 0 or 1 mod
//     4 and not a perfect square).
//

bool
is_quadratic_discriminant(const bigint & D)
{
	debug_handler("quadratic_order", "is_quadratic_discriminant(bigint)");

	long k;
	bool is_disc;

	is_disc = true;

	// testing whether D is = 0,1 (mod 4)
	k = remainder(D, 4);
	if (k < 0)
		k += 4;
	if (k > 1)
		is_disc = false;
	else if (D.is_gt_zero()) {
		// testing whether D is a perfect square
		if (D.is_square())
			is_disc = false;
	}

	return is_disc;
}



//
// quadratic_order::is_imaginary()
//
// Task:
//      returns true if the quadratic order is imaginary (discriminant is less
//      than 0)
//

bool
quadratic_order::is_imaginary() const
{
	debug_handler("quadratic_order", "is_imaginary");

	return Delta.is_lt_zero();
}



//
// quadratic_order::is_real()
//
// Task:
//      returns true if the quadratic order is real (discriminant is greater
//      than 0)
//

bool
quadratic_order::is_real() const
{
	debug_handler("quadratic_order", "is_real");

	return Delta.is_gt_zero();
}



//
// quadratic_order::is_R_computed()
//
// Task:
//      returns true if the regulator has already been computed.
//

bool
quadratic_order::is_R_computed() const
{
	debug_handler("quadratic_order", "is_R_computed");

	long oprec = bigfloat::get_precision();
	bool iscomp;

	bigfloat::set_precision(prec);
	iscomp = (R.is_gt_zero()) || (Delta.is_zero());
	bigfloat::set_precision(oprec);

	return iscomp;
}



//
// quadratic_order::is_h_computed()
//
// Task:
//      returns true if the class number has already been computed.
//

bool
quadratic_order::is_h_computed() const
{
	debug_handler("quadratic_order", "is_h_computed");

	return (h.is_gt_zero()) || (Delta.is_zero());
}



//
// quadratic_order::is_L_computed()
//
// Task:
//      returns true if the L(1) function has already been computed.
//

bool
quadratic_order::is_L_computed() const
{
	debug_handler("quadratic_order", "is_L_computed");

	long oprec = bigfloat::get_precision();
	bool iscomp;

	bigfloat::set_precision(prec);
	iscomp = (L.is_gt_zero()) || (Delta.is_zero());
	bigfloat::set_precision(oprec);

	return iscomp;
}



//
// quadratic_order::is_C_computed()
//
// Task:
//      returns true if the C(Delta) function has already been computed.
//

bool
quadratic_order::is_C_computed() const
{
	debug_handler("quadratic_order", "is_C_computed");

	long oprec = bigfloat::get_precision();
	bool iscomp;

	bigfloat::set_precision(prec);
	iscomp = (Cfunc.is_gt_zero()) || (Delta.is_zero());
	bigfloat::set_precision(oprec);

	return iscomp;
}



//
// quadratic_order::is_CL_computed()
//
// Task:
//      returns true if the class group has already been computed.
//

bool
quadratic_order::is_CL_computed() const
{
	debug_handler("quadratic_order", "is_CL_computed");

	return (CL.size() > 0) || (Delta.is_zero());
}



//
// quadratic_order::is_subexp_computed()
//
// Task:
//      returns true if the class group has already been computed, and it was
//      done with the subexponential algorithm.
//

bool
quadratic_order::is_subexp_computed() const
{
	debug_handler("quadratic_order", "is_subexp_computed");

	return (subused);
}



//<MM>
//
// quadratic_order::is_fundamental_unit_computed()
//
// Task:
//      returns true if the fundamental unit has already been computed.
//

bool quadratic_order
::is_fundamental_unit_computed() const
{
	debug_handler("quadratic_order", "is_fundamental_unit_computed");

	return fund_unit.is_initialized();
}



//<MM>




//
// quadratic_order::bach_bound()
//
// Task:
//      returns Bach's bound on the norms of the prime ideals which generate
//      the class group.
//

int
quadratic_order::bach_bound()
{
	debug_handler("quadratic_order", "bach_bound()");

	return bach_bound(conductor());
}



//
// quadratic_order::bach_bound(const bigint & f)
//
// Task:
//      returns Bach's bound on the norms of the prime ideals which generate
//      the class group, using Table 3 of Bach90.
//

int
quadratic_order::bach_bound(const bigint & f)
{
	debug_handler("quadratic_order", "bach_bound(bigint)");

	bigfloat temp, C1, C2;
	bigint D;
	int BachBound, dec;
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);

	if (f.is_one()) {
		// order is maximal
		dec = decimal_length(Delta*Delta);

		if (dec >= 1000) {
			C1.assign(1.011);
			C2.assign(14.236);
		}
		else if (dec >= 100) {
			C1.assign(1.044);
			C2.assign(2.336);
		}
		else if (dec >= 50) {
			C1.assign(1.065);
			C2.assign(0.439);
		}
		else if (dec >= 20) {
			C1.assign(1.110);
			C2.assign(-1.367);
		}
		else if (dec >= 10) {
			C1.assign(1.166);
			C2.assign(-2.350);
		}
		else if (dec >= 5) {
			C1.assign(1.259);
			C2.assign(-3.077);
		}
		else {
			C1.assign(1.571);
			C2.assign(-3.570);
		}

		temp = C1 * log(bigfloat(Delta*Delta)) + C2;
		square(temp, temp);
		floor(temp, temp);
		temp.intify(BachBound);
	}
	else {
		D = Delta*Delta / (f*f);
		dec = decimal_length(D);

		if (dec >= 1000) {
			C1.assign(1.076);
			C2.assign(4.206);
		}
		else if (dec >= 100) {
			C1.assign(1.287);
			C2.assign(-0.887);
		}
		else if (dec >= 50) {
			C1.assign(1.420);
			C2.assign(-1.701);
		}
		else if (dec >= 20) {
			C1.assign(1.682);
			C2.assign(-2.471);
		}
		else if (dec >= 10) {
			C1.assign(1.974);
			C2.assign(-2.885);
		}
		else if (dec >= 5) {
			C1.assign(2.379);
			C2.assign(-3.189);
		}
		else {
			C1.assign(3.166);
			C2.assign(-3.468);
		}

		temp = C1 * log(bigfloat(D)) + C2;
		square(temp, temp);
		floor(temp, temp);
		temp.intify(BachBound);
	}

	bigfloat::set_precision(oprec);

	return BachBound;
}



//
// swap
//
// Task:
//      swaps Q1 and Q2
//

void
swap(quadratic_order & Q1, quadratic_order & Q2)
{
	debug_handler("quadratic_order", "swap");

	quadratic_order C;

	C.assign(Q1);
	Q1.assign(Q2);
	Q2.assign(C);
}



//
// quadratic_order::conductor()
//
// Task:
//      computes the conductor of the order, i.e., the integer f such that
//      Delta = f^2 Delta_0 for some fundamental discriminant Delta_0.
//      The discriminant is factored if it has not been already.
//

bigint
quadratic_order::conductor()
{
	debug_handler("quadratic_order", "conductor");

	bigint cond, b, temp;
	int e;
	long m4;
	lidia_size_t num_facts, i;

	// set last-used quadratic_order
	qo_l().set_last(this);

	cond.assign_one();
	factor_discriminant();
	num_facts = disc_fact.no_of_comp();
	for (i = 0; i < num_facts; ++i) {
		b = disc_fact.base(i);
		e = disc_fact.exponent(i);
		if ((b > 0) && (e > 1)) {
			if (b == 2) {
				if ((e % 2) == 1)
					shift_right(temp, Delta, e-1);
				else
					shift_right(temp, Delta, e);
				remainder(m4, temp, 4);
				if (m4 < 0)
					m4 += 4;
				if (m4 == 1)
					e >>= 1;
				else
					e = (e-2) >> 1;
			}
			else
				e >>= 1;

			if (e > 0) {
				power(temp, b, e);
				multiply(cond, cond, temp);
			}
		}
	}

	return cond;
}



//
// quadratic_order::is_maximal()
//
// Task:
//      returns true if the quadratic order is maximal (conductor = 1).
//

bool
quadratic_order::is_maximal()
{
	debug_handler("quadratic_order", "is_maximal");

	return (conductor() == bigint(1));
}



//
// quadratic_order::maximize()
//
// Task:
//      maximizes the quadratic order
//

void
quadratic_order::maximize()
{
	debug_handler("quadratic_order", "maximize()");

	bigint cond;

	// set last-used quadratic_order
	qo_l().set_last(this);

	square(cond, conductor());
	assign(Delta / cond);
}



//
// quadratic_order::maximize(quadratic_order)
//
// Task:
//      returns the maximal quadratic order of which this quadratic order is
//      a subset.
//

void
quadratic_order::maximize(quadratic_order & max_ord)
{
	debug_handler("quadratic_order", "maximize()");

	bigint cond;

	// set last-used quadratic_order
	qo_l().set_last(this);

	square(cond, conductor());
	max_ord.assign(Delta / cond);
}



//
// kronecker
//
// Task:
//      computes the kronecker symbol (D/p)
//

int
kronecker(const bigint & D, const bigint & p)
{
	debug_handler("quadratic_order", "kronecker");

	bigint c1;
	int kro;

	// if p > 2, use jacobi symbol
	if (p.compare(bigint(2))) {
		remainder(c1, D, p);
		if (c1.is_lt_zero())
			add(c1, c1, p);
		kro = jacobi(c1, p);
	}
	else {
		remainder(c1, D, 8);
		if (c1.is_lt_zero())
			add(c1, c1, 8);
		if (c1.is_one())
			kro = 1;
		else if (!c1.compare(bigint(5)))
			kro = -1;
		else
			kro = 0;
	}

	return kro;
}



//
// quadratic_order::generate_optimal_Q()
//
// Task:
//      computes the optimal value of Q for computing h* such that
//      h* < h < 2h*
//

long
quadratic_order::generate_optimal_Q()
{
	debug_handler("quadratic_order", "generate_optimal_Q");

	long OQ;
	bigfloat A, l2;
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if ((PL.get_lower_bound() != 2) || (PL.get_upper_bound() < 1000000))
		PL.resize(2, 1000000);

	l2 = log(sqrt(bigfloat(2.0)));
	OQ = PL.get_first_prime(3);
	A = estimate_L1_error(OQ);
	while (A >= l2) {
		OQ = PL.get_next_prime();
		A = estimate_L1_error(OQ);
	}

	bigfloat::set_precision(oprec);

	return OQ;
}



//
// quadratic_order::get_optimal_Q()
//
// Task:
//      returns a value of Q from a pre-computed table which will compute h*
//      such that h* < h < 2h*.
//

long
quadratic_order::get_optimal_Q()
{
	debug_handler("quadratic_order", "get_optimal_Q");

	long Dlog;
	bigfloat temp;

	// set last-used quadratic_order
	qo_l().set_last(this);

	temp = floor(log(bigfloat(abs(Delta))) / std::log(10.0));
	temp.longify(Dlog);

	return OQvals[Dlog / 5];
}



//
// quadratic_order::generate_optimal_Q_cnum()
//
// Task:
//      computes the optimal value of Q required to determine whether the
//      approximation of h gleaned from the analytic class number formula is
//      actually the class number, provided that h <=3.
//

long
quadratic_order::generate_optimal_Q_cnum(const bigfloat & h2,
                                         const bigfloat & t)
{
	debug_handler("quadratic_order", "generate_optimal_Q_cnum");

	long OQ;
	bigfloat A, val;
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if ((PL.get_lower_bound() != 2) || (PL.get_upper_bound() < 1000000))
		PL.resize(2, 1000000);

	val = log((h2 + bigfloat(1.0)) / (h2 + t));
	OQ = PL.get_first_prime(3);
	A = estimate_L1_error(OQ);
	while (A >= val) {
		OQ = PL.get_next_prime();
		A = estimate_L1_error(OQ);
	}

	bigfloat::set_precision(oprec);

	return OQ;
}



//
// quadratic_order::get_optimal_Q_cnum()
//
// Task:
//      returns a value of Q from a pre-computed table which is sufficient
//      to determine whether the approximation of h gleaned from the analytic
//      class number formula is actually the class number, provided that h <=3.
//

long
quadratic_order::get_optimal_Q_cnum()
{
	debug_handler("quadratic_order", "get_optimal_Q_cnum");

	long Dlog;
	bigfloat temp;

	// set last-used quadratic_order
	qo_l().set_last(this);

	temp = floor(log(bigfloat(abs(Delta))) / std::log(10.0));
	temp.longify(Dlog);

	return OQvals_cnum[Dlog / 5];
}



//
// quadratic_order::generate_optimal_Q_cfunc()
//
// Task:
//      computes the optimal value of Q required to compute C(Delta) to 8
//      significant digits.
//

long
quadratic_order::generate_optimal_Q_cfunc()
{
	debug_handler("quadratic_order", "generate_optimal_Q_cfunc");

	long OQ;
	bigfloat B, sb, val, lDelta, lQ, temp1, temp2;
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);

	// set last-used quadratic_order
	qo_l().set_last(this);

	sqrt(val, bigfloat(1.0000002));
	val.divide_by_2();
	inc(val);
	val.assign(log(val));

	lDelta.assign(log(bigfloat(abs(Delta))));

	OQ = 100;
	lQ.assign(log(bigfloat(OQ)));

	multiply(temp1, Pi(), lQ);
	temp1.invert();
	square(temp2, lQ);
	divide(temp2, bigfloat(5.3), temp2);
	add(temp1, temp1, temp2);
	multiply(B, lDelta, temp1);

	divide(temp1, bigfloat(4.0), lQ);
	divide(temp2, bigfloat(1.0), Pi());
	add(B, B, temp1);
	add(B, B, temp2);

	power(temp1, bigfloat(OQ), bigfloat(1.5));
	multiply(temp1, temp1, 9);
	multiply(temp2, lQ, 13);
	add(temp2, temp2, 8);
	divide(sb, temp2, temp1);
	multiply(sb, sb, B);

	square(temp1, bigfloat(OQ));
	divide(temp1, bigfloat(2.0), temp1);
	add(sb, sb, temp1);

	while (sb >= val) {
		if (OQ == 100)
			OQ = 1000;
		else if (OQ < 5000)
			OQ += 1000;
		else
			OQ += 5000;

		lQ.assign(log(bigfloat(OQ)));

		multiply(temp1, Pi(), lQ);
		temp1.invert();
		square(temp2, lQ);
		divide(temp2, bigfloat(5.3), temp2);
		add(temp1, temp1, temp2);
		multiply(B, lDelta, temp1);

		divide(temp1, bigfloat(4.0), lQ);
		divide(temp2, bigfloat(1.0), Pi());
		add(B, B, temp1);
		add(B, B, temp2);

		power(temp1, bigfloat(OQ), bigfloat(1.5));
		multiply(temp1, temp1, 9);
		multiply(temp2, lQ, 13);
		add(temp2, temp2, 8);
		divide(sb, temp2, temp1);
		multiply(sb, sb, B);

		square(temp1, bigfloat(OQ));
		divide(temp1, bigfloat(2.0), temp1);
		add(sb, sb, temp1);
	}

	bigfloat::set_precision(oprec);

	return OQ;
}



//
// quadratic_order::get_optimal_Q_cfunc()
//
// Task:
//      returns a value of Q from a pre-computed list which is sufficient to
//      compute an approximation of C(Delta) accurate to 8 significant digits.
//

long
quadratic_order::get_optimal_Q_cfunc()
{
	debug_handler("quadratic_order", "get_optimal_Q_cfunc");

	long Dlog;
	bigfloat temp;

	// set last-used quadratic_order
	qo_l().set_last(this);

	temp = floor(log(bigfloat(abs(Delta))) / std::log(10.0));
	temp.longify(Dlog);

	return OQvals_cfunc[Dlog / 5];
}



//
// quadratic_order::estimate_C
//
// Task:
//      computes a truncated product approximation of C(D) using primes
//      less than Q.
//

bigfloat
quadratic_order::estimate_C(const long nQ)
{
	debug_handler("quadratic_order", "estimate_C");

	int kron;
	long P, QQ;
	double E;

	// set last-used quadratic_order
	qo_l().set_last(this);

	E = 1.0;

	// generate sufficiently many primes
	QQ = 1000000;
	if (QQ > nQ)
		QQ = nQ+10;
	P = 3;

	if ((QQ > static_cast<long>(PL.get_upper_bound())) || (PL.get_lower_bound() != 2))
		PL.resize(2, QQ);
	P = PL.get_first_prime(P);

	// compute partial product for LQ <= p < QQ
	while ((P< nQ) && (P > 0)) {
		kron = kronecker(Delta, bigint(P));
		E *= 1.0L - (static_cast<double>(kron) / (P - 1.0));

		P = PL.get_next_prime();
		if (P == 0) {
			QQ = PL.get_upper_bound() + 1000000;
			if (QQ > nQ)
				QQ = nQ+10;
			PL.resize(PL.get_upper_bound() + 1, QQ);
			P = PL.get_first_prime();
		}
	}

	return bigfloat(E);
}



//
// quadratic_order::estimate_L
//
// Task:
//      computes a truncated product approximation of L(s,X) using primes
//      less than Q.
//

bigfloat
quadratic_order::estimate_L(const int s, const long nQ)
{
	debug_handler("quadratic_order", "estimate_L");

	int kron, i;
	long m8, P, s2, QQ;
	double E, Ps;
	bigfloat estL;

	// set last-used quadratic_order
	qo_l().set_last(this);

	// compute partial product for p = 2
	m8 = remainder(Delta, 8);
	if (m8 < 0)
		m8 += 8;
	if (m8 == 1)
		kron = 1;
	else if (m8 == 5)
		kron = -1;
	else
		kron = 0;
	s2 = 1;
	s2 <<= s;
	E = std::log(static_cast<double>(s2) / static_cast<double>(s2-kron));

	// generate sufficiently many primes
	QQ = 1000000;
	if (QQ > nQ)
		QQ = nQ+10;
	P = 3;

	if ((QQ > static_cast<long>(PL.get_upper_bound())) || (PL.get_lower_bound() != 2))
		PL.resize(2, QQ);
	P = PL.get_first_prime(P);

	// compute partial product for LQ <= p < QQ
	while ((P< nQ) && (P > 0)) {
		kron = kronecker(Delta, bigint(P));
		Ps = P;
		for (i = 1; i < s; ++i)
			Ps *= P;
		E += std::log(Ps / (Ps - kron));

		P = PL.get_next_prime();
		if (P == 0) {
			QQ = PL.get_upper_bound() + 1000000;
			if (QQ > nQ)
				QQ = nQ+10;
			PL.resize(PL.get_upper_bound() + 1, QQ);
			P = PL.get_first_prime();
		}
	}

	return std::exp(E);
}



//
// quadratic_order::estimate_L1
//
// Task:
//      computes a truncated product approximation of L(1,X) using primes
//      less than 2Q.  This function uses Bach's method of a weighted average
//      of truncated Euler products.
//

bigfloat
quadratic_order::estimate_L1(const long nQ)
{
	debug_handler("quadratic_order", "estimate_L1");

	int kron;
	long Q2, m8, P, QQ;
	double E, C, wt, s;
	register long i;
	bigfloat nFI;

	// set last-used quadratic_order
	qo_l().set_last(this);

	// if nQ = Q, then the estimate has already been computed
	if ((nQ == Q) && (!FI.is_zero())) {
		nFI = FI;
		return(nFI);
	}

	// compute partial product for p = 2
	Q = nQ;
	m8 = remainder(Delta, 8);
	if (m8 < 0)
		m8 += 8;
	if (m8 == 1)
		E = std::log(2.0);
	else if (m8 == 5)
		E = std::log(2.0 / 3.0);
	else
		E = 0.0;

	// compute weights
	Q2 = Q << 1;
	C = 0.0;
	for (i = Q; i <= Q2-1; ++i)
		C += static_cast<double>(i) * std::log(static_cast<double>(i));

	// generate list of primes
	QQ = 1000000;
	if (QQ > Q2)
		QQ = Q2+10;
	P = 3;

	if ((QQ > static_cast<long>(PL.get_upper_bound())) || (PL.get_lower_bound() != 2))
		PL.resize(2, QQ);
	P = PL.get_first_prime(P);

	// compute partial product for p < Q
	while ((P< Q) && (P > 0)) {
		kron = kronecker(Delta, bigint(P));
		E += std::log(static_cast<double>(P) / (P - kron));

		P = PL.get_next_prime();
		if (P == 0) {
			QQ = PL.get_upper_bound() + 1000000;
			if (QQ > Q2)
				QQ = Q2+10;
			PL.resize(PL.get_upper_bound() + 1, QQ);
			P = PL.get_first_prime();
		}
	}

	// computed weighted partial products for Q < p < 2Q
	wt = 1.0;
	s = 0.0;
	for (i = Q; i <= P; ++i)
		s += static_cast<double>(i) * std::log(static_cast<double>(i));
	wt -= s / C;

	while ((P< Q2) && (P > 0)) {
		kron = kronecker(Delta, bigint(P));
		E += wt*std::log(static_cast<double>(P) / (P - kron));

		i = P + 1;
		P = PL.get_next_prime();
		if (P == 0) {
			QQ = PL.get_upper_bound() + 1000000;
			if (QQ > Q2)
				QQ = Q2+10;
			PL.resize(PL.get_upper_bound() + 1, QQ);
			P = PL.get_first_prime();
		}
		s = 0.0;
		for (; i <= P; ++i)
			s += static_cast<double>(i) * std::log(static_cast<double>(i));
		wt -= s / C;
	}

	nFI = exp(bigfloat(E));

	// set static value
	FI = nFI;

	return(nFI);
}



//
// quadratic_order::estimate_L1_error
//
// Task:
//      computes an upper bound of the error in the approximation of L(1,X)
//      computed by the function estimate_L1
//

bigfloat
quadratic_order::estimate_L1_error(const long nQ) const
{
	debug_handler("quadratic_order", "estimate_L1_error");

	bigfloat A, lq;

	lq.assign(log(bigfloat(nQ)));
	if (Q >= 100000)
		A = (6.338*log(bigfloat(abs(Delta))) + 17.031) / (lq * sqrt(bigfloat(nQ)));
	else if (Q >= 50000)
		A = (6.378*log(bigfloat(abs(Delta))) + 17.397) / (lq * sqrt(bigfloat(nQ)));
	else if (Q >= 10000)
		A = (6.510*log(bigfloat(abs(Delta))) + 18.606) / (lq * sqrt(bigfloat(nQ)));
	else if (Q >= 5000)
		A = (6.593*log(bigfloat(abs(Delta))) + 19.321) / (lq * sqrt(bigfloat(nQ)));
	else if (Q >= 1000)
		A = (6.897*log(bigfloat(abs(Delta))) + 21.528) / (lq * sqrt(bigfloat(nQ)));
	else
		A = (7.106*log(bigfloat(abs(Delta))) + 22.845) / (lq * sqrt(bigfloat(nQ)));

	return A;
}



//
// quadratic_order::Lfunction
//
// Task:
//      returns an approximation of L(1,X) computed via the analytic class
//      number formula.
//

bigfloat
quadratic_order::Lfunction()
{
	debug_handler("quadratic_order", "Lfunction");

	long oprec = bigfloat::get_precision();
	bigfloat outL;

	bigfloat::set_precision(prec);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_L_computed()) {
		// L has not been computed - compute it

		// compute class number, if needed
		if (h.is_zero())
			class_number();

		if (is_imaginary()) {
			// L = pi*h / sqrt(Delta)
			multiply(L, Pi(), bigfloat(h));
			divide(L, L, sqrt(bigfloat(-Delta)));
			if (Delta == -4)  L.divide_by_2();
			if (Delta == -3)  divide(L, L, bigfloat(3.0));
		}
		else {
			// L = 2hR / sqrt(Delta)
			multiply(L, R, bigfloat(h));
			L.multiply_by_2();
			divide(L, L, sqrt(bigfloat(Delta)));
		}
	}

	outL.assign(L);
	bigfloat::set_precision(oprec);

	return outL;
}



//
// quadratic_order::LDfunction
//
// Task:
//      returns an approximation of the LD(1) function defined by Shanks
//      (L(1,X) without the 2-factor) computed via the analytic class
//      number formula.
//

bigfloat
quadratic_order::LDfunction()
{
	debug_handler("quadratic_order", "LDfunction");

	bigfloat LD;
	long m8;
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);

	// set last-used quadratic_order
	qo_l().set_last(this);

	m8 = remainder(Delta, 8);
	if (m8 < 0)
		m8 += 8;
	LD.assign(Lfunction());
	if (m8 == 1)
		LD.divide_by_2();
	else if (m8 == 5)
		multiply(LD, LD, bigfloat(1.5));

	bigfloat::set_precision(oprec);

	return LD;
}



//
// quadratic_order::LLI
//
// Task:
//      returns an approximation of the Lower Littlewood Index (defined by
//      Shanks).
//

bigfloat
quadratic_order::LLI()
{
	debug_handler("quadratic_order", "LLI");

	bigfloat LLI, c;
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (Delta.is_zero())
		LLI.assign_zero();
	else {
		c.assign(exp(Euler()) / (Pi() * Pi()));

		if (Delta.is_odd())
			c *= bigfloat(12.0);
		else
			c *= bigfloat(8.0);

		LLI = Lfunction() * c * log(log(bigfloat(abs(Delta))));
	}

	bigfloat::set_precision(oprec);

	return LLI;
}



//
// quadratic_order::LLI_D
//
// Task:
//      returns an approximation of the LLI_D value (defined by Shanks).
//

bigfloat
quadratic_order::LLI_D()
{
	debug_handler("quadratic_order", "LLI_D");

	bigfloat LLI, c;
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (Delta.is_zero())
		LLI.assign_zero();
	else {
		c.assign(bigfloat(8.0) * exp(Euler()) / (Pi() * Pi()));
		LLI = LDfunction() * c * log(log(bigfloat(abs(Delta << 2))));
	}

	bigfloat::set_precision(oprec);

	return LLI;
}



//
// quadratic_order::ULI
//
// Task:
//      returns an approximation of the Upper Littlewood Index (defined by
//      Shanks).
//

bigfloat
quadratic_order::ULI()
{
	debug_handler("quadratic_order", "ULI");

	bigfloat ULI, c;
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (Delta.is_zero())
		ULI.assign_zero();
	else {
		c.assign(exp(Euler()));

		if (Delta.is_odd())
			c.multiply_by_2();

		ULI = Lfunction() / (c * log(log(bigfloat(abs(Delta)))));
	}

	bigfloat::set_precision(oprec);

	return ULI;
}



//
// quadratic_order::ULI_D
//
// Task:
//      returns an approximation of the ULI_D value (defined by Shanks).
//

bigfloat
quadratic_order::ULI_D()
{
	debug_handler("quadratic_order", "ULI_D");

	bigfloat ULI, c;
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (Delta.is_zero())
		ULI.assign_zero();
	else {
		c.assign(exp(Euler()));
		ULI = LDfunction() / (c * log(log(bigfloat(abs(Delta << 2)))));
	}

	bigfloat::set_precision(oprec);

	return ULI;
}



//
// quadratic_order::Cfunction
//
// Task:
//      returns an approximation of C(Delta) accurate to 8 significant digits.
//

bigfloat
quadratic_order::Cfunction()
{
	debug_handler("quadratic_order", "Cfunction");

	bigfloat temp, E4, outC;
	bigint PP;
	int kron;
	long Q2, m8, P, num, QQ;
	double E, E2, E3, C, wt, dP, dP2, s;
	register long i;
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_C_computed()) {
		Q2 = get_optimal_Q_cfunc();

		if (h.is_zero()) {
			Q = Q2 >> 1;

			// compute partial product for p = 2
			m8 = remainder(Delta, 8);
			if (m8 < 0)
				m8 += 8;
			if (m8 == 1) {
				E = std::log(2.0);
				E2 = std::log(4.0 / 3.0);
				Cfunc.assign(2.5);
			}
			else if (m8 == 5) {
				E = std::log(2.0 / 3.0);
				E2 = std::log(4.0 / 5.0);
				Cfunc.assign(0.5);
			}
			else {
				E = 0.0;
				E2 = 0.0;
				Cfunc.assign(15.0/16.0);
			}
			E3 = 0.0;

			// compute weights
			C = 0.0;
			for (i = Q; i <= Q2-1; ++i)
				C += static_cast<double>(i) * std::log(static_cast<double>(i));

			// generate list of primes
			QQ = 1000000;
			if (QQ > Q2)
				QQ = Q2+10;
			P = 3;

			if ((QQ > static_cast<long>(PL.get_upper_bound())) || (PL.get_lower_bound() != 2))
				PL.resize(2, QQ);
			P = PL.get_first_prime(P);

			// compute partial product for p < Q
			while ((P< Q) && (P > 0)) {
				kron = kronecker(Delta, bigint(P));
				E += std::log(static_cast<double>(P) / (P - kron));

				dP = static_cast<double>(P);
				dP2 = dP*dP;
				E2 += std::log(dP2 / (dP2 - kron));

				if (kron == 1) {
					dP2 = (dP - 1.0) * (dP - 1.0);
					dP2 *= dP;
					E3 += std::log(1.0 - (2.0 / dP2));
				}

				P = PL.get_next_prime();
				if (P == 0) {
					QQ = PL.get_upper_bound() + 1000000;
					if (QQ > Q2)
						QQ = Q2+10;
					PL.resize(PL.get_upper_bound() + 1, QQ);
					P = PL.get_first_prime();
				}
			}

			// computed weighted partial products for Q < p < 2Q
			wt = 1.0;
			s = 0.0;
			for (i = Q; i <= P; ++i)
				s += static_cast<double>(i) * std::log(static_cast<double>(i));
			wt -= s / C;
			while ((P< Q2) && (P > 0)) {
				kron = kronecker(Delta, bigint(P));
				E += wt * std::log(static_cast<double>(P) / (P - kron));

				dP = static_cast<double>(P);
				dP2 = dP*dP;
				E2 += std::log(dP2 / (dP2 - kron));

				if (kron == 1) {
					dP2 = (dP - 1.0) * (dP - 1.0);
					dP2 *= dP;
					E3 += std::log(1.0 - (2.0 / dP2));
				}

				i = P + 1;
				P = PL.get_next_prime();
				if (P == 0) {
					QQ = PL.get_upper_bound() + 1000000;
					if (QQ > Q2)
						QQ = Q2+10;
					PL.resize(PL.get_upper_bound() + 1, QQ);
					P = PL.get_first_prime();
				}
				s = 0.0;
				for (; i <= P; ++i)
					s += static_cast<double>(i) * std::log(static_cast<double>(i));
				wt -= s / C;
			}

			FI = exp(bigfloat(E));
			class_number();
		}
		else {
			// compute partial product for p = 2
			m8 = remainder(Delta, 8);
			if (m8 < 0)
				m8 += 8;
			if (m8 == 1) {
				E2 = std::log(4.0 / 3.0);
				Cfunc.assign(2.5);
			}
			else if (m8 == 5) {
				E2 = std::log(4.0 / 5.0);
				Cfunc.assign(0.5);
			}
			else {
				E2 = 0.0;
				Cfunc.assign(15.0/16.0);
			}
			E3 = 0.0;

			// generate list of primes
			QQ = 1000000;
			if (QQ > Q2)
				QQ = Q2+10;
			P = 3;

			if ((QQ > static_cast<long>(PL.get_upper_bound())) || (PL.get_lower_bound() != 2))
				PL.resize(2, QQ);
			P = PL.get_first_prime(P);

			// compute partial product for p < Q
			while ((P< Q2) && (P > 0)) {
				kron = kronecker(Delta, bigint(P));

				dP = static_cast<double>(P);
				dP2 = dP*dP;
				E2 += std::log(dP2 / (dP2 - kron));

				if (kron == 1) {
					dP2 = (dP - 1.0) * (dP - 1.0);
					dP2 *= dP;
					E3 += std::log(1.0 - (2.0 / dP2));
				}

				P = PL.get_next_prime();
				if (P == 0) {
					QQ = PL.get_upper_bound() + 1000000;
					if (QQ > Q2)
						QQ = Q2+10;
					PL.resize(PL.get_upper_bound() + 1, QQ);
					P = PL.get_first_prime();
				}
			}
		}

		// compute C
		if (is_imaginary()) {
			multiply(Cfunc, Cfunc, sqrt(bigfloat(-Delta)));
			power(temp, Pi(), 3);
			divide(temp, temp, bigfloat(h));
			divide(temp, temp, 90);
			multiply(Cfunc, Cfunc, temp);

			// multiply by the number of roots of unity
			if (Delta == -4)  Cfunc.multiply_by_2();
			if (Delta == -3)  multiply(Cfunc, Cfunc, bigfloat(3.0));
		}
		else {
			multiply(Cfunc, Cfunc, sqrt(bigfloat(Delta)));
			power(temp, Pi(), 4);
			divide(temp, temp, bigfloat(h));
			divide(temp, temp, R);
			divide(temp, temp, 180);
			multiply(Cfunc, Cfunc, temp);
		}

		factor_discriminant();
		E4.assign_one();
		num = disc_fact.no_of_comp();
		for (i = 0; i < num; ++i) {
			PP = disc_fact.base(static_cast<lidia_size_t>(i));
			if (PP > 2) {
				power(temp, bigfloat(PP), 4);
				divide(temp, temp-bigfloat(1.0), temp);
				E4 *= temp;
			}
		}

		multiply(Cfunc, Cfunc, E4);
		divide(Cfunc, Cfunc, exp(bigfloat(E2)));
		multiply(Cfunc, Cfunc, exp(bigfloat(E3)));
	}

	outC.assign(Cfunc);
	bigfloat::set_precision(oprec);

	return outC;
}



//
// quadratic_order::Cfunction_bigfloat
//
// Task:
//      returns an approximation of C(Delta) accurate to 8 significant digits
//      using bigfloats
//

bigfloat_int
quadratic_order::Cfunction_bigfloat(int pmult)
{
	debug_handler("quadratic_order", "Cfunction_bigfloat");

	bigfloat_int temp, E4, outC, Cfunc_int;
	bigint PP;
	int kron;
	long Q2, m8, P, num, QQ;
	bigfloat_int E, E2, E3, C, wt, dP, dP2, s;
	register long i;
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec*pmult);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_C_computed()) {
		Q2 = get_optimal_Q_cfunc();

		if (h.is_zero()) {
			Q = Q2 >> 1;

			// compute partial product for p = 2
			m8 = remainder(Delta, 8);
			if (m8 < 0)
				m8 += 8;
			if (m8 == 1) {
				E = log(bigfloat_int(2.0));
				E2.assign(bigint(4), bigint(3));
				Cfunc_int.assign(2.5);
			}
			else if (m8 == 5) {
				E = log(bigfloat_int(2.0) / bigfloat_int(3.0));
				E2.assign(0.8);
				Cfunc_int.assign(0.5);
			}
			else {
				E.assign_zero();
				E2.assign(bigfloat_int(1));
				Cfunc_int.assign(0.9375);
			}
			E3.assign(bigfloat_int(1));

			// compute weights
			C.assign_zero();
			for (i = Q; i <= Q2-1; ++i)
				C += bigfloat_int(i) * log(bigfloat_int(i));

			// generate list of primes
			QQ = 1000000;
			if (QQ > Q2)
				QQ = Q2+10;
			P = 3;

			if ((QQ > static_cast<long>(PL.get_upper_bound())) || (PL.get_lower_bound() != 2))
				PL.resize(2, QQ);
			P = PL.get_first_prime(P);

			// compute partial product for p < Q
			while ((P< Q) && (P > 0)) {
				kron = kronecker(Delta, bigint(P));
				E += log(bigfloat_int(P) / bigfloat_int(P - kron));

				dP.assign(P);
				dP2 = dP*dP;
				multiply(E2, E2, (dP2 / (dP2 - kron)));

				if (kron == 1) {
					dP2 = (dP - bigfloat_int(1)) * (dP - bigfloat_int(1));
					dP2 *= dP;
					multiply(E3, E3, (bigfloat_int(1) - (bigfloat_int(2) / dP2)));
				}

				P = PL.get_next_prime();
				if (P == 0) {
					QQ = PL.get_upper_bound() + 1000000;
					if (QQ > Q2)
						QQ = Q2+10;
					PL.resize(PL.get_upper_bound() + 1, QQ);
					P = PL.get_first_prime();
				}
			}

			// computed weighted partial products for Q < p < 2Q
			wt.assign_one();
			s.assign_zero();
			for (i = Q; i <= P; ++i)
				s += bigfloat_int(i)*log(bigfloat_int(i));
			wt -= s / C;
			while ((P< Q2) && (P > 0)) {
				kron = kronecker(Delta, bigint(P));
				E += wt*log(bigfloat_int(P) / bigfloat_int(P - kron));

				dP.assign(P);
				dP2 = dP*dP;
				multiply(E2, E2, (dP2 / (dP2 - kron)));

				if (kron == 1) {
					dP2 = (dP - bigfloat_int(1)) * (dP - bigfloat_int(1));
					dP2 *= dP;
					multiply(E3, E3, (bigfloat_int(1) - (bigfloat_int(2) / dP2)));
				}

				i = P + 1;
				P = PL.get_next_prime();
				if (P == 0) {
					QQ = PL.get_upper_bound() + 1000000;
					if (QQ > Q2)
						QQ = Q2+10;
					PL.resize(PL.get_upper_bound() + 1, QQ);
					P = PL.get_first_prime();
				}
				s.assign_zero();
				for (; i <= P; ++i)
					s += bigfloat_int(i)*log(bigfloat_int(i));
				wt -= s / C;
			}

			exp(bigfloat_int(E)).approx(FI);
			class_number();
		}
		else {
			// compute partial product for p = 2
			m8 = remainder(Delta, 8);
			if (m8 < 0)
				m8 += 8;
			if (m8 == 1) {
				E2.assign(bigint(4), bigint(3));
				Cfunc_int.assign(2.5);
			}
			else if (m8 == 5) {
				E2.assign(bigfloat_int(0.8));
				Cfunc_int.assign(0.5);
			}
			else {
				E2.assign(bigfloat_int(1));
				Cfunc_int.assign(0.9375);
			}
			E3.assign(bigfloat_int(1));

			// generate list of primes
			QQ = 1000000;
			if (QQ > Q2)
				QQ = Q2+10;
			P = 3;

			if ((QQ > static_cast<long>(PL.get_upper_bound())) || (PL.get_lower_bound() != 2))
				PL.resize(2, QQ);
			P = PL.get_first_prime(P);

			// compute partial product for p < Q
			while ((P< Q2) && (P > 0)) {
				kron = kronecker(Delta, bigint(P));

				dP.assign(P);
				square(dP2, dP);
				multiply(E2, E2, (dP2 / (dP2 - bigfloat_int(kron))));

				if (kron == 1) {
					dP2 = (dP - bigfloat_int(1)) * (dP - bigfloat_int(1));
					dP2 *= dP;
					multiply(E3, E3, (bigfloat_int(1) - (bigfloat_int(2) / dP2)));
				}

				P = PL.get_next_prime();
				if (P == 0) {
					QQ = PL.get_upper_bound() + 1000000;
					if (QQ > Q2)
						QQ = Q2+10;
					PL.resize(PL.get_upper_bound() + 1, QQ);
					P = PL.get_first_prime();
				}
			}
		}

		// compute C
		if (is_imaginary()) {
			multiply(Cfunc_int, Cfunc_int, sqrt(bigfloat_int(-Delta)));
			temp = Pi();
			power(temp, temp, 3);
			divide(temp, temp, bigfloat_int(h));
			divide(temp, temp, 90);
			multiply(Cfunc_int, Cfunc_int, temp);

			// multiply by the number of roots of unity
			if (Delta == -4)  multiply(Cfunc_int, Cfunc_int, bigfloat_int(2));
			if (Delta == -3)  multiply(Cfunc_int, Cfunc_int, bigfloat_int(3));
		}
		else {
			multiply(Cfunc_int, Cfunc_int, sqrt(bigfloat_int(Delta)));
			temp = Pi();
			power(temp, temp, 4);
			divide(temp, temp, bigfloat_int(h));
			divide(temp, temp, R);
			divide(temp, temp, 180);
			multiply(Cfunc_int, Cfunc_int, temp);
		}

		factor_discriminant();
		E4.assign(bigfloat_int(1));
		num = disc_fact.no_of_comp();
		for (i = 0; i < num; ++i) {
			PP = disc_fact.base(static_cast<lidia_size_t>(i));
			if (PP > 2) {
				power(PP, PP, 4);
				temp.assign(PP-1, PP);
				E4 *= temp;
			}
		}
		divide(Cfunc_int, Cfunc_int, E2);
		multiply(Cfunc_int, Cfunc_int, E3);
		multiply(Cfunc_int, Cfunc_int, E4);
	}

	Cfunc_int.approx(Cfunc);
	outC.assign(Cfunc_int);
	bigfloat::set_precision(oprec);

	return outC;
}



//
// quadratic_order::LLI
//
// Task:
//      returns an approximation of the Lower Littlewood Index (defined by
//      Shanks).
//

bigfloat
quadratic_order::LLI(const bigint & nh)
{
	debug_handler("quadratic_order", "LLI(bigint)");

	bigfloat LL;
	bigint oh = h;

	h = nh;
	LL = LLI();
	h = oh;

	return LL;
}



//
// quadratic_order::LLI
//
// Task:
//      returns an approximation of the Lower Littlewood Index (defined by
//      Shanks).
//

bigfloat
quadratic_order::LLI(const bigint & nh, const bigfloat & nR)
{
	debug_handler("quadratic_order", "LLI(bigint, bigfloat)");

	bigfloat LL, oR = R;
	bigint oh = h;

	h = nh;
	R = nR;
	LL = LLI();
	h = oh;
	R = oR;

	return LL;
}



//
// quadratic_order::ULI
//
// Task:
//      returns an approximation of the Lower Littlewood Index (defined by
//      Shanks).
//

bigfloat
quadratic_order::ULI(const bigint & nh)
{
	debug_handler("quadratic_order", "ULI(bigint)");

	bigfloat UL;
	bigint oh = h;

	h = nh;
	UL = ULI();
	h = oh;

	return UL;
}



//
// quadratic_order::ULI
//
// Task:
//      returns an approximation of the Lower Littlewood Index (defined by
//      Shanks).
//

bigfloat
quadratic_order::ULI(const bigint & nh, const bigfloat & nR)
{
	debug_handler("quadratic_order", "ULI(bigint, bigfloat)");

	bigfloat UL, oR = R;
	bigint oh = h;

	h = nh;
	R = nR;
	UL = ULI();
	h = oh;
	R = oR;

	return UL;
}



//
// quadratic_order::Cfunction
//
// Task:
//      returns an approximation of C(Delta) accurate to 8 significant digits.
//

bigfloat
quadratic_order::Cfunction(const bigint & nh,
                           const rational_factorization & dfact,
                           int pmult)
{
	debug_handler("quadratic_order", "Cfunction(bigint)");

	bigfloat CC;
	bigfloat_int CC_int;
	bigint oh = h;
	rational_factorization ofact = disc_fact;

	h = nh;
	disc_fact = dfact;
	if (pmult) {
		CC_int = Cfunction_bigfloat(pmult);
		CC_int.approx(CC);
		std::cout << "-->C in " << CC_int << ", = " << CC << std::endl;
	}
	else
		CC = Cfunction();
	h = oh;
	disc_fact = ofact;

	return CC;
}



//
// quadratic_order::Cfunction
//
// Task:
//      returns an approximation of C(Delta) accurate to 8 significant digits.
//

bigfloat
quadratic_order::Cfunction(const bigint & nh,
                           const bigfloat & nR,
                           const rational_factorization & dfact,
                           int pmult)
{
	debug_handler("quadratic_order", "Cfunction(bigint, bigfloat)");

	bigfloat CC, oR = R;
	bigfloat_int CC_int;
	bigint oh = h;
	rational_factorization ofact = disc_fact;

	h = nh;
	R = nR;
	disc_fact = dfact;
	if (pmult) {
		CC_int = Cfunction_bigfloat(pmult);
		CC_int.approx(CC);
		std::cout << "-->C in " << CC_int << ", = " << CC << std::endl;
	}
	else
		CC = Cfunction();
	h = oh;
	R = oR;
	disc_fact = ofact;

	return CC;
}



//
// quadratic_order::regulator
//
// Task:
//      returns an approximation of the regulator.  The algorithm is chosen
//      based on the size of the discriminant.
//

bigfloat
quadratic_order::regulator()
{
	debug_handler("quadratic_order", "regulator");

	int num;
	bigfloat outR;
	quadratic_order *QO = qi_class::get_current_order_ptr();
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);
	qi_class::set_current_order(*this);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_R_computed()) {
		// R has not been computed yet - compute it

		if (is_imaginary())
			R.assign_one();

		else {
			num = decimal_length(Delta);
			if (num < RBJTB)
				regulator_BJT();
			else if (num < RSHANKSB)
				regulator_shanks();
			else
				class_group_siqs();
		}
	}

	outR.assign(R);
	bigfloat::set_precision(oprec);
	qi_class::set_current_order(*QO);

	return outR;
}



//
// quadratic_order::regulator(long k)
//
// Task:
//      returns an absolute k approximation to the regulator.
//      The algorithm is chosen based on the size of the discriminant.
//

xbigfloat quadratic_order
::regulator(long k)
{
	debug_handler("quadratic_order", "regulator(long k");

	xbigfloat outR;

	if (is_imaginary())
		outR.assign_one();

	else {
		lidia_error_handler("quadratic_order::regulator(long)",
				    "Not yet implemented.");

		if (R_xbig.is_zero() || R_xbig_prec < k+1) {
			this->fundamental_unit();
			R_xbig = fund_unit.get_absolute_Ln_approximation(k);
			R_xbig_prec = k;
			outR = R_xbig;
		}
		else {
			truncate(outR, R_xbig, R_xbig.exponent()+k+1);
		}
	}

	return outR;
}



//
// quadratic_order::class_number
//
// Task:
//      returns the class number.  The algorithm is chosen based on the size
//      of the discriminant.
//

bigint
quadratic_order::class_number()
{
	debug_handler("quadratic_order", "class_number");

	int num;

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_h_computed()) {
		quadratic_order *QO = qi_class::get_current_order_ptr();
		long oprec = bigfloat::get_precision();

		bigfloat::set_precision(prec);
		qi_class::set_current_order(*this);

		// h has not been computed yet - compute it
		num = decimal_length(Delta);

		if (is_real()) {
			if (num < CNBJT_RB)
				class_group_BJT();
			if (num < CNSHANKS_RB)
				class_number_shanks();
			else
				class_group_siqs();
		}
		else {
			if (num < CNBJT_IB)
				class_group_BJT();
			if (num < CNSHANKS_IB)
				class_number_shanks();
			else
				class_group_siqs();
		}

		bigfloat::set_precision(oprec);
		qi_class::set_current_order(*QO);
	}

	return h;
}



//
// quadratic_order::class_group
//
// Task:
//      returns a vector containing the invariants of the class group.  The
//      algorithm is chosen based on the size of the discriminant.
//

base_vector< bigint >
quadratic_order::class_group()
{
	debug_handler("quadratic_order", "class_group");

	int num;

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_CL_computed()) {
		quadratic_order *QO = qi_class::get_current_order_ptr();
		long oprec = bigfloat::get_precision();

		bigfloat::set_precision(prec);
		qi_class::set_current_order(*this);

		// CL has not been computed yet - compute it

		if (h.is_gt_zero()) {
			// class number is known
//      class_group_h();
			h.assign_zero();
			class_group();
		}
		else {
			num = decimal_length(Delta);

			if (is_real()) {
				if (num < CGBJT_RB)
					class_group_BJT();
				else if (num < CGRAND_RB)
					class_group_randexp();
				else
					class_group_siqs();
			}
			else {
				if (num < CGBJT_IB)
					class_group_BJT();
				else if (num < CGSHANKS_IB)
					class_group_shanks();
				else
					class_group_siqs();
			}
		}

		bigfloat::set_precision(oprec);
		qi_class::set_current_order(*QO);
	}

	return CL;
}



//
// quadratic_order::regulator(bool)
//
// Task:
//      returns an approximation of the regulator.  If the parameter is true,
//      then the regulator will be computed using the subexponential algorithm
//      and recomputed if it was previously computed with a different
//      algorithm.
//

bigfloat
quadratic_order::regulator(bool usesub)
{
	debug_handler("quadratic_order", "regulator(bool)");

	bigfloat outR;
	quadratic_order *QO = qi_class::get_current_order_ptr();
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);
	qi_class::set_current_order(*this);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (usesub) {
		if ((!is_R_computed()) || ((is_R_computed()) && (!is_subexp_computed()))) {
			assign(Delta);
			class_group_siqs();
		}
	}
	else
		regulator();

	outR.assign(R);
	bigfloat::set_precision(oprec);
	qi_class::set_current_order(*QO);

	return outR;
}



//
// quadratic_order::class_group(bool)
//
// Task:
//      returns a vector containing the invariants of the class group.  If the
//      parameter is true, then the class group will be computed using the
//      subexponential algorithm, and recomputed if it was previously computed
//      with a different algorithm.
//

base_vector< bigint >
quadratic_order::class_group(bool usesub)
{
	debug_handler("quadratic_order", "class_group(bool)");

	quadratic_order *QO = qi_class::get_current_order_ptr();
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);
	qi_class::set_current_order(*this);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (usesub) {
		if ((!is_CL_computed()) ||
		    ((is_CL_computed()) && (!is_subexp_computed()))) {
			assign(Delta);
			class_group_siqs();
		}
	}
	else
		class_group();

	bigfloat::set_precision(oprec);
	qi_class::set_current_order(*QO);

	return CL;
}



//
// quadratic_order::regulator_BJT
//
// Task:
//      computes an approximation of the regulator using a variation of
//      baby-step giant-step due to Buchmann.
//

bigfloat
quadratic_order::regulator_BJT()
{
	debug_handler("quadratic_order", "regulator_BJT");

	qi_class_real A, UNIT;
	timer t;
	bigfloat outR;
	quadratic_order *QO = qi_class::get_current_order_ptr();
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);
	qi_class::set_current_order(*this);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_R_computed()) {
		if (is_imaginary())
			R.assign_one();
		else {
			if (prin_list.no_of_elements() > 0) {
				A = prin_list.last_entry();
				if ((A.is_one()) && (!R.is_zero())) {
					UNIT.assign_one();
					A.assign(nearest(UNIT, A.get_distance()));
					R.assign(A.get_distance());
				}
			}

			if (R.is_zero()) {
				if (info)
					t.start_timer();

				rBJT();

				if (info) {
					t.stop_timer();
					std::cout << "\nelapsed CPU time (regulator_BJT) = ";
					MyTime(t.user_time());
					std::cout << "\n";
				}
			}

			if (do_verify) {
				if (verify_regulator()) {
					if (info)
						std::cout << "Regulator is correct." << std::endl;
				}
				else
					std::cout << "Regulator is not correct!" << std::endl;
			}
		}
	}

	outR.assign(R);
	bigfloat::set_precision(oprec);
	qi_class::set_current_order(*QO);

	return outR;
}



//
// quadratic_order::regulator_shanks
//
// Task:
//      computes an approximation of the regulator using the O(D^1/5) algorithm
//      due to Mollin and Williams.
//

bigfloat
quadratic_order::regulator_shanks()
{
	debug_handler("quadratic_order", "regulator_shanks");

	qi_class_real A, UNIT;
	timer t;
	bigfloat outR;
	quadratic_order *QO = qi_class::get_current_order_ptr();
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);
	qi_class::set_current_order(*this);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_R_computed()) {
		if (is_imaginary())
			R.assign_one();
		else {
			if (prin_list.no_of_elements() > 0) {
				A = prin_list.last_entry();
				if ((A.is_one()) && (!R.is_zero()))
					UNIT.assign_one();
				A.assign(nearest(UNIT, A.get_distance()));
				R.assign(A.get_distance());
			}

			if (R.is_zero()) {
				if (info)
					t.start_timer();

				rshanks();

				if (info) {
					t.stop_timer();
					std::cout << "\nelapsed CPU time (regulator_shanks) = ";
					MyTime(t.user_time());
					std::cout << "\n";
				}
			}

			if (do_verify) {
				if (verify_regulator()) {
					if (info)
						std::cout << "Regulator is correct." << std::endl;
				}
				else
					std::cout << "Regulator is not correct!" << std::endl;
			}
		}
	}

	outR.assign(R);
	bigfloat::set_precision(oprec);
	qi_class::set_current_order(*QO);

	return outR;
}



//
// quadratic_order::verify_regulator()
//
// Task:
//      verifies whether the regulator is correct.
//

bool
quadratic_order::verify_regulator()
{
	debug_handler("quadratic_order", "verify_regulator");

	bool is_reg;
	timer t;
	quadratic_order *QO = qi_class::get_current_order_ptr();

	if (!is_real())  return true;

	qi_class::set_current_order(*this);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_R_computed())
		regulator();

	if (info)
		t.start_timer();


	if (info > 1) {
		std::cout << "\nVerifying regulator:" << std::endl;
		std::cout << "  R = " << R << std::endl;
	}

	//<MM>
	// Requires that R is an absolute 3-approximation
	// to the regulator.
	//

	is_reg = could_be_regulator_multiple(R);

	//</MM>


	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (verify_regulator) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	qi_class::set_current_order(*QO);

	return is_reg;
}



//
// quadratic_order::class_number_shanks
//
// Task:
//      computes the class number using the O(D^1/5) algorithm due to Shanks.
//

bigint
quadratic_order::class_number_shanks()
{
	debug_handler("quadratic_order", "class_number_shanks");

	timer t;

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_h_computed()) {
		quadratic_order *QO = qi_class::get_current_order_ptr();
		long oprec = bigfloat::get_precision();

		bigfloat::set_precision(prec);
		qi_class::set_current_order(*this);

		if (info)
			t.start_timer();

		if (is_imaginary())
			cns_imag();
		else
			cns_real();

		if (info) {
			t.stop_timer();
			std::cout << "\nelapsed CPU time (class_number_shanks) = ";
			MyTime(t.user_time());
			std::cout << "\n";
		}

		bigfloat::set_precision(oprec);
		qi_class::set_current_order(*QO);
	}

	return h;
}



//
// quadratic_order::class_group_BJT
//
// Task:
//      computes the class group using a variation of baby-step giant-step
//      due to Buchmann, Jacobson, and Teske.
//

base_vector< bigint >
quadratic_order::class_group_BJT()
{
	debug_handler("quadratic_order", "class_group_BJT");

	timer t;

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_CL_computed()) {
		quadratic_order *QO = qi_class::get_current_order_ptr();
		long oprec = bigfloat::get_precision();

		bigfloat::set_precision(prec);
		qi_class::set_current_order(*this);

		if (info)
			t.start_timer();

		if (is_imaginary())
			cgBJT_imag();
		else
			cgBJT_real();

		if (info) {
			t.stop_timer();
			std::cout << "\nelapsed CPU time (class_group_BJT) = ";
			MyTime(t.user_time());
			std::cout << "\n";
		}

		if (do_verify) {
			if (verify_class_group()) {
				if (info)
					std::cout << "Class group is correct." << std::endl;
			}
			else
				std::cout << "Class group is not correct!" << std::endl;
		}

		bigfloat::set_precision(oprec);
		qi_class::set_current_order(*QO);
	}

	return CL;
}



//
// quadratic_order::class_group_shanks
//
// Task:
//      computes the class group using an unpublished variation of Shanks
//      baby-step giant-step algorithm.  It is essentially the same as the
//      class number algorith, except we continue computing the subgroup until
//      we have the entire class group.
//

base_vector< bigint >
quadratic_order::class_group_shanks()
{
	debug_handler("quadratic_order", "class_group_shanks");

	bigint temp;
	bigfloat d1, d2;
	timer t;

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_CL_computed()) {
		quadratic_order *QO = qi_class::get_current_order_ptr();
		long oprec = bigfloat::get_precision();

		bigfloat::set_precision(prec);
		qi_class::set_current_order(*this);

		d1.assign(0);
		d2.assign(0);

		if (info)
			t.start_timer();

		if (is_imaginary())
			cgs_imag(d1, d2, temp);
		else
			cgs_real(d1, temp);

		if (info) {
			t.stop_timer();
			std::cout << "\nelapsed CPU time (class_group_shanks) = ";
			MyTime(t.user_time());
			std::cout << "\n";
		}

		if (do_verify) {
			if (verify_class_group()) {
				if (info)
					std::cout << "Class group is correct." << std::endl;
			}
			else
				std::cout << "Class group is not correct!" << std::endl;
		}

		bigfloat::set_precision(oprec);
		qi_class::set_current_order(*QO);
	}

	return CL;
}



//
// quadratic_order::class_group_h
//
// Task:
//      computes the class group given the class number.
//

base_vector< bigint >
quadratic_order::class_group_h()
{
	debug_handler("quadratic_order", "class_group_h");

	timer t;

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_CL_computed()) {
		quadratic_order *QO = qi_class::get_current_order_ptr();
		long oprec = bigfloat::get_precision();

		bigfloat::set_precision(prec);
		qi_class::set_current_order(*this);

		if (info)
			t.start_timer();

		// compute h if necessary
		if (h.is_zero());
		class_number();

		if (is_imaginary())
			cgh_imag();
		else
			cgh_real();

		if (info) {
			t.stop_timer();
			std::cout << "\nelapsed CPU time (class_group_h) = ";
			MyTime(t.user_time());
			std::cout << "\n";
		}

		if (do_verify) {
			if (verify_class_group()) {
				if (info)
					std::cout << "Class group is correct." << std::endl;
			}
			else
				std::cout << "Class group is not correct!" << std::endl;
		}

		bigfloat::set_precision(oprec);
		qi_class::set_current_order(*QO);
	}

	return CL;
}



//
// quadratic_order::class_group_randexp
//
// Task:
//      computes the class group using a subexponential algorithm with
//      random exponent vectors (Duellmann's variation)
//

base_vector< bigint >
quadratic_order::class_group_randexp()
{
	debug_handler("quadratic_order", "class_group_randexp");

	int decimal_digits;
	timer t;

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_CL_computed()) {
		quadratic_order *QO = qi_class::get_current_order_ptr();
		long oprec = bigfloat::get_precision();

		bigfloat::set_precision(prec);
		qi_class::set_current_order(*this);

		if (info)
			t.start_timer();

		decimal_digits = decimal_length(Delta);
		if (decimal_digits < 7)
			class_group();
		else {
			subused = true;
			cg_randexp();
		}

		if (info) {
			t.stop_timer();
			std::cout << "\nelapsed CPU time (class_group_randexp) = ";
			MyTime(t.user_time());
			std::cout << "\n";
		}

		if (do_verify) {
			if (verify_regulator()) {
				if (info)
					std::cout << "Regulator is correct." << std::endl;
			}
			else
				std::cout << "Regulator is not correct!" << std::endl;

			if (verify_cg(false)) {
				if (info)
					std::cout << "Class group is correct." << std::endl;
			}
			else
				std::cout << "Class group is not correct!" << std::endl;
		}

		bigfloat::set_precision(oprec);
		qi_class::set_current_order(*QO);
	}

	return CL;
}



//
// quadratic_order::class_group_mpqs
//
// Task:
//      computes the class group using an algorithm based on the MPQS
//      factoring algorithm due to Jacobson.
//

base_vector< bigint >
quadratic_order::class_group_mpqs()
{
	debug_handler("quadratic_order", "class_group_mpqs");

	int decimal_digits;
	timer t;

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_CL_computed()) {
		quadratic_order *QO = qi_class::get_current_order_ptr();
		long oprec = bigfloat::get_precision();

		bigfloat::set_precision(prec);
		qi_class::set_current_order(*this);

		if (info)
			t.start_timer();

		decimal_digits = decimal_length(Delta);
		if (decimal_digits < 7)
			class_group();
		else {
			subused = true;
			cg_mpqs(false, false);
		}

		if (info) {
			t.stop_timer();
			std::cout << "\nelapsed CPU time (class_group_mpqs) = ";
			MyTime(t.user_time());
			std::cout << "\n";
		}

		if (do_verify) {
			if (verify_regulator()) {
				if (info)
					std::cout << "Regulator is correct." << std::endl;
			}
			else
				std::cout << "Regulator is not correct!" << std::endl;

			if (verify_class_group()) {
				if (info)
					std::cout << "Class group is correct." << std::endl;
			}
			else
				std::cout << "Class group is not correct!" << std::endl;
		}

		bigfloat::set_precision(oprec);
		qi_class::set_current_order(*QO);
	}

	return CL;
}



//
// quadratic_order::class_group_siqs
//
// Task:
//      computes the class group using an algorithm based on the
//      self-initializing MPQS factoring algorithm due to Jacobson.
//

base_vector< bigint >
quadratic_order::class_group_siqs()
{
	debug_handler("quadratic_order", "class_group_siqs");

	int decimal_digits;
	timer t;

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_CL_computed()) {
		quadratic_order *QO = qi_class::get_current_order_ptr();
		long oprec = bigfloat::get_precision();

		bigfloat::set_precision(prec);
		qi_class::set_current_order(*this);

		if (info)
			t.start_timer();

		decimal_digits = decimal_length(Delta);
		if (decimal_digits < 7)
			class_group();
		else {
			subused = true;
			cg_mpqs(true, false);
		}

		if (info) {
			t.stop_timer();
			std::cout << "\nelapsed CPU time (class_group_siqs) = ";
			MyTime(t.user_time());
			std::cout << "\n";
		}

		if (do_verify) {
			if (verify_regulator()) {
				if (info)
					std::cout << "Regulator is correct." << std::endl;
			}
			else
				std::cout << "Regulator is not correct!" << std::endl;

			if (verify_class_group()) {
				if (info)
					std::cout << "Class group is correct." << std::endl;
			}
			else
				std::cout << "Class group is not correct!" << std::endl;
		}

		bigfloat::set_precision(oprec);
		qi_class::set_current_order(*QO);
	}

	return CL;
}



//
//
// quadratic_order::class_group_lp
//
// Task:
//      computes the class group using large primes an algorithm based on the
//      self-initializing MPQS factoring algorithm due to Jacobson.
//

base_vector< bigint >
quadratic_order::class_group_lp()
{
	debug_handler("quadratic_order", "class_group_lp");

	int decimal_digits;
	timer t;

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (!is_CL_computed()) {
		quadratic_order *QO = qi_class::get_current_order_ptr();
		long oprec = bigfloat::get_precision();

		bigfloat::set_precision(prec);
		qi_class::set_current_order(*this);

		if (info)
			t.start_timer();

		decimal_digits = decimal_length(Delta);
		if (decimal_digits < 7)
			class_group();
		else {
			subused = true;
			cg_mpqs(true, true);
		}

		if (info) {
			t.stop_timer();
			std::cout << "\nelapsed CPU time (class_group_lp) = ";
			MyTime(t.user_time());
			std::cout << "\n";
		}

		if (do_verify) {
			if (verify_regulator()) {
				if (info)
					std::cout << "Regulator is correct." << std::endl;
			}
			else
				std::cout << "Regulator is not correct!" << std::endl;

			if (verify_class_group()) {
				if (info)
					std::cout << "Class group is correct." << std::endl;
			}
			else
				std::cout << "Class group is not correct!" << std::endl;
		}

		bigfloat::set_precision(oprec);
		qi_class::set_current_order(*QO);
	}

	return CL;
}



//
// quadratic_order::count_primes_to_verify(int BachBound)
//
// Task:
//      counts the number of primes not in the factor base but less than
//      Bach's bound.
//

int
quadratic_order::count_primes_to_verify(int BachBound)
{
	debug_handler("quadratic_order", "count_primes_to_verify(int)");

	int nump = 0, p, size_FB;
	qi_class P;
	prime_list PList;
	quadratic_order *QO = qi_class::get_current_order_ptr();

	qi_class::set_current_order(*this);

	if (fact_base.size() == 0) {
		qo_read_parameters(decimal_length(Delta), size_FB);
		qo_get_fb(size_FB);
		fact_base[fact_base.size()-1].get_a().intify(p);
		fact_base.reset();
	}
	else
		fact_base[fact_base.size()-1].get_a().intify(p);

	if (p == 2)
		++p;
	else
		p += 2;
	if (p <= BachBound) {
		PList.resize(p, BachBound+10);
		p = PList.get_first_prime();
		while ((p <= BachBound) && (p > 0)) {
			if (generate_prime_ideal(P, p))
				++nump;
			p = PList.get_next_prime();
		}
	}

	qi_class::set_current_order(*QO);

	return nump;
}



//
// quadratic_order::verify_class_group(int level)
//
// Task:
//      verifies whether the class group is correct using the given level
//

bool
quadratic_order::verify_class_group(int level)
{
	debug_handler("quadratic_order", "verify_class_group(int)");

	int templ;
	bool temp = true;;

	if (level) {
		templ = do_verify;
		do_verify = level;

		if (sieve.is_initialized())
			temp = verify_cg(true);
		else
			temp = verify_cg(false);

		do_verify = templ;
	}

	return temp;
}



//
// quadratic_order::verify_class_group()
//
// Task:
//      verifies whether the class group is correct.
//

bool
quadratic_order::verify_class_group()
{
	debug_handler("quadratic_order", "verify_class_group");

	if (sieve.is_initialized())
		return verify_cg(true);
	else
		return verify_cg(false);
}



//
// quadratic_order::verify_cg()
//
// Task:
//      verifies whether the class group is correct.
//

bool
quadratic_order::verify_cg(bool use_sieve)
{
	debug_handler("quadratic_order", "verify_cg");

	lidia_size_t i, j, numFB, numCL, matcols, matrows;
	int p, BachBound = 0, num, nump = 0, countp = 0;
	bigint ord;
	qi_class A, C, P;
	math_vector< bigint > vec, rvec, nvec, pvec;
	sort_vector< int > D_primes;
	base_vector< qi_class > new_FB;
	matrix< bigint > vmat;
	bool OK;
	timer t, tprime;
	prime_list PList;
	quadratic_number_standard q;
	quadratic_order *QO = qi_class::get_current_order_ptr();
	long oprec = bigfloat::get_precision();

	// FOR TABLE GENERATION
	if (special)
		t.start_timer();

	bigfloat::set_precision(prec);
	qi_class::set_current_order(*this);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (info)
		tprime.start_timer();

	if (info)
		t.start_timer();

	OK = true;

	// verify primes less than Bach bound
	numFB = fact_base.size();
	numCL = CL.size();
	if (numFB == 0) {
		qo_read_parameters(decimal_length(Delta), numFB);
		std::cout << "\nInput factor base size:  ";
		std::cin >> numFB;
		std::cout << std::endl;
		qo_get_fb(numFB);
	}


	if (do_verify > 1) {
		BachBound = bach_bound();

#if 0
		std::cout << "Input bach bound:  ";
		std::cin >> BachBound;
		std::cout << std::endl;
#endif

		nump = count_primes_to_verify(BachBound);

		// FOR TABLE GENERATION
		if (special)
			std::cout << nump << "  " << std::flush;
	}
	else {
		// FOR TABLE GENERATION
		if (special)
			std::cout << "0  " << std::flush;
	}

	if (info > 1) {
		std::cout << "\nVerifying class group (level " << do_verify << "):" << std::endl;
		std::cout << "  CL = " << CL << std::endl;
		std::cout << "  generators = " << gens << std::endl;
	}

	num = decimal_length(Delta);
	if (is_subexp_computed() || !is_CL_computed()) {
		if (info > 1) {
			if (do_verify == 2) {
				std::cout << "  BachBound = " << BachBound << std::endl;
				std::cout << "Representing each prime ideal less than Bach's bound over "
					"the factor base:" << std::endl;
				std::cout << nump << " primes to verify!" << std::endl;
			}
			else if (do_verify > 2) {
				std::cout << "  BachBound = " << BachBound << std::endl;
				std::cout << "Representing each prime ideal less than Bach's bound over "
					"the generators:" << std::endl;
				std::cout << nump+fact_base.size() << " primes to verify!" << std::endl;
			}
		}

		pvec.set_mode(EXPAND);
		nvec.set_mode(EXPAND);
		vec.set_mode(EXPAND);
		vec.set_size(numFB);
		rvec.set_capacity(numCL);
		rvec.set_size(numCL);

		if (do_verify > 2) {
			generators();

			// verify primes in factor base
			for (i = 0; i < numFB; ++i) {
				for (j = 0; j < numFB; ++j)
					vec[j].assign_zero();
				vec[i].assign_one();
				rvec = represent_over_generators(vec);

				C.assign_one();
				for (j = 0; j < numCL; ++j) {
					power(A, gens[j], rvec[j]);
					multiply(C, C, A);
				}

				if (fact_base[i].is_equivalent(C)) {
					if (info > 2) {
						std::cout << "  " << fact_base[i] << " ~ " << rvec << " ~ " << C << std::endl;
					}
				}
				else {
					OK = false;
					if (info > 2) {
						std::cout << "  ERROR:  " << fact_base[i] << " <>" << rvec << " <>";
						std::cout << C << std::endl;
					}
				}
			}
		}

		if (do_verify > 1) {
			matcols = 0;
			vmat.set_representation(matrix_flags::sparse_representation);
			vmat.set_orientation(matrix_flags::column_oriented);

			new_FB.set_mode(EXPAND);
			new_FB = fact_base;

			fact_base[numFB-1].get_a().intify(p);
			p += 2;
			if (p <= BachBound) {
				PList.resize(p, BachBound+10);
			}
			p = PList.get_first_prime();
			while ((p <= BachBound) && (p > 0)) {
				if (generate_prime_ideal(P, p)) {
					TR.touch_files();
					++countp;

					if ((info > 2) && (info <= 4)) {
						if ((countp % 200) == 0) {
							tprime.stop_timer();
							std::cout << countp << "   ";
							MyTime(tprime.user_time());
							std::cout << std::endl;
							tprime.cont_timer();
						}
					}

					if (info > 4) {
						std::cout << "P = " << P << std::flush;
					}

					if (use_sieve) {
						if (do_verify > 2)
							vec = represent_over_FB_sieve(P, q);
						else
							represent_over_FB_sieve(P);
					}
					else
						vec = represent_over_FB(P, q);

					if (info > 4) {
						std::cout << " --- done!" << std::endl;
					}

					if (do_verify > 2) {
						pvec = vec;

						// compute representation over original factor base
						for (j = vec.size()-1; j >= numFB; --j) {
							if (!vec[j].is_zero()) {
								nvec = (math_vector< bigint > ) vmat.get_column_vector(j-numFB);
								multiply(nvec, nvec, vec[j]);
								vec[j].assign_zero();
								nvec.set_size(vec.size());
								add(vec, vec, nvec);

							}
						}
						vec.set_size(numFB);


						// compute representation over generators
						fact_base.set_size(numFB);

						rvec = represent_over_generators(vec);

						fact_base = new_FB;

						C.assign_one();
						for (j = 0; j < numCL; ++j) {
							power(A, gens[j], rvec[j]);
							multiply(C, C, A);
						}

						if (P.is_equivalent(C)) {
							if (info > 2) {
								std::cout << "  " << P << " ~ " << rvec << " ~ " << C << std::endl;
							}
						}
						else {
							OK = false;
							if (info > 2) {
								std::cout << "  ERROR:  " << P << " = " << rvec << " <>" << C;
								std::cout << std::endl;
							}
						}
					}

					// add verified prime to factor base
					if ((fact_base.size() < (numFB*4)) && (fact_base.size() < 8000)
					    && (num > 28)) {
						new_FB[new_FB.size()] = P;
						fact_base[fact_base.size()] = P;

						if (do_verify > 2) {
							matrows = pvec.size();
							vmat.resize(matrows, matcols+1);
							vmat.sto_column_vector(pvec, matrows, matcols);
							++matcols;
						}
					}
				}
				p = PList.get_next_prime();
			}

			fact_base.set_size(numFB);
		}
	}

	if (do_verify > 0 && is_CL_computed()) {
		generators();

		// verify orders of generators
		if (info > 1)
			std::cout << "Checking orders of generators:" << std::endl;

		for (i = 0; i < numCL; ++i) {
			A.assign(gens[i]);
			ord.assign(A.order_in_CL());
			if (ord == CL[i]) {
				if (info > 2) {
					std::cout << "  ord(G[" << i << "]) = " << ord << " = " << CL[i] << std::endl;
				}
			}
			else {
				OK = false;
				if (info > 2) {
					std::cout << "  ERROR:  ord(G[" << i << "]) = " << ord << " <>" << CL[i];
					std::cout << std::endl;
				}
			}
		}
	}

	if (info) {
		tprime.stop_timer();
		std::cout << "\nelapsed CPU time (verify_class_group) = ";
		MyTime(tprime.user_time());
		std::cout << "\n";
	}

	bigfloat::set_precision(oprec);
	qi_class::set_current_order(*QO);

	// FOR TABLE GENERATION
	if (special) {
		t.stop_timer();
		std::cout << t.user_time() << std::endl;
	}

	return OK;
}



//
// quadratic_order::is_in_lattice(qi_class_real)
//
// Task:
//      determines whether the factorization of A over the factorbase used
//      by the subexponential algorithm is in the lattice of vectors
//      representing the principal class.
//

bool
quadratic_order::is_in_lattice(const qi_class_real & AA)
{
	debug_handler("quadratic_order", "is_in_lattice(qi_class_real)");

	math_vector< bigint > vec;
	bool isin;
	quadratic_number_standard q;
	quadratic_order *QO = qi_class::get_current_order_ptr();

	// set last-used quadratic_order
	qo_l().set_last(this);

	qi_class::set_current_order(*this);

	vec = represent_over_FB_sieve(AA, q);
	isin = is_in_lattice(vec);

	qi_class::set_current_order(*QO);

	return isin;
}



//
// quadratic_order::is_in_lattice(math_vector)
//
// Task:
//      determines whether the factorization of A over the factorbase used
//      by the subexponential algorithm is in the lattice of vectors
//      representing the principal class.
//

bool
quadratic_order::is_in_lattice(math_vector< bigint > & vec)
{
	debug_handler("quadratic_order", "is_in_lattice(math_vector)");

	// set last-used quadratic_order
	qo_l().set_last(this);

	math_vector< bigint > x;
	bool is_in_lat;

	is_in_lat = RHNF.solve_hnf(vec, x);
	return is_in_lat;
}



//
// quadratic_order::is_in_lattice(qi_class_real, bigfloat)
//
// Task:
//      determines whether the factorization of A over the factorbase used
//      by the subexponential algorithm is in the lattice of vectors
//      representing the principal class.  If so, the distance from (1) is
//      also returned.
//

bool
quadratic_order::is_in_lattice(const qi_class_real & AA, bigfloat & dist)
{
	debug_handler("quadratic_order", "is_in_lattice(qi_class_real, bigfloat)");

	lidia_size_t i;
	bigint max, DET2;
	bigfloat ldist, rd2, a, b, rtemp1, rtemp2, maxlog;
	qi_class_real UNIT, C, D;
	math_vector< bigint > vec, x;
	math_vector< bigfloat > logs;
	bool in_lat;
	quadratic_order *QO = qi_class::get_current_order_ptr();
	long oprec = bigfloat::get_precision(), xprec2;
	quadratic_number_standard q;

	xbigfloat ldist_x;
	timer t;

	// set last-used quadratic_order
	qo_l().set_last(this);

	qi_class::set_current_order(*this);

	bigfloat::set_precision(prec);

	logs.set_capacity(minima_qn.size());
	logs.set_size(minima_qn.size());

	// FOR TABLE GENERATION
	if (special2)
		t.start_timer();

	// find dependency over relation matrix
	vec = represent_over_FB_sieve(AA, q);

	in_lat = RHNF.solve_hnf(vec, x);

	if (in_lat) {
		if (info > 3) {
			C.assign_one();
			qi_class_real PP;
			for (i = 0; i < vec.size(); ++i)
				if (!vec[i].is_zero()) {
					power_real(PP, static_cast<qi_class_real>(fact_base[i]), vec[i]);
					C *= PP;
				}
			std::cout << "A = " << AA << std::endl;
			std::cout << "C = " << C << std::endl;
		}

		DET2 = 1;
		math_vector< bigint > w;
		w.set_mode(EXPAND);

		// compute vector
		x = TR.multiply(vec);

		// remove common factors with DET
		if (DET != 1) {
			bigint G = DET;
			for (i = 0; i < x.size(); ++i)
				G = gcd(G, x[i]);
			if (info > 3) {
				std::cout << "DET = " << DET << std::endl;
				std::cout << "G = " << G << std::endl;
				std::cout << "DET/G = " << DET/G << std::endl;
			}
			if (DET.is_lt_zero())
				G.negate();
			if (abs(G) > 1)
				divide(x, x, G);
			divide(DET2, DET, G);
		}

		matrix< bigint > B;
		B.set_storage_mode(matrix_flags::sparse_representation);
		B.resize(AMAT.get_no_of_rows(), AMAT.get_no_of_columns());
		for (lidia_size_t ii1 = 0; ii1 < AMAT.get_no_of_rows(); ++ii1)
			for (lidia_size_t ii2 = 0; ii2 < AMAT.get_no_of_columns(); ++ii2)
				if (AMAT.member(ii1, ii2))
					B.sto(ii1, ii2, bigint(AMAT.member(ii1, ii2)));

		multiply(w, B, x);

		if (info > 3) {
			std::cout << "DET2 = " << DET2 << std::endl;
			std::cout << "Ax = " << w << std::endl;
			std::cout << "vec = " << vec << std::endl;
			if (DET2 != 1)
				std::cout << "Ax/DET2 = " << w / DET2 << std::endl;
		}

		if ((w/DET2) != vec) {
			subtract(w, w, vec);
			if (info > 3)
				std::cout << "Ax - vec = " << w << std::endl;

			if (!RHNF.solve_hnf(w, w)) {
				std::cout << "ERROR!  Couldn't solve Hx = w!!!" << std::endl;
				std::cout << "w = " << w << std::endl;
			}

			w = TR.multiply(w);
			subtract(x, x, w);

			multiply(w, B, x);
			if (info > 3) {
				std::cout << "Ax = " << w << std::endl;
				std::cout << "vec = " << vec << std::endl;
				if (DET2 != 1)
					std::cout << "Ax/DET2 = " << w / DET2 << std::endl;
			}
		}

		// FOR TABLE GENERATION
		if (special2) {
			t.stop_timer();
			std::cout << t.user_time() << std::endl;
			t.start_timer();
		}


		// Markus' Method
		bigint q;
		base_vector< xbigfloat > lbw, lbr;
		base_vector< long > plbw, plbr;
		long br, lowregbound, e_tmp;
		xbigfloat Rx, ldist_x, LL;
		math_vector< bigint > fud;
		bigfloat tmp;

		if (info > 3) {
			max = x[0];
			for (i = 1; i < x.size(); ++i)
				if (x[i] > max)
					max = x[i];
			std::cout << "MAX in x:  " << decimal_length(max) << " decimal digits" << std::endl;
		}


		//<MM>
		//
		// NOTE: The reduction modulo the regulator is a hack that must be
		//       removed by using a not yet existing function
		//       quadratic_number_power_product::short_principal_ideal_generator.
		//
		//       !! HERE, we assume that the basis of "fund_unit" is "minima_qn" !!
		//
		//</MM>

//    const math_vector<bigint> & fund_unit_exp = fund_unit.get_exponents();
		math_vector< bigint > fund_unit_exp = fund_unit.get_exponents();
		fund_unit_exp.set_size(minima_qn.size());

		if (DET2 == 1) {
			Rx.assign(R_xbig * fu_mult);
			br = Rx.b_value();

			// compute distance modulo R
			reduce_modulo_regulator(q, minima_qn, x, lbw, plbw, minima_qn, fund_unit_exp, lbr, plbr, br);

			// compute new vector
			if (info > 3) {
				std::cout << "q = " << decimal_length(q) << " decimal digits" << std::endl;
			}
			x = x - q*fund_unit_exp;
		}
		else {
			Rx.assign(R_xbig * fu_mult * DET2);
			br = Rx.b_value();

			fud = DET2*fund_unit_exp;

			// compute distance modulo R
			reduce_modulo_regulator(q, minima_qn, x, lbw, plbw, minima_qn, fud, lbr, plbr, br);

			// compute new vector
			if (info > 3) {
				std::cout << "q = " << decimal_length(q) << " decimal digits" << std::endl;
			}
			x = x - q*fud;
		}

		if (info > 3) {
			max = x[0];
			for (i = 1; i < x.size(); ++i)
				if (x[i] > max)
					max = x[i];
			std::cout << "MAX in x:  " << decimal_length(max) << " decimal digits" << std::endl;
		}


		lowregbound = lower_regulator_bound(LL, Delta, h, 0);

		if (DET2 == 1) {
			xprec2 = static_cast<long>(std::ceil(static_cast<double>(bigfloat::get_precision()) *
							     std::log(10.0L) / std::log(2.0L)));

			relative_Ln_approximation (ldist_x, minima_qn, x, xprec2, lowregbound);
		}
		else {
			max = abs(x[0]);
			for (i = 1; i < x.size(); ++i)
				if (abs(x[i]) > max)
					max.assign(abs(x[i]));

			xprec2 = static_cast<long>(std::ceil(static_cast<double>(4*decimal_length(max) +
										 5 + decimal_length(DET2)) *
							     std::log(10.0L) / std::log(2.0L)));

			relative_Ln_approximation (ldist_x, minima_qn, x, xprec2, lowregbound);

			std::cout << ldist_x << std::endl;
			ldist.assign(ldist_x.get_mantissa());
			e_tmp = ldist_x.get_exponent()-ldist_x.get_mantissa().bit_length();
			if (e_tmp > 0)
				shift_left(ldist, ldist, e_tmp);
			else
				shift_right(ldist, ldist, -e_tmp);

			std::cout << "ldist = " << ldist << std::endl;

			divide(ldist_x, ldist_x, xbigfloat(DET2), xprec2);

			std::cout << ldist_x << std::endl;
		}


		// store in bigfloat
		ldist.assign(ldist_x.get_mantissa());
		e_tmp = ldist_x.get_exponent()-ldist_x.get_mantissa().bit_length();
		if (e_tmp > 0)
			shift_left(ldist, ldist, e_tmp);
		else
			shift_right(ldist, ldist, -e_tmp);

		if (info > 3) {
			std::cout << "ldist = " << ldist << std::endl;
		}

		// final reduction mod R
		while (ldist < 0)
			add(ldist, ldist, R);
		while (ldist > R)
			subtract(ldist, ldist, R);

		if (info > 3)
			std::cout << "ldist = " << ldist << std::endl;

		// FOR TABLE GENERATION
		if (special2) {
			t.stop_timer();
			std::cout << t.user_time() << std::endl;
		}

		// convert to standard formulation
		UNIT.assign_one();
		dist = AA.convert_distance(UNIT, ldist);
	}
	else {
		// FOR TABLE GENERATION
		if (special2) {
			t.stop_timer();
			std::cout << t.user_time() << std::endl;
			std::cout << "0" << std::endl;
		}

		dist.assign(R);
	}

	qi_class::set_current_order(*QO);
	bigfloat::set_precision(oprec);

	return in_lat;
}



//
// quadratic_order::factor_h
//
// Task:
//      computes the class number and factors it.
//

rational_factorization
quadratic_order::factor_h()
{
	debug_handler("quadratic_order", "factor_h");

	// set last-used quadratic_order
	qo_l().set_last(this);

	// compute class number, if needed
	if (!is_h_computed())
		class_number();

	// factor h, if necessary
	if (h.is_gt_zero()) {
		if (hfact != h) {
			hfact.assign(h);
			hfact.factor();
		}
	}

	return hfact;
}



//
// quadratic_order::exponent
//
// Task:
//      computes the exponent of the class group.
//

bigint
quadratic_order::exponent()
{
	debug_handler("quadratic_order", "exponent");

	bigint expo;

	// set last-used quadratic_order
	qo_l().set_last(this);

	// compute class group, if needed
	if (!is_CL_computed())
		class_group();

	// return largest non-cyclic factor
	if (!is_CL_computed())
		expo.assign_zero();
	else
		expo.assign(CL[CL.size()-1]);

	return expo;
}



//
// quadratic_order::p_rank
//
// Task:
//      computes the p_rank of the class group.
//

int
quadratic_order::p_rank(const bigint & p)
{
	debug_handler("quadratic_order", "p_rank");

	int prank;
	lidia_size_t i, rank;
	bigint rem;

	// set last-used quadratic_order
	qo_l().set_last(this);

	// compute class group, if needed
	if (!is_CL_computed())
		class_group();

	if (h.is_zero())
		prank = 0;
	else {
		// compute p-rank
		rank = CL.size();
		i = 0;
		remainder(rem, CL[i], p);
		while ((!rem.is_zero()) && (i < rank)) {
			++i;
			if (i < rank)
				remainder(rem, CL[i], p);
		}
		prank = rank - i;
	}

	return prank;
}



//
// quadratic_order::generators
//
// Task:
//      computes a set of generators of the class group.
//

base_vector< qi_class >
quadratic_order::generators()
{
	debug_handler("quadratic_order", "generators");

	qi_class G, A;
	lidia_size_t i, j, numCL, numFB;
	quadratic_order *QO = qi_class::get_current_order_ptr();

	// set last-used quadratic_order
	qo_l().set_last(this);

	qi_class::set_current_order(*this);

	// compute class group, if needed
	if (!is_CL_computed())
		class_group();

	if ((gens.size() == 0) && (!Delta.is_zero())) {
		if (h.is_one()) {
			G.assign_one();
			gens[0] = G;
		}
		else if (fact_base.size() == 1)
			gens[0].assign(fact_base[0]);
		else {
			// reduce transformation matrix elements and compute inverse
			set_transformations();

			// columns of UINV are the exponents vectors for the factor base elements
			// which correspond to each generator

			numCL = CL.size();
			numFB = contributors.size();
			for (i = 0; i < numCL; ++i) {
				G.assign_one();
				for (j = 0; j < numFB; ++j) {
					power(A, fact_base[contributors[j]], UINV.member(j, i));
					G *= A;
				}
				gens[i] = G;
			}
		}
	}

	qi_class::set_current_order(*QO);

	return gens;
}



//
// quadratic_order::factor_discriminant
//
// Task:
//      factors the discriminant by computing the class group and utiliziting
//      its 2-Sylow subgroup.
//

rational_factorization
quadratic_order::factor_discriminant()
{
	debug_handler("quadratic_order", "factor_discriminant");

	timer t;
	quadratic_order *QO = qi_class::get_current_order_ptr();
	long oprec = bigfloat::get_precision();

	bigfloat::set_precision(prec);
	qi_class::set_current_order(*this);

	// set last-used quadratic_order
	qo_l().set_last(this);

	if (info)
		t.start_timer();

	if (!Delta.is_zero()) {
		if (is_imaginary())
			factor_imag();
		else
			factor_real();
	}

	if (info) {
		t.stop_timer();
		std::cout << "\nelapsed CPU time (factor) = ";
		MyTime(t.user_time());
		std::cout << "\n";
	}

	bigfloat::set_precision(oprec);
	qi_class::set_current_order(*QO);

	return disc_fact;
}



//<MM>

//
// ::fundamental_unit
//
// Returns the fundamental unit of the order.
//
//

const quadratic_number_power_product & quadratic_order
::fundamental_unit ()
{
	lidia_error_handler("quadratic_order::fundamental_unit()",
			    "Not yet implemented.");

	if (!is_fundamental_unit_computed()) {

		// fundamental unit has not been computed yet - compute it

		if (is_imaginary())
			fund_unit.assign_one(*this);

		else {

			int num;
			quadratic_order *QO = qi_class::get_current_order_ptr();
			long oprec = bigfloat::get_precision();

			bigfloat::set_precision(prec);
			qi_class::set_current_order(*this);

			// set last-used quadratic_order
			qo_l().set_last(this);

			num = decimal_length(Delta);
			if (num < RBJTB)
				regulator_BJT();
			else if (num < RSHANKSB)
				regulator_shanks();
			else
				class_group_siqs();

			bigfloat::set_precision(oprec);
			qi_class::set_current_order(*QO);

			// The result of class_group_siqs should be
			// a compact representation of the fundamental
			// unit or an absolute 4 approximation to the
			// regulator. Same applies for the regulator
			// functions.
			//
			if (!is_fundamental_unit_computed()) {
				fund_unit.assign(*this, R, 1);
			}
		}
	}

	return fund_unit;
}



//
// ::could_be_regulator_multiple
//
// If the function returns false, then l is not an absolute
// 3 approximation to the regulator.
//
bool quadratic_order
::could_be_regulator_multiple (const xbigfloat & l)
{
	debug_handler ("quadratic_order",
		       "could_be_regulator_multiple()");

	// Determine the accuracy, 2^{-k+2} < ln2.
	//
	long k = 3;

	// Determine minimum in O close to l
	//
	quadratic_ideal J, K;
	quadratic_number_power_product beta;
	xbigfloat b;

	K.assign(*this);
	J.assign(*this);
	J.order_close(beta, b, l, k+1);

	if (info > 1)
		std::cout << "  nearest from (1) = " << J << std::endl;

	// Found O again ?
	//
	if (J != K) {
		quadratic_ideal A1, A2;
		A1 = J; A1.rho();
		A2 = J; A2.inverse_rho();
		if (K != A1 && K != A2)
			return false;
	}
	return true;
}



//</MM>




//
// operator >>
//
// Task:
//      inputs a quadratic_order from the std::istream in.
//

std::istream & operator >> (std::istream & in, quadratic_order & QO)
{
	debug_handler("quadratic_order", "operator >>");

	bigint D;

	in >> D;
	QO.assign(D);

	return in;
}



//
// operator <<
//
// Task:
//      outputs a quadratic_order to the std::ostream out.
//

std::ostream & operator << (std::ostream & out, const quadratic_order & QO)
{
	debug_handler("quadratic_order", "operator << ");

	int num;
	long oprec = bigfloat::get_precision();

	num = decimal_length(QO.Delta);
	bigfloat::set_precision(QO.prec);

	out << "Quadratic Order:" << std::endl;
	out << "   Delta = " << QO.Delta << " (" << num << ")" << std::endl;
	if (!QO.Delta.is_zero()) {
		if (QO.disc_fact.no_of_comp() > 0)
			out << " = " << QO.disc_fact << std::endl;
		if (QO.is_R_computed())
			out << "   R = " << QO.R << std::endl;
		if (QO.is_h_computed())
			out << "   h = " << QO.h << std::endl;
		if (QO.is_CL_computed())
			out << "   CL = " << QO.CL << std::endl;
		if (QO.gens.size() > 0)
			out << "   generators = " << QO.gens << std::endl;
		if (QO.is_L_computed())
			out << "   L(1, X) = " << QO.L << std::endl;
		if (QO.is_C_computed()) {
			bigfloat::set_precision(8);
			out << "   C(Delta) = " << QO.Cfunc << std::endl;
		}
	}

	bigfloat::set_precision(oprec);

	return out;
}



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
#endif
