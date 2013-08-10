// -*- C++ -*-
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
//	Author	: Werner Backes (WB), Thorsten Lauer (TL)
//	Changes	: See CVS log
//
//==============================================================================================


#ifndef LIDIA_LATTICE_GENSYS_H_GUARD_
#define LIDIA_LATTICE_GENSYS_H_GUARD_



#ifndef HEADBANGER

//
// include files needed
//



#ifndef LIDIA_BI_LATTICE_BASIS_H_GUARD_
# include	"LiDIA/lattices/bi_lattice_basis.h"
#endif
#ifndef LIDIA_BF_LATTICE_BASIS_H_GUARD_
# include	"LiDIA/lattices/bf_lattice_basis.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



//
// computing mode defines
//

const sdigit double_mode = 0x00000001; // double = bigfloat
const sdigit bigint_mode = 0x00000010;
const sdigit bigfloat_mode = 0x00000100;

//
// IO mode defines
//

const sdigit own_io_mode = 0x00000001;
const sdigit pari_io_mode = 0x00000010;
const sdigit maple_io_mode = 0x00000100;
const sdigit mathematica_io_mode = 0x00001000;

//
// Class lattice_gensys
//

class lattice_gensys
{
public :
	sdigit             		comp_mode;
	sdigit             		io_mode;

	timer                       CT;

	bigint_lattice_basis*	lbi;
	bigfloat_lattice_basis*	lbf;

public :

	//
	// Constructors / Destructor
	//
	lattice_gensys(lidia_size_t, lidia_size_t);
	lattice_gensys(lidia_size_t, lidia_size_t, const double**);
	lattice_gensys(lidia_size_t, lidia_size_t, const bigint**);
	lattice_gensys(lidia_size_t, lidia_size_t, const bigfloat**);
	lattice_gensys(const lattice_gensys&);
	lattice_gensys(const math_matrix< bigint > &);
	lattice_gensys(const math_matrix< bigfloat > &);
	virtual ~lattice_gensys();

	//
	// Input / Output
	//
	friend std::istream& operator >> (std::istream&, lattice_gensys&);
	friend std::ostream& operator << (std::ostream&, const lattice_gensys&);

	//
	// Assignments
	//
	lattice_gensys& operator = (const lattice_gensys&);
	lattice_gensys& operator = (const math_matrix< bigint > &);
	lattice_gensys& operator = (const math_matrix< bigfloat > &);

	void assign(const lattice_gensys&);
	void assign(const math_matrix< bigint > &);
	void assign(const math_matrix< bigfloat > &);

	friend void assign(lattice_gensys&, const lattice_gensys&);
	friend void assign(lattice_gensys&, const math_matrix< bigint > &);
	friend void assign(lattice_gensys&, const math_matrix< bigfloat > &);

	//
	// Element Operations
	//
	void member(lidia_size_t, lidia_size_t, double&);
	void member(lidia_size_t, lidia_size_t, bigint&);
	void member(lidia_size_t, lidia_size_t, bigfloat&);

	void sto(lidia_size_t, lidia_size_t, const double&);
	void sto(lidia_size_t, lidia_size_t, const bigint&);
	void sto(lidia_size_t, lidia_size_t, const bigfloat&);

	void member(lidia_size_t&, lidia_size_t&, double**&);
	void member(lidia_size_t&, lidia_size_t&, bigint**&);
	void member(lidia_size_t&, lidia_size_t&, bigfloat**&);

	void member(math_matrix< bigint > &);
	void member(math_matrix< bigfloat > &);

	void sto(lidia_size_t, lidia_size_t, const double**);
	void sto(lidia_size_t, lidia_size_t, const bigint**);
	void sto(lidia_size_t, lidia_size_t, const bigfloat**);

	lidia_size_t get_no_of_rows();
	lidia_size_t get_no_of_columns();

	//
	// Type Checking
	//
	bool check_double();
	bool check_bigint();

	//
	// Modes Operations
	//
	void set_computing_mode(sdigit);
	sdigit get_computing_mode();

	void set_io_mode(sdigit);
	sdigit get_io_mode();

	//
	// Quality of Reduction
	//
	void set_precision(sdigit) { }
	sdigit get_computed_precision()
	{
		return 0;
	}
	sdigit get_read_precision()
	{
		return 0;
	}

	void set_reduction_parameter(double);
	void set_reduction_parameter(sdigit, sdigit);
	double get_reduction_parameter();
	void get_reduction_parameter(sdigit&, sdigit&);

	//
	// algorithm information
	//
	sdigit get_no_of_reduction_steps();
	sdigit get_no_of_swaps();
	sdigit get_no_of_corrections();

	//
	// Time needed
	//
	timer& get_computing_time();

	//
	// Algorithms
	//

	//
	// Buchmann - Kessler
	//
	void lin_gen_system(lattice_gensys&, lidia_size_t&);
        friend lattice_gensys lin_gen_system(const lattice_gensys& L,
					     lidia_size_t& rank);

	void lll_gensys(lidia_size_t&);
	void lll_gensys(const lattice_gensys&, lidia_size_t&);
        friend lattice_gensys lll_gensys(const lattice_gensys& L, lidia_size_t& rank);

	void lll_trans_gensys(lattice_gensys&, lidia_size_t&);
	friend lattice_gensys lll_trans_gensys(const lattice_gensys& L,
					       lattice_gensys& T,
					       lidia_size_t& rank);

	void mlll(base_vector< bigint > &);
	void mlll(const lattice_gensys&, base_vector< bigint > &);
        friend lattice_gensys mlll(const lattice_gensys& L, base_vector< bigint > & v);

	void mlll(bigint*&);
	void mlll(const lattice_gensys&, bigint*&);
        friend lattice_gensys mlll(const lattice_gensys& L, bigint*& v);

	void compute_basis(const lattice_gensys&, const lattice_gensys&);

	//
	// friend functions
	//
	friend void mlll(math_matrix< bigint > &, const math_matrix< bigint > &, base_vector< bigint > &);
	friend void mlll(math_matrix< bigfloat > &, const math_matrix< bigfloat > &, base_vector< bigint > &);
	friend void mlll(math_matrix< bigfloat > &, const math_matrix< bigfloat > &, base_vector< bigfloat > &);

	friend void mlll(math_matrix< bigint > &, const math_matrix< bigint > &, bigint*&);
	friend void mlll(math_matrix< bigfloat > &, const math_matrix< bigfloat > &, bigint*&);
	friend void mlll(math_matrix< bigfloat > &, const math_matrix< bigfloat > &, bigfloat*&);

	friend void lll_gensys(math_matrix< bigfloat > &, const math_matrix< bigfloat > &, lidia_size_t&);
	friend void lll_gensys(math_matrix< bigint > &, const math_matrix< bigint > &, lidia_size_t&);

	friend void lll_trans_gensys(math_matrix< bigint > &, math_matrix< bigint > &, lidia_size_t&);
	friend void lll_trans_gensys(math_matrix< bigint > &, math_matrix< bigfloat > &, lidia_size_t&);
	friend void lll_trans_gensys(math_matrix< bigfloat > &, math_matrix< bigfloat > &, lidia_size_t&);

	friend void lin_gen_system(math_matrix< bigint > &, const math_matrix< bigint > &, lidia_size_t&);
	friend void lin_gen_system(math_matrix< bigint > &, const math_matrix< bigfloat > &, lidia_size_t&);
	friend void lin_gen_system(math_matrix< bigfloat > &, const math_matrix< bigfloat > &, lidia_size_t&);


};

// friend functions

//
// Input / Output
//
std::istream& operator >> (std::istream&, lattice_gensys&);
std::ostream& operator << (std::ostream&, const lattice_gensys&);

//
// Assignments
//

void assign(lattice_gensys&, const lattice_gensys&);
void assign(lattice_gensys&, const math_matrix< bigint > &);
void assign(lattice_gensys&, const math_matrix< bigfloat > &);

//
// Algorithms
//

//
// Buchmann - Kessler
//
inline
lattice_gensys lin_gen_system(const lattice_gensys& L,
			      lidia_size_t& rank) {
    lattice_gensys LL(L);
    lattice_gensys T(1, 1);
    LL.lin_gen_system(T, rank);
    return(T);
}

inline
lattice_gensys lll_gensys(const lattice_gensys& L, lidia_size_t& rank) {
    lattice_gensys TL(L);
    TL.lll_gensys(rank);
    return (TL);
}

inline
lattice_gensys lll_trans_gensys(const lattice_gensys& L,
				lattice_gensys& T,
				lidia_size_t& rank) {
    lattice_gensys LT(L);
    LT.lll_gensys(T, rank);
    return (LT);
}

inline
lattice_gensys mlll(const lattice_gensys& L,
		    base_vector< bigint > & v) {
    lattice_gensys TL(L);
    TL.mlll(v);
    return(TL);
}

inline
lattice_gensys mlll(const lattice_gensys& L, bigint*& v) {
    lattice_gensys TL(L);
    TL.mlll(v);
    return(TL);
}

void mlll(math_matrix< bigint > &, const math_matrix< bigint > &,
	  base_vector< bigint > &);
void mlll(math_matrix< bigfloat > &, const math_matrix< bigfloat > &,
	  base_vector< bigint > &);
void mlll(math_matrix< bigfloat > &, const math_matrix< bigfloat > &,
	  base_vector< bigfloat > &);
void mlll(math_matrix< bigint > &, const math_matrix< bigint > &,
	  bigint*&);
void mlll(math_matrix< bigfloat > &, const math_matrix< bigfloat > &,
	  bigint*&);
void mlll(math_matrix< bigfloat > &, const math_matrix< bigfloat > &,
	  bigfloat*&);

void lll_gensys(math_matrix< bigfloat > &,
		const math_matrix< bigfloat > &, lidia_size_t&);
void lll_gensys(math_matrix< bigint > &, const math_matrix< bigint > &,
		lidia_size_t&);

void lll_trans_gensys(math_matrix< bigint > &, math_matrix< bigint > &,
		      lidia_size_t&);
void lll_trans_gensys(math_matrix< bigint > &,
		      math_matrix< bigfloat > &, lidia_size_t&);
void lll_trans_gensys(math_matrix< bigfloat > &,
		      math_matrix< bigfloat > &, lidia_size_t&);

void lin_gen_system(math_matrix< bigint > &,
		    const math_matrix< bigint > &, lidia_size_t&);
void lin_gen_system(math_matrix< bigint > &,
		    const math_matrix< bigfloat > &, lidia_size_t&);
void lin_gen_system(math_matrix< bigfloat > &,
		    const math_matrix< bigfloat > &, lidia_size_t&);


#endif  // HEADBANGER



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#endif	// LIDIA_LATTICE_GENSYS_H_GUARD_
