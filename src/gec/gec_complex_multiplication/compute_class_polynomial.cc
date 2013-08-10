// -*- C++ -*-
//==============================================================================================
//
//	This file is part of LiDIA --- a library for computational number theory
//
//	Copyright (c) 1994--2002 the LiDIA Group.  All rights reserved.
//
//	See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//
//  File    : compute_class_polynomial.cc
//  Author  : Harald Baier (HB)
//	Changes	: See CVS log
//
//==============================================================================================

#ifdef HAVE_CONFIG_H
# include       "config.h"
#endif
# include        "LiDIA/gec_complex_multiplication.h"
#include         "LiDIA/modular_functions.h"
#ifdef VERBOSE
# include        <cassert>
#endif


#ifdef LIDIA_NAMESPACE
namespace LiDIA {
#endif

    namespace {

//
// fundamental functions:
// * zeta(n) returns the primitve n-th root of unity e^( 2\pi*i / n )
// * tau(Q) returns the imaginary quadratic integer corresponding
//   to the quadratic form Q
//
	bigcomplex 
	zeta( int n )
	{
	    return( exp( Pi() * 2 * bigcomplex(0,1) / n ) );
	}

	bigcomplex
	tau( const quadratic_form & Q )
	{
	    // Q = (a,b,c)  ==>  tau_Q = ( -b + i*sqrt(|\Delta|) / (2a)
	    bigcomplex tau_Q( -Q.get_b(), sqrt( bigfloat( abs(Q.discriminant()) ) ) );
	    tau_Q /= 2 * Q.get_a();
	
	    return( tau_Q );
	}

// 
// gamma_2_star function as introduced by [Atkin/Morain] or [Lay/Zimmer]
//
	bigcomplex 
	gamma_2_star( const quadratic_form & Q )
	{
	    bigint residue = 
		(Q.get_b()*(Q.get_a()-Q.get_c()+Q.get_a()*Q.get_a()*Q.get_c())) % 3;
  
	    return( power( zeta( 3 ), residue) * gamma_2( tau(Q), true ) );
	}
  
//
// function weber_f_star as introduced by [Yui/Zagier]:
// weber_f_star is a class invariant provided 
// delta = 1 mod 8, delta not divisible by 3
// see PhD, Theorem 2.2.19
//
	bigcomplex
	weber_f_star( const quadratic_form & Q )
	{
	    // remark that since in this case delta = 1 mod 8,
	    // we have either 2|a or 2|c

	    // sign stores a possible factor -1, depending on the parity of
	    //  (\Delta-1)/8
	    int sign = 1; 
	    if ( ( ( Q.discriminant() - 1 ) / 8 ).is_odd() )
		sign = -1;

	    bigint exponent, residue; // residue of exponent modulo 48
	    bigcomplex result; // result stores the result

	    if( (Q.get_a() ).is_odd() ) // it follows that c is even
	    {
		exponent = 
		    Q.get_b()*(Q.get_a()-Q.get_c()+Q.get_a()*Q.get_a()*Q.get_c());
		residue = exponent % 48;

		result =  sign * power(zeta(48),residue) * weber_f2(tau(Q), true );

		return( result );
	    }
  
	    if( Q.get_c().is_odd() ) // it follows that a is even
	    {
		exponent = 
		    Q.get_b()*(Q.get_a()-Q.get_c()-Q.get_a()*Q.get_c()*Q.get_c());
		residue = exponent % 48;

		result = 
		    sign * power(zeta(48),residue) * weber_f1(tau(Q), true );

		return( result );
	    }

	    exponent = 
		Q.get_b()*(Q.get_a()-Q.get_c()-Q.get_a()*Q.get_c()*Q.get_c());
	    residue = exponent % 48;

	    result = power(zeta(48),residue) * weber_f(tau(Q), true );

	    return( result );
	}

    } // anonymous namespace

//
// the function class_invariant( Q ) returns an approximation
// of the complex number corresponding to the quadratic form Q
// and the generation_mode
//
    bigcomplex 
    gec_complex_multiplication::class_invariant( const quadratic_form & Q )
    {
	if ( generation_mode == 4 )
	    return( weber_f_star( Q ) );
  
	if ( generation_mode == 3 )
	    return( gamma_2_star( Q ) );

	return( modular_j( tau( Q ), true ) );
    }

//
// the standard algorithm for multiplying two monic polynomials
//
    void
    gec_complex_multiplication::
    product_real( base_vector<bigfloat> & c,
		  const base_vector<bigfloat> & a,
		  const base_vector<bigfloat> & b )
    {
	lidia_size_t deg_a = a.size() - 1;
	lidia_size_t deg_b = b.size() - 1;
	lidia_size_t i, j, kk, deg_ab = deg_a + deg_b;

	bigfloat temp, temp_a;
	base_vector<bigfloat> c_temp( deg_ab + 1, FIXED );

	for( i = 0; i <= deg_ab; i++ )
	    c_temp[ i ].assign_zero();

	for( i = 0; i <= deg_a; i++ )
	{
	    temp_a.assign( a[i] );
	    for( j = 0, kk = i; j <= deg_b; j++, kk++ )
	    {
		multiply( temp, temp_a, b[j] );
		add( c_temp[kk], c_temp[kk], temp);
	    }
	}
	
	c.set_capacity( deg_ab + 1 );
	for( i = 0; i <= deg_ab; i++ )
	    swap( c[i], c_temp[i] ) ;
    }

//
// the Karatsuba variant for multiplying two monic polynomials
// of arbitrary degree (see polynomialHybrid of PhD and
// algorithm productPolynomialKaratsuba (8.2) of PhD)
//
// It is assumed that degree(f) >= degree(g). Otherwise,
// the error_handler is invoked. Furthermore, it is assumed result != f.
//
    void
    gec_complex_multiplication::
    product_karatsuba( base_vector<bigfloat> & result,
		       const base_vector<bigfloat> & f,
		       const base_vector<bigfloat> & g)
    {
	// assign the degrees
	lidia_size_t k = f.size() - 1;
	lidia_size_t l = g.size() - 1;
	lidia_size_t diff = k - l;

	// if degree(f) < degree(g), abort
	if( diff < 0 )
	    lidia_error_handler( "gec_complex_multiplication",
				 "product_karatsuba: degree(f) < degree(g)" );

	// to allow result = g
	base_vector<bigfloat> local_g(k + 1, FIXED);

	// initialize local_g with g and fill the rest with zeros
	for( lidia_size_t kk = 0; kk <= l; kk++ )
	    local_g[kk+diff].assign( g[kk] );
	for( lidia_size_t kk = 0; kk < diff; kk++ )
	    local_g[kk].assign_zero();

	//
	// for the rest: look at algorithm productPolynomialKaratsuba (8.2) of PhD
	//
	base_vector<bigfloat> d( k, FIXED );
	lidia_size_t kk, ll, n, j, upper_bound_j, counter_1, counter_2;
	bigfloat factor_1, factor_2, tmp, current_coefficient;

	for( kk = 0; kk < diff; kk++ )
	    d[kk].assign_zero();
	for( kk = diff; kk < k; kk++ )
	    multiply( d[ kk ], f[ kk ], local_g[ kk ] );

	for( n = diff; n < k; n++ )
	{
	    upper_bound_j = ( n - 1 ) / 2;
	    current_coefficient.assign_zero();
		
	    for( j = 0; j <= upper_bound_j; j++)
	    {
		add( factor_1, f[ j ], f[ n - j ] );
		add( factor_2, local_g[ j ], local_g[ n - j ] );
		multiply( tmp, factor_1, factor_2 );

		subtract( tmp, tmp, d[ j ] );
		subtract( tmp, tmp, d[ n - j ] );
			
		add( current_coefficient, current_coefficient, tmp );
	    }

	    if( n % 2 == 0 )
		add( result[ n - diff ], current_coefficient, d[ n/2 ] );
	    else
		result[n-diff].assign( current_coefficient );
	}

	for( n = k; n <= 2*k - 3; n++ )
	{
	    upper_bound_j = ( n - 1 ) / 2;
	    current_coefficient.assign_zero();
		
	    for( j = n - k + 1; j <= upper_bound_j; j++)
	    {
		add( factor_1, f[ j ], f[ n - j ] );
		add( factor_2, local_g[ j ], local_g[ n - j ] );
		multiply( tmp, factor_1, factor_2 );

		subtract( tmp, tmp, d[ j ] );
		subtract( tmp, tmp, d[ n - j ] );
			
		add( current_coefficient, current_coefficient, tmp );
	    }

	    if( n % 2 == 0 )
		add( current_coefficient, current_coefficient, d[ n/2 ] );

	    add( current_coefficient, current_coefficient, f[ n - k ] );
	    add( current_coefficient, current_coefficient, local_g[ n - k ] );
	    result[n-diff].assign( current_coefficient );
	}

	result[2*k-2-diff] = f[ k - 2 ] + d[k-1] + local_g[ k-2 ];
	result[2*k-1-diff] = f[ k - 1 ] + local_g[ k - 1 ];
	result[2*k-diff].assign_one();
    }

//
// compute the class polynomial
//
    void
    gec_complex_multiplication::
    compute_class_polynomial()
    {
	// test if all values are set
	if( generation_mode == 0 )
	{
	    class_group.set_mode( EXPAND );
	    class_group = compute_class_group( delta ); // see chapter 6 of PhD
	    h = class_group.size();

	    set_generation_mode( 0 );
	    set_complex_precision( 0 );
	}

	// if class number is low take the standard algorithm
	if( h < 50 )
	    compute_class_polynomial_real();
	// otherwise use Karatsuba variant
	else
	{	
	    if( complex_precision < 350 )
		compute_class_polynomial_karatsuba_low_precision();
	    else
		compute_class_polynomial_karatsuba_high_precision();
	}
    }

//
// this implementation corresponds to algorithm
// computeClassPolynomialReal of PhD (Algorithm 8.4)
//
    void
    gec_complex_multiplication::compute_class_polynomial_real()
    {
	bigfloat::set_precision( complex_precision );

	//
	// polynomial minimal_polynomial_approximation stores 
	// the real approximation of the class polynomial
	//
	polynomial<bigfloat> minimal_polynomial_approximation;
	minimal_polynomial_approximation.assign_one();

	// linear_polynomial stores a real root of the class polynomial
	polynomial<bigfloat> linear_polynomial;
	linear_polynomial.assign_x();

	// quadratic_polynomial stores the minimal polynomial
	// of a non-real root of the class polynomial
	polynomial<bigfloat> quadratic_polynomial;
	quadratic_polynomial.set_degree( 2 );
	quadratic_polynomial[ 2 ].assign( bigfloat( 1 ) );

	// to store the current class invariant and some real linear coefficient
	bigcomplex current_class_invariant;
	bigfloat linear_coefficient;

	for( lidia_size_t i = 0; i < h; i++ )
	{
	    // first case: c is not real ==> pol = X^2 - 2*Re(c) X + |c|^2
	    if ( (i < h - 1) && (class_group[i] == - class_group[ i + 1 ]) )
	    {
		// compute the class invariant c
		current_class_invariant.assign( class_invariant(class_group[i]) );
			
		// compute the linear coefficient  -2*Re(c)
		linear_coefficient.assign( current_class_invariant.real() );
		shift_left( linear_coefficient, linear_coefficient, 1 );
		linear_coefficient.negate();

		// assign minimal polynomial
		quadratic_polynomial[ 1 ].assign( linear_coefficient );
		quadratic_polynomial[ 0 ].assign( norm(current_class_invariant) );

		// multiply (X^2 - 2*Re(c) X + |c|^2) with accumulated polynomial
		multiply( minimal_polynomial_approximation, 
			  minimal_polynomial_approximation,
			  quadratic_polynomial );

		// the next quadratic form is obsolete in this case
		i++;
	    }
	    // otherwise: c is real
	    else
	    {
		// compute the class invariant c
		current_class_invariant.assign( class_invariant(class_group[i]) );

		// compute the linear polynomial X - c
		linear_polynomial[0].assign( current_class_invariant.real() );
		linear_polynomial[0].negate();

		// multiply (X - c) with accumulated polynomial
		multiply( minimal_polynomial_approximation, 
			  minimal_polynomial_approximation,
			  linear_polynomial );
	    }
	}

	// 
	// assign the minimal polynomial by rounding the coefficients
	// of minimal_polynomial_approximation
	//
	class_polynomial.set_degree( h );
	class_polynomial[ h ] = 1; // minimal polynomial is monic

	for( lidia_size_t i = 0; i < h; i++ )
	    minimal_polynomial_approximation[i].bigintify( class_polynomial[i] );
    }


//
// this implementation corresponds to algorithm
// computeClassPolynomialKaratsuba of PhD (Algorithm 8.5):
// we only make use of it if the precision is high, that is
// if the precision is at least 350.
//
    void
    gec_complex_multiplication::
    compute_class_polynomial_karatsuba_high_precision()
    {
	bigint count_mult( 0 );
	lidia_size_t degree_of_pol_1, degree_of_pol_2;

	bigfloat::set_precision( complex_precision );

	//
	// bigfloat d_1 stores absolute term of linear factor of an ideal class
	// of order at most 2
	//
	bigfloat d_1;
	d_1.assign_zero();
	
	// to store the coefficients of the polynomials of degree 2
	base_vector<bigfloat> degree_2( h, FIXED );
	lidia_size_t counter_degree_2 = 0;

	//
	// polynomial minimal_polynomial_approximation stores 
	// the real approximation of the class polynomial
	//	
	base_vector<bigfloat> minimal_polynomial_approximation(0, EXPAND);
	base_vector<bigfloat> tmp_pol(0, EXPAND);

	bigcomplex current_class_invariant; // current class invariant
	bigfloat current_coefficient; // current real coefficient

	//
	// compute the class invariants and polynomials of degree 2
	// lines 5 to 15 of PhD-algorithm
	for( lidia_size_t i = 0; i < h; i++ )
	{
	    // the class invariant c is not real
	    if ( (i < h - 1 ) && ( class_group[i] == - class_group[ i + 1 ] ) )
	    {
		// compute class invariant c
		current_class_invariant.assign( class_invariant(class_group[i]) );

		// assign |c|^2 to current stack entry
		degree_2[ counter_degree_2++ ].assign( 
		    norm( current_class_invariant ) );

		count_mult += 2;

		// compute  -2*Re(c)
		current_coefficient.assign( current_class_invariant.real() );
		shift_left( current_coefficient, 
			    current_coefficient, 1 );
		current_coefficient.negate();
		// assign  -2*Re(c)  to current stack entry
		degree_2[ counter_degree_2++ ].assign( current_coefficient );

		// the next quadratic form is obsolete in this case
		i++;
	    }
	    // the class invariant c is real
	    else
	    {
		// compute class invariant c
		current_class_invariant.assign( class_invariant(class_group[i]) );

		current_coefficient.assign( current_class_invariant.real() );
		current_coefficient.negate();

		// if d_1 == 0 (T=1) then store c in d_1
		if( d_1.is_zero() )
		    d_1.assign( current_coefficient );
		// otherwise compute polynomial X^2 + (d_1-c)X + d_1*c
		// and store coefficients on the stack
		else
		{
		    // compute d_1*c and store it on the stack
		    multiply( degree_2[ counter_degree_2++ ],
			      current_coefficient,
			      d_1 );
				
		    inc( count_mult );

		    // compute d_1-c and store it on the stack
		    add( degree_2[ counter_degree_2++ ],
			 current_coefficient,
			 d_1 );

		    // set d_1 = 0 to store that no linear polynomial is set
		    d_1.assign_zero();
		}
	    }
	}

	// if h is odd assign the remaining class invariant to 
	// approximation of minimal polynomial
	if( h & 1 )
	{
	    minimal_polynomial_approximation[1].assign_one();
	    minimal_polynomial_approximation[0].assign( d_1 );
	}

	//
	// do the computations for degree 4
	//

	// the number of coefficients
	lidia_size_t l = counter_degree_2 / 4;
	// base_vector to store them
	base_vector<bigfloat> old_vector( 4*l, FIXED );

	// variables for storing length and position of both arrays
	lidia_size_t length_old = 0;
	lidia_size_t position_old = 0, position_old_bis = 2;
	lidia_size_t counter_old = 1, counter_old_bis = 3;
	lidia_size_t ii, jj, kk, upper_bound_jj;

	bigfloat factor_1, factor_2, tmp;

	// compute the coefficients of the polynomials of degree 4
	// and store them in the base_vector old_vector
	for( kk = 0; kk < l; kk++ )
	{
	    multiply( old_vector[ position_old ], degree_2[ position_old ], 
		      degree_2[ position_old_bis ] );	  
	    inc( count_mult );

	    multiply( d_1, degree_2[ counter_old ], 
		      degree_2[ counter_old_bis ] );
	    inc( count_mult );

	    add( factor_1, degree_2[ position_old ], 
		 degree_2[ counter_old ] );
	    add( factor_2, degree_2[ position_old_bis ], 
		 degree_2[ counter_old_bis ] );
	    multiply( tmp, factor_1, factor_2 );
	    inc( count_mult );

	    subtract( current_coefficient, tmp, d_1 );
	    subtract( old_vector[ counter_old ], current_coefficient, 
		      old_vector[ position_old ] );

	    add( current_coefficient, degree_2[ position_old ], 
		 degree_2[ position_old_bis ] );
	    add( old_vector[ position_old_bis ], current_coefficient, d_1 );

	    add( old_vector[ counter_old_bis ], degree_2[ counter_old ], 
		 degree_2[ counter_old_bis ] );

	    position_old += 4;
	    counter_old += 4;
	    position_old_bis += 4;
	    counter_old_bis += 4;
	}

	if( position_old != counter_degree_2 )
	{
	    if( minimal_polynomial_approximation.size() != 0 )
	    {
		tmp_pol[ 2 ].assign_one();
		tmp_pol[ 1 ].assign( degree_2[ counter_degree_2 - 1 ] );
		tmp_pol[ 0 ].assign( degree_2[ counter_degree_2 - 2 ] );
		  
		product_real( minimal_polynomial_approximation,
			      minimal_polynomial_approximation,
			      tmp_pol );
	    }
	    else
	    {
		minimal_polynomial_approximation[ 2 ].assign_one();
		minimal_polynomial_approximation[ 1 ].assign( degree_2[ counter_degree_2 - 1 ] );
		minimal_polynomial_approximation[ 0 ].assign( degree_2[ counter_degree_2 - 2 ] );
	    }
	}

	// free the memory of degree_2
	degree_2.kill();

	//
	// do the further computations (if necessary)
	// for polynomials of degree 2^(j+1), where j>=2

	length_old = position_old;

	// vector to store intermediate products
	base_vector<bigfloat> d( 0, EXPAND );
	// to store coefficients of new polynomials (vector Q in PhD)
	base_vector < bigfloat > new_vector( 0, EXPAND );

	// some more counters
	lidia_size_t counter_new, new_degree;
	lidia_size_t current_degree = 4;
	lidia_size_t count_single_addends;

	while( true )
	{
	    new_degree = 2 * current_degree;
	    l = length_old / new_degree;
	  
	    new_vector.set_capacity( l * new_degree  );

	    position_old = 0;
	    position_old_bis = current_degree;
	    counter_new = 0;

	    //
	    // set capacity of d
	    //
	    d.set_capacity( current_degree );

	    for( kk = 0; kk < l; kk++ )
	    {
		// initialize d
		for( ii = 0; ii < current_degree; ii++ )
		{
		    multiply( d[ ii ], old_vector[ position_old + ii ], 
			      old_vector[ position_old_bis + ii ] );
		    inc( count_mult );
		}

		new_vector[ counter_new++ ].assign( d[ 0 ] );

		for( ii = 1; ii < current_degree; ii++ )
		{
		    upper_bound_jj = ii / 2;

		    if( ! ( ii & 1 ) )
			current_coefficient.assign( d[ upper_bound_jj ] );
		    else
		    {
			current_coefficient.assign_zero();
			upper_bound_jj++;
		    }

		    for( jj = 0; jj < upper_bound_jj; jj++ )
		    {
			add( factor_1, old_vector[ position_old + jj ], 
			     old_vector[ position_old + ii - jj ] );
			add( factor_2, old_vector[ position_old_bis + jj ], 
			     old_vector[ position_old_bis + ii - jj ] );
			multiply( tmp, factor_1, factor_2 );
			inc( count_mult );
			subtract( tmp, tmp, d[ jj ] );
			subtract( tmp, tmp, d[ ii - jj ] );

			add( current_coefficient, current_coefficient, tmp );
		    }

		    new_vector[ counter_new++ ].assign( current_coefficient );
		}

		count_single_addends = 0;

		for( ii = current_degree; ii <= new_degree - 3; ii++ )
		{
		    upper_bound_jj = ii / 2;

		    if( ! ( ii & 1 ) )
			current_coefficient.assign( d[ upper_bound_jj ] );
		    else
		    {
			current_coefficient.assign_zero();
			upper_bound_jj++;
		    }

		    for( jj = ii - current_degree + 1; 
			 jj < upper_bound_jj; jj++ )
		    {
			add( factor_1, old_vector[ position_old + jj ], 
			     old_vector[ position_old + ii - jj ] );
			add( factor_2, old_vector[ position_old_bis + jj ], 
			     old_vector[ position_old_bis + ii - jj ] );
			multiply( tmp, factor_1, factor_2 );
			inc( count_mult );
			subtract( tmp, tmp, d[ jj ] );
			subtract( tmp, tmp, d[ ii - jj ] );

			add( current_coefficient, current_coefficient, tmp );
		    }

		    add( current_coefficient, current_coefficient,
			 old_vector[ position_old + count_single_addends ] );
		    add( current_coefficient, current_coefficient,
			 old_vector[ position_old_bis + count_single_addends ] );

		    new_vector[ counter_new++ ].assign( current_coefficient );
		    count_single_addends++;
		}

		current_coefficient.assign_zero();
		add( current_coefficient, current_coefficient,
		     old_vector[ position_old + count_single_addends ] );
		add( current_coefficient, current_coefficient,
		     old_vector[ position_old_bis + count_single_addends ] );
		add( current_coefficient, current_coefficient, d[ ii / 2 ] );
		new_vector[ counter_new++ ].assign( current_coefficient );

		count_single_addends++;
		current_coefficient.assign_zero();
		add( current_coefficient, current_coefficient,
		     old_vector[ position_old + count_single_addends ] );
		add( current_coefficient, current_coefficient,
		     old_vector[ position_old_bis + count_single_addends ] );
		new_vector[ counter_new++ ].assign( current_coefficient );
		  
		position_old += new_degree;
		position_old_bis += new_degree;
	    }

	    if( counter_new < length_old )
	    {
		if( minimal_polynomial_approximation.size() != 0 )
		{
		    tmp_pol[ current_degree ].assign_one();
		    for( kk = 0; kk < current_degree; kk++ )
			tmp_pol[ kk ].assign( old_vector[ counter_new++ ] );
			  
		    if( current_degree < 
			( minimal_polynomial_approximation.size() - 1 ) / 2 )
			product_real( minimal_polynomial_approximation,
				      minimal_polynomial_approximation,
				      tmp_pol );
		    else
			product_karatsuba( minimal_polynomial_approximation,
					   tmp_pol,
					   minimal_polynomial_approximation );
		}
		else
		{
		    minimal_polynomial_approximation[ current_degree ].assign_one();
		    for( kk = 0; kk < current_degree; kk++ )
			minimal_polynomial_approximation[kk].assign( old_vector[ counter_new++ ] );
		}
	    }

	    current_degree = new_degree;

	    if( current_degree > h/2 )
		break;

	    old_vector.assign( new_vector );
	    length_old = old_vector.size();
	}

	if( h > current_degree )
	{
	    tmp_pol[ current_degree ].assign_one();
	    for( kk = 0; kk < current_degree; kk++ )
		tmp_pol[ kk ].assign( new_vector[ kk ] );
	  
	    if( current_degree < 
		( minimal_polynomial_approximation.size() - 1 ) / 2 )
		product_real( minimal_polynomial_approximation,
			      minimal_polynomial_approximation,
			      tmp_pol );
	    else
		product_karatsuba( minimal_polynomial_approximation,
				   tmp_pol, minimal_polynomial_approximation );
	}
	else
	    for( kk = 0; kk < h; kk++ )
		minimal_polynomial_approximation[ kk ].assign( new_vector[ kk ] );

	// 
	// assign the minimal polynomial by rounding the coefficients
	// of minimal_polynomial_approximation
	//
	class_polynomial.set_degree( h );
	class_polynomial[ h ] = 1; // minimal polynomial is monic

	for( kk = 0; kk < h; kk++ )
	    minimal_polynomial_approximation[kk].bigintify( class_polynomial[kk] );
    }



//
// this implementation corresponds to algorithm
// computeClassPolynomialKaratsuba of PhD (Algorithm 8.5)
// we only make use of this algorithm in case of low precision,
// that is if the precision is at most 349. The difference
// to the high precision case is simply that product_hybrid
// from PhD is ignored and always product_real is invoked.
//
    void
    gec_complex_multiplication::
    compute_class_polynomial_karatsuba_low_precision()
    {
	bigint count_mult( 0 );
	lidia_size_t degree_of_pol_1, degree_of_pol_2;

	bigfloat::set_precision( complex_precision );

	//
	// bigfloat d_1 stores absolute term of linear factor of an ideal class
	// of order at most 2
	//
	bigfloat d_1;
	d_1.assign_zero();
	
	// to store the coefficients of the polynomials of degree 2
	base_vector<bigfloat> degree_2( h, FIXED );
	lidia_size_t counter_degree_2 = 0;

	//
	// polynomial minimal_polynomial_approximation stores 
	// the real approximation of the class polynomial
	//	
	polynomial<bigfloat> minimal_polynomial_approximation;
	minimal_polynomial_approximation.assign_one();

	bigcomplex current_class_invariant; // current class invariant
	bigfloat current_coefficient; // current real coefficient
	
	//
	// compute the class invariants and polynomials of degree 2
	// lines 5 to 15 of PhD-algorithm
	for( lidia_size_t i = 0; i < h; i++ )
	{
	    // the class invariant c is not real
	    if ( (i < h - 1 ) && ( class_group[i] == - class_group[ i + 1 ] ) )
	    {
		// compute class invariant c
		current_class_invariant.assign( class_invariant(class_group[i]) );

		// assign |c|^2 to current stack entry
		degree_2[ counter_degree_2++ ].assign( 
		    norm( current_class_invariant ) );

		count_mult += 2;

		// compute  -2*Re(c)
		current_coefficient.assign( current_class_invariant.real() );
		shift_left( current_coefficient, 
			    current_coefficient, 1 );
		current_coefficient.negate();
		// assign  -2*Re(c)  to current stack entry
		degree_2[ counter_degree_2++ ].assign( current_coefficient );

		// the next quadratic form is obsolete in this case
		i++;
	    }
	    // the class invariant c is real
	    else
	    {
		// compute class invariant c
		current_class_invariant.assign( class_invariant(class_group[i]) );

		current_coefficient.assign( current_class_invariant.real() );
		current_coefficient.negate();

		// if d_1 == 0 (T=1) then store c in d_1
		if( d_1.is_zero() )
		    d_1.assign( current_coefficient );
		// otherwise compute polynomial X^2 + (d_1-c)X + d_1*c
		// and store coefficients on the stack
		else
		{
		    // compute d_1*c and store it on the stack
		    multiply( degree_2[ counter_degree_2++ ],
			      current_coefficient,
			      d_1 );
				
		    inc( count_mult );

		    // compute d_1-c and store it on the stack
		    add( degree_2[ counter_degree_2++ ],
			 current_coefficient,
			 d_1 );

		    // set d_1 = 0 to store that no linear polynomial is set
		    d_1.assign_zero();
		}
	    }
	}

	// if h is odd assign the remaining class invariant to 
	// approximation of minimal polynomial
	if( h & 1 )
	{
	    minimal_polynomial_approximation.assign_x();
	    minimal_polynomial_approximation[ 0 ].assign( d_1 );
	}

	//
	// do the computations for degree 4
	//

	// the number of coefficients
	lidia_size_t l = counter_degree_2 / 4;
	// base_vector to store them
	base_vector<bigfloat> old_vector( 4*l, FIXED );

	// variables for storing length and position of both arrays
	lidia_size_t length_old = 0;
	lidia_size_t position_old = 0, position_old_bis = 2;
	lidia_size_t counter_old = 1, counter_old_bis = 3;
	lidia_size_t ii, jj, kk, upper_bound_jj;

	bigfloat factor_1, factor_2, tmp;

	// compute the coefficients of the polynomials of degree 4
	// and store them in the base_vector old_vector
	for( kk = 0; kk < l; kk++ )
	{
	    multiply( old_vector[ position_old ], degree_2[ position_old ], 
		      degree_2[ position_old_bis ] );	  
	    inc( count_mult );

	    multiply( d_1, degree_2[ counter_old ], 
		      degree_2[ counter_old_bis ] );
	    inc( count_mult );

	    add( factor_1, degree_2[ position_old ], 
		 degree_2[ counter_old ] );
	    add( factor_2, degree_2[ position_old_bis ], 
		 degree_2[ counter_old_bis ] );
	    multiply( tmp, factor_1, factor_2 );
	    inc( count_mult );

	    subtract( current_coefficient, tmp, d_1 );
	    subtract( old_vector[ counter_old ], current_coefficient, 
		      old_vector[ position_old ] );

	    add( current_coefficient, degree_2[ position_old ], 
		 degree_2[ position_old_bis ] );
	    add( old_vector[ position_old_bis ], current_coefficient, d_1 );

	    add( old_vector[ counter_old_bis ], degree_2[ counter_old ], 
		 degree_2[ counter_old_bis ] );

	    position_old += 4;
	    counter_old += 4;
	    position_old_bis += 4;
	    counter_old_bis += 4;
	}

	polynomial<bigfloat> tmp_pol;
  
	if( position_old != counter_degree_2 )
	{
	    tmp_pol.set_degree( 2 );
	    tmp_pol[ 2 ].assign_one();
	    tmp_pol[ 1 ].assign( degree_2[ counter_degree_2 - 1 ] );
	    tmp_pol[ 0 ].assign( degree_2[ counter_degree_2 - 2 ] );

	    degree_of_pol_1 = minimal_polynomial_approximation.degree();
	    degree_of_pol_2 = 2;
	    count_mult += degree_of_pol_1 * degree_of_pol_2;

	    multiply( minimal_polynomial_approximation,
		      minimal_polynomial_approximation,
		      tmp_pol );

	}

	// free the memory of degree_2
	degree_2.kill();

	//
	// do the further computations (if necessary)
	// for polynomials of degree 2^(j+1), where j>=2

	length_old = position_old;

	// vector to store intermediate products
	base_vector<bigfloat> d( 0, EXPAND );
	// to store coefficients of new polynomials (vector Q in PhD)
	base_vector < bigfloat > new_vector( 0, EXPAND );

	// some more counters
	lidia_size_t counter_new, new_degree;
	lidia_size_t current_degree = 4;
	lidia_size_t count_single_addends;

	while( true )
	{
	    new_degree = 2 * current_degree;
	    l = length_old / new_degree;
	  
	    new_vector.set_capacity( l * new_degree  );

	    position_old = 0;
	    position_old_bis = current_degree;
	    counter_new = 0;

	    //
	    // set capacity of d
	    //
	    d.set_capacity( current_degree );

	    for( kk = 0; kk < l; kk++ )
	    {
		// initialize d
		for( ii = 0; ii < current_degree; ii++ )
		{
		    multiply( d[ ii ], old_vector[ position_old + ii ], 
			      old_vector[ position_old_bis + ii ] );
		    inc( count_mult );
		}

		new_vector[ counter_new++ ].assign( d[ 0 ] );

		for( ii = 1; ii < current_degree; ii++ )
		{
		    upper_bound_jj = ii / 2;

		    if( ! ( ii & 1 ) )
			current_coefficient.assign( d[ upper_bound_jj ] );
		    else
		    {
			current_coefficient.assign_zero();
			upper_bound_jj++;
		    }

		    for( jj = 0; jj < upper_bound_jj; jj++ )
		    {
			add( factor_1, old_vector[ position_old + jj ], 
			     old_vector[ position_old + ii - jj ] );
			add( factor_2, old_vector[ position_old_bis + jj ], 
			     old_vector[ position_old_bis + ii - jj ] );
			multiply( tmp, factor_1, factor_2 );
			inc( count_mult );
			subtract( tmp, tmp, d[ jj ] );
			subtract( tmp, tmp, d[ ii - jj ] );

			add( current_coefficient, current_coefficient, tmp );
		    }

		    new_vector[ counter_new++ ].assign( current_coefficient );
		}

		count_single_addends = 0;

		for( ii = current_degree; ii <= new_degree - 3; ii++ )
		{
		    upper_bound_jj = ii / 2;

		    if( ! ( ii & 1 ) )
			current_coefficient.assign( d[ upper_bound_jj ] );
		    else
		    {
			current_coefficient.assign_zero();
			upper_bound_jj++;
		    }

		    for( jj = ii - current_degree + 1; 
			 jj < upper_bound_jj; jj++ )
		    {
			add( factor_1, old_vector[ position_old + jj ], 
			     old_vector[ position_old + ii - jj ] );
			add( factor_2, old_vector[ position_old_bis + jj ], 
			     old_vector[ position_old_bis + ii - jj ] );
			multiply( tmp, factor_1, factor_2 );
			inc( count_mult );
			subtract( tmp, tmp, d[ jj ] );
			subtract( tmp, tmp, d[ ii - jj ] );

			add( current_coefficient, current_coefficient, tmp );
		    }

		    add( current_coefficient, current_coefficient,
			 old_vector[ position_old + count_single_addends ] );
		    add( current_coefficient, current_coefficient,
			 old_vector[ position_old_bis + count_single_addends ] );

		    new_vector[ counter_new++ ].assign( current_coefficient );
		    count_single_addends++;
		}

		current_coefficient.assign_zero();
		add( current_coefficient, current_coefficient,
		     old_vector[ position_old + count_single_addends ] );
		add( current_coefficient, current_coefficient,
		     old_vector[ position_old_bis + count_single_addends ] );
		add( current_coefficient, current_coefficient, d[ ii / 2 ] );
		new_vector[ counter_new++ ].assign( current_coefficient );

		count_single_addends++;
		current_coefficient.assign_zero();
		add( current_coefficient, current_coefficient,
		     old_vector[ position_old + count_single_addends ] );
		add( current_coefficient, current_coefficient,
		     old_vector[ position_old_bis + count_single_addends ] );
		new_vector[ counter_new++ ].assign( current_coefficient );
		  
		position_old += new_degree;
		position_old_bis += new_degree;
	    }

	    if( counter_new < length_old )
	    {
		tmp_pol.set_degree( current_degree );
		tmp_pol[ current_degree ].assign_one();
	  
		for( kk = 0; kk < current_degree; kk++ )
		    tmp_pol[ kk ].assign( old_vector[ counter_new++ ] );
			  
		degree_of_pol_1 = minimal_polynomial_approximation.degree();
		degree_of_pol_2 = current_degree;
		count_mult += degree_of_pol_1 * degree_of_pol_2;

		multiply( minimal_polynomial_approximation,
			  minimal_polynomial_approximation,
			  tmp_pol );
	    }

	    current_degree = new_degree;

	    if( current_degree > h/2 )
		break;

	    old_vector.assign( new_vector );
	    length_old = old_vector.size();
	}

	tmp_pol.set_degree( current_degree );
	tmp_pol[ current_degree ].assign_one();
	
	for( kk = 0; kk < current_degree; kk++ )
	    tmp_pol[ kk ].assign( new_vector[ kk ] );
	
	degree_of_pol_1 = minimal_polynomial_approximation.degree();
	degree_of_pol_2 = current_degree;
	count_mult += degree_of_pol_1 * degree_of_pol_2;

	multiply( minimal_polynomial_approximation,
		  minimal_polynomial_approximation,
		  tmp_pol );

	// 
	// assign the minimal polynomial by rounding the coefficients
	// of minimal_polynomial_approximation
	//
	class_polynomial.set_degree( h );
	class_polynomial[ h ] = 1; // minimal polynomial is monic

	for( kk = 0; kk < h; kk++ )
	    minimal_polynomial_approximation[kk].bigintify( class_polynomial[kk] );
    }



#ifdef LIDIA_NAMESPACE
}	// end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif
