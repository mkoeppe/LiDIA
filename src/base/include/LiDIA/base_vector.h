// -*- C++ -*-
//==============================================================================================
//
//      This file is part of LiDIA --- a library for computational number theory
//
//      Copyright (c) 1994--2001 the LiDIA Group.  All rights reserved.
//
//      See http://www.informatik.tu-darmstadt.de/TI/LiDIA/
//
//----------------------------------------------------------------------------------------------
//
//      $Id$
//
//      Author  : Frank Lehmann (FL), Markus Maurer (MM),
//                Thomas Papanikolaou (TP), Patrick Theobald (PT)
//      Changes : See CVS log
//
//==============================================================================================


#ifndef LIDIA_BASE_VECTOR_H_GUARD_
#define LIDIA_BASE_VECTOR_H_GUARD_



#ifndef LIDIA_LIDIA_H_GUARD_
# include       "LiDIA/LiDIA.h"
#endif
#ifndef LIDIA_VECTOR_REPRESENTATION_H_GUARD_
# include       "LiDIA/vector_representation.h"
#endif



#ifdef LIDIA_NAMESPACE
namespace LiDIA {
# define IN_NAMESPACE_LIDIA
#endif



class crt;
template< class T > class matrix;
template< class T > class base_matrix;
template< class T > class ring_matrix;
template< class T > class field_matrix;
template< class T > class dense_power_series;



template< class T >
class base_vector : public vector_representation< T >
{
private:

        friend class crt;
        friend class base_matrix< T >;
        friend class ring_matrix< T >;
        friend class field_matrix< T >;
        friend class matrix< T >;
        friend class dense_power_series< T >;



        //
        // constructors
        //

public:

        base_vector();
        explicit base_vector(const vector_flags & md);
        explicit base_vector(lidia_size_t all);
        base_vector(lidia_size_t, const vector_flags &);
        base_vector(lidia_size_t, lidia_size_t);
        base_vector(lidia_size_t, lidia_size_t, const vector_flags &);
        base_vector(const base_vector< T > &);
        base_vector(const base_vector< T > &, const vector_flags &);

#ifndef HEADBANGER
        base_vector(const T *, lidia_size_t);
        base_vector(const T *, lidia_size_t, const vector_flags &);
#endif

        //
        // destructor
        //

public:

        virtual ~base_vector();

        //
        // Input / Output
        //

        void write(std::ostream &) const;
        void read(std::istream &);



        //
        // BEGIN: access functions
        //

        //
        // capacity
        //

public:

        lidia_size_t capacity() const
        {
                debug_handler ("base_vector< T >", "capacity()");
                return this->allocated;
        }

        lidia_size_t get_capacity() const
        {
                debug_handler ("base_vector< T >", "get_capacity()");
                return this->allocated;
        }

        void set_capacity(lidia_size_t all);

        void reset()
        {
                debug_handler ("base_vector< T >", "reset()");
                this->kill();
        }

        void kill();

        //
        // size
        //

public:

        lidia_size_t size() const
        {
                debug_handler ("base_vector< T >", "size()");
                return this->length;
        }

        lidia_size_t get_size() const
        {
                debug_handler ("base_vector< T >", "get_size()");
                return this->length;
        }

        void set_size(lidia_size_t len);

        //
        // exp_ratio
        //

public:

        float exp_ratio() const
        {
                debug_handler ("base_vector< T >", "exp_ratio()");
                return this->bitfield.dyn_exp_ratio;
        }

        float get_exp_ratio() const
        {
                debug_handler ("base_vector< T >", "get_exp_ratio()");
                return this->bitfield.dyn_exp_ratio;
        }

        void set_exp_ratio (float ratio)
        {
                debug_handler ("base_vector< T >", "set_exp_ratio(float)");
                this->bitfield.dyn_exp_ratio = ((ratio >= (float) 1.0) ? ratio : (float) 1.0);
        }

        //
        // info_mode
        //

public:

        unsigned long get_info_mode() const
        {
                debug_handler ("base_vector< T >", "get_info_mode()");
                return this->bitfield.get_info_mode();
        }

        void set_info_mode(char md)
        {
                debug_handler ("base_vector< T >", "set_info_mode(char)");
                this->bitfield.set_info_mode(md);
        }

        unsigned long get_mode() const
        {
                debug_handler ("base_vector< T >", "get_mode()");
                return this->bitfield.get_info_mode();
        }

        void set_mode(char md)
        {
                debug_handler ("base_vector< T >", "set_mode(char)");
                this->bitfield.set_info_mode(md);
        }

        void set_mode(const vector_flags &md)
        {
                debug_handler ("base_vector< T >", "set_mode(const vector_flags&)");
                this->bitfield = md;
        }

        //
        // member
        //

public:

        T & operator[] (lidia_size_t);

        const T & operator[] (lidia_size_t i) const
        {
                debug_handler ("base_vector< T >", "operator[](lidia_size_t)");
                return this->member(i);
        }

        const T & member(lidia_size_t) const;

        //
        // value array
        //

public:

        void set_data(const T *, lidia_size_t);

        T *get_data() const;

public:

        T * get_data_address() const
        {
                debug_handler ("base_vector< T >", "get_data_address()");
                return this->value;
        }                         // returns a pointer to the

#ifndef HEADBANGER                        // internal elts of the vector
        T * set_data_address(T *d)
        {
                debug_handler ("base_vector< T >", "set_data_address(T*)");
                T *tmp = this->value;
                this->value = d;
                return tmp;
        }
#endif

        //
        // END: access functions
        //

        //
        // assignment
        //

public:

        base_vector< T > & operator = (const base_vector< T > &v)
        {
                debug_handler ("base_vector< T >", "operator = (const base_vector< T > &)");
                this->assign(v);
                return *this;
        }

        void assign(const base_vector< T > &);

        void assign(lidia_size_t, const base_vector< T > &, lidia_size_t, lidia_size_t);

        void fill(const T&);

        //
        // reverse
        //

        void reverse();
        void reverse(const base_vector< T > &);


        //
        // swap functions
        //

        void swap(base_vector< T > &);



        //
        // concat functions
        //

        void concat(const base_vector< T > &, const base_vector< T > &);



        //
        // shift functions
        //

        void shift_left(lidia_size_t, lidia_size_t);
        void shift_right(lidia_size_t, lidia_size_t);



        //
        // remove functions
        //

        void remove_from(lidia_size_t, lidia_size_t l = 1);

        //
        // insert functions
        //

        void insert_at(const T &, lidia_size_t);

};



//
// I/O
//

template< class T >
inline std::ostream &
operator << (std::ostream & s, const base_vector< T > & x)
{
        debug_handler ("base_vector< T >", "operator << (std::ostream&, "
                       "const base_vector< T > &)");
        x.write(s);
        return s;
}



template< class T >
inline std::istream &
operator >> (std::istream & s, base_vector< T > & x)
{
        debug_handler ("base_vector< T >", "operator >> (std::istream&, "
                       "base_vector< T > &)");
        x.read(s);
        return s;
}



//
// swap functions
//

template< class T >
inline void
swap(base_vector< T > & a, base_vector< T > & b)
{
        debug_handler ("base_vector< T >", "swap(base_vector< T > &, "
                       "base_vector< T > &)");
        a.swap(b);
}



#ifdef LIDIA_NAMESPACE
}       // end of namespace LiDIA
# undef IN_NAMESPACE_LIDIA
#endif



#ifdef LIDIA_INCLUDE_CC
# include       "LiDIA/base_vector.cc"
#endif



#endif  // LIDIA_BASE_VECTOR_H_GUARD_
