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
//      Author  : Patrick Theobald (PT)
//      Changes : See CVS log
//
//==============================================================================================


#ifndef LIDIA_BASE_VECTOR_CC_GUARD_
#define LIDIA_BASE_VECTOR_CC_GUARD_



#include "LiDIA/base_vector.h"
#include "LiDIA/arith.inl"
#include "LiDIA/base/array_functions.h"
#include "LiDIA/precondition_error.h"


#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
namespace LiDIA {
# endif
#endif



//
// debug defines / error defines
//

extern const char *PRT;
extern const char *vector_error_msg[];

#define DV_BV LDBL_VECTOR         // Debug value
#define DM_BV "base_vector"       // Debug message / Error message
#define LBV_ERROR vector_error_msg

//
// debug level
//
//   0 : constructors / destructor
//   1 : input / output
//   2 : access functions
//   3 : divide / concat / assign
//   4 : swap functions
//   5 : structur functions
//   6 : assignments
//   7 : reverse functions
//   8 : stream handling
//   9 : remove / insert
//

//
// constructors
//

template< class T >
base_vector< T >::base_vector ()
        : vector_representation< T > ()
{
        debug_handler_l(DM_BV,
                        "base_vector()", DV_BV);

        this->length = this->allocated = 0;
        this->value = NULL;
}



template < class T >
base_vector< T >::base_vector (const vector_flags & md)
{
        debug_handler_l(DM_BV,
                        "base_vector(vector_flags &)", DV_BV);

        this->length = this->allocated = 0;
        this->bitfield = md;
        this->value = NULL;
}



template< class T >
base_vector< T >::base_vector (lidia_size_t all)
{
        debug_handler_l(DM_BV,
                        "base_vector(lidia_size_t)", DV_BV);

        if (all < 0)
                precondition_error_handler(all, "all", "all > 0",
                                    "base_vector< T >::base_vector(lidia_size_t all)",
                                    DM_BV, LBV_ERROR[0]);
        this->length = 0;
        this->allocated = all;

        this->value = new T[all];
        memory_handler(this->value, DM_BV, "base_vector :: "
                       "out of memory");
}



template< class T >
base_vector< T >::base_vector (lidia_size_t all, const vector_flags & md)
{
        debug_handler_l(DM_BV, "in constructor "
                        "base_vector(lidia_size_t, const vector_flags &)", DV_BV);

        this->length = 0;
        this->allocated = all;
        this->bitfield = md;

        this->value = NULL;
        if (all > 0) {
                this->value = new T[all];
                memory_handler(this->value, DM_BV, "base_vector :: "
                               "out of memory");
        }
        else
                if (all < 0)
                        precondition_error_handler(all, "all", "all > 0",
                                            "base_vector< T >::"
                                            "base_vector(lidia_size_t all, const vector_flags & md)",
                                            DM_BV, LBV_ERROR[0]);
}



template< class T >
base_vector< T >::base_vector (lidia_size_t all, lidia_size_t len)
{
        debug_handler_l(DM_BV, "in constructor "
                        "base_vector(lidia_size_t, lidia_size_t)", DV_BV);

        this->length = (len < all) ? len : all;
        this->allocated = all;

        this->value = NULL;
        if (all > 0) {
                this->value = new T[all];
                memory_handler(this->value, DM_BV, "base_vector :: "
                               "out of memory");
        }
        else
                if (all < 0)
                        precondition_error_handler(all, "all", "all > 0",
                                            "base_vector< T >::"
                                            "base_vector(lidia_size_t all, lidia_size_t len)",
                                            DM_BV, LBV_ERROR[0]);
}



template< class T >
base_vector< T >::base_vector (lidia_size_t all, lidia_size_t len, const vector_flags & md)
{
        debug_handler_l(DM_BV, "in constructor "
                        "base_vector(lidia_size_t, lidia_size_t, const vector_flags &)", DV_BV);

        this->length = (len < all) ? len : all;
        this->allocated = all;

        this->bitfield = md;

        this->value = NULL;
        if (all > 0) {
                this->value = new T[all];
                memory_handler(this->value, DM_BV, "base_vector :: "
                               "out of memory");
        }
        else
                if (all < 0)
                        precondition_error_handler(all, "all", "all > 0",
                                            "base_vector< T >::"
                                            "base_vector(lidia_size_t all, lidia_size_t len, vector_flags & md)",
                                            DM_BV, LBV_ERROR[0]);
}



template< class T >
base_vector< T >::base_vector (const base_vector< T > &v)
{
        debug_handler_l(DM_BV, "in constructor "
                        "base_vector(base_vector< T > &)", DV_BV);

        this->length = v.length;
        this->allocated = v.length;

        this->bitfield = v.bitfield;

        if (this->allocated > 0) {
                this->value = new T[this->allocated];
                memory_handler(this->value, DM_BV, "base_vector :: "
                               "out of memory");

                copy_data(this->value, v.value, this->length);
        }
        else
                this->value = NULL;
}



template< class T >
base_vector< T >::base_vector (const base_vector< T > & v, const vector_flags & md)
{
        debug_handler_l(DM_BV, "in constructor "
                        "base_vector(base_vector< T > &, const vector_flags &)", DV_BV);

        this->length = v.length;
        this->allocated = v.length;

        this->bitfield = md;
        if (this->allocated > 0) {
                this->value = new T[this->allocated];
                memory_handler(this->value, DM_BV, "base_vector :: "
                               "out of memory");

                copy_data(this->value, v.value, this->length);
        }
        else
                this->value = NULL;
}



template< class T >
base_vector< T >::base_vector (const T *v, lidia_size_t len)
{
        debug_handler_l(DM_BV, "in constructor "
                        "base_vector(const T *, lidia_size_t)", DV_BV);

        this->length = len;
        this->allocated = len;

        this->value = NULL;
        if (this->allocated > 0) {
                this->value = new T[this->allocated];
                memory_handler(this->value, DM_BV, "base_vector :: "
                               "out of memory");

                LiDIA::copy_data(static_cast<T *>(this->value), static_cast<const T *>(v), this->length);
        }
        else
                if (this->allocated < 0)
                        precondition_error_handler(this->allocated, "allocated", "allocated > 0",
                                            "base_vector< T >::"
                                            "base_vector(const T *v, lidia_size_t len)",
                                            DM_BV, LBV_ERROR[6]);
}



template< class T >
base_vector< T >::base_vector (const T *v, lidia_size_t len, const vector_flags & md)
{
        debug_handler_l(DM_BV, "in constructor "
                        "base_vector(const T *, lidia_size_t, vector_flags &)", DV_BV);

        this->length = len;
        this->allocated = len;
        this->bitfield = md;

        this->value = NULL;

        if (this->allocated > 0) {
                this->value = new T[this->allocated];
                memory_handler(this->value, DM_BV, "base_vector :: "
                               "out of memory");

                LiDIA::copy_data(static_cast<T *>(this->value), static_cast<const T *>(v), this->length);
        }
        else
                if (this->allocated < 0)
                        precondition_error_handler(this->allocated, "allocated", "allocated > 0",
                                            "base_vector< T >::"
                                            "base_vector(const T *v, lidia_size_t len, vector_flags & md)",
                                            DM_BV, LBV_ERROR[6]);
}



//
// destructor
//

template< class T >
base_vector< T >::~base_vector()
{
        debug_handler_l(DM_BV, "in destructor "
                        "~base_vector()", DV_BV);

        if (this->allocated)
                delete[] this->value;
}



//
//  Input / Output
//

template< class T >
void
base_vector< T >::read (std::istream & in)
{
        debug_handler_l(DM_BV, "in member - function "
                        "read(std::istream &)", DV_BV + 1);

        register lidia_size_t n = 0;
        lidia_size_t sz = 4;

        register lidia_size_t i;

        if (in.flags() & std::ios::binary) {
                // binary-mode for stream 'in'
                in >> sz;

                set_capacity(sz);
                this->length = sz;

                for (i = 0; i < this->length; i++)
                        in >> this->value[i];
                return;
        }


        char c;

        T *vbuf = new T[sz];
        memory_handler(vbuf, DM_BV, "read :: "
                       "out of memory");
        in >> std::ws >> c;
        if (c != '[')
                lidia_error_handler(DM_BV, "read :: "
                                    "[ expected");
        in >> std::ws >> c;
        while (c != ']') {
                in.putback(c);
                in >> vbuf[n];
                n++;

                if (n == sz) {
                        debug_handler_l(DM_BV, "read :: "
                                        "doubling input size", DV_BV);

                        lidia_size_t osz = sz;

                        sz *= 2;

                        T *nbuf = new T[sz];
                        memory_handler(nbuf, DM_BV, "read :: "
                                       "out of memory");

                        copy_data(nbuf, vbuf, osz);
                        delete[] vbuf;
                        vbuf = nbuf;
                }
                in >> std::ws >> c;
        }

        set_capacity(n);
        this->length = n;
        copy_data(this->value, vbuf, n);
        delete[] vbuf;
}



template< class T >
void
base_vector< T >::write (std::ostream & out) const
{
        debug_handler_l (DM_BV, "in member - function "
                         "write(std::ostream &)", DV_BV + 1);

        register lidia_size_t i;

        if (out.flags() & std::ios::binary) {
                // binary-mode for stream 'out'
                out << this->length;

                for (i = 0; i < this->length; i++)
                        out << this->value[i] << " " << std::flush;
                return;
        }

        out << "[ ";

        for (i = 0; i < this->length; i++)
                out << this->value[i] << " ";

        out << "]" << std::flush;
}



//
// structur functions: reading and modifying the variables of the vector
//

template< class T >
void
base_vector< T >::set_capacity (lidia_size_t all)
{
        debug_handler_l(DM_BV, "in member - function "
                        "set_capacity(lidia_size_t)", DV_BV + 5);

        if (all != this->allocated)
                if (all == 0) {
                        delete[] this->value;
                        this->value = NULL;

                        this->length = 0;
                        this->allocated = 0;
                }
                else
                        if (all > 0) {
                                T *tmp = new T[all];
                                memory_handler(tmp, DM_BV , "set_capacity :: "
                                               "out of memory");

                                this->length = (all < this->length) ? all : this->length;
                                copy_data(tmp, this->value, this->length);

                                if (this->value != NULL)
                                        delete[] this->value;

                                this->value = tmp;
                                this->allocated = all;
                        }
                        else
                                precondition_error_handler(all, "all", "all > 0",
                                                    "void base_vector< T >::"
                                                    "set_capacity(lidia_size_t all)",
                                                    DM_BV, LBV_ERROR[0]);
}



template< class T >
void
base_vector< T >::kill ()
{
        debug_handler_l(DM_BV, "in member - function "
                        "kill()", DV_BV + 5);

        if (this->allocated) {
                delete[] this->value;
                this->value = NULL;

                this->length = 0;
                this->allocated = 0;
        }
}



template< class T >
void
base_vector< T >::set_size (lidia_size_t len)
{
        debug_handler_l(DM_BV, "in member - function "
                        "set_size(lidia_size_t)", DV_BV + 5);

        if (len < 0)
                precondition_error_handler(len, "len", "len > 0",
                                    "void base_vector< T >::"
                                    "set_size(lidia_size_t len)",
                                    DM_BV, LBV_ERROR[0]);
        else
                if (len <= this->allocated)
                        this->length = len;
                else
                        if (this->bitfield.info_mode == vector_flags::expand) {
                                //  realloc 'value' to provide access to expanding vectors
                                set_capacity(static_cast<lidia_size_t>(this->bitfield.dyn_exp_ratio * len));
                                this->length = len;
                        }
                        else
                                precondition_error_handler(this->bitfield.info_mode, "info_mode", "mode == vector_flags::expand",
                                                    "void base_vector< T >::"
                                                    "set_size(lidia_size_t len)",
                                                    DM_BV, LBV_ERROR[7]);
}



//
// access functions
//

template< class T >
T &
base_vector< T >::operator [](lidia_size_t i)
{
        debug_handler_l(DM_BV, "in operator "
                        "[] (lidia_size_t)", DV_BV + 2);

        if (i < 0)
                precondition_error_handler(i, "i", "i >= 0",
                                    "T & base_vector< T >::"
                                    "operator [](lidia_size_t i)",
                                    DM_BV, LBV_ERROR[3]);
        else
                if (i < this->allocated)
                        this->length = (i + 1 > this->length) ? (i + 1) : this->length;
                else
                        if (this->bitfield.info_mode == vector_flags::expand) {
                                set_capacity(static_cast<lidia_size_t>(this->bitfield.dyn_exp_ratio * (i + 1)));
                                this->length = i + 1;
                        }
                        else
                                precondition_error_handler(i, "i", "i < allocated",
                                                    this->bitfield.info_mode, "info_mode", "info_mode == vector_flags::expand",
                                                    "T & base_vector< T >::"
                                                    "operator [](lidia_size_t i)",
                                                    DM_BV, LBV_ERROR[3]);
        return this->value[i];
}



#if 0
template < class T >
const T &
base_vector< T >::operator [] (lidia_size_t i) const
{
        debug_handler_l(DM_BV, "in operator "
                        "[] (lidia_size_t) const", DV_BV + 2);

        if (i< 0 || i >= this->length)
                precondition_error_handler(i, "i", "0 <= i < length",
                                    "const T & base_vector< T >::"
                                    "operator [](lidia_size_t i) const",
                                    DM_BV, LBV_ERROR[3]);
        return this->value[i];
}
#endif



template< class T >
const T &
base_vector< T >::member(lidia_size_t i) const
{
        debug_handler_l(DM_BV, "in meber - function "
                        "member(lidia_size_t)", DV_BV + 2);

        if (i< 0 || i >= this->length)
                precondition_error_handler(i, "i", "0 <= i < length",
                                    "T & base_vector< T >::"
                                    "member(lidia_size_t i) const",
                                    DM_BV, LBV_ERROR[3]);
        return this->value[i];
}



template< class T >
void
base_vector< T >::set_data (const T *d, lidia_size_t l)
{
        debug_handler_l(DM_BV, "in member - function "
                        "set_data(const T *, lidia_size_t)", DV_BV + 2);

        set_capacity(l);
        this->length = l;

        register lidia_size_t i;
        for (i = 0; i < l; i++)
                this->value[i] = d[i];
}



template< class T >
T *
base_vector< T >::get_data () const
{
        debug_handler_l(DM_BV, "in base_vector :: "
                        "get_data()", DV_BV + 2);

        T *d = new T[this->length];
        memory_handler(d, DM_BV, "get_data :: "
                       "memory exhausted");

        register lidia_size_t i;
        for (i = 0; i < this->length; i++)
                d[i] = this->value[i];
        return d;
}



//
// assignment
//

#if 0
template < class T >
base_vector< T > &
base_vector< T >::operator = (const base_vector< T > & v)
{
        debug_handler_l(DM_BV, "in operator "
                        " = (const base_vector< T > &)", DV_BV + 6);

        assign(v);
        return *this;
}
#endif



template< class T >
void
base_vector< T >::assign (lidia_size_t at, const base_vector< T > & v, lidia_size_t from, lidia_size_t to)
{
        debug_handler_l(DM_BV, "in member - function "
                        "assign(lidia_size_t, base_vector, lidia_size_t, lidia_size_t)", DV_BV + 3);

        register lidia_size_t i, j;

        if (at< 0 || from < 0 || from >= v.length || to< 0 || to >= v.length || to < from)
                lidia_error_handler(DM_BV, "assign :: "
                                    "invalid indices");

        register lidia_size_t n = to - from + 1; // this value is positive
        register lidia_size_t new_len = (this->length < at + n) ? (at + n) : this->length;

        if (this != &v) {
                // v aliases current instance
                set_size(new_len);
                for (i = from, j = at; i <= to; i++, j++)
                        this->value[j] = v.value[i];
        }
        else {
                set_size(new_len);
                T *tmp = new T[n];
                memory_handler(tmp, DM_BV, "assign :: "
                               "out of memory");

                for (i = from, j = 0; i <= to; i++, j++)
                        tmp[j] = v.value[i];

                for (i = 0, j = at; i < n; i++, j++)
                        this->value[j] = tmp[i];

                delete[] tmp;
        }
}



template< class T >
void
base_vector< T >::assign (const base_vector< T > &v)
{
        debug_handler_l(DM_BV, "in member - function "
                        "assign(const base_vector< T > &)", DV_BV + 6);

        if (this != &v) {
                set_capacity(v.length);
                this->length = v.length;

                copy_data(this->value, v.value, this->length);

                this->sort_dir = v.sort_dir;
                this->el_cmp = v.el_cmp;
        }
}


template< class T >
void
base_vector< T >::fill(const T& x) 
{
        for(lidia_size_t i = 0; i < this->length; ++i)
        {
                this->value[i] = x;
        }
}


//
// reverse function
//

template< class T >
void
base_vector< T >::reverse ()
{
        debug_handler_l(DM_BV, "in member - function "
                        "reverse()", DV_BV + 7);

        register lidia_size_t i, j;

        for (i = 0, j = this->length - 1; i < j; i++, j--)
                LiDIA::swap(this->value[i], this->value[j]);

}



template< class T >
void
base_vector< T >::reverse (const base_vector< T > & b)
{
        debug_handler_l(DM_BV, "in member - function "
                        "reverse(const base_vector< T > &)", DV_BV + 7);

        register lidia_size_t i, j;

        if (this == &b)
                for (i = 0, j = this->length - 1; i < j; i++, j--)
                        LiDIA::swap(this->value[i], this->value[j]);
        else {
                set_capacity(b.length);
                this->length = b.length;

                for (i = 0; i < this->length; i++)
                        this->value[i] = b.value[this->length - i - 1];
        }
}



//
// swap functions
//

template< class T >
void
base_vector< T >::swap (base_vector< T > & b)
{
        debug_handler_l(DM_BV, "in member - function "
                        "swap(base_vector< T > &)", DV_BV + 4);


        int (*e) (const T & a, const T & b);

        LiDIA::swap(this->length, b.length);
        LiDIA::swap(this->allocated, b.allocated);
        this->bitfield.swap(b.bitfield);

        T *d = this->value;
        this->value = b.value;
        b.value = d;

        char s = this->sort_dir;
        this->sort_dir = b.sort_dir;
        b.sort_dir = s;

        e = this->el_cmp;
        this->el_cmp = b.el_cmp;
        b.el_cmp = e;
}



//
// concat function
//

template< class T >
void
base_vector< T >::concat (const base_vector< T > &b, const base_vector< T > &c)
{
        debug_handler_l(DM_BV, "in member - function "
                        "concat(base_vector< T > &, const base_vector< T > &)", DV_BV + 3);

        register lidia_size_t i, j, l, oldlength;

        if (this != &b && this != &c) {
                l = b.length + c.length;
                set_capacity(l);
                this->length = l;

                for (i = 0; i < b.length; i++)
                        this->value[i] = b.value[i];

                j = b.length;

                for (i = 0; i < c.length; i++ , j++)
                        this->value[j] = c.value[i];
        }
        else {
                if (this == &b) {
                        j = this->length;
                        oldlength = c.length;
                        l = this->length + c.length;
                        set_capacity(l);

                        for (i = 0; i < oldlength; i++, j++)
                                this->value[j] = c.value[i];
                }
                else {
                        base_vector< T > tmp(b.length + c.length, vector_flags::fixed);

                        for (i = 0; i < b.length; i++)
                                tmp.value[i] = b.value[i];

                        j = b.length;

                        for (i = 0; i < c.length; i++, j++)
                                tmp.value[j] = c.value[i];

                        // exchange pointers instead of copying
                        T *p = this->value;
                        this->value = tmp.value;
                        tmp.value = p;

                        tmp.allocated = this->allocated;
                }
        }
        this->allocated = this->length = b.length + c.length;
}



//
// shift functions
//

template< class T >
void
base_vector< T >::shift_left (lidia_size_t pos, lidia_size_t num)
{
        debug_handler_l(DM_BV, "in member - function "
                        "shift_left(lidia_size_t, lidia_size_t)", DV_BV + 3);

        register lidia_size_t i;
        register lidia_size_t old_len = this->length;

        if (pos < num)
                precondition_error_handler (pos, "pos", "pos >= num",
                                     "void base_vector< T >::"
                                     "shift_left(lidia_size_t pos, lidia_size_t num)",
                                     DM_BV, LBV_ERROR[3]);
        else {
                set_size(old_len - num);
                for (i = pos; i < old_len; i++)
                        this->value[i - num] = this->value[i];
        }
}



template< class T >
void
base_vector< T >::shift_right (lidia_size_t pos, lidia_size_t num)
{
        debug_handler_l(DM_BV, "in member - function "
                        "shift_right(lidia_size_t, lidia_size_t)", DV_BV + 3);

        register lidia_size_t i;
        register lidia_size_t old_len = this->length;

        if (pos < 0)
                precondition_error_handler (pos, "pos", "pos >= 0",
                                     "void base_vector< T >::"
                                     "shift_right(lidia_size_t pos, lidia_size_t num)",
                                     DM_BV, LBV_ERROR[3]);
        else {
                set_size(old_len + num);
                for (i = old_len - 1; i >= pos; i--)
                        this->value[i + num] = this->value[i];
        }
}



//
// remove function
//

template< class T >
void
base_vector< T >::remove_from (lidia_size_t pos, lidia_size_t l)
{
        debug_handler_l(DM_BV, "in member - function "
                        "remove_from(lidia_size_t, lidia_size_t)", DV_BV + 9);

        register lidia_size_t i;

        if (l <= 0 || pos< 0 || pos + l > this->length)
                precondition_error_handler (l, "l", "0 < l",
                                     pos, "pos", "0 <= pos and pos + l <= length",
                                     "void sort_vector< T >::"
                                     "remove_from(lidia_size_t pos, lidia_size_t l)",
                                     DM_BV, LBV_ERROR[12]);

        for (i = pos; i < this->length - l; i++)
                this->value[i] = this->value[i + l];
        set_size(this->length - l);
}



//
// insert function
//

template< class T >
void
base_vector< T >::insert_at (const T & x, lidia_size_t pos)
{
        debug_handler_l(DM_BV, "in member - function "
                        "insert_at (const T &, lidia_size_t)", DV_BV + 9);

        register lidia_size_t i, l = this->length;

        if (pos < 0)
                precondition_error_handler(pos, "pos", "pos >= 0",
                                    "void base_vector< T >::"
                                    "insert_at(const T & x, lidia_size_t pos)",
                                    DM_BV , LBV_ERROR[1]);

        set_size((l+1) > (pos + 1) ? l + 1 : pos + 1);

        // if size has been increased, insert at position 'pos'
        for (i = this->length - 1; i > pos; i--)
                this->value[i] = this->value[i - 1];
        this->value[pos] = x;
}



#undef DV_BV
#undef DM_BV
#undef LBV_ERROR



#ifdef LIDIA_NAMESPACE
# ifndef IN_NAMESPACE_LIDIA
}       // end of namespace LiDIA
# endif
#endif



#endif  // LIDIA_BASE_VECTOR_CC_GUARD_
