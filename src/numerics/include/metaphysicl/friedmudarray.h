//-----------------------------------------------------------------------bl-
//--------------------------------------------------------------------------
//
// MetaPhysicL - A metaprogramming library for physics calculations
//
// Copyright (C) 2013 The PECOS Development Team
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the Version 2.1 GNU Lesser General
// Public License as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc. 51 Franklin Street, Fifth Floor,
// Boston, MA  02110-1301  USA
//
//-----------------------------------------------------------------------el-
//
// $Id: core.h 37197 2013-02-21 05:49:09Z roystgnr $
//
//--------------------------------------------------------------------------


#ifndef METAPHYSICL_FRIEDMUDARRAY_H
#define METAPHYSICL_FRIEDMUDARRAY_H

#include <algorithm>
#include <functional>
#include <stdexcept>
#include <ostream>

#include "metaphysicl/compare_types.h"
#include "metaphysicl/ct_set.h"
#include "metaphysicl/metaphysicl_asserts.h"
#include "metaphysicl/raw_type.h"
#include "metaphysicl/sparsenumberutils.h"
#include "metaphysicl/testable.h"

namespace MetaPhysicL {

// Forward declarations

// Data type T, index type I
template <typename T, typename I>
class FriedmudArray;

// Helper structs

template<typename I1, typename I2, typename S, typename T, bool reverseorder>
struct DotType<FriedmudArray<S,I1>,
               FriedmudArray<T,I2>, reverseorder> {
  typedef
    FriedmudArray
      <typename DotType<S,T,reverseorder>::supertype,
       typename CompareTypes<I1, I2>::supertype>
      supertype;
};

template<typename I1, typename I2, typename S, typename T, bool reverseorder>
struct OuterProductType<FriedmudArray<S, I1>,
                        FriedmudArray<T, I2>, reverseorder> {
  typedef
    FriedmudArray
      <typename OuterProductType<S,T,reverseorder>::supertype,
       typename CompareTypes<I1, I2>::supertype>
      supertype;
};

template<typename S, typename I>
struct SumType<FriedmudArray<S, I> > {
  typedef FriedmudArray<typename SumType<S>::supertype, I> supertype;
};


template <typename T, typename I>
class FriedmudArray :
  public safe_bool<FriedmudArray<T,I> >
{
public:
  typedef T value_type;

  template <unsigned int i>
  struct entry_type {
    typedef value_type type;
  };

  typedef I index_value_type;

  template <typename T2>
  struct rebind {
    typedef FriedmudArray<T2, I> other;
  };

  std::size_t size() const
    { _data.size(); }

  void resize(std::size_t s)
    {}

  FriedmudArray() {}

  FriedmudArray(const T& val) {
    // This makes no sense unless val is 0!
#ifndef NDEBUG
    if (val)
      throw std::domain_error("Cannot initialize FriedmudArray with non-zero scalar");
#endif
  }

  template <typename T2>
  FriedmudArray(const T2& val) {
    // This makes no sense unless val is 0!
#ifndef NDEBUG
    if (val)
      throw std::domain_error("Cannot initialize FriedmudArray with non-zero scalar");
#endif
  }

  template <typename T2, typename I2>
  FriedmudArray(FriedmudArray<T2, I2> src)
    { std::copy(src._data.begin(), src._data.end(), _data.begin()); }

  T& raw_at(unsigned int i)
    { return _data[i]; }

  const T& raw_at(unsigned int i) const
    { return _data[i]; }

  // FIXME: these encapsulation violations are necessary for std::pow
  // until I can figure out the right friend declaration.
  const std::map<I,T>& nude_data() const
    { return _data; }

  std::map<I,T>& nude_data()
    { return _data; }

  /*
  const std::map<I>& nude_indices() const
    { return _indices; }

  std::vector<I>& nude_indices()
    { return _indices; }
  */
/*
  std::size_t runtime_index_of(index_value_type i) const
    {
      typename std::vector<I>::const_iterator it =
        std::lower_bound(_indices.begin(), _indices.end(), i);
      metaphysicl_assert(it != _indices.end());
      std::size_t offset = it - _indices.begin();
      metaphysicl_assert_equal_to(_indices[offset], i);
      return offset;
    }
*/

  T& operator[](index_value_type i)
    { return _data[i]; }

  const T& operator[](index_value_type i) const
    { return _data[i]; }

  template <unsigned int i>
  typename entry_type<i>::type& get() {
    return _data[i];
  }

  template <unsigned int i>
  const typename entry_type<i>::type& get() const {
    return _data[i];
  }

  value_type& insert(unsigned int i)
  {
    return _data[i];
  }

  template <unsigned int i>
  typename entry_type<i>::type& insert() {
    return this->insert(i);
  }

  template <unsigned int i, typename T2>
  void set(const T2& val) {
    _data[i] = val;
  }

  bool boolean_test() const {
    typename std::map<I,T>::iterator it  = _data.begin();
    typename std::map<I,T>::iterator end = _data.end();

    for (; it != end; ++it)
      if (it->second)
        return true;

    return false;
  }

  FriedmudArray<T,I> operator- () const {
    FriedmudArray<I,T> returnval;

    typename std::map<I,T>::iterator it  = _data.begin();
    typename std::map<I,T>::iterator end = _data.end();

    for (; it != end; ++it)
      returnval[it->first] = -it->second;

    return returnval;
  }

  /*
  // Since this is a dynamically allocated sparsity pattern, we can
  // increase it as needed to support e.g. operator+=
  template <typename I2>
  void sparsity_union (std::vector<I2> new_indices) {
    typename std::vector<I>::iterator index_it = _indices.begin();
    typename std::vector<I2>::const_iterator index2_it = new_indices.begin();

    typedef typename CompareTypes<I,I2>::supertype max_index_type;
    max_index_type unseen_indices = 0;

    I maxI = std::numeric_limits<I>::max();

    while (index2_it != new_indices.end()) {
      I idx1 = (index_it == _indices.end()) ? maxI : *index_it;
      I2 idx2 = *index2_it;

      while (idx1 < idx2) {
        ++index_it;
        idx1 = (index_it == _indices.end()) ? maxI : *index_it;
      }

      while ((idx1 == idx2) &&
             (idx1 != maxI)) {
        ++index_it;
        idx1 = (index_it == _indices.end()) ? maxI : *index_it;
        ++index2_it;
        idx2 = (index2_it == new_indices.end()) ? maxI : *index2_it;
      }

      while (idx2 < idx1) {
        ++unseen_indices;
        ++index2_it;
        if (index2_it == new_indices.end())
          break;
        idx2 = *index2_it;
      }
    }

    // The common case is cheap
    if (!unseen_indices)
      return;

    std::vector<T> merged_data(_data.size() + unseen_indices);
    std::vector<I> merged_indices(_indices.size() + unseen_indices);

    typename std::vector<T>::iterator md_it = merged_data.begin();
    typename std::vector<I>::iterator mi_it = merged_indices.begin();

    typename std::vector<T>::iterator d_it = _data.begin();
    typename std::vector<I>::iterator i_it = _indices.begin();
    typename std::vector<I2>::const_iterator i2_it = new_indices.begin();

    for (; mi_it != merged_indices.end(); ++md_it, ++mi_it) {
      if ((i_it == _indices.end()) ||
          ((i2_it != new_indices.end()) &&
           (*i2_it < *i_it))) {
        *mi_it = *i2_it;
        ++i2_it;
      } else {
        if ((i2_it != new_indices.end()) &&
            (*i2_it == *i_it))
          ++i2_it;
        metaphysicl_assert(d_it < _data.end());
        metaphysicl_assert(md_it < merged_data.end());
        *md_it = *d_it;
        *mi_it = *i_it;
        ++d_it;
        ++i_it;
      }
    }

    metaphysicl_assert(i_it  == _indices.end());
    metaphysicl_assert(i2_it == new_indices.end());
    metaphysicl_assert(d_it  == _data.end());
    metaphysicl_assert(md_it == merged_data.end());

    _indices.swap(merged_indices);
    _data.swap(merged_data);
  }
  */


  /*
  // Since this is a dynamically allocated sparsity pattern, we can
  // decrease it when possible for efficiency
  template <typename I2>
  void sparsity_intersection (std::vector<I2> new_indices) {
    typename std::vector<I>::iterator index_it = _indices.begin();
    typename std::vector<I2>::const_iterator index2_it = new_indices.begin();

    typedef typename CompareTypes<I,I2>::supertype max_index_type;
    max_index_type shared_indices = 0;

    I maxI = std::numeric_limits<I>::max();

    while (index2_it != new_indices.end()) {
      I idx1 = (index_it == _indices.end()) ? maxI : *index_it;
      I2 idx2 = *index2_it;

      while (idx1 < idx2) {
        ++index_it;
        idx1 = (index_it == _indices.end()) ? maxI : *index_it;
      }

      while ((idx1 == idx2) &&
             (idx1 != maxI)) {
        ++index_it;
        idx1 = (index_it == _indices.end()) ? maxI : *index_it;
        ++index2_it;
        idx2 = (index2_it == new_indices.end()) ? maxI : *index2_it;
        ++shared_indices;
      }

      while (idx2 < idx1) {
        ++index2_it;
        if (index2_it == new_indices.end())
          break;
        idx2 = *index2_it;
      }
    }

    std::vector<T> merged_data(shared_indices);
    std::vector<I> merged_indices(shared_indices);

    typename std::vector<T>::iterator md_it = merged_data.begin();
    typename std::vector<I>::iterator mi_it = merged_indices.begin();

    typename std::vector<T>::iterator d_it = _data.begin();
    typename std::vector<I>::iterator i_it = _indices.begin();
    typename std::vector<I2>::const_iterator i2_it = new_indices.begin();

    for (; i_it != _indices.end() && i2_it != new_indices.end();
         ++md_it, ++mi_it, ++d_it, ++i_it, ++i2_it) {
      while (*i2_it < *i_it) {
        ++i2_it;
        if (i2_it == new_indices.end())
          break;
      }
      if (i2_it == new_indices.end())
        break;
      while (*i2_it > *i_it) {
          ++i_it;
        if (i_it == _indices.end())
          break;
      }
      if (i_it == _indices.end())
        break;

      *md_it = *d_it;
      *mi_it = *i_it;
    }

    _indices.swap(merged_indices);
    _data.swap(merged_data);
  }
  */


  // Not defineable since !0 != 0
  // FriedmudArray<T,I> operator! () const;


  template <typename T2, typename I2>
  FriedmudArray<T,I>&
    operator+= (const FriedmudArray<T2,I2>& a)
  {
    typename std::map<I,T>::const_iterator it  = a._data.begin();
    typename std::map<I,T>::const_iterator end = a._data.end();

    for (; it != end; ++it)
      _data[it->first] += it->second;

    return *this;
  }


  template <typename T2, typename I2>
  FriedmudArray<T,I>&
    operator-= (const FriedmudArray<T2,I2>& a) {
    typename std::map<I,T>::const_iterator it  = a._data.begin();
    typename std::map<I,T>::const_iterator end = a._data.end();

    for (; it != end; ++it)
      _data[it->first] -= it->second;

    return *this;
  }


  template <typename T2, typename I2>
  FriedmudArray<T,I>&
    operator*= (const FriedmudArray<T2,I2>& a) {
    typename std::map<I,T>::const_iterator it  = a._data.begin();
    typename std::map<I,T>::const_iterator end = a._data.end();

    for (; it != end; ++it)
      _data[it->first] *= it->second;

    return *this;
  }


  template <typename T2, typename I2>
  FriedmudArray<T,I>&
    operator/= (const FriedmudArray<T2,I2>& a) {
    typename std::map<I,T>::iterator it  = a._data.begin();
    typename std::map<I,T>::iterator end = a._data.end();

    for (; it != end; ++it)
      _data[it->first] /= it->second;

    return *this;
  }


  template <typename T2>
  FriedmudArray<T,I>& operator*= (const T2& a) {
    typename std::map<I,T>::iterator it  = _data.begin();
    typename std::map<I,T>::iterator end = _data.end();

    for (; it != end; ++it)
      it->second *= a;

    return *this;
  }

  template <typename T2>
  FriedmudArray<T,I>& operator/= (const T2& a) {
    typename std::map<I,T>::iterator it  = _data.begin();
    typename std::map<I,T>::iterator end = _data.end();

    for (; it != end; ++it)
      it->second /= a;

    return *this;
  }

  template <typename T2, typename I2>
  FriedmudArray
    <typename DotType<T,T2>::supertype,
     typename CompareTypes<I, I2>::supertype>
  dot (const FriedmudArray<T2,I2>& a) const
  {
    typedef typename DotType<T,T2>::supertype TS;
    typedef typename CompareTypes<I, I2>::supertype IS;

    FriedmudArray<TS, IS> returnval;

    // FIXME
    metaphysicl_not_implemented();

    return returnval;
  }

  template <typename T2, typename I2>
  FriedmudArray<
    typename OuterProductType<T,T2>::supertype,
    typename CompareTypes<I, I2>::supertype>
  outerproduct (const FriedmudArray<T2, I2>& a) const
  {
    typedef typename OuterProductType<T,T2>::supertype TS;
    typedef typename CompareTypes<I, I2>::supertype IS;
    FriedmudArray<TS, IS> returnval;

    // FIXME
    metaphysicl_not_implemented();

    return returnval;
  }

private:

  std::map<I, T> _data;
};


//
// Non-member functions
//

template <unsigned int N,
          unsigned int index1=0, typename Data1=void,
          unsigned int index2=0, typename Data2=void,
          unsigned int index3=0, typename Data3=void,
          unsigned int index4=0, typename Data4=void,
          unsigned int index5=0, typename Data5=void,
          unsigned int index6=0, typename Data6=void,
          unsigned int index7=0, typename Data7=void,
          unsigned int index8=0, typename Data8=void>
struct FriedmudArrayOf
{
  typedef
  typename SymmetricCompareTypes<Data1,
    typename SymmetricCompareTypes<Data2,
      typename SymmetricCompareTypes<Data3,
        typename SymmetricCompareTypes<Data4,
          typename SymmetricCompareTypes<Data5,
            typename SymmetricCompareTypes<Data6,
              typename SymmetricCompareTypes<Data7,Data8>::supertype
            >::supertype
          >::supertype
        >::supertype
      >::supertype
    >::supertype
  >::supertype supertype;

  typedef FriedmudArray<supertype, unsigned int> type;
};



template <std::size_t N, unsigned int index, typename T>
struct FriedmudArrayUnitVector
{
  typedef FriedmudArray<T, unsigned int> type;

  static const type value() {
    type returnval;
    returnval[index] = 1;
    return returnval;
  }
};


template <std::size_t N, typename T>
struct FriedmudArrayFullVector
{
  typedef FriedmudArray<T,unsigned int> type;

  static const type value() {
    type returnval;
    for (unsigned int i=0; i != N; ++i)
      returnval[i] = 1;
    return returnval;
  }
};



template <typename T, typename I, typename I2>
inline
FriedmudArray<FriedmudArray<T, I>, I2>
transpose(const FriedmudArray<FriedmudArray<T, I2>, I>& a)
{
  FriedmudArray<FriedmudArray<T, I>, I2> returnval;

  metaphysicl_not_implemented();

  return returnval;
}


template <typename T, typename I>
FriedmudArray<typename SumType<T>::supertype, I>
sum (const FriedmudArray<T, I> &a)
{
  metaphysicl_not_implemented();
}



#define FriedmudArray_op_ab(opname, atype, btype, functorname) \
template <typename T, typename T2, typename I, typename I2> \
inline \
typename Symmetric##functorname##Type<atype,btype>::supertype \
operator opname (const atype& a, const btype& b) \
{ \
  typedef typename Symmetric##functorname##Type<atype,btype>::supertype type; \
  type returnval = a; \
  returnval opname##= b; \
  return returnval; \
}

#define FriedmudArray_op(opname, functorname) \
FriedmudArray_op_ab(opname, FriedmudArray<T MacroComma I>, FriedmudArray<T2 MacroComma I2>, functorname)

FriedmudArray_op(+, Plus)       // Union)
FriedmudArray_op(-, Minus)      // Union)
FriedmudArray_op(*, Multiplies) // Intersection)
FriedmudArray_op(/, Divides)    // First)

// Let's also allow scalar times vector.
// Scalar plus vector, etc. remain undefined in the sparse context.

template <typename T, typename T2, typename I>
inline
typename MultipliesType<FriedmudArray<T2,I>,T,true>::supertype
operator * (const T& a, const FriedmudArray<T2,I>& b)
{
  typename MultipliesType<FriedmudArray<T2,I>,T,true>::supertype returnval;

  typename std::map<I,T>::iterator it  = b.nude_data().begin();
  typename std::map<I,T>::iterator end = b.nude_data().end();

  for (; it != end; ++it)
    returnval[it->first] = a * it->second;

  return returnval;
}

template <typename T, typename T2, typename I>
inline
typename MultipliesType<FriedmudArray<T,I>,T2>::supertype
operator * (const FriedmudArray<T,I>& a, const T2& b)
{
  typename MultipliesType<FriedmudArray<T,I>,T2>::supertype returnval;

  typename std::map<I,T>::iterator it  = a.nude_data().begin();
  typename std::map<I,T>::iterator end = a.nude_data().end();

  for (; it != end; ++it)
    returnval[it->first] = b * it->second;

  return returnval;
}

template <typename T, typename T2, typename I>
inline
typename DividesType<FriedmudArray<T,I>,T2>::supertype
operator / (const FriedmudArray<T,I>& a, const T2& b)
{
  typename DividesType<FriedmudArray<T,I>,T2>::supertype returnval;

  typename std::map<I,T>::iterator it  = a.nude_data().begin();
  typename std::map<I,T>::iterator end = a.nude_data().end();

  for (; it != end; ++it)
    returnval[it->first] = it->second / b;

  return returnval;
}


#define FriedmudArray_operator_binary(opname, functorname) \
template <typename T, typename T2, typename I, typename I2> \
inline \
FriedmudArray<bool, typename CompareTypes<I,I2>::supertype> \
operator opname (const FriedmudArray<T,I>& a, \
                 const FriedmudArray<T2,I2>& b) \
{ \
    metaphysicl_not_implemented(); \
} \
template <typename T, typename T2, typename I> \
inline \
FriedmudArray<bool, I> \
operator opname (const FriedmudArray<T, I>& a, const T2& b) \
{ \
    metaphysicl_not_implemented(); \
} \
template <typename T, typename T2, typename I> \
inline \
FriedmudArray<bool, I> \
operator opname (const T& a, const FriedmudArray<T2,I>& b) \
{ \
    metaphysicl_not_implemented(); \
}

// NOTE: unary functions for which 0-op-0 is true are undefined compile-time
// errors, because there's no efficient way to have them make sense in
// the sparse context.

FriedmudArray_operator_binary(<, less)
// FriedmudArray_operator_binary(<=)
FriedmudArray_operator_binary(>, greater)
// FriedmudArray_operator_binary(>=)
// FriedmudArray_operator_binary(==)
FriedmudArray_operator_binary(!=, not_equal_to)

// FIXME - make && an intersection rather than a union for efficiency
FriedmudArray_operator_binary(&&, logical_and)
FriedmudArray_operator_binary(||, logical_or)

template <typename T, typename I>
inline
std::ostream&
operator<< (std::ostream& output, const FriedmudArray<T, I>& a)
{
  // Enclose the entire output in braces
  output << '{';

  typename std::map<I,T>::iterator it  = a.nude_data().begin();
  typename std::map<I,T>::iterator end = a.nude_data().end();

  for (; it != end; ++it)
      output << ", (" << it->first << ',' << it->second << ')';

  output << '}';
  return output;
}


// CompareTypes, RawType, ValueType specializations

#define FriedmudArray_comparisons(templatename, settype) \
template<typename T, typename I, bool reverseorder> \
struct templatename<FriedmudArray<T,I>, FriedmudArray<T,I>, reverseorder> { \
  typedef FriedmudArray<T,I> supertype; \
}; \
 \
template<typename T, typename T2, typename I, typename I2, bool reverseorder> \
struct templatename<FriedmudArray<T,I>, FriedmudArray<T2,I2>, reverseorder> { \
  typedef FriedmudArray<typename Symmetric##templatename<T, T2, reverseorder>::supertype, \
                            typename CompareTypes<I,I2>::supertype> supertype; \
}; \
 \
template<typename T, typename T2, typename I, bool reverseorder> \
struct templatename<FriedmudArray<T, I>, T2, reverseorder, \
                    typename boostcopy::enable_if<BuiltinTraits<T2> >::type> { \
  typedef FriedmudArray<typename Symmetric##templatename<T, T2, reverseorder>::supertype, I> supertype; \
}

FriedmudArray_comparisons(CompareTypes, Union);
FriedmudArray_comparisons(PlusType, Union);
FriedmudArray_comparisons(MinusType, Union);
FriedmudArray_comparisons(MultipliesType, Intersection);
FriedmudArray_comparisons(DividesType, First);
FriedmudArray_comparisons(AndType, Intersection);
FriedmudArray_comparisons(OrType, Union);


template <typename T, typename I>
struct RawType<FriedmudArray<T, I> >
{
  typedef FriedmudArray<typename RawType<T>::value_type, I> value_type;

  static value_type value(const FriedmudArray<T, I>& a)
    {
      value_type returnval;

      typename std::map<I,T>::iterator it  = a.nude_data().begin();
      typename std::map<I,T>::iterator end = a.nude_data().end();

      for (; it != end; ++it)
        returnval[it->first] = RawType<T>::value(it->second);
      return returnval;
    }
};

template <typename T, typename I>
struct ValueType<FriedmudArray<T, I> >
{
  typedef typename ValueType<T>::type type;
};

} // namespace MetaPhysicL


namespace std {

using MetaPhysicL::CompareTypes;
using MetaPhysicL::FriedmudArray;
using MetaPhysicL::SymmetricCompareTypes;

#define FriedmudArray_std_unary(funcname) \
template <typename T, typename I> \
inline \
FriedmudArray<T, I> \
funcname (FriedmudArray<T, I> a) \
{ \
  typename std::map<I,T>::iterator it  = a.nude_data.begin(); \
  typename std::map<I,T>::iterator end = a.nude_data.end(); \
  \
  for (; it != end; ++it) \
    it->second = std::funcname(it->second); \
 \
  return a; \
}


#define FriedmudArray_std_binary_union(funcname) \
template <typename T, typename T2, typename I, typename I2> \
inline \
FriedmudArray<typename SymmetricCompareTypes<T,T2>::supertype, \
                         typename CompareTypes<I,I2>::supertype> \
funcname (const FriedmudArray<T, I>& a, \
          const FriedmudArray<T2, I2>& b) \
{ \
    metaphysicl_not_implemented(); \
} \
 \
template <typename T, typename T2, typename I> \
inline \
FriedmudArray<typename SymmetricCompareTypes<T,T2>::supertype, I> \
funcname (const FriedmudArray<T, I>& a, const T2& b) \
{ \
    metaphysicl_not_implemented(); \
} \
 \
template <typename T, typename T2, typename I> \
inline \
FriedmudArray<typename SymmetricCompareTypes<T,T2>::supertype, I> \
funcname (const T& a, const FriedmudArray<T2, I>& b) \
{ \
    metaphysicl_not_implemented(); \
}


// Pow needs its own specialization, both to avoid being confused by
// pow<T1,T2> and because pow(x,0) isn't 0.
template <typename T, typename T2, typename I>
inline
FriedmudArray<typename SymmetricCompareTypes<T,T2>::supertype, I>
pow (const FriedmudArray<T, I>& a, const T2& b)
{
  metaphysicl_not_implemented();
}


// NOTE: unary functions for which f(0) != 0 are undefined compile-time
// errors, because there's no efficient way to have them make sense in
// the sparse context.

// FriedmudArray_std_binary(pow) // separate definition
// FriedmudArray_std_unary(exp)
// FriedmudArray_std_unary(log)
// FriedmudArray_std_unary(log10)
FriedmudArray_std_unary(sin)
// FriedmudArray_std_unary(cos)
FriedmudArray_std_unary(tan)
FriedmudArray_std_unary(asin)
// FriedmudArray_std_unary(acos)
FriedmudArray_std_unary(atan)
FriedmudArray_std_binary_union(atan2)
FriedmudArray_std_unary(sinh)
// FriedmudArray_std_unary(cosh)
FriedmudArray_std_unary(tanh)
FriedmudArray_std_unary(sqrt)
FriedmudArray_std_unary(abs)
FriedmudArray_std_unary(fabs)
FriedmudArray_std_binary_union(max)
FriedmudArray_std_binary_union(min)
FriedmudArray_std_unary(ceil)
FriedmudArray_std_unary(floor)
FriedmudArray_std_binary_union(fmod) // TODO: optimize this


template <typename T, typename I>
class numeric_limits<FriedmudArray<T, I> > :
  public MetaPhysicL::raw_numeric_limits<FriedmudArray<T, I>, T> {};

} // namespace std


#endif // METAPHYSICL_FRIEDMUDARRAY_H
