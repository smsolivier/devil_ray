// Copyright 2019 Lawrence Livermore National Security, LLC and other
// Devil Ray Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef DRAY_DOF_ACCESS_HPP
#define DRAY_DOF_ACCESS_HPP

#include <dray/exports.hpp>
#include <dray/types.hpp>
#include <dray/vec.hpp>

namespace dray
{

//
// SharedDofPtr - support for double indirection  val = dof_array[ele_offsets[dof_idx]];
//
template <typename DofT> struct SharedDofPtr
{
  const int32 *m_offset_ptr; // Points to element dof map, [dof_idx]-->offset
  const DofT *m_dof_ptr; // Beginning of dof data array, i.e. offset==0.

  // Iterator offset dereference operator.
  DRAY_EXEC const DofT &operator[] (const int32 i) const
  {
    return m_dof_ptr[m_offset_ptr[i]];
  }

  // Iterator offset operator.
  DRAY_EXEC SharedDofPtr operator+ (const int32 &i) const
  {
    return { m_offset_ptr + i, m_dof_ptr };
  }

  // Iterator pre-increment operator.
  DRAY_EXEC SharedDofPtr &operator++ ()
  {
    ++m_offset_ptr;
    return *this;
  }

  // Iterator dereference operator.
  DRAY_EXEC const DofT &operator* () const
  {
    return m_dof_ptr[*m_offset_ptr];
  }
};

} //namespace dray

#endif // DRAY_ELEMENT_HPP