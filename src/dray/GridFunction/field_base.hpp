// Copyright 2019 Lawrence Livermore National Security, LLC and other
// Devil Ray Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#ifndef DRAY_FIELD_BASE_HPP
#define DRAY_FIELD_BASE_HPP

#include <string>
#include <vector>
#include <dray/range.hpp>

namespace dray
{

class FieldBase
{
protected:
  std::string m_name;
public:
  virtual ~FieldBase() {}

  std::string name() const
  {
    return m_name;
  }
  void name(const std::string &name)
  {
    m_name = name;
  }

  virtual std::vector<Range> range() const = 0;
  virtual int32 order() const = 0;
  virtual std::string type_name() const = 0;
};

} // namespace dray

#endif // DRAY_FIELD_BASE_HPP
