// Copyright 2019 Lawrence Livermore National Security, LLC and other
// Devil Ray Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <dray/GridFunction/low_order_field.hpp>
#include <dray/policies.hpp>

namespace dray
{

LowOrderField::LowOrderField(Array<Float> values, Assoc assoc)
  : m_assoc(assoc),
    m_values(values)
{
  const Float *values_ptr = m_values.get_device_ptr_const();
  const int32 size = m_values.size();

  RAJA::ReduceMin<reduce_policy, Float> xmin (infinity<Float>());
  RAJA::ReduceMax<reduce_policy, Float> xmax (neg_infinity<Float>());

  RAJA::forall<for_policy> (RAJA::RangeSegment (0, size), [=] DRAY_LAMBDA (int32 ii)
  {
    const Float value = values_ptr[ii];
    xmin.min (value);
    xmax.max (value);
  });
  m_range.include (xmin.get ());
  m_range.include (xmax.get ());
}

LowOrderField::~LowOrderField()
{

}

LowOrderField::Assoc
LowOrderField::assoc() const
{
  return m_assoc;
}

std::vector<Range> LowOrderField::range() const
{
  std::vector<Range> ranges;
  ranges.push_back(m_range);
  return ranges;
}

int32 LowOrderField::order() const
{
  if(m_assoc == Assoc::Vertex)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

std::string LowOrderField::type_name() const
{
  std::string name = "low_order_";
  if(m_assoc == Assoc::Vertex)
  {
    name += "vertex";
  }
  else
  {
    name += "element";
  }
  return name;
}

} // namespace dray
