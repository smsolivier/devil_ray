#ifndef DRAY_COMM_DATA_HPP
#define DRAY_COMM_DATA_HPP

#include <vector>
#include <iostream>

namespace dray
{

class MsgCommData
{
  public:
    MsgCommData()
    {
      m_rank= -1;
    }

    MsgCommData(int r, const std::vector<int> &m)
    {
      m_rank = r;
      m_message = m;
    }

    MsgCommData(const MsgCommData &d)
    {
      m_rank = d.m_rank;
      m_message = d.m_message;
    }

    MsgCommData &operator=(const MsgCommData &d)
    {
      m_rank = d.m_rank;
      m_message = d.m_message;
      return *this;
    }

    int m_rank;
    std::vector<int> m_message;
};

} //namespace dray
#endif
