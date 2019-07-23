#ifndef DRAY_INTEGER_UTILS_HPP
#define DRAY_INTEGER_UTILS_HPP

#include <dray/exports.hpp>
#include <dray/types.hpp>

namespace dray
{
template <int32 dim>
class MultinomialCoeff;

using BinomialCoeff = MultinomialCoeff<1>;

template <int32 dim>
class MultinomialCoeff
{
  // Invariant: m_val = MultinomialCoeff(n; i, j, k).
  // Invariant: i+j+k = n.
  protected:
    int32 m_val;
    int32 m_n;
    int32 m_ijk[dim + 1];

  public:
    // Before using MultinomialCoeff, call construct(n).
    DRAY_EXEC void construct(int32 n)
    {
      constexpr int32 full_place = dim;
      m_val = 1;
      m_n = n;
      for (int32 d = 0; d <= dim; d++)
        m_ijk[d] = 0;
      m_ijk[full_place] = n;
    }

    // Getters.
    DRAY_EXEC int32        get_val() const { return m_val; }
    DRAY_EXEC int32        get_n()   const { return m_n; }
    DRAY_EXEC const int32 *get_ijk() const { return m_ijk; }

    // slice_over() - Advance to next coefficient along a direction.
    //                Be careful not to slide off Pascal's simplex.
    DRAY_EXEC int32 slide_over(int32 inc_place)
    {
      constexpr int32 dec_place = dim;
      //       n!              n!         k
      // ---------------  =  -------  *  ---
      // (i+1)! M (k-1)!     i! M k!     i+1
      /// if (m_ijk[dec_place])
      m_val *= m_ijk[dec_place];
      m_ijk[dec_place]--;

      m_ijk[inc_place]++;
      if (m_ijk[inc_place])
        m_val /= m_ijk[inc_place];
      return m_val;
    }

    // swap_places() - The multinomial coefficient is symmetric in i, j, k.
    DRAY_EXEC void swap_places(int32 place1, int32 place2 = dim)
    {
      const int32 s = m_ijk[place2];
      m_ijk[place2] = m_ijk[place1];
      m_ijk[place1] = s;
    }
};

}//namespace dray

#endif//DRAY_INTEGER_UTILS_HPP