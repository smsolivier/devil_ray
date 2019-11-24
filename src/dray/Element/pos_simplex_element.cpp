#include <dray/Element/element.hpp>

namespace dray
{

// Template instantiations.
template class TriRefSpace<2u>;
template class TriRefSpace<3u>;

// Template instantiations for general-order triangle/tetrahedral elements.
template class Element_impl<2u, 1, ElemType::Tri, Order::General>;
template class Element_impl<2u, 3, ElemType::Tri, Order::General>;
template class Element_impl<3u, 1, ElemType::Tri, Order::General>;
template class Element_impl<3u, 3, ElemType::Tri, Order::General>;
// If fixed-order implementations are needed as well, add instantiations for them here.


} // namespace dray
