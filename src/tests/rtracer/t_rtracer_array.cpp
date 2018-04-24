#include "gtest/gtest.h"
#include <rtracer/array.hpp>

TEST(rtracer_array, rtracer_array_basic)
{
  rtracer::Array<int> int_array;
  int_array.resize(2);
  int *host = int_array.get_host_ptr();
  host[0] = 0; 
  host[1] = 1; 

  int *dev = int_array.get_device_ptr();

  ASSERT_EQ(dev[0], 0);
  ASSERT_EQ(dev[1], 1);
  
  dev[0] = 1;
  dev[1] = 2;

  int *host2 = int_array.get_host_ptr();
  ASSERT_EQ(host2[0], 1);
  ASSERT_EQ(host2[1], 2);
  
}
