// Copyright 2019 Lawrence Livermore National Security, LLC and other
// Devil Ray Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "t_utils.hpp"
#include "test_config.h"
#include "gtest/gtest.h"
#include <dray/io/blueprint_reader.hpp>
#include <dray/filters/path_lengths.hpp>
#include <conduit_blueprint.hpp>
#include <conduit_relay.hpp>

#define EXAMPLE_MESH_SIDE_DIM 10

TEST (dray_slice, dray_slice)
{
  std::string output_path = prepare_output_dir ();
  std::string output_file =
  conduit::utils::join_file_path (output_path, "uniform");
  remove_test_image (output_file);

  conduit::Node data;
  conduit::blueprint::mesh::examples::braid("uniform",
                                            EXAMPLE_MESH_SIDE_DIM,
                                            EXAMPLE_MESH_SIDE_DIM,
                                            EXAMPLE_MESH_SIDE_DIM,
                                            data);
  //data.print();
  dray::DataSet dataset = dray::BlueprintReader::blueprint_to_dray(data);
  conduit::relay::io::save(data,"uniform", "hdf5");

  dray::PathLengths pl;
  pl.resolution(100,100);
  pl.absorption_field("radial");
  pl.emission_field("radial");
  pl.execute(dataset);

}
