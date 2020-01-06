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

#define EXAMPLE_MESH_SIDE_DIM 2

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

  dray::PathLengths pl;
  pl.execute(dataset);

}
