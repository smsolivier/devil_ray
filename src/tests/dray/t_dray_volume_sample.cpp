// Copyright 2019 Lawrence Livermore National Security, LLC and other
// Devil Ray Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include "t_utils.hpp"
#include "test_config.h"
#include "gtest/gtest.h"
#include <dray/io/blueprint_reader.hpp>
#include <dray/rendering/camera.hpp>
#include <dray/rendering/volume_sampler.hpp>
#include <dray/utils/appstats.hpp>

TEST (dray_sampler, volume_sampler)
{
  std::string root_file = std::string (DATA_DIR) + "laghos_tg.cycle_000350.root";
  std::string output_path = prepare_output_dir ();
  std::string output_file =
  conduit::utils::join_file_path (output_path, "tg_sampling");

//  std::string root_file = std::string (DATA_DIR) + "impeller_p2_000000.root";
//  std::string output_path = prepare_output_dir ();
//  std::string output_file =
//  conduit::utils::join_file_path (output_path, "impeller_faces");

  remove_test_image (output_file);

  dray::Collection dataset = dray::BlueprintReader::load (root_file);

  dray::ColorTable color_table ("Spectral");

  // Camera
  const int c_width  = 1000;
  const int c_height = 1000;

  dray::Camera camera;
  camera.set_width (c_width);
  camera.set_height (c_height);
  camera.reset_to_bounds (dataset.bounds());
  dray::Array<dray::Ray> rays;
  camera.create_rays(rays);

  dray::DataSet domain = dataset.domain(0);
  dray::VolumeSampler sampler(domain);
  sampler.sample(rays);

}
