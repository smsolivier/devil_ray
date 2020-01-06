// Copyright 2019 Lawrence Livermore National Security, LLC and other
// Devil Ray Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <dray/error.hpp>
#include <dray/io/blueprint_reader.hpp>
#include <dray/mfem2dray.hpp>
#include <dray/derived_topology.hpp>
#include <dray/uniform_topology.hpp>
#include <dray/GridFunction/field.hpp>
#include <dray/GridFunction/low_order_field.hpp>
#include <dray/utils/data_logger.hpp>

#include <mfem/fem/conduitdatacollection.hpp>
// conduit includes
#include <conduit.hpp>
#include <conduit_blueprint.hpp>
#include <conduit_relay.hpp>

#include <fstream>

using namespace conduit;

namespace dray
{
namespace detail
{

bool is_high_order(const conduit::Node &dom)
{
  if(dom.has_path("fields"))
  {
    const conduit::Node &fields = dom["fields"];
    const int num_fields= fields.number_of_children();
    for(int t = 0; t < num_fields; ++t)
    {
      const conduit::Node &field = fields.child(t);
      if(field.has_path("basis")) return true;
    }
 }
  return false;
}

std::string append_cycle (const std::string &base, const int cycle)
{
  std::ostringstream oss;

  char fmt_buff[64];
  snprintf (fmt_buff, sizeof (fmt_buff), "%06d", cycle);
  oss.str ("");
  oss << base << "_" << std::string (fmt_buff);
  return oss.str ();
}


class BlueprintTreePathGenerator
{
  public:
  BlueprintTreePathGenerator (const std::string &file_pattern,
                              const std::string &tree_pattern,
                              int num_files,
                              int num_trees,
                              const std::string &protocol,
                              const Node &mesh_index)
  : m_file_pattern (file_pattern), m_tree_pattern (tree_pattern),
    m_num_files (num_files), m_num_trees (num_trees), m_protocol (protocol),
    m_mesh_index (mesh_index)
  {
  }

  //-------------------------------------------------------------------//
  ~BlueprintTreePathGenerator ()
  {
  }

  //-------------------------------------------------------------------//
  std::string Expand (const std::string pattern, int idx) const
  {
    //
    // Note: This currently only handles format strings:
    // "%05d" "%06d" "%07d"
    //

    std::size_t pattern_idx = pattern.find ("%05d");

    if (pattern_idx != std::string::npos)
    {
      char buff[16];
      snprintf (buff, 16, "%05d", idx);
      std::string res = pattern;
      res.replace (pattern_idx, 4, std::string (buff));
      return res;
    }

    pattern_idx = pattern.find ("%06d");

    if (pattern_idx != std::string::npos)
    {
      char buff[16];
      snprintf (buff, 16, "%06d", idx);
      std::string res = pattern;
      res.replace (pattern_idx, 4, std::string (buff));
      return res;
    }

    pattern_idx = pattern.find ("%07d");

    if (pattern_idx != std::string::npos)
    {
      char buff[16];
      snprintf (buff, 16, "%07d", idx);
      std::string res = pattern;
      res.replace (pattern_idx, 4, std::string (buff));
      return res;
    }
    return pattern;
  }


  //-------------------------------------------------------------------//
  std::string GenerateFilePath (int tree_id) const
  {
    // for now, we only support 1 tree per file.
    int file_id = tree_id;
    return Expand (m_file_pattern, file_id);
  }

  //-------------------------------------------------------------------//
  std::string GenerateTreePath (int tree_id) const
  {
    // the tree path should always end in a /
    std::string res = Expand (m_tree_pattern, tree_id);
    if ((res.size () > 0) && (res[res.size () - 1] != '/'))
    {
      res += "/";
    }
    return res;
  }

  private:
  std::string m_file_pattern;
  std::string m_tree_pattern;
  int m_num_files;
  int m_num_trees;
  std::string m_protocol;
  Node m_mesh_index;
};

void relay_blueprint_mesh_read (const Node &options, Node &data)
{
  std::string full_root_fname = options["root_file"].as_string ();

  // std::cout<<"ROOOT "<<root_fname<<"\n";
  // read the root file, it can be either json or hdf5

  // assume hdf5, but check for json file
  std::string root_protocol = "hdf5";
  char buff[6] = { 0, 0, 0, 0, 0, 0};

  // heuristic, if json, we expect to see "{" in the first 5 chars of the file.
  std::ifstream ifs;
  ifs.open (full_root_fname.c_str ());
  if (!ifs.is_open ())
  {
    DRAY_ERROR ("failed to open relay root file: " + full_root_fname);
  }
  // std::cout<<"OPEN\n";
  ifs.read ((char *)buff, 5);
  ifs.close ();

  std::string test_str (buff);

  if (test_str.find ("{") != std::string::npos)
  {
    root_protocol = "json";
  }

  std::cout<<"OPEN2 "<<root_protocol<<"\n";
  Node root_node;
  relay::io::load (full_root_fname, root_protocol, root_node);

  // std::cout<<"OPEN2\n";

  if (!root_node.has_child ("file_pattern"))
  {
    DRAY_ERROR ("Root file missing 'file_pattern'");
  }

  if (!root_node.has_child ("blueprint_index"))
  {
    DRAY_ERROR ("Root file missing 'blueprint_index'");
  }

  NodeConstIterator itr = root_node["blueprint_index"].children ();
  Node verify_info;
  // TODO, for now lets verify the first mesh index

  const Node &mesh_index = itr.next ();

  if (!blueprint::mesh::index::verify (mesh_index, verify_info[itr.name ()]))
  {
    std::cout << "Mesh Blueprint index verify failed" << std::endl
              << verify_info.to_json () << "\n";
  }

  std::string data_protocol = "hdf5";

  if (root_node.has_child ("protocol"))
  {
    data_protocol = root_node["protocol/name"].as_string ();
  }

  // read the first mesh (all domains ...)

  int num_domains = root_node["number_of_trees"].to_int ();
  if (num_domains != 1)
  {
    DRAY_ERROR ("only supports single domain");
  }

  BlueprintTreePathGenerator gen (root_node["file_pattern"].as_string (),
                                  root_node["tree_pattern"].as_string (),
                                  root_node["number_of_files"].to_int (),
                                  num_domains, data_protocol, mesh_index);

  std::ostringstream oss;

  char domain_fmt_buff[64];
  int domain_id = 0;
  snprintf (domain_fmt_buff, sizeof (domain_fmt_buff), "%06d", domain_id);
  oss.str ("");
  oss << "domain_" << std::string (domain_fmt_buff);

  std::string current, next;
  utils::rsplit_file_path (full_root_fname, current, next);
  std::string domain_file = utils::join_path (next, gen.GenerateFilePath (domain_id));
  relay::io::load (domain_file, data_protocol, data);
}

Array<Float> fill_array(const conduit::Node &values)
{
  Array<Float> res;

  if(!values.dtype().is_float32() &&
     !values.dtype().is_float64())
  {
    return res;
  }

  const int32 size = values.dtype().number_of_elements();
  res.resize(size);
  Float *res_ptr = res.get_host_ptr();

  if(values.dtype().is_float32())
  {
    const float32 *values_ptr = values.value();

    for(int32 i = 0; i < size; ++i)
    {
      res_ptr[i] = static_cast<Float>(values_ptr[i]);
    }
  }
  else
  {
    const float64 *values_ptr = values.value();

    for(int32 i = 0; i < size; ++i)
    {
      res_ptr[i] = static_cast<Float>(values_ptr[i]);
    }
  }

  return res;
}

void uniform_low_order_fields(const conduit::Node &n_dataset, DataSet &dataset)
{
  // we are assuming that this is uniform
  if(n_dataset.has_child("fields"))
  {
    // add all of the fields:
    NodeConstIterator itr = n_dataset["fields"].children();
    while(itr.has_next())
    {
      const Node &n_field = itr.next();
      std::string field_name = itr.name();

      const int num_children = n_field["values"].number_of_children();

      if(n_field["values"].number_of_children() == 0 )
      {
        Array<Float> values = fill_array(n_field["values"]);
        if(values.size() == 0)
        {
          std::cout<<"skipping non-floating point field '"<<field_name<<"'\n";
        }

        std::string assoc_str = n_field["association"].as_string();
        LowOrderField::Assoc assoc;
        if(assoc_str == "vertex")
        {
          assoc = LowOrderField::Assoc::Vertex;
        }
        else
        {
          assoc = LowOrderField::Assoc::Element;
        }

        std::shared_ptr<LowOrderField> field
          = std::make_shared<LowOrderField>(values, assoc);
        field->name(field_name);
        dataset.add_field(field);
      }

      if(n_field["values"].number_of_children() == 3 )
      {
        std::cout<<"skipping vector field\n";
      }
    } //while
  } // if has fields
}

DataSet low_order(const conduit::Node &n_dataset)
{
  const int num_topos = n_dataset["topologies"].number_of_children();
  if(num_topos != 1)
  {
    DRAY_ERROR("Only a single topology is supported");
  }
  const conduit::Node &topo = n_dataset["topologies"].child(0);

  if(topo["type"].as_string() != "uniform")
  {
    DRAY_ERROR("Only uniform topology implemented");
  }
  const std::string cname = topo["coordset"].as_string();
  const conduit::Node coords = n_dataset["coordsets/"+cname];

  const Node &n_dims = coords["dims"];

  int dims_i = n_dims["i"].to_int();
  int dims_j = n_dims["j"].to_int();
  int dims_k = 1;
  bool is_2d = true;

  // check for 3d
  if(n_dims.has_path("k"))
  {
    dims_k = n_dims["k"].to_int();
    is_2d = false;
  }

  Float origin_x = 0.0f;
  Float origin_y = 0.0f;
  Float origin_z = 0.0f;

  if(coords.has_path("origin"))
  {
    const Node &n_origin = coords["origin"];

    if(n_origin.has_child("x"))
    {
      origin_x = n_origin["x"].to_float32();
    }

    if(n_origin.has_child("y"))
    {
      origin_y = n_origin["y"].to_float32();
    }

    if(n_origin.has_child("z"))
    {
      origin_z = n_origin["z"].to_float32();
    }
  }

  Float spacing_x = 1.0f;
  Float spacing_y = 1.0f;
  Float spacing_z = 1.0f;

  if(coords.has_path("spacing"))
  {
    const Node &n_spacing = coords["spacing"];

    if(n_spacing.has_path("dx"))
    {
        spacing_x = n_spacing["dx"].to_float32();
    }

    if(n_spacing.has_path("dy"))
    {
        spacing_y = n_spacing["dy"].to_float32();
    }

    if(n_spacing.has_path("dz"))
    {
        spacing_z = n_spacing["dz"].to_float32();
    }
  }

  Vec<Float,3> spacing{spacing_x, spacing_y, spacing_z};
  Vec<Float,3> origin{origin_x, origin_y, origin_z};
  Vec<int32,3> dims{dims_i, dims_j, dims_k};

  std::shared_ptr<UniformTopology> utopo
    = std::make_shared<UniformTopology>(spacing, origin, dims);

  DataSet dataset(utopo);
  uniform_low_order_fields(n_dataset, dataset);
  return dataset;
}
//-----------------------------------------------------------------------------

template <typename T>
DataSet bp2dray (const conduit::Node &n_dataset)
{
  bool high_order = is_high_order(n_dataset);
  if(high_order)
  {
    std::cout<<"HO\n";
  }
  else
  {
    return low_order(n_dataset);
  }


  using MeshElemT = MeshElem<3u, Quad, General>;
  using FieldElemT = FieldOn<MeshElemT, 1u>;

  mfem::Mesh *mfem_mesh_ptr = mfem::ConduitDataCollection::BlueprintMeshToMesh (n_dataset);

  mfem::Geometry::Type geom_type = mfem_mesh_ptr->GetElementBaseGeometry(0);
  if(geom_type != mfem::Geometry::CUBE)
  {
    DRAY_ERROR("Only hex imports implemented");
  }

  mfem_mesh_ptr->GetNodes ();
  int space_p;

  GridFunction<3> space_data = import_mesh (*mfem_mesh_ptr, space_p);

  Mesh<MeshElemT> mesh (space_data, space_p);

  std::shared_ptr<HexTopology> topo = std::make_shared<HexTopology>(mesh);
  DataSet dataset(topo);

  NodeConstIterator itr = n_dataset["fields"].children ();

  std::string nodes_gf_name = "";
  std::string topo_name = "main";

  if (n_dataset["topologies"].number_of_children () == 0)
  {
    // this should not happen if verify is called before
    DRAY_ERROR ("Blueprint dataset has no topologies");
  }
  else
  {
    std::vector<std::string> names = n_dataset["topologies"].child_names ();
    topo_name = names[0];
  }

  const Node &n_topo = n_dataset["topologies/" + topo_name];
  if (n_topo.has_child ("grid_function"))
  {
    nodes_gf_name = n_topo["grid_function"].as_string ();
  }

  while (itr.has_next ())
  {
    const Node &n_field = itr.next ();
    std::string field_name = itr.name ();

    // skip mesh nodes gf since they are already processed
    // skip attribute fields, they aren't grid functions
    if (field_name != nodes_gf_name && field_name.find ("_attribute") == std::string::npos)
    {
      mfem::GridFunction *grid_ptr =
      mfem::ConduitDataCollection::BlueprintFieldToGridFunction (mfem_mesh_ptr, n_field);
      const mfem::FiniteElementSpace *fespace = grid_ptr->FESpace ();
      const int32 P = fespace->GetOrder (0);
      if (P == 0)
      {
        DRAY_INFO ("Field has unsupported order " << P);
        continue;
      }
      const int components = grid_ptr->VectorDim ();
      if (components == 1)
      {
        std::cout<<"Field "<<field_name<<"\n";

        int field_p;
        GridFunction<1> field_data = import_grid_function<1> (*grid_ptr, field_p);
        Field<FieldElemT> field (field_data, field_p, field_name);

        std::shared_ptr<Field<FieldElemT>> ffield
          = std::make_shared<Field<FieldElemT>>(field);
        dataset.add_field(ffield);
      }
      else if (components == 3)
      {
        Field<FieldElemT> field_x =
          import_vector_field_component<MeshElemT> (*grid_ptr, 0);
        field_x.name(field_name + "_x");

        Field<FieldElemT> field_y =
          import_vector_field_component<MeshElemT> (*grid_ptr, 1);
        field_y.name(field_name + "_y");

        Field<FieldElemT> field_z =
          import_vector_field_component<MeshElemT> (*grid_ptr, 2);
        field_z.name(field_name + "_z");

        std::shared_ptr<Field<FieldElemT>> ffield_x
          = std::make_shared<Field<FieldElemT>>(field_x);
        dataset.add_field(ffield_x);

        std::shared_ptr<Field<FieldElemT>> ffield_y
          = std::make_shared<Field<FieldElemT>>(field_y);
        dataset.add_field(ffield_y);

        std::shared_ptr<Field<FieldElemT>> ffield_z
          = std::make_shared<Field<FieldElemT>>(field_z);
        dataset.add_field(ffield_z);
      }
      else
      {
        DRAY_INFO ("Import field: number of components = " << components << " not supported");
      }
      delete grid_ptr;
      DRAY_INFO ("Imported field name " << field_name);
    }
  }
  delete mfem_mesh_ptr;
  return dataset;
}

DataSet load_bp (const std::string &root_file)
{
  Node options, data;
  options["root_file"] = root_file;
  detail::relay_blueprint_mesh_read (options, data);
  return bp2dray<Float> (data);
}

} // namespace detail

DataSet BlueprintReader::load (const std::string &root_file)
{
  return detail::load_bp (root_file);
}

DataSet BlueprintReader::load (const std::string &root_file, const int cycle)
{
  std::string full_root = detail::append_cycle (root_file, cycle) + ".root";
  return detail::load_bp (full_root);
}

DataSet
BlueprintReader::blueprint_to_dray (const conduit::Node &n_dataset)
{
  return detail::bp2dray<Float> (n_dataset);
}

} // namespace dray
