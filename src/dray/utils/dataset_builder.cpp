// Copyright 2019 Lawrence Livermore National Security, LLC and other
// Devil Ray Developers. See the top-level COPYRIGHT file for details.
//
// SPDX-License-Identifier: (BSD-3-Clause)

#include <dray/utils/dataset_builder.hpp>


namespace dray
{
  //
  // HexRecord definitions
  //

  /** HexRecord() : Keeps consistent ordering from input. */
  HexRecord::HexRecord(const std::map<std::string, int32> &scalar_vidx,
                       const std::map<std::string, int32> &scalar_eidx)
    : HexRecord(scalar_vidx, scalar_eidx, {}, {})
  {
  }

  HexRecord::HexRecord(const std::map<std::string, int32> &scalar_vidx,
                       const std::map<std::string, int32> &scalar_eidx,
                       const std::map<std::string, int32> &vector_vidx,
                       const std::map<std::string, int32> &vector_eidx )
    : m_coord_data_initd(false),
      m_coord_data(),
      m_scalar_vidx(scalar_vidx),
      m_scalar_vdata_initd(scalar_vidx.size(), false),
      m_scalar_vdata(scalar_vidx.size()),
      m_scalar_vname(scalar_vidx.size()),
      m_scalar_eidx(scalar_eidx),
      m_scalar_edata_initd(scalar_eidx.size(), false),
      m_scalar_edata(scalar_eidx.size()),
      m_scalar_ename(scalar_eidx.size()),

      m_vector_vidx(vector_vidx),
      m_vector_vdata_initd(vector_vidx.size(), false),
      m_vector_vdata(vector_vidx.size()),
      m_vector_vname(vector_vidx.size()),
      m_vector_eidx(vector_eidx),
      m_vector_edata_initd(vector_eidx.size(), false),
      m_vector_edata(vector_eidx.size()),
      m_vector_ename(vector_eidx.size())

  {
    for (const auto &name_idx : scalar_vidx)
      m_scalar_vname[name_idx.second] = name_idx.first;
    for (const auto &name_idx : scalar_eidx)
      m_scalar_ename[name_idx.second] = name_idx.first;
    for (const auto &name_idx : vector_vidx)
      m_vector_vname[name_idx.second] = name_idx.first;
    for (const auto &name_idx : vector_eidx)
      m_vector_ename[name_idx.second] = name_idx.first;
  }

  /** is_initd_self() */
  bool HexRecord::is_initd_self() const
  {
    return (m_coord_data_initd
            && *std::min_element(m_scalar_vdata_initd.begin(), m_scalar_vdata_initd.end())
            && *std::min_element(m_scalar_edata_initd.begin(), m_scalar_edata_initd.end())
            && *std::min_element(m_vector_vdata_initd.begin(), m_vector_vdata_initd.end())
            && *std::min_element(m_vector_edata_initd.begin(), m_vector_edata_initd.end())
            );
  }

  /** is_initd_extern() */
  bool HexRecord::is_initd_extern(const std::map<std::string, int32> &scalar_vidx,
                                  const std::map<std::string, int32> &scalar_eidx,
                                  const std::map<std::string, int32> &vector_vidx,
                                  const std::map<std::string, int32> &vector_eidx ) const
  {
    bool initd = true;
    initd &= m_coord_data_initd;
    for (const auto & name_idx : scalar_vidx)
      initd &= m_scalar_vdata_initd[m_scalar_vidx.at(name_idx.first)];
    for (const auto & name_idx : scalar_eidx)
      initd &= m_scalar_edata_initd[m_scalar_eidx.at(name_idx.first)];
    for (const auto & name_idx : vector_vidx)
      initd &= m_vector_vdata_initd[m_vector_vidx.at(name_idx.first)];
    for (const auto & name_idx : vector_eidx)
      initd &= m_vector_edata_initd[m_vector_eidx.at(name_idx.first)];
    return initd;
  }

  /** print_uninitd_coords */
  void HexRecord::print_uninitd_coords(bool println) const
  {
    const char *RED = "\u001b[31m";
    const char *NRM = "\u001b[0m";

    if (!m_coord_data_initd)
      printf("%sCoordinate data is uninitialized.%c%s", RED, (println ? '\n' : ' '), NRM);
  }

  /** print_uninitd_fields */
  void HexRecord::print_uninitd_fields(bool println) const
  {
    const char *RED = "\u001b[31m";
    const char *NRM = "\u001b[0m";

    const char end = (println ? '\n' : ' ');
    for (int32 idx = 0; idx < m_scalar_vname.size(); ++idx)
      if (!m_scalar_vdata_initd[idx])
        printf("%sField data (vert) '%s' is uninitialized.%c%s", RED, m_scalar_vname[idx].c_str(), end, NRM);
    for (int32 idx = 0; idx < m_scalar_ename.size(); ++idx)
      if (!m_scalar_edata_initd[idx])
        printf("%sField data (elem) '%s' is uninitialized.%c%s", RED, m_scalar_ename[idx].c_str(), end, NRM);
    for (int32 idx = 0; idx < m_vector_vname.size(); ++idx)
      if (!m_vector_vdata_initd[idx])
        printf("%sField data (vert) '%s' is uninitialized.%c%s", RED, m_vector_vname[idx].c_str(), end, NRM);
    for (int32 idx = 0; idx < m_vector_ename.size(); ++idx)
      if (!m_vector_edata_initd[idx])
        printf("%sField data (elem) '%s' is uninitialized.%c%s", RED, m_vector_ename[idx].c_str(), end, NRM);
  }

  /** reset_extern() */
  void HexRecord::reset_extern(const std::map<std::string, int32> &scalar_vidx,
                               const std::map<std::string, int32> &scalar_eidx,
                               const std::map<std::string, int32> &vector_vidx,
                               const std::map<std::string, int32> &vector_eidx )
  {
    m_coord_data_initd = false;
    for (const auto & name_idx : scalar_vidx)
      m_scalar_vdata_initd[m_scalar_vidx.at(name_idx.first)] = false;
    for (const auto & name_idx : scalar_eidx)
      m_scalar_edata_initd[m_scalar_eidx.at(name_idx.first)] = false;
    for (const auto & name_idx : vector_vidx)
      m_vector_vdata_initd[m_vector_vidx.at(name_idx.first)] = false;
    for (const auto & name_idx : vector_eidx)
      m_vector_edata_initd[m_vector_eidx.at(name_idx.first)] = false;
  }

  /** reset_all() */
  void HexRecord::reset_all()
  {
    m_coord_data_initd = false;
    m_scalar_vdata_initd.clear();
    m_scalar_vdata_initd.resize(m_scalar_vidx.size(), false);
    m_scalar_edata_initd.clear();
    m_scalar_edata_initd.resize(m_scalar_eidx.size(), false);
    m_vector_vdata_initd.clear();
    m_vector_vdata_initd.resize(m_vector_vidx.size(), false);
    m_vector_edata_initd.clear();
    m_vector_edata_initd.resize(m_vector_eidx.size(), false);
  }

  /** reuse_all() */
  void HexRecord::reuse_all()
  {
    m_coord_data_initd = true;
    m_scalar_vdata_initd.clear();
    m_scalar_vdata_initd.resize(m_scalar_vidx.size(), true);
    m_scalar_edata_initd.clear();
    m_scalar_edata_initd.resize(m_scalar_eidx.size(), true);
    m_vector_vdata_initd.clear();
    m_vector_vdata_initd.resize(m_vector_vidx.size(), true);
    m_vector_edata_initd.clear();
    m_vector_edata_initd.resize(m_vector_eidx.size(), true);
  }

  /** coords() */
  const HexRecord::CoordT & HexRecord::coords() const
  {
    return m_coord_data;
  }

  /** coords() */
  void HexRecord::coords(const CoordT &coord_data)
  {
    m_coord_data = coord_data;
    m_coord_data_initd = true;
  }

  /** scalar_vdata() */
  const HexRecord::VScalarT & HexRecord::scalar_vdata(const std::string &fname) const
  {
    return m_scalar_vdata[m_scalar_vidx.at(fname)];
  }

  /** scalar_vdata() */
  void HexRecord::scalar_vdata(const std::string &fname, const VScalarT &vdata)
  {
    const int32 idx = m_scalar_vidx.at(fname);
    m_scalar_vdata[idx] = vdata;
    m_scalar_vdata_initd[idx] = true;
  }

  /** scalar_edata() */
  const HexRecord::EScalarT & HexRecord::scalar_edata(const std::string &fname) const
  {
    return m_scalar_edata[m_scalar_eidx.at(fname)];
  }

  /** scalar_edata() */
  void HexRecord::scalar_edata(const std::string &fname, const EScalarT &edata)
  {
    const int32 idx = m_scalar_eidx.at(fname);
    m_scalar_edata[idx] = edata;
    m_scalar_edata_initd[idx] = true;
  }

  /** vector_vdata() */
  const HexRecord::VVectorT & HexRecord::vector_vdata(const std::string &fname) const
  {
    return m_vector_vdata[m_vector_vidx.at(fname)];
  }

  /** vector_vdata() */
  void HexRecord::vector_vdata(const std::string &fname, const VVectorT &vdata)
  {
    const int32 idx = m_vector_vidx.at(fname);
    m_vector_vdata[idx] = vdata;
    m_vector_vdata_initd[idx] = true;
  }

  /** vector_edata() */
  const HexRecord::EVectorT & HexRecord::vector_edata(const std::string &fname) const
  {
    return m_vector_edata[m_vector_eidx.at(fname)];
  }

  /** vector_edata() */
  void HexRecord::vector_edata(const std::string &fname, const EVectorT &edata)
  {
    const int32 idx = m_vector_eidx.at(fname);
    m_vector_edata[idx] = edata;
    m_vector_edata_initd[idx] = true;
  }




  //
  // DataSetBuilder definitions
  //

  /** DataSetBuilder() */
  DataSetBuilder::DataSetBuilder(ShapeMode shape_mode,
                                 const std::vector<std::string> &scalar_vnames,
                                 const std::vector<std::string> &scalar_enames,
                                 const std::vector<std::string> &vector_vnames,
                                 const std::vector<std::string> &vector_enames )
    : m_shape_mode(shape_mode),
      m_num_elems(0),
      m_coord_data(),
      m_scalar_vdata(scalar_vnames.size()),
      m_scalar_edata(scalar_enames.size()),
      m_vector_vdata(vector_vnames.size()),
      m_vector_edata(vector_enames.size())
  {
    int32 idx;

    idx = 0;
    for (const std::string &fname : scalar_vnames)
      m_scalar_vidx[fname] = idx++;

    idx = 0;
    for (const std::string &fname : scalar_enames)
      m_scalar_eidx[fname] = idx++;

    idx = 0;
    for (const std::string &fname : vector_vnames)
      m_vector_vidx[fname] = idx++;

    idx = 0;
    for (const std::string &fname : vector_enames)
      m_vector_eidx[fname] = idx++;
  }

  /** new_empty_hex_record() */
  HexRecord DataSetBuilder::new_empty_hex_record() const
  {
    if (m_shape_mode != Hex)
      throw std::logic_error("Cannot call new_empty_hex_record() on a non-Hex DataSetBuilder.");
    return HexRecord(m_scalar_vidx, m_scalar_eidx, m_vector_vidx, m_vector_eidx);
  }

  /** add_hex_record() : Copies all registered data fields, then flags them as uninitialized. */
  void DataSetBuilder::add_hex_record(HexRecord &record)
  {
    using VScalarT = HexRecord::VScalarT;
    using EScalarT = HexRecord::EScalarT;
    using VVectorT = HexRecord::VVectorT;
    using EVectorT = HexRecord::EVectorT;
    using CoordT   = HexRecord::CoordT;

    if (m_shape_mode != Hex)
      throw std::logic_error("Cannot call add_hex_record() on a non-Hex DataSetBuilder.");

    if (!record.is_initd_extern(m_scalar_vidx, m_scalar_eidx, m_vector_vidx, m_vector_eidx))
    {
      record.print_uninitd_coords();
      record.print_uninitd_fields();
      throw std::logic_error("Attempt to add to DataSetBuilder, but record is missing fields.");
    }

    constexpr int32 verts_per_elem = 8;
    const int32 vtk_2_lex[8] = {0, 1, 3, 2,  4, 5, 7, 6};

    m_num_elems++;

    const CoordT &cdata = record.coords();
    for (int32 j = 0; j < verts_per_elem; ++j)
      m_coord_data.push_back(cdata.m_data[vtk_2_lex[j]]);

    for (const auto &name_idx : m_scalar_vidx)
    {
      const std::string &fname = name_idx.first;
      const int32 fidx = name_idx.second;
      const VScalarT &fdata = record.scalar_vdata(fname);
      for (int32 j = 0; j < verts_per_elem; ++j)
        m_scalar_vdata[fidx].push_back(fdata.m_data[vtk_2_lex[j]]);
    }

    for (const auto &name_idx : m_scalar_eidx)
    {
      const std::string &fname = name_idx.first;
      const int32 fidx = name_idx.second;
      const EScalarT &fdata = record.scalar_edata(fname);
      m_scalar_edata[fidx].push_back(fdata.m_data[0]);
    }

    for (const auto &name_idx : m_vector_vidx)
    {
      const std::string &fname = name_idx.first;
      const int32 fidx = name_idx.second;
      const VVectorT &fdata = record.vector_vdata(fname);
      for (int32 j = 0; j < verts_per_elem; ++j)
        m_vector_vdata[fidx].push_back(fdata.m_data[vtk_2_lex[j]]);
    }

    for (const auto &name_idx : m_vector_eidx)
    {
      const std::string &fname = name_idx.first;
      const int32 fidx = name_idx.second;
      const EVectorT &fdata = record.vector_edata(fname);
      m_vector_edata[fidx].push_back(fdata.m_data[0]);
    }

    record.reset_extern(m_scalar_vidx, m_scalar_eidx, m_vector_vidx, m_vector_eidx);
  }


  /** to_blueprint() */
  void DataSetBuilder::to_blueprint(conduit::Node &bp_dataset,
                                    const std::string &coordset_name) const
  {
    /*
     * https://llnl-conduit.readthedocs.io/en/latest/blueprint_mesh.html#outputting-meshes-for-visualization
     */

    const int32 n_elems = m_num_elems;
    const int32 n_verts = m_coord_data.size();

    //
    // Init node.
    //
    bp_dataset.reset();
    bp_dataset["state/time"] = (float64) 0.0f;
    bp_dataset["state/cycle"] = (uint64) 0;

    conduit::Node &coordset = bp_dataset["coordsets/" + coordset_name];
    conduit::Node &topo = bp_dataset["topologies/mesh"];
    conduit::Node &fields = bp_dataset["fields"];

    //
    // Coordset.
    //
    coordset["type"] = "explicit";
    conduit::Node &coord_vals = coordset["values"];
    coordset["values/x"].set(conduit::DataType::float64(n_verts));
    coordset["values/y"].set(conduit::DataType::float64(n_verts));
    coordset["values/z"].set(conduit::DataType::float64(n_verts));
    float64 *x_vals = coordset["values/x"].value();
    float64 *y_vals = coordset["values/y"].value();
    float64 *z_vals = coordset["values/z"].value();
    for (int32 vidx = 0; vidx < n_verts; ++vidx)
    {
      x_vals[vidx] = (float64) m_coord_data[vidx][0];
      y_vals[vidx] = (float64) m_coord_data[vidx][1];
      z_vals[vidx] = (float64) m_coord_data[vidx][2];
    }

    //
    // Topology.
    //
    topo["type"] = "unstructured";
    topo["coordset"] = coordset_name;
    topo["elements/shape"] = "hex";
    topo["elements/connectivity"].set(conduit::DataType::int32(n_verts));
    int32 * conn = topo["elements/connectivity"].value();
    std::iota(conn, conn + n_verts, 0);

    //
    // Fields.
    //
    for (const auto &name_idx : m_scalar_vidx)  // Scalar vertex fields.
    {
      const std::string &fname = name_idx.first;
      const int32 fidx = name_idx.second;

      conduit::Node &field = fields[fname];
      field["association"] = "vertex";
      field["type"] = "scalar";
      field["topology"] = "mesh";
      field["values"].set(conduit::DataType::float64(n_verts));

      float64 *out_vals = field["values"].value();
      const std::vector<Vec<Float, 1>> &in_field_data = m_scalar_vdata[fidx];
      for (int32 i = 0; i < in_field_data.size(); ++i)
        out_vals[i] = in_field_data[i][0];
    }

    for (const auto &name_idx : m_scalar_eidx)  // Scalar element fields.
    {
      const std::string &fname = name_idx.first;
      const int32 fidx = name_idx.second;

      conduit::Node &field = fields[fname];
      field["association"] = "element";
      field["type"] = "scalar";
      field["topology"] = "mesh";
      field["values"].set(conduit::DataType::float64(n_elems));

      float64 *out_vals = field["values"].value();
      const std::vector<Vec<Float, 1>> &in_field_data = m_scalar_edata[fidx];
      for (int32 i = 0; i < in_field_data.size(); ++i)
        out_vals[i] = in_field_data[i][0];
    }


    const std::string tangent_names[3] = {"u", "v", "w"};

    for (const auto &name_idx : m_vector_vidx)  // Vector vertex fields.
    {
      const std::string &fname = name_idx.first;
      const int32 fidx = name_idx.second;

      constexpr int32 ncomp = 3;
      conduit::Node &field = fields[fname];
      field["association"] = "vertex";
      field["type"] = "vector";
      field["topology"] = "mesh";
      field["values"];

      float64 * out_vals[ncomp];
      for (int32 d = 0; d < ncomp; ++d)
      {
        field["values"][tangent_names[d]].set(conduit::DataType::float64(n_verts));
        out_vals[d] = field["values"][tangent_names[d]].value();
      }

      const std::vector<Vec<Float, 3>> &in_field_data = m_vector_vdata[fidx];
      for (int32 i = 0; i < in_field_data.size(); ++i)
        for (int32 d = 0; d < ncomp; ++d)
          out_vals[d][i] = in_field_data[i][d];
    }

    for (const auto &name_idx : m_vector_eidx)  // Vector element fields.
    {
      const std::string &fname = name_idx.first;
      const int32 fidx = name_idx.second;

      constexpr int32 ncomp = 3;
      conduit::Node &field = fields[fname];
      field["association"] = "element";
      field["type"] = "vector";
      field["topology"] = "mesh";
      field["values"];

      float64 * out_vals[ncomp];
      for (int32 d = 0; d < ncomp; ++d)
      {
        field["values"][tangent_names[d]].set(conduit::DataType::float64(n_verts));
        out_vals[d] = field["values"][tangent_names[d]].value();
      }

      const std::vector<Vec<Float, 3>> &in_field_data = m_vector_edata[fidx];
      for (int32 i = 0; i < in_field_data.size(); ++i)
        for (int32 d = 0; d < ncomp; ++d)
          out_vals[d][i] = in_field_data[i][d];
    }

  }

}//namespace dray
