


//=====================================================================
/// Helper function to build the mesh; assumed to live in namespace
/// where dimensions of mesh are defined
//=====================================================================
template<class ELEMENT>
RefineableTriangleMesh<ELEMENT>* build_the_mesh(const double& 
                                                uniform_element_area)
{

  double left_edge_x  = -domain_width / 2.0;
  double right_edge_x =  domain_width / 2.0;
  
  // The boundary is bounded by six distinct boundaries, each
  // represented by its own polyline
  Vector<TriangleMeshCurveSection*> boundary_polyline_pt(5);
  TriangleMeshClosedCurve* closed_curve_pt = 0;
 
  // Boundary 0: Inflow boundary
  Vector<Vector <double> > bound_0(2);
  for(unsigned ipoint=0; ipoint<2; ipoint++)
  {
    // Resize the vector
    bound_0[ipoint].resize(2);
  }
  bound_0[0][0] = left_edge_x;
  bound_0[0][1] = 0.0;
  bound_0[1][0] = left_edge_x;
  bound_0[1][1] = domain_height;
 
  unsigned boundary_id = Inflow_boundary_id;
  boundary_polyline_pt[0] = new TriangleMeshPolyLine(bound_0, boundary_id);
 
  // Boundary 1: Top no-slip wall
  Vector<Vector <double> > bound_1(2);
  for(unsigned ipoint=0; ipoint<2; ipoint++)
  {
    // Resize the vector
    bound_1[ipoint].resize(2);
  }
  bound_1[0][0] = left_edge_x;
  bound_1[0][1] = domain_height;
  bound_1[1][0] = 0.0;
  bound_1[1][1] = domain_height;
 
  boundary_id = No_slip_boundary_id;
  boundary_polyline_pt[1] = new TriangleMeshPolyLine(bound_1, boundary_id);
 
  // Boundary 2: Top exit boundary
  Vector<Vector <double> > bound_2(2);
  for(unsigned ipoint=0; ipoint<2; ipoint++)
  {
    // Resize the vector
    bound_2[ipoint].resize(2);
  }
  bound_2[0][0] = 0.0;
  bound_2[0][1] = domain_height;
  bound_2[1][0] = right_edge_x;
  bound_2[1][1] = domain_height;
 
  boundary_id = Top_slip_boundary_id;
  boundary_polyline_pt[2] = new TriangleMeshPolyLine(bound_2, boundary_id);
 

  // Boundary 3: Outflow boundary
  Vector<Vector <double> > bound_3(2);
  for(unsigned ipoint=0; ipoint<2; ipoint++)
  {
    // Resize the vector
    bound_3[ipoint].resize(2);
  }
  bound_3[0][0] = right_edge_x;
  bound_3[0][1] = domain_height;
  bound_3[1][0] = right_edge_x;
  bound_3[1][1] = 0.0;

  boundary_id = Outflow_boundary_id;
  boundary_polyline_pt[3] = new TriangleMeshPolyLine(bound_3, boundary_id);
 

  // Boundary 4: bottom boundary
  Vector<Vector <double> > bound_4(2);
  for(unsigned ipoint=0; ipoint<2; ipoint++)
  {
    // Resize the vector
    bound_4[ipoint].resize(2);
  }
  bound_4[0][0] = right_edge_x;
  bound_4[0][1] = 0.0;
  bound_4[1][0] = left_edge_x;
  bound_4[1][1] = 0.0;

  boundary_id = Bottom_boundary_id;
  boundary_polyline_pt[4] = new TriangleMeshPolyLine(bound_4, boundary_id);
 
  
  // Create the triangle mesh polygon for outer boundary
  //----------------------------------------------------
  TriangleMeshPolygon *outer_polygon = new TriangleMeshPolygon(boundary_polyline_pt);
 
  // Enable redistribution of polylines
  outer_polygon -> enable_redistribution_of_segments_between_polylines();
 
  // Set the pointer
  closed_curve_pt = outer_polygon;
 
  // Now build the mesh
  //===================
 
  // Use the TriangleMeshParameters object for helping on the manage of the
  // TriangleMesh parameters
  TriangleMeshParameters triangle_mesh_parameters(closed_curve_pt);
 
  // Specify the maximum area element
  triangle_mesh_parameters.element_area() = uniform_element_area;
 
  // Specify the internal open boundaries
  /* triangle_mesh_parameters.internal_open_curves_pt() = inner_open_boundaries_pt; */

  // Identify the upper and lower halves of the enriched regions
  /* triangle_mesh_parameters.add_region_coordinates(Enriched_region_upper_id, region1_point); */
  /* triangle_mesh_parameters.add_region_coordinates(Enriched_region_lower_id, region2_point); */

  // set the high-resolution target element area for the enriched regions
  /* triangle_mesh_parameters.set_target_area_for_region(Enriched_region_upper_id, High_res_element_area); */
  /* triangle_mesh_parameters.set_target_area_for_region(Enriched_region_lower_id, High_res_element_area); */
  
  // Create the mesh
  RefineableTriangleMesh<ELEMENT>* Bulk_mesh_pt = 
    new RefineableTriangleMesh<ELEMENT>(triangle_mesh_parameters);
 
  return Bulk_mesh_pt;
}
