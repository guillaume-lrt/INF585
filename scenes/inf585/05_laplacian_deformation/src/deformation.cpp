#include "deformation.hpp"

using namespace vcl;

void build_matrix(linear_system_structure& linear_system, constraint_structure const& constraints, mesh const& shape, buffer<vec3> const& initial_position, buffer<buffer<unsigned int> > const& one_ring)
{
    // TO DO: Build and fill the matrix M and the rhs (M q_x/y/z = rhs_x/y/z)
    //
    //  linear_system.rhs_x.resize(...)
    //  for all vertices
    //      linear_system.rhs_x[k] = ....
    //  
    //  Fill the matrix:
    //  linear_system.M.resize(...), linear_system.M.reserve(...)
    //  for all non-zero elements:
    //      linear_system.M.coeffRef(i,j) = ...
    // 
}


void update_deformation(linear_system_structure& linear_system, constraint_structure const& constraints, vcl::mesh& shape, vcl::mesh_drawable& visual, vcl::buffer<vcl::vec3> const& initial_position, vcl::buffer<vcl::buffer<unsigned int> > const& one_ring)
{
    // TO DO: Update the RHS with new constraints and Solve the system
    //
    // For all vertex and contraints
    //   linear_system.rhs_x[k] = ... (+ rhs_y, rhs_z)
    // 
    // linear_system.solver.solve(linear_system.rhs_x) ... or linear_system.solver.solveWithGuess(linear_system.rhs_x, linear_system.guess_x), etc.
}
