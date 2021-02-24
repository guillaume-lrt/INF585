#include "simulation.hpp"
//#include "boundary.hpp"

using namespace vcl;



void divergence_free(grid_2D<vec2>& new_velocity, grid_2D<vec2> const& velocity, grid_2D<float>& divergence, grid_2D<float>& gradient_field)
{
    // v = projection of v0 on divergence free vector field
    //
    // v : Final vector field to be filled
    // v0: Initial vector field (non divergence free)
    // divergence: temporary buffer used to compute the divergence of v0
    // gradient_field: temporary buffer used to compute v = v0 - nabla(gradient_field)


    // TO do:
    // 1. Compute divergence of v0
    size_t N_max = 15;
    int const Nx = int(velocity.dimension.x);
    int const Ny = int(velocity.dimension.y);

    // compute div(v0)
    for (int x = 1; x < Nx - 1; x++) {
        for (int y = 1; y < Ny - 1; y++) {
            divergence(x, y) = 0.5f*(velocity(x + 1, y).x- velocity(x - 1, y).x + velocity(x, y + 1).y - velocity(x, y - 1).y);
        }
    }

    // 2. Compute gradient_field such that nabla(gradient_field)^2 = div(v0)
    for (size_t i = 0; i < N_max; ++i) {
        for (int x = 1; x < Nx - 1; x++) {
            for (int y = 1; y < Ny - 1; y++) {
                gradient_field(x, y) = 0.25f* (gradient_field(x + 1, y) + gradient_field(x - 1, y) + gradient_field(x, y - 1) + gradient_field(x, y + 1) - divergence(x,y));
            }
        }
        //if (boundary == copy)
        set_boundary(gradient_field);
        //else
        //    set_boundary_reflective(new_field);
    }

    // 3. Compute v = v0 - nabla(gradient_field)
    for (int x = 1; x < Nx - 1; x++) {
        for (int y = 1; y < Ny - 1; y++) {
            vec2 n_grad = { gradient_field(x + 1, y) - gradient_field(x - 1, y), gradient_field(x, y + 1) - gradient_field(x, y - 1) };
            new_velocity(x, y) = velocity(x,y) - 0.5f * n_grad;
        }
    }
}