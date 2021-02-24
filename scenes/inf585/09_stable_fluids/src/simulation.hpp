#pragma once

#include "vcl/vcl.hpp"
#include "boundary.hpp"

enum density_type_structure {density_color, density_texture, view_velocity_curl} ;


void divergence_free(vcl::grid_2D<vcl::vec2>& new_velocity, vcl::grid_2D<vcl::vec2> const& velocity, vcl::grid_2D<float>& divergence, vcl::grid_2D<float>& gradient_field);

template <typename T> void diffuse(vcl::grid_2D<T>& new_field, vcl::grid_2D<T> const& field_reference, float mu, float dt, boundary_condition boundary);
template <typename T> void advect(vcl::grid_2D<T>& new_value, vcl::grid_2D<T> const& value_reference, vcl::grid_2D<vcl::vec2> const& velocity, float dt);


template <typename T>
void diffuse(vcl::grid_2D<T>& new_field, vcl::grid_2D<T> const& field_reference, float mu, float dt, boundary_condition boundary)
{
    using namespace vcl;
    // Compute diffusion on f
    //  Use f as current value, f_prev as previous value
    //  The function is generic in order to handle f as being either a velocity (T=vec2), or a color density (T=vec3)

    // TO do:
    //  Update new_field in solving the diffusion equation for the time step dt
    //  Use Gauss-Seidel iterations over 10 to 15 times
    //   For each iteration, you may reset the boundary conditions with the following code:
    //   	if(boundary==copy)
	//		    set_boundary(new_field);
	//       else
	//		    set_boundary_reflective(new_field);

    size_t N_max = 15;
    int const Nx = int(field_reference.dimension.x);
    int const Ny = int(field_reference.dimension.y);

    float const a = mu * dt;
    float const c = 1.f / (1.f + 4.f * a);

    for (size_t i = 0; i < N_max; ++i) {
        for (int x = 1; x < Nx-1; x++) {
            for (int y = 1; y < Ny - 1; y++) {
                new_field(x, y) = c * (field_reference(x, y) + a * (new_field(x + 1, y) + new_field(x - 1, y) + new_field(x, y - 1) + new_field(x, y + 1)));
            }
        }
        if (boundary == copy)
            set_boundary(new_field);
        else
            set_boundary_reflective(new_field);
    }


}



template <typename T>
void advect(vcl::grid_2D<T>& new_value, vcl::grid_2D<T> const& value_reference, vcl::grid_2D<vcl::vec2> const& velocity, float dt)
{
    using namespace vcl;
    // Compute advection of value along the velocity v, given its previous state value_prev
    int const N = int(new_value.dimension.x);
	
	for(int x=1; x<N-1; ++x) {
        for(int y=1; y<N-1; ++y) {
            vec2 const p  = {float(x),float(y)};

            // retrieve an averaged velocity value around the position p to avoid grid artifacts
			vec2 const v = 0.25f*(interpolation_bilinear(velocity, x+0.5f, y+0.5f)+
				interpolation_bilinear(velocity, x-0.5f, y+0.5f)+
				interpolation_bilinear(velocity, x+0.5f, y-0.5f)+
				interpolation_bilinear(velocity, x-0.5f, y-0.5f));

            // Back tracing
			vec2 const p_back = p-dt*v;

            // Retrieve back-traced value and set it to the new value
			new_value(x,y) = interpolation_bilinear(value_reference, p_back.x, p_back.y);

        }
    }
}




