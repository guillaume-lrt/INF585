#include "deformation.hpp"

//#include "vcl/math/rotation/rotation.hpp"

using namespace vcl;

void apply_deformation(mesh& shape, // The position of shape are the one to be deformed
	vec2 const& tr,                 // Input gesture of the user in the 2D-screen coordinates - tr must be converted into a transformation applied to the positions of shape
	buffer<vec3> const& position_before_deformation,  // Initial reference position before the deformation
	buffer<vec3> const& normal_before_deformation,    // Initial reference normals before the deformation
	gui_widget const& widget,                         // Current values of the GUI widget
	picking_parameters const& picking,                // Information on the picking point
	rotation const& camera_orientation)               // Current camera orientation - allows to convert the 2D-screen coordinates into 3D coordinates
{
	float const r = widget.falloff; // radius of influence of the deformation
	size_t const N = shape.position.size();
	for (size_t k = 0; k < N; ++k)
	{
		vec3& p_shape = shape.position[k];                             // position to deform
		vec3 const& p_shape_original = position_before_deformation[k]; // reference position before deformation
		vec3 const& p_clicked = picking.p_clicked;                     // 3D position of picked position
		vec3 const& n_clicked = picking.n_clicked;                     // normal of the surface (before deformation) at the picked position

		float const dist = norm( p_clicked - p_shape_original );       // distance between the picked position and the vertex before deformation

		vec3 translation = camera_orientation * vec3(tr, 0.0f);
		float const wi = std::exp(-(dist * dist) / (r * r));

		// TO DO: Implement the deformation models
		// **************************************************************** //
		// ...
		if (widget.deformer_type == deform_translate) // Case of translation
		{
			// Hint: You can convert the 2D translation in screen space into a 3D translation in the view plane in multiplying 
			//       camera_orientation * (tr.x, tr.y, 0)

			//std::cout << widget.deformer_direction << std::endl;
			if (widget.deformer_direction == direction_surface_normal)
				translation = (translation * n_clicked) * n_clicked;


			// Fake deformation (linear translation in the screen space) 
			//   the following lines should be modified to get the expected smooth deformation
			p_shape = p_shape_original + wi * translation;

		}
		if (widget.deformer_type == deform_twist)
		{
			// Deformation to implement
			translation = camera_orientation.matrix_col_z();
			if (widget.deformer_direction == direction_surface_normal)
				translation = (translation * n_clicked) * n_clicked;
			translation = normalize(translation);
			//std::cout << translation << std::endl;
			
			auto R = rotation::axis_angle_to_matrix(translation, wi*tr.x);
			p_shape = R * (p_shape_original - p_clicked) + p_clicked;
		}
		if (widget.deformer_type == deform_scale)
		{
			// Deformation to implement"
			float s = 1 + wi * tr.x;
			//std::cout << s << std::endl;
			p_shape = s * (p_shape_original - p_clicked) + p_clicked;
		}

	}


}