#include "skinning.hpp"

namespace vcl
{
	void normalize_weights(buffer<buffer<float>>& weights)
	{
		size_t const N = weights.size();
		for (size_t k = 0; k < N; ++k) {
			float s = 0.0f;
			//std::cout << weights[k].size() << std::endl;
			for(float w : weights[k]) s += w;
			assert_vcl_no_msg(s>1e-5f);
			for(float& w : weights[k]) w /= s;
		}
	}


	// Linear Blend Skinning
	void skinning_LBS_compute(
		buffer<vec3>& position_skinned,  // position to deform
		buffer<vec3>& normal_skinned,    // normal to deform
		buffer<affine_rt> const& skeleton_current,    // rigid transforms for skeleton joints in current pose
		buffer<affine_rt> const& skeleton_rest_pose,  // rigid transforms of skeleton joints in rest pose
		buffer<vec3> const& position_rest_pose,       // vertex positions of the mesh in rest pose
		buffer<vec3> const& normal_rest_pose,         // normal coordinates of the mesh in rest pose
		rig_structure const& rig)                     // information of the skinning weights (joints and weights associated to a vertex)
	{
		size_t const N_vertex = position_rest_pose.size();
		size_t const N_joint = skeleton_current.size();

		// Sanity check on sizes of buffers
		assert_vcl_no_msg(position_skinned.size()==N_vertex);
		assert_vcl_no_msg(normal_skinned.size()==N_vertex);
		assert_vcl_no_msg(normal_rest_pose.size()==N_vertex);
		assert_vcl_no_msg(skeleton_rest_pose.size()==N_joint);
		assert_vcl_no_msg(rig.joint.size()==N_vertex);
		assert_vcl_no_msg(rig.weight.size()==N_vertex);

		//std::cout << N_vertex << std::endl << N_joint << std::endl << rig.joint[10] << std::endl << rig.weight[10] << std::endl;

		// To do
		//   Compute the Linear Blend Skinning ...

		for (int i = 0; i < N_vertex; i++) {
			mat3 res_rotate;
			//res_rotate.fill(0);
			//std::cout << res_rotate << std::endl;
			vec3 res_transf = { 0.,0.,0. };
			mat4 M;
			int t = 0;
			for (auto& j : rig.joint[i]) {
				//if (N_joint > 5) std::cout << rig.joint[i][t] << std::endl << rig.joint[i] << std::endl << rig.weight[i] << std::endl << std::endl;
				mat4 M_temp = skeleton_current[j].matrix() * inverse(skeleton_rest_pose[j]).matrix();
				//res_rotate = res_rotate + rig.weight[i][t] * M.rotate.matrix();
				//res_transf += rig.weight[i][t] * M.translate;
				M += rig.weight[i][t] * M_temp;
				t += 1;
			}
	
			/*affine_rt res_aff = affine_rt(rotation(res_rotate),res_transf);
			position_skinned[i] = res_aff* position_rest_pose[i];
			normal_skinned[i] = res_aff * normal_rest_pose[i];*/
			position_skinned[i] = (M * vec4(position_rest_pose[i], 1.0f)).xyz();
			normal_skinned[i] = (M * vec4(normal_rest_pose[i], 1.0f)).xyz();
		}

	}

}