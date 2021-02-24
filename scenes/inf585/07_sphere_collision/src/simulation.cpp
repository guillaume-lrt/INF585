#include "simulation.hpp"

using namespace vcl;




void simulate(std::vector<particle_structure>& particles, float dt_true)
{
	vec3 const g = {0,0,-9.81f};
	size_t const N = particles.size();
	buffer<vec3> cube_sides = { { 0.,0.,-1. }, { 0.,0.,1. }, { 0.,-1.,0. }, { 0.,1.,0. }, { -1.,0.,0. }, { 1.,0.,0. } };
	buffer<vec3> normals = { { 0.,0., 1. }, { 0.,0., -1. }, { 0.,1.,0. }, { 0.,-1.,0. }, { 1.,0.,0. }, { -1.,0.,0. } };
	float const epsilon = 1e-5f;
	float alpha, beta;

	size_t const N_substep = 5;
	float const dt = dt_true / N_substep;
	for (size_t k_substep = 0; k_substep < N_substep; ++k_substep)
	{
		size_t const N = particles.size();
		for (size_t k = 0; k < N; ++k)
		{
			particle_structure& particle = particles[k];

			vec3 const f = particle.m * g;

			particle.v = (1 - 0.9f * dt) * particle.v + dt * f;
			particle.p = particle.p + dt * particle.v;
		}

		for (size_t k = 0; k < N; ++k) {
			for (size_t kp = k+1; kp < N; ++kp) {
				particle_structure& particle = particles[k];

				particle_structure& particle_2 = particles[kp];

				vec3& p = particle.p;
				vec3& p2 = particle_2.p;

				float& r = particle.r;
				float& r2 = particle_2.r;

				vec3& v = particle.v;
				vec3& v2 = particle_2.v;

				alpha = 0.95f;
				beta = 0.9f;
				float gamma = 0.83f;

				if (norm(p - p2) < r + r2) {
					//std::cout << "Collision between spheres" << std::endl;
					//float relative_speed = norm(particle.v) / norm(particle_2.v + 0.00001f);
					float relative_speed = norm(v - v2);
					vec3 u = (p - p2) / norm(p - p2);
					if (relative_speed > 0.2f) {
						//vec3 v_new;
						//vec3 v_new_2;
						if (particle.m != particle_2.m) {
							float const j = 2.f * particle.m * particle_2.m / (particle.m + particle_2.m) * dot(v2 - v, u);
							vec3 J = j * u;
							v = alpha * v + beta * J / particle.m;
							v2 = alpha * v2 + beta * J / particle_2.m;
						}
						else {
							float const j = dot(v - v2, u);
							v = v - alpha * j * u;
							v2 = v2 + alpha * j * u;
							//v_new = alpha * particle.v + beta * (dot(particle_2.v - particle.v, u) * u);
							//v_new_2 = alpha * particle_2.v - beta * (dot(particle_2.v - particle.v, u) * u);
						}
						//particle.v = v_new;
						//particle_2.v = v_new_2;
					}
					else {
						v = gamma * v;
						v2 = gamma * v2;
					}
					float d = r + r2 - norm(p - p2);
					p = p + (d / 2.f + epsilon) * u;
					p2 = p2 - (d / 2.f + epsilon) * u;
				}
			}
		}
		// collision sphere plane
		//#pragma omp parallel for
		for (int k = 0; k < N; ++k)
		{
			particle_structure& particle = particles[k];

			//vec3 const f = particle.m * g;

			//particle.v = (1 - 0.9f * dt) * particle.v + dt * f;
			//particle.p = particle.p + dt * particle.v;

			alpha = 0.9f;
			beta = 0.95f;

			for (size_t i = 0; i < cube_sides.size(); i++) {
				vec3 a = cube_sides[i];
				vec3 n = normals[i];
				float const detection = dot(particle.p - a, n);
				if (detection <= 0)
					std::cout << "is outside" << std::endl;
				if (detection <= particle.r - epsilon) {
					//std::cout << "Collision with ground" << std::endl;
					vec3 v_perp = dot(particle.v, n) * n;
					vec3 v_para = particle.v - v_perp;
					particle.v = alpha * v_para - beta * v_perp;
					float d = particle.r - detection;		// distance from the plane
					particle.p = particle.p + d * n;			// position at the exact point of penetration
					break;
				}
			}
		}
	}
}
