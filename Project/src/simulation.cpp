#include "simulation.hpp"

static buffer<uint3> connectivity_grid(size_t Nu, size_t Nv)
{
	buffer<uint3> connectivity;
	for (size_t ku = 0; ku < Nu - 1; ++ku) {
		for (size_t kv = 0; kv < Nv - 1; ++kv) {
			unsigned int k00 = static_cast<unsigned int>(kv + Nv * ku);
			unsigned int k10 = static_cast<unsigned int>(kv + 1 + Nv * ku);
			unsigned int k01 = static_cast<unsigned int>(kv + Nv * (ku + 1));
			unsigned int k11 = static_cast<unsigned int>(kv + 1 + Nv * (ku + 1));

			connectivity.push_back(uint3{ k00, k10, k11 });
			connectivity.push_back(uint3{ k00, k11, k01 });
		}
	}
	return connectivity;
}

mesh mesh_primitive_half_cylinder(float radius, vec3 const& p0, vec3 const& p1, int Nu, int Nv, bool is_closed)
{
	vec3 const p01 = p1 - p0;
	float const L = norm(p01);
	assert_vcl(L > 1e-6f, "Cylinder has 0 length");

	vec3 const dir = p01 / L;
	rotation const R = rotation_between_vector({ 0,0,1 }, dir);

	mesh shape;
	for (size_t ku = 0; ku < size_t(Nu); ++ku) {
		for (size_t kv = 0; kv < size_t(Nv); ++kv) {
			float const u = ku / (Nu - 1.0f);
			float const v = kv / (Nv - 1.0f);

			//float const theta = 2 * pi * v;
			float const theta = pi * v + pi / 2.f;

			// cylinder oriented along local z-axis
			vec3 const q = { radius * std::cos(theta), radius * std::sin(theta), L * u };

			// rotate and translate to p1
			vec3 const p = R * q + p0;

			// normal
			vec3 const n = R * vec3{ std::cos(theta), std::sin(theta), 0 };
			// uv
			vec2 const uv = { u,v };

			shape.position.push_back(p);
			shape.normal.push_back(n);
			shape.uv.push_back(uv);
		}
	}

	shape.connectivity = connectivity_grid(Nu, Nv);

	if (is_closed) {
		shape.push_back(mesh_primitive_disc(radius, p0, dir, Nv).flip_connectivity());
		shape.push_back(mesh_primitive_disc(radius, p1, dir, Nv));
	}

	shape.fill_empty_field();

	return shape;
}

mesh mesh_primitive_cubic_grid_without_top(vec3 const& p000, vec3 const& p100, vec3 const& p110, vec3 const& p010, vec3 const& p001, vec3 const& p101, vec3 const& p111, vec3 const& p011, int Nx, int Ny, int Nz)
{
	assert_vcl(Nx >= 2 && Ny >= 2 && Nz >= 2, "Nx, Ny, Nz must be > 2");

	mesh shape;
	shape.push_back(mesh_primitive_grid(p000, p100, p101, p001, Nx, Nz));
	shape.push_back(mesh_primitive_grid(p100, p110, p111, p101, Ny, Nz));
	shape.push_back(mesh_primitive_grid(p110, p010, p011, p111, Nx, Nz));
	shape.push_back(mesh_primitive_grid(p010, p000, p001, p011, Ny, Nz));
	shape.push_back(mesh_primitive_grid(p100, p000, p010, p110, Nx, Ny));

	return shape;
}


template <class T>
void sphere_object(T c, particle_structure& particle, float alpha, float beta) {
	int s = c.positions().size();
	for (size_t i = 0; i < s; i++) {
		vec3 a = c.positions()[i];
		vec3 n = c.normals()[i];
		float const detection = dot(particle.p - a, n);
		if (detection <= 0)
			std::cout << "is outside" << std::endl;
		if (detection <= particle.r - epsilon) {
			vec3 v_perp = dot(particle.v, n) * n;
			vec3 v_para = particle.v - v_perp;
			particle.v = alpha * v_para - beta * v_perp;
			float d = particle.r - detection;		// distance from the plane
			particle.p = particle.p + d * n;			// position at the exact point of penetration
			break;
		}
	}
}

bool sphere_object(Asset ass, particle_structure& particle, float alpha, float beta) {
	int s = ass.faces().size();
	//#pragma omp parallel for
	for (int i = 0; i < s; i++) {
		vec3 a = ass.faces()[i];
		vec3 n = ass.faces_normal()[i];
		float const detection = dot(particle.p - a, n);
		//if (detection <= 0)
			//std::cout << "is outside" << std::endl;
		if (-epsilon <= detection && detection <= particle.r - epsilon) {
			float d = particle.r - detection;		// distance from the plane
			vec3 P = particle.p + d * n;			// position at the exact point of penetration

			// Compute vectors        
			vec3 A = ass.faces_vertices()[i][0];
			vec3 C = ass.faces_vertices()[i][2];
			vec3 B = ass.faces_vertices()[i][1];
			//vec3 A = { 0.f,0.f,0.f };
			//vec3 B = { 0.f,1.f,0.f };
			//vec3 C = { 0.f,0.f,1.f };
			//vec3 P = { -5.f,-.3f,.3f };
			vec3 v0 = C - A;
			vec3 v1 = B - A;
			vec3 v2 = P - A;

			// Compute dot products
			float dot00 = dot(v0, v0);
			float dot01 = dot(v0, v1);
			float dot02 = dot(v0, v2);
			float dot11 = dot(v1, v1);
			float dot12 = dot(v1, v2);

			// Compute barycentric coordinates
			float invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
			float u = (dot11 * dot02 - dot01 * dot12) * invDenom;
			float v = (dot00 * dot12 - dot01 * dot02) * invDenom;

			// Check if point is in triangle
			//cout << "is in triangle with u = " << u << "and v = " << v << endl;
			if ((u >= 0) && (v >= 0) && (u + v <= 1)) {
				//cout << -epsilon << ", " << detection << endl;
				particle.p = P;
				vec3 v_perp = dot(particle.v, n) * n;
				vec3 v_para = particle.v - v_perp;
				particle.v = alpha * v_para - beta * v_perp;
				//break;
				return true;
			}
		}
	}
	return false;
}

void simulate(Scene scene, std::vector<particle_structure>& particles, float dt_true)
{
	vec3 g = {0,0,-9.81f};
	size_t const N = particles.size();
	//buffer<vec3> cube_sides = { { 0.,0.,-1. }, { 0.,0.,1. }, { 0.,-1.,0. }, { 0.,1.,0. }, { -1.,0.,0. }, { 1.,0.,0. } };
	//buffer<vec3> normals = { { 0.,0., 1. }, { 0.,0., -1. }, { 0.,1.,0. }, { 0.,-1.,0. }, { 1.,0.,0. }, { -1.,0.,0. } };
	//Cylinder c1 = scene.cylinders()[0];
	float alpha, beta;

	size_t const N_substep = 10;
	float const dt = dt_true / N_substep;
	for (size_t k_substep = 0; k_substep < N_substep; ++k_substep)
	{
		size_t const N = particles.size();
		for (int k = 0; k < N; ++k)
		{
			particle_structure& particle = particles[k];
			auto p = particle.p;

			if (particle.gravity == 1)
				if (-2.f < p.x && p.x < -1.7f)
					if (-0.3f < p.y && p.y < 0.3f)
						if (-1.7f < p.z && p.z < -1.f)
							particle.gravity = -1;

			vec3 const f = particle.m * particle.gravity * g;

			particle.v = (1 - 0.9f * dt) * particle.v + dt * f;
			particle.p = particle.p + dt * particle.v;

		}

		// collision between spheres
		if (k_substep == 0) {
			for (size_t k = 0; k < N; ++k) {
				for (size_t kp = k + 1; kp < N; ++kp) {
					particle_structure& particle = particles[k];

					particle_structure& particle_2 = particles[kp];

					vec3& p = particle.p;
					vec3& p2 = particle_2.p;

					float& r = particle.r;
					float& r2 = particle_2.r;

					vec3& v = particle.v;
					vec3& v2 = particle_2.v;

					alpha = 0.99f;
					beta = 0.99f;
					float gamma = 0.83f;

					if (norm(p - p2) < r + r2) {
						//std::cout << "Collision between spheres" << std::endl;
						//float relative_speed = norm(particle.v) / norm(particle_2.v + 0.00001f);
						float relative_speed = norm(v - v2);
						vec3 u = (p - p2) / norm(p - p2);
						if (relative_speed > 0.1f) {
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
		}

		// collision sphere plane
		//#pragma omp parallel for
		for (int k = 0; k < N; ++k)
		{
			particle_structure& particle = particles[k];
			alpha = 1.f;
			beta = 0.2f;
			bool is_intersection = false;
			//if (k_substep == 0)
			//std::cout << particle.p << std::endl;

			if (norm(particle.v) > 0.2f || k_substep == 0) {
				//cout << norm(particle.v) << ", " << epsilon << endl;
				for (auto& c : scene.cylinders()) {

					vec3 A = c.p0(); vec3 B = c.p1();
					vec3 AB = B - A;
					vec3 AC = particle.p - A;
					float n_AB = norm(AB);
					float proj = dot(AC, AB) / n_AB; // norm(AC) * cos_angle;

					vec3 K = A + proj * AB / n_AB;

					float KC = norm(particle.p - K);  // = ||KC||

					//std::cout << proj << ", " << proj_2 << ", " << n_AB << std::endl;
					//std::cout << K << ", " << K_2 << std::endl << std::endl;
					if ((proj >= 0) && (proj <= n_AB)) {
						if (KC < c.radius()) {
							sphere_object(c, particle, alpha, beta);
							is_intersection = true;
							break;
						}
					}
				}
				if (!is_intersection) {
					for (auto& c : scene.cubes()) {
						bool is_inside = false;
						auto p = particle.p;
						auto vertex = c.vertex();
						if (p.y >= vertex[0].y - epsilon && p.y <= vertex[1].y + epsilon)			// check the sphere is inside the cube
							if (p.x >= vertex[0].x - epsilon && p.x <= vertex[3].x + epsilon)
								if (p.z >= vertex[0].z - epsilon && p.z <= vertex[4].z + epsilon) {
									sphere_object(c, particle, alpha, beta);
									is_intersection = true;
									break;
								}
					}
				}
				if (!is_intersection) {
					int count = 0;
					//cout << particle.p << endl;
					for (auto& a : scene.assets()) {
						//Asset a = scene.assets()[0];
						count++;
						if (sphere_object(a, particle, alpha, beta))   // found an intersection
							//cout << count << endl;
							break;
					}
				}
			}
		}
	}
}
