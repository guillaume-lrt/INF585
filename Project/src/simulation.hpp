#pragma once

#include "vcl/vcl.hpp"

using namespace vcl;

struct particle_structure
{
    vcl::vec3 p; // Position
    vcl::vec3 v; // Speed

    vcl::vec3 c; // Color
    float r;     // Radius
    float m;     // mass
};

class Cube {
    public:
        inline Cube(vec3 p0 = { 0.,0.,0. },
            vec3 p1 = { 1.,0.,0. },
            vec3 p3 = { 0.,1.,0. },
            vec3 p4 = { 0.,0.,1. }, int N = 3) {
            vec3 p2 = p1 + p3 - p0;
            vec3 p7 = p3 + p4 - p0;
            vec3 p5 = p1 + p4 - p0;
            vec3 p6 = p2 + p4 - p0;
            m_vertex = { p0,p1,p2,p3,p4,p5,p6,p7 };
            m_mesh = mesh_primitive_cubic_grid(p0, p1, p2, p3, p4, p5, p6, p7, N, N, N);
        }

        inline ~Cube(){}

        inline const buffer<vec3>& vertex() const { return m_vertex; }
        inline buffer<vec3>& vertex() { return m_vertex; }

        inline const mesh& c_mesh() const { return m_mesh; }
        inline mesh& c_mesh() { return m_mesh; }

    private:
        buffer<vec3> m_vertex;
        //buffer<vec3> sides;
        buffer<vec3> m_normals;
        mesh m_mesh;
};

class Cylinder {
    public:
        int N_sample_height;
        int N_sample_circ;
        float m_height;
        bool is_closed;

        inline Cylinder(float radius = 1.f, vec3 p0 = { 0.,0.,0. }, vec3 p1 = { 1.,1.,1. },
                        int Sample_height = 20, int Sample_circ = 20, bool closed = false) : 
                        m_radius(radius), m_p0(p0), m_p1(p1), N_sample_height(Sample_height), N_sample_circ(Sample_circ), is_closed(closed){
            m_height = norm(p1 - p0);
            m_normals = {}; m_positions = {};
            //N_sample_length = int(m_height / (2 * 3.14f * radius) * (N_sample_circ - 1) + 1.5f);
            //N_sample_length = int(m_height * 10 + 2.f);
            //N_sample_height = int(m_height * 10 + 2.f);
            //std::cout << N_sample_length << std::endl;            
        }

        inline ~Cylinder(){}

        inline void update_mesh() {
            //N_sample_length = int(m_height / (2 * 3.14f * m_radius) * (N_sample_circ - 1) + 1.5f);
            //N_sample_height = int(m_height * 10 + 2.f);
            m_mesh = mesh_primitive_cylinder(m_radius, m_p0, m_p1, N_sample_height, N_sample_circ, is_closed);
            m_mesh.compute_normal();
            buffer<vec3> m_positions_temp = m_mesh.position;
            buffer<vec3> m_normals_temp = m_mesh.normal;
            int N = m_positions_temp.size();
            m_normals = {}; m_positions = {};

            for (int i = 0; i < N; i++) {
                bool is_equal = false;          // check if pos[i] is diff from all other points
                //for (int j = i + 1; j < N; j++) {
                //    if (norm(m_positions_temp[i] - m_positions_temp[j]) < 0.001f) {
                //        is_equal = true;
                //        break;
                //    }
                //}
                if (!is_equal) {            // no other similar points
                    m_normals.push_back(-m_normals_temp[i]);
                    m_positions.push_back(m_positions_temp[i]);
                }
            }
        }

        inline const float& radius() const { return m_radius; }
        inline float& radius() { return m_radius; }

        inline const vec3& p0() const { return m_p0; }
        inline vec3& p0() { return m_p0; }

        inline const vec3& p1() const { return m_p1; }
        inline vec3& p1() { return m_p1; }
        inline void height(float new_height) {
            m_height = new_height;
            vec3 v = normalize(m_p1 - m_p0);
            m_p1 = m_p0 + v * m_height;
        }

        inline const mesh& c_mesh() const { return m_mesh; }
        inline mesh& c_mesh() { 
            if (m_normals.size() == 0)
                update_mesh();
            return m_mesh; 
        }

        inline const buffer<vec3>& positions() const { return m_positions; }
        inline buffer<vec3>& positions() { return m_positions; }

        inline const buffer<vec3>& normals() const { return m_normals; }
        inline buffer<vec3>& normals() { return m_normals; } 

    private:
        float m_radius;
        vec3 m_p0;
        vec3 m_p1;
        mesh m_mesh;
        buffer<vec3> m_positions;
        buffer<vec3> m_normals;
};

class Scene {
    public:
        inline Scene(){
            m_cubes = {};
            Cylinder c1 = Cylinder();
            m_cylinders = {c1};
        }
        inline ~Scene(){}

        inline const buffer<Cube*>& cubes() const { return m_cubes; }
        inline buffer<Cube*>& cubes() { return m_cubes; }

        inline const buffer<Cylinder>& cylinders() const { return m_cylinders; }
        inline buffer<Cylinder>& cylinders() { return m_cylinders; }

    private:
        buffer<Cube*> m_cubes;
        buffer<Cylinder> m_cylinders;
};

void simulate(Scene scene, std::vector<particle_structure>& particles, float dt);


