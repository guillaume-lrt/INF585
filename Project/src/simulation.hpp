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
        bool is_open;  // if it's open at the top

        inline Cube(vec3 p0 = { 0.,0.,0. },
            vec3 p1 = { 1.,0.,0. },
            vec3 p3 = { 0.,1.,0. },
            vec3 p4 = { 0.,0.,1. }, bool open = true,  int N = 2) {
            vec3 p2 = p1 + p3 - p0;
            vec3 p7 = p3 + p4 - p0;
            vec3 p5 = p1 + p4 - p0;
            vec3 p6 = p2 + p4 - p0;
            is_open = open;
            m_vertex = { p0,p1,p2,p3,p4,p5,p6,p7 };
            m_mesh = mesh_primitive_cubic_grid(p0, p1, p2, p3, p4, p5, p6, p7, N, N, N);
            m_mesh.compute_normal();
            //m_positions = m_mesh.position;
            //m_normals = m_mesh.normal;
            buffer<vec3> m_positions_temp = m_mesh.position;
            buffer<vec3> m_normals_temp = m_mesh.normal;
            int K = m_positions_temp.size();
            for (int i = 0; i < K; i++) {
                if (m_normals_temp[i].z != -1) {           
                    m_normals.push_back(m_normals_temp[i]);
                    m_positions.push_back(m_positions_temp[i]);
                }
            }
        }

        inline ~Cube(){}

        inline const buffer<vec3>& vertex() const { return m_vertex; }
        inline buffer<vec3>& vertex() { return m_vertex; }

        inline const buffer<vec3>& positions() const { return m_positions; }
        inline buffer<vec3>& positions() { return m_positions; }

        inline const buffer<vec3>& normals() const { return m_normals; }
        inline buffer<vec3>& normals() { return m_normals; }

        inline const mesh& get_mesh() const { return m_mesh; }
        inline mesh& get_mesh() { return m_mesh; }

    private:
        buffer<vec3> m_vertex;
        buffer<vec3> m_positions;
        buffer<vec3> m_normals;
        mesh m_mesh;
};

class Cylinder {
    public:
        int N_sample_height;
        int N_sample_circ;
        float m_height;
        bool is_closed;
        bool is_half;

        inline Cylinder(float radius = .2f, vec3 p0 = { 0.,0.,0. }, vec3 p1 = { 1.,1.,1. },
                        int Sample_height = 3, int Sample_circ = 20, bool closed = false, bool half = false) : 
                        m_radius(radius), m_p0(p0), m_p1(p1), N_sample_height(Sample_height), N_sample_circ(Sample_circ), is_closed(closed), is_half(half){
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
            if (is_half)
                m_mesh = mesh_primitive_half_cylinder(m_radius, m_p0, m_p1, N_sample_height, N_sample_circ, is_closed);
            else
                m_mesh = mesh_primitive_cylinder(m_radius, m_p0, m_p1, N_sample_height, N_sample_circ, is_closed);
            m_mesh.compute_normal();
            m_normals = -m_mesh.normal;
            m_positions = m_mesh.position;
            //buffer<vec3> m_positions_temp = m_mesh.position;
            //buffer<vec3> m_normals_temp = m_mesh.normal;
            //buffer<vec2> uv_temp = {};
            //buffer<vec3> color_temp = {};
            //int N = m_positions_temp.size();
            //m_normals = {}; m_positions = {};

            //for (int i = 0; i < N; i++) {
            //    bool is_equal = false;          // check if pos[i] is diff from all other points
            //    if (!is_half) {
            //        m_normals.push_back(-m_normals_temp[i]);
            //        m_positions.push_back(m_positions_temp[i]);
            //    }
            //    else if (m_normals_temp[i].z <= 0){
            //        m_normals.push_back(-m_normals_temp[i]);
            //        m_positions.push_back(m_positions_temp[i]);
            //        uv_temp.push_back(m_mesh.uv[i]);
            //        color_temp.push_back(m_mesh.color[i]);
            //    }
            //}
            //m_mesh.position = m_positions;
            //m_mesh.normal = m_normals;
            //m_mesh.uv = uv_temp;
            //m_mesh.color = color_temp;
            //printf("test");
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

        inline const mesh& get_mesh() const { return m_mesh; }
        inline mesh& get_mesh() { 
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

class Asset {
    public:
        vec3 p0;
        float scale;
        float rotation;
        

        inline Asset(std::string path = "null", vec3 position = { 0.f,0.f,0.f }, float rescale = 0.18f, float angle = pi/2.f): 
            m_path(path), p0(position), scale(rescale), rotation(angle){}
        inline ~Asset(){}

        inline void update_mesh() {
            m_mesh = mesh_load_file_obj(m_path);
            auto& pos = m_mesh.position;
            pos *= scale;
            vec3 n = { 0.f,0.f,1.f };   // assume only rotation around z axis
            size_t N = pos.size();
            for (size_t i = 0; i < N;i++)
                pos[i] = cos(rotation) * pos[i] + sin(rotation) * cross(n, pos[i]) + (1 - cos(rotation))*dot(pos[i], n) * n;

            m_mesh.position += p0;

            m_mesh.compute_normal();
            m_normals = -m_mesh.normal;
            m_positions = m_mesh.position;
        }

        inline const mesh& get_mesh() const { return m_mesh; }
        inline mesh& get_mesh() {
            if (m_mesh.position.size() == 0)
                update_mesh();
            return m_mesh;
        }

        inline const buffer<vec3>& positions() const { return m_positions; }
        inline buffer<vec3>& positions() { return m_positions; }

        inline const buffer<vec3>& normals() const { return m_normals; }
        inline buffer<vec3>& normals() { return m_normals; }

    private:
        std::string m_path;
        mesh m_mesh;
        buffer<vec3> m_positions;
        buffer<vec3> m_normals;
};

class Scene {
    public:
        inline Scene(){
            Cube cube1 = Cube({ -2.,-1.,0. }, { -2.,1.,0. }, { 0.,-1.,0. }, {-2.,-1.,0.5});
            m_cubes = {cube1};
            Cylinder c1 = Cylinder(); Cylinder c2 = Cylinder(); Cylinder c3 = Cylinder();
            m_cylinders = {c1,c2,c3};

            Asset c_in_1 = Asset("assets/cylinder_in.obj");
            m_assets = { c_in_1 };
        }
        inline ~Scene(){}

        inline const buffer<Cube>& cubes() const { return m_cubes; }
        inline buffer<Cube>& cubes() { return m_cubes; }

        inline const buffer<Cylinder>& cylinders() const { return m_cylinders; }
        inline buffer<Cylinder>& cylinders() { return m_cylinders; }

        inline const buffer<Asset>& assets() const { return m_assets; }
        inline buffer<Asset>& assets() { return m_assets; }

    private:
        buffer<Cube> m_cubes;
        buffer<Cylinder> m_cylinders;
        buffer<Asset> m_assets;
};

void simulate(Scene scene, std::vector<particle_structure>& particles, float dt);
template <class T>
void sphere_object(T c, particle_structure& particle, float alpha, float beta);
//void sphere_object(Cube c, particle_structure& particle, float alpha, float beta); 
//void sphere_object(Asset c, particle_structure& particle, float alpha, float beta);


