#pragma once

#include "vcl/vcl.hpp"
#include "vcl/shape/mesh/primitive/mesh_primitive.hpp"

float const epsilon = 1e-4f;
using namespace vcl;
enum class Type { CYLINDER_IN, CYLINDER_OUT, SPIRAL, CUBE_OUT, CYLINDER_END };
using namespace std;

static buffer<uint3> connectivity_grid(size_t Nu, size_t Nv);
mesh mesh_primitive_half_cylinder(float radius, vec3 const& p0, vec3 const& p1, int Nu, int Nv, bool is_closed);
mesh mesh_primitive_cubic_grid_without_top(vec3 const& p000, vec3 const& p100, vec3 const& p110, vec3 const& p010, vec3 const& p001, vec3 const& p101, vec3 const& p111, vec3 const& p011, int Nx, int Ny, int Nz);

struct particle_structure
{
    vcl::vec3 p; // Position
    vcl::vec3 v; // Speed

    vcl::vec3 c; // Color
    float r;     // Radius
    float m;     // mass

    int gravity = 1; // -+ 1 to have inverse gravity
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
            m_mesh = mesh_primitive_cubic_grid_without_top(p0, p1, p2, p3, p4, p5, p6, p7, N, N, N);
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
        vec3 p_in;          // save the coordinates of in and out points
        vec3 p_out; 
        float scale;
        float rotation;
        bool flip_normals;
        Type type;       // if it's a spiral, box...

        inline Asset(std::string path = "null", vec3 position = { 0.f,0.f,0.f }, float rescale = 0.18f, float angle = pi / 2.f, bool flip = false, Type obj_type = Type::CYLINDER_IN) :
            m_path(path), p0(position), scale(rescale), rotation(angle), flip_normals(flip), type(obj_type) {}
        inline ~Asset(){}

        inline void update_mesh() {
            m_mesh = mesh_load_file_obj(m_path);
            auto& pos = m_mesh.position;

            vec3 p_in_temp;
            vec3 p_out_temp;

            // relative position of in and out points for each obj
            switch (type) {
                case  Type::SPIRAL:
                    scale = 1.5f;
                    p_in_temp = { -.95f,1.f,-.05f };
                    p_out_temp = { 0.f,0.f,-.58f };
                    break;

                case  Type::CYLINDER_OUT:
                    p_in_temp = { 0.f,0.f,1.6f };
                    p_out_temp = { 0.f,2.25f,-.25f };
                    break;

                case  Type::CYLINDER_IN:
                    p_in_temp = { 0.f,2.3f,0.05f };
                    p_out_temp = { 0.f,0.f,-1.60f };
                    break;

                case  Type::CUBE_OUT:
                    p_out_temp = { 0.f,2.35f,-.7f };
                    break;

                case  Type::CYLINDER_END:
                    p_in_temp = { 0.f,2.3f,0.05f };
                    p_out_temp = { 0.f,0.f,1.6f };
                    break;
            }

            

            vec3 n = { 0.f,0.f,1.f };   // assume only rotation around z axis
            size_t N = pos.size();
            for (size_t i = 0; i < N; i++)
                Rotation(pos[i], n, rotation);
            //pos[i] = cos(rotation) * pos[i] + sin(rotation) * cross(n, pos[i]) + (1 - cos(rotation))*dot(pos[i], n) * n;
            Rotation(p_in_temp, n, rotation);
            Rotation(p_out_temp, n, rotation);

            pos *= scale;
            p_in_temp *= scale;
            p_out_temp *= scale;

            vec3 translation = p0;
            if (!is_equal(p_in, vec3{ 0.f,0.f,0.f }))
                translation = p_in - p_in_temp;
            if (!is_equal(p_out, vec3{ 0.f,0.f,0.f }))
                translation = p_out - p_out_temp;

            pos += translation;
            p_in = p_in_temp + translation;
            p_out = p_out_temp + translation;
            p0 = translation;

            m_mesh.compute_normal();
            if (flip_normals) {
                m_mesh.normal *= -1;
            }
            m_normals = m_mesh.normal;
            m_positions = m_mesh.position;

            m_faces = {}; m_faces_normal = {}; m_face_vertices = {};

            N = m_mesh.connectivity.size();
            int a, b, c;
            for (size_t i = 0; i < N; i++) {
                a = m_mesh.connectivity[i][0]; b = m_mesh.connectivity[i][1]; c = m_mesh.connectivity[i][2];
                m_faces.push_back((m_mesh.position[a] + m_mesh.position[b] + m_mesh.position[c]) / 3.f);             // average of all vertices
                m_faces_normal.push_back((m_mesh.normal[a] + m_mesh.normal[b] + m_mesh.normal[c]) / 3.f);
                m_face_vertices.push_back({ m_mesh.position[a], m_mesh.position[b], m_mesh.position[c] });
            }
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

        inline const buffer<vec3>& faces() const { return m_faces; }
        inline buffer<vec3>& faces() { return m_faces; }

        inline const buffer<vec3>& faces_normal () const { return m_faces_normal; }
        inline buffer<vec3>& faces_normal() { return m_faces_normal; }

        inline const buffer<buffer_stack3<vec3>>& faces_vertices() const { return m_face_vertices; }
        inline buffer<buffer_stack3<vec3>>& faces_vertices() { return m_face_vertices; }

    private:
        std::string m_path;
        mesh m_mesh;
        buffer<vec3> m_positions;
        buffer<vec3> m_normals;
        buffer<vec3> m_faces;
        buffer<vec3> m_faces_normal;
        buffer<buffer_stack3<vec3>> m_face_vertices;

        inline void Rotation(vec3& v, vec3 n, float rotation) {
            v = cos(rotation) * v + sin(rotation) * cross(n, v) + (1 - cos(rotation)) * dot(v, n) * n;
        }
};

class Scene {
    public:
        inline Scene(){
            Cube cube1 = Cube({ -2.,-2.,0. }, { -2.,1.,0. }, { 0.,-2.,0. }, {-2.,-2,0.5});
            m_cubes = {};
            Cylinder c1 = Cylinder(); Cylinder c2 = Cylinder(); Cylinder c3 = Cylinder(); Cylinder c4 = Cylinder();
            m_cylinders = {c1,c2,c3, c4};

            Asset c1_in= Asset("assets/cylinder_in.obj");
            Asset c2_out = Asset("assets/cylinder_out_ter.obj");
            Asset spiral_1 = Asset("assets/spiral_1.obj");
            Asset cube_out = Asset("assets/box_out_bis.obj");
            Asset c3_out = Asset("assets/cylinder_out_ter.obj");
            Asset c4_in_end = Asset("assets/cylinder_in_end.obj");
            m_assets = { c1_in, c2_out, spiral_1, c3_out, cube_out, c4_in_end };
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
bool sphere_object(Asset a, particle_structure& particle, float alpha, float beta); 
//void sphere_object(Asset c, particle_structure& particle, float alpha, float beta);


