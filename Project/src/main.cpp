#include "vcl/vcl.hpp"
#include <iostream>

#include "simulation.hpp"


using namespace vcl;

struct gui_parameters {
	bool display_frame = true;
	bool add_sphere = false;
	bool wireframe = false;
	bool solid = true;
	bool show_normals = false;
};

struct user_interaction_parameters {
	vec2 mouse_prev;
	timer_fps fps_record;
	mesh_drawable global_frame;
	gui_parameters gui;
	bool cursor_on_gui;
	
};
user_interaction_parameters user;

struct scene_environment
{
	camera_around_center camera;
	mat4 projection;
	vec3 light;
};
scene_environment scene_env;

mesh_drawable sphere;

Scene scene = Scene();
//Cube cube1;

segments_drawable cube_wireframe;

timer_event_periodic timer(1.f);
std::vector<particle_structure> particles;

//enum { CYLINDER, SPIRAL_IN, SPIRAL_OUT, CUBE_OUT };

void mouse_move_callback(GLFWwindow* window, double xpos, double ypos);
void window_size_callback(GLFWwindow* window, int width, int height);

void initialize_data();
void display_scene();
void display_interface();
void emit_particle();


int main(int, char* argv[])
{

	std::cout << "Run " << argv[0] << std::endl;

	int const width = 1280, height = 1024;
	GLFWwindow* window = create_window(width, height);
	window_size_callback(window, width, height);
	std::cout << opengl_info_display() << std::endl;;

	imgui_init(window);
	glfwSetCursorPosCallback(window, mouse_move_callback);
	glfwSetWindowSizeCallback(window, window_size_callback);
	
	std::cout<<"Initialize data ..."<<std::endl;
	initialize_data();


	std::cout<<"Start animation loop ..."<<std::endl;
	user.fps_record.start();
	timer.start();
	glEnable(GL_DEPTH_TEST);
	while (!glfwWindowShouldClose(window))
	{
		//scene_env.light = scene_env.camera.position();
		scene_env.light = { 10.f,10.f,10.f };
		user.fps_record.update();
		timer.update();
		//cout << scene_env.camera.center_of_rotation << ", " << scene_env.camera.distance_to_center << endl;
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		glClear(GL_DEPTH_BUFFER_BIT);
		imgui_create_frame();
		if(user.fps_record.event) {
			std::string const title = "VCL Display - "+str(user.fps_record.fps)+" fps";
			glfwSetWindowTitle(window, title.c_str());
		}

		ImGui::Begin("GUI",NULL,ImGuiWindowFlags_AlwaysAutoResize);
		user.cursor_on_gui = ImGui::IsAnyWindowFocused();

		if(user.gui.display_frame) draw(user.global_frame, scene_env);

		emit_particle();
		display_interface();
		float const dt = 0.01f * timer.scale;
		//std::cout << dt << std::endl;
		simulate(scene, particles, dt);
		//simulate(particles, dt);
		display_scene();


		ImGui::End();
		imgui_render_frame(window);
		glfwSwapBuffers(window);
		glfwPollEvents();
	}

	imgui_cleanup();
	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}


void emit_particle()
{
	// Emit particle with random velocity
	//  Assume first that all particles have the same radius and mass
	static buffer<vec3> const color_lut = {{1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1}};
	if (timer.event && user.gui.add_sphere) {
		float const theta = rand_interval(0, 2*pi);
		vec3 const v = vec3(.3f*std::cos(theta), .3f*std::sin(theta), 0.f);
		//vec3 const v = vec3(0.f,0.f, 0.0f);

		particle_structure particle;
		particle.p = {0.1,0,6.};
		//particle.p = { 1.2,0.,4.2 };
		particle.r = 0.16f;
		particle.c = color_lut[int(rand_interval()*color_lut.size())];
		particle.v = v;
		particle.m = 7.f; //

		particles.push_back(particle);
	}
}


void initialize_data()
{
	GLuint const shader_mesh = opengl_create_shader_program(opengl_shader_preset("mesh_vertex"), opengl_shader_preset("mesh_fragment"));
	GLuint const shader_uniform_color = opengl_create_shader_program(opengl_shader_preset("single_color_vertex"), opengl_shader_preset("single_color_fragment"));
	GLuint const texture_white = opengl_texture_to_gpu(image_raw{1,1,image_color_type::rgba,{255,255,255,255}});
	mesh_drawable::default_shader = shader_mesh;
	mesh_drawable::default_texture = texture_white;
	curve_drawable::default_shader = shader_uniform_color;
	segments_drawable::default_shader = shader_uniform_color;
	
	user.global_frame = mesh_drawable(mesh_primitive_frame());
	user.gui.display_frame = true;
	scene_env.camera.distance_to_center = 13.147f;
	scene_env.camera.look_at({1.25f,8.2f,10.2f}, {0.78f,-1.46f,1.4f}, {0,0,1});

	sphere = mesh_drawable(mesh_primitive_sphere());

	Cylinder &c1 = scene.cylinders()[0];
	c1.p1() = { -1.6,0.,5.9f };
	c1.p0() = { 4.f,0.,5.2f };
	c1.is_half = true;
	c1.update_mesh();

	Asset& c1_in = scene.assets()[0];
	c1_in.p_in = c1.p0();
	c1_in.type = Type::CYLINDER_IN;
	c1_in.update_mesh();

	Cylinder& c2 = scene.cylinders()[1];
	c2.p0() = c1_in.p_out;
	c2.p1() = c2.p0();
	c2.p1().z = 3.;
	c2.update_mesh();

	Asset& c2_out = scene.assets()[1];
	c2_out.p_in = c2.p1();
	c2_out.type = Type::CYLINDER_OUT;
	c2_out.update_mesh();

	Cylinder& c3 = scene.cylinders()[2];
	c3.p0() = c2_out.p_out;

	c3.p1() = c3.p0();
	c3.p1() += {-2.f, 0.f, -0.4f};
	c3.is_half = true;
	c3.update_mesh();

	Asset& spiral_1 = scene.assets()[2];
	spiral_1.p_in = c3.p1();
	spiral_1.rotation = -pi / 2.f;
	spiral_1.type = Type::SPIRAL;
	spiral_1.update_mesh();

	Asset& c3_out = scene.assets()[3];
	c3_out.p_in = spiral_1.p_out;
	c3_out.type = Type::CYLINDER_OUT;
	c3_out.rotation = 0.f;
	c3_out.update_mesh();

	Asset& cube_out = scene.assets()[4];
	cube_out.rotation = pi / 2.f;
	cube_out.p0.x = 1.f;
	cube_out.p0.z += -1.f;
	cube_out.type = Type::CUBE_OUT;
	cube_out.update_mesh();

	Cylinder& c4 = scene.cylinders()[3];
	c4.p0() = cube_out.p_out;
	c4.p1() = c4.p0();
	c4.p1() += {-2.f, 0.f, -0.4f};
	c4.is_half = true;
	c4.update_mesh();

	Asset& c4_end = scene.assets()[5];
	c4_end.p_in = c4.p1();
	c4_end.type = Type::CYLINDER_END;
	c4_end.rotation = -pi/2.f;
	c4_end.update_mesh();

	Cylinder& c5 = scene.cylinders()[4];
	c5.p0() = c4_end.p_out;
	c5.p1() = c5.p0();
	c5.p1().z = 5.8f;
	c5.update_mesh();

	Cylinder& c6 = scene.cylinders()[5];
	c6.p1() = { -.5f,0.,6.1f };
	c6.p0() = { -2.2f,0.,6.f };
	c6.N_sample_height = 4;
	c6.N_sample_circ = 30;
	c6.is_half = true;
	c6.update_mesh();
}

void display_scene()
{
	size_t const N = particles.size();
	for(size_t k=0; k<N; ++k)
	{
		particle_structure const& particle = particles[k];
		sphere.shading.color = particle.c;
		sphere.transform.translate = particle.p;
		sphere.transform.scale = particle.r;

		draw(sphere, scene_env);
	}
	
	//sphere.transform.translate = scene.assets()[5].p_in;
	//sphere.transform.scale = .1f;
	//draw(sphere, scene_env);
	//sphere.transform.translate = scene.assets()[5].p_out;
	//sphere.transform.scale = .1f;
	//draw(sphere, scene_env);

	// draw color on vertices to debug
	Cylinder c = scene.cylinders()[2];
	//auto pos = c.get_mesh().position;
	//for (int i = 0; i < pos.size(); i++) {
	//	float t = i / float(pos.size());
	//	sphere.transform.translate = pos[i];
	//	sphere.transform.scale = 0.008f;
	//	vec3 v = t * vec3(1.f, 0.f, 0.f) + (1 - t) * vec3(0.f, 1.f, 0.f);
	//	//std::cout << v << std::endl;
	//	sphere.shading.color = v;
	//	draw(sphere, scene_env);
	//}

	mesh cylinder_mesh;// = c.c_mesh();
	mesh_drawable cylinder_mesh_drawable; // = mesh_drawable(cylinder_mesh);

	for (auto& c : scene.cylinders()) {
		cylinder_mesh = c.get_mesh();
		cylinder_mesh_drawable = mesh_drawable(cylinder_mesh);
		if (user.gui.solid)
			draw(cylinder_mesh_drawable, scene_env);

		if (user.gui.wireframe)
			draw_wireframe(cylinder_mesh_drawable, scene_env, { 0,0,0 });
	}

	mesh cube_mesh = scene.cubes()[0].get_mesh();
	mesh_drawable cube_mesh_drawable = mesh_drawable(cube_mesh);
	if (user.gui.solid) {
		draw(cube_mesh_drawable, scene_env);
	}
	if (user.gui.wireframe) {
		draw_wireframe(cube_mesh_drawable, scene_env, { 0,0,0 });
	}

	for (auto& c : scene.assets()) {
		cylinder_mesh = c.get_mesh();
		cylinder_mesh_drawable = mesh_drawable(cylinder_mesh);
		if (user.gui.solid)
			draw(cylinder_mesh_drawable, scene_env);

		if (user.gui.wireframe)
			draw_wireframe(cylinder_mesh_drawable, scene_env, { 0,0,0 });
	}

	auto faces = scene.assets()[2].faces();
	auto faces_normal = scene.assets()[2].faces_normal();
	// show normals
	if (user.gui.show_normals)
		for (int i = 0; i < faces.size(); i++) {
			//auto mm = mesh_primitive_arrow(cyl_pos[i], cyl_pos[i] + cylinder_in.normal[i]/10.f,0.01f,2.f,1.25f);
			if (norm(faces_normal[i]) > 0.f) {
				auto mm = mesh_primitive_arrow(faces[i], faces[i] + faces_normal[i] / 10.f, 0.01f, 2.f, 1.25f);
				//auto mm = mesh_primitive_cylinder(0.05f, cyl_pos[i], cyl_pos[i]+cylinder_in.normal[i]);
				auto mm_d = mesh_drawable(mm);
				draw(mm_d, scene_env);
			}
		}	
}


void display_interface()
{
	ImGui::Checkbox("Frame", &user.gui.display_frame);
	ImGui::Checkbox("Wireframe", &user.gui.wireframe);
	ImGui::Checkbox("Solid", &user.gui.solid);
	ImGui::Checkbox("Show normals", &user.gui.show_normals);
	ImGui::SliderFloat("Time scale", &timer.scale, 0.05f, 2.0f, "%.2f s");
    ImGui::SliderFloat("Interval create sphere", &timer.event_period, 0.05f, 5.0f, "%.2f s");
    ImGui::Checkbox("Add sphere", &user.gui.add_sphere);
}


void window_size_callback(GLFWwindow* , int width, int height)
{
	glViewport(0, 0, width, height);
	float const aspect = width / static_cast<float>(height);
	scene_env.projection = projection_perspective(50.0f*pi/180.0f, aspect, 0.1f, 100.0f);
}


void mouse_move_callback(GLFWwindow* window, double xpos, double ypos)
{
	vec2 const  p1 = glfw_get_mouse_cursor(window, xpos, ypos);
	vec2 const& p0 = user.mouse_prev;
	glfw_state state = glfw_current_state(window);

	auto& camera = scene_env.camera;
	if(!user.cursor_on_gui){
		if(state.mouse_click_left && !state.key_ctrl)
			scene_env.camera.manipulator_rotate_trackball(p0, p1);
		if(state.mouse_click_left && state.key_ctrl)
			camera.manipulator_translate_in_plane(p1-p0);
		if(state.mouse_click_right)
			camera.manipulator_scale_distance_to_center( (p1-p0).y );
	}

	user.mouse_prev = p1;
}

void opengl_uniform(GLuint shader, scene_environment const& current_scene)
{
	opengl_uniform(shader, "projection", current_scene.projection);
	opengl_uniform(shader, "view", scene_env.camera.matrix_view());
	opengl_uniform(shader, "light", scene_env.light, false);
}



