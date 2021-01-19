/**
	Objective: Reproduce the scene with the rotating bubbles going out of the pot and the smoke appearance.

	Bubble part
	----------
	- Complete the function create_new_bubble 
	   - possibly add new parameters to particle_bubble structure
	- Complete the function compute_bubble_position

	Smoke/billboard part
	-----------
	- Complete the function create_new_sprite
	- Complete the function compute_billboard_position and the display procedure of the billboard in the function display_scene().

*/


#include "vcl/vcl.hpp"
#include <iostream>

using namespace vcl;


struct user_interaction_parameters {
	vec2 mouse_prev;
	bool cursor_on_gui;
	bool display_frame = true;
	bool display_transparent_billboard = true;
	timer_fps fps_record;
};
user_interaction_parameters user;

struct scene_environment
{
	camera_around_center camera;
	mat4 projection;
	vec3 light;
};
scene_environment scene;


struct particle_bubble
{
	vec3 p0;
	float t0;
	vec3 color;
	float radius;
	float RotationRadius;   // radius around the center of rotation
	float offset;   // in [0,2pi] so that the rotation not alwats start at 0
	//  Add parameters you need to store for the bubbles
};

struct particle_billboard
{
	vec3 p0;
	float t0;
	vec3 v0;
};



void mouse_move_callback(GLFWwindow* window, double xpos, double ypos);
void window_size_callback(GLFWwindow* window, int width, int height);

void initialize_data();
void display_scene();
particle_bubble create_new_bubble(float t);
vec3 compute_bubble_position(particle_bubble const& bubble, float t_current);
particle_billboard create_new_billboard(float t);
vec3 compute_billboard_position(particle_billboard const& billboard, float t_current);
template <typename T> void remove_old_particles(std::vector<T>& particles, float t_current, float t_max);

float compute_billboard_alpha(particle_billboard const& billboard, float t_current);


// Visual elements of the scene
mesh_drawable global_frame;
mesh_drawable cooking_pot;
mesh_drawable spoon;
mesh_drawable liquid_surface;
mesh_drawable sphere; // used to display the bubbles
mesh_drawable quad;   // used to display the sprites



// Particles and their timer
std::vector<particle_bubble> bubbles;
timer_event_periodic timer_bubble(0.15f);
	
std::vector<particle_billboard> billboards;
timer_event_periodic timer_billboard(0.05f);

float x2;
float y2;



int main(int, char* argv[])
{

	std::cout << "Run " << argv[0] << std::endl;
	scene.projection = projection_perspective(50*pi/180.0f, 1.0f, 0.1f, 100.0f);

	GLFWwindow* window = create_window(1280,1024);
	window_size_callback(window, 1280, 1024);
	std::cout << opengl_info_display() << std::endl;;

	imgui_init(window);
	glfwSetCursorPosCallback(window, mouse_move_callback);
	glfwSetWindowSizeCallback(window, window_size_callback);
	
	std::cout<<"Initialize data ..."<<std::endl;
	initialize_data();


	std::cout<<"Start animation loop ..."<<std::endl;
	user.fps_record.start();
	glEnable(GL_DEPTH_TEST);
	while (!glfwWindowShouldClose(window))
	{
		scene.light = scene.camera.position();
		user.fps_record.update();
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		glClear(GL_DEPTH_BUFFER_BIT);
		imgui_create_frame();
		ImGui::Begin("GUI",NULL,ImGuiWindowFlags_AlwaysAutoResize);
		user.cursor_on_gui = ImGui::IsAnyWindowFocused();
		ImGui::Checkbox("Display frame", &user.display_frame);
		ImGui::Checkbox("Transparent billboard", &user.display_transparent_billboard);
		ImGui::SliderFloat("Bubble emission rate", &timer_bubble.event_period, 0.01f, 1.0f, "%.2f");
		ImGui::SliderFloat("Billboard emission rate", &timer_billboard.event_period, 0.005f, 1.0f, "%.3f");
		ImGui::SliderFloat("Billboard alpha x0", &x2, 0.1f, 3.0f);			// x2 is the half life of the billboard
		ImGui::SliderFloat("Billboard alpha y0", &y2, 0.1f, 1.0f);			// max alpha value

		if(user.fps_record.event) {
			std::string const title = "VCL Display - "+str(user.fps_record.fps)+" fps";
			glfwSetWindowTitle(window, title.c_str());
		}

		if(user.display_frame)
			draw(global_frame, scene);

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



void initialize_data()
{
	// Initialize basic shaders
	GLuint const shader_mesh = opengl_create_shader_program(opengl_shader_preset("mesh_vertex"), opengl_shader_preset("mesh_fragment"));
	GLuint const texture_white = opengl_texture_to_gpu(image_raw{1,1,image_color_type::rgba,{255,255,255,255}});
	mesh_drawable::default_shader = shader_mesh;
	mesh_drawable::default_texture = texture_white;

	// Global frame to visualize the global coordinates
	global_frame = mesh_drawable(mesh_primitive_frame());

	// Load 3D model of pot
	//   Note: The mesh should be accessible from the directory assets/cauldron.obj 
	//   (check the position of your executable if the loading fails)
	cooking_pot = mesh_drawable(mesh_load_file_obj("assets/cauldron.obj"));
	cooking_pot.shading.color = {0.9f, 0.8f, 0.6f};
	cooking_pot.transform.translate = {-0.1f, -0.3f, 0.0f};
	cooking_pot.transform.scale = 0.43f;

	// Load 3D model of spoon
	spoon = mesh_drawable(mesh_load_file_obj("assets/spoon.obj"));
	spoon.shading.color = {0.9f, 0.8f, 0.6f};
    spoon.transform.translate = {-0.1f, -0.3f, 0};
    spoon.transform.scale     = 0.43f;

	// Create the flat surface representing the liquid in the pot
	liquid_surface = mesh_drawable(mesh_primitive_disc(0.73f,{0,0,0},{0,1,0},60));
    liquid_surface.shading.color   = {0.5f,0.6f,0.8f};
    liquid_surface.shading.phong = {0.7f, 0.3f, 0.0f, 128};

	// Sphere used to display the bubbles
	sphere = mesh_drawable(mesh_primitive_sphere(1.0f));

	// Billboard texture and associated quadrange
	GLuint const texture_billboard = opengl_texture_to_gpu(image_load_png("assets/smoke.png"));
	float const L = 0.35f; //size of the quad
	quad = mesh_drawable(mesh_primitive_quadrangle({-L,-L,0},{L,-L,0},{L,L,0},{-L,L,0}));
	//std::cout << quad.shading.alpha << std::endl;
	quad.texture = texture_billboard;

	x2 = 1.8f; y2 = 0.3f;
}

void display_scene() 
{
	// Display the static elements: pot, spoon and surface
	draw(cooking_pot, scene);
	draw(spoon, scene);
	draw(liquid_surface, scene);

	// Create new bubble if needed
	timer_bubble.update();
	if (timer_bubble.event)
		bubbles.push_back( create_new_bubble(timer_bubble.t) );

	// Create new billboards if needed
	timer_billboard.update();
	if(timer_billboard.event)
		billboards.push_back( create_new_billboard(timer_billboard.t) );
		

	// Display bubbles
	for(size_t k = 0; k < bubbles.size(); ++k)
	{
		vec3 const p = compute_bubble_position(bubbles[k], timer_bubble.t);
		sphere.transform.translate = p;
		sphere.shading.color = bubbles[k].color;
		sphere.transform.scale = bubbles[k].radius;
		draw(sphere, scene);
	}


	// Enable transparency using alpha blending (if display_transparent_billboard is true)	
	if(user.display_transparent_billboard){
		glEnable(GL_BLEND);
		glDepthMask(false);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else
		glDisable(GL_BLEND);

	for(size_t k = 0; k < billboards.size(); ++k)
	{
		// Display the sprites here
		//  Note : to make the sprite constantly facing the camera set
		//       quad.transform.rotation = scene.camera.orientation();
		// ...
		quad.transform.rotate = scene.camera.orientation();
		quad.transform.translate = compute_billboard_position(billboards[k], timer_billboard.t);
		//std::cout << "billboard display: " << compute_billboard_position(billboards[k], timer_billboard.t) << std::endl;
		quad.shading.alpha = compute_billboard_alpha(billboards[k], timer_billboard.t);
		draw(quad, scene);
	}
	glDepthMask(true);
		
	remove_old_particles(bubbles, timer_bubble.t, 3.0f);
	remove_old_particles(billboards, timer_billboard.t, 5.0f);
}



particle_bubble create_new_bubble(float t)
{
	particle_bubble bubble;
	bubble.t0 = t;

	float const theta = rand_interval(0.0f, 2*pi);
    float const radius = rand_interval(0.0f, 0.7f);
	float const dist_to_border = 0.8f - radius;
    bubble.p0     = radius*vec3(std::cos(theta), 0.25, std::sin(theta));
	bubble.radius = rand_interval(0.03f,0.06f);
	bubble.color  = {0.5f+rand_interval(0,0.2f),0.6f+rand_interval(0,0.2f),1.0f-rand_interval(0,0.2f)};
	bubble.offset = rand_interval(0.0f, 2 * pi);
	bubble.RotationRadius = rand_interval(dist_to_border / 3.f, dist_to_border);

	return bubble;
}
vec3 compute_bubble_position(particle_bubble const& bubble, float t_current)
{
	float const t = t_current - bubble.t0;

	vec3 p = { bubble.RotationRadius * std::sin(3*t + bubble.offset), t, bubble.RotationRadius * std::cos(3*t + bubble.offset) };
	p = p + bubble.p0;
	return p;
}

particle_billboard create_new_billboard(float t)
{
	particle_billboard billboard;
	billboard.t0 = t;

	float const theta = rand_interval(0, 2 * pi);
	billboard.p0 = { 0,0,0 };
	billboard.v0 = { .5f * std::sin(theta),.5f,.5f * std::cos(theta)};

	return billboard;
}
vec3 compute_billboard_position(particle_billboard const& billboard, float t_current)
{
	vec3 const g = { 0,-0.4f,0 }; // gravity constant

	float t = t_current - billboard.t0;

	vec3 const p = 0.5f * g * t * t + billboard.v0 * t + billboard.p0; // currently only models the first parabola

	return p;
}

float compute_billboard_alpha(particle_billboard const& billboard, float t_current) {
	// compute the alpha channel as a parabola defined with x2 and y2 (coordinates of the apex)
	float t = t_current - billboard.t0;
	//std::cout << x2 << ", "  << y2 << std::endl;
	bool parabola = true;
	if parabola{
		float const b = 2 * y2 / x2;
		float const a = -b / (2 * x2);
		return a * t * t + b * t;
	}
	else {
		float const alpha = t / (2*x2);
		auto temp = y2*(1 - alpha) * std::sqrt(alpha);
		return std::max(temp, 0);
	}
	return 0.0f;
}

// Generic function allowing to remove particles with a lifetime greater than t_max
template <typename T>
void remove_old_particles(std::vector<T>& particles, float t_current, float t_max)
{
	for (auto it = particles.begin(); it != particles.end();)
	{
		if (t_current - it->t0 > t_max)
			it = particles.erase(it);
		if(it!=particles.end())
			++it;
	}
}




void window_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
	float const aspect = width / static_cast<float>(height);
	scene.projection = projection_perspective(50.0f*pi/180.0f, aspect, 0.1f, 100.0f);
}


void mouse_move_callback(GLFWwindow* window, double xpos, double ypos)
{
	vec2 const  p1 = glfw_get_mouse_cursor(window, xpos, ypos);
	vec2 const& p0 = user.mouse_prev;

	glfw_state state = glfw_current_state(window);

	auto& camera = scene.camera;
	if(!user.cursor_on_gui){
		if(state.mouse_click_left && !state.key_ctrl)
			scene.camera.manipulator_rotate_trackball(p0, p1);
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
	opengl_uniform(shader, "view", scene.camera.matrix_view());
	opengl_uniform(shader, "light", scene.light, false);
}


