#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width = width;
    state.image_height = height;
    state.image_color = 0;
    state.image_depth = 0;
    //std::cout << "TODO: allocate and initialize state.image_color and state.image_depth." << std::endl;

    state.image_color = new pixel[width * height];

    for (int i = 0; i < width * height; ++i)
        state.image_color[i] = make_pixel(0, 0, 0);

    state.image_depth = new float[width * height];

    for (int i = 0; i < width * height; ++i)
        state.image_depth[i] = 1;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{

    int numTriangles = state.num_vertices / 3.0;

    	switch(type) {
		case render_type::triangle:
		{
			const data_geometry *d_GeometryArr[3];
        		data_geometry g[3];
        		data_vertex v[3];
        		int num_triangles = state.num_vertices / 3;
		
        		int k = 0;

			for (int i = 0; i < num_triangles; ++i)
			{
			    for (int j = 0; j < 3; ++j)
			    {
				v[j].data = &state.vertex_data[k];
				g[j].data = v[j].data;
				state.vertex_shader(v[j], g[j], state.uniform_data);
				d_GeometryArr[j] = &g[j];
				k += state.floats_per_vertex;
			    }

			    rasterize_triangle(state, d_GeometryArr);
			}
			break;
		}

		case render_type::indexed:
		{
			//
		}
		
		case render_type::fan:
		{
			//
		}
		
		case render_type::strip:
		{
			//
		}
	}
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    data_geometry in2; 
    if (face == 6)
    {
        rasterize_triangle(state, in);
        return;
    }
    if(face == 0) {

    }
    if(face == 1) {

    }
    if(face == 2) {

    }
    if(face == 3) {

    }
    if(face == 4) {

    }
    if(face == 5) {

    }
    // std::cout << "TODO: implement clipping. (The current code passes the triangle through without clipping them.)" << std::endl;
    clip_triangle(state, in, face + 1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    vec3 w, x, y, z;
    vec2 v0, v1, v2;
    
    auto width = state.image_width;
    auto height = state.image_height;

    w = {in[0]->gl_Position[3], in[1]->gl_Position[3], in[2]->gl_Position[3]};
    x = {in[0]->gl_Position[0] / w[0], in[1]->gl_Position[0] / w[1], in[2]->gl_Position[0] / w[2]};
    y = {in[0]->gl_Position[1] / w[0], in[1]->gl_Position[1] / w[1], in[2]->gl_Position[1] / w[2]};
    z = {in[0]->gl_Position[2] / w[0], in[1]->gl_Position[2] / w[1], in[2]->gl_Position[2] / w[2]};

    v0 = {(float)((width / 2) * x[0]) + (float)((width / 2) - 0.5),
               (float)((height / 2) * y[0]) + (float)((height / 2) - 0.5)};
    v1 = {(float)((width / 2) * x[1]) + (float)((width / 2) - 0.5),
               (float)((height / 2) * y[1]) + (float)((height / 2) - 0.5)};
    v2 = {(float)((width / 2) * x[2]) + (float)((width / 2) - 0.5),
               (float)((height / 2) * y[2]) + (float)((height / 2) - 0.5)};

    vec3 area = {(v2[0] - v0[0]) * (v1[1] - v0[1]) - (v2[1] - v0[1]) * (v1[0] - v0[0]), 0, 0};

    for (int i = 0; i < width; ++i)
        for (int j = 0; j < height; ++j)
        {
            vec2 p = {(float)i, (float)j};
            area[1] = (p[0] - v1[0]) * (v2[1] - v1[1]) - (p[1] - v1[1]) * (v2[0] - v1[0]);
            area[2] = (p[0] - v2[0]) * (v0[1] - v2[1]) - (p[1] - v2[1]) * (v0[0] - v2[0]);
            float alpha = area[1] / area[0];
            float beta = area[2] / area[0];
            float gamma = 1 - alpha - beta;
            float a_prime, b_prime, g_prime;

            if (alpha >= 0 && beta >= 0 && gamma >= 0)
            {
                float zDepth = (alpha * z[0]) + (beta * z[1]) + (gamma * z[2]);

                if (state.image_depth[i + j * width] > zDepth)
                {

                    data_fragment frag;
                    frag.data = new float[state.floats_per_vertex];
                    data_output d_output;

                    for (int k = 0; k < state.floats_per_vertex; ++k)
                        if (state.interp_rules[k] == interp_type::flat) {
                            frag.data[k] = in[0]->data[k];
			}
                        else if (state.interp_rules[k] == interp_type::noperspective) {
                            frag.data[k] = alpha * in[0]->data[k] + beta * in[1]->data[k] + gamma * in[2]->data[k];
			}
                        else if (state.interp_rules[k] == interp_type::smooth)
                        {
                            float c = (alpha / w[0]) + (beta / w[1]) + (gamma / w[2]);
                            a_prime = (alpha / w[0]) / c;
                            b_prime = (beta / w[1]) / c;
                            g_prime = (gamma / w[2]) / c;
                            frag.data[k] = (a_prime * in[0]->data[k]) + (b_prime * in[1]->data[k]) + (g_prime * in[2]->data[k]);
                        }

                    state.fragment_shader(frag, d_output, state.uniform_data);

                    float r = d_output.output_color[0] * 255;
                    float g = d_output.output_color[1] * 255;
                    float b = d_output.output_color[2] * 255;

                    state.image_color[i + j * width] = make_pixel(r, g, b);
                    state.image_depth[i + j * width] = zDepth;
               }
          }
    }
}

