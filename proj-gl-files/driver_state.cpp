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
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;

    int size = state.image_width * state.image_height;
    state.image_color = new pixel[size];
    for(int i = 0; i < size; ++i) {
	state.image_color[i] = make_pixel(0, 0, 0);
    }
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
    //std::cout<<"TODO: implement rendering."<<std::endl;
    switch(type) {
		case render_type::triangle:
		{
			const data_geometry *vertexArray[3];
			data_geometry temp[3];
			data_vertex a, b, c;
	
			a.data = &state.vertex_data[0];
			b.data = &state.vertex_data[1 * state.floats_per_vertex];
			c.data = &state.vertex_data[2 * state.floats_per_vertex];
	
			temp[0].data = a.data;
			temp[1].data = b.data;
			temp[2].data = c.data;
	
			state.vertex_shader(a, temp[0], state.uniform_data);
			state.vertex_shader(b, temp[1], state.uniform_data);
			state.vertex_shader(c, temp[2], state.uniform_data);
	
			for(int i = 0; i < 3; i++) {
				vertexArray[i] = &temp[i];
			}	
			
			rasterize_triangle(state, vertexArray);
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
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    //std::cout<<"TODO: implement rasterization"<<std::endl;
    
    float aX = 0, aY = 0, bX = 0, bY = 0, cX = 0, cY = 0;
    float alpha = 0, beta = 0, gamma = 0;
    float alpha0 = 0, alpha1 = 0, alpha2 = 0, beta0 = 0, beta1 = 0, beta2 = 0, gamma0 = 0, gamma1 = 0, gamma2 = 0;
    
    auto w = state.image_width;
    auto h = state.image_height;
    
    aX = (w/2) * (*in)[0].gl_Position[0] + ((w/2) - (1/2));
    bX = (w/2) * (*in)[1].gl_Position[0] + ((w/2) - (1/2));
    cX = (w/2) * (*in)[2].gl_Position[0] + ((w/2) - (1/2));
    aY = (h/2) * (*in)[0].gl_Position[1] + ((h/2) - (1/2));
    bY = (h/2) * (*in)[1].gl_Position[1] + ((h/2) - (1/2));
    cY = (h/2) * (*in)[2].gl_Position[1] + ((h/2) - (1/2));
    
    float area = 0.5f * (((bX * cY) - (cX * bY)) - ((aX * cY) - (cX * aY)) - ((aX * bY) - (bX * aY)));
    
    alpha0 = (bX * cY - cX * bY);
    alpha1 = (bY - cY);
    alpha2 = (cX - bX);
    
    beta0 = (cX * aY - aX * cY);
    beta1 = (cY - aY);
    beta2 = (aX - cX);
	
    gamma0 = (aX * bY - bX * aY);
    gamma1 = (aY - bY);
    gamma2 = (bX - aX);
    
    for(auto j = 0; j < state.image_height; ++j) {
		for(auto i = 0; i < state.image_width; ++i) {
			alpha = (0.5f * (alpha0 + alpha1 * i + alpha2 * j))/area;
			beta = (0.5f * (beta0 + beta1 * i + beta2 * j))/area;
			gamma = (0.5f * (gamma0 + gamma1 * i + gamma2 * j))/area;
			if(alpha >= 0 && beta >= 0 && gamma >= 0) {
				state.image_color[i + j * w] = make_pixel(255.0, 255.0, 255.0);
			}
		}
     }
    
}

