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

    int numTriangles = state.num_vertices / 3.0;

    	switch(type) {
		case render_type::triangle:
		{
			const data_geometry *vertexArray[3];
			data_geometry temp[3];
			data_vertex a, b, c;
			int index = 0;
		
			for(int i = 0; i < numTriangles; i++) {
				a.data = &state.vertex_data[index];
				index += state.floats_per_vertex;
				b.data = &state.vertex_data[index];
				index += state.floats_per_vertex;
				c.data = &state.vertex_data[index];
				index += state.floats_per_vertex;
			
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
    
    float aX = 0, aY = 0, bX = 0, bY = 0, cX = 0, cY = 0;
    float alpha = 0, beta = 0, gamma = 0;
    
    auto w = state.image_width;
    auto h = state.image_height;
    
    aX = (w/2) * in[0]->gl_Position[0] / in[0]->gl_Position[3] + ((w/2) - (1/2));
    aY = (h/2) * in[0]->gl_Position[1] / in[0]->gl_Position[3] + ((h/2) - (1/2));
    bX = (w/2) * in[1]->gl_Position[0] / in[1]->gl_Position[3] + ((w/2) - (1/2));
    bY = (h/2) * in[1]->gl_Position[1] / in[1]->gl_Position[3] + ((h/2) - (1/2));
    cX = (w/2) * in[2]->gl_Position[0] / in[2]->gl_Position[3] + ((w/2) - (1/2));
    cY = (h/2) * in[2]->gl_Position[1] / in[2]->gl_Position[3] + ((h/2) - (1/2));
    
    for(auto i = 0; i < w; i++) {
	float area = (cX - aX) * (bY - aY) - (cY - aY) * (bX - aX);
	for(auto j = 0; j < h; j++) {
		alpha = (i - bX) * (cY - bY) - (j - bY) * (cX - bX); 
		beta = (i - cX) * (aY - cY) - (j - cY) * (aX - cX);
		gamma = (i - aX) * (bY - aY) - (j - aY) * (bX - aX);

		alpha /= area;
		beta /= area;
		gamma /= area;
		
		float fragData[MAX_FLOATS_PER_VERTEX];
		data_fragment fragment{fragData};

		if(alpha >= 0 && beta >= 0 && gamma >= 0) {
			for(int k = 0; k < state.floats_per_vertex; k++) {
				switch(state.interp_rules[k]) {
					case interp_type::flat:{
						fragment.data[k] = in[0]->data[k];
						break;
					}
				
					case interp_type::smooth:{
						float const_k = (alpha / in[0]->gl_Position[3]) + (beta / in[1]->gl_Position[3]) + (gamma / in[2]->gl_Position[3]);
						float alpha1 = (alpha / in[0]->gl_Position[3]) / const_k;
						float beta1 = (beta / in[1]->gl_Position[3]) / const_k;
						float gamma1 = (gamma / in[2]->gl_Position[3]) / const_k;
						fragment.data[k] = alpha1 * in[0]->data[k] + beta1 * in[1]->data[k] + gamma1 * in[2]->data[k];
						break;
					}		

					case interp_type::noperspective:{
						fragment.data[k] = alpha * in[0]->data[k] + beta * in[1]->data[k] + gamma * in[2]->data[k];
						break;
					}		
			
					case interp_type::invalid:{	
						break;	
					}
				}
    			}

			data_output out;
			state.fragment_shader(fragment, out, state.uniform_data);
			state.image_color[i + j * w] = make_pixel(out.output_color[0] * 255.0, out.output_color[1] * 255.0, out.output_color[2] * 255.0);
				
		}
	}
    }
}
