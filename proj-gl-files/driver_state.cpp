#include "driver_state.h"
#include <cstring>
#include <limits>

using namespace std;

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete[] image_color;
    delete[] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state &state, int width, int height)
{
    state.image_width = width;
    state.image_height = height;
    state.image_color = 0;
    state.image_depth = 0;
    // std::cout << "TODO: allocate and initialize state.image_color and state.image_depth." << std::endl;

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
void render(driver_state &state, render_type type)
{
    // std::cout << "TODO: implement rendering." << std::endl;
    int numTriangles = state.num_vertices / 3;
    switch (type)
    {
    case render_type::triangle:
    {
        const data_geometry *d_GeometryArr[3];
        data_geometry g[3];
        data_vertex v[3];

        int k = 0;

        for (int i = 0; i < numTriangles; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                v[j].data = &state.vertex_data[k];
                g[j].data = v[j].data;
                state.vertex_shader(v[j], g[j], state.uniform_data);
                d_GeometryArr[j] = &g[j];
                k += state.floats_per_vertex;
            }
            clip_triangle(state, d_GeometryArr, 0);
        }
        break;
    }
    case render_type::invalid:
        break;
    case render_type::indexed:
    {
      const data_geometry *d_GeometryArr[3];
      data_geometry g[3];
      data_vertex v[3];

      for (int i = 0; i < state.num_triangles*3; i +=3)
      {
        for(int j = 0; j < 3; ++j) {
          v[j].data = &state.vertex_data[state.index_data[i + j]*state.floats_per_vertex];
          g[j].data = v[j].data;
          state.vertex_shader(v[j], g[j], state.uniform_data);
          d_GeometryArr[j] = &g[j];
        }
        clip_triangle(state, d_GeometryArr, 0);
      }
      break;
    }
    case render_type::fan:
    {
      const data_geometry *d_GeometryArr[3];
      data_geometry g[3];
      data_vertex v[3];

      int k = state.floats_per_vertex;

      for (int i = 0; i < state.num_vertices; ++i)
      {
          for (int j = 0; j < 3; ++j)
          {
            if(j == 0) {
              v[j].data = &state.vertex_data[0];
            }
            else {
              v[j].data = &state.vertex_data[k];
              if(j == 1) {
                k += state.floats_per_vertex;
              }
            }
            g[j].data = v[j].data;
            state.vertex_shader(v[j], g[j], state.uniform_data);
            d_GeometryArr[j] = &g[j];
          }
          clip_triangle(state, d_GeometryArr, 0);
      }
      break;
  	}
    case render_type::strip:
		{
      const data_geometry *d_GeometryArr[3];
      data_geometry g[3];
      data_vertex v[3];

      int k = 0;

      for (int i = 0; i < state.num_vertices - 2; ++i)
      {
          if(k != 0) {k -= (2*state.floats_per_vertex);}
          for (int j = 0; j < 3; ++j)
          {
              v[j].data = &state.vertex_data[k];
              g[j].data = v[j].data;
              state.vertex_shader(v[j], g[j], state.uniform_data);
              d_GeometryArr[j] = &g[j];
              k += state.floats_per_vertex;
          }
          clip_triangle(state, d_GeometryArr, 0);
      }
			break;
		}
    default:
        break;
    }
}

void interpolate(driver_state &state, data_geometry &new_dg, const data_geometry dg[2], const data_geometry *in[3], const int &v_num, const float &beta)
{
    for (int i = 0; i < state.floats_per_vertex; ++i)
        switch (state.interp_rules[i])
        {
        case interp_type::flat:
            new_dg.data[i] = in[v_num]->data[i];
            break;
        case interp_type::smooth:
            new_dg.data[i] = beta * dg[0].data[i] + (1 - beta) * dg[1].data[i];
            break;
        case interp_type::noperspective:
        {
            float alpha = beta * dg[0].gl_Position[3] / (beta * dg[0].gl_Position[3] + (1 - beta) * dg[1].gl_Position[3]);
            new_dg.data[i] = alpha * dg[0].data[i] + (1 - alpha) * dg[1].data[i];
            break;
        }
        default:
            break;
        }
}

float compute_beta(const vec4 &s_v, const vec4 &e_v, const int &coord)
{
    return (-e_v[3] - e_v[coord]) / (s_v[coord] + s_v[3] - e_v[3] - e_v[coord]);
}

vec4 compute_intersection(const vec4 &s_v, const vec4 &e_v, const float &beta)
{
    return beta * s_v + (1 - beta) * e_v;
}

void clip_triangle(driver_state &state, const data_geometry *in[3], int face)
{
    if (face == 1)
    {
        rasterize_triangle(state, in);
        return;
    }

    vec4 a = in[0]->gl_Position, b = in[1]->gl_Position, c = in[2]->gl_Position;
    const data_geometry *_in[3] = {in[0], in[1], in[2]};
    data_geometry new_dg[3];
    data_geometry dg[2];
    float beta0, beta1;
    vec4 point0, point1;

    if (a[2] < -a[3] && b[2] < -b[3] && c[2] < -c[3])
        return;
    else
        if (a[2] < -a[3] && b[2] >= -b[3] && c[2] >= -c[3])
        {
            beta0 = compute_beta(a, b, 2);
            beta1 = compute_beta(c, a, 2);
            point0 = compute_intersection(a, b, beta0);
            point1 = compute_intersection(c, a, beta1);

            new_dg[0].data = new float[state.floats_per_vertex];
            new_dg[0].gl_Position = point1;
            new_dg[1] = *in[1];
            new_dg[2] = *in[2];
            dg[0] = *in[2];
            dg[1] = *in[0];
            interpolate(state, new_dg[0], dg, in, 0, beta1);
            _in[0] = &new_dg[0];
            _in[1] = &new_dg[1];
            _in[2] = &new_dg[2];

            clip_triangle(state, _in, face + 1);

            new_dg[0].data = new float[state.floats_per_vertex];
            new_dg[0].gl_Position = point0;
            new_dg[1] = *in[1];
            new_dg[2].data = new float[state.floats_per_vertex];
            new_dg[2].gl_Position = point1;
            dg[0] = *in[0];
            dg[1] = *in[1];
            interpolate(state, new_dg[0], dg, in, 0, beta0);
            dg[0] = *in[2];
            dg[1] = *in[0];
            interpolate(state, new_dg[2], dg, in, 2, beta1);
            _in[0] = &new_dg[0];
            _in[1] = &new_dg[1];
            _in[2] = &new_dg[2];
        }

        clip_triangle(state, _in, face + 1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state &state, const data_geometry *in[3])
{
    //std::cout << "TODO: implement rasterization" << std::endl;
    vec3 w = {in[0]->gl_Position[3], in[1]->gl_Position[3], in[2]->gl_Position[3]};
    vec3 x = {in[0]->gl_Position[0] / w[0], in[1]->gl_Position[0] / w[1], in[2]->gl_Position[0] / w[2]};
    vec3 y = {in[0]->gl_Position[1] / w[0], in[1]->gl_Position[1] / w[1], in[2]->gl_Position[1] / w[2]};
    vec3 z = {in[0]->gl_Position[2] / w[0], in[1]->gl_Position[2] / w[1], in[2]->gl_Position[2] / w[2]};

    vec2 v0 = {(float)((state.image_width / 2) * x[0]) + (float)((state.image_width / 2) - 0.5),
               (float)((state.image_height / 2) * y[0]) + (float)((state.image_height / 2) - 0.5)};
    vec2 v1 = {(float)((state.image_width / 2) * x[1]) + (float)((state.image_width / 2) - 0.5),
               (float)((state.image_height / 2) * y[1]) + (float)((state.image_height / 2) - 0.5)};
    vec2 v2 = {(float)((state.image_width / 2) * x[2]) + (float)((state.image_width / 2) - 0.5),
               (float)((state.image_height / 2) * y[2]) + (float)((state.image_height / 2) - 0.5)};

    vec3 area = {(v2[0] - v0[0]) * (v1[1] - v0[1]) - (v2[1] - v0[1]) * (v1[0] - v0[0]), 0, 0};

    for (int i = 0; i < state.image_width; ++i)
        for (int j = 0; j < state.image_height; ++j)
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

                if (state.image_depth[i + j * state.image_width] > zDepth)
                {
                    data_fragment d_frag;
                    d_frag.data = new float[state.floats_per_vertex];
                    data_output d_output;

                    for (int k = 0; k < state.floats_per_vertex; ++k)
                        if (state.interp_rules[k] == interp_type::flat)
                            d_frag.data[k] = in[0]->data[k];
                        else if (state.interp_rules[k] == interp_type::noperspective)
                            d_frag.data[k] = alpha * in[0]->data[k] + beta * in[1]->data[k] + gamma * in[2]->data[k];
                        else if (state.interp_rules[k] == interp_type::smooth)
                        {
                            float c = (alpha / w[0]) + (beta / w[1]) + (gamma / w[2]);
                            a_prime = (alpha / w[0]) / c;
                            b_prime = (beta / w[1]) / c;
                            g_prime = (gamma / w[2]) / c;
                            d_frag.data[k] = (a_prime * in[0]->data[k]) + (b_prime * in[1]->data[k]) + (g_prime * in[2]->data[k]);
                        }

                    state.fragment_shader(d_frag, d_output, state.uniform_data);

                    float r = d_output.output_color[0] * 255;
                    float g = d_output.output_color[1] * 255;
                    float b = d_output.output_color[2] * 255;

                    state.image_color[i + j * state.image_width] = make_pixel(r, g, b);
                    state.image_depth[i + j * state.image_width] = zDepth;
                }
            }
        }
}
