//This header file includes the definition of shader.
//And also a function "Triangle_Traversal" used to test if a pixel in a traiangle.

#ifndef __SHADER_H__
#define __SHADER_H__

#include"tgaimage.h"
#include"geometry.h"
#include"model.h"
#include"scene.h"
#include<cmath>

class IShader {
public:
    virtual ~IShader();
    virtual vec4 vertex(int iface, int nthvert, Scene &scene) = 0;
    virtual bool fragment(vec3 bar, TGAColor &color) = 0;
};

void triangle(vec4 *pts, IShader &shader, TGAImage &image, TGAImage &zbuffer);

class GouraudShader : public IShader {
    vec3 varying_intensity; // written by vertex shader, read by fragment shader
    ~GouraudShader() {}

    virtual vec4 vertex(int iface, int nthvert, Scene* scene, Model* model) {
    	vec3 ver3 = model->get_ver(iface, nthvert);
    	vec4 ver4 = ver3.extend_4();
    	vec4 gl_Vertex = scene->View*scene->Projection*scene->Modelview*ver4;
    	vec3 normal = model->get_nor(iface, nthvert);
    	varying_intensity[nthvert] = std::max(0, dot(normal, scene->light_dir));
    	return gl_Vertex;
    }

    virtual bool fragment(vec3 bar, TGAColor &color) {
        float intensity = dot(varying_intensity, bar);   // interpolate intensity for the current pixel
        color = TGAColor(255, 255, 255)*intensity; // well duh
        return false;                              // no, we do not discard this pixel
    }
};

vec3 barycentric(vec3 A, vec3 B, vec3 C, vec3 P);
void triangle(vec4 *gl_pts, IShader &shader, TGAImage &image, TGAImage &zbuffer);

#endif