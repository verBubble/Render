#include <vector>
#include <cmath>
#include <iostream>
#include "tgaimage.h"
#include "model.h"
#include "geometry.h"

const TGAColor white = TGAColor(255, 255, 255, 255);
const TGAColor red   = TGAColor(255, 0,   0,   255);
Model *model = NULL;
const int width  = 800;
const int height = 800;
Vec3f dir = Vec3f(0,0,-1);
float z_buffer[width][height];
TGAImage Texture;

void line(int x0, int y0, int x1, int y1, TGAImage &image, TGAColor color) { 
    bool steep = false; 
    if (std::abs(x0-x1)<std::abs(y0-y1)) { 
        std::swap(x0, y0); 
        std::swap(x1, y1); 
        steep = true; 
    } 
    if (x0>x1) { 
        std::swap(x0, x1); 
        std::swap(y0, y1); 
    } 
    int dx = x1-x0; 
    int dy = y1-y0; 
    int derror2 = std::abs(dy)*2; //notice it's abs !!! 
    int error2 = 0; 
    int y = y0; 
    for (int x=x0; x<=x1; x++) { 
        if (steep) { 
            image.set(y, x, color); 
        } else { 
            image.set(x, y, color); 
        } 
        error2 += derror2; 
        if (error2 > dx) { 
            y += (y1>y0?1:-1); 
            error2 -= dx*2; 
        } 
    } 
} 

inline 
float edgeFunction(const float* a, const float* b, const Vec3f &c) 
{ return (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0]); } 

inline 
float edgeFunction(const Vec3f &a, const Vec3f &b, const Vec3f &c) 
{ return (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0]); } 

inline 
float edgeFunction(const float* a, const float* b, const float* c) 
{ return (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0]); } 

//this function is used for coloring inside a triangle, via colors stored in vertex.
//maybe vertexs' attribute can be packed in a class(will be optimizing ...)
void triangle(Vec3f t0, Vec3f t1, Vec3f t2, TGAImage &image, Vec3f color0, Vec3f color1, Vec3f color2,
    Vec2f vt0, Vec2f vt1, Vec2f vt2){
    int w = image.get_width();
    int h = image.get_height();

    float v0[3];
    float v1[3];
    float v2[3];

    int TexW = Texture.get_width();
    int TexH = Texture.get_height();

    TGAColor vtc0 = Texture.get(int(vt0[0]*TexW), int(vt0[1]*TexH));
    TGAColor vtc1 = Texture.get(int(vt1[0]*TexW), int(vt1[1]*TexH));
    TGAColor vtc2 = Texture.get(int(vt2[0]*TexW), int(vt2[1]*TexH));

    //std::cout<<"vtc0"<<vtc0.get_r()<<std::endl;

/*    //project triangle into screen
    v0[0] = t0[0]/t0[2]; v0[1] = t0[1]/t0[2];
    v1[0] = t1[0]/t1[2]; v1[1] = t1[1]/t1[2];
    v2[0] = t2[0]/t2[2]; v2[1] = t2[1]/t2[2];
*/
    //convert from screen space to NDC then raster
    v0[0] = (1+t0[0])*0.5*w; v0[1] = (1+t0[1])*0.5*h;
    v1[0] = (1+t1[0])*0.5*w; v1[1] = (1+t1[1])*0.5*h;
    v2[0] = (1+t2[0])*0.5*w; v2[1] = (1+t2[1])*0.5*h;

/*  code about color via lighting
    Vec3f n = (t2 - t0)^(t1 - t0);
    n.normalize();
    float intensity = n*dir;

    if(intensity<0) return;
    float c0[3] = {intensity, intensity, intensity};
    float c1[3] = {intensity, intensity, intensity};
    float c2[3] = {intensity, intensity, intensity};
*/

    float c0[3] = {(float)vtc0.get_r(), (float)vtc0.get_g(), (float)vtc0.get_b()};
    float c1[3] = {(float)vtc1.get_r(), (float)vtc1.get_g(), (float)vtc1.get_b()};
    float c2[3] = {(float)vtc2.get_r(), (float)vtc2.get_g(), (float)vtc2.get_b()};

    c0[0] /= t0[2]; c0[1] /= t0[2]; c0[2] /= t0[2];
    c1[0] /= t1[2]; c1[1] /= t1[2]; c1[2] /= t1[2];
    c2[0] /= t2[2]; c2[1] /= t2[2]; c2[2] /= t2[2];


/*
   //divide vertex attribute by z-coordinate
    float c0[3] = {color0[0]/t0[2], color0[1]/t0[2], color0[2]/t0[2]};
    float c1[3] = {color1[0]/t1[2], color1[1]/t1[2], color1[2]/t1[2]};
    float c2[3] = {color2[0]/t2[2], color2[1]/t2[2], color2[2]/t2[2]};
*/

/* 
    float c0[3] = {color0[0], color0[1], color0[2]};
    float c1[3] = {color1[0], color1[1], color1[2]};
    float c2[3] = {color2[0], color2[1], color2[2]};
*/
    //pre-compute z-coordinate
    float z1 = (float)1./(float)t0[2]; 
    float z2 = (float)1./(float)t1[2]; 
    float z3 = (float)1./(float)t2[2];

    float area = edgeFunction(v0, v1, v2);

    for(int j=0; j<h; j++){
        for(int i=0; i<w; i++){
            Vec3f p = Vec3f(i, j, 0);
            float w0 = edgeFunction(v2, v1, p); 
            float w1 = edgeFunction(v0, v2, p); 
            float w2 = edgeFunction(v1, v0, p); 

            if(w0>=0 && w1>=0 && w2>=0){

                w0 /= area;
                w1 /= area;
                w2 /= area;

                float oneOverZ = z1 * w0 + z2 * w1 + z3 * w2;
                float z = 1/oneOverZ;

                if(z_buffer[i][j]>z){
                    z_buffer[i][j]=z; 
                    float r = w0*c0[0] + w1*c1[0] + w2*c2[0];
                    float g = w0*c0[1] + w1*c1[1] + w2*c2[1];
                    float b = w0*c0[2] + w1*c1[2] + w2*c2[2];
                    r *= z;
                    g *= z;
                    b *= z;

                    TGAColor c = TGAColor(r, g, b, 255);
                    image.set(i, j, c);   
                }        
            }
        }
    }
}


int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    for(int j=0;j<height;j++){
        for(int i=0;i<width;i++){
            z_buffer[i][j]=0x7fffffff;
        }
    }

    Texture.read_tga_file("african_head_diffuse.tga");
    Texture.flip_vertically();
    TGAImage image(width, height, TGAImage::RGB);
    for (int i=0; i<model->nfaces(); i++) {
        std::vector<int> face = model->face(i);
        std::vector<int> text = model->texture(i);
        Vec3f v0 = model->vert(face[0]);
        Vec3f v1 = model->vert(face[1]);
        Vec3f v2 = model->vert(face[2]);
        Vec2f vt0 = model->vt(text[0]);
        Vec2f vt1 = model->vt(text[1]);
        Vec2f vt2 = model->vt(text[2]);

        Vec3f c = Vec3f(rand()/double(RAND_MAX), rand()/double(RAND_MAX), rand()/double(RAND_MAX));
        Vec3f c1 = Vec3f(rand()/double(RAND_MAX), rand()/double(RAND_MAX), rand()/double(RAND_MAX));
        Vec3f c2 = Vec3f(rand()/double(RAND_MAX), rand()/double(RAND_MAX), rand()/double(RAND_MAX));


        triangle(v0, v1, v2, image, c, c1, c2, vt0, vt1, vt2);

    }

    std::vector<int> face = model->face(5);
    Vec3f v1 = model->vert(face[1]);
    std::vector<int> text = model->texture(5);
   // std::cout<<text[1]<<std::endl;
    Vec2f vt = model->vt(text[2]);
    std::cout<<"test_texture_ver"<<vt<<std::endl;
    std::cout<<"test_ver"<<v1<<std::endl;
    image.flip_vertically();
    image.write_tga_file("output.tga");
    delete model;
    return 0;
}

