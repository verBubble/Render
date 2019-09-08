#include<vector>
#include<iostream>

#include"geometry.h"
#include"model.h"
#include"scene.h"
#include"shader.h"
#include"tgaimage.h"

Model *model     = NULL;
const int width  = 800;
const int height = 800;

vec3 light_dir(1,1,1);
vec3      eye(0,-1,3);
vec3    center(0,0,0);
vec3        up(0,1,0);

int main(int argc, char** argv) {
    if (2==argc) {
        model = new Model(argv[1]);
    } else {
        model = new Model("obj/african_head.obj");
    }

    Scene* scene;
    scene = new Scene(up, center, eye, light_dir, width, height);

    GouraudShader shader;
    TGAImage image  (width, height, TGAImage::RGB);
    TGAImage zbuffer(width, height, TGAImage::GRAYSCALE);

    for (int i=0; i<model->nfaces(); i++) {
        Vec4f screen_coords[3];
        for (int j=0; j<3; j++) {
            screen_coords[j] = shader.vertex(i, j, scene, model);
        }
        triangle(screen_coords, shader, image, zbuffer);
    }

    image.  flip_vertically(); // to place the origin in the bottom left corner of the image
    zbuffer.flip_vertically();
    image.  write_tga_file("output.tga");
    zbuffer.write_tga_file("zbuffer.tga");

    delete model;
    return 0;
}