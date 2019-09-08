#include<cmath>
#include<limits>
#include"shader.h"
#include"geometry.h"
#include"model.h"
#include"scene.h"
#include"tgaimage.h"


vec3 barycentric(vec3 A, vec3 B, vec3 C, vec3 P){
	vec3 s1(B[0]-A[0], C[0]-A[0], A[0]-P[0]);
	vec3 s2(B[1]-A[1], C[1]-A[1], A[1]-P[1]);
	vec3 res = cross(s1, s2);
	if(std::abs(res[2])>1e-2) return vec3(1-(res[0]+res[1])/res[2], res[0]/res[2], res[1]/res[2]);
	else return vec3(-1,1,1);
}

void triangle(vec4 *gl_pts, IShader &shader, TGAImage &image, TGAImage &zbuffer) {
    vec3 bboxmin( std::numeric_limits<float>::max(),  std::numeric_limits<float>::max(), 0);
    vec3 bboxmax(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), 0);
    for (int i=0; i<3; i++) {
        for (int j=0; j<2; j++) {
            bboxmin[j] = std::min(bboxmin[j], gl_pts[i][j]/gl_pts[i][3]);
            bboxmax[j] = std::max(bboxmax[j], gl_pts[i][j]/gl_pts[i][3]);
        }
    }
    vec3 A(gl_pts[0][0], gl_pts[0][1], 0);
    vec3 B(gl_pts[1][0], gl_pts[1][1], 0);
    vec3 C(gl_pts[2][0], gl_pts[2][1], 0);
    TGAColor color;
    for (int i=bboxmin.x(); i<=bboxmax.x(); i++) {
        for (int j=bboxmin.y(); j<=bboxmax.y(); j++) {
        	vec3 P(i, j, 0);
            vec3 c = barycentric(A, B, C, P);
            float z = gl_pts[0][2]*c.x() + gl_pts[1][2]*c.y() + gl_pts[2][2]*c.z();
            float w = gl_pts[0][3]*c.x() + gl_pts[1][3]*c.y() + gl_pts[2][3]*c.z();
            int frag_depth = std::max(0, std::min(255, int(z/w+.5)));
            if (c.x()<0 || c.y()<0 || c.z()<0 || zbuffer.get(i, j)[0]>frag_depth) continue;
            bool discard = shader.fragment(c, color);
            if (!discard) {
                zbuffer.set(i, j, TGAColor(frag_depth));
                image.set(i, j, color);
            }
        }
    }
}
