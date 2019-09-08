//This header file includes the class model, which is defined to operate the .obj file, 
//and write data to data structure (triangle, vertexAttribute).

#ifndef __MODEL_H__
#define __MODEL_H__

#include <vector>
#include "geometry.h"

class Model {
private:
	std::vector<vec3> verts_;
	std::vector<vec3> vts_;   //just ignore the 3th data.
	std::vector<vec3> vns_;
	std::vector<std::vector<int> > faces_;
	std::vector<std::vector<int> > vertextures_;
	std::vector<std::vector<int> > vernormals_;
public:
	Model(const char *filename);
	~Model();
	int nverts();
	int nfaces();
	vec3 vert(int i);
	vec3 vt(int i);
	vec3 vn(int i);
	std::vector<int> face(int idx);
	std::vector<int> texture(int idx);  //the ith triangle's texture.
	std::vector<int> normal(int idx);

	vec3 get_ver(int iface, int idx);
	vec3 get_nor(int iface, int idx);
	vec3 get_text(int iface, int idx);
};

#endif //__MODEL_H__