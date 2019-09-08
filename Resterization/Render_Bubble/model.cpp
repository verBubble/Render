#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

Model::Model(const char *filename){
    std::ifstream in;
    in.open (filename, std::ifstream::in);
    if (in.fail()) return;
    std::string line;
    while (!in.eof()) {
        std::getline(in, line);
        std::istringstream iss(line.c_str());
        char trash;
        if (!line.compare(0, 2, "v ")) {
            iss >> trash;
            Vec3 v;
            float t0,t1,t2;
            iss>>t0>>t1>>t2;
            v.set_x(t0);
            v.set_y(t1);
            v.set_z(t2);
            verts_.push_back(v);
        } else if (!line.compare(0, 2, "f ")) {
            std::vector<int> f;
            std::vector<int> t;  //texture
            std::vector<int> n;  //normal
            int nidx, idx, tidx;
            iss >> trash;
            while (iss >> idx >> trash >> tidx >> trash >> nidx) {
                idx--; // in wavefront obj all indices start at 1, not zero
                tidx--;
                nidx--;
                f.push_back(idx);
                t.push_back(tidx);
                n.push_back(nidx);
            }
            faces_.push_back(f);
            vertextures_.push_back(t);
            vernormals_.push_back(n);
        } else if (!line.compare(0, 3, "vt ")){
            iss>>trash>>trash;
            Vec3 vt;
            float temp,temp1;
            iss>>temp>>temp1>>trash;
            vt.set_x(temp);
            vt.set_y(temp1);
            vts_.push_back(vt);
        } else if (!line.compare(0, 3, "vn ")){
            iss>>trash>>trash;
            Vec3 vn;
            float temp,temp1,temp2;
            iss>>temp>>temp1>>temp2>>trash;
            vn.set_x(temp);
            vn.set_y(temp1);
            vn.set_z(temp2);
            vns_.push_back(vn);
        }
    }
    std::cerr << "# v# " << verts_.size() << " f# "  << faces_.size() << std::endl;
}

int Model::nverts() {
    return (int)verts_.size();
}

int Model::nfaces() {
    return (int)faces_.size();
}

std::vector<int> Model::face(int idx) {
    return faces_[idx];
}

Vec3 Model::vert(int i) {
    return verts_[i];
}

std::vector<int> Model::texture(int idx) {
    return vertextures_[idx];
}

std::vector<int> Model::normal(int idx) {
    return vernormals_[idx];
}

Vec3 Model::vt(int i) {
    return vts_[i];
}

Vec3 Model::vn(int i) {
    return vns_[i];
}

vec3 Model::get_ver(int iface, int idx){
    return vert(face(iface)[idx]);
}

vec3 Model::get_nor(int iface, int idx){
    return vn(normal(iface)[idx]);
}

vec3 Model::get_text(int iface, int idx){
    return vt(texture(iface)[idx]);
}