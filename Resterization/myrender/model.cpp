#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include "model.h"

Model::Model(const char *filename) : verts_(), vts_(), faces_(), vertextures_() {
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
            Vec3f v;
            for (int i=0;i<3;i++) iss >> v.raw[i];
            verts_.push_back(v);
        } else if (!line.compare(0, 2, "f ")) {
            std::vector<int> f;
            std::vector<int> t;  //texture
            int itrash, idx, tidx;
            iss >> trash;
            while (iss >> idx >> trash >> tidx >> trash >> itrash) {
                idx--; // in wavefront obj all indices start at 1, not zero
                tidx--;
                f.push_back(idx);
                t.push_back(tidx);
            }
            faces_.push_back(f);
            vertextures_.push_back(t);
        } else if (!line.compare(0, 3, "vt ")){
            iss>>trash>>trash;
            Vec2f vt;
            iss>>vt.raw[0]>>vt.raw[1]>>trash;
            vts_.push_back(vt);
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