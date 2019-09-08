//This header file includes the definition of Vec2, Vec3, Matrix.

#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__
#include<math.h>
#include<iostream>
#include<stdlib.h>

class vec4 {
public:
	vec4() {}
	vec4(float e0,float e1,float e2,float e3){
		e[0]=e0;
		e[1]=e1;
		e[2]=e2;
		e[3]=e3;
	}

	inline float x() const {return e[0];}
	inline float y() const {return e[1];}
	inline float z() const {return e[2];}
	inline float h() const {return e[3];}

	inline void set_x(float x) {e[0] = x;}
	inline void set_y(float y) {e[1] = y;}
	inline void set_z(float z) {e[2] = z;}
	inline void set_h(float h) {e[3] = h;}

	inline void set(int idx, float x){
		e[idx] = x;
	}

	inline const vec4& operator+() const {return *this;}
	inline vec4 operator-() const {return vec4(-e[0],-e[1],-e[2],-e[3]);}
	inline float operator[](int i) const {return e[i];}
	inline float& operator[](int i) {return e[i];}

	inline vec4& operator+=(const vec4 &v2);
	inline vec4& operator-=(const vec4 &v2);
	inline vec4& operator*=(const vec4 &v2);
	inline vec4& operator/=(const vec4 &v2);
	inline vec4& operator*=(const float t);
	inline vec4& operator/=(const float t);
	inline void  operator=(const vec4 v);

	inline float length() const{
		return sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2] + e[3]*e[3]);}
	inline float squared_length() const {
		return e[0]*e[0] + e[1]*e[1] + e[2]*e[2] + e[3]*e[3];}
	inline void make_unit_vector();
	inline void standardization();

	float e[4];
};


class vec3 {
public:
	vec3() {}
	vec3(float e0,float e1,float e2){
		e[0]=e0;
		e[1]=e1;
		e[2]=e2;
	}

	inline float x() const {return e[0];}
	inline float y() const {return e[1];}
	inline float z() const {return e[2];}
	inline float r() const {return e[0];}
	inline float g() const {return e[1];}
	inline float b() const {return e[2];}

	inline void set_x(float x) {e[0] = x;}
	inline void set_y(float y) {e[1] = y;}
	inline void set_z(float z) {e[2] = z;}

	inline const vec3& operator+() const {return *this;}
	inline vec3 operator-() const {return vec3(-e[0],-e[1],-e[2]);}
	inline float operator[](int i) const {return e[i];}
	inline float& operator[](int i) {return e[i];}

	inline vec3& operator+=(const vec3 &v2);
	inline vec3& operator-=(const vec3 &v2);
	inline vec3& operator*=(const vec3 &v2);
	inline vec3& operator/=(const vec3 &v2);
	inline vec3& operator*=(const float t);
	inline vec3& operator/=(const float t);
	inline void  operator=(const vec3 v);

	inline float length() const{
		return sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);}
	inline float squared_length() const {
		return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];}
	inline void make_unit_vector();
	inline vec4 extend_4(){
		return vec4(e[0],e[1],e[2],1);
	}

	float e[3];
};


class Matrix44{
public:
	float matrix[4][4];

	Matrix44(){
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				if(i==j) matrix[i][j]=1;
				else matrix[i][j]=0;
			}
		}
	}

	void set(int i, int j, int x){
		matrix[i][j] = x;
	}

	float* operator[](int idx) {
		return matrix[idx];
	}

	void operator=(Matrix44 mat){
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				matrix[i][j] = mat[i][j];
			}
		}
		return;
	}
	
	void identity(){
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				if(i==j) matrix[i][j]=1;
				else matrix[i][j]=0;
			}
		}
	}
};

//operations about vec3
inline vec3 operator+(const vec3 &v1, const vec3 &v2);
inline vec3 operator-(const vec3 &v1, const vec3 &v2);
inline vec3 operator*(const vec3 &v1, const vec3 &v2);
inline vec3 operator/(const vec3 &v1, const vec3 &v2);
inline vec3 operator*(float t,const vec3 &v);
inline vec3 operator/(vec3 v, float t);
inline vec3 operator*(const vec3 &v, float t);
inline float distance(const vec3 &v1, const vec3 &v2);
inline float dot(const vec3 &v1, const vec3 &v2);
inline vec3 cross(const vec3 &v1,const vec3 &v2);
inline vec3 unit_vector(vec3 v);

//operations about vec4
inline vec4 operator+(const vec4 &v1, const vec4 &v2);
inline vec4 operator-(const vec4 &v1, const vec4 &v2);
inline vec4 operator*(const vec4 &v1, const vec4 &v2);
inline vec4 operator/(const vec4 &v1, const vec4 &v2);
inline vec4 operator*(float t,const vec4 &v);
inline vec4 operator/(vec4 v, float t);
inline vec4 operator*(const vec4 &v, float t);
inline float dot(const vec4 &v1, const vec4 &v2);
inline vec4 unit_vector(vec4 v);

//operations about multiply
vec4 operator*( Matrix44 mat, const vec4 v);
Matrix44 operator*( Matrix44 mat1, Matrix44 mat2);

#endif