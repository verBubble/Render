#include<math.h>
#include<stdlib.h>
#include<iostream>
#include"geometry.h"

//operations about vec3

inline void vec3::make_unit_vector(){
	float k=1.0 /sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);
	e[0]*=k; e[1]*=k;e[2]*=k;
}

inline vec3 operator+(const vec3 &v1, const vec3 &v2){
	return vec3(v1.e[0] + v2.e[0],v1.e[1]+v2.e[1], v1.e[2]+v2.e[2]);
}

inline vec3 operator-(const vec3 &v1, const vec3 &v2){
	return vec3(v1.e[0]-v2.e[0],v1.e[1]-v2.e[1], v1.e[2]-v2.e[2]);
}

inline vec3 operator*(const vec3 &v1, const vec3 &v2){
	return vec3(v1.e[0]*v2.e[0],v1.e[1]*v2.e[1], v1.e[2]*v2.e[2]);
}

inline vec3 operator/(const vec3 &v1, const vec3 &v2){
	return vec3(v1.e[0]/v2.e[0],v1.e[1]/v2.e[1], v1.e[2]/v2.e[2]);
}

inline vec3 operator*(float t,const vec3 &v){
	return vec3(t*v.e[0],t*v.e[1],t*v.e[2]);
}

inline vec3 operator/(vec3 v, float t){
	return vec3(v.e[0]/t, v.e[1]/t, v.e[2]/t);
}

inline vec3 operator*(const vec3 &v, float t){
	return vec3(t*v.e[0],t*v.e[1],t*v.e[2]);
}

inline float distance(const vec3 &v1, const vec3 &v2){
	float t1 = v1[0] - v2[0];
	float t2 = v1[1] - v2[1];
	float t3 = v1[2] - v2[2];
	return sqrt(t1*t1+t2*t2+t3*t3);
}


inline float dot(const vec3 &v1, const vec3 &v2){
	return v1.e[0]*v2.e[0] + v1.e[1]*v2.e[1] + v1.e[2]*v2.e[2];
}

inline vec3 cross(const vec3 &v1,const vec3 &v2){
	return vec3( (v1.e[1]*v2.e[2] - v1.e[2]*v2.e[1]),
				-(v1.e[0]*v2.e[2] - v1.e[2]*v2.e[0]),
				(v1.e[0]*v2.e[1] - v1.e[1]*v2.e[0]) );
}

inline vec3& vec3::operator+=(const vec3 &v){
	e[0] += v.e[0];
	e[1] += v.e[1];
	e[2] += v.e[2];
	return *this;
}

inline vec3& vec3::operator*=(const vec3 &v){
	e[0] *= v.e[0];
	e[1] *= v.e[1];
	e[2] *= v.e[2];
	return *this;
}

inline vec3& vec3::operator/=(const vec3 &v){
	e[0] /= v.e[0];
	e[1] /= v.e[1];
	e[2] /= v.e[2];
	return *this;
}

inline vec3& vec3::operator-=(const vec3 &v){
	e[0] -= v.e[0];
	e[1] -= v.e[1];
	e[2] -= v.e[2];
	return *this;
}

inline vec3& vec3::operator*=(const float t){
	e[0] *= t;
	e[1] *= t;
	e[2] *= t;
	return *this;
}

inline vec3& vec3::operator/=(const float t){
	float k= 1.0/t;
	e[0] *= k;
	e[1] *= k;
	e[2] *= k;
	return *this;
}

inline void vec3::operator=(const vec3 v){
	set_x(v[0]);
	set_y(v[1]);
	set_z(v[2]);
	return;
}

inline vec3 unit_vector(vec3 v){
	return v / v.length();
}

//operations about vec4

inline void vec4::make_unit_vector(){
	float k=1.0 /sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2] + e[3]*e[3]);
	e[0]*=k; e[1]*=k;e[2]*=k; e[3]*=k;
}

inline vec4 operator+(const vec4 &v1, const vec4 &v2){
	return vec4(v1.e[0] + v2.e[0],v1.e[1]+v2.e[1], v1.e[2]+v2.e[2], v1.e[3]+v2.e[3]);
}

inline vec4 operator-(const vec4 &v1, const vec4 &v2){
	return vec4(v1.e[0]-v2.e[0],v1.e[1]-v2.e[1], v1.e[2]-v2.e[2], v1.e[3]-v2.e[3]);
}

inline vec4 operator*(const vec4 &v1, const vec4 &v2){
	return vec4(v1.e[0]*v2.e[0],v1.e[1]*v2.e[1], v1.e[2]*v2.e[2], v1.e[3]*v2.e[3]);
}

inline vec4 operator/(const vec4 &v1, const vec4 &v2){
	return vec4(v1.e[0]/v2.e[0],v1.e[1]/v2.e[1], v1.e[2]/v2.e[2], v1.e[3]/v2.e[3]);
}

inline vec4 operator*(float t,const vec4 &v){
	return vec4(t*v.e[0],t*v.e[1],t*v.e[2],t*v.e[3]);
}

inline vec4 operator/(vec4 v, float t){
	return vec4(v.e[0]/t, v.e[1]/t, v.e[2]/t, v.e[3]/t);
}

inline vec4 operator*(const vec4 &v, float t){
	return vec4(t*v.e[0],t*v.e[1],t*v.e[2],t*v.e[3]);
}

inline float dot(const vec4 &v1, const vec4 &v2){
	return v1.e[0]*v2.e[0] + v1.e[1]*v2.e[1] + v1.e[2]*v2.e[2] + v1.e[3]*v2.e[3];
}

inline vec4& vec4::operator+=(const vec4 &v){
	e[0] += v.e[0];
	e[1] += v.e[1];
	e[2] += v.e[2];
	e[3] += v.e[3];
	return *this;
}

inline vec4& vec4::operator*=(const vec4 &v){
	e[0] *= v.e[0];
	e[1] *= v.e[1];
	e[2] *= v.e[2];
	e[3] *= v.e[3];
	return *this;
}

inline vec4& vec4::operator/=(const vec4 &v){
	e[0] /= v.e[0];
	e[1] /= v.e[1];
	e[2] /= v.e[2];
	e[3] /= v.e[3];
	return *this;
}

inline vec4& vec4::operator-=(const vec4 &v){
	e[0] -= v.e[0];
	e[1] -= v.e[1];
	e[2] -= v.e[2];
	e[3] -= v.e[3];
	return *this;
}

inline vec4& vec4::operator*=(const float t){
	e[0] *= t;
	e[1] *= t;
	e[2] *= t;
	e[3] *= t;
	return *this;
}

inline vec4& vec4::operator/=(const float t){
	float k= 1.0/t;
	e[0] *= k;
	e[1] *= k;
	e[2] *= k;
	e[3] *= k;
	return *this;
}

inline void vec4::operator=(const vec4 v){
	set_x(v[0]);
	set_y(v[1]);
	set_z(v[2]);
	set_h(v[3]);
	return;
}

inline void vec4::standardization(){
	e[0]/=e[3];
	e[1]/=e[3];
	e[2]/=e[3];
	e[3] = 1;
}

inline vec4 unit_vector(vec4 v){
	return v / v.length();
}

//operations about multiply
vec4 operator*(Matrix44 mat, const vec4 v){
	vec4 res;
	for(int i=0;i<4;i++){
		float temp = 0;
		for(int j=0;j<4;j++){
			temp += mat[i][j]*v[j];
		}
		res.set(i,temp);
	}
	return res;
}

Matrix44 operator*(Matrix44 mat1, Matrix44 mat2){
	Matrix44 res;
	for(int i=0;i<4;i++){
		for(int j=0;j<4;j++){
			float temp = 0;
			for(int k=0;k<4;k++){
				temp += mat1[i][k] * mat2[k][j];
			}
			res.set(i,j,temp);
		}
	}
	return res;
}