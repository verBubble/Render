//This file includes the definition of scene, such as modelview, projection, view and so on.

#ifndef __SCENE_H__
#define __SCENE_H__

#include"geometry.h"

class Scene{
public:
	vec3 up;
	vec3 center;
	vec3 eye;
	vec3 light_dir;
	int width;
	int height;

	Matrix44 Modelview;
	Matrix44 Projection;
	Matrix44 View;

	Scene(vec3 Up, vec3 Center, vec3 Eye, vec3 light_Dir, int w, int h){
		light_dir = light_Dir;
		light_dir.make_unit_vector();
		
		up = Up;
		eye = Eye;
		center = Center;
		width = w;
		height = h;

		//modelview
		Matrix44 temp1;
		Matrix44 temp2;
		vec3 c = unit_vector(eye-center);
    	vec3 a = unit_vector(cross(up,c));
    	vec3 b = unit_vector(cross(c,a));
    	for (int i=0; i<3; i++) {
        	temp1.set(0,i,a[i]);
        	temp1.set(1,i,b[i]);
        	temp1.set(2,i,c[i]);
        	temp2.set(i,3,-center[i]);
   		}
   		Modelview = temp1*temp2;

   		//projection
   		Projection[3][2] = -1.0/(distance(eye,center));

   		//view
   		float x = (float)w/8;
   		float y = (float)h/8;
   		float w1 = (float)w*3/4;
   		float h1 = (float)h*3/4;
   		View.set(0,0, w1/2);
   		View.set(1,1, h1/2);
   		View.set(2,2, (float)255.0/2);
   		View.set(0,3, x+w1/2);
   		View.set(1,3, y+h1/2);
   		View.set(2,3, (float)255.0/2);
	}


	vec3 get_lightDir(){
		return light_dir;
	}
};

#endif