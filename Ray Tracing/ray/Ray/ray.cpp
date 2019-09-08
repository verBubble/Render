#include <iostream>
#include<fstream>
#include<limits>
#include<cstdlib>
#include<time.h>
#include "Ray.h"
using namespace std;

float MAXFLOAT = numeric_limits<float>::max();

vec3 random_in_unit_sphere(){
	vec3 p;
	do{
		vec3 temp=vec3(rand()/double(RAND_MAX),rand()/double(RAND_MAX),rand()/double(RAND_MAX));
		p=2.0*temp -vec3(1,1,1);  //p�ķ�ΧΪ-1~1 
	}while(p.squared_length()>=1.0);
	return p;
}

class lambertian : public material {
public:
	lambertian(texture* a):albedo(a) {}
	virtual bool scatter(const ray& r_in,const hit_record& rec,vec3& attenuation,ray& scattered) const{
		vec3 target = rec.p+rec.normal+random_in_unit_sphere();
		scattered = ray(rec.p,target-rec.p,r_in.time());  //��ɢ����p������ķ����������������һ������ 
		attenuation = albedo->value(0,0,rec.p);  //˥�� 
		return true;   //�������κι��� 
	}
	texture* albedo;
};

vec3 reflect(const vec3& v, const vec3& n){
	return v-2*dot(v,n)*n;
}

class metal:public material{
public:
	metal(const vec3& a,float f) : albedo(a) {if(f<1) fuzz=f; else fuzz =1;}
	virtual bool scatter(const ray& r_in,const hit_record& rec, vec3& attenuation, ray& scattered) const{
		vec3 reflected =reflect(unit_vector(r_in.direction()),rec.normal);  //��¼������ߵķ��� 
		scattered=ray(rec.p,reflected+fuzz*random_in_unit_sphere(),r_in.time());  //��ɢ����P�㣬����Ϊ�������
		attenuation =albedo; //˥�� 
		return (dot(scattered.direction(),rec.normal)>0);   //���Ϊ��������ô˵���ù��߱����գ����ô���
		                                                   //���Ϊ����˵�����߱����� 
	}
	vec3 albedo;
	float fuzz;    //ģ������ 
};

//define refract

bool refract(const vec3& v,const vec3& n,float ni_over_nt, vec3& refracted){
	vec3 uv = unit_vector(v);
	float dt = dot(uv,n);  //cos theta1 
	float discriminant = 1.0 -ni_over_nt*ni_over_nt*(1-dt*dt);  //cos theta2 ��ƽ��(��������)
	if(discriminant > 0){   //������ֵ�Ǹ�������ʾ��ȫ���䣬��û�з������� 
		refracted = ni_over_nt*(uv - n*dt) - n*sqrt(discriminant);
		return true;
	}
	else{
		return false;
	}
}

//define schlick  
//ʹ��schlick��������������ͷ���ı��� 

float schlick(float cosine, float ref_idx){
	float r0 = (1-ref_idx) / (1+ref_idx);
	r0 = r0*r0;
	return r0+(1-r0)*pow((1-cosine),5);
}

//define dielectric 

class dielectric : public material{
public:
	dielectric(float ri):ref_idx(ri) {}
	virtual bool scatter(const ray& r_in,const hit_record& rec,vec3& attenuation, ray& scattered) const {
		vec3 outward_normal;
		vec3 reflected = reflect(r_in.direction(),rec.normal);
		float ni_over_nt;
		attenuation = vec3(1.0,1.0,1.0);
		vec3 refracted;
		float reflect_prob;
		float cosine;
		
		//���ȴ�����ߵ����������������״�ڲ������ⲿ 
		
		if(dot(r_in.direction(),rec.normal)>0){   //��ֵΪ��˵������������״�ڲ� 
			outward_normal = -rec.normal;    //��Ҫ��ת������ 
			ni_over_nt = ref_idx;      //�����ʾͱ�ʾΪ��������������(ref_idx)���Կ�����������(Լ����1) 
			cosine = ref_idx * dot(r_in.direction(),rec.normal) / r_in.direction().length();  
			                         //��������뷨����֮��нǵ�cosֵ�����üӸ���  
		}
		else{   //�������������״�ⲿ 
			outward_normal = rec.normal;  //���������ֲ��� 
			ni_over_nt = 1.0 / ref_idx;   //��������Ҫ�ǿ�����������(1)���������������  
			cosine = -dot(r_in.direction(),rec.normal) / r_in.direction().length();
			                //����cosֵʱҪע�⸺�� 
		}
		
		//Ԥ�������֮�����Կ�ʼ����������߷���Ĳ���
		 
		if(refract(r_in.direction(),outward_normal,ni_over_nt,refracted)){  //����ɲ����Բ������� 
			reflect_prob = schlick(cosine,ref_idx);  //������ԵĻ����ȼ��������ĸ��� 
		}
		else{
			scattered=ray(rec.p,reflected);   //������ܲ������䣬��ֱ�ӽ��з��� 
			reflect_prob=1.0;   //����ĸ�����Ϊ1 
		}
		
		//�����������ģ����ʵ��� 
		if(rand()/(double)RAND_MAX < reflect_prob ){
			scattered = ray(rec.p,reflected);     //���С��������ʾ����� 
		}
		else{
			scattered = ray(rec.p,refracted);   //������ھͷ��� 
		}
		return true;
	}
		
	float ref_idx; 
};

//������sphere��֮��˺���������Ҫ 
float hit_sphere(const vec3& center, float radius,const ray& r){
	vec3 oc = r.origin() - center;
	float a = dot(r.direction(),r.direction());
	float b = 2.0*dot(oc,r.direction());
	float c = dot(oc,oc) - radius*radius;
	float discriminant =b*b - 4*a*c;
	if(discriminant <0){
		return -1.0;
	}
	else{
		return (-b-sqrt(discriminant))/(2.0*a);
	}
}

class flip_normals:public hitable{
	public:
		flip_normals(hitable *p):ptr(p) {}
		virtual bool hit(const ray& r,float t_min,float t_max,hit_record&  rec) const {
			if(ptr->hit(r,t_min,t_max,rec)){
				rec.normal = -rec.normal;
				return true;
			}
			else return false;
		}
		
		virtual bool bounding_box(float t0,float t1,aabb& box) const {
			return ptr->bounding_box(t0,t1,box);
		}
		hitable *ptr;
};


hitable *cornell_box(){
	hitable** list = new hitable*[6];
	int i=0;
	material *red = new lambertian(new constant_texture(vec3(0.65,0.05,0.05)));
	material *white = new lambertian(new constant_texture(vec3(0.73,0.73,0.73)));
	material *green = new lambertian(new constant_texture(vec3(0.12,0.45,0.15)));
	material *light = new diffuse_light(new constant_texture(vec3(15,15,15)));
	list[i++] = new flip_normals(new yz_rect(0,555,0,555,555,green));
	list[i++] = new yz_rect(0,555,0,555,0,red);
	list[i++] = new xz_rect(213,343,227,332,554,light);
	list[i++] = new flip_normals(new xz_rect(0,555,0,555,555,white));
	list[i++] = new xz_rect(0,555,0,555,0,white);
	list[i++] = new flip_normals(new xy_rect(0,555,0,555,555,white));
	return new hitable_list(list,i);
}

vec3 color(const ray& r,hitable *world,int depth) {
	hit_record rec;
	if(world->hit(r,0.001,MAXFLOAT,rec)){  //������hitable������ʹ��0.001��Ҫ���Ե��������������Ĺ��� 
		ray scattered;
		vec3 attenuation;
		vec3 emitted = rec.mat_ptr->emitted(rec.u,rec.v,rec.p);
		
		if(depth <50 && rec.mat_ptr->scatter(r,rec,attenuation,scattered)){
			return emitted + attenuation*color(scattered,world,depth+1);    //�ݹ�˥������ 
		}
		else{
			return emitted;
		}
	}
	//������߲������κζ�����������յ���ɫ���� 
	else{
		return vec3(0,0,0); 
	} 
}


hitable *random_scene(){
	texture *pertext = new noise_texture(4);
	hitable **list = new hitable*[4];
	list[0] = new sphere(vec3(0,-1000,0),1000,new lambertian(pertext));
	list[1] = new sphere(vec3(0,2,0),2,new lambertian(pertext));
	list[2] = new sphere(vec3(0,7,0),2,new diffuse_light(new constant_texture(vec3(4,4,4))));
	list[3] = new xy_rect(3,5,1,3,-2,new diffuse_light(new constant_texture(vec3(4,4,4))));
	return new hitable_list(list,4);
}


int main(){
	int nx=800;
	int ny=800;
	int ns=100;
	
	ofstream fout;
	fout.open("final.ppm");
	fout<<"P3\n"<<nx<<" "<<ny<<"\n255\n";
	float R = cos(M_PI/4);
	hitable *world = cornell_box();

	vec3 lookfrom(278,278,-800);
	vec3 lookat(278,278,0);
	float disk_to_focus = 10.0;
	float aperture = 0.0;
	camera cam(lookfrom,lookat,vec3(0,1,0),40,float(nx)/float(ny),aperture,disk_to_focus,0.0,1.0);

	for(int j=ny-1;j>=0;j--){
		for(int i=0;i<nx;i++){
			vec3 col(0,0,0);
			for(int s=0;s<ns;s++){
				float u=float(i+rand()/double(RAND_MAX)) /float(nx);  //��200�ϵı��� 
				float v=float(j+rand()/double(RAND_MAX)) /float(ny);  //��100�ϵı���������˵��λ�� 
		                                                           //ѭ��ȷ��ÿһ�����ص��rgbֵ 
		                                                           //���ֵ������ΧΪ0~1 
				ray r=cam.get_ray(u,v);  //��ȡ��ǰ�㵽������������������� 
				vec3 p = r.point_at_parameter(2.0);  
				col += color(r,world,0);  //��ǰ���ߵ�����֮���������ɫ 
			}
			col/=float(ns);  //���ѭ�����ȡƽ��ֵ 
			col = vec3(sqrt(col[0]),sqrt(col[1]),sqrt(col[2]));  //ʹ��ɫ�䵭 
			int ir=int(255.99*col[0]);
			int ig=int(255.99*col[1]);
			int ib=int(255.99*col[2]);   //���䵽0~255 
			fout<<ir<<" "<<ig<<" "<<ib<<"\n";
			cout<<i<<" "<<j<<endl;
		}
	}
	return 0;
}































/*int main()
{
	ofstream fout;
	int nx=200;
	int ny=100;
	int ns=100;
	//srand(time(0));
	fout.open("ray.ppm");
	fout<<"P3\n"<<nx<<" "<<ny<<"\n255\n";
	
	
	hitable *list[5];  //list��hitable��ָ������� list��hitable** 
	list[0] = new sphere(vec3(0,0,-1),0.5,new lambertian(vec3(0.1,0.2,0.5)));
	list[1] = new sphere(vec3(0,-100.5,-1),100,new lambertian(vec3(0.8,0.8,0.0)));
	list[2] = new sphere(vec3(1,0,-1),0.5,new metal(vec3(0.8,0.6,0.2),0.3));
	list[3] = new sphere(vec3(-1,0,-1),0.5,new dielectric(1.5));
	list[4] = new sphere(vec3(-1,0,-1),-0.45,new dielectric(1.5));
	hitable *world = new hitable_list(list,5);
	
	vec3 lookfrom(3,3,2);
	vec3 lookat(0,0,-1);
	float disk_to_focus = (lookfrom-lookat).length();
	float aperture = 2.0;
	camera cam(lookfrom,lookat,vec3(0,1,0),20,float(nx)/float(ny),aperture,disk_to_focus);
	
	
	for(int j=ny-1;j>=0;j--){
		for(int i=0;i<nx;i++){
			vec3 col(0,0,0);
			for(int s=0;s<ns;s++){
				float u=float(i+rand()/double(RAND_MAX)) /float(nx);  //��200�ϵı��� 
				float v=float(j+rand()/double(RAND_MAX)) /float(ny);  //��100�ϵı���������˵��λ�� 
		                                                           //ѭ��ȷ��ÿһ�����ص��rgbֵ 
		                                                           //���ֵ������ΧΪ0~1 
				ray r=cam.get_ray(u,v);  //��ȡ��ǰ�㵽������������������� 
				vec3 p = r.point_at_parameter(2.0);  
				col += color(r,world,0);  //��ǰ���ߵ�����֮���������ɫ 
			}
			col/=float(ns);  //���ѭ�����ȡƽ��ֵ 
			col = vec3(sqrt(col[0]),sqrt(col[1]),sqrt(col[2]));  //ʹ��ɫ�䵭 
			int ir=int(255.99*col[0]);
			int ig=int(255.99*col[1]);
			int ib=int(255.99*col[2]);   //���䵽0~255 
			fout<<ir<<" "<<ig<<" "<<ib<<"\n";
		}
	}
	return 0;
}

*/
