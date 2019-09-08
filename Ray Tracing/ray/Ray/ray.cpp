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
		p=2.0*temp -vec3(1,1,1);  //p的范围为-1~1 
	}while(p.squared_length()>=1.0);
	return p;
}

class lambertian : public material {
public:
	lambertian(texture* a):albedo(a) {}
	virtual bool scatter(const ray& r_in,const hit_record& rec,vec3& attenuation,ray& scattered) const{
		vec3 target = rec.p+rec.normal+random_in_unit_sphere();
		scattered = ray(rec.p,target-rec.p,r_in.time());  //分散，从p点出发的法向量加上随机化的一个向量 
		attenuation = albedo->value(0,0,rec.p);  //衰减 
		return true;   //不吸收任何光线 
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
		vec3 reflected =reflect(unit_vector(r_in.direction()),rec.normal);  //记录反射光线的方向 
		scattered=ray(rec.p,reflected+fuzz*random_in_unit_sphere(),r_in.time());  //分散，过P点，方向为反射光线
		attenuation =albedo; //衰减 
		return (dot(scattered.direction(),rec.normal)>0);   //如果为负数，那么说明该光线被吸收，则不用处理；
		                                                   //如果为正，说明光线被反射 
	}
	vec3 albedo;
	float fuzz;    //模糊参数 
};

//define refract

bool refract(const vec3& v,const vec3& n,float ni_over_nt, vec3& refracted){
	vec3 uv = unit_vector(v);
	float dt = dot(uv,n);  //cos theta1 
	float discriminant = 1.0 -ni_over_nt*ni_over_nt*(1-dt*dt);  //cos theta2 的平方(几何意义)
	if(discriminant > 0){   //如果这个值是负数，表示内全反射，即没有发生折射 
		refracted = ni_over_nt*(uv - n*dt) - n*sqrt(discriminant);
		return true;
	}
	else{
		return false;
	}
}

//define schlick  
//使用schlick近似来计算折射和反射的比例 

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
		
		//首先处理光线的情况，光线来自形状内部或者外部 
		
		if(dot(r_in.direction(),rec.normal)>0){   //此值为正说明光线来自形状内部 
			outward_normal = -rec.normal;    //需要反转法向量 
			ni_over_nt = ref_idx;      //折射率就表示为这个物体的折射率(ref_idx)除以空气的折射率(约等于1) 
			cosine = ref_idx * dot(r_in.direction(),rec.normal) / r_in.direction().length();  
			                         //求出光线与法向量之间夹角的cos值，不用加负号  
		}
		else{   //如果光线来自形状外部 
			outward_normal = rec.normal;  //法向量保持不变 
			ni_over_nt = 1.0 / ref_idx;   //折射率需要是空气的折射率(1)除以物体的折射率  
			cosine = -dot(r_in.direction(),rec.normal) / r_in.direction().length();
			                //计算cos值时要注意负号 
		}
		
		//预处理结束之后便可以开始进行折射或者反射的操作
		 
		if(refract(r_in.direction(),outward_normal,ni_over_nt,refracted)){  //计算可不可以产生折射 
			reflect_prob = schlick(cosine,ref_idx);  //如果可以的话就先计算出反射的概率 
		}
		else{
			scattered=ray(rec.p,reflected);   //如果不能产生折射，就直接进行反射 
			reflect_prob=1.0;   //反射的概率设为1 
		}
		
		//利用随机数来模拟真实情况 
		if(rand()/(double)RAND_MAX < reflect_prob ){
			scattered = ray(rec.p,reflected);     //如果小于这个概率就折射 
		}
		else{
			scattered = ray(rec.p,refracted);   //如果大于就反射 
		}
		return true;
	}
		
	float ref_idx; 
};

//加入了sphere类之后此函数不再需要 
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
	if(world->hit(r,0.001,MAXFLOAT,rec)){  //调用完hitable函数，使用0.001是要忽略掉和物体擦肩而过的光线 
		ray scattered;
		vec3 attenuation;
		vec3 emitted = rec.mat_ptr->emitted(rec.u,rec.v,rec.p);
		
		if(depth <50 && rec.mat_ptr->scatter(r,rec,attenuation,scattered)){
			return emitted + attenuation*color(scattered,world,depth+1);    //递归衰减光线 
		}
		else{
			return emitted;
		}
	}
	//如果光线不经过任何东西，反射天空的颜色即可 
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
				float u=float(i+rand()/double(RAND_MAX)) /float(nx);  //在200上的比例 
				float v=float(j+rand()/double(RAND_MAX)) /float(ny);  //在100上的比例，或者说是位置 
		                                                           //循环确定每一个像素点的rgb值 
		                                                           //随机值产生范围为0~1 
				ray r=cam.get_ray(u,v);  //获取当前点到照相机的向量，即光线 
				vec3 p = r.point_at_parameter(2.0);  
				col += color(r,world,0);  //当前光线到世界之后产生的颜色 
			}
			col/=float(ns);  //多次循环最后取平均值 
			col = vec3(sqrt(col[0]),sqrt(col[1]),sqrt(col[2]));  //使颜色变淡 
			int ir=int(255.99*col[0]);
			int ig=int(255.99*col[1]);
			int ib=int(255.99*col[2]);   //扩充到0~255 
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
	
	
	hitable *list[5];  //list是hitable的指针的数组 list是hitable** 
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
				float u=float(i+rand()/double(RAND_MAX)) /float(nx);  //在200上的比例 
				float v=float(j+rand()/double(RAND_MAX)) /float(ny);  //在100上的比例，或者说是位置 
		                                                           //循环确定每一个像素点的rgb值 
		                                                           //随机值产生范围为0~1 
				ray r=cam.get_ray(u,v);  //获取当前点到照相机的向量，即光线 
				vec3 p = r.point_at_parameter(2.0);  
				col += color(r,world,0);  //当前光线到世界之后产生的颜色 
			}
			col/=float(ns);  //多次循环最后取平均值 
			col = vec3(sqrt(col[0]),sqrt(col[1]),sqrt(col[2]));  //使颜色变淡 
			int ir=int(255.99*col[0]);
			int ig=int(255.99*col[1]);
			int ib=int(255.99*col[2]);   //扩充到0~255 
			fout<<ir<<" "<<ig<<" "<<ib<<"\n";
		}
	}
	return 0;
}

*/
