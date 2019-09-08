#include<math.h>
#include<stdlib.h>
#include<iostream>
using namespace std;

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

	inline float length() const{
		return sqrt(e[0]*e[0] + e[1]*e[1] + e[2]*e[2]);}
	inline float squared_length() const {
		return e[0]*e[0] + e[1]*e[1] + e[2]*e[2];}
	inline void make_unit_vector();

	float e[3];
};

inline istream& operator>>(istream &is,vec3 &t){
	is>>t.e[0]>>t.e[1]>>t.e[2];
	return is;
}

inline ostream& operator<<(ostream &os,const vec3 &t){
	os<<t.e[0]<<" "<<t.e[1]<<" "<<t.e[2];
	return os;
}

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

inline vec3 unit_vector(vec3 v){
	return v / v.length();
}

//define ray 

class ray
{
public:
	ray() {}
	ray(const vec3& a,const vec3& b,float ti = 0.0) {A=a; B=b;_time=ti;}  //add time
	vec3 origin() const {return A;}
	vec3 direction() const {return B;}
	float time() const {return _time;}  //new
	vec3 point_at_parameter(float t) const {return A+t*B;}

	vec3 A;
	vec3 B;
	float _time;  //new
};

//更新：添加了参数时间 

//define fitable
class material;
 
//add u,v for a plain
 
struct hit_record {
	float t;   //光线与hitable交点的t 
	float u;   //用于平面，指光线与平面交点的水平方向的占比 
	float v;   //用于平面，指光线与平面交点的竖直方向的占比 
	vec3 p;   //交点 
	vec3 normal;  //法线 
	material *mat_ptr;
};

//define aabb

inline float ffmin(float a,float b) {
	return a<b ? a:b;
} 

inline float ffmax(float a, float b){
	return a>b ? a:b;
}

class aabb{
	public:
		vec3 _min;
		vec3 _max;
		
		aabb(){}
		aabb(const vec3& a, const vec3& b){
			_min=a;
			_max=b;
		}
		
		vec3 min() const {return _min;}
		vec3 max() const {return _max;}
		
		bool hit(const ray& r,float tmin,float tmax) const {
			for(int a=0;a<3;a++){
				float t0 = ffmin((_min[a]-r.origin()[a])/r.direction()[a],
				(_max[a] - r.origin()[a])/r.direction()[a]);
				float t1 = ffmax((_min[a]-r.origin()[a])/r.direction()[a],
				(_max[a] - r.origin()[a])/r.direction()[a]);
				
				tmin = ffmin(t0,tmin);
				tmax = ffmax(t1,tmax);
				if(tmax <= tmin) return false;
			}
			return true;
		}
};

//define materials

class material {
public:
	virtual bool scatter(const ray& r_in, const hit_record& rec, vec3& attenuation, ray& scattered) const=0;
	virtual vec3 emitted(float u,float v,const vec3& p) const {return vec3(0,0,0);}
};

//define hitable
//add bounding_box

class hitable{
public:
	virtual bool hit(const ray& r,float t_min, float t_max, hit_record& rec) const = 0;
	virtual bool bounding_box(float t0,float t1, aabb& box) const = 0;
};

//define sphere

class sphere:public hitable{
public:
	sphere() {}
	sphere(vec3 cen, float r,material *m) :center(cen), radius(r),mat_ptr(m) {};
	virtual bool hit(const ray& r,float tmin,float tmax,hit_record& rec) const;
	virtual bool bounding_box(float t0,float t1,aabb& box) const;
	vec3 center;
	float radius;
	material *mat_ptr;
};

bool sphere::hit(const ray& r,float t_min,float t_max,hit_record& rec) const {
	vec3 oc = r.origin() -center;
	float a = dot(r.direction(),r.direction());
	float b = dot(oc,r.direction());  //注意这里没有乘2
	float c = dot(oc,oc) - radius*radius;
	float discriminant =b*b - a*c;  //所以这里不需要减去4ac
	if(discriminant > 0) {
		float temp = (-b-sqrt(b*b-a*c))/a;   //这里分子和分母都约去了2
		if(temp<t_max && temp>t_min){
			rec.t =temp;
			rec.p =r.point_at_parameter(rec.t);
			rec.normal=(rec.p -center)/radius; //圆心指向被撞击点的单位向量
			rec.mat_ptr=mat_ptr;
			return true;
		}
		temp = (-b+sqrt(b*b-a*c))/a;
		if(temp<t_max && temp>t_min){
			rec.t =temp;
			rec.p =r.point_at_parameter(rec.t);
			rec.normal=(rec.p -center)/radius; //圆心指向被撞击点的单位向量
			rec.mat_ptr=mat_ptr;
			return true;
		}
	}
	return false;
}

//add funtion bounding_box

bool sphere::bounding_box(float t0,float t1,aabb& box) const{
	box = aabb(center-vec3(radius,radius,radius),center + vec3(radius,radius,radius));
	return true;
}

//define hitable_list 

class hitable_list: public hitable{
public:
	hitable_list(){}
	hitable_list(hitable **l,int n) {list =l; list_size =n;}
	virtual bool hit(const ray& r,float tmin,float tmax,hit_record& rec) const;
	virtual bool bounding_box(float t0,float t1,aabb& box) const;
	hitable **list;
	int list_size;
};

//add bounding_box funtion for hitable_list

aabb surrounding_box(aabb box0,aabb box1);

bool hitable_list::bounding_box(float t0,float t1,aabb& box) const {
	if(list_size<1) return false;
	aabb temp_box;
	bool first_true = list[0]->bounding_box(t0,t1,temp_box);
	if(!first_true) return false;
	else{
		box = temp_box;
	}
	
	for(int i=1;i<list_size;i++){
		if(list[0]->bounding_box(t0,t1,temp_box)) box=surrounding_box(box,temp_box);
		else return false;
	}
	return true;
}

bool hitable_list::hit(const ray& r,float t_min, float t_max, hit_record& rec) const {
	hit_record temp_rec;
	bool hit_anything =false;
	double closest_so_far = t_max;
	for(int i=0;i<list_size;i++){
		if(list[i]->hit(r,t_min,closest_so_far,temp_rec)){
			hit_anything = true;
			closest_so_far = temp_rec.t;
			rec = temp_rec;
		}
	}
	return hit_anything;
}

//random_in_unit_disk
 
vec3 random_in_unit_disk(){
	vec3 p;
	do{
		vec3 temp = vec3(rand()/(double)RAND_MAX,rand()/(double)RAND_MAX,0);
		p=2.0*temp - vec3(1,1,0);
	}while(dot(p,p) >= 1.0);
	return p;
} 

//define camera 

class camera{
public:
	//new: add t0 and t1
	camera(vec3 lookfrom, vec3 lookat, vec3 vup,float vfov,float aspect,float aperture, float focus_disk,
	float t0, float t1){
		// vfov代表的是欧拉角，aperture是光圈，aspect是画布的长宽比 
		time0 =t0;
		time1 =t1;
		lens_radius = aperture/2;
		float theta = vfov*M_PI/180;
		float half_height = tan(theta/2);
		float half_width = aspect * half_height;
 
		origin=lookfrom;    //定义照相机的原点位置，即光线必须经过的地方 	
		w=unit_vector(lookfrom-lookat);
		u=unit_vector(cross(vup,w));
		v=cross(w,u);	   //建立正交坐标系 
		
		lower_left_corner = origin - focus_disk*half_width*u - focus_disk*half_height*v - focus_disk*w;
		horizontal=2*half_width*focus_disk*u;   //水平方向的长度 
		vertical=2*half_height*focus_disk*v;  //竖直方向的长度，以上三个是定义画布的位置和大小
	}

	// add time
	ray get_ray(float s,float t){
		vec3 rd = lens_radius*random_in_unit_disk();
		vec3 offset = u*rd.x() + v*rd.y();
		float time = time0 + (rand()/(double)RAND_MAX)*(time1-time0);
		return ray(origin + offset, lower_left_corner + s*horizontal + t*vertical - origin -offset,time); 
	}
	
	vec3 lower_left_corner;
	vec3 horizontal;
	vec3 vertical;
	vec3 origin;
	vec3 u,v,w;
	float lens_radius;
	float time0,time1;  //add time,time0 and time1 are open time and close time
};

// add moving sphere

class moving_sphere:public hitable{
public:
	moving_sphere() {}
	moving_sphere(vec3 cen0, vec3 cen1, float t0, float t1,float r,material *m) :
	center0(cen0),center1(cen1),time0(t0),time1(t1), radius(r),mat_ptr(m) {};  //add cen0 for t0 and cen1 for t1
	virtual bool hit(const ray& r,float tmin,float tmax,hit_record& rec) const;
	virtual bool bounding_box(float t0,float t1,aabb& box) const;
	vec3 center(float time) const;  //new
	vec3 center0,center1;  // add center
	float radius;
	float time0,time1;  //and time
	material *mat_ptr;
};

//add to know the center for a exact time 

vec3 moving_sphere::center(float time) const {
	return center0 + ((time-time0)/(time1-time0))*(center1-center0);
}

//replace center with center(r.time())

bool moving_sphere::hit(const ray& r,float t_min,float t_max,hit_record& rec) const {
	vec3 oc = r.origin() -center(r.time());
	float a = dot(r.direction(),r.direction());
	float b = dot(oc,r.direction());  //注意这里没有乘2
	float c = dot(oc,oc) - radius*radius;
	float discriminant =b*b - a*c;  //所以这里不需要减去4ac
	if(discriminant > 0) {
		float temp = (-b-sqrt(b*b-a*c))/a;   //这里分子和分母都约去了2
		if(temp<t_max && temp>t_min){
			rec.t =temp;
			rec.p =r.point_at_parameter(rec.t);
			rec.normal=(rec.p -center(r.time()))/radius; //圆心指向被撞击点的单位向量
			rec.mat_ptr=mat_ptr;
			return true;
		}
		temp = (-b+sqrt(b*b-a*c))/a;
		if(temp<t_max && temp>t_min){
			rec.t =temp;
			rec.p =r.point_at_parameter(rec.t);
			rec.normal=(rec.p - center(r.time()))/radius; //圆心指向被撞击点的单位向量
			rec.mat_ptr=mat_ptr;
			return true;
		}
	}
	return false;
}

//add bounding_box funtion for moving sphere

aabb surrounding_box(aabb box0,aabb box1){
	vec3 small( fmin(box0.min().x(),box1.min().x()),
				fmin(box0.min().y(),box1.min().y()),
				fmin(box0.min().z(),box1.min().z()) );
	vec3 big( fmax(box0.max().x(),box1.max().x()),
				fmax(box0.max().y(),box1.max().y()),
				fmax(box0.max().z(),box1.max().z()) );
	return aabb(small, big);				
}

bool moving_sphere::bounding_box(float t0,float t1,aabb& box) const {
	aabb box0(center(t0)-vec3(radius,radius,radius),center(t0)+vec3(radius,radius,radius));
	aabb box1(center(t1)-vec3(radius,radius,radius),center(t1)+vec3(radius,radius,radius));
	box = surrounding_box(box0,box1);
	return true;
}


//define bvh_node

class bvh_node : public hitable{
	public:
		bvh_node(){}
		bvh_node(hitable **l,int n,float time0,float time1);
		virtual bool hit(const ray& r,float tmin,float tmax,hit_record& rec) const;
		virtual bool bounding_box(float t0, float t1,aabb& box) const;
		hitable *left;
		hitable *right;
		aabb box;
};

bool bvh_node::bounding_box(float t0,float t1, aabb &b) const {
	b = box;
	return true;
}

int box_x_compare (const void*a,const void* b){
	aabb box_left,box_right;
	hitable *ah = *(hitable**) a;
	hitable *bh = *(hitable**) b;
	if(!ah->bounding_box(0,0,box_left) || !bh->bounding_box(0,0,box_right))
				std::cerr<<"no bounding box in bvh_node constructor\n";
	if(box_left.min().x() - box_right.min().x() < 0.0) 
				return -1;
	else
		return 1;
}

int box_y_compare (const void*a,const void* b){
	aabb box_left,box_right;
	hitable *ah = *(hitable**) a;
	hitable *bh = *(hitable**) b;
	if(!ah->bounding_box(0,0,box_left) || !bh->bounding_box(0,0,box_right))
				std::cerr<<"no bounding box in bvh_node constructor\n";
	if(box_left.min().y() - box_right.min().y() < 0.0) 
				return -1;
	else
		return 1;
}

int box_z_compare (const void*a,const void* b){
	aabb box_left,box_right;
	hitable *ah = *(hitable**) a;
	hitable *bh = *(hitable**) b;
	if(!ah->bounding_box(0,0,box_left) || !bh->bounding_box(0,0,box_right))
				std::cerr<<"no bounding box in bvh_node constructor\n";
	if(box_left.min().z() - box_right.min().z() < 0.0) 
				return -1;
	else
		return 1;
}

bool bvh_node::hit(const ray& r,float t_min,float t_max,hit_record& rec) const {
	if(box.hit(r,t_min,t_max)){
		hit_record left_rec,right_rec;
		bool hit_left=left->hit(r,t_min,t_max,left_rec);
		bool hit_right=right->hit(r,t_min,t_max,right_rec);
		if(hit_left && hit_right){
			if(left_rec.t<right_rec.t)
				rec = left_rec;
			else rec = right_rec;
			return true;
		}
		else if(hit_left) {
			rec = left_rec;
			return true;
		}
		else if(hit_right){
			rec = right_rec;
			return true;
		}
		else return false;
	}
	else return false;
}

bvh_node::bvh_node(hitable** l,int n,float time0,float time1){
	int axis = int(3*(rand()/RAND_MAX));
	if(axis == 0){
		qsort(l,n,sizeof(hitable*),box_x_compare);
	}
	else if (axis ==1) qsort(l,n,sizeof(hitable*),box_y_compare);
	else qsort(l,n,sizeof(hitable*),box_z_compare);
	if(n==1) left = right =l[0];
	else if(n==2) {
		left = l[0];
		right = l[1];
	}
	else {
		left = new bvh_node(l,n/2,time0,time1);
		right = new bvh_node(l+n/2,n-n/2,time0,time1);
	}
	aabb box_left,box_right;
	if(!left->bounding_box(time0,time1,box_left) || !right->bounding_box(time0,time1,box_right))
	  std::cerr<<"no bounding box in bvh_node consructor\n";
	box=surrounding_box(box_left,box_right);
}


//define texture

class texture{
	public:
		virtual vec3 value(float u,float v,const vec3& p) const =0;
};


class constant_texture:public texture{
	public:
		constant_texture() {}
		constant_texture(vec3 c) :color(c) {}
		virtual vec3 value(float u,float x,const vec3& p) const {
			return color;
		}
		vec3 color;
};

class checker_texture :public texture{
	public:
		checker_texture() {}
		checker_texture(texture *t0,texture *t1) : even(t0),odd(t1) {}
		virtual vec3 value(float u,float v,const vec3& p)const {
			float sines =sin(10*p.x())*sin(10*p.y())*sin(10*p.z());
			if(sines < 0 ) return odd->value(u,v,p);
			else return even->value(u,v,p);		
		}
		
		texture *odd;
		texture *even;
};


//define perlin noise
 
inline float perlin_interp(vec3 c[2][2][2], float u, float v, float w) {  //传过来的参数分别是三个维度上要求的点的坐标x与floor(x)的差 
    float uu = u*u*u*(u*(u*6-15)+10);    //c是梯度 
    float vv = v*v*v*(v*(v*6-15)+10);
    float ww = w*w*w*(w*(w*6-15)+10);  //平滑函数 
    float accum = 0;
    for (int i=0; i < 2; i++)
        for (int j=0; j < 2; j++)
            for (int k=0; k < 2; k++) {
                vec3 weight_v(u-i, v-j, w-k);   //循环遍历从每一个顶点到要求的点P的坐标 
                accum += (i*uu + (1-i)*(1-uu))*
                    (j*vv + (1-j)*(1-vv))*
                    (k*ww + (1-k)*(1-ww))*dot(c[i][j][k], weight_v);
            }
    return accum/8;
}

class perlin {
    public:
        float noise(const vec3& p) const {  //给某点坐标，记录该点的坐标值 
            float u = p.x() - floor(p.x());
            float v = p.y() - floor(p.y());
            float w = p.z() - floor(p.z());  //计算tx，ty，tz 
            int i = floor(p.x());
            int j = floor(p.y());
            int k = floor(p.z());
            vec3 c[2][2][2];
            for (int di=0; di < 2; di++)
                for (int dj=0; dj < 2; dj++)
                    for (int dk=0; dk < 2; dk++){
                    	int x = (i+di) & 255;
                    	int y = (j+dj) & 255;
                    	int z = (k+dk) & 255;
                    	int temp = perm_x[perm_x[perm_x[x]+y]+z];
                    	c[di][dj][dk] = ranvec[temp];
					}
          //c[di][dj][dk] = ranvec[perm_y[perm_x[(i+di) & 255] + ((j+dj) & 255)] +( (k+dk) & 255)]; //将一个格子中的8个顶点中的向量取出来 
            return perlin_interp(c, u, v, w);  //进行interpolation 
        }
        float turb(const vec3& p, int depth=7) const {
            float accum = 0;
            vec3 temp_p = p;
            float weight = 1.0;
            for (int i = 0; i < depth; i++) {
                accum += weight*noise(temp_p);
                weight *= 0.5;
                temp_p *= 2;
            }
            return fabs(accum);
        }
        static vec3 *ranvec;
        static int *perm_x;
        static int *perm_y;
        static int *perm_z;
};

static vec3* perlin_generate() {  
    vec3 * p = new vec3[256];
    for ( int i = 0; i < 256; ++i )
        p[i] = unit_vector(vec3(-1 + 2*(rand()/RAND_MAX), -1 + 2*(rand()/RAND_MAX), -1 + 2*(rand()/RAND_MAX)));
    return p;
}

//std::uniform_int_distribution distributionInt;
//aotu diceInt = std::bind(distributionInt,generator);

void permute(int *p, int n) {
    for (int i = n-1; i > 0; i--) {
        int target = int((rand()/RAND_MAX) * 255);
        int tmp = p[i];
        p[i] = p[target];
        p[target] = tmp;
    }
    for(int i=n;i<2*n;i++) p[i]=p[i-n];
    return;
}



static int* perlin_generate_perm() {  
    int * p = new int[512];
    for (int i = 0; i < 256; i++)
        p[i] = i;
    permute(p, 256);
    return p;
}

vec3 *perlin::ranvec = perlin_generate();   //随机生成最初的三维向量的数组 哈希表 
int *perlin::perm_x = perlin_generate_perm();  //先按照顺序生成数组，然后打乱，从前面那个函数生成的表中取向量 

class noise_texture: public texture{
	public:
		noise_texture() {}
		noise_texture(float sc):scale(sc) {}
		virtual vec3 value(float u, float v,const vec3& p) const {
			return vec3(1,1,1)*0.5*(1 + sin(scale*p.x() + 5*noise.turb(scale*p))) ;
		}
		perlin noise;	
		float scale;
};


//define light

class diffuse_light:public material{
	public:
		diffuse_light(texture *a):emit(a) {}
		virtual bool scatter(const ray& r_in,const hit_record& rec,vec3& attenuation,ray& scattered) const {
			return false;
		}
		virtual vec3 emitted(float u,float v,const vec3& p) const {
			return emit->value(u,v,p);
		}
		texture *emit;
}; 

//define x-y plain

class xy_rect:public hitable{
	public:
		xy_rect() {}
		xy_rect(float _x0,float _x1,float _y0,float _y1,float _k,material *mat):x0(_x0),x1(_x1),y0(_y0),y1(_y1),k(_k),
		mp(mat) {} ;
		virtual bool hit(const ray& r,float t0,float t1,hit_record& rec) const;
		virtual bool bounding_box(float t0,float t1,aabb& box) const {
			box = aabb(vec3(x0,y0,k-0.0001),vec3(x1,y1,k+0.0001));
			return true;
		}
		material *mp;
		float x0,x1,y0,y1,k;
};

bool xy_rect::hit(const ray& r,float t0,float t1,hit_record& rec) const {
	float t=(k-r.origin().z())/r.direction().z();
	if(t<t0 || t>t1) return false;
	float x=r.origin().x() + t*r.direction().x();
	float y=r.origin().y() + t*r.direction().y();
	if(x<x0 || x>x1 ||y<y0 || y>y1) return false;
	rec.u = (x-x0)/(x1-x0);
	rec.v = (y-y0)/(y1-y0);
	rec.t = t;
	rec.mat_ptr = mp;
	rec.p = r.point_at_parameter(t);
	rec.normal = vec3(0,0,1);
	return true;
}

class xz_rect: public hitable  {
    public:
        xz_rect() {}
        xz_rect(float _x0, float _x1, float _z0, float _z1, float _k, material *mat) : x0(_x0), x1(_x1), z0(_z0), z1(_z1), k(_k), mp(mat) {};
        virtual bool hit(const ray& r, float t0, float t1, hit_record& rec) const;
        virtual bool bounding_box(float t0, float t1, aabb& box) const {
               box =  aabb(vec3(x0,k-0.0001,z0), vec3(x1, k+0.0001, z1));
               return true; }
        material  *mp;
        float x0, x1, z0, z1, k;
};

class yz_rect: public hitable  {
    public:
        yz_rect() {}
        yz_rect(float _y0, float _y1, float _z0, float _z1, float _k, material *mat) : y0(_y0), y1(_y1), z0(_z0), z1(_z1), k(_k), mp(mat) {};
        virtual bool hit(const ray& r, float t0, float t1, hit_record& rec) const;
        virtual bool bounding_box(float t0, float t1, aabb& box) const {
               box =  aabb(vec3(k-0.0001, y0, z0), vec3(k+0.0001, y1, z1));
               return true; }
        material  *mp;
        float y0, y1, z0, z1, k;
};

bool xz_rect::hit(const ray& r, float t0, float t1, hit_record& rec) const {
    float t = (k-r.origin().y()) / r.direction().y();
    if (t < t0 || t > t1)
        return false;
    float x = r.origin().x() + t*r.direction().x();
    float z = r.origin().z() + t*r.direction().z();
    if (x < x0 || x > x1 || z < z0 || z > z1) 
        return false;
    rec.u = (x-x0)/(x1-x0);
    rec.v = (z-z0)/(z1-z0); 
    rec.t = t;
    rec.mat_ptr = mp;
    rec.p = r.point_at_parameter(t);
    rec.normal = vec3(0, 1, 0);
    return true;
}

bool yz_rect::hit(const ray& r, float t0, float t1, hit_record& rec) const {
    float t = (k-r.origin().x()) / r.direction().x();
    if (t < t0 || t > t1)
        return false;
    float y = r.origin().y() + t*r.direction().y();
    float z = r.origin().z() + t*r.direction().z();
    if (y < y0 || y > y1 || z < z0 || z > z1) 
        return false;
    rec.u = (y-y0)/(y1-y0);
    rec.v = (z-z0)/(z1-z0); 
    rec.t = t;
    rec.mat_ptr = mp;
    rec.p = r.point_at_parameter(t);
    rec.normal = vec3(1, 0, 0);
    return true;
}


