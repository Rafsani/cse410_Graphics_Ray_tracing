#include <bits/stdc++.h>
#include <windows.h>
#include <GL/glut.h>
#include<stdlib.h>
using namespace std;
extern int recursion_level;
class pixelColor
{
public:
    double R,G,B;
    void set_color(double r,double g, double b)
    {
        R=r;
        G=g;
        B=b;
    }
};

struct point
{
    double x,y,z;

    point operator+(const point& p)
    {
        point res;
        res.x = x + p.x;
        res.y = y + p.y;
        res.z = z + p.z;
        return res;
    }

    point operator-(const point& p)
    {
        point res;
        res.x = x - p.x;
        res.y = y - p.y;
        res.z = z - p.z;
        return res;
    }

    point operator*(const double& scl)
    {
        point res;
        res.x = x*scl;
        res.y = y*scl;
        res.z = z*scl;
        return res;
    }

    point operator*(const point& p)
    {
        point res;
        res.x = y * p.z - z * p.y;
        res.y = z * p.x - p.z * x;
        res.z = x * p.y - p.x * y;
        return res;
    }

    point operator=(const point& p)
    {
        x = p.x;
        y = p.y;
        z = p.z;
        return *this;

    }

};


extern point pos,l,r,u;

double Dot_mult(point p1,point p2)
{
    return p1.x*p2.x + p1.y*p2.y + p1.z* p2.z;
}

point normalize_point (point arg)
{
    double val = sqrt(arg.x*arg.x + arg.y*arg.y + arg.z*arg.z);

    point p;
    p.x = arg.x/val;
    p.y = arg.y/val;
    p.z = arg.z/val;

    return p;
}



class Ray
{

public:
    point start;
    point dir_vect;

    Ray(point st, point dir)
    {
        start = st;
        dir_vect = normalize_point(dir);
    }
};

class LightSource
{
public:
    point p;
    LightSource() = default;
    double R,G,B;
    LightSource(double x, double y, double z)
    {
        p.x = x;
        p.y = y;
        p.z = z;
    }
    void set_color(double r, double g, double b)
    {
        R=r;
        G=g;
        B=b;
    }


    void draw()
    {
        glPushMatrix();
        glTranslated(p.x, p.y, p.z);
        drawSphere(2.0, 10, 10,R,G,B);
        glPopMatrix();
    }
};


class object
{
public:
    // should have x, y, z
    //double height, width, length ;
    double color[3] ;
    double coEfficients[4] ;// reflection coefficients
    int shine ;// exponent term of specular component
    object() {};
    virtual void draw() {} ;
    virtual double getTvalue(Ray &ray) {};
    virtual point normalvect(point intersect) { };
    void intersectmethod(Ray &ray,double t, pixelColor color,int level) {};
    void set_properties(double R,double G,double B,double ambient,double diffuse,double speculer,double reflection,double shines)
    {
        color[0]  =R;
        color[1]=G;
        color[2] = B;
        coEfficients[0] = ambient;
        coEfficients[1] = diffuse;
        coEfficients[2] = speculer;
        coEfficients[3] = reflection;
        shine = shines;
    }


    void intersect_method(Ray &ray,double t, pixelColor &p,int level);





};


std::vector <object*> objects;
std::vector<LightSource> lights;

void object::intersect_method(Ray &ray,double t, pixelColor &p,int level)
{
    point intersecting_point = ray.start +(ray.dir_vect * t)  ;
    point normal_atIntercectPoint = normalvect(intersecting_point);

    point reflectedRay = ray.dir_vect - normal_atIntercectPoint * (2.0 * Dot_mult(ray.dir_vect, normal_atIntercectPoint));
    reflectedRay = normalize_point(reflectedRay);

    p.set_color(color[0] *coEfficients[0],color[1]*coEfficients[0],color[2]*coEfficients[0]);

    for(int i=0; i<lights.size(); i++)
    {
        point lightToInsectDir = lights[i].p - intersecting_point;
        lightToInsectDir = normalize_point(lightToInsectDir);
        double lambert = 0.0,phong = 0.0;
        Ray L(intersecting_point + (lightToInsectDir * 0.01), lightToInsectDir);
        // point Vvect = normalize_point(intersecting_point * (-1.0));
        point Vvect = normalize_point( lightToInsectDir - intersecting_point);
       //point Vvect = normalize_point( intersecting_point - lightToInsectDir);
        point Rvect = normalize_point(L.dir_vect - normal_atIntercectPoint *(Dot_mult(L.dir_vect,normal_atIntercectPoint)*2.0) );
        bool obscured = false;

        for(int j=0; j<objects.size(); j++)
        {
            double isobscured = objects[j]->getTvalue(L);
            if(isobscured > 0)
            {
                obscured = true;
                break;
            }
        }


        if(!obscured)
        {
            lambert = coEfficients[1] * Dot_mult(L.dir_vect,normal_atIntercectPoint);
            phong = coEfficients[2] * pow(Dot_mult(Rvect,Vvect),shine);
            lambert = max(0.0,lambert);
            phong = max(0.0,phong);
        }

        p.R +=  ( phong*lights[i].R + lambert*lights[i].R) * color[0];
        p.G +=  (  phong * lights[i].G+ lambert * lights[i].G) *color[1];
        p.B += (  phong * lights[i].B + lambert * lights[i].B) * color[2];

        //cout << p.R << ' ' << p.G << ' ' << p.B << endl;

    }
    p.set_color(max(0.0,min(p.R,1.0)), max(0.0,min(p.G,1.0)), max(0.0,min(p.B,1.0)));

}


class Sphere: public object
{
public:
    point reference_point;
    double radius;
    Sphere()
    {
        reference_point.x = 0;
        reference_point.y = 0;
        reference_point.z = 0;
        radius = 0;
    }
    Sphere(double x,double y, double z, double rad)
    {
        reference_point.x = x;
        reference_point.y = y;
        reference_point.z = z;
        radius = rad;
    }

    point normalvect(point intersect)
    {
        point normal = intersect - reference_point;
        return normalize_point(normal);
        //return normal;
    }

    void draw()
    {
        glPushMatrix();
        {

            glTranslatef(reference_point.x,reference_point.y,reference_point.z);
            drawSphere(radius,50,50,color[0],color[1],color[2]);

        }
        glPopMatrix();


    }


    double getTvalue(Ray &ray)
    {
        point Rd = ray.dir_vect;
        point Ro = ray.start - reference_point;

        double discriminant = 4* Dot_mult(Rd,Ro) * Dot_mult(Rd,Ro) - 4* Dot_mult(Ro,Ro) + radius*radius;
        if(discriminant < 0 )
            return -1;
        double disc = sqrt(discriminant);
        double t1 = (- 2 * Dot_mult(Rd,Ro) + disc)/2;
        double t2 = (- 2 * Dot_mult(Rd,Ro) - disc)/2;
        //std::cout << "t = " << t1 << " " << t2 << std::endl;
        if(t1>t2)
            return t2;
        else
            return t1;
    }
};


class Triangle: public object
{
public:
    point points[3];
    double radius;
    Triangle(double x1,double y1, double z1,double x2,double y2, double z2,double x3,double y3, double z3)
    {
        points[0].x = x1;
        points[0].y = y1;
        points[0].z = z1;
        points[1].x = x2;
        points[1].y = y2;
        points[1].z = z2;
        points[2].x = x3;
        points[2].y = y3;
        points[2].z = z3;
    }

    void draw()
    {
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_TRIANGLES);
        {
            glVertex3f(points[0].x, points[0].y, points[0].z);
            glVertex3f(points[1].x, points[1].y, points[1].z);
            glVertex3f(points[2].x, points[2].y, points[2].z);
        }
        glEnd();
    }
    double getTvalue(Ray &ray)
    {
        const double EPSILON = 0.0000001;
        point edge1, edge2, h, s, q;
        double a,f,u,v;
        edge1 = points[1] - points[0];
        edge2 = points[2] - points[0];
        h = ray.dir_vect * edge2;
        a = Dot_mult(edge1,h);
        if (a > -EPSILON && a < EPSILON)
            return -1;    // This ray is parallel to this triangle.
        f = 1.0/a;
        s = ray.start - points[0];
        u = f * Dot_mult(s,h);
        if (u < 0.0 || u > 1.0)
            return -1;
        q = s*edge1;
        v = f * Dot_mult(ray.dir_vect,q);
        if (v < 0.0 || u + v > 1.0)
            return -1;
        // At this stage we can compute t to find out where the intersection point is on the line.
        float t = f * Dot_mult(edge2,q);
        if (t > EPSILON) // ray intersection
        {
            //outIntersectionPoint = rayOrigin + rayVector * t;
            return t;
        }
        else // This means that there is a line intersection but not a ray intersection.
            return -1;
    }

     point normalvect(point intersect)
    {
        point AB = points[1] - points[0];
        point AC = points[2] - points[0];
        point V = AB * AC;
        point P;
        return normalize_point(V);

    }
};


class General: public object
{
public:
    double  A, B, C, D, E, F, G, H, I, J;
    point   cube_reference_point;
    double length, width, height;
    General(double a,double b, double c, double d, double e, double f, double g, double h, double i, double j, double l, double w, double hh)
    {
        length = l ;
        width = w;
        height = hh;
        A=a;
        B=b;
        C=c;
        D=d;
        E=e;
        F=f;
        G=g;
        H=h;
        I=i;
        J=j;
    }
    void set_cube_reference_point(double a,double b, double c)
    {
        cube_reference_point.x = a;
        cube_reference_point.y = b;
        cube_reference_point.z = c;
    }

    void draw()
    {

    }


    bool cutoff(point intersect1)
    {
        if(intersect1.x <  cube_reference_point.x && length > 0 )
        {
            return true;
        }
        else if(intersect1.x >  (cube_reference_point.x + length)  && length > 0)
        {
            return true;
        }
        else if(intersect1.y <  cube_reference_point.y && width > 0)
        {
            return true;
        }
        else if(intersect1.y >  (cube_reference_point.y + width)  && width > 0)
        {
            return true;
        }
        else if(intersect1.z <  cube_reference_point.z  && height > 0)
        {
            return true;
        }
        else if(intersect1.z >  (cube_reference_point.z + height) && height > 0 )
        {
            return true;
        }

        return false;

    }

    double getTvalue(Ray &ray)
    {
        /// F(x, y, z) = Ax2 + By2 + Cz2 + Dxy+ Exz + Fyz + Gx + Hy + Iz + J = 0
        /// Aqt2 + Bqt + Cq = 0
        /// Aq = Axd2 + Byd2 + Czd2 + Dxdyd + Exdzd + Fydzd

        double a = (A*ray.dir_vect.x*ray.dir_vect.x) + (B * ray.dir_vect.y * ray.dir_vect.y) + (C * ray.dir_vect.z * ray.dir_vect.z ) + (D*ray.dir_vect.x*ray.dir_vect.y)
                    + (E * ray.dir_vect.x * ray.dir_vect.z ) + (F*ray.dir_vect.y * ray.dir_vect.z);

        double b = (2*A* ray.start.x * ray.dir_vect.x) + (2*B* ray.start.y * ray.dir_vect.y)
                    + (2*C* ray.start.z * ray.dir_vect.z) + D * (ray.start.x * ray.dir_vect.y + ray.start.y * ray.dir_vect.x)
                    + E * ( ray.start.x * ray.dir_vect.z + ray.start.z * ray.dir_vect.x )
                    + F * ( ray.start.y * ray.dir_vect.z + ray.start.z * ray.dir_vect.y )
                    + G * ray.dir_vect.x + H*ray.dir_vect.y + I* ray.dir_vect.z ;

        double c = (A* ray.start.x * ray.start.x) + (B* ray.start.y *  ray.start.y ) + (C *  ray.start.z *  ray.start.z )
                    + D* ray.start.x * ray.start.y + E*  ray.start.x * ray.start.z + F *  ray.start.y *  ray.start.z
                    + G*  ray.start.x + H* ray.start.y  + I *  ray.start.z  + J;


        double discriminent = b*b - 4*a*c;



        if( discriminent < 0) return -1;

        double t0 = (-b - sqrt(discriminent))*1.0 / (2*a);
        double t1 = (-b + sqrt(discriminent))*1.0 / (2*a);
        point intersect1 = ray.start + ray.dir_vect * t0;
        point intersect2 = ray.start + ray.dir_vect * t1;
        bool in_1  =  cutoff(intersect1);
        bool in_2 = cutoff(intersect2);


        if(in_1 && in_2) return -1;
        else if(in_1) return t1;
        else if(in_2) return t0;
        else return min(t0,t1);




    }

    point normalvect(point intersect)
    {
        point p;
        p.x = 2*A*intersect.x + D*intersect.y + E*intersect.z + G;
        p.y = 2*B*intersect.y + D*intersect.x + F*intersect.z + H;
        p.z = 2*C*intersect.z + E*intersect.x + F*intersect.y + I;
        return normalize_point(p);
      }
};



bool is_outside_range(double value, double min_val, double max_val)
{
    return (value < min_val) || (value > max_val);
}

class Floor: public object
{
public:
    point reference_point;
    double floorWidth,tileWidth;
    int no_of_tiles;
    Floor(double floor_Width, double tilesize)
    {
        floorWidth = floor_Width;
        reference_point.x = -floor_Width/2;
        reference_point.y = -floor_Width/2;
        reference_point.z = 0;
        tileWidth = tilesize;
        no_of_tiles = floor_Width/tilesize;
        set_properties(0,0,0,0.4, 0.2, 0.1, 0.3,5);
    }

    void draw()
    {
        //std::cout << "floor";
        glBegin(GL_QUADS);
        {
            for (int i = 0; i < no_of_tiles; i++)
                for (int j = 0; j < no_of_tiles; j++)
                {
                    bool c = (i + j) % 2;
                    glColor3f(c, c, c);

                    glVertex3f(reference_point.x + tileWidth * i, reference_point.y + tileWidth * j, reference_point.z);
                    glVertex3f(reference_point.x + tileWidth * (i + 1), reference_point.y + tileWidth * j, reference_point.z);
                    glVertex3f(reference_point.x + tileWidth * (i + 1), reference_point.y + tileWidth * (j + 1), reference_point.z);
                    glVertex3f(reference_point.x + tileWidth * i, reference_point.y + tileWidth * (j + 1), reference_point.z);
                }
        }
        glEnd();
    }
    double getTvalue(Ray &ray)
    {


        point normal = normalvect(reference_point);
        double up = Dot_mult(normal,ray.start);
        double down = Dot_mult(normal,ray.dir_vect);
        double t = (-1.0) * up/down;

        point intersecting_point = ray.start + ray.dir_vect * t;
        if(intersecting_point.x < reference_point.x || intersecting_point.x > -reference_point.x || intersecting_point.y < reference_point.y || intersecting_point.y > -reference_point.y)
        {
            return -1;
        }

        int x = (intersecting_point.x - reference_point.x) / tileWidth;
        int y = (intersecting_point.y - reference_point.y) / tileWidth;

        if((x + y) % 2 == 0)
        {
            color[0] = color[1] = color[2] = 0;
        }
        else
        {
            color[0] = color[1] = color[2] = 1;
        }


        return t;
    }




    point normalvect(point intersect)
    {
        point p;
        p.x=0;
        p.y=0;
        p.z=1;
        return p;
    }
};




