#pragma once
#include <string.h>
//#include <ctype.h>
//#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <math.h>
#include <float.h> //FLT_EPSILON, DBL_EPSILON
#define loopi(start_l,end_l) for ( int i=start_l;i<end_l;++i )

struct vector3
{
	double x, y, z;
};

struct vec3f
{
	double x, y, z;

	inline vec3f(void) {}

	//inline vec3f operator =( vector3 a )
	// { vec3f b ; b.x = a.x; b.y = a.y; b.z = a.z; return b;}

	inline vec3f(vector3 a)
	{
		x = a.x; y = a.y; z = a.z;
	}

	inline vec3f(const double X, const double Y, const double Z)
	{
		x = X; y = Y; z = Z;
	}

	inline vec3f operator + (const vec3f& a) const
	{
		return vec3f(x + a.x, y + a.y, z + a.z);
	}

	inline vec3f operator += (const vec3f& a) const
	{
		return vec3f(x + a.x, y + a.y, z + a.z);
	}

	inline vec3f operator * (const double a) const
	{
		return vec3f(x * a, y * a, z * a);
	}

	inline vec3f operator * (const vec3f a) const
	{
		return vec3f(x * a.x, y * a.y, z * a.z);
	}

	inline vec3f v3() const
	{
		return vec3f(x, y, z);
	}

	inline vec3f operator = (const vector3 a)
	{
		x = a.x;y = a.y;z = a.z;return *this;
	}

	inline vec3f operator = (const vec3f a)
	{
		x = a.x;y = a.y;z = a.z;return *this;
	}

	inline vec3f operator / (const vec3f a) const
	{
		return vec3f(x / a.x, y / a.y, z / a.z);
	}

	inline vec3f operator - (const vec3f& a) const
	{
		return vec3f(x - a.x, y - a.y, z - a.z);
	}

	inline vec3f operator / (const double a) const
	{
		return vec3f(x / a, y / a, z / a);
	}

	inline double dot(const vec3f& a) const
	{
		return a.x*x + a.y*y + a.z*z;
	}

	inline vec3f cross(const vec3f& a, const vec3f& b)
	{
		x = a.y * b.z - a.z * b.y;
		y = a.z * b.x - a.x * b.z;
		z = a.x * b.y - a.y * b.x;
		return *this;
	}

	inline double angle(const vec3f& v)
	{
		vec3f a = v, b = *this;
		double dot = v.x*x + v.y*y + v.z*z;
		double len = a.length() * b.length();
		if (len == 0)len = 0.00001f;
		double input = dot / len;
		if (input < -1) input = -1;
		if (input > 1) input = 1;
		return (double)acos(input);
	}

	inline double angle2(const vec3f& v, const vec3f& w)
	{
		vec3f a = v, b = *this;
		double dot = a.x*b.x + a.y*b.y + a.z*b.z;
		double len = a.length() * b.length();
		if (len == 0)len = 1;

		vec3f plane; plane.cross(b, w);

		if (plane.x * a.x + plane.y * a.y + plane.z * a.z > 0)
			return (double)-acos(dot / len);

		return (double)acos(dot / len);
	}

	inline vec3f rot_x(double a)
	{
		double yy = cos(a) * y + sin(a) * z;
		double zz = cos(a) * z - sin(a) * y;
		y = yy; z = zz;
		return *this;
	}
	inline vec3f rot_y(double a)
	{
		double xx = cos(-a) * x + sin(-a) * z;
		double zz = cos(-a) * z - sin(-a) * x;
		x = xx; z = zz;
		return *this;
	}
	inline void clamp(double min, double max)
	{
		if (x < min) x = min;
		if (y < min) y = min;
		if (z < min) z = min;
		if (x > max) x = max;
		if (y > max) y = max;
		if (z > max) z = max;
	}
	inline vec3f rot_z(double a)
	{
		double yy = cos(a) * y + sin(a) * x;
		double xx = cos(a) * x - sin(a) * y;
		y = yy; x = xx;
		return *this;
	}
	inline vec3f invert()
	{
		x = -x;y = -y;z = -z;return *this;
	}
	inline vec3f frac()
	{
		return vec3f(
			x - double(int(x)),
			y - double(int(y)),
			z - double(int(z))
		);
	}

	inline vec3f integer()
	{
		return vec3f(
			double(int(x)),
			double(int(y)),
			double(int(z))
		);
	}

	inline double length() const
	{
		return (double)sqrt(x*x + y * y + z * z);
	}

	inline vec3f normalize(double desired_length = 1)
	{
		double square = sqrt(x*x + y * y + z * z);
		/*
		if (square <= 0.00001f )
		{
		x=1;y=0;z=0;
		return *this;
		}*/
		//double len = desired_length / square;
		x /= square;y /= square;z /= square;

		return *this;
	}
};

class SymetricMatrix {
public:

	// Constructor

	SymetricMatrix(double c = 0) { loopi(0, 10) m[i] = c; }

	SymetricMatrix(double m11, double m12, double m13, double m14,
		double m22, double m23, double m24,
		double m33, double m34,
		double m44) {
		m[0] = m11;  m[1] = m12;  m[2] = m13;  m[3] = m14;
		m[4] = m22;  m[5] = m23;  m[6] = m24;
		m[7] = m33;  m[8] = m34;
		m[9] = m44;
	}

	// Make plane

	SymetricMatrix(double a, double b, double c, double d)
	{
		m[0] = a * a;  m[1] = a * b;  m[2] = a * c;  m[3] = a * d;
		m[4] = b * b;  m[5] = b * c;  m[6] = b * d;
		m[7] = c * c; m[8] = c * d;
		m[9] = d * d;
	}

	double operator[](int c) const { return m[c]; }

	// Determinant

	double det(int a11, int a12, int a13,
		int a21, int a22, int a23,
		int a31, int a32, int a33)
	{
		double det = m[a11] * m[a22] * m[a33] + m[a13] * m[a21] * m[a32] + m[a12] * m[a23] * m[a31]
			- m[a13] * m[a22] * m[a31] - m[a11] * m[a23] * m[a32] - m[a12] * m[a21] * m[a33];
		return det;
	}

	const SymetricMatrix operator+(const SymetricMatrix& n) const
	{
		return SymetricMatrix(m[0] + n[0], m[1] + n[1], m[2] + n[2], m[3] + n[3],
			m[4] + n[4], m[5] + n[5], m[6] + n[6],
			m[7] + n[7], m[8] + n[8],
			m[9] + n[9]);
	}

	SymetricMatrix& operator+=(const SymetricMatrix& n)
	{
		m[0] += n[0];   m[1] += n[1];   m[2] += n[2];   m[3] += n[3];
		m[4] += n[4];   m[5] += n[5];   m[6] += n[6];   m[7] += n[7];
		m[8] += n[8];   m[9] += n[9];
		return *this;
	}

	double m[10];
};
///////////////////////////////////////////

class Simplify
{
	// Global Variables & Strctures
public:
	struct Triangle { int v[3];double err[4];int deleted, dirty;vec3f n; };
	struct Vertex { vec3f p;int tstart, tcount;SymetricMatrix q;int border; };
	struct Ref { int tid, tvertex; };
	std::vector<Triangle> triangles;
	std::vector<Vertex> vertices;
	std::vector<Ref> refs;
	// Helper functions
public:

	double vertex_error(SymetricMatrix q, double x, double y, double z);
	double calculate_error(int id_v1, int id_v2, vec3f &p_result);
	bool flipped(vec3f p, int i0, int i1, Vertex &v0, Vertex &v1, std::vector<int> &deleted);
	void update_triangles(int i0, Vertex &v, std::vector<int> &deleted, int &deleted_triangles);
	void update_mesh(int iteration);
	void compact_mesh();
	double min(double v1, double v2);
	void simplify_mesh(int target_count, double agressiveness = 7, bool verbose = false);
	void simplify_mesh_lossless(bool verbose = false);
	void write_obj(const char* filename);
	void load_obj(const char* filename);
	void AddVertex(double x, double y, double z);
	void AddTriangle(int v0, int v1, int v2);
	int vCount();
	int triCount();
	void GetVertex(int index, vec3f &a);
	void GetTriangle(int index, std::vector<int> &v);
};
