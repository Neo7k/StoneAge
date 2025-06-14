#pragma once

#include <math.h>

#include <algorithm>

#include "Common.h"

struct i2
{
	int x;
	int y;

	constexpr i2() = default;
	constexpr i2(int in_x, int in_y)
		: x(in_x), y(in_y) {}

	constexpr i2 operator * (int i) const { return i2{ x * i, y * i }; }
	constexpr i2 operator + (i2 other) const { return i2{ x + other.x, y + other.y }; }
};

struct v4
{
	float x;
	float y;
	float z;
	float w;

	v4() = default;
	v4(float in_x, float in_y, float in_z, float in_w)
		: x(in_x), y(in_y), z(in_z), w(in_w) {}
	v4(float in_x, float in_y)
		: x(in_x), y(in_y), z(0.0f), w(1.0f) {}
	v4(float in_x, float in_y, float in_z)
		: x(in_x), y(in_y), z(in_z), w(1.0f) {}
	
	v4 operator * (i2 other) const { return v4{x * other.x, y * other.y, z, w};	}
	v4 operator * (v4 other) const { return v4{x * other.x, y * other.y, z * other.z, w * other.w}; }
	v4 operator * (float other) const { return v4{x * other, y * other, z * other, w * other}; }
	v4 operator += (v4 other) { x += other.x; y += other.y; z += other.z; w += other.w; return *this; }
	v4 operator + (v4 other) const { return v4{x + other.x, y + other.y, z + other.z, w + other.w}; }
	v4 operator - (v4 other) const { return v4{x - other.x, y - other.y, z - other.z, w - other.w}; }

	static v4 Zero()
	{
		return v4{0.0f, 0.0f, 0.0f, 0.0f};
	}

	float Dot(v4 v) const
	{
		return x * v.x + y * v.y + z * v.z + w * v.w;
	}

	v4 Cross(v4 v) const
	{
		return {y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x, 0.0f};
	}

	v4 Normalized() const
	{
		return *this * (1.0f / sqrt(Dot(*this)));
	}

	friend std::ostream& operator << (std::ostream& os, v4 v)
	{
		os << "{" << v.x << ", " << v.y << ", " << v.z << ", " << v.w << "}";
		return os;
	}
};

struct b4
{
	uchar x;
	uchar y;
	uchar z;
	uchar w;

	b4() = default;
	b4(uchar value) : x(value), y(value), z(value), w(value) {}
	b4(uchar in_x, uchar in_y, uchar in_z, uchar in_w) : x(in_x), y(in_y), z(in_z), w(in_z) {}
	b4(v4 in) 
		: x((uchar)std::clamp(in.x, 0.0f, 255.0f))
		, y((uchar)std::clamp(in.y, 0.0f, 255.0f))
		, z((uchar)std::clamp(in.z, 0.0f, 255.0f))
		, w((uchar)std::clamp(in.w, 0.0f, 255.0f)) {}
	void operator += (b4 other) { x += other.x; y += other.y; z += other.z; w += other.w; }
};

i2 Floor_i2(v4 vec)
{
	return i2{(int)vec.x, (int)vec.y};
}

i2 Ceil_i2(v4 vec)
{
	return i2{(int)ceilf(vec.x), (int)ceilf(vec.y)};
}

i2 Clamp(i2 v, i2 min, i2 max)
{
	return {std::clamp(v.x, min.x, max.x), std::clamp(v.y, min.y, max.y)};
}

b4 Clamp(b4 v, uchar min, uchar max)
{
	return {std::clamp(v.x, min, max), std::clamp(v.y, min, max), std::clamp(v.z, min, max), std::clamp(v.w, min, max)};
}

struct Mtx
{
	Mtx() = default;
	Mtx(v4 l1, v4 l2)
	{
		lines[0] = l1;
		lines[1] = l2;
		lines[2] = v4::Zero();
		lines[3] = v4{0.0f, 0.0f};
	}
	Mtx(v4 l1, v4 l2, v4 l3)
	{
		lines[0] = l1;
		lines[1] = l2;
		lines[2] = l3;
		lines[3] = v4{0.0f, 0.0f};
	}
	Mtx(v4 l1, v4 l2, v4 l3, v4 l4)
	{
		lines[0] = l1;
		lines[1] = l2;
		lines[2] = l3;
		lines[3] = l4;
	}

	v4 operator * (v4 v) const
	{
		return v4{lines[0].Dot(v),
							lines[1].Dot(v),
							lines[2].Dot(v),
							lines[3].Dot(v)};
	}

	Mtx operator * (const Mtx& o) const
	{
		return
		{
			{
				lines[0].Dot({o.lines[0].x, o.lines[1].x, o.lines[2].x, o.lines[3].x}),
				lines[0].Dot({o.lines[0].y, o.lines[1].y, o.lines[2].y, o.lines[3].y}),
				lines[0].Dot({o.lines[0].z, o.lines[1].z, o.lines[2].z, o.lines[3].z}),
				lines[0].Dot({o.lines[0].w, o.lines[1].w, o.lines[2].w, o.lines[3].w})
			},
			{
				lines[1].Dot({o.lines[0].x, o.lines[1].x, o.lines[2].x, o.lines[3].x}),
				lines[1].Dot({o.lines[0].y, o.lines[1].y, o.lines[2].y, o.lines[3].y}),
				lines[1].Dot({o.lines[0].z, o.lines[1].z, o.lines[2].z, o.lines[3].z}),
				lines[1].Dot({o.lines[0].w, o.lines[1].w, o.lines[2].w, o.lines[3].w})
			},
			{
				lines[2].Dot({o.lines[0].x, o.lines[1].x, o.lines[2].x, o.lines[3].x}),
				lines[2].Dot({o.lines[0].y, o.lines[1].y, o.lines[2].y, o.lines[3].y}),
				lines[2].Dot({o.lines[0].z, o.lines[1].z, o.lines[2].z, o.lines[3].z}),
				lines[2].Dot({o.lines[0].w, o.lines[1].w, o.lines[2].w, o.lines[3].w})
			},
			{
				lines[3].Dot({o.lines[0].x, o.lines[1].x, o.lines[2].x, o.lines[3].x}),
				lines[3].Dot({o.lines[0].y, o.lines[1].y, o.lines[2].y, o.lines[3].y}),
				lines[3].Dot({o.lines[0].z, o.lines[1].z, o.lines[2].z, o.lines[3].z}),
				lines[3].Dot({o.lines[0].w, o.lines[1].w, o.lines[2].w, o.lines[3].w})
			}
		};
	}
	v4 lines[4];
};
