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
	
	v4 operator * (i2 other) { return v4{x * other.x, y * other.y, z, w};	}
	v4 operator * (v4 other) { return v4{x * other.x, y * other.y, z * other.z, w * other.w}; }

	static v4 Zero()
	{
		return v4{0.0f, 0.0f, 0.0f, 0.0f};
	}

	float Dot(v4 v) const
	{
		return x * v.x + y * v.y + z * v.z + w * v.w;
	}

};

struct b4
{
	uchar x;
	uchar y;
	uchar z;
	uchar w;

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

	v4 operator * (v4 v) const
	{
		return v4{lines[0].Dot(v),
							lines[1].Dot(v),
							lines[2].Dot(v),
							lines[3].Dot(v)};
	}

	v4 lines[4];
};
