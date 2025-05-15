#pragma once

#include <math.h>

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
};

i2 Floor_i2(v4 vec)
{
	return i2{(int)vec.x, (int)vec.y};
}

i2 Ceil_i2(v4 vec)
{
	return i2{(int)ceilf(vec.x), (int)ceilf(vec.y)};
}

struct b4
{
	uchar x;
	uchar y;
	uchar z;
	uchar w;

	void operator += (b4 other) { x += other.x; y += other.y; z += other.z; w += other.w; }
};

