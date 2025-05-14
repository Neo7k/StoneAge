#pragma once

#include <math.h>

#include "Common.h"

struct i2
{
	int x;
	int y;

	constexpr i2 operator * (int i) const
	{
		return i2{ x * i, y * i };
	}

	constexpr i2 operator + (i2 other) const
	{
		return i2{ x + other.x, y + other.y };
	}
};

struct v4
{
	float x;
	float y;
	float z;
	float w;
};

struct b4
{
	uchar x;
	uchar y;
	uchar z;
	uchar w;
};

