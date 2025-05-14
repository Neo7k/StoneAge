#include <iostream>

#include <assert.h>
#include <math.h>
#include <string.h>

#include "Common.h"
#include "Math.h"
#include "GameWindow.h"

template<i2 size>
struct Image
{
	constexpr i2 GetSize() const { return size; }
	void Clear() { memset(pix, 0, sizeof(pix)); }
	void Set(i2 coords, b4 color) { pix[coords.y][coords.x] = color; }

	b4 pix[size.y][size.x];
};

struct Quad
{
	v4 GetMin() const
	{
		v4 m = verts[0];
		for (int i = 1; i < 4; ++i)
		{
			m.x = verts[i].x < m.x ? verts[i].x : m.x;
			m.y = verts[i].y < m.y ? verts[i].y : m.y;
		}
		return m;
	}

	v4 GetMax() const
	{
		v4 m = verts[0];
		for (int i = 1; i < 4; ++i)
		{
			m.x = verts[i].x > m.x ? verts[i].x : m.x;
			m.y = verts[i].y > m.y ? verts[i].y : m.y;
		}
		return m;
	}

	v4 verts[4];
};

int main()
{
	Image<{32, 32}> frame;

	GameWindow window(frame, 32);

	auto&& frame_fn =
		[&]
		()
		{
			frame.Clear();

			Quad q{ v4{0.25f, 0.25f}, v4{0.75f, 0.25f},
							v4{0.25f, 0.75f}, v4{0.75f, 0.75f}};

			i2 top_left = Floor_i2(q.GetMin() * frame.GetSize());
			i2 bottom_right = Ceil_i2(q.GetMax() * frame.GetSize());
			for (int y = top_left.y; y < bottom_right.y; ++y)
				for (int x = top_left.x; x < bottom_right.x; ++x)
					frame.Set({x, y}, {255, 0, 0, 0});
		};
	
	window.Run(frame_fn);
	
	return 0;
}
