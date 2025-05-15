#include <iostream>
#include <thread>
#include <chrono>

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
	void Add(i2 coords, b4 color) { pix[coords.y][coords.x] += color; }

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
	frame.Clear();

	GameWindow window(frame, 32);

	auto&& frame_fn =
		[&]
		()
		{
			Quad q{ v4{0.45f, 0.25f}, v4{0.25f, 0.75f},
							v4{0.75f, 0.55f}, v4{0.95f, 0.05f}};

			static int x = 0;
			static int y = 0;
			static int pass = 0;
			++x;
			if (x >= frame.GetSize().x)
			{
				++y;
				x = 0;
			}

			if (y >= frame.GetSize().y)
			{
				x = 0;
				y = 0;
				++pass;
			}

			if (pass > 3)
				return;

				{
					v4 p = v4{x + 0.5f, y + 0.5f} * v4{1.0f / frame.GetSize().x, 1.0f / frame.GetSize().y};
					float signs[] =
					{
						(q.verts[1].x - q.verts[0].x) * (p.y - q.verts[0].y) - (q.verts[1].y - q.verts[0].y) * (p.x - q.verts[0].x),
						(q.verts[2].x - q.verts[1].x) * (p.y - q.verts[1].y) - (q.verts[2].y - q.verts[1].y) * (p.x - q.verts[1].x),
						(q.verts[3].x - q.verts[2].x) * (p.y - q.verts[2].y) - (q.verts[3].y - q.verts[2].y) * (p.x - q.verts[2].x),
						(q.verts[0].x - q.verts[3].x) * (p.y - q.verts[3].y) - (q.verts[0].y - q.verts[3].y) * (p.x - q.verts[3].x),
					};

					if (pass == 0)
						frame.Add({x, y}, {signs[0] < 0 ? (uchar)127u : (uchar)64u, 0, 0, 0});
					if (pass == 1)
						frame.Add({x, y}, {0, signs[1] < 0 ? (uchar)127u : (uchar)64u, 0, 0});
					if (pass == 2)
						frame.Add({x, y}, {0, 0, signs[2] < 0 ? (uchar)127u : (uchar)64u, 0});
					if (pass == 3)
						frame.Add({x, y}, signs[3] < 0 ? b4{127, 127, 127, 0} : b4{64, 64, 64, 0}); 
				}
		};
	
	window.Run(frame_fn);
	
	return 0;
}
