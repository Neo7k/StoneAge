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

	float GetEdgeValue(v4 p, int edge) const
	{
		int n = edge;
		int n_1 = n + 1 <= 3 ? n + 1 : 0;
		return (verts[n_1].y - verts[n].y) * (p.x - verts[n].x) - (verts[n_1].x - verts[n].x) * (p.y - verts[n].y);
	}

	friend Quad operator * (Mtx m, const Quad& q)
	{
		return Quad{m * q.verts[0],
								m * q.verts[1],
								m * q.verts[2],
								m * q.verts[3]};
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
			frame.Clear();
			Quad q{ v4{0.45f, 0.25f}, v4{0.25f, 0.75f},
							v4{0.75f, 0.55f}, v4{0.95f, 0.05f}};

			static float alpha = 0.0f;
			alpha += 0.02f;
			Mtx rotation{v4{cos(alpha), -sin(alpha), 0.0f, 0.0f},
									v4{sin(alpha), cos(alpha), 0.0f, 0.0f}};
			q = rotation * q;

			Mtx to_screen_ctr{v4{0.5f, 0.0f, 0.0f, 0.5f},
												v4{0.0f, 0.5f, 0.0f, 0.5f}};
			q = to_screen_ctr * q;

			i2 top_left = Floor_i2(q.GetMin() * frame.GetSize());
			i2 bottom_right = Ceil_i2(q.GetMax() * frame.GetSize());
			top_left = Clamp(top_left, {0, 0}, frame.GetSize());
			bottom_right = Clamp(bottom_right, {0, 0}, frame.GetSize());
			v4 frame_size_inv{1.0f / frame.GetSize().x, 1.0f / frame.GetSize().y};
			for (int y = top_left.y; y < bottom_right.y; ++y)
				for (int x = top_left.x; x < bottom_right.x; ++x)
				{
					v4 p = v4{x + 0.5f, y + 0.5f} * frame_size_inv; 
					if (q.GetEdgeValue(p, 0) > 0 &&
							q.GetEdgeValue(p, 1) > 0 &&
							q.GetEdgeValue(p, 2) > 0 &&
							q.GetEdgeValue(p, 3) > 0)
					{
						frame.Set({x, y}, {255, 0, 0, 0}); 
					}
				}
		};
	
	window.Run(frame_fn);
	
	return 0;
}
