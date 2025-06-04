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
	Quad() = default;
	Quad(v4 vert1, v4 vert2, v4 vert3, v4 vert4)
		: verts{vert1, vert2, vert3, vert4} {}
	Quad(v4 vert1, v4 vert2, v4 vert3, v4 vert4, v4 in_color)
		: verts{vert1, vert2, vert3, vert4}, color(in_color) {}

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
								m * q.verts[3],
								q.color};
	}

	v4 verts[4];
	v4 color{1.0f, 0.0f};
};

int main()
{
	Image<{128, 128}> frame;
	frame.Clear();

	GameWindow window(frame, 8);

	std::vector<Quad> quads;

	

	auto&& frame_fn =
		[&]
		()
		{
			frame.Clear();
			static float t = 0.0f;
			t += 0.005;
			quads.clear();
			for (int i = 0; i < 4; ++i)
			for (int j = 0; j < 2; ++j)
			{
				Mtx mtx
				{
					v4{sinf(t), 0.0f, cosf(t), -0.5f + j},
					v4{0.0f, 1.0f, 0.0f, 0.0f},
					v4{-cosf(t), 0.0f, sinf(t), 4.0f - (float)i},
				};
				quads.emplace_back(mtx*Quad
				{ 
					v4{-0.15f, 0.5f, 0.15f},
					v4{-0.15f, -0.5f, 0.15f},
					v4{-0.15f, -0.5f, -0.15f},
					v4{-0.15f, 0.5f, -0.15f},
					v4{0.0f, 0.0f, 1.0f}
				});
				quads.emplace_back(mtx*Quad
				{ 
					v4{0.15f, 0.5f, -0.15f},
					v4{0.15f, -0.5f, -0.15f},
					v4{0.15f, -0.5f, 0.15f},
					v4{0.15f, 0.5f, 0.15f},
					v4{0.0f, 0.0f, 1.0f}
				});
				quads.emplace_back(mtx*Quad
				{ 
					v4{-0.15f, 0.5f, -0.15f},
					v4{-0.15f, -0.5f, -0.15f},
					v4{0.15f, -0.5f, -0.15f},
					v4{0.15f, 0.5f, -0.15f},
					v4{1.0f, 0.0f, 0.0f}
				});
				quads.emplace_back(mtx*Quad
				{ 
					v4{0.15f, 0.5f, 0.15f},
					v4{0.15f, -0.5f, 0.15f},
					v4{-0.15f, -0.5f, 0.15f},
					v4{-0.15f, 0.5f, 0.15f},
					v4{1.0f, 0.0f, 0.0f}
				});
			}
			Mtx view
			{
				{1.0f, 0.0f, 0.0f, 0.0f},
				{0.0f, 1.0f, 0.0f, 0.0f},
				{0.0f, 0.0f, 1.0f, -t},
				{0.0f, 0.0f, 0.0f, 0.0f}
			};
			Mtx projection 
			{
				{1.0f, 0.0f, 0.0f, 0.0f},
				{0.0f, -1.0f, 0.0f, 0.0f},
				{0.0f, 0.0f, 1.0f, 0.0f},
				{0.0f, 0.0f, 1.0f, 0.0f} 
			};
			Mtx view_proj = projection * view; 
			for (Quad q : quads)
			{
				q = view_proj * q;
				if (q.verts[0].z <= 0.0f &&
						q.verts[1].z <= 0.0f &&
						q.verts[2].z <= 0.0f &&
						q.verts[3].z <= 0.0f)
				{
					continue;
				}
				for (v4& vert : q.verts)
					if (vert.w != 0.0f)
						vert = vert * (1.0f / vert.w); 

				Mtx to_screen_ctr
				{
					v4{0.5f, 0.0f, 0.0f, 0.5f},
					v4{0.0f, 0.5f, 0.0f, 0.5f}
				};
				q = to_screen_ctr * q;

				i2 top_left = Floor_i2(q.GetMin() * frame.GetSize());
				i2 bottom_right = Ceil_i2(q.GetMax() * frame.GetSize());
				top_left = Clamp(top_left, {0, 0}, frame.GetSize());
				bottom_right = Clamp(bottom_right, {0, 0}, frame.GetSize());
				v4 frame_size_inv{1.0f / frame.GetSize().x, 1.0f / frame.GetSize().y};
				for (int y = top_left.y; y < bottom_right.y; ++y)
					for (int x = top_left.x; x < bottom_right.x; ++x)
					{
						const float ms[][2] = {{0.05f, 0.25f}, {0.45f, 0.05f}, {0.55f, 0.95f}, {0.95f, 0.75f}};
						v4 pixel_color = v4::Zero();
						for (int i = 0; i < 4; ++i)
						{
							v4 p = v4{x + ms[i][0], y + ms[i][1]} * frame_size_inv; 
							float e[4] {q.GetEdgeValue(p, 0),
													q.GetEdgeValue(p, 1),
													q.GetEdgeValue(p, 2),
													q.GetEdgeValue(p, 3)};
							if (e[0] > 0 && e[1] > 0 && e[2] > 0 && e[3] > 0)
							{
								pixel_color = v4{q.color.x, q.color.y, q.color.z, pixel_color.w + 0.25f};
							}
							else if (e[0] < 0 && e[1] < 0 && e[2] < 0 && e[3] < 0)
							{
							}
						}
						b4 fp = frame.pix[y][x];
						v4 frame_pixel = v4{(float)fp.x, (float)fp.y, (float)fp.z, (float)fp.w} * (1.0f / 255.0f);
						frame_pixel = pixel_color * pixel_color.w + frame_pixel * (1.0f - pixel_color.w);
						frame.Set({x, y}, b4(frame_pixel * 255.0f));
					}
			}
		};
	
	window.Run(frame_fn);
	
	return 0;
}
