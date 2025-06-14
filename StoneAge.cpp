#include <iostream>
#include <thread>
#include <chrono>

#include <assert.h>
#include <math.h>
#include <string.h>

#include "Common.h"
#include "Math.h"
#include "GameWindow.h"

constexpr int Key_W = 25;
constexpr int Key_A = 38;
constexpr int Key_S = 39;
constexpr int Key_D = 40;
constexpr int Key_Left = 113;
constexpr int Key_Right = 114;


template<typename T, i2 size>
struct Image
{
	constexpr i2 GetSize() const { return size; }
	void Clear() { memset(pix, 0, sizeof(pix)); }
	void Clear(T value) { for (int i = 0; i < size.y; ++i) for (int j = 0; j < size.x; ++j) pix[i][j] = value; }
	void Set(i2 coords, T color) { pix[coords.y][coords.x] = color; }
	void Add(i2 coords, T color) { pix[coords.y][coords.x] += color; }
	T Get(i2 coords) const { return pix[coords.y][coords.x]; }

	T pix[size.y][size.x];
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

	v4 GetNormal() const
	{
		return (verts[1] - verts[0]).Cross(verts[2] - verts[0]);
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

struct Input 
{
	Input() { memset(keys, 0, sizeof(keys)); memset(keys_time, 0, sizeof(keys_time)); }
	void Update(float dt) 
	{
		for (int i = 0; i < 256; ++i)
			if (keys[i])
				keys_time[i] += dt;
	}
	bool keys[256];
	float keys_time[256];
};

struct Camera
{
	v4 GetForward() const { return {-sinf(angle), 0.0f, cosf(angle), 0.0f}; }
	v4 GetRight() const { return {cosf(angle), 0.0f, sinf(angle), 0.0f}; }
	v4 pos = v4::Zero();
	float angle = 0.0f;
};

int main()
{
	Image<b4, {128, 128}> frame;
	Image<float, frame.GetSize()> depth;

	GameWindow window(frame, 8);

	Input input;
	auto&& key_fn = [&input](int key, bool is_pressed)
	{
		input.keys[key] = is_pressed;
		if (!is_pressed)
			input.keys_time[key] = 0.0f;
	};

	auto frame_start = std::chrono::steady_clock::now();

	Camera camera;
	std::vector<Quad> quads;

	auto&& frame_fn =
		[&]
		()
		{
			auto now = std::chrono::steady_clock::now();
			std::chrono::duration<float> dt_dur = now - frame_start;
			frame_start = now;
			float dt = dt_dur.count();

			input.Update(dt);
			static float t = 0.0f;
			t += dt;
			t = 0.0f;
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

			const float camera_speed = 1.0f;
			camera.pos += camera.GetForward() * ((float)input.keys[Key_W] * camera_speed * dt);
			camera.pos += camera.GetForward() * ((float)input.keys[Key_S] * camera_speed * dt * -1.0f);
			camera.pos += camera.GetRight() * ((float)input.keys[Key_D] * camera_speed * dt);
			camera.pos += camera.GetRight() * ((float)input.keys[Key_A] * camera_speed * dt * -1.0f);
			const float camera_rotation_speed = 1.75f;
			camera.angle += camera_rotation_speed * dt * (float)input.keys[Key_Left];
			camera.angle -= camera_rotation_speed * dt * (float)input.keys[Key_Right];

			Mtx view
			{
				{cosf(camera.angle), 0.0f, sinf(camera.angle), -camera.pos.x * cos(camera.angle) - camera.pos.z * sin(camera.angle)},
				{0.0f, 1.0f, 0.0f, 															-camera.pos.y																											},
				{-sinf(camera.angle), 0.0f, cosf(camera.angle), camera.pos.x * sin(camera.angle) - camera.pos.z * cos(camera.angle)},
				{0.0f, 0.0f, 0.0f, 1.0f}
			};
			Mtx projection 
			{
				{1.0f, 0.0f, 0.0f, 0.0f},
				{0.0f, -1.0f, 0.0f, 0.0f},
				{0.0f, 0.0f, 1.0f, 0.0f},
				{0.0f, 0.0f, 1.0f, 0.0f} 
			};
			Mtx view_proj = projection * view; 
			frame.Clear();
			const float far = 100.0f;
			depth.Clear(far);
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
					{
						vert = vert * v4{1.0f / vert.w, 1.0f / vert.w, 1.0f, 1.0f / vert.w};
					}
				
				Mtx to_screen_ctr
				{
					v4{0.5f, 0.0f, 0.0f, 0.5f},
					v4{0.0f, 0.5f, 0.0f, 0.5f},
					v4{0.0f, 0.0f, 1.0f, 0.0f}
				};
				q = to_screen_ctr * q;
				v4 normal = q.GetNormal();
				if (normal.z >= 0.0f)
					continue;
				
				i2 top_left = Floor_i2(q.GetMin() * frame.GetSize());
				i2 bottom_right = Ceil_i2(q.GetMax() * frame.GetSize());
				top_left = Clamp(top_left, {0, 0}, frame.GetSize());
				bottom_right = Clamp(bottom_right, {0, 0}, frame.GetSize());
				v4 frame_size_inv{1.0f / frame.GetSize().x, 1.0f / frame.GetSize().y};
				for (int y = top_left.y; y < bottom_right.y; ++y)
					for (int x = top_left.x; x < bottom_right.x; ++x)
					{
						float frame_depth = depth.Get({x, y});
						float pixel_depth = (q.verts[0].Dot(normal) - (v4{x + 0.5f, y + 0.5f} * frame_size_inv).Dot(normal)) / normal.z;
						if (pixel_depth > frame_depth)
							continue;

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
								depth.Set({x, y}, pixel_depth);
							}
						}
						
						b4 fp = frame.Get({x, y});
						v4 frame_pixel = v4{(float)fp.x, (float)fp.y, (float)fp.z, (float)fp.w} * (1.0f / 255.0f);
						frame_pixel = pixel_color * pixel_color.w + frame_pixel * (1.0f - pixel_color.w);
						frame.Set({x, y}, b4(frame_pixel * 255.0f));
					}
			}
		};
	
	window.Run(frame_fn, key_fn);
	
	return 0;
}
