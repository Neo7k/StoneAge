#include <iostream>
#include <thread>
#include <chrono>

#include <assert.h>
#include <math.h>
#include <string.h>

#include "Common.h"
#include "Math.h"
#include "GameWindow.h"

#define CLOG(condition, ...) \
	if (condition) \
		std::cout << std::format(__VA_ARGS__) << std::endl;

constexpr int Key_W = 25;
constexpr int Key_A = 38;
constexpr int Key_S = 39;
constexpr int Key_D = 40;
constexpr int Key_Left = 113;
constexpr int Key_Right = 114;

constexpr float Near = 0.1f;
constexpr float Far = 20.0f;

template<typename T, i2 size>
struct Image
{
	constexpr i2 GetSize() const { return size; }
	void Clear() { memset(pix, 0, sizeof(pix)); }
	void Clear(T value) { for (int i = 0; i < size.y; ++i) for (int j = 0; j < size.x; ++j) pix[i][j] = value; }
	void Set(i2 coords, T color) { pix[coords.y][coords.x] = color; }
	void Add(i2 coords, T color) { pix[coords.y][coords.x] += color; }
	T Get(i2 coords) const { return pix[coords.y][coords.x]; }
	T& At(i2 coords) { return pix[coords.y][coords.x]; }
	template<typename CopyFunc>
	void CopyFrom(auto& image, CopyFunc&& copy_fn)
	{
		for (int i = 0; i < size.y; ++i)
			for (int j = 0; j < size.x; ++j)
				pix[i][j] = copy_fn(image.pix[i][j]);
	}

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

	std::array<v4, 3> GetEdgeVertsExcept(int index) const
	{
		if (index == 0)
			return std::array<v4, 3>{verts[1], verts[3], verts[2]};
		if (index == 1)
			return std::array<v4, 3>{verts[2], verts[0], verts[3]};
		if (index == 2)
			return std::array<v4, 3>{verts[3], verts[1], verts[0]};
		if (index == 3)
			return std::array<v4, 3>{verts[0], verts[2], verts[1]};
		return {};
	}

	v4 GetNormal() const
	{
		v4 n1 = (verts[1] - verts[0]).Cross(verts[2] - verts[0]);
		v4 n2 = (verts[2] - verts[0]).Cross(verts[3] - verts[0]);
		return (fabs(n1.Dot(v4{1.0f, 1.0f, 1.0f, 1.0f})) > fabs(n2.Dot(v4{1.0f, 1.0f, 1.0f, 1.0f}))) ? n1 : n2;
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
	Input() { memset(keys, 0, sizeof(keys)); memset(keys_time, 0, sizeof(keys_time)); memset(keys_has_been_pressed, 0, sizeof(keys_has_been_pressed)); }
	void Update(float dt) 
	{
		for (int i = 0; i < 256; ++i)
			if (keys[i])
				keys_time[i] += dt;
	}
	bool ConsumeKey(int key)
	{
		if (keys_has_been_pressed[key] & 0x2)
		{
			keys_has_been_pressed[key] &= 0x1;
			return true;
		}
		return false;
	}
	bool keys[256];
	uchar keys_has_been_pressed[256];
	float keys_time[256];
};

struct Camera
{
	v4 GetForward() const { return {-sinf(angle), 0.0f, cosf(angle), 0.0f}; }
	v4 GetRight() const { return {cosf(angle), 0.0f, sinf(angle), 0.0f}; }

	void Update(const Input& input, float dt)
	{
		pos += GetForward() * ((float)input.keys[Key_W] - (float)input.keys[Key_S]) * speed * dt;
		pos += GetRight() * ((float)input.keys[Key_D] - (float)input.keys[Key_A]) * speed * dt;
		angle += rotation_speed * dt * ((float)input.keys[Key_Left] - (float)input.keys[Key_Right]);
	}

	Mtx GetViewMatrix() const
	{
		return
			{
				{cosf(angle), 0.0f, sinf(angle), -pos.x * cos(angle) - pos.z * sin(angle)},
				{0.0f, 1.0f, 0.0f, 										-pos.y														 },
				{-sinf(angle), 0.0f, cosf(angle), pos.x * sin(angle) - pos.z * cos(angle)},
				{0.0f, 0.0f, 0.0f, 1.0f}
			};
	}

	v4 pos = v4::Zero();
	float angle = 0.0f;

	float speed = 1.0f;
	float rotation_speed = 1.15f;
};

void GenerateScene(std::vector<Quad>& quads, float dt)
{
	static float t = 0.0f;
	t += dt;
	quads.clear();

	for (int i = 0; i < 6; ++i)
	for (int j = 0; j < 6; ++j)
	{
		float angle = t + i + j * 0.5;
		Mtx mtx
		{
			v4{cosf(angle), 0.0f, -sinf(angle), -0.0f + j},
			v4{0.0f, 1.0f, 0.0f, 0.0f},
			v4{sinf(angle), 0.0f, cosf(angle), 1.0f - (float)i},
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
			v4{0.0f, 1.0f, 0.0f}
		});
		quads.emplace_back(mtx*Quad
		{ 
			v4{-0.15f, 0.5f, -0.15f},
			v4{-0.15f, -0.5f, -0.15f},
			v4{0.15f, -0.5f, -0.15f},
			v4{0.15f, 0.5f, -0.15f},
			v4{1.0f, 1.0f, 0.0f}
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
	quads.emplace_back(Quad
	{
		v4{-5.0f, -0.5f, 5.0f},
		v4{-5.0f, -0.5f, -5.0f},
		v4{5.0f, -0.5f, -5.0f},
		v4{5.0f, -0.5f, 5.0f},
		v4{0.5f, 0.5f, 0.5f}
	});
}


enum class TransformResult
{
	Discard,
	Ok,
	Split
};

TransformResult TransformQuad(Quad& q, const Mtx& view_proj, Quad& split_quad)
{
	if (q.verts[0].z < 0.0f &&
			q.verts[1].z < 0.0f &&
			q.verts[2].z < 0.0f &&
			q.verts[3].z < 0.0f)
	{
		return TransformResult::Discard;
	}
	Quad qcopy = q;
	for (int i = 0; i < 4; ++i)
	{
		v4& vert = q.verts[i];
		float depth = vert.z;
		if (depth < 0.0f)
		{
			auto others = qcopy.GetEdgeVertsExcept(i);
			using Intersection = std::optional<v4>;
			auto&& f = [depth](v4 a, v4 b) -> Intersection
			{
				v4 e = b - a;
				constexpr float Epsilon = 1e-6f;
				if (fabs(e.z) < Epsilon)
					return {};

				v4 c = {a.x - depth * e.x / e.z, a.y - depth * e.y / e.z, 0.0f, Near};
				if ((c - b).Dot(c - a) < 0.0f)
					return c;
				return {};
			};
			Intersection intersects[2] = {f(vert, others[0]), f(vert, others[1])};
			int num = std::ranges::count(intersects, true, &Intersection::has_value);
			if (num == 0)
			{
				vert = f(vert, others[2]).value();
			}
			if (num == 1)
			{
				vert = intersects[0] ? intersects[0].value() : intersects[1].value();
			}
			if (num == 2)
			{
				split_quad = qcopy;
				v4 middle = f(vert, others[2]).value();
				vert = intersects[0].value();
				q.verts[(i + 3) % 4] = middle;
				split_quad.verts[i] = intersects[1].value(); 
				split_quad.verts[(i + 1) % 4] = middle;
				return TransformResult::Split;
			}
		}
	}
	return TransformResult::Ok;
}

void PerspectiveTransformQuad(Quad& q)
{
	for (v4& vert : q.verts)
	{
		vert /= vert.w;
	}
	
	Mtx to_screen_ctr
	{
		v4{0.5f, 0.0f, 0.0f, 0.5f},
		v4{0.0f, 0.5f, 0.0f, 0.5f},
		v4{0.0f, 0.0f, 1.0f, 0.0f}
	};
	q = to_screen_ctr * q;
}

void RasterizeQuad(Quad& q, auto& frame, auto& depth)
{
	v4 normal = q.GetNormal();
	if (normal.z >= 0.0f)
		return;

	i2 top_left = Floor_i2(q.GetMin() * frame.GetSize());
	i2 bottom_right = Ceil_i2(q.GetMax() * frame.GetSize());
	top_left = Clamp(top_left, {0, 0}, frame.GetSize());
	bottom_right = Clamp(bottom_right, {0, 0}, frame.GetSize());
	v4 frame_size_inv{1.0f / frame.GetSize().x, 1.0f / frame.GetSize().y};
	for (int y = top_left.y; y < bottom_right.y; ++y)
		for (int x = top_left.x; x < bottom_right.x; ++x)
		{
			v4& frame_depth = depth.At({x, y}); // v4 is for 4 float MSAA values
			
			const float ms[][2] = {{0.05f, 0.25f}, {0.45f, 0.05f}, {0.55f, 0.95f}, {0.95f, 0.75f}};
			auto& pixel_color = frame.At({x, y});
			for (int i = 0; i < 4; ++i)
			{
				v4 p = v4{x + ms[i][0], y + ms[i][1]} * frame_size_inv; 
				float e[4] {q.GetEdgeValue(p, 0),
										q.GetEdgeValue(p, 1),
										q.GetEdgeValue(p, 2),
										q.GetEdgeValue(p, 3)};
				if (e[0] > 0 && e[1] > 0 && e[2] > 0 && e[3] > 0)
				{
					float pixel_depth = (q.verts[0].Dot(normal) - (p).Dot(normal)) / normal.z;
					if (pixel_depth < frame_depth[i])
					{
						pixel_color[i] = q.color;
						frame_depth[i] = pixel_depth;
					}
				}
			}
		}
}

int main()
{
	Image<b4, {128, 128}> frame;
	Image<std::array<v4, 4>, frame.GetSize()> frame_ms; // 4 colors for multisampling
	Image<v4, frame.GetSize()> depth;

	constexpr int scale_factor = 8;
	GameWindow window(frame, scale_factor);

	Input input;
	auto&& key_fn = [&input](int key, bool is_pressed)
	{
		input.keys[key] = is_pressed;
		if (!is_pressed)
		{
			input.keys_time[key] = 0.0f;
			input.keys_has_been_pressed[key] &= 0x2;
		}
		else if (input.keys_has_been_pressed[key] == 0u)
			input.keys_has_been_pressed[key] = 3;
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

			camera.Update(input, dt);
			
			GenerateScene(quads, dt);
			
			Mtx view = camera.GetViewMatrix();
			
			Mtx projection 
			{
				{1.0f, 0.0f, 0.0f, 0.0f},
				{0.0f, -1.0f, 0.0f, 0.0f},
				{0.0f, 0.0f, Far / (Far - Near), -Far * Near / (Far - Near)},
				{0.0f, 0.0f, 1.0f, 0.0f} 
			};
			Mtx view_proj = projection * view; 
			frame.Clear();
			frame_ms.Clear();
			depth.Clear(v4{Far});
			for (Quad q : quads)
			{
				Quad split_quad;
				q = view_proj * q;
				auto transform_result = TransformQuad(q, view_proj, split_quad);
				if (transform_result != TransformResult::Discard)
				{
					PerspectiveTransformQuad(q);
					RasterizeQuad(q, frame_ms, depth);
					if (transform_result == TransformResult::Split)
					{
						PerspectiveTransformQuad(split_quad);
						RasterizeQuad(split_quad, frame_ms, depth);
					}
				}
			}
			frame.CopyFrom(frame_ms, 
					[](const std::array<v4, 4>& pix) 
					{return b4((pix[0] + pix[1] + pix[2] + pix[3]) * (255.0f / 4.0f));});
		};
	
	window.Run(frame_fn, key_fn);
	
	return 0;
}
