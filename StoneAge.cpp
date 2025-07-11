#include <assert.h>
#include <math.h>
#include <string.h>

#include <algorithm>
#include <iostream>
#include <thread>
#include <chrono>
#include <thread>
#include <barrier>
#include <vector>
#include <ranges>
#include <numeric>
#include <fstream>
#include <random>

using uchar = unsigned char;
using uint = unsigned int;
namespace view = std::ranges::views;

#include "Math.h"
#include "GameWindow.h"

#define LOG(...) \
		std::cout << std::format(__VA_ARGS__) << std::endl;
#define CLOG(condition, ...) \
	if (condition) \
		std::cout << std::format(__VA_ARGS__) << std::endl;

constexpr int Key_W = 25;
constexpr int Key_A = 38;
constexpr int Key_S = 39;
constexpr int Key_D = 40;
constexpr int Key_Left = 113;
constexpr int Key_Right = 114;
constexpr int Key_Space = 65;

constexpr float Near = 0.1f;
constexpr float Far = 20.0f;

bool ltf = false;

template<typename T, i2 size>
struct Image
{
	Image()
	{
		pix = new T[size.x * size.y];
	}
	~Image()
	{
		delete[] pix;
	}

	constexpr i2 GetSize() const { return size; }
	void Clear() { memset(pix, 0, size.x * size.y * sizeof(T)); }
	void Clear(uint thread_id, uint thread_num) 
	{ 
		uint num_pix = size.x * size.y;
		uint chunk_size = DivCeil(num_pix, thread_num);
		chunk_size = std::min(chunk_size, num_pix - chunk_size * thread_id);
		memset(pix + chunk_size * thread_id, 0, chunk_size * sizeof(T)); 
	}
	void Clear(T value) { for (int i = 0; i < size.x * size.y; ++i) pix[i] = value; }
	void Clear(T value, uint thread_id, uint thread_num) 
	{ 
		uint num_pix = size.x * size.y;
		int chunk_size = (num_pix + thread_num - 1) / thread_num;
		for (auto chunk_view : view::iota(0u, num_pix) | view::chunk(chunk_size) | view::drop(thread_id) | view::take(1))
			for (uint i : chunk_view)
				pix[i] = value; 
	}
	void Set(i2 coords, T color) { pix[coords.y * size.x + coords.x] = color; }
	void Add(i2 coords, T color) { pix[coords.y * size.x + coords.x] += color; }
	T Get(i2 coords) const { return pix[coords.y * size.x + coords.x]; }
	T& At(i2 coords) { return pix[coords.y * size.x + coords.x]; }
	template<typename CopyFunc>
	void CopyFrom(auto& image, CopyFunc&& copy_fn)
	{
		for (int i = 0; i < size.x * size.y; ++i)
				pix[i] = copy_fn(image.pix[i]);
	}
	template<typename CopyFunc>
	void CopyFrom(auto& image, CopyFunc&& copy_fn, int thread_id, int thread_num)
	{
		uint num_pix = size.x * size.y;
		int chunk_size = (num_pix + thread_num - 1) / thread_num;
		for (auto chunk_view : view::iota(0u, num_pix) | view::chunk(chunk_size) | view::drop(thread_id) | view::take(1))
			for (uint i : chunk_view)
				pix[i] = copy_fn(image.pix[i]);
	}

	T* pix;
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
	Input() { memset(keys, 0, sizeof(keys)); memset(keys_time, 0, sizeof(keys_time)); memset(keys_has_been_pressed, 0, sizeof(keys_has_been_pressed)); memset(keys_up, 0, sizeof(keys_up)); }
	void Update(float dt) 
	{
		for (int i = 0; i < 256; ++i)
			if (keys[i])
				keys_time[i] += dt;
	}
	void EndFrame()
	{
		for (int i = 0; i < 256; ++i)
			if (keys_up[i])
				keys_time[i] = 0.0f;
		memset(keys_up, 0, sizeof(keys_up));
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
	bool keys_up[256];
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

	v4 GetPos2d() const { return {pos.x, 0.0f, pos.z}; }

	static constexpr float h = 1.0f;
	static constexpr float speed = 1.0f;
	static constexpr float rotation_speed = 1.15f;
	v4 pos = v4{0.0f, h, 0.0f, 1.0f};
	float angle = 0.0f;
};


using Rng = std::mt19937;
using RngF = std::uniform_real_distribution<float>;
using RngI = std::uniform_int_distribution<int>;

struct Model
{
	Model(const char* filename)
	{
		std::ifstream f(filename);
		v4 verts[4] = { {1.0f}, {1.0f}, {1.0f}, {1.0f} };
		v4 color = {1.0f};
		while (f)
		{
			for (int i = 0; i < 4; ++i)
				f >> verts[i].x >> verts[i].y >> verts[i].z;
			f >> color.x >> color.y >> color.z;

			std::string s;
			std::getline(f, s);

			Quad q
			{
				verts[0], verts[1], verts[2], verts[3],
				color
			};

			quads.push_back(q);
		}
	}

	Model& Transform(const Mtx& mtx)
	{
		for (auto& q : quads)
			q = mtx * q;
		return *this;
	}

	Model& Light(v4 light_dir)
	{
		light_dir = light_dir.Normalized();
		for (auto& q : quads)
			q.color = q.color * std::max(q.GetNormal().Normalized().Dot(light_dir), 0.3f);
		return *this;
	}

	Model& Colorize(v4 min, v4 max, Rng& rng)
	{
		for (auto& q : quads)
			if (q.color.x < 0.0f)
				q.color = v4{RngF(min.x, max.x)(rng), RngF(min.y, max.y)(rng), RngF(min.z, max.z)(rng)};
		return *this;
	}

	void CopyTo(std::vector<Quad>& in_quads) const
	{
		in_quads.insert(in_quads.end(), quads.begin(), quads.end());
	}

	std::vector<Quad> quads;
};

struct Map
{
	enum class Cell
	{
		Empty,
		Cactus
	};

	virtual ~Map() = default;
	virtual int GetSize() const = 0; 
	virtual Cell& At(i2 c) = 0;
	virtual const Cell& At(i2 c) const = 0;
};

template<int size>
struct MapT : public Map
{
	
	MapT()
	{
		memset(map, 0, sizeof(map));
		constexpr int map_sz = size;
		constexpr int map_sz2 = map_sz / 2;
		constexpr float ms2f = (float)map_sz2;
		constexpr int num_cactuses = 50;

		Rng rng(2);
		RngI cpick(-map_sz2 + 1, map_sz2 - 1);
		RngF rpick(0.0f, M_PI * 2.0f);
		for (int i = 0; i < num_cactuses; ++i)
		{
			int j = 0;
			while (j++ < map_sz * map_sz)
			{
				i2 coords{cpick(rng), cpick(rng)};
				if (At(coords) == Cell::Empty)
				{
					At(coords) = Cell::Cactus;
					break;
				}
			}
		}
	}

	virtual int GetSize() const override { return size; }
	virtual Cell& At(i2 c) override { c = c + i2{size/2}; if (c.x < 0 || c.x >= size || c.y < 0 || c.y >= size) return invalid; return map[c.y][c.x]; }
	virtual const Cell& At(i2 c) const override { c = c + i2{size/2}; if (c.x < 0 || c.x >= size || c.y < 0 || c.y >= size) return invalid; return map[c.y][c.x]; }

	Cell map[size][size];
	Cell invalid = Cell::Empty;
};

struct Box
{
	Box() = default;
	Box(float w, float h, float d)
	{
		float hw = w / 2.0f; float hh = h / 2.0f; float hd = d / 2.0f;
		quads[0] = Quad{{-hw,-hh,hd},{-hw,-hh,-hd},{-hw,hh,-hd},{-hw,hh,hd}};
		quads[1] = Quad{{hw,-hh,hd},{hw,hh,hd},{hw,hh,-hd},{hw,-hh,-hd}};
		quads[2] = Quad{{-hw,-hh,-hd},{hw,-hh,-hd},{hw,hh,-hd},{-hw,hh,-hd}};
		quads[3] = Quad{{-hw,-hh,hd},{-hw,hh,hd},{hw,hh,hd},{hw,-hh,hd}};
		quads[4] = Quad{{-hw,hh,-hd},{hw,hh,-hd},{hw,hh,hd},{-hw,hh,hd}};
		quads[5] = Quad{{-hw,-hh,-hd},{-hw,-hh,hd},{hw,-hh,hd},{hw,-hh,-hd}};
	}

	void Transform(const Mtx& mtx)
	{
		for (Quad& q : quads)
			q = mtx * q;
	}

	bool Contains(v4 p) const
	{
		for (const Quad& q : quads)
			if (q.GetNormal().Dot(p - q.verts[0]) < 0.0f)
				return false;

		return true;
	}

	Quad quads[6];
};

struct Player;
struct Dino
{
	enum State
	{
		Wander,
		Chase,
		Flee
	};

	Dino(const Model& in_model_wl, const Model& in_model_wr, Rng& in_rng, const auto& in_map)
		: models{in_model_wl, in_model_wr}, rng(in_rng), map(&in_map)
		, bbox{{1.5f, 0.95f, 0.5f},{1.1f, 1.1f, 0.5f}}
	{
		v4 min{0.3f, 0.5f, 0.0f}, max{0.9f, 0.9f, 0.1f};
		auto rngc = rng;
		models[0].Colorize(min, max, rngc);
		rngc = rng;
		models[1].Colorize(min, max, rngc);
		RngI mp(-map->GetSize() / 2, map->GetSize() / 2);
		pos = v4{(float)mp(rng), 0.0f, (float)mp(rng)};
	}

	void Update(float dt, Player& player)
	{
		if (push.Len2() > 0.0f)
		{
			push.y -= 2.5f * dt;
			pos += push * dt;
			if (pos.y < 0.0f)
			{
				pos.y = 0.0f;
				push = v4::Zero();
			}
			assert(pos.w == 1.0f);
			UpdateTransform();
			return;
		}

		if (state == State::Wander)
		{
			UpdateWander(dt, player);
		}
		if (state == State::Chase)
		{
			UpdateChase(dt, player);
		}
		if (state == State::Flee)
		{
			UpdateFlee(dt, player);
		}
		
		if (pos.Dist2(target_pos) > 1.5f)
		{
			angle = atan2f(target_pos.z - pos.z, target_pos.x - pos.x);
			pos += GetForward() * speed * dt; 
			UpdateTransform();
			i2 ipos{(int)pos.x, (int)pos.z};
			if (map->At(ipos) == Map::Cell::Cactus)
			{
				if (pos.Dist2(v4{ipos.x + 0.5f, 0.0f, ipos.y + 0.5f}) < 1.5f)
					Push(GetForward() * speed + v4(0.0f, 2.5f, 0.0f, 0.0f));
			}

			if ((lr_t += dt) > lr_period)
			{
				lr_t = 0.0f;
				model_id = !model_id;
			}
		}
	}
	
	void UpdateWander(float dt, Player& player);
	void UpdateChase(float dt, Player& player);
	void UpdateFlee(float dt, Player& player);

	void UpdateTransform()
	{
		transform = Mtx::Translate(pos) * Mtx::RotateY(angle);
		memcpy(bboxt, bbox, sizeof(bbox));
		bboxt[0].Transform(transform * Mtx::Translate({-0.45f, 0.80f, 0.0f}));
		bboxt[1].Transform(transform * Mtx::Translate({0.40f, 1.55f, 0.0f}));
	}

	void Push(v4 in_push)
	{
		push = in_push;
	}

	void Hit(float s)
	{
		if (state == State::Wander)
		{
			aggression = 1.0f;
			state = State::Chase;
		}
		else if (state == State::Chase)
		{
			aggression -= s;
			if (aggression < 0.0f)
			{
				state = State::Flee;
				flee_time = 8.0f;
			}
		}
		else if (state == State::Flee)
		{
			flee_time += 2.0f;
		}
	}

	v4 GetForward() const { return {cosf(angle), 0.0f, sinf(angle), 0.0f}; }

	void CopyTo(std::vector<Quad>& in_quads) const
	{
		Model(models[model_id]).Transform(transform).Light({0.5f, -0.5f, 0.5f}).CopyTo(in_quads);
#define DRAW_HITBOX 0
#if DRAW_HITBOX
		for (int i = 0; i < 2; ++i)
		for (auto& q : bboxt[i].quads)
			in_quads.push_back(q);
#endif
	}

	v4 pos{0.0f, 0.0f, 0.0f};
	float angle = 0.0f;
	float speed = 1.0f;
	float wander_time = 1.0f;
	State state = State::Wander;
	float aggression = 0.0f;
	float flee_time = 0.0f;
	float lr_period = 0.5f;
	float lr_t = 0.0f;
	v4 target_pos;
	uint model_id = 0u;
	Model models[2];
	Rng rng;
	const Map* map = nullptr;
	v4 push = v4::Zero();
	Box bbox[2];
	Box bboxt[2];
	Mtx transform;
};

struct Rock
{
	Rock(const Model& in_model, Rng* in_rng)
		: model(&in_model), rng(in_rng) {}

	void Update(float dt, std::vector<Dino>& dinos)
	{
		if (!in_game || sleeps)
			return;

		constexpr float floor_y = 0.1f;
		if (vel.Len2() < sleep_vel2 && fabs(pos.y - floor_y) < 0.1f)
		{
			sleeps = true;
			return;
		}

		vel -= v4{0.0f, g, 0.0f, 0.0f} * dt;
		pos += vel * dt;
		angle_vel -= 0.15f * angle_vel * dt; 
		angle_x += angle_vel * dt;
		angle_z += angle_vel * 0.5f * dt;
		if (pos.y < floor_y)
		{
			bounce_num++;
			pos.y = floor_y;
			vel /= bounce;
			vel.y = -vel.y;
			if (bounce_num = 1)
				angle_vel = angle_speed * 2.0f;
		}
		else 
		{
			for (auto& dino : dinos)
				if (dino_hit != &dino && pos.Dist2(dino.pos) < 8.5f)
				{
					if (dino.bboxt[0].Contains(pos) || dino.bboxt[1].Contains(pos))
					{
						float hit_s = 0.05f * sqrtf(vel.Len2());
						v4 hit = vel * hit_s;
						hit.y = 3.5f * hit_s;
						dino.Push(hit);
						dino.Hit(hit_s);
						dino_hit = &dino;
						vel = Mtx::RotateY(RngF(M_PI / 2.0f, M_PI * 3/4)(*rng)) * vel * 0.2f;
						vel.y = RngF(1.0f, 3.0f)(*rng);
						break;
					}
				}
		}
	}

	void Throw(v4 p, v4 v)
	{
		in_game = true;
		sleeps = false;
		pos = p;
		vel = v;
		angle_vel = angle_speed;
		bounce_num = 0;
		dino_hit = nullptr;
	}

	void CopyTo(std::vector<Quad>& in_quads) const
	{
		if (in_game)
			Model(*model).Transform(Mtx::Translate(pos) * Mtx::RotateX(angle_x) * Mtx::RotateZ(angle_z)).Light({0.5f, -0.5f, 0.5f}).CopyTo(in_quads);
	}

	v4 pos = v4::Zero();
	v4 vel = v4::Zero();
	const Model* model = nullptr;
	const Dino* dino_hit = nullptr;
	Rng* rng = nullptr;
	float angle_x = 0.0f;
	float angle_z = 0.0f;
	float angle_speed = 3.0f;
	float angle_vel = 3.0f;
	float g = 2.0f;
	float bounce = 2.0f;
	float sleep_vel2 = 0.5f;
	int bounce_num = 0;
	bool in_game = false;
	bool sleeps = true;
};

struct Player
{
	Player(Input& in_input, const Model& rock_model)
		: input(&in_input), rocks(8, Rock(rock_model, &rng))
	{}

	void Update(float dt, std::vector<Dino>& dinos)
	{
		if (push.Len2() > 0.0f)
		{
			push.y -= 2.5f * dt;
			camera.pos += push * dt;
			if (camera.pos.y < camera.h)
			{
				camera.pos.y = camera.h;
				push = v4::Zero();
			}
		}
		else
		{
			camera.Update(*input, dt);
			if ((rock_cd -= dt) < 0.0f && input->keys_up[Key_Space])
			{
				rock_cd = 0.5f;
				rock_id = (rock_id + 1) % rocks.size();
				v4 throw_v = camera.GetForward() * 5.5f * std::clamp(input->keys_time[Key_Space] * 2.0f, 0.5f, 1.5f);
				throw_v.y = std::clamp(input->keys_time[Key_Space] * 2.0f, 0.0f, 2.5f);
				rocks[rock_id].Throw(camera.pos + camera.GetRight() * 0.25f, throw_v);
			}
		}
		for (auto& rock : rocks)
			rock.Update(dt, dinos);
	}

	void Push(v4 push_v)
	{
		push = push_v;
	}
	
	void CopyTo(std::vector<Quad>& in_quads) const
	{
		for (auto& rock : rocks)
			rock.CopyTo(in_quads);
	}

	Camera camera;
	Input* input = nullptr;
	Rng rng;
	std::vector<Rock> rocks;
	int rock_id = 0;
	float rock_cd = 0.0f;
	v4 push = v4::Zero();
};

void Dino::UpdateChase(float dt, Player& player)
{
	if (player.camera.pos.y != player.camera.h)
		return;

	target_pos = player.camera.GetPos2d();
	if (target_pos.Dist2(pos) < 2.0f)
	{
		player.Push((target_pos - pos).Normalized() * 2.0f + v4{0.0f, 2.0f, 0.0f, 0.0f});
	}
}

void Dino::UpdateFlee(float dt, Player& player)
{
	if ((flee_time -= dt) < 0.0f)
	{
		state = State::Wander;
		aggression = 0.0f;
		return;
	}

	v4 player_pos = player.camera.GetPos2d();
	v4 way_out = pos - player_pos;
	way_out = way_out * 15.0f;
	target_pos = player_pos + way_out;
	if (pos.Dist2(player_pos) > 225.0f)
	{
		state = State::Wander;
		aggression = 0.0f;
	}
}

void Dino::UpdateWander(float dt, Player& player)
{
	if ((wander_time -= dt) < 0.0f)
	{
		wander_time = RngF(1.5f, 8.0f)(rng);
		RngF wr(-5.0f, 5.0f);
		target_pos = pos + v4{wr(rng), 0.0f, wr(rng)};
	}
}

void GenerateStaticScene(std::vector<Quad>& quads, const Map& map, const Model& cactus)
{
	quads.clear();

	Rng rng;
	int ms2 = map.GetSize() / 2 - 1;
	float ms2f = map.GetSize() / 2.0f;
	for (int i = -ms2; i < ms2; ++i)
	for (int j = -ms2; j < ms2; ++j)
	{
		if (map.At({i, j}) == Map::Cell::Cactus)
			Model(cactus).Transform(Mtx::Translate({i + 0.5f, 0.0f, j + 0.5f}) * Mtx::RotateY(RngF(0.0f, M_PI * 2)(rng)) * Mtx::Scale(0.5f)).Colorize({0.2f, 0.5f, 0.0f}, {0.6f, 0.9f, 0.1f}, rng).Light({0.5f, -0.5f, 0.5f}).CopyTo(quads);
	}
}

enum class TransformResult
{
	Discard,
	Ok,
	Split
};

TransformResult TransformQuad(Quad& q, const Mtx& view_proj, Quad& split_quad)
{
	q = view_proj * q;
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
				if ((c - b).Dot(c - a) <= 0.0f)
					return c;
				return {};
			};
			Intersection intersects[2] = {f(vert, others[0]), f(vert, others[1])};
			int num = std::ranges::count(intersects, true, &Intersection::has_value);
			if (num == 0)
			{
				auto cross = f(vert, others[2]);
				//assert(cross);
				if (cross)
					vert = cross.value();
				else
					return TransformResult::Discard;
			}
			if (num == 1)
			{
				vert = intersects[0] ? intersects[0].value() : intersects[1].value();
			}
			if (num == 2)
			{
				split_quad = qcopy;
				auto middle = f(vert, others[2]);
				//assert(middle);
				if (!middle)
					return TransformResult::Discard;

				vert = intersects[0].value();
				q.verts[(i + 3) % 4] = middle.value();
				split_quad.verts[i] = intersects[1].value(); 
				split_quad.verts[(i + 1) % 4] = middle.value();
				return TransformResult::Split;
			}
		}
	}
	
	return TransformResult::Ok;
}

void PerspectiveTransformQuad(Quad& q)
{
	for (v4& vert : q.verts)
		vert /= vert.w;
	
	Mtx to_screen_ctr
	{
		v4{0.5f, 0.0f, 0.0f, 0.5f},
		v4{0.0f, 0.5f, 0.0f, 0.5f},
		v4{0.0f, 0.0f, 1.0f, 0.0f}
	};
	q = to_screen_ctr * q;
}

void RasterizeQuad(Quad& q, auto& frame, auto& depth, i2 tiles)
{
	v4 normal = q.GetNormal();

	int tile = tiles.x;
	int num_tiles = tiles.y;
	int ntsqrt = (int)sqrtf(num_tiles);
	i2 ntiles{num_tiles / ntsqrt, ntsqrt};
	i2 tile_xy{tile % ntiles.x, tile / ntiles.x};

	i2 top_left = Floor_i2(q.GetMin() * frame.GetSize());
	i2 bottom_right = Ceil_i2(q.GetMax() * frame.GetSize());

	top_left = Clamp(top_left, (frame.GetSize() * tile_xy) / ntiles, (frame.GetSize() * (tile_xy + i2{1, 1})) / ntiles);
	bottom_right = Clamp(bottom_right, (frame.GetSize() * tile_xy) / ntiles, (frame.GetSize() * (tile_xy + i2{1, 1})) / ntiles);

	v4 frame_size_inv{1.0f / frame.GetSize().x, 1.0f / frame.GetSize().y};
	for (int y = top_left.y; y < bottom_right.y; ++y)
		for (int x = top_left.x; x < bottom_right.x; ++x)
		{
			v4& frame_depth = depth.At({x, y}); // v4 is for 4 float MSAA values
			
			constexpr float ms[][2] = {{0.05f, 0.25f}, {0.45f, 0.05f}, {0.55f, 0.95f}, {0.95f, 0.75f}};
			auto& pixel_color = frame.At({x, y});
			for (int i = 0; i < 4; ++i)
			{
				v4 p = v4{x + ms[i][0], y + ms[i][1]} * frame_size_inv; 
				float e[4] {q.GetEdgeValue(p, 0),
										q.GetEdgeValue(p, 1),
										q.GetEdgeValue(p, 2),
										q.GetEdgeValue(p, 3)};
				if (e[0] >= 0.0f && e[1] >= 0.0f && e[2] >= 0.0f && e[3] >= 0.0f)
				{
					float pixel_depth = (q.verts[0].Dot(normal) - p.Dot(normal)) / normal.z;
					if (pixel_depth < frame_depth[i])
					{
						pixel_color[i] = q.color;
						frame_depth[i] = pixel_depth;
					}
				}
			}
		}
}

template<int num_frames>
struct PerfTracker
{
	PerfTracker()
	{
		times.fill(0.0f);	
	}
	void Update(float dt)
	{
		static uint frame_index = 0;
		++frame_index;
		float dt_ms = dt * 1000.0f;
		times[frame_index % num_frames] = dt_ms;
		auto times_cpy = times;
		std::sort(times_cpy.begin(), times_cpy.end());
		float mean = times_cpy[1 + (num_frames / 2)];
		float ratio = mean / last;
		if (ratio < 0.7f || ratio > 1.3f)
		{
			LOG("frame time: {} ms", mean);
			last = mean;
		}
	}

	std::array<float, num_frames> times;
	float last = 0.0f;
};

struct MeasureTime
{
	MeasureTime(const char* in_comment)
		: comment(in_comment)
	{
		time_point = std::chrono::high_resolution_clock::now();
	}
	MeasureTime(const char* in_comment, int in_thread_id, int in_thread_count)
		: comment(in_comment)
		, thread_id(in_thread_id)
		, thread_count(in_thread_count)
	{
		time_points[thread_id] = std::chrono::high_resolution_clock::now();
	}

	~MeasureTime()
	{
		if (thread_count == 0)
		{
			std::chrono::duration<float, std::milli> duration = std::chrono::high_resolution_clock::now() - time_point;
			CLOG(ltf, "{} = {}", comment, duration.count());
			return;
		}

		finish_time_points[thread_id] = std::chrono::high_resolution_clock::now();
		arrived_threads_count.fetch_add(1);
		if (arrived_threads_count.load() == thread_count)
		{
			std::stringstream ss;
			ss << comment << " ";
			for (int i = 0; i < thread_count; ++i)
			{
				std::chrono::duration<float, std::milli> duration = finish_time_points[i] - time_points[i];
				ss << duration.count() << " | ";
			}
			ss << std::endl;
			CLOG(ltf, "{}", ss.str());
			arrived_threads_count.store(0);
		}
	}
	const char* comment = nullptr;
	int thread_id = 0;
	int thread_count = 0;
	std::chrono::high_resolution_clock::time_point time_point;
	static std::atomic_int arrived_threads_count;
	static std::chrono::high_resolution_clock::time_point time_points[64];
	static std::chrono::high_resolution_clock::time_point finish_time_points[64];
};
std::chrono::high_resolution_clock::time_point MeasureTime::time_points[64];
std::chrono::high_resolution_clock::time_point MeasureTime::finish_time_points[64];
std::atomic_int MeasureTime::arrived_threads_count{0};

struct Jobs
{
	Jobs()
		: size(std::thread::hardware_concurrency())
	{}

	template<typename ParallelFunc>
	void Run(ParallelFunc&& parallel_func)
	{
		threads.reserve(size);
		
		for (int i = 0; i < size; ++i)
			threads.emplace_back([i, this, parallel_func]() { while (!done.test()) parallel_func(i); } );
	}

	template<typename OnEndFunc>
	void End(OnEndFunc&& on_end_func)
	{
		done.test_and_set();
		on_end_func();
	}

	~Jobs()
	{
		for (auto& thread : threads)
			thread.join();
	}

	int GetSize() const { return size; }

	int size = 0;
	std::vector<std::thread> threads;
	std::atomic_flag done;
};


int main()
{
	Image<b4, {1024, 1024}> frame;
	Image<std::array<v4, 4>, frame.GetSize()> frame_ms; // 4 colors for multisampling
	Image<v4, frame.GetSize()> depth;

	constexpr int scale_factor = 1;
	GameWindow window(frame, scale_factor);

	Input input;
	auto&& key_fn = [&input](int key, bool is_pressed)
	{
		input.keys[key] = is_pressed;
		if (!is_pressed)
		{
			input.keys_has_been_pressed[key] &= 0x2;
			input.keys_up[key] = true;
		}
		else if (input.keys_has_been_pressed[key] == 0u)
			input.keys_has_been_pressed[key] = 3;
	};

	auto frame_start = std::chrono::steady_clock::now();

	PerfTracker<5> perf;

	Jobs jobs;
	std::vector<Quad> quads;
	std::vector<uint> valid_qids;
	std::vector<Quad> splinters;
	struct ThreadData
	{
		std::vector<uint> valid_qids;
		std::vector<Quad> splinters;
	};
	std::vector<ThreadData> thread_data(jobs.GetSize());

	Mtx projection 
	{
		{1.0f, 0.0f, 0.0f, 0.0f},
		{0.0f, -1.0f, 0.0f, 0.0f},
		{0.0f, 0.0f, Far / (Far - Near), -Far * Near / (Far - Near)},
		{0.0f, 0.0f, 1.0f, 0.0f} 
	};
	Mtx view_proj;

	MapT<40> map;
	Model cactus("cactus.qobj");
	std::vector<Quad> static_scene;
	GenerateStaticScene(static_scene, map, cactus);
	Model dino_model_wl("dino_wl.qobj");
	Model dino_model_wr("dino_wr.qobj");
	constexpr int dinos_num = 16;
	std::vector<Dino> dinos;
	for (int i = 0; i < dinos_num; ++i)
	{
		Rng rng(i);
		dinos.emplace_back(dino_model_wl, dino_model_wr, rng, map);
	}
	Model rock_model("rock.qobj");
	Player player(input, rock_model);

	std::barrier frame_start_barrier(jobs.GetSize() + 1, [&]()
	{
		if (jobs.done.test())
			return;

		auto now = std::chrono::steady_clock::now();
		std::chrono::duration<float> dt_dur = now - frame_start;
		frame_start = now;
		float dt = dt_dur.count();

		MeasureTime mt("frame start");
		perf.Update(dt);

		input.Update(dt);
		ltf = input.ConsumeKey(24);
		player.Update(dt, dinos);
		view_proj = projection * player.camera.GetViewMatrix(); 

		quads.clear();
		quads.insert(quads.end(), static_scene.begin(), static_scene.end());
		for (Dino& dino : dinos)
		{
			dino.Update(dt, player);
			dino.CopyTo(quads);
		}
		player.CopyTo(quads);
		// Floor
		float ms2f = map.GetSize() / 2.0f;
		quads.emplace_back(Mtx::Translate(player.camera.GetPos2d()) * Quad
		{
			v4{-ms2f, 0.0f, ms2f},
			v4{-ms2f, 0.0f, -ms2f},
			v4{ms2f, 0.0f, -ms2f},
			v4{ms2f, 0.0f, ms2f},
			v4{0.45f, 0.3f, 0.10f} * 1.5f
		});

		input.EndFrame();
	});

	std::barrier collect_quads_barrier(jobs.GetSize() + 1, [&]()
	{
		MeasureTime mt("collect quads");
		uint n_vq = std::accumulate(thread_data.begin(), thread_data.end(), 0u, [](uint v, auto& s) {return v + s.valid_qids.size();});
		uint n_s = std::accumulate(thread_data.begin(), thread_data.end(), 0u, [](uint v, auto& s) {return v + s.splinters.size();});
		valid_qids.clear();
		splinters.clear();
		valid_qids.reserve(n_vq);
		splinters.reserve(n_s);
		for (auto& td : thread_data)
		{
			valid_qids.insert(valid_qids.end(), td.valid_qids.begin(), td.valid_qids.end());
			splinters.insert(splinters.end(), td.splinters.begin(), td.splinters.end());
		}
	});

	std::barrier copy_image_barrier(jobs.GetSize() + 1);
	std::barrier end_frame_barrier(jobs.GetSize() + 1);

	jobs.Run([&](int thread_id)
	{
		frame_start_barrier.arrive_and_wait();
		if (jobs.done.test())
			return;

		int nthreads = jobs.GetSize();
		int chunk_size = (quads.size() + nthreads - 1) / nthreads;

		// Can do without barriers here because the code below (until the nearest barrier) doesn't touch this data
		{
			frame.Clear(thread_id, nthreads);
			const v4 sky_blue = v4{0.0f, 0.63f, 0.9f};
			frame_ms.Clear({sky_blue, sky_blue, sky_blue, sky_blue}, thread_id, nthreads);
			depth.Clear(v4{Far}, thread_id, nthreads);
		}

		{
			MeasureTime mt("transform quads", thread_id, nthreads);
			thread_data[thread_id].valid_qids.clear();
			thread_data[thread_id].splinters.clear();
			for (auto chunk_view : view::iota(0, (int)quads.size()) | view::chunk(chunk_size) | view::drop(thread_id) | view::take(1))
			{
				for (uint qid : chunk_view)
				{
					auto& quad = quads[qid];
					
					Quad split_quad;
					auto transform_result = TransformQuad(quad, view_proj, split_quad);
					if (transform_result == TransformResult::Discard)
						continue;

					PerspectiveTransformQuad(quad);
					if (quad.GetNormal().z > 0.0f)
						continue;

					thread_data[thread_id].valid_qids.push_back(qid);
					if (transform_result == TransformResult::Split)
					{
						PerspectiveTransformQuad(split_quad);
						thread_data[thread_id].splinters.push_back(split_quad);
					}
				}
			}
		}

		collect_quads_barrier.arrive_and_wait();

		{
			MeasureTime mt("rasterize quads", thread_id, nthreads);
			for (uint qid : valid_qids)
				RasterizeQuad(quads[qid], frame_ms, depth, {thread_id, nthreads});
			for (Quad& quad : splinters)
				RasterizeQuad(quad, frame_ms, depth, {thread_id, nthreads});
		}

		copy_image_barrier.arrive_and_wait();
		{
			MeasureTime mt("copy image", thread_id, nthreads);
			frame.CopyFrom(frame_ms, 
			[](const std::array<v4, 4>& pix) 
			{return b4((pix[0] + pix[1] + pix[2] + pix[3]) * (255.0f / 4.0f));},
			thread_id, nthreads);
		}

		end_frame_barrier.arrive_and_wait();
	});


	auto&& frame_fn = [&]()
	{
		MeasureTime mt("total frame");
		frame_start_barrier.arrive_and_wait();
		collect_quads_barrier.arrive_and_wait();
		copy_image_barrier.arrive_and_wait();
		end_frame_barrier.arrive_and_wait();
	};
	
	window.Run(frame_fn, key_fn);
	jobs.End([&]() { frame_start_barrier.arrive_and_wait(); });
	
	return 0;
}

