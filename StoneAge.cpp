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

int main()
{
	Image<{32, 32}> frame;

	GameWindow window(frame, 32);

	auto&& frame_fn =
		[&]
		()
		{
			frame.Clear();

			i2 top_left{8, 8};
			for (int x = 0; x < 16; ++x)
				for (int y = 0; y < 16; ++y)
					frame.Set(top_left + i2{x, y}, {255, 0, 0, 0});
		};
	
	window.Run(frame_fn);
	
	return 0;
}
