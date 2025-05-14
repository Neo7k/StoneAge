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

	b4 pix[size.y][size.x];
};

int main()
{
	Image<{8, 8}> frame;

	GameWindow window(frame, 64);

	auto&& frame_fn =
		[&]
		()
		{
			frame.Clear();

			frame.pix[1][1] = {255, 0, 0, 0};
			frame.pix[2][2] = {255, 127, 0, 0};
			frame.pix[3][3] = {255, 255, 0, 0};
			frame.pix[4][4] = {0, 255, 0, 0};
			frame.pix[5][5] = {0, 0, 255, 0};
			frame.pix[6][6] = {127, 0, 255, 0};

		};
	
	window.Run(frame_fn);
	
	return 0;
}
