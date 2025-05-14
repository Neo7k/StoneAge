#pragma once

#include <X11/Xlib.h>
#include <X11/extensions/XShm.h>
#include <sys/ipc.h>
#include <sys/shm.h>

#include "Math.h"

constexpr int Key_Escape = 9;

template<typename TImage>
struct GameWindow
{
	GameWindow(TImage& in_frame, int in_scale_factor)
		: frame(in_frame)
		, scale_factor(in_scale_factor)
	{
		xdisplay = XOpenDisplay(nullptr);
		assert(xdisplay);

		screen = DefaultScreen(xdisplay);
		XSetWindowAttributes attributes = 
		{
			.background_pixel = BlackPixel(xdisplay, screen),
			.event_mask = ExposureMask | KeyPressMask
		};

		window_size = in_frame.GetSize() * scale_factor;
		xwindow = XCreateWindow(xdisplay,
				DefaultRootWindow(xdisplay),
				0, 0, window_size.x, window_size.y, 0,
				DefaultDepth(xdisplay, screen),
				InputOutput,
				DefaultVisual(xdisplay, screen),
				CWBackPixel | CWEventMask,
				&attributes);
		XMapWindow(xdisplay, xwindow);

		wm_delete_window = XInternAtom(xdisplay, "WM_DELETE_WINDOW", false);
		XSetWMProtocols(xdisplay, xwindow, &wm_delete_window, 1);

		xframe = XShmCreateImage(xdisplay,
				DefaultVisual(xdisplay, screen), DefaultDepth(xdisplay, screen),
				ZPixmap, nullptr, &shminfo, window_size.x, window_size.y);
		shminfo.shmid = shmget(IPC_PRIVATE, xframe->bytes_per_line * xframe->height, IPC_CREAT | 0777);
		xframe->data = reinterpret_cast<char*>(shmat(shminfo.shmid, 0, 0));
		shminfo.shmaddr = xframe->data;
		shminfo.readOnly = false;
		XShmAttach(xdisplay, &shminfo);
		XSync(xdisplay, false);
	}

	~GameWindow()
	{
		XUnmapWindow(xdisplay, xwindow);
		XDestroyWindow(xdisplay, xwindow);
		XCloseDisplay(xdisplay);
	}

	template<typename TFrameFn>
	void Run(TFrameFn&& frame_fn)
	{
		while (!done)
		{
			while (XPending(xdisplay))
			{
				XEvent event;
				XNextEvent(xdisplay, &event);
				switch (event.type)
				{
					case KeyPress:
						if (event.xkey.keycode == Key_Escape)
							done = true;
						break;
					case ClientMessage:
						if (event.xclient.data.l[0] == wm_delete_window)
							done = true;
						break;
				}
			}

			frame_fn();
			
			i2 frame_size = frame.GetSize();
			for (int i = 0; i < frame_size.y; ++i)
				for (int j = 0; j < frame_size.x; ++j)
					for (int k = 0; k < scale_factor; ++k)
						for (int l = 0; l < scale_factor; ++l)
						{
							int si = (i * scale_factor + k) * frame_size.x * scale_factor * 4;
							int sj = (j * scale_factor + l) * 4;
							xframe->data[si + sj + 0] = frame.pix[i][j].z;
							xframe->data[si + sj + 1] = frame.pix[i][j].y;
							xframe->data[si + sj + 2] = frame.pix[i][j].x;
							xframe->data[si + sj + 3] = frame.pix[i][j].w;
						}

			
			XShmPutImage(xdisplay, xwindow, DefaultGC(xdisplay, screen), xframe,
					0, 0, 0, 0,
					window_size.x, window_size.y, false);
			XFlush(xdisplay);
		}		
	}

	TImage& frame;
	int scale_factor = 1;
	i2 window_size;
	Display* xdisplay = nullptr;
	Window xwindow;
	XImage* xframe = nullptr;
	bool done = false;
	Atom wm_delete_window;
	int screen = 0;
	XShmSegmentInfo shminfo;
};


