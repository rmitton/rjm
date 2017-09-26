// rjm_texbleed.h
//
// Fills in the color of pixels where alpha==0.
//
// To generate the implementation, place this define in exactly one source
// file before including the header:
// #define TEXBLEED_IMPLEMENTATION


// This is free and unencumbered software released into the public domain.
// 
// Anyone is free to copy, modify, publish, use, compile, sell, or
// distribute this software, either in source code form or as a compiled
// binary, for any purpose, commercial or non-commercial, and by any
// means.
// 
// In jurisdictions that recognize copyright laws, the author or authors
// of this software dedicate any and all copyright interest in the
// software to the public domain. We make this dedication for the benefit
// of the public at large and to the detriment of our heirs and
// successors. We intend this dedication to be an overt act of
// relinquishment in perpetuity of all present and future rights to this
// software under copyright law.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
// OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
// ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
// 
// For more information, please refer to <http://unlicense.org/>


#ifndef __RJM_TEXBLEED_H__
#define __RJM_TEXBLEED_H__

// Given an RGBA texture of w*h, finds all pixels
// where alpha==0 and fills in a suitable RGB color for it.
//
// ac = index of the byte within the pixel that contains the alpha channel.
// pixstride = size of one pixel, in bytes.
// rowstride = size of one row of pixels, in bytes.
void rjm_texbleed(unsigned char *pixels, int w, int h, int ac, int pixstride, int rowstride);


//--- Implementation follows ----------------------------------------------

#ifdef TEXBLEED_IMPLEMENTATION

// We search for pixels with alpha greater than this 
// to bleed outwards from.
#define BLEED_THRESHOLD	128

typedef struct { int dx, dy; } TbPoint;

static void bleedcompare(TbPoint *p, int gstride, int offsetx, int offsety)
{
	TbPoint other = p[offsety*gstride + offsetx];
	other.dx += offsetx;
	other.dy += offsety;

	int odist = other.dx*other.dx + other.dy*other.dy;
	int pdist = p->dx*p->dx + p->dy*p->dy;
	if (odist < pdist)
		*p = other;
}

void rjm_texbleed(unsigned char *pixels, int w, int h, int ac, int pixstride, int rowstride)
{
	int gstride = w + 2;
	int cellcount = gstride*(h+2);
	TbPoint *storage = (TbPoint *)malloc(cellcount*sizeof(TbPoint));
	TbPoint *grid = storage + gstride + 1;

	// Initialize to empty.
	for (int n=0;n<cellcount;n++)
		storage[n].dx = storage[n].dy = 9999;
	
	// Fill in the solid pixels.
	int any = 0;
	for (int y=0;y<h;y++)
	{
		for (int x=0;x<w;x++)
		{
			unsigned char *pix = pixels + y*rowstride + x*pixstride;
			if (pix[ac] > BLEED_THRESHOLD) {
				TbPoint *p = &grid[y*gstride+x];
				p->dx = 0;
				p->dy = 0;
				any = 1;
			}
		}
	}

	if (any)
	{
		// Distance field sweep - Pass 0
		for (int y=0;y<h;y++)
		{
			for (int x=0;x<w;x++)
			{
				TbPoint *p = &grid[y*gstride + x];
				bleedcompare(p, gstride, -1,  0);
				bleedcompare(p, gstride,  0, -1);
				bleedcompare(p, gstride, -1, -1);
				bleedcompare(p, gstride,  1, -1);
			}

			for (int x=w-1;x>=0;x--)
			{
				TbPoint *p = &grid[y*gstride + x];
				bleedcompare(p, gstride, 1, 0);
			}
		}

		// Distance field sweep - Pass 1
		for (int y=h-1;y>=0;y--)
		{
			for (int x=w-1;x>=0;x--)
			{
				TbPoint *p = &grid[y*gstride + x];
				bleedcompare(p, gstride,  1, 0);
				bleedcompare(p, gstride,  0, 1);
				bleedcompare(p, gstride, -1, 1);
				bleedcompare(p, gstride,  1, 1);
			}

			for (int x=0;x<w;x++)
			{
				TbPoint *p = &grid[y*gstride + x];
				bleedcompare(p, gstride, -1, 0);
			}
		}

		// Read back the nearest pixels.
		for (int y=0;y<h;y++)
		{
			for (int x=0;x<w;x++)
			{
				TbPoint *p = &grid[y*gstride+x];
				int sx = x + p->dx;
				int sy = y + p->dy;

				// Copy the RGB over.
				unsigned char *src = pixels + sy*rowstride + sx*pixstride;
				unsigned char *dst = pixels + y*rowstride + x*pixstride;
				if (dst[ac] == 0) {
					for (int n=0;n<pixstride;n++)
						dst[n] = src[n];
					dst[ac] = 0;
				}
			}
		}
	}

	free(storage);
}

#endif // TEXBLEED_IMPLEMENTATION
#endif // __RJM_TEXBLEED_H__
