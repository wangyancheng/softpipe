#ifndef __RASTER_H__
#define __RASTER_H__

#include "tgaimage.h"
#include "geometry.h"
#define RASTERMODE 0

extern void rs_line(int x0, int y0, int x1, int y1, TGAImage *image, TGAColor color);
void rs_triangle(Vec2i *pts, TGAImage *image, TGAColor color);

#endif //__RASTER_H__