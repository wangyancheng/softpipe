#include "context.h"
#include "state.h"
#include "tile_cache.h"
#include <stdio.h>
#include <stdlib.h>

#ifndef UNUSED
#define UNUSED __attribute__((unused))
#endif

#define MALLOC_STRUCT(T)   (struct T *) malloc(sizeof(struct T))

#define MALLOC(_size)  malloc(_size)

#pragma pack(push,2)
struct bmp_file_header {
   uint16_t bfType;
   uint32_t bfSize;
   uint16_t bfReserved1;
   uint16_t bfReserved2;
   uint32_t bfOffBits;
};
#pragma pack(pop)
struct bmp_info_header {
   uint32_t biSize;
   int32_t biWidth;
   int32_t biHeight;
   uint16_t biPlanes;
   uint16_t biBitCount;
   uint32_t biCompression;
   uint32_t biSizeImage;
   int32_t biXPelsPerMeter;
   int32_t biYPelsPerMeter;
   uint32_t biClrUsed;
   uint32_t biClrImportant;
};

struct bmp_rgb_quad {
   uint8_t rgbBlue;
   uint8_t rgbGreen;
   uint8_t rgbRed;
   uint8_t rgbAlpha;
};

/**
 * Convert float in [0,1] to ubyte in [0,255] with clamping.
 */
static inline ubyte
float_to_ubyte(float f)
{
   /* return 0 for NaN too */
   if (!(f > 0.0f)) {
      return (ubyte) 0;
   }
   else if (f >= 1.0f) {
      return (ubyte) 255;
   }
   else {
      union fi tmp;
      tmp.f = f;
      tmp.f = tmp.f * (255.0f/256.0f) + 32768.0f;
      return (ubyte) tmp.i;
   }
}

void
debug_dump_float_rgba_bmp(const char *filename,
                          unsigned width, unsigned height,
                          float *rgba, unsigned stride)
{
   FILE *stream;
   struct bmp_file_header bmfh;
   struct bmp_info_header bmih;
   unsigned x, y;

   if (!rgba)
      goto error1;

   bmfh.bfType = 0x4d42;
   bmfh.bfSize = 14 + 40 + height*width*4;
   bmfh.bfReserved1 = 0;
   bmfh.bfReserved2 = 0;
   bmfh.bfOffBits = 14 + 40;

   bmih.biSize = 40;
   bmih.biWidth = width;
   bmih.biHeight = height;
   bmih.biPlanes = 1;
   bmih.biBitCount = 32;
   bmih.biCompression = 0;
   bmih.biSizeImage = height*width*4;
   bmih.biXPelsPerMeter = 0;
   bmih.biYPelsPerMeter = 0;
   bmih.biClrUsed = 0;
   bmih.biClrImportant = 0;

   stream = fopen(filename, "wb");
   if (!stream)
      goto error1;

   fwrite(&bmfh, 14, 1, stream);
   fwrite(&bmih, 40, 1, stream);

   y = height;
   while (y--) {
      float *ptr = rgba + (stride * y * 4);
      for (x = 0; x < width; ++x) {
         struct bmp_rgb_quad pixel;
         pixel.rgbRed   = float_to_ubyte(ptr[x*4 + 0]);
         pixel.rgbGreen = float_to_ubyte(ptr[x*4 + 1]);
         pixel.rgbBlue  = float_to_ubyte(ptr[x*4 + 2]);
         pixel.rgbAlpha = float_to_ubyte(ptr[x*4 + 3]);
         fwrite(&pixel, 1, 4, stream);
      }
   }

   fclose(stream);
error1:
   ;
}


void
debug_dump_transfer_bmp(UNUSED struct softpipe_context *pipe,
                        const char *filename,
                        struct pipe_transfer *transfer, UNUSED void *ptr)
{
   float *rgba;

   if (!transfer)
      goto error1;

   rgba = MALLOC(transfer->box.width *
		 transfer->box.height *
		 transfer->box.depth *
		 4*sizeof(float));
   if (!rgba)
      goto error1;

    pipe_get_tile_rgba(transfer, ptr, 0, 0,
                      transfer->box.width, transfer->box.height,
                      rgba);

   debug_dump_float_rgba_bmp(filename,
                             transfer->box.width, transfer->box.height,
                             rgba, transfer->box.width);

   FREE(rgba);
error1:
   ;
}

void
debug_dump_surface_bmp(struct softpipe_context *pipe,
                       const char *filename,
                       struct pipe_surface *surface)
{
   struct pipe_transfer *transfer;
   struct pipe_resource *texture = surface->texture;
   void *ptr;

   ptr = texture->userdata;
   if(1) {
      transfer = MALLOC_STRUCT(pipe_transfer);
      transfer->resource = surface->texture;
      transfer->level = surface->u.tex.level;
      transfer->usage = PIPE_TRANSFER_READ;
      transfer->box.x = 0;
      transfer->box.y = 0;
      transfer->box.z = 0;
      transfer->box.width = surface->width;
      transfer->box.height = surface->height;
      transfer->box.depth = 1;
      transfer->stride = (surface->width + 4) * 4;
      transfer->layer_stride = 0;
      }

   debug_dump_transfer_bmp(pipe, filename, transfer, ptr);

   FREE(transfer);
}