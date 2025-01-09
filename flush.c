#include "context.h"
#include "flush.h"
#include "debug_img.h"


void
softpipe_flush( struct softpipe_context *softpipe,
                unsigned flags)
{
    unsigned int i;

    #if 0
    draw_flush(softpipe->draw);
    #endif

    #if 0
    if (flags & SP_FLUSH_TEXTURE_CACHE) {
        unsigned sh;

        for (sh = 0; sh < ARRAY_SIZE(softpipe->tex_cache); sh++) {
            for (i = 0; i < softpipe->num_sampler_views[sh]; i++) {
                sp_flush_tex_tile_cache(softpipe->tex_cache[sh][i]);
            }
        }
    }
   #endif

   /* 如果这是一个交换缓冲区，则只需刷新颜色缓冲区。
    *
    * zbuffer 的改变不会被丢弃，而是被保存在缓存中
    * 希望稍后的清除操作能将它们消除。
    */
   for (i = 0; i < softpipe->framebuffer.nr_cbufs; i++)
      if (softpipe->cbuf_cache[i])
         sp_flush_tile_cache(softpipe->cbuf_cache[i]);
   #if 0
   if (softpipe->zsbuf_cache)
      sp_flush_tile_cache(softpipe->zsbuf_cache);
    #endif
   softpipe->dirty_render_cache = false;

   /* Enable to dump BMPs of the color/depth buffers each frame */
#if 1
    if (flags & PIPE_FLUSH_END_OF_FRAME) {
        static unsigned frame_no = 1;
        static char filename[256];
        util_snprintf(filename, sizeof(filename), "cbuf_%u.bmp", frame_no);
        debug_dump_surface_bmp(softpipe, filename, softpipe->framebuffer.cbufs[0]);
        #if 0
        util_snprintf(filename, sizeof(filename), "zsbuf_%u.bmp", frame_no);
        debug_dump_surface_bmp(softpipe, filename, softpipe->framebuffer.zsbuf);
        #endif
        ++frame_no;
    }
#endif

    #if 0
    if (fence)
        *fence = (void*)(intptr_t)1;
    #endif
}