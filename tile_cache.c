#include "tile_cache.h"
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <string.h>

#ifndef UNUSED
#define UNUSED __attribute__((unused))
#endif



#define MALLOC_STRUCT(T)   (struct T *) malloc(sizeof(struct T))
#define MALLOC(_size)  malloc(_size)
#define CALLOC(_count, _size) calloc(_count, _size)
/* Compute the size of an array */
#ifndef ARRAY_SIZE
#  define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))
#endif

static inline int addr_to_clear_pos(union tile_address addr)
{
   int pos;
   pos = addr.bits.layer * (MAX_WIDTH / TILE_SIZE) * (MAX_HEIGHT / TILE_SIZE);
   pos += addr.bits.y * (MAX_WIDTH / TILE_SIZE);
   pos += addr.bits.x;
   return pos;
}
/**
 * Is the tile at (x,y) in cleared state?
 */
static inline uint
is_clear_flag_set(const uint *bitvec, union tile_address addr, unsigned max)
{
   int pos, bit;
   pos = addr_to_clear_pos(addr);
   assert(pos / 32 < (int)max);
   bit = bitvec[pos / 32] & (1 << (pos & 31));
   return bit;
}

/**
 * Mark the tile at (x,y) as not cleared.
 */
static inline void
clear_clear_flag(uint *bitvec, union tile_address addr, unsigned max)
{
   int pos;
   pos = addr_to_clear_pos(addr);
   assert((unsigned)pos / 32 < max);
   bitvec[pos / 32] &= ~(1 << (pos & 31));
}

/**
 * 创建tile cache，包括多个tile数据
 * 
 * \param pipe -[in] 管线上下文  //softpipe_context
 * \return 返回tile cache
 */
struct softpipe_tile_cache *
sp_create_tile_cache( struct softpipe_context *pipe )
{
   struct softpipe_tile_cache *tc;
   uint pos;
   MAYBE_UNUSED int maxTexSize;
   int maxLevels;

   /* sanity checking: max sure MAX_WIDTH/HEIGHT >= largest texture image */
   maxLevels = SP_MAX_TEXTURE_2D_LEVELS;     //15
   maxTexSize = 1 << (maxLevels - 1);          //2^14=16384
   assert(MAX_WIDTH >= maxTexSize);

   STATIC_ASSERT(sizeof(union tile_address) == 4);

   STATIC_ASSERT((TILE_SIZE << TILE_ADDR_BITS) >= MAX_WIDTH);

   tc = CALLOC_STRUCT( softpipe_tile_cache );
   if (tc) {
      tc->pipe = pipe;        //pipe
      for (pos = 0; pos < ARRAY_SIZE(tc->tile_addrs); pos++) {
         tc->tile_addrs[pos].bits.invalid = 1;        //tile地址无效
      }
      tc->last_tile_addr.bits.invalid = 1;         //最近tile地址无效

      /* 这种分配使我们能够保证分配
       * 失败以后永远不会致命
       */
      tc->tile = MALLOC_STRUCT( softpipe_cached_tile );     //创建用于清除的tile
      if (!tc->tile)
      {
         FREE(tc);
         return NULL;
      }

      /* xxx此代码可防止valgrind警告使用非专业化
       * 在渲染之前不会清除表面的程序中的内存。
       * 但是，它在其他情况下破裂了（例如
       * progs/tests/drawbuffers，请参见错误24402）。
       */
#if 0
      /* set flags to indicate all the tiles are cleared */
      memset(tc->clear_flags, 255, sizeof(tc->clear_flags));
#endif
   }
   return tc;
}

/**
 * 销毁tile cache
 * 
 * \param tc -[in] tile cache
 */
void
sp_destroy_tile_cache(struct softpipe_tile_cache *tc)
{
   if (tc) {
      uint pos;

      for (pos = 0; pos < ARRAY_SIZE(tc->entries); pos++) {
         /*assert(tc->entries[pos].x < 0);*/
         FREE( tc->entries[pos] );     //释放tile数据
      }
      FREE( tc->tile );    //释放用于清除的tile

      if (tc->num_maps) {
         //int i;
         FREE(tc->transfer_map);
         FREE(tc->clear_flags);
      }

      FREE( tc );
   }
}

/**
* 指定framebufer surface到tile cache
*
* \param tc -[in] 图块，用于暂存颜色，深度等
* \param ps -[in] framebuffer的surface
*/
void
sp_tile_cache_set_surface(struct softpipe_tile_cache *tc,
                          struct pipe_surface *ps)
{
   //struct softpipe_context *pipe = tc->pipe;
   int i;

   //释放原来已映射的surface
   if (tc->num_maps) {
      if (ps == tc->surface)
         return;

      for (i = 0; i < tc->num_maps; i++) {
         tc->transfer_map[i] = NULL;
      }
      FREE(tc->transfer_map);
      tc->num_maps = 0;

      FREE(tc->clear_flags);
      tc->clear_flags_size = 0;
   }

   //设置图块的surface指向framebuffer的surface
   tc->surface = ps;

   if (ps) {

      tc->num_maps = ps->u.tex.last_layer - ps->u.tex.first_layer + 1;     //映射的数量
      tc->transfer = CALLOC(tc->num_maps, sizeof(struct pipe_transfer *));
      tc->transfer_map = CALLOC(tc->num_maps, sizeof(void *));

      tc->clear_flags_size = (MAX_WIDTH / TILE_SIZE) * (MAX_HEIGHT / TILE_SIZE) * tc->num_maps / 32 * sizeof(uint);     //清除标志大小
      tc->clear_flags = CALLOC(1, tc->clear_flags_size);    //清除标志

      if (ps->texture->target != PIPE_BUFFER) {
         for (i = 0; i < tc->num_maps; i++) {
            if(i >= PIPE_TEXTURE_BUFFER_NUMBER){
               break;
            }
		 	//2. 图块的transfer_map指向framebuffer
            if(1) {
            tc->transfer[i] = CALLOC(1, sizeof(struct pipe_transfer));
            struct pipe_transfer *transfer = tc->transfer[i];
            transfer->resource = tc->surface->texture;
            transfer->level = tc->surface->u.tex.level;
            transfer->usage = 1027;
            transfer->box.x = 0;
            transfer->box.y = 0;
            transfer->box.z = 0;
            transfer->box.width = tc->surface->width;
            transfer->box.height = tc->surface->height;
            transfer->box.depth = 1;
            transfer->stride = (tc->surface->width + 4) * 4;
            transfer->layer_stride = 0;
            }
            
            tc->transfer_map[i] = ps->texture->userdata + i * PIPE_TEXTURE_BUFFER_SIZE;
         }
      }
      else {
         /* can't render to buffers */
         assert(0);
      }

      //临时设置为非深度模板
      //tc->depth_stencil = util_format_is_depth_or_stencil(ps->format);
      tc->depth_stencil = false;
   }
}

/**
 * 获得tile cache的表面
 * 
 * \param tc -[in] tile cache
 * \return 返回tile cache的表面
 */
struct pipe_surface *
sp_tile_cache_get_surface(struct softpipe_tile_cache *tc)
{
   return tc->surface;
}

/**
 * 夹子瓷砖针对传输昏暗。
 *
 * XXX: 这只是剪辑宽度和高度！
 *
 * \return TRUE if tile is totally clipped, FALSE otherwise
 */
static inline bool
u_clip_tile(uint x, uint y, uint *w, uint *h, const struct pipe_box *box)
{
   if ((int) x >= box->width)
      return true;
   if ((int) y >= box->height)
      return true;
   if ((int) (x + *w) > box->width)
      *w = box->width - x;
   if ((int) (y + *h) > box->height)
      *h = box->height - y;
   return false;
}

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

static inline void
util_format_b8g8r8x8_unorm_pack_rgba_float(uint8_t *dst_row, unsigned dst_stride, const float *src_row, unsigned src_stride, unsigned width, unsigned height)
{
   unsigned x, y;
   for(y = 0; y < height; y += 1) {
      const float *src = src_row;
      uint8_t *dst = dst_row;
      for(x = 0; x < width; x += 1) {
#ifdef PIPE_ARCH_BIG_ENDIAN
         uint32_t value = 0;
         value |= (float_to_ubyte(src[2])) << 24;
         value |= ((float_to_ubyte(src[1])) & 0xff) << 16;
         value |= ((float_to_ubyte(src[0])) & 0xff) << 8;
         *(uint32_t *)dst = value;
#else
         uint32_t value = 0;
         value |= (float_to_ubyte(src[2])) & 0xff;
         value |= ((float_to_ubyte(src[1])) & 0xff) << 8;
         value |= ((float_to_ubyte(src[0])) & 0xff) << 16;
         *(uint32_t *)dst = value;
#endif
         src += 4;
         dst += 4;
      }
      dst_row += dst_stride;
      src_row += src_stride/sizeof(*src_row);
   }
}

/**
 * 按照制定格式拷贝数据
 * 
 * \param format -[in] 数据格式
 * \param src -[in] 源数据
 * \param src_stride -[in] 源数据步幅
 * \param dst -[in] 目标数据
 * \param dst_stride -[in] 目标数据步幅
 * \param x -[in] x坐标
 * \param y -[in] y坐标
 * \param w -[in] 宽度
 * \param h -[in] 高度
 */
void
util_format_write_4f(UNUSED enum pipe_format format,
                     const float *src, unsigned src_stride,
                     void *dst, unsigned dst_stride,
                     unsigned x, unsigned y, unsigned w, unsigned h)
{
   struct util_format_description *format_desc;
   uint8_t *dst_row;
   const float *src_row;

   //暂时注释掉
   #if 0
   format_desc = util_format_description(format);
   #endif
   format_desc = MALLOC_STRUCT(util_format_description);
   format_desc->format = PIPE_FORMAT_B8G8R8X8_UNORM;
   format_desc->name = "PIPE_FORMAT_B8G8R8X8_UNORM";
   format_desc->short_name = "b8g8r8x8_unorm";
   format_desc->block.width = 1;
   format_desc->block.height = 1;
   format_desc->block.bits = 32;
   format_desc->layout = UTIL_FORMAT_LAYOUT_PLAIN; 
   format_desc->nr_channels = 4;
   format_desc->is_array = 1;
   format_desc->is_bitmask = 1;
   format_desc->is_mixed = 0;
   for(int i = 0; i < 4; i++){
      format_desc->channel[i].type = 1;
      format_desc->channel[i].normalized = 1;
      format_desc->channel[i].pure_integer = 0;
      format_desc->channel[i].size = 8;
      format_desc->channel[i].shift = 0;
   }
   format_desc->swizzle[0] = 2;
   format_desc->swizzle[1] = 1;
   format_desc->swizzle[2] = 0;
   format_desc->swizzle[3] = 5;
   format_desc->colorspace = UTIL_FORMAT_COLORSPACE_RGB;


   assert(x % format_desc->block.width == 0);
   assert(y % format_desc->block.height == 0);

   dst_row = (uint8_t *)dst + y*dst_stride + x*(format_desc->block.bits/8);
   src_row = src;

   util_format_b8g8r8x8_unorm_pack_rgba_float(dst_row, dst_stride, src_row, src_stride, w, h);
}

/**
 * Copy 2D rect from one place to another.
 * Position and sizes are in pixels.
 * src_stride may be negative to do vertical flip of pixels from source.
 * /param 
 */
void
util_copy_rect(ubyte * dst,
               UNUSED enum pipe_format format,
               unsigned dst_stride,
               unsigned dst_x,
               unsigned dst_y,
               unsigned width,
               unsigned height,
               const ubyte * src,
               int src_stride,
               unsigned src_x,
               unsigned src_y)
{
   unsigned i;
   int src_stride_pos = src_stride < 0 ? -src_stride : src_stride;
   //暂时注释掉
   #if 0
   int blocksize = util_format_get_blocksize(format);		//单个像素字节数
   int blockwidth = util_format_get_blockwidth(format);
   int blockheight = util_format_get_blockheight(format);
   #endif
   int blocksize = 4;		//单个像素字节数
   int blockwidth = 1;
   int blockheight = 1;

   assert(blocksize > 0);
   assert(blockwidth > 0);
   assert(blockheight > 0);

   dst_x /= blockwidth;
   dst_y /= blockheight;
   width = (width + blockwidth - 1)/blockwidth;
   height = (height + blockheight - 1)/blockheight;
   src_x /= blockwidth;
   src_y /= blockheight;

   dst += dst_x * blocksize;
   src += src_x * blocksize;
   dst += dst_y * dst_stride;
   src += src_y * src_stride_pos;
   width *= blocksize;

   if (width == dst_stride && width == (unsigned)src_stride)
      memcpy(dst, src, height * width);
   else {
      for (i = 0; i < height; i++) {
         memcpy(dst, src, width);
         dst += dst_stride;
         src += src_stride;
      }
   }
}

/**
 * 将原始像素的原始块从用户存储器移到传输对象。
 */
void
pipe_put_tile_raw(struct pipe_transfer *pt,
                  void *dst,
                  uint x, uint y, uint w, uint h,
                  const void *src, int src_stride)
{
   enum pipe_format format = pt->resource->format;

   if (src_stride == 0){
      //暂时注释掉
      #if 0
      src_stride = util_format_get_stride(format, w);
      #endif
      src_stride = 256;
   }

   if (u_clip_tile(x, y, &w, &h, &pt->box))
      return;

   util_copy_rect(dst, format, pt->stride, x, y, w, h, src, src_stride, 0, 0);
}

/**
 * 把tile cache刷新到framebuffer
 * 
 * \param pt -[in] pipe transfer
 * \param dst -[in] 目标(framebuffer)
 * \param x -[in] x坐标(窗口)
 * \param y -[in] y坐标（窗口
 * \param w -[in] 宽度（tile）
 * \param h -[in] 高度（tile）
 * \param format -[in] 数据格式
 * \param p -[in] 源数据（tile）
 */
void
pipe_put_tile_rgba_format(struct pipe_transfer *pt,
                          void *dst,
                          uint x, uint y, uint w, uint h,
                          enum pipe_format format,
                          const float *p)
{
   unsigned src_stride = w * 4;
   void *packed;

   //夹子瓷砖针对传输昏暗
   if (u_clip_tile(x, y, &w, &h, &pt->box))
      return;

   //暂时注释掉
   //packed = MALLOC(util_format_get_nblocks(format, w, h) * util_format_get_blocksize(format));
   packed = MALLOC(4096 * 4);

   if (!packed)
      return;

   switch (format) {
   case PIPE_FORMAT_Z16_UNORM:
      /*z16_put_tile_rgba((ushort *) packed, w, h, p, src_stride);*/
      break;
   case PIPE_FORMAT_Z32_UNORM:
      /*z32_put_tile_rgba((unsigned *) packed, w, h, p, src_stride);*/
      break;
   case PIPE_FORMAT_Z24_UNORM_S8_UINT:
   case PIPE_FORMAT_Z24X8_UNORM:
      /*s8z24_put_tile_rgba((unsigned *) packed, w, h, p, src_stride);*/
      break;
   case PIPE_FORMAT_S8_UINT_Z24_UNORM:
   case PIPE_FORMAT_X8Z24_UNORM:
      /*z24s8_put_tile_rgba((unsigned *) packed, w, h, p, src_stride);*/
      break;
   case PIPE_FORMAT_Z32_FLOAT:
      /*z32f_put_tile_rgba((unsigned *) packed, w, h, p, src_stride);*/
      break;
   case PIPE_FORMAT_Z32_FLOAT_S8X24_UINT:
      /*z32f_s8x24_put_tile_rgba((unsigned *) packed, w, h, p, src_stride);*/
      break;
   default:
      //暂时注释掉
      /*
      util_format_write_4f(format,
                           p, src_stride * sizeof(float),
                           packed, util_format_get_stride(format, w),
                           0, 0, w, h);
      */
     util_format_write_4f(format,
                     p, src_stride * sizeof(float),
                     packed, 256,
                     0, 0, w, h);
   }

   pipe_put_tile_raw(pt, dst, x, y, w, h, packed, 0);

   FREE(packed);
}


/**
 * Set pixels in a tile to the given clear color/value, float.
 * 
 * \param tile -[in] 用于清除的tile
 * \param format -[in] 数据格式
 * \param clear_value -[in] 清除颜色值
 */
static void
clear_tile_rgba(struct softpipe_cached_tile *tile,
                UNUSED enum pipe_format format,
                const union pipe_color_union *clear_value)
{
   if (clear_value->f[0] == 0.0 &&
       clear_value->f[1] == 0.0 &&
       clear_value->f[2] == 0.0 &&
       clear_value->f[3] == 0.0) {
      memset(tile->data.color, 0, sizeof(tile->data.color));
   }
   else {
      uint i, j;

      //暂时注释掉
      #if 0
      if (util_format_is_pure_uint(format)) {
         for (i = 0; i < TILE_SIZE; i++) {
            for (j = 0; j < TILE_SIZE; j++) {
               tile->data.colorui128[i][j][0] = clear_value->ui[0];
               tile->data.colorui128[i][j][1] = clear_value->ui[1];
               tile->data.colorui128[i][j][2] = clear_value->ui[2];
               tile->data.colorui128[i][j][3] = clear_value->ui[3];
            }
         }
      } else if (util_format_is_pure_sint(format)) {
         for (i = 0; i < TILE_SIZE; i++) {
            for (j = 0; j < TILE_SIZE; j++) {
               tile->data.colori128[i][j][0] = clear_value->i[0];
               tile->data.colori128[i][j][1] = clear_value->i[1];
               tile->data.colori128[i][j][2] = clear_value->i[2];
               tile->data.colori128[i][j][3] = clear_value->i[3];
            }
         }
      } else {
      #endif
         for (i = 0; i < TILE_SIZE; i++) {
            for (j = 0; j < TILE_SIZE; j++) {
               tile->data.color[i][j][0] = clear_value->f[0];
               tile->data.color[i][j][1] = clear_value->f[1];
               tile->data.color[i][j][2] = clear_value->f[2];
               tile->data.color[i][j][3] = clear_value->f[3];
            }
         }
      //}
   }
}

/**
 * Actually clear the tiles which were flagged as being in a clear state.
 * 
 * \param tc -[in] tile cache
 * \param layer -[in] 层索引
 */
static void
sp_tile_cache_flush_clear(struct softpipe_tile_cache *tc, int layer)
{
   struct pipe_transfer *pt = tc->transfer[layer];
   const uint w = tc->transfer[layer]->box.width;
   const uint h = tc->transfer[layer]->box.height;
   uint x, y;
   uint numCleared = 0;

   assert(pt->resource);

   /* clear the scratch tile to the clear value */
   if (tc->depth_stencil) {
      //暂时注释掉
      #if 0
      clear_tile(tc->tile, pt->resource->format, tc->clear_val);
      #endif
   } else {
      clear_tile_rgba(tc->tile, pt->resource->format, &tc->clear_color);
   }

   /* push the tile to all positions marked as clear */
   for (y = 0; y < h; y += TILE_SIZE) {
      for (x = 0; x < w; x += TILE_SIZE) {
         union tile_address addr = tile_address(x, y, layer);

         if (is_clear_flag_set(tc->clear_flags, addr, tc->clear_flags_size)) {
            /* write the scratch tile to the surface */
            if (tc->depth_stencil) {
               //暂时注释掉
               #if 0
               pipe_put_tile_raw(pt, tc->transfer_map[layer],
                                 x, y, TILE_SIZE, TILE_SIZE,
                                 tc->tile->data.any, 0/*STRIDE*/);
               #endif
            }
            else {
               //暂时注释掉
               #if 0
               if (util_format_is_pure_uint(tc->surface->format)) {
                  pipe_put_tile_ui_format(pt, tc->transfer_map[layer],
                                          x, y, TILE_SIZE, TILE_SIZE,
                                          tc->surface->format,
                                          (unsigned *) tc->tile->data.colorui128);
               } else if (util_format_is_pure_sint(tc->surface->format)) {
                  pipe_put_tile_i_format(pt, tc->transfer_map[layer],
                                         x, y, TILE_SIZE, TILE_SIZE,
                                         tc->surface->format,
                                         (int *) tc->tile->data.colori128);
               } else {
               #endif
                  pipe_put_tile_rgba_format(pt, tc->transfer_map[layer],
                                            x, y, TILE_SIZE, TILE_SIZE,
                                            tc->surface->format,
                                            (float *) tc->tile->data.color);
               //}
            }
            numCleared++;
         }
      }
   }


#if 0
   debug_printf("num cleared: %u\n", numCleared);
#endif
}

/**
 * 刷新tile数据
 * \param tc -[in] tile cache
 * \param pos -[in] tile位置
 */
static void
sp_flush_tile(struct softpipe_tile_cache* tc, unsigned pos)
{
   int layer = tc->tile_addrs[pos].bits.layer;
   if (!tc->tile_addrs[pos].bits.invalid) {     //当前tile有效
      if (tc->depth_stencil) {      //深度模板
         //暂时注释掉
         #if(0)
         pipe_put_tile_raw(tc->transfer[layer], tc->transfer_map[layer],
                           tc->tile_addrs[pos].bits.x * TILE_SIZE,
                           tc->tile_addrs[pos].bits.y * TILE_SIZE,
                           TILE_SIZE, TILE_SIZE,
                           tc->entries[pos]->data.depth32, 0/*STRIDE*/);
         #
         #endif
      }
      else {
         //暂时注释掉
         #if 0
         if (util_format_is_pure_uint(tc->surface->format)) {
            pipe_put_tile_ui_format(tc->transfer[layer], tc->transfer_map[layer],
                                    tc->tile_addrs[pos].bits.x * TILE_SIZE,
                                    tc->tile_addrs[pos].bits.y * TILE_SIZE,
                                    TILE_SIZE, TILE_SIZE,
                                    tc->surface->format,
                                    (unsigned *) tc->entries[pos]->data.colorui128);
         } else if (util_format_is_pure_sint(tc->surface->format)) {
            pipe_put_tile_i_format(tc->transfer[layer], tc->transfer_map[layer],
                                   tc->tile_addrs[pos].bits.x * TILE_SIZE,
                                   tc->tile_addrs[pos].bits.y * TILE_SIZE,
                                   TILE_SIZE, TILE_SIZE,
                                   tc->surface->format,
                                   (int *) tc->entries[pos]->data.colori128);
         } else
         #endif 
         {
            pipe_put_tile_rgba_format(tc->transfer[layer], tc->transfer_map[layer],
                                      tc->tile_addrs[pos].bits.x * TILE_SIZE,
                                      tc->tile_addrs[pos].bits.y * TILE_SIZE,
                                      TILE_SIZE, TILE_SIZE,
                                      tc->surface->format,
                                      (float *) tc->entries[pos]->data.color);
         }
      }
      tc->tile_addrs[pos].bits.invalid = 1;  /* mark as empty */
   }
}

/**
 * 分配一个tile数据
 * 
 * \param tc -[in] tile cache
 * \return 返回一个tile
 */
static struct softpipe_cached_tile *
sp_alloc_tile(struct softpipe_tile_cache *tc)
{
   struct softpipe_cached_tile * tile = MALLOC_STRUCT(softpipe_cached_tile);    //tile数据
   //分配失败
   if (!tile)
   {
      //在这种情况下，窃取现有的瓷砖
      if (!tc->tile)        //用于清除的tile不是空
      {
         unsigned pos;
         for (pos = 0; pos < ARRAY_SIZE(tc->entries); ++pos) {
            if (!tc->entries[pos])      //当前位置没数据
               continue;

            sp_flush_tile(tc, pos);       //把当前tile刷新到framebuffer里
            tc->tile = tc->entries[pos];        //用于清除的tile为当前tile
            tc->entries[pos] = NULL;
            break;
         }

         /* this should never happen */
         if (!tc->tile)
            abort();
      }

      tile = tc->tile;      //返回用于清除的tile
      tc->tile = NULL;

      tc->last_tile_addr.bits.invalid = 1;      //最近tile地址无效
   }
   return tile;
}

/**
* 刷新 tile cache：将所有脏块写到transfer
*
* \param tc -[in] tile cache用于暂存图块数据
*/
void
sp_flush_tile_cache(struct softpipe_tile_cache *tc)
{
   int inuse = 0, pos;
   int i;
   if (tc->num_maps) {
      //1. 把数据写到transfer
      for (pos = 0; (unsigned)pos < ARRAY_SIZE(tc->entries); pos++) {
         struct softpipe_cached_tile *tile = tc->entries[pos];
         if (!tile)
         {
            assert(tc->tile_addrs[pos].bits.invalid);
            continue;
         }
         sp_flush_tile(tc, pos);
         ++inuse;
      }

      //2. 分配用于清除的tile
      if (!tc->tile)
         tc->tile = sp_alloc_tile(tc);

      //3. 清除
      for (i = 0; i < tc->num_maps; i++)
         sp_tile_cache_flush_clear(tc, i);
      /* reset all clear flags to zero */
      memset(tc->clear_flags, 0, tc->clear_flags_size);

      tc->last_tile_addr.bits.invalid = 1;
   }

#if 0
   debug_printf("flushed tiles in use: %d\n", inuse);
#endif
}

/**
 * 将像素的原始块从传输对象移至用户内存。
 * 
 * \param pt -[in] pipe transfer
 * \param src -[in] 源数据(framebuffer)
 * \param x -[in] x坐标
 * \param y -[in] y坐标
 * \param w -[in] 宽度
 * \param h -[in] 高度
 * \param dst -[in] 目标数据
 * \param dst_stride -[in] 目标数据步幅
 */
void
pipe_get_tile_raw(struct pipe_transfer *pt,
                  const void *src,
                  uint x, uint y, uint w, uint h,
                  void *dst, int dst_stride)
{
   if (dst_stride == 0){
      //暂时打个桩
      //dst_stride = util_format_get_stride(pt->resource->format, w);
      dst_stride = w * 4;
   }
   
   
   if (u_clip_tile(x, y, &w, &h, &pt->box))
      return;

   util_copy_rect(dst, pt->resource->format, dst_stride, 0, 0, w, h, src, pt->stride, x, y);
}

/**
 * Convert ubyte to float in [0, 1].
 */
static inline float
ubyte_to_float(ubyte ub)
{
   return (float) ub * (1.0f / 255.0f);
}

static inline void
util_format_b8g8r8x8_unorm_unpack_rgba_float(float *dst_row, unsigned dst_stride, const uint8_t *src_row, unsigned src_stride, unsigned width, unsigned height)
{
   unsigned x, y;
   for(y = 0; y < height; y += 1) {
      float *dst = dst_row;
      const uint8_t *src = src_row;
      for(x = 0; x < width; x += 1) {
#ifdef PIPE_ARCH_BIG_ENDIAN
         uint32_t value = *(const uint32_t *)src;
         uint32_t b;
         uint32_t g;
         uint32_t r;
         b = value >> 24;
         g = (value >> 16) & 0xff;
         r = (value >> 8) & 0xff;
         dst[0] = ubyte_to_float(r); /* r */
         dst[1] = ubyte_to_float(g); /* g */
         dst[2] = ubyte_to_float(b); /* b */
         dst[3] = 1; /* a */
#else
         uint32_t value = *(const uint32_t *)src;
         uint32_t b;
         uint32_t g;
         uint32_t r;
         b = (value) & 0xff;
         g = (value >> 8) & 0xff;
         r = (value >> 16) & 0xff;
         dst[0] = ubyte_to_float(r); /* r */
         dst[1] = ubyte_to_float(g); /* g */
         dst[2] = ubyte_to_float(b); /* b */
         dst[3] = 1; /* a */
#endif
         src += 4;
         dst += 4;
      }
      src_row += src_stride;
      dst_row += dst_stride/sizeof(*dst_row);
   }
}

void
util_format_read_4f(UNUSED enum pipe_format format,
                    float *dst, unsigned dst_stride,
                    const void *src, unsigned src_stride,
                    unsigned x, unsigned y, unsigned w, unsigned h)
{
   struct util_format_description *format_desc;
   const uint8_t *src_row;
   float *dst_row;

   //暂时注释掉
   #if 0
   format_desc = util_format_description(format);
   #endif
   format_desc = MALLOC_STRUCT(util_format_description);
   format_desc->format = PIPE_FORMAT_B8G8R8X8_UNORM;
   format_desc->name = "PIPE_FORMAT_B8G8R8X8_UNORM";
   format_desc->short_name = "b8g8r8x8_unorm";
   format_desc->block.width = 1;
   format_desc->block.height = 1;
   format_desc->block.bits = 32;
   format_desc->layout = UTIL_FORMAT_LAYOUT_PLAIN; 
   format_desc->nr_channels = 4;
   format_desc->is_array = 1;
   format_desc->is_bitmask = 1;
   format_desc->is_mixed = 0;
   format_desc->swizzle[0] = 2;
   format_desc->swizzle[1] = 1;
   format_desc->swizzle[2] = 0;
   format_desc->swizzle[3] = 5;

   assert(x % format_desc->block.width == 0);
   assert(y % format_desc->block.height == 0);

   src_row = (const uint8_t *)src + y*src_stride + x*(format_desc->block.bits/8);
   dst_row = dst;

   util_format_b8g8r8x8_unorm_unpack_rgba_float(dst_row, dst_stride, src_row, src_stride, w, h);
}


void
pipe_tile_raw_to_rgba(enum pipe_format format,
                      const void *src,
                      uint w, uint h,
                      float *dst, unsigned dst_stride)
{
   switch (format) {
   case PIPE_FORMAT_Z16_UNORM:
      //z16_get_tile_rgba((ushort *) src, w, h, dst, dst_stride);
      break;
   case PIPE_FORMAT_Z32_UNORM:
      //z32_get_tile_rgba((unsigned *) src, w, h, dst, dst_stride);
      break;
   case PIPE_FORMAT_Z24_UNORM_S8_UINT:
   case PIPE_FORMAT_Z24X8_UNORM:
      //s8z24_get_tile_rgba((unsigned *) src, w, h, dst, dst_stride);
      break;
   case PIPE_FORMAT_S8_UINT:
      //s8_get_tile_rgba((unsigned char *) src, w, h, dst, dst_stride);
      break;
   case PIPE_FORMAT_X24S8_UINT:
      //s8x24_get_tile_rgba((unsigned *) src, w, h, dst, dst_stride);
      break;
   case PIPE_FORMAT_S8_UINT_Z24_UNORM:
   case PIPE_FORMAT_X8Z24_UNORM:
      //z24s8_get_tile_rgba((unsigned *) src, w, h, dst, dst_stride);
      break;
   case PIPE_FORMAT_S8X24_UINT:
      //x24s8_get_tile_rgba((unsigned *) src, w, h, dst, dst_stride);
      break;
   case PIPE_FORMAT_Z32_FLOAT:
      //z32f_get_tile_rgba((float *) src, w, h, dst, dst_stride);
      break;
   case PIPE_FORMAT_Z32_FLOAT_S8X24_UINT:
      //z32f_x24s8_get_tile_rgba((float *) src, w, h, dst, dst_stride);
      break;
   case PIPE_FORMAT_X32_S8X24_UINT:
      //x32_s8_get_tile_rgba((unsigned *) src, w, h, dst, dst_stride);
      break;
   default:
      /*
      util_format_read_4f(format,
                          dst, dst_stride * sizeof(float),
                          src, util_format_get_stride(format, w),
                          0, 0, w, h);
      */
     util_format_read_4f(format,
                    dst, dst_stride * sizeof(float),
                    src, 4*w,
                    0, 0, w, h);
   }
}

/**
 * 获取一个tile
 * 
 * \param pt -[in] pipe transfer
 * \param src -[in] 源数据(framebuffer)
 * \param x -[in] x坐标
 * \param y -[in] y坐标
 * \param w -[in] 宽度
 * \param h -[in] 高度
 * \param format -[in] 数据格式
 * \param p -[in] 目标数据(tile)
 */
void
pipe_get_tile_rgba_format(struct pipe_transfer *pt,
                          const void *src,
                          uint x, uint y, uint w, uint h,
                          enum pipe_format format,
                          float *p)
{
   unsigned dst_stride = w * 4;
   void *packed;

   if (u_clip_tile(x, y, &w, &h, &pt->box)) {
      return;
   }

   //暂时注释掉
   //packed = MALLOC(util_format_get_nblocks(format, w, h) * util_format_get_blocksize(format));
   packed = MALLOC(w * h * 4);
   if (!packed) {
      return;
   }

   if (format == PIPE_FORMAT_UYVY || format == PIPE_FORMAT_YUYV) {
      assert((x & 1) == 0);
   }

   pipe_get_tile_raw(pt, src, x, y, w, h, packed, 0);

   pipe_tile_raw_to_rgba(format, packed, w, h, p, dst_stride);

   FREE(packed);
}


void
pipe_get_tile_rgba(struct pipe_transfer *pt,
                   const void *src,
                   uint x, uint y, uint w, uint h,
                   float *p)
{
   pipe_get_tile_rgba_format(pt, src, x, y, w, h, pt->resource->format, p);
}

/**
 * 从缓存中获取一个图块。
 * 如果缓存中没有空闲的图块，则返回一个新的图块。
 * 
 * \param tc -[in] tile cache
 * \param addr -[in] tile地址
 * \return 返回一个图块
 * 
 */
struct softpipe_cached_tile *
sp_find_cached_tile(struct softpipe_tile_cache *tc,
                    union tile_address addr )
{
   struct pipe_transfer *pt;
   
   //获得tile地址
   const int pos = CACHE_POS(addr.bits.x,
                             addr.bits.y, addr.bits.layer);      //每行5个，每列10个
   struct softpipe_cached_tile *tile = tc->entries[pos];        //tile
   int layer;
   //如果当前位置没有tile数据，则分配一个tile
   if (!tile) {
      tile = sp_alloc_tile(tc);
      tc->entries[pos] = tile;
   }

   //如果当前tile地址和最近tile地址不一样，则把最近tile刷新到framebuffer里
   if (addr.value != tc->tile_addrs[pos].value) {

      layer = tc->tile_addrs[pos].bits.layer;         //层索引
      if (tc->tile_addrs[pos].bits.invalid == 0) {
         /* put dirty tile back in framebuffer */
         if (tc->depth_stencil) {
            //暂时注释掉
            #if 0
            pipe_put_tile_raw(tc->transfer[layer], tc->transfer_map[layer],
                              tc->tile_addrs[pos].bits.x * TILE_SIZE,
                              tc->tile_addrs[pos].bits.y * TILE_SIZE,
                              TILE_SIZE, TILE_SIZE,
                              tile->data.depth32, 0/*STRIDE*/);
            #endif
         }
         else {
            #if 0
            if (util_format_is_pure_uint(tc->surface->format)) {
               pipe_put_tile_ui_format(tc->transfer[layer], tc->transfer_map[layer],
                                      tc->tile_addrs[pos].bits.x * TILE_SIZE,
                                      tc->tile_addrs[pos].bits.y * TILE_SIZE,
                                      TILE_SIZE, TILE_SIZE,
                                      tc->surface->format,
                                      (unsigned *) tile->data.colorui128);
            } else if (util_format_is_pure_sint(tc->surface->format)) {
               pipe_put_tile_i_format(tc->transfer[layer], tc->transfer_map[layer],
                                      tc->tile_addrs[pos].bits.x * TILE_SIZE,
                                      tc->tile_addrs[pos].bits.y * TILE_SIZE,
                                      TILE_SIZE, TILE_SIZE,
                                      tc->surface->format,
                                      (int *) tile->data.colori128);
            } else
            #endif 
            {
               pipe_put_tile_rgba_format(tc->transfer[layer], tc->transfer_map[layer],
                                         tc->tile_addrs[pos].bits.x * TILE_SIZE,
                                         tc->tile_addrs[pos].bits.y * TILE_SIZE,
                                         TILE_SIZE, TILE_SIZE,
                                         tc->surface->format,
                                         (float *) tile->data.color);
            }
         }
      }

      tc->tile_addrs[pos] = addr;      //

      layer = tc->tile_addrs[pos].bits.layer;
      pt = tc->transfer[layer];
      assert(pt->resource);

      //
      if (is_clear_flag_set(tc->clear_flags, addr, tc->clear_flags_size)) {
         /* don't get tile from framebuffer, just clear it */
         if (tc->depth_stencil) {
            #if 0
            clear_tile(tile, pt->resource->format, tc->clear_val);
            #endif
         }
         else {
            clear_tile_rgba(tile, pt->resource->format, &tc->clear_color);
         }
         clear_clear_flag(tc->clear_flags, addr, tc->clear_flags_size);
      }
      else {
         /* get new tile data from transfer */
         if (tc->depth_stencil) {
            #if 0
            pipe_get_tile_raw(tc->transfer[layer], tc->transfer_map[layer],
                              tc->tile_addrs[pos].bits.x * TILE_SIZE,
                              tc->tile_addrs[pos].bits.y * TILE_SIZE,
                              TILE_SIZE, TILE_SIZE,
                              tile->data.depth32, 0/*STRIDE*/);
            #endif
         }
         else {
            #if 0
            if (util_format_is_pure_uint(tc->surface->format)) {
               pipe_get_tile_ui_format(tc->transfer[layer], tc->transfer_map[layer],
                                         tc->tile_addrs[pos].bits.x * TILE_SIZE,
                                         tc->tile_addrs[pos].bits.y * TILE_SIZE,
                                         TILE_SIZE, TILE_SIZE,
                                         tc->surface->format,
                                         (unsigned *) tile->data.colorui128);
            } else if (util_format_is_pure_sint(tc->surface->format)) {
               pipe_get_tile_i_format(tc->transfer[layer], tc->transfer_map[layer],
                                         tc->tile_addrs[pos].bits.x * TILE_SIZE,
                                         tc->tile_addrs[pos].bits.y * TILE_SIZE,
                                         TILE_SIZE, TILE_SIZE,
                                         tc->surface->format,
                                         (int *) tile->data.colori128);
            } else
            #endif 
            {
               pipe_get_tile_rgba_format(tc->transfer[layer], tc->transfer_map[layer],
                                         tc->tile_addrs[pos].bits.x * TILE_SIZE,
                                         tc->tile_addrs[pos].bits.y * TILE_SIZE,
                                         TILE_SIZE, TILE_SIZE,
                                         tc->surface->format,
                                         (float *) tile->data.color);
            }
         }
      }
   }

   tc->last_tile = tile;
   tc->last_tile_addr = addr;
   return tile;
}

/**
 * 当整个表面被清除为某个值时，我们可以避免
 * 获取上面的图块。
 * 保存颜色并为屏幕的每个图块设置一个“清除标志”。
 * 
 * \param tc -[in] tile cache
 * \param color -[in] 清除颜色(用于清除color buffer)
 * \param clearValue -[in] 清除值（用于清除深度模板buffer）
 */
void
sp_tile_cache_clear(struct softpipe_tile_cache *tc,
                    const union pipe_color_union *color,
                    uint64_t clearValue)
{
   uint pos;

   tc->clear_color = *color;

   tc->clear_val = clearValue;

   /* 设置标志以指示所有瓷砖已清除 */
      memset(tc->clear_flags, 255, tc->clear_flags_size);

   for (pos = 0; pos < ARRAY_SIZE(tc->tile_addrs); pos++) {
      tc->tile_addrs[pos].bits.invalid = 1;
   }
   tc->last_tile_addr.bits.invalid = 1;
}
