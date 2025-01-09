#ifndef __TILE_CACHE_H__
#define __TILE_CACHE_H__

#include <stdint.h>
#include <sys/types.h>
#include <stdbool.h>
#include "context.h"
#include "state.h"

typedef unsigned int       uint;
typedef unsigned short     ushort;
typedef unsigned char      ubyte;

struct softpipe_tile_cache;

#define MAYBE_UNUSED __attribute__((unused))

/**
 * Tile size (width and height). This needs to be a power of two.
 */
#define TILE_SIZE_LOG2 6
#define TILE_SIZE (1 << TILE_SIZE_LOG2)      //64
#define SP_MAX_TEXTURE_2D_LEVELS 15  /* 16K x 16K */
#define TILE_ADDR_BITS (SP_MAX_TEXTURE_2D_LEVELS - 1 - TILE_SIZE_LOG2)     //15-1-6=8
/** Max surface size */
#define MAX_WIDTH (1 << (SP_MAX_TEXTURE_2D_LEVELS - 1))
#define MAX_HEIGHT (1 << (SP_MAX_TEXTURE_2D_LEVELS - 1))

/*
 * 静态（编译时）断言。
 * 基本上，使用COND对数组进行维数。如果COND为假/零，数组大小将为-1，我们将得到一个编译错误
 */
#define STATIC_ASSERT(COND) \
   do { \
      (void) sizeof(char [1 - 2*!(COND)]); \
   } while (0)

#define CALLOC_STRUCT(T)   (struct T *) calloc(1, sizeof(struct T))
#define FREE(x) free(x)

/**
 * tile地址
 */
union tile_address {
   struct {
      unsigned x:TILE_ADDR_BITS;     /* 16K / TILE_SIZE */
      unsigned y:TILE_ADDR_BITS;     /* 16K / TILE_SIZE */
      unsigned invalid:1;           //是否无效
      unsigned layer:8;             //层索引
      unsigned pad:7;               //填充
   } bits;
   unsigned value;
};

/**
 * 一个cached tile数据
 */
struct softpipe_cached_tile
{
   union {
      float color[TILE_SIZE][TILE_SIZE][4];
      uint color32[TILE_SIZE][TILE_SIZE];
      uint depth32[TILE_SIZE][TILE_SIZE];
      ushort depth16[TILE_SIZE][TILE_SIZE];
      ubyte stencil8[TILE_SIZE][TILE_SIZE];
      uint colorui128[TILE_SIZE][TILE_SIZE][4];
      int colori128[TILE_SIZE][TILE_SIZE][4];
      uint64_t depth64[TILE_SIZE][TILE_SIZE];
      ubyte any[1];
   } data;
};

/**
 * 指向每个r,g,b,a的fiu数组的并集
 */
union pipe_color_union
{
   float f[4];
   int i[4];
   unsigned int ui[4];
};

#define NUM_ENTRIES 50     //tile缓存的大小

/**
 *tile cache，包括tile地址和tile数据
 */
struct softpipe_tile_cache
{
   struct softpipe_context *pipe;   //上下文
   struct pipe_surface *surface;    //我们正在缓存的表面
   struct pipe_transfer **transfer;   //指向framebuffer的transfer
   void **transfer_map;   //指向framebuffer
   int num_maps;

   union tile_address tile_addrs[NUM_ENTRIES];        //tile地址
   struct softpipe_cached_tile *entries[NUM_ENTRIES];    //tile数据
   uint *clear_flags;      //
   uint clear_flags_size;
   union pipe_color_union clear_color;    //清除颜色(color buffer)
   uint64_t clear_val;        //清除值(depth buffer)
   bool depth_stencil;        //Is the surface a depth/stencil format? 

   struct softpipe_cached_tile *tile;  //用于清除的tile

   union tile_address last_tile_addr;  //最近一次tile地址
   struct softpipe_cached_tile *last_tile;  //最近检索到的tile
};

extern struct softpipe_tile_cache *
sp_create_tile_cache( struct softpipe_context *pipe );         //创建tile cache

extern void
sp_destroy_tile_cache(struct softpipe_tile_cache *tc);         //销毁tile cache

extern void
sp_tile_cache_set_surface(struct softpipe_tile_cache *tc,      //设置tile cache的surface
                          struct pipe_surface *sps);

extern struct pipe_surface *
sp_tile_cache_get_surface(struct softpipe_tile_cache *tc);     //获得tile cache的surface

extern void
sp_flush_tile_cache(struct softpipe_tile_cache *tc);        //刷新tile cache

extern void
sp_tile_cache_clear(struct softpipe_tile_cache *tc,         //清除tile cache
                    const union pipe_color_union *color,
                    uint64_t clearValue);

extern struct softpipe_cached_tile *
sp_find_cached_tile(struct softpipe_tile_cache *tc,         //查找cached tile
                    union tile_address addr );

extern void
pipe_get_tile_rgba(struct pipe_transfer *pt,
                   const void *src,
                   uint x, uint y, uint w, uint h,
                   float *p);
/**
 * 根据给出的窗口坐标获得tile地址
 * 
 * \param x -[in]    窗口坐标x
 * \param y -[in]    窗口坐标y
 * \param layer -[in]  层索引
 */
static inline union tile_address
tile_address( unsigned x,
              unsigned y, unsigned layer )
{
   union tile_address addr;

   addr.value = 0;
   addr.bits.x = x / TILE_SIZE;
   addr.bits.y = y / TILE_SIZE;
   addr.bits.layer = layer;
   return addr;
}

/**
 *  根据给出的窗口坐标获得tile，如果与上次查找匹配，则快速检索图块。
 * 
 * \param tc -[in]   tile cache
 * \param x -[in]    窗口坐标x
 * \param y -[in]    窗口坐标y
 * \param layer -[in]  层索引
 * \return  返回tile
 */
static inline struct softpipe_cached_tile *
sp_get_cached_tile(struct softpipe_tile_cache *tc,
                   int x, int y, int layer )
{
   union tile_address addr = tile_address( x, y, layer );   //tile地址

   if (tc->last_tile_addr.value == addr.value)
      return tc->last_tile;

   return sp_find_cached_tile( tc, addr );
}

/**
 * 返回包含 win pos (x,y) 的图块在缓存中的位置。
 * 我们目前使用直接映射缓存，因此这就像一个黑客密钥。
 * 在某些时候，我们应该研究一些更复杂的东西，比如
 * LRU 替换策略。
 */
#define CACHE_POS(x, y, l)                        \
   (((x) + (y) * 5 + (l) * 10) % NUM_ENTRIES)

#endif