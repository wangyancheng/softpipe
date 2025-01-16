#ifndef __SETUP_H__
#define __SETUP_H__

#include <stdbool.h>
#include "context.h"
#include "quad.h"
#include "state.h"

#define DEBUG_VERTS 1
#define DEBUG_FRAGS 1


/**
 * 每批处理的最大四边形数量（2x2 像素块）。
 * 由于我们依赖一些 32 位位掩码（每个四边形两位），因此不能任意增加该数量。
 */
#define MAX_QUADS 16

 /** 
 * 基元类型
 */
enum pipe_prim_type {
   PIPE_PRIM_POINTS,			//对应 GL_POINTS
   PIPE_PRIM_LINES,				//对应 GL_LINES
   PIPE_PRIM_LINE_LOOP,			//对应 GL_LINE_STRIP
   PIPE_PRIM_LINE_STRIP,
   PIPE_PRIM_TRIANGLES,
   PIPE_PRIM_TRIANGLE_STRIP,
   PIPE_PRIM_TRIANGLE_FAN,
   PIPE_PRIM_QUADS,
   PIPE_PRIM_QUAD_STRIP,
   PIPE_PRIM_POLYGON,
   PIPE_PRIM_LINES_ADJACENCY,
   PIPE_PRIM_LINE_STRIP_ADJACENCY,
   PIPE_PRIM_TRIANGLES_ADJACENCY,
   PIPE_PRIM_TRIANGLE_STRIP_ADJACENCY,
   PIPE_PRIM_PATCHES,
   PIPE_PRIM_MAX,
};

/**
* 属性插值模式
* 参考：https://www.khronos.org/opengl/wiki/Type_Qualifier_(GLSL)
* flat: 该值不会被插值。赋予片段着色器的值是来自该图元的激发顶点的值。
* noperspective: 该值将在窗口空间中线性插入。这通常不是您想要的，但它可以有其用途。
* smooth: 该值将以透视校正方式插入。如果不存在限定符，则这是默认值。
*/
enum sp_interp_mode {
   SP_INTERP_POS,       //用于点精灵
   SP_INTERP_CONSTANT,  //常量插值
   SP_INTERP_LINEAR,    //线性插值
   SP_INTERP_PERSPECTIVE   //透视插值
};

/**
* 插值系数
*/
struct tgsi_interp_coef
{
   float a0[TGSI_NUM_CHANNELS];	   //三角形在起始点的深度值类的值
   float dadx[TGSI_NUM_CHANNELS];   //表示当x坐标增加1个单位时，深度值的变化量
   float dady[TGSI_NUM_CHANNELS];   //表示当y坐标增加1个单位时，深度值的变化量
};

/**
 * Triangle edge info
 */
struct edge {
   float dx;		/**< X(v1) - X(v0), used only during setup */
   float dy;		/**< Y(v1) - Y(v0), used only during setup */
   float dxdy;		/**< dx/dy */
   float sx, sy;	/**< first sample point coord */
   int lines;		/**< number of lines on this edge */
};

/**
 * setup上下文
 */
struct setup_context {
   struct softpipe_context *softpipe;  //管线上下文

   struct quad_header quad[MAX_QUADS]; //四边形
   struct quad_header *quad_ptrs[MAX_QUADS];
   
   /* 顶点只是构成每个属性的一系列浮子
    * 转动。目前固定在4个浮子上，但应随时间变化。
    * Codegen将有助于应对。
   */
   const float (*vmax)[4];          //顶点最大值
   const float (*vmid)[4];          //顶点中间值
   const float (*vmin)[4];          //顶点最小值
   const float (*vprovoke)[4];		//provoking vertex（激发顶点）

   struct edge ebot;                //底边
   struct edge etop;                //顶边
   struct edge emaj;                //主要边
   
   float oneoverarea;               //面积的倒数
   int facing;                      //面 Front (0) , back (1)

   float pixel_offset;              //像素偏移
   unsigned max_layer;              //最大层

   struct tgsi_interp_coef coef[PIPE_MAX_SHADER_INPUTS];
   struct tgsi_interp_coef posCoef;  //插值系数/* For Z, W */

      struct {
      int left[2];   /**< [0] = row0, [1] = row1 */
      int right[2];
      int y;
   } span;        //跨度


   #if DEBUG_FRAGS
   unsigned int numFragsEmitted;  /**< per primitive */
   unsigned int numFragsWritten;  /**< per primitive */
   #endif

   unsigned cull_face;		/* which faces cull */
   unsigned nr_vertex_attrs;  //顶点属性数量
};

union fi {
   float f;
   int32_t i;
   uint32_t ui;
};

/**
 * 是否是无穷大或NaN
 */
static inline bool
util_is_inf_or_nan(float x)
{
   union fi tmp;
   tmp.f = x;
   return (tmp.ui & 0x7f800000) == 0x7f800000;
}

static inline unsigned
sp_clamp_viewport_idx(int idx)
{
   return (PIPE_MAX_VIEWPORTS > idx && idx >= 0) ? idx : 0;
}

/**
 * Return number of bits set in n.
 */
static inline unsigned
util_bitcount(unsigned n)
{
#if defined(HAVE___BUILTIN_POPCOUNT)
   return __builtin_popcount(n);
#else
   /* K&R classic bitcount.
    *
    * For each iteration, clear the LSB from the bitfield.
    * Requires only one iteration per set bit, instead of
    * one iteration per bit less than highest set bit.
    */
   unsigned bits;
   for (bits = 0; n; bits++) {
      n &= n - 1;
   }
   return bits;
#endif
}

#define MIN2( A, B )   ( (A)<(B) ? (A) : (B) )
#define MAX2( A, B )   ( (A)>(B) ? (A) : (B) )

extern void setup_test(void);

#endif