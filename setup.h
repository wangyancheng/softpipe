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
   float a0[TGSI_NUM_CHANNELS];	   //provoking vertex, in an xyzw layout 
   float dadx[TGSI_NUM_CHANNELS];   //x坐标偏移
   float dady[TGSI_NUM_CHANNELS];   //y坐标偏移
};

/**
 * setup上下文
 */
struct setup_context {
   int facing;                      //面 Front (0) , back (1)
   struct softpipe_context *softpipe;  //管线上下文
   struct quad_header quad[MAX_QUADS]; //四边形
   const float (*vprovoke)[4];		//provoking vertex（激发顶点）
   struct tgsi_interp_coef coef[PIPE_MAX_SHADER_INPUTS];
   struct tgsi_interp_coef posCoef;  //插值系数/* For Z, W */
   unsigned max_layer;              //支持的最大层数
   unsigned nr_vertex_attrs;	      //顶点属性数量
   #if DEBUG_FRAGS
   unsigned int numFragsEmitted;  /**< per primitive */
   unsigned int numFragsWritten;  /**< per primitive */
#endif
};

union fi {
   float f;
   int32_t i;
   uint32_t ui;
};

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