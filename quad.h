#ifndef __QUAD_H__
#define __QUAD_H__

#include <stdint.h>
#include "state.h"

#define TGSI_QUAD_SIZE    4  /* 4 pixel/quad */
//#define PIPE_MAX_COLOR_BUFS 1    //颜色buffer数量
#define TGSI_NUM_CHANNELS 4      //RGBA

/* 光栅化器生成 2x2 四边形片段并将它们提供给 them to
 * 当前 fp_machine（见下文）。
 * 请记住，Y=0=顶部，Y 沿窗口向下增加。
 */
#define QUAD_TOP_LEFT     0
#define QUAD_TOP_RIGHT    1
#define QUAD_BOTTOM_LEFT  2
#define QUAD_BOTTOM_RIGHT 3

#define MASK_TOP_LEFT     (1 << QUAD_TOP_LEFT)
#define MASK_TOP_RIGHT    (1 << QUAD_TOP_RIGHT)
#define MASK_BOTTOM_LEFT  (1 << QUAD_BOTTOM_LEFT)
#define MASK_BOTTOM_RIGHT (1 << QUAD_BOTTOM_RIGHT)
#define MASK_ALL          0xf


/**
 * Quad 阶段输入(pos, coverage, front/back face, etc)
 */
struct quad_header_input
{
   int x0, y0;                  //窗口坐标
   unsigned layer;              //层索引
   unsigned viewport_index;     //视口索引
   float coverage[TGSI_QUAD_SIZE];  //抗锯齿的片段覆盖
   unsigned facing:1;         //面 Front (0) , back (1)
   unsigned prim:2;           //基元类型< QUAD_PRIM_POINT, LINE, TRI >
};

/**
 * Quad 阶段输入/输出。
 */
struct quad_header_inout
{
   unsigned mask:4;             //覆盖值
};

/**
 * Quad 阶段输出(color & depth)
 */
struct quad_header_output
{
   /** colors in SOA format (rrrr, gggg, bbbb, aaaa) */
   float color[PIPE_MAX_COLOR_BUFS][TGSI_NUM_CHANNELS][TGSI_QUAD_SIZE]; //颜色
   float depth[TGSI_QUAD_SIZE];        //深度
   uint8_t stencil[TGSI_QUAD_SIZE];    //模板
};

/**
 * 对我们需要了解的有关 2x2 像素块的所有内容进行编码。
 * 使用“Channel-Serial”或“SoA”布局。
 */
struct quad_header {
   struct quad_header_input input;
   struct quad_header_inout inout;
   struct quad_header_output output;

   /* Redundant/duplicated:
    */
   //const struct tgsi_interp_coef *posCoef;
   //const struct tgsi_interp_coef *coef;
};

#endif