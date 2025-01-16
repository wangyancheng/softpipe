#include <stdbool.h>

#ifndef UNUSED
#define UNUSED __attribute__((unused))
#endif
#include <assert.h>
#include <stdio.h>
#include "setup.h"
#include <stdlib.h>
#include "quad_depth_test.h"
#include "quad_blend.h"
#include "quad.h"
#include <string.h>
#include "flush.h"
#include <math.h>

#if DEBUG_VERTS
static void
print_vertex(const struct setup_context *setup,
             const float (*v)[4])
{
   int i;
   printf("   Vertex: (%p)\n", (void *) v);
   for (i = 0; (unsigned)i < setup->nr_vertex_attrs; i++) {
      printf("     %d: %f %f %f %f\n",  i,
              v[i][0], v[i][1], v[i][2], v[i][3]);
      if (util_is_inf_or_nan(v[i][0])) {
         printf("   NaN!\n");
      }
   }
}
#endif

/**
 * 给定 X 或 Y 坐标，返回其所属的块/四边形(block/quad)坐标。
 */
static inline int
block(int x)
{
   return x & ~(2-1);
}

static inline int
block_x(int x)
{
   return x & ~(16-1);
}

#define CLAMP( X, MIN, MAX )  ( (X)>(MIN) ? ((X)>(MAX) ? (MAX) : (X)) : (MIN) )
/**
 * Render a horizontal span of quads
 */
static void
flush_spans(struct setup_context *setup)
{
   const int step = MAX_QUADS;
   const int xleft0 = setup->span.left[0];
   const int xleft1 = setup->span.left[1];
   const int xright0 = setup->span.right[0];
   const int xright1 = setup->span.right[1];
   struct quad_stage *pipe = setup->softpipe->quad.first;

   const int minleft = block_x(MIN2(xleft0, xleft1));
   const int maxright = MAX2(xright0, xright1);
   int x;

   /* process quads in horizontal chunks of 16 */
   for (x = minleft; x < maxright; x += step) {
      unsigned skip_left0 = CLAMP(xleft0 - x, 0, step);
      unsigned skip_left1 = CLAMP(xleft1 - x, 0, step);
      unsigned skip_right0 = CLAMP(x + step - xright0, 0, step);
      unsigned skip_right1 = CLAMP(x + step - xright1, 0, step);
      unsigned lx = x;
      unsigned q = 0;

      unsigned skipmask_left0 = (1U << skip_left0) - 1U;
      unsigned skipmask_left1 = (1U << skip_left1) - 1U;

      /* These calculations fail when step == 32 and skip_right == 0.
       */
      unsigned skipmask_right0 = ~0U << (unsigned)(step - skip_right0);
      unsigned skipmask_right1 = ~0U << (unsigned)(step - skip_right1);

      unsigned mask0 = ~skipmask_left0 & ~skipmask_right0;
      unsigned mask1 = ~skipmask_left1 & ~skipmask_right1;

      if (mask0 | mask1) {
         do {
            unsigned quadmask = (mask0 & 3) | ((mask1 & 3) << 2);
            if (quadmask) {
               setup->quad[q].input.x0 = lx;
               setup->quad[q].input.y0 = setup->span.y;
               setup->quad[q].input.facing = setup->facing;
               setup->quad[q].inout.mask = quadmask;
               setup->quad_ptrs[q] = &setup->quad[q];
               q++;
#if DEBUG_FRAGS
               setup->numFragsEmitted += util_bitcount(quadmask);
#endif
            }
            mask0 >>= 2;
            mask1 >>= 2;
            lx += 2;
         } while (mask0 | mask1);

         //pipe->run( pipe, setup->quad_ptrs, q );
         shade_quads(pipe, setup->quad_ptrs, q);
      }
   }


   setup->span.y = 0;
   setup->span.right[0] = 0;
   setup->span.right[1] = 0;
   setup->span.left[0] = 1000000;     /* greater than right[0] */
   setup->span.left[1] = 1000000;     /* greater than right[1] */
}

/**
* 计算常量系数（GL_FLAT shading）的a0
* 值来自vertex[slot][i]
* 结果将放入setup->coef[slot].a0[i]
*
* \param vertSlot -[in] 哪个属性槽
* \param i -[in] 槽的哪个组件（ 0...3)
*/
static void
const_coeff(struct setup_context *setup,
            struct tgsi_interp_coef *coef,
            unsigned int vertSlot, unsigned int i)
{
   assert(i <= 3);

   coef->dadx[i] = 0;
   coef->dady[i] = 0;

   /* need provoking vertex info!
    */
   coef->a0[i] = setup->vprovoke[vertSlot][i];
}

/**
 * 透视插值计算
 * 
 * 
 */
static void
point_persp_coeff(UNUSED const struct setup_context *setup,
                  const float (*vert)[4],
                  struct tgsi_interp_coef *coef,
                  unsigned int vertSlot, unsigned int i)
{
   assert(i <= 3);
   coef->dadx[i] = 0.0F;
   coef->dady[i] = 0.0F;
   coef->a0[i] = vert[vertSlot][i] * vert[0][3];
}

/**
 * Special coefficient setup for gl_FragCoord.
 * X and Y are trivial, though Y may have to be inverted for OpenGL.
 * Z and W are copied from posCoef which should have already been computed.
 * We could do a bit less work if we'd examine gl_FragCoord's swizzle mask.
 */
static void
setup_fragcoord_coeff(struct setup_context *setup, unsigned int slot)
{
   const struct tgsi_shader_info *fsInfo = &setup->softpipe->fs_variant->info;
   bool origin_lower_left =
         fsInfo->properties[TGSI_PROPERTY_FS_COORD_ORIGIN];
   bool pixel_center_integer =
         fsInfo->properties[TGSI_PROPERTY_FS_COORD_PIXEL_CENTER];

   /*X*/
   setup->coef[slot].a0[0] = pixel_center_integer ? 0.0f : 0.5f;
   setup->coef[slot].dadx[0] = 1.0f;
   setup->coef[slot].dady[0] = 0.0f;
   /*Y*/
   setup->coef[slot].a0[1] =
		   (origin_lower_left ? setup->softpipe->framebuffer.height-1 : 0)
		   + (pixel_center_integer ? 0.0f : 0.5f);
   setup->coef[slot].dadx[1] = 0.0f;
   setup->coef[slot].dady[1] = origin_lower_left ? -1.0f : 1.0f;
   /*Z*/
   setup->coef[slot].a0[2] = setup->posCoef.a0[2];
   setup->coef[slot].dadx[2] = setup->posCoef.dadx[2];
   setup->coef[slot].dady[2] = setup->posCoef.dady[2];
   /*W*/
   setup->coef[slot].a0[3] = setup->posCoef.a0[3];
   setup->coef[slot].dadx[3] = setup->posCoef.dadx[3];
   setup->coef[slot].dady[3] = setup->posCoef.dady[3];
}

/**
* 根据scissor和surfce bounds 裁剪
*
* \param setup -[in/out] setup上下文
* \param quad - [in] 要渲染的四变形
*/
static inline void
quad_clip(struct setup_context *setup, struct quad_header *quad)
{
   unsigned viewport_index = quad[0].input.viewport_index;
   const struct pipe_scissor_state *cliprect = &setup->softpipe->cliprect[viewport_index];
   const int minx = (int) cliprect->minx;
   const int maxx = (int) cliprect->maxx;
   const int miny = (int) cliprect->miny;
   const int maxy = (int) cliprect->maxy;

   if (quad->input.x0 >= maxx ||
       quad->input.y0 >= maxy ||
       quad->input.x0 + 1 < minx ||
       quad->input.y0 + 1 < miny) {
      /* totally clipped */
      quad->inout.mask = 0x0;
      return;
   }
   if (quad->input.x0 < minx)
      quad->inout.mask &= (MASK_BOTTOM_RIGHT | MASK_TOP_RIGHT);
   if (quad->input.y0 < miny)
      quad->inout.mask &= (MASK_BOTTOM_LEFT | MASK_BOTTOM_RIGHT);
   if (quad->input.x0 == maxx - 1)
      quad->inout.mask &= (MASK_BOTTOM_LEFT | MASK_TOP_LEFT);
   if (quad->input.y0 == maxy - 1)
      quad->inout.mask &= (MASK_TOP_LEFT | MASK_TOP_RIGHT);
}

/**
 * 发射一个quad（传递到下一阶段）并进行剪辑。
 * 
 * \param setup -[in/out] setup上下文
 * \param quad -[in] 要渲染的四变形
 */
static inline void
clip_emit_quad(struct setup_context *setup, struct quad_header *quad)
{
   quad_clip(setup, quad);

   if (quad->inout.mask) {
    struct softpipe_context *sp = setup->softpipe;

#if DEBUG_FRAGS
      setup->numFragsEmitted += util_bitcount(quad->inout.mask);
#endif

      shade_quads(sp->quad.first, &quad, 1);
   }
}

/**
 * 计算A0，dadx和dady，以线性插值系数，
 * 为了一行。
 * V [0]和V [1]分别为Vmin和Vmax。
 */
static void
line_linear_coeff(const struct setup_context *setup,
                  struct tgsi_interp_coef *coef,
                  uint i,
                  const float v[2])
{
   const float da = v[1] - v[0];
   const float dadx = da * setup->emaj.dx * setup->oneoverarea;
   const float dady = da * setup->emaj.dy * setup->oneoverarea;
   coef->dadx[i] = dadx;
   coef->dady[i] = dady;
   coef->a0[i] = (v[0] -
                  (dadx * (setup->vmin[0][0] - setup->pixel_offset) +
                   dady * (setup->vmin[0][1] - setup->pixel_offset)));
}

/* 如果启用了V0，将圆柱形包裹应用于V0，V1坐标。
 * 输入坐标必须在[0，1]范围内，否则结果不确定。
 */
static void
line_apply_cylindrical_wrap(float v0,
                            float v1,
                            uint cylindrical_wrap,
                            float output[2])
{
   if (cylindrical_wrap) {
      float delta;

      delta = v1 - v0;
      if (delta > 0.5f) {
         v0 += 1.0f;
      }
      else if (delta < -0.5f) {
         v1 += 1.0f;
      }
   }

   output[0] = v0;
   output[1] = v1;
}

/**
 * Compute a0, dadx and dady for a perspective-corrected interpolant,
 * for a line.
 * v[0] and v[1] are vmin and vmax, respectively.
 */
static void
line_persp_coeff(const struct setup_context *setup,
                 struct tgsi_interp_coef *coef,
                 uint i,
                 const float v[2])
{
   const float a0 = v[0] * setup->vmin[0][3];
   const float a1 = v[1] * setup->vmax[0][3];
   const float da = a1 - a0;
   const float dadx = da * setup->emaj.dx * setup->oneoverarea;
   const float dady = da * setup->emaj.dy * setup->oneoverarea;
   coef->dadx[i] = dadx;
   coef->dady[i] = dady;
   coef->a0[i] = (a0 -
                  (dadx * (setup->vmin[0][0] - setup->pixel_offset) +
                   dady * (setup->vmin[0][1] - setup->pixel_offset)));
}

/**
 * 计算设置 - > coef [] array dadx，dady，a0值。
 * 必须在设置 - > vmin后调用VMAX初始化。
 * 
 * \param setup -[in/out] setup上下文
 * \param coef -[in/out] 插值系数
 * \param i -[in] 组件索引
 * \param v -[in] 顶点数据
 */
static bool
setup_line_coefficients(struct setup_context *setup,
                        const float (*v0)[4],
                        const float (*v1)[4])
{
   struct softpipe_context *softpipe = setup->softpipe;
   const struct tgsi_shader_info *fsInfo = &setup->softpipe->fs_variant->info;  //片段着色器信息
   const struct sp_setup_info *sinfo = &softpipe->setup_info;   //顶点设置信息
   uint fragSlot;
   float area;
   float v[2];          //线段的两个顶点

   assert(sinfo->valid);

   /* use setup->vmin, vmax to point to vertices */
   if (softpipe->rasterizer->flatshade_first)
      setup->vprovoke = v0;
   else
      setup->vprovoke = v1;
   setup->vmin = v0;
   setup->vmax = v1;

   setup->emaj.dx = setup->vmax[0][0] - setup->vmin[0][0];      //线段的x长度
   setup->emaj.dy = setup->vmax[0][1] - setup->vmin[0][1];      //线段的y长度

   /* 注意：这实际上不是面积，而是与面积成比例的。 */
   area = setup->emaj.dx * setup->emaj.dx + setup->emaj.dy * setup->emaj.dy;
   if (area == 0.0f || util_is_inf_or_nan(area))
      return false;
   setup->oneoverarea = 1.0f / area;

   /* z 和 w 通过线性插值完成：    */
   v[0] = setup->vmin[0][2];            //线段的z坐标
   v[1] = setup->vmax[0][2];
   line_linear_coeff(setup, &setup->posCoef, 2, v);

   v[0] = setup->vmin[0][3];
   v[1] = setup->vmax[0][3];
   line_linear_coeff(setup, &setup->posCoef, 3, v);

   /* 为所有剩余属性设置插值：
    */
   for (fragSlot = 0; fragSlot < fsInfo->num_inputs; fragSlot++) {
      const uint vertSlot = sinfo->attrib[fragSlot].src_index;
      uint j;

      switch (sinfo->attrib[fragSlot].interp) {
      case SP_INTERP_CONSTANT:              //常量插值
         for (j = 0; j < TGSI_NUM_CHANNELS; j++)
            const_coeff(setup, &setup->coef[fragSlot], vertSlot, j);
         break;
      case SP_INTERP_LINEAR:                 //线性插值
         for (j = 0; j < TGSI_NUM_CHANNELS; j++) {
            line_apply_cylindrical_wrap(setup->vmin[vertSlot][j],
                                        setup->vmax[vertSlot][j],
                                        fsInfo->input_cylindrical_wrap[fragSlot] & (1 << j),
                                        v);
            line_linear_coeff(setup, &setup->coef[fragSlot], j, v);
         }
         break;
      case SP_INTERP_PERSPECTIVE:            //透视插值
         for (j = 0; j < TGSI_NUM_CHANNELS; j++) {
            line_apply_cylindrical_wrap(setup->vmin[vertSlot][j],
                                        setup->vmax[vertSlot][j],
                                        fsInfo->input_cylindrical_wrap[fragSlot] & (1 << j),
                                        v);
            line_persp_coeff(setup, &setup->coef[fragSlot], j, v);
         }
         break;
      case SP_INTERP_POS:                   //点精灵，纹理坐标插值
         setup_fragcoord_coeff(setup, fragSlot);
         break;
      default:
         assert(0);
      }

      if (fsInfo->input_semantic_name[fragSlot] == TGSI_SEMANTIC_FACE) {
         /* convert 0 to 1.0 and 1 to -1.0 */
         setup->coef[fragSlot].a0[0] = setup->facing * -2.0f + 1.0f;
         setup->coef[fragSlot].dadx[0] = 0.0;
         setup->coef[fragSlot].dady[0] = 0.0;
      }
   }
   return true;
}

/**
 * 在线段中绘制一个像素。
 * 
 * \param setup -[in/out] setup上下文
 * \param x -[in] 线的x坐标
 * \param y -[in] 线的y坐标
 */
static inline void
plot(struct setup_context *setup, int x, int y)
{
   const int iy = y & 1;
   const int ix = x & 1;
   const int quadX = x - ix;
   const int quadY = y - iy;
   const int mask = (1 << ix) << (2 * iy);

   if (quadX != setup->quad[0].input.x0 ||
       quadY != setup->quad[0].input.y0)
   {
      /* flush prev quad, start new quad */

      if (setup->quad[0].input.x0 != -1)
         clip_emit_quad(setup, &setup->quad[0]);

      setup->quad[0].input.x0 = quadX;
      setup->quad[0].input.y0 = quadY;
      setup->quad[0].inout.mask = 0x0;
   }

   setup->quad[0].inout.mask |= mask;
}

/**
 * 线光栅化，单位线宽度。线宽大于1被渲染成两个三角形。
 * 
 * \param setup -[in/out] setup上下文
 * \param v0 -[in] 线的第一个顶点
 * \param v1 -[in] 线的第二个顶点
 */
void
sp_setup_line(struct setup_context *setup,
              const float (*v0)[4],
              const float (*v1)[4])
{
   int x0 = (int) v0[0][0];
   int x1 = (int) v1[0][0];
   int y0 = (int) v0[0][1];
   int y1 = (int) v1[0][1];
   int dx = x1 - x0;
   int dy = y1 - y0;
   int xstep, ystep;
   uint layer = 0;
   unsigned viewport_index = 0;

#if DEBUG_VERTS
   printf("Setup line:\n");
   print_vertex(setup, v0);
   print_vertex(setup, v1);
#endif

   if (setup->softpipe->rasterizer->rasterizer_discard)
      return;

   if (dx == 0 && dy == 0)
      return;
    
   if (!setup_line_coefficients(setup, v0, v1))
      return;

   assert(v0[0][0] < 1.0e9);
   assert(v0[0][1] < 1.0e9);
   assert(v1[0][0] < 1.0e9);
   assert(v1[0][1] < 1.0e9);

   if (dx < 0) {
      dx = -dx;   /* make positive */
      xstep = -1;
   }
   else {
      xstep = 1;
   }

   if (dy < 0) {
      dy = -dy;   /* make positive */
      ystep = -1;
   }
   else {
      ystep = 1;
   }

   assert(dx >= 0);
   assert(dy >= 0);
   assert(setup->softpipe->reduced_prim == PIPE_PRIM_LINES);

   setup->quad[0].input.x0 = setup->quad[0].input.y0 = -1;
   setup->quad[0].inout.mask = 0x0;
   if (setup->softpipe->layer_slot > 0) {
      layer = *(unsigned *)setup->vprovoke[setup->softpipe->layer_slot];
      layer = MIN2(layer, setup->max_layer);
   }
   setup->quad[0].input.layer = layer;

   if (setup->softpipe->viewport_index_slot > 0) {
      unsigned *udata = (unsigned*)setup->vprovoke[setup->softpipe->viewport_index_slot];
      viewport_index = sp_clamp_viewport_idx(*udata);
   }
   setup->quad[0].input.viewport_index = viewport_index;

   /* XXX temporary: set coverage to 1.0 so the line appears
    * if AA mode happens to be enabled.
    */
   setup->quad[0].input.coverage[0] =
   setup->quad[0].input.coverage[1] =
   setup->quad[0].input.coverage[2] =
   setup->quad[0].input.coverage[3] = 1.0;

   if (dx > dy) {
      /*** X-major line ***/
      int i;
      const int errorInc = dy + dy;
      int error = errorInc - dx;
      const int errorDec = error - dx;

      for (i = 0; i < dx; i++) {
         plot(setup, x0, y0);

         x0 += xstep;
         if (error < 0) {
            error += errorInc;
         }
         else {
            error += errorDec;
            y0 += ystep;
         }
      }
   }
   else {
      /*** Y-major line ***/
      int i;
      const int errorInc = dx + dx;
      int error = errorInc - dy;
      const int errorDec = error - dy;

      for (i = 0; i < dy; i++) {
         plot(setup, x0, y0);

         y0 += ystep;
         if (error < 0) {
            error += errorInc;
         }
         else {
            error += errorDec;
            x0 += xstep;
         }
      }
   }

   /* draw final quad */
   if (setup->quad[0].inout.mask) {
      clip_emit_quad(setup, &setup->quad[0]);
   }
}


/**
 * setup点
 * \param setup -[in/out] setup上下文
 * \param v0 -[in] 定点数据
 */
void
setup_point(struct setup_context *setup,
               const float (*v0)[4])
{
    struct softpipe_context *softpipe = setup->softpipe;		//softpipe上下文
    const struct tgsi_shader_info *fsInfo = &setup->softpipe->fs_variant->info;	//片段着色器信息
    const int sizeAttr = setup->softpipe->psize_slot;		//顶点着色器输出槽是否包括点大小
    const float size							//点的大小
        = sizeAttr > 0 ? v0[sizeAttr][0]
        : setup->softpipe->rasterizer->point_size;
    const float halfSize = 0.5F * size;
    const bool round = (bool) setup->softpipe->rasterizer->point_smooth;   //点抗锯齿开启后，点渲染成圆形（正常为正方形）
    const float x = v0[0][0];  /* Note: data[0] is always position */    //点的屏幕坐标位置
    const float y = v0[0][1];
    const struct sp_setup_info *sinfo = &softpipe->setup_info;    //  顶点的设置信息，主要包括插值类型
    unsigned int fragSlot;
    unsigned int layer = 0;
    unsigned viewport_index = 0;
#if DEBUG_VERTS
    printf("Setup point:\n");
    print_vertex(setup, v0);
#endif

    assert(sinfo->valid);

    if (setup->softpipe->rasterizer->rasterizer_discard)
        return;

    assert(setup->softpipe->reduced_prim == PIPE_PRIM_POINTS);

    //1. 设置quad[0].input layer和viewport_index参数
    if (setup->softpipe->layer_slot > 0) {
        layer = *(unsigned *)v0[setup->softpipe->layer_slot];
        layer = MIN2(layer, setup->max_layer);
    }
    setup->quad[0].input.layer = layer;

    if (setup->softpipe->viewport_index_slot > 0) {
        unsigned *udata = (unsigned*)v0[setup->softpipe->viewport_index_slot];
        viewport_index = sp_clamp_viewport_idx(*udata);
    }
    setup->quad[0].input.viewport_index = viewport_index;

    /**
    * 对于点，所有插值都是常量值。
    * 但是，对于点精灵，我们需要适当的设置纹理坐标。
    * XXX：哪些系数是纹理坐标？
    * 我们可以将点精灵做成纹理四边形......
    *
    * KW：我们不知道哪些系数是纹理坐标-最终，每个属性使用哪种插值模式的选择应该由片段程序决定，
    * 使用包含插值模式作为参数的每个属性声明语句。因此，要么片段程序必须针对点精灵与正常点
    * 行为进行调整，要么必须定义一种特殊的插值模式来匹配点精灵所需的行为。
    * 但是-后者不是普通硬件的功能，因此可能应该在此基础上排除。
    */
    //2. 计算setup->posCoef（插值系数）
    
    setup->vprovoke = v0;     //provoking vertex（激发顶点）

    // setup Z, W 
    const_coeff(setup, &setup->posCoef, 0, 2);
    const_coeff(setup, &setup->posCoef, 0, 3);

    for (fragSlot = 0; fragSlot < fsInfo->num_inputs; fragSlot++) {
        const unsigned int vertSlot = sinfo->attrib[fragSlot].src_index;
        unsigned int j;

        switch (sinfo->attrib[fragSlot].interp) {
        case SP_INTERP_CONSTANT:    //常量插值
            // fall-through 
        case SP_INTERP_LINEAR:		//线性插值
            for (j = 0; j < TGSI_NUM_CHANNELS; j++)
                const_coeff(setup, &setup->coef[fragSlot], vertSlot, j);
            break;
        case SP_INTERP_PERSPECTIVE:			//透视插值
            for (j = 0; j < TGSI_NUM_CHANNELS; j++)
                point_persp_coeff(setup, setup->vprovoke,
                                    &setup->coef[fragSlot], vertSlot, j);
            break;
        case SP_INTERP_POS:		//点精灵，纹理坐标插值
            setup_fragcoord_coeff(setup, fragSlot);
            break;
        default:
            assert(0);
        }

        if (fsInfo->input_semantic_name[fragSlot] == TGSI_SEMANTIC_FACE) {
            // convert 0 to 1.0 and 1 to -1.0 
            setup->coef[fragSlot].a0[0] = setup->facing * -2.0f + 1.0f;
            setup->coef[fragSlot].dadx[0] = 0.0;
            setup->coef[fragSlot].dady[0] = 0.0;
        }
    }
    

    //3.1 渲染成一个点
    if (halfSize <= 0.5 && !round) {
        /* special case for 1-pixel points */
        const int ix = ((int) x) & 1;
        const int iy = ((int) y) & 1;
        setup->quad[0].input.x0 = (int) x - ix;
        setup->quad[0].input.y0 = (int) y - iy;
        setup->quad[0].inout.mask = (1 << ix) << (2 * iy);        //???
        clip_emit_quad(setup, &setup->quad[0]);
    }
    else {
        //3.2 渲染成圆
        if (round) {
            /* rounded points */
            const int ixmin = block((int) (x - halfSize));
            const int ixmax = block((int) (x + halfSize));
            const int iymin = block((int) (y - halfSize));
            const int iymax = block((int) (y + halfSize));
            const float rmin = halfSize - 0.7071F;  /* 0.7071 = sqrt(2)/2 */
            const float rmax = halfSize + 0.7071F;
            const float rmin2 = MAX2(0.0F, rmin * rmin);
            const float rmax2 = rmax * rmax;
            const float cscale = 1.0F / (rmax2 - rmin2);
            int ix, iy;

            for (iy = iymin; iy <= iymax; iy += 2) {
                for (ix = ixmin; ix <= ixmax; ix += 2) {
                    float dx, dy, dist2, cover;

                    setup->quad[0].inout.mask = 0x0;    //圆外mask是0x0

                    dx = (ix + 0.5f) - x;
                    dy = (iy + 0.5f) - y;
                    dist2 = dx * dx + dy * dy;
                    if (dist2 <= rmax2) {       //在圆内
                        cover = 1.0F - (dist2 - rmin2) * cscale;
                        setup->quad[0].input.coverage[QUAD_TOP_LEFT] = MIN2(cover, 1.0f);
                        setup->quad[0].inout.mask |= MASK_TOP_LEFT;
                    }

                    dx = (ix + 1.5f) - x;
                    dy = (iy + 0.5f) - y;
                    dist2 = dx * dx + dy * dy;
                    if (dist2 <= rmax2) {
                        cover = 1.0F - (dist2 - rmin2) * cscale;
                        setup->quad[0].input.coverage[QUAD_TOP_RIGHT] = MIN2(cover, 1.0f);
                        setup->quad[0].inout.mask |= MASK_TOP_RIGHT;
                    }

                    dx = (ix + 0.5f) - x;
                    dy = (iy + 1.5f) - y;
                    dist2 = dx * dx + dy * dy;
                    if (dist2 <= rmax2) {
                        cover = 1.0F - (dist2 - rmin2) * cscale;
                        setup->quad[0].input.coverage[QUAD_BOTTOM_LEFT] = MIN2(cover, 1.0f);
                        setup->quad[0].inout.mask |= MASK_BOTTOM_LEFT;
                    }

                    dx = (ix + 1.5f) - x;
                    dy = (iy + 1.5f) - y;
                    dist2 = dx * dx + dy * dy;
                    if (dist2 <= rmax2) {
                        cover = 1.0F - (dist2 - rmin2) * cscale;
                        setup->quad[0].input.coverage[QUAD_BOTTOM_RIGHT] = MIN2(cover, 1.0f);
                        setup->quad[0].inout.mask |= MASK_BOTTOM_RIGHT;
                    }

                    if (setup->quad[0].inout.mask) {
                        setup->quad[0].input.x0 = ix;
                        setup->quad[0].input.y0 = iy;
                        clip_emit_quad(setup, &setup->quad[0]);
                    }
                }
            }
        }
        //3.3 渲染成正方形
        else {
            /* square points */
            const int xmin = (int) (x + 0.75 - halfSize);
            const int ymin = (int) (y + 0.25 - halfSize);
            const int xmax = xmin + (int) size;
            const int ymax = ymin + (int) size;
            /* XXX 现在可以将剪刀应用于 xmin、ymin、xmax、ymax*/
            const int ixmin = block(xmin);
            const int ixmax = block(xmax - 1);
            const int iymin = block(ymin);
            const int iymax = block(ymax - 1);
            int ix, iy;

            if(1){
                printf("(%f, %f) -> X:%d..%d Y:%d..%d\n", x, y, xmin, xmax,ymin,ymax);
            }

            for (iy = iymin; iy <= iymax; iy += 2) {
                unsigned int rowMask = 0xf;
                if (iy < ymin) {
                    /* above the top edge */
                    rowMask &= (MASK_BOTTOM_LEFT | MASK_BOTTOM_RIGHT);
                }
                if (iy + 1 >= ymax) {
                    /* below the bottom edge */
                    rowMask &= (MASK_TOP_LEFT | MASK_TOP_RIGHT);
                }

                for (ix = ixmin; ix <= ixmax; ix += 2) {
                    unsigned int mask = rowMask;

                    if (ix < xmin) {
                        /* fragment is past left edge of point, turn off left bits */
                        mask &= (MASK_BOTTOM_RIGHT | MASK_TOP_RIGHT);
                    }
                    if (ix + 1 >= xmax) {
                        /* past the right edge */
                        mask &= (MASK_BOTTOM_LEFT | MASK_TOP_LEFT);
                    }

                    setup->quad[0].inout.mask = mask;
                    setup->quad[0].input.x0 = ix;
                    setup->quad[0].input.y0 = iy;
                    clip_emit_quad(setup, &setup->quad[0]);
                }
            }
        }
    }
}

/**
 * 计算a0，dadx和dady，以线性插值系数，
 * 用于三角形。
 * V[0]，V[1]和V[2]分别为VMIN，VMID和VMAX。
 * 
 * \param setup -[in/out] setup上下文
 * \param coef -[in/out] 插值系数
 * \param i -[in] 组件索引
 * \param v -[in] 顶点数据
 */
static void
tri_linear_coeff(struct setup_context *setup,
                 struct tgsi_interp_coef *coef,
                 uint i,
                 const float v[3])
{
   float botda = v[1] - v[0];       //底边的深度差值
   float majda = v[2] - v[0];       //主边的深度差值
   float a = setup->ebot.dy * majda - botda * setup->emaj.dy;     //
   float b = setup->emaj.dx * botda - majda * setup->ebot.dx;     //
   float dadx = a * setup->oneoverarea;      //计算dadx
   float dady = b * setup->oneoverarea;      //计算dady

   assert(i <= 3);

   coef->dadx[i] = dadx;
   coef->dady[i] = dady;

   /* 
    * 计算 a0 作为在 (0,0) 处采样得到的片段的值。
    * 考虑到我们希望在像素中心 (pixel_offset, pixel_offset) 处进行采样。
    * 这种方法对于 dadx 或 dady 值非常大的三角形来说效果不佳，
    * 因为这将导致从 a0 中减去并重新加上一个非常大的数字，
    * 这意味着我们将丢失 a0 中的大量小数位精度。
    * 解决这个问题的方法是将 a0 定义为靠近 vmin 的某个像素中心的采样值。
    * 我稍后将切换到这种方法。
    */
   coef->a0[i] = (v[0] -
                  (dadx * (setup->vmin[0][0] - setup->pixel_offset) +
                   dady * (setup->vmin[0][1] - setup->pixel_offset)));
}

/* Apply cylindrical wrapping to v0, v1, v2 coordinates, if enabled.
 * Input coordinates must be in [0, 1] range, otherwise results are undefined.
 * Some combinations of coordinates produce invalid results,
 * but this behaviour is acceptable.
 */
static void
tri_apply_cylindrical_wrap(float v0,
                           float v1,
                           float v2,
                           uint cylindrical_wrap,
                           float output[3])
{
   if (cylindrical_wrap) {
      float delta;

      delta = v1 - v0;
      if (delta > 0.5f) {
         v0 += 1.0f;
      }
      else if (delta < -0.5f) {
         v1 += 1.0f;
      }

      delta = v2 - v1;
      if (delta > 0.5f) {
         v1 += 1.0f;
      }
      else if (delta < -0.5f) {
         v2 += 1.0f;
      }

      delta = v0 - v2;
      if (delta > 0.5f) {
         v2 += 1.0f;
      }
      else if (delta < -0.5f) {
         v0 += 1.0f;
      }
   }

   output[0] = v0;
   output[1] = v1;
   output[2] = v2;
}

/**
 * Compute a0, dadx and dady for a perspective-corrected interpolant,
 * for a triangle.
 * We basically multiply the vertex value by 1/w before computing
 * the plane coefficients (a0, dadx, dady).
 * Later, when we compute the value at a particular fragment position we'll
 * divide the interpolated value by the interpolated W at that fragment.
 * v[0], v[1] and v[2] are vmin, vmid and vmax, respectively.
 */
static void
tri_persp_coeff(struct setup_context *setup,
                struct tgsi_interp_coef *coef,
                uint i,
                const float v[3])
{
   /* premultiply by 1/w  (v[0][3] is always W):
    */
   float mina = v[0] * setup->vmin[0][3];
   float mida = v[1] * setup->vmid[0][3];
   float maxa = v[2] * setup->vmax[0][3];
   float botda = mida - mina;
   float majda = maxa - mina;
   float a = setup->ebot.dy * majda - botda * setup->emaj.dy;
   float b = setup->emaj.dx * botda - majda * setup->ebot.dx;
   float dadx = a * setup->oneoverarea;
   float dady = b * setup->oneoverarea;

   assert(i <= 3);

   coef->dadx[i] = dadx;
   coef->dady[i] = dady;
   coef->a0[i] = (mina -
                  (dadx * (setup->vmin[0][0] - setup->pixel_offset) +
                   dady * (setup->vmin[0][1] - setup->pixel_offset)));
}


/**
 * 计算 setup->coef[] 数组中的 dadx、dady、a0 值。
 * 必须在 setup->vmin、vmid、vmax、vprovoke 初始化之后调用。
 * 
 * \param setup -[in/out] setup上下文
 */
static void
setup_tri_coefficients(struct setup_context *setup)
{
   struct softpipe_context *softpipe = setup->softpipe;
   const struct tgsi_shader_info *fsInfo = &setup->softpipe->fs_variant->info;
   const struct sp_setup_info *sinfo = &softpipe->setup_info;
   uint fragSlot;
   float v[3];

   assert(sinfo->valid);

   /* Z和W是通过线性插值完成的：
    */
   v[0] = setup->vmin[0][2];     //顶点最小值的z坐标
   v[1] = setup->vmid[0][2];     //顶点中间值的z坐标
   v[2] = setup->vmax[0][2];     //顶点最大值的z坐标
   tri_linear_coeff(setup, &setup->posCoef, 2, v);

   v[0] = setup->vmin[0][3];     //顶点最小值的w坐标
   v[1] = setup->vmid[0][3];     //顶点中间值的w坐标
   v[2] = setup->vmax[0][3];     //顶点最大值的w坐标
   tri_linear_coeff(setup, &setup->posCoef, 3, v);

   /* 所有其余属性的设置插值:
    */
   for (fragSlot = 0; fragSlot < fsInfo->num_inputs; fragSlot++) {
      const uint vertSlot = sinfo->attrib[fragSlot].src_index;
      uint j;

      switch (sinfo->attrib[fragSlot].interp) {
      //常量插值
      case SP_INTERP_CONSTANT:
         for (j = 0; j < TGSI_NUM_CHANNELS; j++) {
            const_coeff(setup, &setup->coef[fragSlot], vertSlot, j);
         }
         break;
      //线性插值
      case SP_INTERP_LINEAR:
         for (j = 0; j < TGSI_NUM_CHANNELS; j++) {
            tri_apply_cylindrical_wrap(setup->vmin[vertSlot][j],
                                       setup->vmid[vertSlot][j],
                                       setup->vmax[vertSlot][j],
                                       fsInfo->input_cylindrical_wrap[fragSlot] & (1 << j),
                                       v);
            tri_linear_coeff(setup, &setup->coef[fragSlot], j, v);
         }
         break;
      case SP_INTERP_PERSPECTIVE:
         for (j = 0; j < TGSI_NUM_CHANNELS; j++) {
            tri_apply_cylindrical_wrap(setup->vmin[vertSlot][j],
                                       setup->vmid[vertSlot][j],
                                       setup->vmax[vertSlot][j],
                                       fsInfo->input_cylindrical_wrap[fragSlot] & (1 << j),
                                       v);
            tri_persp_coeff(setup, &setup->coef[fragSlot], j, v);
         }
         break;
      //
      case SP_INTERP_POS:
         setup_fragcoord_coeff(setup, fragSlot);
         break;
      default:
         assert(0);
      }

      if (fsInfo->input_semantic_name[fragSlot] == TGSI_SEMANTIC_FACE) {
         /* convert 0 to 1.0 and 1 to -1.0 */
         setup->coef[fragSlot].a0[0] = setup->facing * -2.0f + 1.0f;
         setup->coef[fragSlot].dadx[0] = 0.0;
         setup->coef[fragSlot].dady[0] = 0.0;
      }

      if (0) {
         for (j = 0; j < TGSI_NUM_CHANNELS; j++) {
            printf("attr[%d].%c: a0:%f dx:%f dy:%f\n",
                         fragSlot, "xyzw"[j],
                         setup->coef[fragSlot].a0[j],
                         setup->coef[fragSlot].dadx[j],
                         setup->coef[fragSlot].dady[j]);
         }
      }
   }
}


static void
setup_tri_edges(struct setup_context *setup)
{
   float vmin_x = setup->vmin[0][0] + setup->pixel_offset;
   float vmid_x = setup->vmid[0][0] + setup->pixel_offset;

   float vmin_y = setup->vmin[0][1] - setup->pixel_offset;
   float vmid_y = setup->vmid[0][1] - setup->pixel_offset;
   float vmax_y = setup->vmax[0][1] - setup->pixel_offset;

   setup->emaj.sy = ceilf(vmin_y);
   setup->emaj.lines = (int) ceilf(vmax_y - setup->emaj.sy);
   setup->emaj.dxdy = setup->emaj.dy ? setup->emaj.dx / setup->emaj.dy : .0f;
   setup->emaj.sx = vmin_x + (setup->emaj.sy - vmin_y) * setup->emaj.dxdy;

   setup->etop.sy = ceilf(vmid_y);
   setup->etop.lines = (int) ceilf(vmax_y - setup->etop.sy);
   setup->etop.dxdy = setup->etop.dy ? setup->etop.dx / setup->etop.dy : .0f;
   setup->etop.sx = vmid_x + (setup->etop.sy - vmid_y) * setup->etop.dxdy;

   setup->ebot.sy = ceilf(vmin_y);
   setup->ebot.lines = (int) ceilf(vmid_y - setup->ebot.sy);
   setup->ebot.dxdy = setup->ebot.dy ? setup->ebot.dx / setup->ebot.dy : .0f;
   setup->ebot.sx = vmin_x + (setup->ebot.sy - vmin_y) * setup->ebot.dxdy;
}


/**
 * 渲染三角形的上半部分或下半部分。
 * 这里也应用了剪裁/裁剪矩形。
 */
static void
subtriangle(struct setup_context *setup,
            struct edge *eleft,
            struct edge *eright,
            int lines,
            unsigned viewport_index)
{
   const struct pipe_scissor_state *cliprect = &setup->softpipe->cliprect[viewport_index];
   const int minx = (int) cliprect->minx;
   const int maxx = (int) cliprect->maxx;
   const int miny = (int) cliprect->miny;
   const int maxy = (int) cliprect->maxy;
   int y, start_y, finish_y;
   int sy = (int)eleft->sy;

   assert((int)eleft->sy == (int) eright->sy);
   assert(lines >= 0);

   /* clip top/bottom */
   start_y = sy;
   if (start_y < miny)
      start_y = miny;

   finish_y = sy + lines;
   if (finish_y > maxy)
      finish_y = maxy;

   start_y -= sy;
   finish_y -= sy;

   /*
   debug_printf("%s %d %d\n", __FUNCTION__, start_y, finish_y);
   */

   for (y = start_y; y < finish_y; y++) {

      /* avoid accumulating adds as floats don't have the precision to
       * accurately iterate large triangle edges that way.  luckily we
       * can just multiply these days.
       *
       * this is all drowned out by the attribute interpolation anyway.
       */
      int left = (int)(eleft->sx + y * eleft->dxdy);
      int right = (int)(eright->sx + y * eright->dxdy);

      /* clip left/right */
      if (left < minx)
         left = minx;
      if (right > maxx)
         right = maxx;

      if (left < right) {
         int _y = sy + y;
         if (block(_y) != setup->span.y) {
            flush_spans(setup);
            setup->span.y = block(_y);
         }

         setup->span.left[_y&1] = left;
         setup->span.right[_y&1] = right;
      }
   }


   /* save the values so that emaj can be restarted:
    */
   eleft->sx += lines * eleft->dxdy;
   eright->sx += lines * eright->dxdy;
   eleft->sy += lines;
   eright->sy += lines;
}

/**
 * 计算三角形的面积。
 * 通过叉乘计算三角形的面积。
 * 三角形的面积等于叉乘的z值的一半。
 * 如果面积为正，则三角形为逆时针方向。
 * 如果面积为负，则三角形为顺时针方向。
 *
 * \param v0 -[in] 三角形的第一个顶点
 * \param v1 -[in] 三角形的第二个顶点
 * \param v2 -[in] 三角形的第三个顶点
 * \return 三角形的面积
 */
static float
calc_det(const float (*v0)[4],
         const float (*v1)[4],
         const float (*v2)[4])
{
   /* edge vectors e = v0 - v2, f = v1 - v2 */
   const float ex = v0[0][0] - v2[0][0];     //v0的x坐标 - v2的x坐标
   const float ey = v0[0][1] - v2[0][1];     //v0的y坐标 - v2的y坐标
   const float fx = v1[0][0] - v2[0][0];     //v1的x坐标 - v2的x坐标
   const float fy = v1[0][1] - v2[0][1];     //v1的y坐标 - v2的y坐标

   /* det = cross(e,f).z */
   return ex * fy - ey * fx;        //叉乘
}

/**
**
 * 将顶点从上到下排序，设置三角形
 * 边缘字段（Ebot，Emaj，Etop）。
 * 
 * \param setup -[in/out] setup上下文
 * \param det -[in] 三角形的面积
 * \param v0 -[in] 三角形的第一个顶点
 * \param v1 -[in] 三角形的第二个顶点
 * \param v2 -[in] 三角形的第三个顶点
 * 
 * \return FALSE if coords are inf/nan (cull the tri), TRUE otherwise
 */
static bool
setup_sort_vertices(struct setup_context *setup,
                    float det,
                    const float (*v0)[4],
                    const float (*v1)[4],
                    const float (*v2)[4])
{
   //1. 设置激发顶点
   if (setup->softpipe->rasterizer->flatshade_first)
      setup->vprovoke = v0;         //激发顶点
   else
      setup->vprovoke = v2;

   //2. 确定顶点的底部到最高顺序 
   {
      float y0 = v0[0][1];
      float y1 = v1[0][1];
      float y2 = v2[0][1];
      if (y0 <= y1) {
	 if (y1 <= y2) {
	    /* y0<=y1<=y2 */
	    setup->vmin = v0;         //最小的顶点
	    setup->vmid = v1;         //中间的顶点
	    setup->vmax = v2;         //最大的顶点
	 }
	 else if (y2 <= y0) {
	    /* y2<=y0<=y1 */
	    setup->vmin = v2;
	    setup->vmid = v0;
	    setup->vmax = v1;
	 }
	 else {
	    /* y0<=y2<=y1 */
	    setup->vmin = v0;
	    setup->vmid = v2;
	    setup->vmax = v1;
	 }
      }
      else {
	 if (y0 <= y2) {
	    /* y1<=y0<=y2 */
	    setup->vmin = v1;
	    setup->vmid = v0;
	    setup->vmax = v2;
	 }
	 else if (y2 <= y1) {
	    /* y2<=y1<=y0 */
	    setup->vmin = v2;
	    setup->vmid = v1;
	    setup->vmax = v0;
	 }
	 else {
	    /* y1<=y2<=y0 */
	    setup->vmin = v1;
	    setup->vmid = v2;
	    setup->vmax = v0;
	 }
      }
   }

   //3. 计算三角形的边长度
   setup->ebot.dx = setup->vmid[0][0] - setup->vmin[0][0];     //底部边的x长度
   setup->ebot.dy = setup->vmid[0][1] - setup->vmin[0][1];     //底部边的y长度
   setup->emaj.dx = setup->vmax[0][0] - setup->vmin[0][0];     //主要边的x长度
   setup->emaj.dy = setup->vmax[0][1] - setup->vmin[0][1];     //主要边的y长度
   setup->etop.dx = setup->vmax[0][0] - setup->vmid[0][0];     //顶部边的x长度
   setup->etop.dy = setup->vmax[0][1] - setup->vmid[0][1];     //顶部边的y长度

   /*
    * 计算三角形的面积。
    * 
    * 后续使用 1/面积 来计算属性的偏导数（用于计算重心坐标）。
    * 三角形的面积将与 prim->det 相同，但符号可能不同，具体取决于顶点排序的方式。
    * 为了确定图元是正面还是背面，我们使用 prim->det 的值，因为它的符号是正确的。
    */
   //4. 计算三角形的面积
   {
      const float area = (setup->emaj.dx * setup->ebot.dy -    //三角形的面积
			    setup->ebot.dx * setup->emaj.dy);

      setup->oneoverarea = 1.0f / area;      //面积的倒数

      /*
      debug_printf("%s one-over-area %f  area %f  det %f\n",
                   __FUNCTION__, setup->oneoverarea, area, det );
      */
      //如果面积为inf或nan，说明在这种情况下，三角形是不可见的
      if (util_is_inf_or_nan(setup->oneoverarea))
         return false;
   }

   /* 
    * 我们需要知道这是一个正面还是背面三角形，用于：
    * GLSL 的 gl_FrontFacing 片段属性（bool 类型）
    * 双面模板测试
    * 其中：
    * 0 = 正面
    * 1 = 背面
    */
   //5. 设置三角形的正面还是背面
   setup->facing = 
      ((det < 0.0) ^ 
       (setup->softpipe->rasterizer->front_ccw));

   {
      unsigned face = setup->facing == 0 ? PIPE_FACE_FRONT : PIPE_FACE_BACK;
      //6. 根据glCullFace剔除三角形
      if (face & setup->cull_face)
	 return false;
   }


   /* 准备像素偏移用于光栅化：
    * - 对于 GL： 像素中心为 (0.5, 0.5)。
    * - 对于其他 API： 假设像素中心为 (0.0, 0.0)。
    */
   //7. 设置像素偏移
   if (setup->softpipe->rasterizer->half_pixel_center) {
      setup->pixel_offset = 0.5f;
   } else {
      setup->pixel_offset = 0.0f;
   }

   return true;
}

/**
 * Do setup for triangle rasterization, then render the triangle.
 * 
 * \param setup -[in/out] setup上下文
 * \param v0 -[in] 三角形的第一个顶点
 * \param v1 -[in] 三角形的第二个顶点
 * \param v2 -[in] 三角形的第三个顶点
 */
void
sp_setup_tri(struct setup_context *setup,
             const float (*v0)[4],
             const float (*v1)[4],
             const float (*v2)[4])
{
   float det;           //三角形的面积
   uint layer = 0;
   unsigned viewport_index = 0;
#if DEBUG_VERTS
   printf("Setup triangle:\n");
   print_vertex(setup, v0);
   print_vertex(setup, v1);
   print_vertex(setup, v2);
#endif

   if (setup->softpipe->rasterizer->rasterizer_discard)
      return;
   
   //1. 计算三角形的面积
   det = calc_det(v0, v1, v2);      //计算三角形的面积，用于判断三角形的正面还是背面
   /*
   debug_printf("%s\n", __FUNCTION__ );
   */

#if DEBUG_FRAGS
   setup->numFragsEmitted = 0;
   setup->numFragsWritten = 0;
#endif

   //2. 三角形排序等一些操作
   if (!setup_sort_vertices( setup, det, v0, v1, v2 ))
      return;

   //3. 计算三角形的插值系数
   setup_tri_coefficients( setup );    //计算三角形的插值系数

   //4. 设置三角形的边
   setup_tri_edges( setup );           //设置三角形的边

   assert(setup->softpipe->reduced_prim == PIPE_PRIM_TRIANGLES);

   setup->span.y = 0;            //初始化span的y坐标
   setup->span.right[0] = 0;     //初始化span的右边界
   setup->span.right[1] = 0;     //初始化span的右边界

   //5. 甚至layer
   /*   setup->span.z_mode = tri_z_mode( setup->ctx ); */
   if (setup->softpipe->layer_slot > 0) {
      layer = *(unsigned *)setup->vprovoke[setup->softpipe->layer_slot];
      layer = MIN2(layer, setup->max_layer);
   }
   setup->quad[0].input.layer = layer;

   //6. 设置viewport_index
   if (setup->softpipe->viewport_index_slot > 0) {
      unsigned *udata = (unsigned*)v0[setup->softpipe->viewport_index_slot];
      viewport_index = sp_clamp_viewport_idx(*udata);
   }
   setup->quad[0].input.viewport_index = viewport_index;

   /*   init_constant_attribs( setup ); */

   if (setup->oneoverarea < 0.0) {
      /* 主边在左:
       */
      subtriangle(setup, &setup->emaj, &setup->ebot, setup->ebot.lines, viewport_index);
      subtriangle(setup, &setup->emaj, &setup->etop, setup->etop.lines, viewport_index);
   }
   else {
      /* 主边在右:
       */
      subtriangle(setup, &setup->ebot, &setup->emaj, setup->ebot.lines, viewport_index);
      subtriangle(setup, &setup->etop, &setup->emaj, setup->etop.lines, viewport_index);
   }

   flush_spans( setup );

   if (setup->softpipe->active_statistics_queries) {
      setup->softpipe->pipeline_statistics.c_primitives++;
   }

#if DEBUG_FRAGS
   printf("Tri: %u frags emitted, %u written\n",
          setup->numFragsEmitted,
          setup->numFragsWritten);
#endif
}


#define MALLOC_STRUCT(T)   (struct T *) malloc(sizeof(struct T))
#define MALLOC(_size)  malloc(_size)

void setup_test(void)
{
    int screen_w = 300;
    int screen_h = 300;

    struct setup_context setup;
    setup.facing = 0;               //面向
    setup.softpipe = MALLOC_STRUCT(softpipe_context);       //管线上下文
    struct softpipe_context *softpipe = setup.softpipe;
    softpipe->rasterizer= MALLOC_STRUCT(pipe_rasterizer_state);     //光栅化器
    struct pipe_rasterizer_state *rasterizer = softpipe->rasterizer;
    rasterizer->flatshade = 0;
    rasterizer->light_twoside = 0;
    rasterizer->clamp_vertex_color = 0;
    rasterizer->clamp_fragment_color = 0;
    rasterizer->front_ccw = 1;
    rasterizer->cull_face = 0;
    rasterizer->fill_front = 0;
    rasterizer->fill_back = 0;
    rasterizer->offset_point = 0;
    rasterizer->offset_line = 0;
    rasterizer->offset_tri = 0;
    rasterizer->scissor = 0;
    rasterizer->poly_smooth = 0;
    rasterizer->poly_stipple_enable = 0;
    rasterizer->point_smooth = 0;
    rasterizer->sprite_coord_enable = 0;
    rasterizer->point_quad_rasterization = 0;
    rasterizer->point_tri_clip = 0;
    rasterizer->point_size_per_vertex = 0;
    rasterizer->multisample = 0;
    rasterizer->force_persample_interp = 0;
    rasterizer->line_smooth = 0;
    rasterizer->line_stipple_enable = 0;
    rasterizer->line_last_pixel = 0;
    rasterizer->conservative_raster_mode = 0;
    rasterizer->flatshade_first = 0;
    rasterizer->half_pixel_center = 1;
    rasterizer->bottom_edge_rule = 1;
    rasterizer->subpixel_precision_x = 0;
    rasterizer->subpixel_precision_y = 0;
    rasterizer->rasterizer_discard = 0;
    rasterizer->tile_raster_order_fixed = 0;
    rasterizer->tile_raster_order_increasing_x = 0;
    rasterizer->tile_raster_order_increasing_y = 0;
    rasterizer->depth_clip_near = 1;
    rasterizer->depth_clip_far = 1;
    rasterizer->clip_halfz = 0;
    rasterizer->offset_units_unscaled = 0;
    rasterizer->clip_plane_enable = 0;
    rasterizer->line_stipple_factor = 0;
    rasterizer->line_stipple_pattern = 65535;
    rasterizer->sprite_coord_enable = 0;
    rasterizer->line_width = 1.0;  
    rasterizer->point_size = 1.0;
    rasterizer->offset_units = 0.0;
    rasterizer->offset_scale = 0.0;
    rasterizer->offset_clamp = 0.0;
    rasterizer->conservative_raster_dilate = 0.0;
    for(int i; i < PIPE_MAX_VIEWPORTS; i++){
        softpipe->cliprect[i].minx = 0;
        softpipe->cliprect[i].miny = 0;
        softpipe->cliprect[i].maxx = screen_w;
        softpipe->cliprect[i].maxy = screen_h;
    }

    softpipe->blend = MALLOC_STRUCT(pipe_blend_state);      //混合状态
    struct pipe_blend_state *blend = softpipe->blend;
    blend->independent_blend_enable = 0;
    blend->logicop_enable = 0;
    blend->logicop_func = 0;
    blend->dither = 1;
    blend->alpha_to_coverage = 0;
    blend->alpha_to_one = 0;
    for(int i = 0; i < PIPE_MAX_COLOR_BUFS; i++){
        if(i == 0){
            blend->rt[i].blend_enable = 0;
            blend->rt[i].rgb_func = 0;
            blend->rt[i].rgb_src_factor = 0;
            blend->rt[i].rgb_dst_factor = 0;
            blend->rt[i].alpha_func = 0;
            blend->rt[i].alpha_src_factor = 0;
            blend->rt[i].alpha_dst_factor = 0;
            blend->rt[i].colormask = 0xf;
        }else{
            blend->rt[i].blend_enable = 0;
            blend->rt[i].rgb_func = 0;
            blend->rt[i].rgb_src_factor = 0;
            blend->rt[i].rgb_dst_factor = 0;
            blend->rt[i].alpha_func = 0;
            blend->rt[i].alpha_src_factor = 0;
            blend->rt[i].alpha_dst_factor = 0;
            blend->rt[i].colormask = 0x0;
        }
    }

    softpipe->framebuffer.width = screen_w;         //framebuffer
    softpipe->framebuffer.height = screen_h;
    softpipe->framebuffer.layers = 0;
    softpipe->framebuffer.samples = 0;
    softpipe->framebuffer.nr_cbufs = PIPE_MAX_COLOR_BUFS;
    for(int i = 0; i < softpipe->framebuffer.nr_cbufs; i++){
        softpipe->framebuffer.cbufs[i] = MALLOC_STRUCT(pipe_surface);
        struct pipe_surface *surface = softpipe->framebuffer.cbufs[i];
        surface->format = PIPE_FORMAT_B8G8R8X8_UNORM;
        surface->writable = 0;
        surface->width = screen_w;
        surface->height = screen_h;
        surface->texture = MALLOC_STRUCT(pipe_resource);
        struct pipe_resource *texture = surface->texture;
        texture->width0 = screen_w;
        texture->height0 = screen_h;
        texture->array_size = 1;
        texture->format = PIPE_FORMAT_B8G8R8X8_UNORM;
        texture->target = PIPE_TEXTURE_2D;   
        texture->last_level = 0;
        texture->nr_samples = 0;
        texture->nr_storage_samples = 0;
        texture->usage = 0;
        texture->bind = 138;
        texture->flags = 0;
        texture->userdata = MALLOC(PIPE_TEXTURE_BUFFER_SIZE * PIPE_TEXTURE_BUFFER_NUMBER);
        memset(texture->userdata, 0, PIPE_TEXTURE_BUFFER_SIZE * PIPE_TEXTURE_BUFFER_NUMBER);
        surface->u.tex.level = 0;
        surface->u.tex.first_layer = 0;
        surface->u.tex.last_layer = 0;
    }

    softpipe->depth_stencil = MALLOC_STRUCT(pipe_depth_stencil_alpha_state);        //深度模板状态
    struct pipe_depth_stencil_alpha_state *depth_stencil = softpipe->depth_stencil;
    depth_stencil->depth.enabled = 0;
    depth_stencil->depth.writemask = 0;
    depth_stencil->depth.func = 0;
    depth_stencil->depth.bounds_test = 0;
    depth_stencil->depth.bounds_min = 0.0;
    depth_stencil->depth.bounds_max = 0.0;
    for(int i = 0; i < 2; i++){
        depth_stencil->stencil[i].enabled = 0;
        depth_stencil->stencil[i].func = 0;
        depth_stencil->stencil[i].fail_op = 0;
        depth_stencil->stencil[i].zfail_op = 0;
        depth_stencil->stencil[i].zpass_op = 0;
        depth_stencil->stencil[i].valuemask = 0;
        depth_stencil->stencil[i].writemask = 0;
    }
    depth_stencil->alpha.enabled = 0;
    depth_stencil->alpha.func = 0;
    depth_stencil->alpha.ref_value = 0.0;

    softpipe->fs_variant = MALLOC_STRUCT(sp_fragment_shader_variant);        //片段着色器变体
    struct sp_fragment_shader_variant *fs_variant = softpipe->fs_variant;
    fs_variant->info.num_inputs = 0;
    fs_variant->info.num_outputs = 0;
    fs_variant->info.properties[0] = 0;
    fs_variant->info.properties[1] = 0;
    fs_variant->info.properties[2] = 0;
    fs_variant->info.properties[3] = 0;
    fs_variant->info.properties[4] = 0;
    fs_variant->info.properties[5] = 1;
    fs_variant->info.properties[6] = 0;
    fs_variant->info.properties[7] = 0;
    fs_variant->info.properties[8] = 1;
    for (int i = 9; i < TGSI_PROPERTY_COUNT; i ++){
        fs_variant->info.properties[i] = 0;
    }
    for( int i = 0; i < PIPE_MAX_SHADER_INPUTS; i++){
        fs_variant->info.input_semantic_name[i] = 0;
    }
    fs_variant->info.writes_z = 0;

    softpipe->setup_info.valid = 1;     //顶点格式
    for(int i = 0; i < PIPE_MAX_SHADER_OUTPUTS; i++){
        softpipe->setup_info.attrib[i].src_index = 0;
        softpipe->setup_info.attrib[i].interp = 0;
    }
    
    softpipe->psize_slot = -1;          //点大小槽
    softpipe->reduced_prim = PIPE_PRIM_POINTS;
    softpipe->viewport_index_slot = -1;
    softpipe->layer_slot = -1;
    softpipe->early_depth = false;
    softpipe->active_statistics_queries = 0;
    softpipe->occlusion_count = 0;      //遮挡计数
    softpipe->active_query_count = 0;   //活动查询计数
    softpipe->pipeline_statistics.ia_vertices = 0;
    softpipe->pipeline_statistics.ia_primitives = 0;
    softpipe->pipeline_statistics.vs_invocations = 0;
    softpipe->pipeline_statistics.gs_invocations = 0;
    softpipe->pipeline_statistics.gs_primitives = 0;
    softpipe->pipeline_statistics.c_invocations = 0;
    softpipe->pipeline_statistics.c_primitives = 0;
    softpipe->pipeline_statistics.ps_invocations = 0;
    softpipe->pipeline_statistics.hs_invocations = 0;
    softpipe->pipeline_statistics.ds_invocations = 0;
    softpipe->pipeline_statistics.cs_invocations = 0;
    softpipe->quad.first = MALLOC_STRUCT(quad_stage);      //quad
    softpipe->quad.first->softpipe = softpipe;
    softpipe->quad.first->run = shade_quads;
    softpipe->quad.shade = MALLOC_STRUCT(quad_stage);
    softpipe->quad.shade->softpipe = softpipe;
    softpipe->quad.shade->run = shade_quads;
    softpipe->quad.depth_test = MALLOC_STRUCT(quad_stage);
    softpipe->quad.depth_test->softpipe = softpipe;
    softpipe->quad.depth_test->run = choose_depth_test;
    softpipe->quad.blend = MALLOC_STRUCT(quad_stage);
    softpipe->quad.blend->softpipe = softpipe;
    softpipe->quad.blend->run = choose_blend_quad;
    softpipe->quad.pstipple = MALLOC_STRUCT(quad_stage);
    softpipe->quad.pstipple->softpipe = softpipe;
    softpipe->quad.pstipple->run = NULL;        //stipple_quad

    for(int i = 0; i< PIPE_MAX_COLOR_BUFS; i ++){
        softpipe->cbuf_cache[i] = sp_create_tile_cache(softpipe);
        sp_tile_cache_set_surface(softpipe->cbuf_cache[i], softpipe->framebuffer.cbufs[i]);
        #if 0
        softpipe->cbuf_cache[i] = MALLOC_STRUCT(softpipe_tile_cache);
        struct softpipe_tile_cache *cbuf_cache = softpipe->cbuf_cache[i];
        cbuf_cache->pipe = softpipe;
        cbuf_cache->surface = softpipe->framebuffer.cbufs[i];
        cbuf_cache->transfer = (struct pipe_transfer **)malloc(sizeof(struct pipe_transfer *) * TRANSFER_MAX);
        for(int j = 0; j < TRANSFER_MAX; j++){
            cbuf_cache->transfer[j] = MALLOC_STRUCT(pipe_transfer);
            struct pipe_transfer *transfer = cbuf_cache->transfer[i];
            transfer->resource = cbuf_cache->surface->texture;
            transfer->level = cbuf_cache->surface->u.tex.level;
            transfer->usage = 1027;
            transfer->box.x = 0;
            transfer->box.y = 0;
            transfer->box.z = 0;
            transfer->box.width = cbuf_cache->surface->width;
            transfer->box.height = cbuf_cache->surface->height;
            transfer->box.depth = 1;
            transfer->stride = 0;
            transfer->layer_stride = 0;
        }
        cbuf_cache->transfer_map = NULL;
        cbuf_cache->num_maps = 1;
        #endif
    }

    softpipe->dirty_render_cache = 0;


    for(int i = 0; i < MAX_QUADS; i ++){
        memset(&setup.quad[i], 0, sizeof(struct quad_header));
    }

    for(int i = 0; i < PIPE_MAX_SHADER_INPUTS; i ++){
        setup.coef[i].a0[0] = 0.0;
        setup.coef[i].a0[1] = 0.0;
        setup.coef[i].a0[2] = 0.0;
        setup.coef[i].a0[3] = 0.0;
        setup.coef[i].dadx[0] = 0.0;
        setup.coef[i].dadx[1] = 0.0;
        setup.coef[i].dadx[2] = 0.0;
        setup.coef[i].dadx[3] = 0.0;
        setup.coef[i].dady[0] = 0.0;
        setup.coef[i].dady[1] = 0.0;
        setup.coef[i].dady[2] = 0.0;
        setup.coef[i].dady[3] = 0.0;
    }

    
    setup.posCoef.a0[0] = 0.0;          //插值系数
    setup.posCoef.a0[1] = 0.0;
    setup.posCoef.a0[2] = 0.0;
    setup.posCoef.a0[3] = 0.0;
    setup.posCoef.dadx[0] = 0.0;
    setup.posCoef.dadx[1] = 0.0;
    setup.posCoef.dadx[2] = 0.0;
    setup.posCoef.dadx[3] = 0.0;
    setup.posCoef.dady[0] = 0.0;
    setup.posCoef.dady[1] = 0.0;
    setup.posCoef.dady[2] = 0.0;
    setup.posCoef.dady[3] = 0.0;


    setup.max_layer = 0;
    setup.nr_vertex_attrs = 1;

    #if 0
    float v0[4][4] = {
        {75.0, 75.0, 0.5, 1.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0}
    };
    setup_point(&setup, v0);
    float v1[4][4] = {
    {225.0, 75.0, 0.5, 1.0},
    {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0}
    };
    softpipe->rasterizer->point_size = 50;
    setup_point(&setup, v1);
    float v2[4][4] = {
    {75.0, 225.0, 0.5, 1.0},
    {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0}
    };
    softpipe->rasterizer->point_size = 50;
    softpipe->rasterizer->point_smooth = true;
    setup_point(&setup, v2);
    #endif
    #if 0
    float v0[4][4] = {
        {75.0, 150.0, 0.5, 1.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0}
    };
    float v1[4][4] = {
        {225.0, 150.0, 0.5, 1.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 0.0}
    };
    setup.softpipe->reduced_prim = PIPE_PRIM_LINES;
    sp_setup_line(&setup, v0, v1);
    float v2[4][4] = {
    {150.0, 75.0, 0.5, 1.0},
    {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0}
    };
    float v3[4][4] = {
    {150.0, 225.0, 0.5, 1.0},
    {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0}
    };
    sp_setup_line(&setup, v2, v3);
    sp_setup_line(&setup, v0, v3);
    sp_setup_line(&setup, v1, v2);
    #endif
    float v0[4][4] = {
    {100.0, 100.0, 0.5, 1.0},
    {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0},
    {0.0, 0.0, 0.0, 0.0}
    };
    float v1[4][4] = {
   {200.0, 100.0, 0.5, 1.0},
   {0.0, 0.0, 0.0, 0.0},
   {0.0, 0.0, 0.0, 0.0},
   {0.0, 0.0, 0.0, 0.0}
   };
   float v2[4][4] = {
   {150.0, 200.0, 0.5, 1.0},
   {0.0, 0.0, 0.0, 0.0},
   {0.0, 0.0, 0.0, 0.0},
   {0.0, 0.0, 0.0, 0.0}
   };
   setup.softpipe->reduced_prim = PIPE_PRIM_TRIANGLES;
    sp_setup_tri(&setup, v0, v1, v2);

    softpipe_flush(softpipe, PIPE_FLUSH_END_OF_FRAME);      //刷新管线
}

