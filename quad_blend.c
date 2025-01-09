#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include "quad.h"
#include "quad_fs.h"
#include "state.h"
#include "tile_cache.h"
#include "quad_blend.h"

/** cast wrapper */
static inline struct blend_quad_stage *
blend_quad_stage(struct quad_stage *stage)
{
   return (struct blend_quad_stage *) stage;
}

/**
 * Clamp X to [MIN, MAX].
 * This is a macro to allow float, int, uint, etc. types.
 * We arbitrarily turn NaN into MIN.
 */
#define CLAMP( X, MIN, MAX )  ( (X)>(MIN) ? ((X)>(MAX) ? (MAX) : (X)) : (MIN) )

/**
 * 将颜色重新基准化到[0,1]范围内
 */
static void
clamp_colors(float (*quadColor)[4])
{
   unsigned i, j;

   for (i = 0; i < 4; i++) {
      for (j = 0; j < TGSI_QUAD_SIZE; j++) {
         quadColor[i][j] = CLAMP(quadColor[i][j], 0.0F, 1.0F);
      }
   }
}

/**
 * If we're drawing to a luminance, luminance/alpha or intensity surface
 * we have to adjust (rebase) the fragment/quad colors before writing them
 * to the tile cache.  The tile cache always stores RGBA colors but if
 * we're caching a L/A surface (for example) we need to be sure that R=G=B
 * so that subsequent reads from the surface cache appear to return L/A
 * values.
 * The piglit fbo-blending-formats test will exercise this.
 */
static void
rebase_colors(enum format base_format, float (*quadColor)[4])
{
   unsigned i;

   switch (base_format) {
   case RGB:
      for (i = 0; i < 4; i++) {
         /* A = 1 */
         quadColor[3][i] = 1.0F;
      }
      break;
   case LUMINANCE:
      for (i = 0; i < 4; i++) {
         /* B = G = R */
         quadColor[2][i] = quadColor[1][i] = quadColor[0][i];
         /* A = 1 */
         quadColor[3][i] = 1.0F;
      }
      break;
   case LUMINANCE_ALPHA:
      for (i = 0; i < 4; i++) {
         /* B = G = R */
         quadColor[2][i] = quadColor[1][i] = quadColor[0][i];
      }
      break;
   case INTENSITY:
      for (i = 0; i < 4; i++) {
         /* A = B = G = R */
         quadColor[3][i] = quadColor[2][i] = quadColor[1][i] = quadColor[0][i];
      }
      break;
   default:
      ; /* nothing */
   }
}

/**
* 仅仅将quad颜色复制到帧缓冲区图块（尊重写入掩码），
* 如果需要，在写入颜色后，Clamping将会处理
*
* \param qs -[in] quad的状态及方法
* \param quads -[in] quad数据
* \param nr -[in] quad个数
*/
static void
single_output_color(struct quad_stage *qs,
                    struct quad_header *quads[],
                    unsigned nr)
{
   struct blend_quad_stage *bqs = blend_quad_stage(qs);
   unsigned int i, j, q;

   //1. 获得图块
   struct softpipe_cached_tile *tile
      = sp_get_cached_tile(qs->softpipe->cbuf_cache[0],
                           quads[0]->input.x0,
                           quads[0]->input.y0, quads[0]->input.layer);

   for (q = 0; q < nr; q++) {
      struct quad_header *quad = quads[q];
      float (*quadColor)[4] = quad->output.color[0];        //颜色
      const int itx = (quad->input.x0 & (TILE_SIZE-1));     //tile的x坐标
      const int ity = (quad->input.y0 & (TILE_SIZE-1));     //tile的y坐标

      //2. clamp color
      if (qs->softpipe->rasterizer->clamp_fragment_color)
         clamp_colors(quadColor);

      //3. 重新设置颜色
      rebase_colors(bqs->base_format[0], quadColor);

      //4. 把颜色拷贝到tile
      for (j = 0; j < TGSI_QUAD_SIZE; j++) {
         if (quad->inout.mask & (1 << j)) {
            int x = itx + (j & 1);
            int y = ity + (j >> 1);
            for (i = 0; i < 4; i++) { /* loop over color chans */
               tile->data.color[y][x][i] = quadColor[i][j];
            }
         }
      }
      if(1){
         union tile_address addr = tile_address(quads[0]->input.x0, 
                                    quads[0]->input.y0, quads[0]->input.layer);
         const int pos = CACHE_POS(addr.bits.x,
                          addr.bits.y, addr.bits.layer);
         char filename[100];
         sprintf(filename, "%d.tile", pos);
         FILE *fp = fopen(filename, "w");
         for(int y = 0; y < TILE_SIZE; y++){
            for(int x = 0; x < TILE_SIZE; x++){
               fprintf(fp, "%f,%f,%f,%f;",tile->data.color[y][x][0], tile->data.color[y][x][1], tile->data.color[y][x][2], tile->data.color[y][x][3]);
            }
            fprintf(fp, "\n");
         }
         fclose(fp);
      }
   }
}

#define MALLOC_STRUCT(T)   (struct T *) malloc(sizeof(struct T))
#define MALLOC(_size)  malloc(_size)
/**
* 混合测试
* \param qs -[in] quad处理的状态和方法
* \param quads -[in] quad数据
* \param nr -[in] quad个数
*/
void
choose_blend_quad(struct quad_stage *qs,
                  struct quad_header *quads[],
                  unsigned nr)
{
   //struct blend_quad_stage *bqs = blend_quad_stage(qs);
   struct blend_quad_stage *bqs = MALLOC_STRUCT(blend_quad_stage);
   bqs->base = *qs;
   struct softpipe_context *softpipe = qs->softpipe;
   const struct pipe_blend_state *blend = softpipe->blend;
   unsigned i;

   //qs->run = blend_fallback;

   if (softpipe->framebuffer.nr_cbufs == 0) {
      //qs->run = blend_noop;
   }
   else if (!softpipe->blend->logicop_enable &&
            softpipe->blend->rt[0].colormask == 0xf &&
            softpipe->framebuffer.nr_cbufs == 1)
   {
      if (softpipe->framebuffer.cbufs[0] == NULL) {
         //qs->run = blend_noop;
      }
      else if (!blend->rt[0].blend_enable) {
         qs->run = single_output_color;
      }
      else if (blend->rt[0].rgb_src_factor == blend->rt[0].alpha_src_factor &&
               blend->rt[0].rgb_dst_factor == blend->rt[0].alpha_dst_factor &&
               blend->rt[0].rgb_func == blend->rt[0].alpha_func)
      {
         if (blend->rt[0].alpha_func == PIPE_BLEND_ADD) {
            if (blend->rt[0].rgb_src_factor == PIPE_BLENDFACTOR_ONE &&
                blend->rt[0].rgb_dst_factor == PIPE_BLENDFACTOR_ONE) {
               //qs->run = blend_single_add_one_one;
            }
            else if (blend->rt[0].rgb_src_factor == PIPE_BLENDFACTOR_SRC_ALPHA &&
                blend->rt[0].rgb_dst_factor == PIPE_BLENDFACTOR_INV_SRC_ALPHA){
               //qs->run = blend_single_add_src_alpha_inv_src_alpha;
            }

         }
      }
   }

   /**
    * 对于每个颜色缓冲区，确定缓冲区是否具有目标 alpha 以及是否需要颜色钳位。
    */
   for (i = 0; i < softpipe->framebuffer.nr_cbufs; i++) {
      if (softpipe->framebuffer.cbufs[i]) {
         #if 0
         const enum pipe_format format = softpipe->framebuffer.cbufs[i]->format;
         const struct util_format_description *desc =
            util_format_description(format);
         // assuming all or no color channels are normalized: 
         bqs->clamp[i] = desc->channel[0].normalized;
         bqs->format_type[i] = desc->channel[0].type;

         if (util_format_is_intensity(format))
            bqs->base_format[i] = INTENSITY;
         else if (util_format_is_luminance(format))
            bqs->base_format[i] = LUMINANCE;
         else if (util_format_is_luminance_alpha(format))
            bqs->base_format[i] = LUMINANCE_ALPHA;
         else if (!util_format_has_alpha(format))
            bqs->base_format[i] = RGB;
         else
            bqs->base_format[i] = RGBA;
         
         #endif
         bqs->clamp[i] = true;
         bqs->base_format[i] = RGB;
         bqs->format_type[i] = UTIL_FORMAT_TYPE_UNSIGNED;
      }
   }

   qs->run(&(bqs->base), quads, nr);
   free(bqs);
}
