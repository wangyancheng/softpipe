#include <stdbool.h>

#ifndef UNUSED
#define UNUSED __attribute__((unused))
#endif

#include "quad_fs.h"
#include "quad.h"
#include "quad_blend.h"

/**
 * General depth/stencil test function.  Used when there's no fast-path.
 */
static void
depth_test_quads_fallback(UNUSED struct quad_stage *qs,
                          UNUSED struct quad_header *quads[],
                          UNUSED unsigned nr)
{

}

static void
depth_noop(struct quad_stage *qs,
           struct quad_header *quads[],
           unsigned nr)
{
   choose_blend_quad(qs, quads, nr);
}

/**
* 深度测试
* \param qs -[in] quad处理阶段的环境和方法
* \param quads -[in] quad数据
* \param nr -[in] quad个数
*/
void
choose_depth_test(struct quad_stage *qs,
                  struct quad_header *quads[],
                  unsigned nr)
{
   const struct tgsi_shader_info *fsInfo = &qs->softpipe->fs_variant->info;

   bool interp_depth = !fsInfo->writes_z || qs->softpipe->early_depth;    //深度插值

   bool alpha = qs->softpipe->depth_stencil->alpha.enabled;    //alpha测试使能

   bool depth = qs->softpipe->depth_stencil->depth.enabled;    //深度测试使能

   unsigned depthfunc = qs->softpipe->depth_stencil->depth.func;   //深度测试功能

   bool stencil = qs->softpipe->depth_stencil->stencil[0].enabled;   //模板测试使能

   bool depthwrite = qs->softpipe->depth_stencil->depth.writemask;   //允许深度buffer写位

   bool occlusion = qs->softpipe->active_query_count;     //遮挡查询使能

   bool clipped = !qs->softpipe->rasterizer->depth_clip_near;    //深度截取使能

   if(!qs->softpipe->framebuffer.zsbuf)
      depth = depthwrite = stencil = false;

   /* default */
   qs->run = depth_test_quads_fallback;

   /* look for special cases */
   if (!alpha &&
       !depth &&
       !occlusion &&
       !clipped &&
       !stencil) {
      qs->run = depth_noop;
   }
   else if (!alpha &&
            interp_depth &&
            depth &&
            depthwrite &&
            !occlusion &&
            !clipped &&
            !stencil)
   {
      if (1) {
         switch (depthfunc) {
         case PIPE_FUNC_NEVER:
            qs->run = depth_test_quads_fallback;
            break;
         case PIPE_FUNC_LESS:
            //qs->run = depth_interp_z16_less_write;
            break;
         case PIPE_FUNC_EQUAL:
            //qs->run = depth_interp_z16_equal_write;
            break;
         case PIPE_FUNC_LEQUAL:
            //qs->run = depth_interp_z16_lequal_write;
            break;
         case PIPE_FUNC_GREATER:
            //qs->run = depth_interp_z16_greater_write;
            break;
         case PIPE_FUNC_NOTEQUAL:
            //qs->run = depth_interp_z16_notequal_write;
            break;
         case PIPE_FUNC_GEQUAL:
            //qs->run = depth_interp_z16_gequal_write;
            break;
         case PIPE_FUNC_ALWAYS:
            //qs->run = depth_interp_z16_always_write;
            break;
         default:
            qs->run = depth_test_quads_fallback;
            break;
         }
      }
   }

   /* next quad/fragment stage */
   qs->run( qs, quads, nr );
}