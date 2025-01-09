#include <stdbool.h>
#include <assert.h>
#include "setup.h"
#include "quad_depth_test.h"

/**
* 为 quad 中的四个片段执行片段着色器。
*
* \param qs -[in] quad处理阶段
* \param quad -[in] quad数据
*/
static inline bool
shade_quad(struct quad_stage *qs, struct quad_header *quad)
{
   struct softpipe_context *softpipe = qs->softpipe;

   //1. 统计像素着色器调用的次数
   if (softpipe->active_statistics_queries) {
      softpipe->pipeline_statistics.ps_invocations +=
         util_bitcount(quad->inout.mask);
   }

   /* run shader */
   //2. 运行fragement shader
   if(1){
      for(int i = 0; i < PIPE_MAX_COLOR_BUFS; i++){
         for(int j = 0; j < TGSI_QUAD_SIZE; j++){
            quad->output.color[i][0][j] = 1;
            quad->output.color[i][1][j] = 0;
            quad->output.color[i][2][j] = 0;
            quad->output.color[i][3][j] = 1;
         }
      }
   }
    return true;
}

/**
 * quad的coverage计算
 * 
 * \param qs -[in] quad处理阶段
 * \param quad -[in] quad数据
 */
static void
coverage_quad(struct quad_stage *qs, struct quad_header *quad)
{
   struct softpipe_context *softpipe = qs->softpipe;
   unsigned int cbuf;

   /* loop over colorbuffer outputs */
   for (cbuf = 0; cbuf < softpipe->framebuffer.nr_cbufs; cbuf++) {
      float (*quadColor)[4] = quad->output.color[cbuf];
      unsigned j;
      for (j = 0; j < TGSI_QUAD_SIZE; j++) {
         assert(quad->input.coverage[j] >= 0.0);
         assert(quad->input.coverage[j] <= 1.0);
         quadColor[3][j] *= quad->input.coverage[j];
      }
   }
}


/**
 * 在quad中，把位置渲染为颜色、深度、模板值
 *
 * \param qs -[in] quad处理阶段
 * \param quads -[in] quads数据
 * \param nr -[in] quad的数量
 */
void
shade_quads(struct quad_stage *qs,
            struct quad_header *quads[],
            unsigned nr)
{
    unsigned i, nr_quads = 0;

    for (i = 0; i < nr; i++) {
        /**
        * 只有当所有碎片都被杀死_并且它不是输出列表中的第一个四边形时，
        * 才会从输出列表中省略该四边形。在（优化的）深度测试代码中，
        * 第一个四元组是特殊的：四元组的 Z 坐标是相对于列表中的第一个四元组逐步插值的。
        * 对于多通道算法，我们需要在每个通道中生成完全相同的 Z 值。
        * 如果插值从不同的四边形开始，我们可以得到相同 (x,y) 的不同 Z 值。
        */
        //1. 单个quad处理
        if (!shade_quad(qs, quads[i]) && i > 0)
            continue; /* quad totally culled/killed */

	    //2. coverage处理
        if (/*do_coverage*/ 0)
            coverage_quad(qs, quads[i] );

        quads[nr_quads++] = quads[i];
    }

    //3. 运行片段测试
    if (nr_quads){
            choose_depth_test(qs, quads, nr_quads);
    }
        
    
}