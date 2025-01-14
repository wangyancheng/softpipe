#ifndef __CONTEXT_H__
#define __CONTEXT_H__

#include <stdint.h>
#include "state.h"
#include "setup.h"
#include "quad_fs.h"
#include "tile_cache.h"

typedef unsigned char ubyte;

enum tgsi_property_name {
   TGSI_PROPERTY_GS_INPUT_PRIM,
   TGSI_PROPERTY_GS_OUTPUT_PRIM,
   TGSI_PROPERTY_GS_MAX_OUTPUT_VERTICES,
   TGSI_PROPERTY_FS_COORD_ORIGIN,
   TGSI_PROPERTY_FS_COORD_PIXEL_CENTER,
   TGSI_PROPERTY_FS_COLOR0_WRITES_ALL_CBUFS,
   TGSI_PROPERTY_FS_DEPTH_LAYOUT,
   TGSI_PROPERTY_VS_PROHIBIT_UCPS,
   TGSI_PROPERTY_GS_INVOCATIONS,
   TGSI_PROPERTY_VS_WINDOW_SPACE_POSITION,
   TGSI_PROPERTY_TCS_VERTICES_OUT,
   TGSI_PROPERTY_TES_PRIM_MODE,
   TGSI_PROPERTY_TES_SPACING,
   TGSI_PROPERTY_TES_VERTEX_ORDER_CW,
   TGSI_PROPERTY_TES_POINT_MODE,
   TGSI_PROPERTY_NUM_CLIPDIST_ENABLED,
   TGSI_PROPERTY_NUM_CULLDIST_ENABLED,
   TGSI_PROPERTY_FS_EARLY_DEPTH_STENCIL,
   TGSI_PROPERTY_FS_POST_DEPTH_COVERAGE,
   TGSI_PROPERTY_NEXT_SHADER,
   TGSI_PROPERTY_CS_FIXED_BLOCK_WIDTH,
   TGSI_PROPERTY_CS_FIXED_BLOCK_HEIGHT,
   TGSI_PROPERTY_CS_FIXED_BLOCK_DEPTH,
   TGSI_PROPERTY_MUL_ZERO_WINS,
   TGSI_PROPERTY_COUNT,
};

enum tgsi_semantic {
   TGSI_SEMANTIC_POSITION,
   TGSI_SEMANTIC_COLOR,
   TGSI_SEMANTIC_BCOLOR,       /**< back-face color */
   TGSI_SEMANTIC_FOG,
   TGSI_SEMANTIC_PSIZE,
   TGSI_SEMANTIC_GENERIC,
   TGSI_SEMANTIC_NORMAL,
   TGSI_SEMANTIC_FACE,
   TGSI_SEMANTIC_EDGEFLAG,
   TGSI_SEMANTIC_PRIMID,
   TGSI_SEMANTIC_INSTANCEID,  /**< doesn't include start_instance */
   TGSI_SEMANTIC_VERTEXID,
   TGSI_SEMANTIC_STENCIL,
   TGSI_SEMANTIC_CLIPDIST,
   TGSI_SEMANTIC_CLIPVERTEX,
   TGSI_SEMANTIC_GRID_SIZE,   /**< grid size in blocks */
   TGSI_SEMANTIC_BLOCK_ID,    /**< id of the current block */
   TGSI_SEMANTIC_BLOCK_SIZE,  /**< block size in threads */
   TGSI_SEMANTIC_THREAD_ID,   /**< block-relative id of the current thread */
   TGSI_SEMANTIC_TEXCOORD,    /**< texture or sprite coordinates */
   TGSI_SEMANTIC_PCOORD,      /**< point sprite coordinate */
   TGSI_SEMANTIC_VIEWPORT_INDEX,  /**< viewport index */
   TGSI_SEMANTIC_LAYER,       /**< layer (rendertarget index) */
   TGSI_SEMANTIC_SAMPLEID,
   TGSI_SEMANTIC_SAMPLEPOS,
   TGSI_SEMANTIC_SAMPLEMASK,
   TGSI_SEMANTIC_INVOCATIONID,
   TGSI_SEMANTIC_VERTEXID_NOBASE,
   TGSI_SEMANTIC_BASEVERTEX,
   TGSI_SEMANTIC_PATCH,       /**< generic per-patch semantic */
   TGSI_SEMANTIC_TESSCOORD,   /**< coordinate being processed by tess */
   TGSI_SEMANTIC_TESSOUTER,   /**< outer tessellation levels */
   TGSI_SEMANTIC_TESSINNER,   /**< inner tessellation levels */
   TGSI_SEMANTIC_VERTICESIN,  /**< number of input vertices */
   TGSI_SEMANTIC_HELPER_INVOCATION,  /**< current invocation is helper */
   TGSI_SEMANTIC_BASEINSTANCE,
   TGSI_SEMANTIC_DRAWID,
   TGSI_SEMANTIC_WORK_DIM,    /**< opencl get_work_dim value */
   TGSI_SEMANTIC_SUBGROUP_SIZE,
   TGSI_SEMANTIC_SUBGROUP_INVOCATION,
   TGSI_SEMANTIC_SUBGROUP_EQ_MASK,
   TGSI_SEMANTIC_SUBGROUP_GE_MASK,
   TGSI_SEMANTIC_SUBGROUP_GT_MASK,
   TGSI_SEMANTIC_SUBGROUP_LE_MASK,
   TGSI_SEMANTIC_SUBGROUP_LT_MASK,
   TGSI_SEMANTIC_COUNT,       /**< number of semantic values */
};

/**
 * Shader摘要信息
 */
struct tgsi_shader_info
{
    ubyte num_inputs;       //输入槽(slot)数量
    ubyte num_outputs;      //输出槽(slot)数量
    unsigned properties[TGSI_PROPERTY_COUNT]; /* index with TGSI_PROPERTY_ */
    ubyte input_semantic_name[PIPE_MAX_SHADER_INPUTS]; /**< TGSI_SEMANTIC_x */
    bool writes_z;      // does fragment shader write Z value? 
    ubyte input_cylindrical_wrap[PIPE_MAX_SHADER_INPUTS];      //输入槽的圆柱包裹
};

/**
 * fragment shader变量
 */
struct sp_fragment_shader_variant
{
    struct tgsi_shader_info info;       //Shader信息
};

/**
 * setup信息
 */
struct sp_setup_info {
   unsigned valid;            //有效位
   struct {
      unsigned interp:8;      //插值类型 < SP_INTERP_X >
      int src_index:8;        //槽索引
   } attrib[PIPE_MAX_SHADER_OUTPUTS];
};

/**
 * Query result for PIPE_QUERY_PIPELINE_STATISTICS.
 */
struct pipe_query_data_pipeline_statistics
{
   uint64_t ia_vertices;    /**< Num vertices read by the vertex fetcher. */
   uint64_t ia_primitives;  /**< Num primitives read by the vertex fetcher. */
   uint64_t vs_invocations; /**< Num vertex shader invocations. */
   uint64_t gs_invocations; /**< Num geometry shader invocations. */
   uint64_t gs_primitives;  /**< Num primitives output by a geometry shader. */
   uint64_t c_invocations;  /**< Num primitives sent to the rasterizer. */
   uint64_t c_primitives;   /**< Num primitives that were rendered. */
   uint64_t ps_invocations; /**< Num pixel shader invocations. */
   uint64_t hs_invocations; /**< Num hull shader invocations. */
   uint64_t ds_invocations; /**< Num domain shader invocations. */
   uint64_t cs_invocations; /**< Num compute shader invocations. */
};




/**
 * 上下文
 */
struct softpipe_context {
   struct pipe_rasterizer_state *rasterizer;
   struct pipe_scissor_state cliprect[PIPE_MAX_VIEWPORTS];
   struct pipe_blend_state *blend;     //混合
   struct pipe_framebuffer_state framebuffer;		//framebuffer
   struct pipe_depth_stencil_alpha_state *depth_stencil;
   struct sp_fragment_shader_variant *fs_variant;  //Fragment shader变量
   struct sp_setup_info setup_info;       //顶点格式
   int8_t psize_slot;          //vertex shader输出槽（output slot）包含点大小（point size)
   unsigned reduced_prim;           //精简后的基元类型
   int8_t viewport_index_slot;     //vertex shader 输出槽包含视口索引
   int8_t layer_slot;              //vertex shader 输出槽包含层索引
   bool early_depth;                //是否启用早期深度测试
   unsigned active_statistics_queries;    //激活静态查询
   /* Counter for occlusion queries.  Note this supports overlapping
    * queries.
    */
   uint64_t occlusion_count;
   unsigned active_query_count;        //遮挡查询使能
   struct pipe_query_data_pipeline_statistics pipeline_statistics;   //查询结果
    /** 软件quad渲染管线*/
   struct {
      struct quad_stage *shade;
      struct quad_stage *depth_test;
      struct quad_stage *blend;
      struct quad_stage *pstipple;
      struct quad_stage *first; /**< points to one of the above stages */
   } quad;
   struct softpipe_tile_cache *cbuf_cache[PIPE_MAX_COLOR_BUFS];		//tile cache
   struct softpipe_tile_cache *zsbuf_cache;

   bool dirty_render_cache;

};

#endif