#ifndef __QUAD_BLEND_H__
#define __QUAD_BLEND_H__

#include "quad_fs.h"
#include "quad.h"

enum format
{
   RGBA,
   RGB,
   LUMINANCE,
   LUMINANCE_ALPHA,
   INTENSITY
};

enum util_format_type {
   UTIL_FORMAT_TYPE_VOID = 0,		
   UTIL_FORMAT_TYPE_UNSIGNED = 1,	
   UTIL_FORMAT_TYPE_SIGNED = 2,		
   UTIL_FORMAT_TYPE_FIXED = 3,		
   UTIL_FORMAT_TYPE_FLOAT = 4		
};

/** Subclass of quad_stage */
struct blend_quad_stage
{
   struct quad_stage base;
   bool clamp[PIPE_MAX_COLOR_BUFS];  /**< clamp colors to [0,1]? */
   enum format base_format[PIPE_MAX_COLOR_BUFS];      //颜色分量
   enum util_format_type format_type[PIPE_MAX_COLOR_BUFS];     //颜色类型
};


extern void
choose_blend_quad(struct quad_stage *qs,
                  struct quad_header *quads[],
                  unsigned nr);

#endif