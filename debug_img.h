#ifndef __DEBUG_IMG_H__
#define __DEBUG_IMG_H__

#include "context.h"
#include "state.h"

extern void
debug_dump_surface_bmp(struct softpipe_context *pipe,
                       const char *filename,
                       struct pipe_surface *surface);


#endif