#ifndef __QUAD_FS_H__
#define __QUAD_FS_H__

#include "context.h"
#include "quad.h"

struct quad_stage {
   struct softpipe_context *softpipe;

   /** the stage action */
   void (*run)(struct quad_stage *qs, struct quad_header *quad[], unsigned nr);
};

extern void shade_quads(struct quad_stage *qs, struct quad_header *quads[], unsigned nr);
#endif