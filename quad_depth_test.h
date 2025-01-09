#ifndef __QUAD_DEPTH_TEST_H__
#define __QUAD_DEPTH_TEST_H__

#include "quad_fs.h"
#include "quad.h"

extern void choose_depth_test(struct quad_stage *qs, struct quad_header *quads[], unsigned nr);

#endif