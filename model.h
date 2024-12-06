#ifndef __MODEL_H__
#define __MODEL_H__

#include "geometry.h"
#include "vec.h"

typedef struct{
	int v;
	int vt;
	int vn;
}FaceIndices;
typedef struct{
	FaceIndices vertex[3];
}Face;
typedef struct{
	Vec3f *verts_;
	Face *faces_;
}Model;

extern Model Model_init(const char *filename);
extern void Model_release(Model *model);
extern int Model_nverts(Model *model);
extern int Model_nfaces(Model *model);
extern Vec3f Model_vert(Model *model, int i);
extern Face Model_face(Model *model, int idx);


#endif //__MODEL_H__
