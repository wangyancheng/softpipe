#include <string.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "vec.h"
#include "model.h"

Model Model_init(const char *filename){
    Model model;
    model.faces_ = vector_create();
    model.verts_ = vector_create();

    FILE *fp = fopen(filename, "r");
    if(fp == NULL){
        printf("can't open file %s\n", filename);
        return model;
    }
    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        if (!strncmp(line, "v ", 2)) {
            Vec3f v;
            sscanf(line+2, "%f %f %f", &v.raw[0],&v.raw[1],&v.raw[2]); 
            vector_add(&model.verts_, v);
        } else if (!strncmp(line, "f ", 2)) {
            unsigned index;
            char c0, c1;
            index = 2;
            char str[255];
            unsigned strlen = 0;
            int indices;
            int vertex=0;
            int which=0;   //0:V,1:VT,2:VN
            Face f;
            c0 = line[index];
            while (c0 != '\n' && c0 != '\r' && index < 200){
                strlen = 0;
                memset(str, 0, sizeof(str));
                c0 = line[index];
                while(c0 >= '0' && c0 <= '9'){
                    if(strlen > 200)
                        break;
                    str[strlen++] =  c0;
                    c0 = line[++index];
                }
                indices = atoi(str) - 1;
                if(which == 0){
                    f.vertex[vertex].v = indices;
                }
                else if(which == 1){
                    f.vertex[vertex].vt = indices;
                }
                else if(which == 2){
                    f.vertex[vertex].vn = indices;
                }
                c0 = line[index++];
                c1 = line[index];
                if(c0 == ' '){
                    which=0;        //V
                    vertex++;
                    if(vertex>3)
                        break;
                    continue;
                }
                else if(c0 == '/' && c1 != '/'){
                    which++;    //VT
                    continue;
                }
                else if(c0 == '/' && c1 == '/'){
                    index++;
                    which += 2;    //VN
                    continue;
                }
            }
            vector_add(&model.faces_, f);
        }
    }
    fclose(fp);
    printf("v:%u f:%u \n", (unsigned)vector_size(model.verts_), (unsigned)vector_size(model.faces_));
    return model;
}

void Model_release(Model *model) {
    vector_free(model->verts_);
    vector_free(model->faces_);
}

int Model_nverts(Model *model) {
    return (int)vector_size(model->verts_);
}

int Model_nfaces(Model *model) {
    return (int)vector_size(model->faces_);
}

Face Model_face(Model *model, int idx) {
    return model->faces_[idx];
}

Vec3f Model_vert(Model *model, int i) {
    return model->verts_[i];
}

