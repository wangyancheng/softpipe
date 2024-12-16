#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "tgaimage.h"
#include "geometry.h"

void swapi(int *a, int *b){
    int t;
    t = *a;
    *a = *b;
    *b = t;
}

void swapv2i(Vec2i *a, Vec2i *b){
    Vec2i t;
    t.x = a->x;
    t.y = a->y;
    a->x = b->x;
    a->y = b->y;
    b->x = t.x;
    b->y = t.y;
}

int max(int a, int b){
    return a > b? a : b;
}

int min(int a, int b){
    return a < b? a: b;
}

void rs_line(int x0, int y0, int x1, int y1, TGAImage *image, TGAColor color) {
    bool steep = false;
    if (abs(x0-x1)<abs(y0-y1)) {
        swapi(&x0, &y0);
        swapi(&x1, &y1);
        steep = true;
    }
    if (x0>x1) {
        swapi(&x0, &x1);
        swapi(&y0, &y1);
    }
    int dx = x1-x0;
    int dy = y1-y0;
    int derror2 = abs(dy)*2;
    int error2 = 0;
    int y = y0;
    for (int x=x0; x<=x1; x++) {
        if (steep) {
            tgai_set(image, y, x, color);
        } else {
            tgai_set(image, x, y, color);
        }
        error2 += derror2;

        if (error2>dx) {
            y += (y1>y0?1:-1);
            error2 -= dx*2;
        }
    }
}
#if !RASTERMODE
void rs_triangle(Vec2i *pts, TGAImage *image, TGAColor color) {
    if (pts[0].y==pts[1].y && pts[0].y==pts[2].y) return; // i dont care about degenerate triangles
    if (pts[0].y>pts[1].y) swapv2i(&pts[0], &pts[1]);
    if (pts[0].y>pts[2].y) swapv2i(&pts[0], &pts[2]);
    if (pts[1].y>pts[2].y) swapv2i(&pts[1], &pts[2]);
    int total_height = pts[2].y-pts[0].y;
    for (int i=0; i<total_height; i++) {
        bool second_half = i>pts[1].y-pts[0].y || pts[1].y==pts[0].y;
        int segment_height = second_half ? pts[2].y-pts[1].y : pts[1].y-pts[0].y;
        float alpha = (float)i/total_height;
        float beta  = (float)(i-(second_half ? pts[1].y-pts[0].y : 0))/segment_height; // be careful: with above conditions no division by zero here
        Vec2i A;
        A.x = pts[0].x + (pts[2].x-pts[0].x)*alpha;
        A.y = pts[0].y + (pts[2].y-pts[0].y)*alpha;
        Vec2i B;
        B.x = second_half ? pts[1].x + (pts[2].x-pts[1].x)*beta : pts[0].x + (pts[1].x-pts[0].x)*beta;
        B.y = second_half ? pts[1].y + (pts[2].y-pts[1].y)*beta : pts[0].y + (pts[1].y-pts[0].y)*beta;
        if (A.x>B.x) swapv2i(&A, &B);
        for (int j=A.x; j<=B.x; j++) {
            tgai_set(image, j, pts[0].y+i, color); // attention, due to int casts t0.y+i != A.y
        }
    }
}
#endif

#if RASTERMODE
Vec3f barycentric(Vec2i *pts, Vec2i P) {
    Vec3f u1 = Vec3f_init(pts[2].y-pts[0].y, pts[1].y-pts[0].y, pts[0].y-P.y);
    Vec3f u2 = Vec3f_init(pts[2].x-pts[0].x, pts[1].x-pts[0].x, pts[0].x-P.x);
    Vec3f u = Vec3f_cross(u2, u1);
    /* `pts` and `P` has integer value as coordinates
       so `abs(u[2])` < 1 means `u[2]` is 0, that means
       triangle is degenerate, in this case return something with negative coordinates */
    Vec3f r = Vec3f_init(-1, 1, 1);
    if (abs(u.z)<1) return r;
    r = Vec3f_init(1.f-(u.x+u.y)/u.z, u.y/u.z, u.x/u.z);
    return r; 
} 
 
void rs_triangle(Vec2i *pts, TGAImage *image, TGAColor color) { 
    Vec2i bboxmin = Vec2i_init(tgai_get_width(image)-1, tgai_get_height(image)-1);
    Vec2i bboxmax = Vec2i_init(0, 0);
    Vec2i clamp = Vec2i_init(tgai_get_width(image)-1, tgai_get_height(image)-1);
    for (int i=0; i<3; i++) { 
        bboxmin.x = max(0, min(bboxmin.x, pts[i].x));
	bboxmin.y = max(0, min(bboxmin.y, pts[i].y));
	bboxmax.x = min(clamp.x, max(bboxmax.x, pts[i].x));
	bboxmax.y = min(clamp.y, max(bboxmax.y, pts[i].y));
    } 
    Vec2i P; 
    for (P.x=bboxmin.x; P.x<=bboxmax.x; P.x++) { 
        for (P.y=bboxmin.y; P.y<=bboxmax.y; P.y++) { 
            Vec3f bc_screen  = barycentric(pts, P); 
            if (bc_screen.x<0 || bc_screen.y<0 || bc_screen.z<0) continue; 
            tgai_set(image, P.x, P.y, color); 
        } 
    } 
}
#endif