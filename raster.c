#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "tgaimage.h"

void swap(int *a, int *b){
    int t;
    t = *a;
    *a = *b;
    *b = t;
}

void rs_line(int x0, int y0, int x1, int y1, TGAImage *image, TGAColor color) {
    bool steep = false;
    if (abs(x0-x1)<abs(y0-y1)) {
        swap(&x0, &y0);
        swap(&x1, &y1);
        steep = true;
    }
    if (x0>x1) {
        swap(&x0, &x1);
        swap(&y0, &y1);
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