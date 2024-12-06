#include "tgaimage.h"
#include "raster.h"
#include "model.h"


int main(int argc, char** argv) {
	Model model;
	if (2==argc) {
		model = Model_init(argv[1]);
	} else {
		model = Model_init("obj/african_head.obj");
	}
	const int width  = 800;
	const int height = 800;
	TGAImage image = tgai_create(width, height, RGB);
	TGAColor white = tgacolor_rgba(255, 255, 255, 255);
	//TGAColor red = tgacolor_rgba(255, 0, 0, 255);

	for (int i=0; i<Model_nfaces(&model); i++) {
    	Face face = Model_face(&model, i);
    	for (int j=0; j<3; j++) {
        	Vec3f v0 = Model_vert(&model, face.vertex[j].v);
        	Vec3f v1 = Model_vert(&model, face.vertex[(j+1)%3].v);
        	int x0 = (v0.x+1.)*width/2.;
        	int y0 = (v0.y+1.)*height/2.;
        	int x1 = (v1.x+1.)*width/2.;
        	int y1 = (v1.y+1.)*height/2.;
        	rs_line(x0, y0, x1, y1, &image, white);
    	}
	}
	Model_release(&model);
	tgai_flip_vertically(&image); // i want to have the origin at the left bottom corner of the image
	tgai_write_tga_file(&image, "output.tga", true);
	return 0;
}

