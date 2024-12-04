#include "tgaimage.h"


int main(int argc, char** argv) {
	//tgai_test();
	//const TGAColor white = tgacolor_rgba(255, 255, 255, 255);
	const TGAColor red   = tgacolor_rgba(255, 0,   0,   255);
	TGAImage image = tgai_create(100, 100, RGB);
	tgai_set(image, 52, 41, red);
	tgai_flip_vertically(image); // i want to have the origin at the left bottom corner of the image
	tgai_write_tga_file(image, "output.tga", true);
	return 0;
}

