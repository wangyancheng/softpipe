#ifndef __IMAGE_H__
#define __IMAGE_H__

#include <stdbool.h>

typedef struct{
	union {
		struct {
			unsigned char b, g, r, a;
		};
		unsigned char raw[4]; 
		unsigned int val;
	};
	int bytespp;
}TGAColor;

extern TGAColor tgacolor_rgba(unsigned char R, unsigned char G, unsigned char B, unsigned char A);
extern TGAColor tgacolor_val(int v, int bpp);
extern TGAColor tgacolor_raw(const unsigned char *p, int bpp);

enum Format {
	GRAYSCALE=1, RGB=3, RGBA=4
};
#pragma pack(push,1)
typedef struct{
	char idlength;           //Length of the image ID field
	char colormaptype;		//Whether a color map is included
	char datatypecode;		//Image type. 10: run-length encoded true-color image
	short colormaporigin;   //index of first color map entry that is included in the file 
	short colormaplength;	//number of entries of the color map that are included in the file
	char colormapdepth;		//number of bits per color map entry
	short x_origin;			//absolute coordinate of lower-left corner for displays where origin is at the lower left
	short y_origin;			//as for X-origin
	short width;			//width in pixels
	short height;			//height in pixels
	char  bitsperpixel;		//bits per pixel
	char  imagedescriptor;	//bits 3–0 give the alpha channel depth, bits 5–4 give pixel ordering
							/*Bit 4 of the image descriptor byte indicates right-to-left pixel ordering if set. 
							Bit 5 indicates an ordering of top-to-bottom. Otherwise, pixels are stored in bottom-to-top, left-to-right order.*/
}TGA_Header;
#pragma pack(pop)
typedef struct{
	unsigned char* data;
	int width;
	int height;
	int bytespp;
}TGAImage;
extern TGAImage tgai_create(int w, int h, int bpp);
extern void tgai_release(TGAImage img);
extern void tgai_copy(TGAImage destin, TGAImage source);
extern int tgai_get_bytespp(TGAImage img);
extern int tgai_get_width(TGAImage img);
extern int tgai_get_height(TGAImage img);
extern unsigned char *tgai_buffer(TGAImage img);
extern void tgai_clear(TGAImage img);
//extern bool tgai_scale(TGAImage img, int w, int h);
extern TGAColor tgai_get(TGAImage img, int x, int y);
extern bool tgai_set(TGAImage img, int x, int y, TGAColor c);
extern bool tgai_flip_vertically(TGAImage img);
extern bool tgai_flip_horizontally(TGAImage img);
extern bool tgai_read_file(TGAImage img, const char *filename);
extern bool tgai_write_tga_file(TGAImage img, const char *filename, bool rle);
extern void tgai_test(void);

#endif //__IMAGE_H__
