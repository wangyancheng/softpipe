#include <string.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "tgaimage.h"
#include "tile_cache.h"


TGAColor tgacolor_rgba(unsigned char R, unsigned char G, unsigned char B, unsigned char A)
{
	TGAColor color;

	color.r = R;
	color.g = G;
	color.b = B;
	color.a = A;
	color.bytespp = 4;

	return color;
}

TGAColor tgacolor_val(int v, int bpp)
{
	TGAColor color;

	color.val = v;
	color.bytespp = bpp;

	return color;
}

TGAColor tgacolor_raw(const unsigned char *p, int bpp)
{
	TGAColor color;

	color.val = 0;
	for (int i=0; i<bpp; i++) {
		color.raw[i] = p[i];
	}
	color.bytespp = bpp;

	return color;
}

TGAImage tgai_create(int w, int h, int bpp)
{
	TGAImage img;
	unsigned long nbytes = w*h*bpp;

	img.data = malloc(nbytes);
	memset(img.data, 0, nbytes);
	img.width = w;
	img.height = h;
	img.bytespp = bpp;

	return img;
}

void tgai_release(TGAImage *img)
{
	free(img->data);
	img->data = NULL;
}

void tgai_copy(TGAImage *destin, TGAImage *source)
{
	if(destin->data != NULL)
		free(destin->data);
	
	destin->width = source->width;
	destin->height = source->height;
	destin->bytespp = source->bytespp;
	unsigned long nbytes = destin->width*destin->height*destin->bytespp;
	destin->data = malloc(nbytes);
	memcpy(destin->data, source->data, nbytes);
}

int tgai_get_bytespp(TGAImage *img) {
	return img->bytespp;
}

int tgai_get_width(TGAImage *img) {
	return img->width;
}

int tgai_get_height(TGAImage *img) {
	return img->height;
}

unsigned char *tgai_buffer(TGAImage *img) {
	return img->data;
}

void tgai_clear(TGAImage *img) {
	memset((void *)img->data, 0, img->width*img->height*img->bytespp);
}

#if 0
bool tgai_scale(TGAImage *img, int w, int h) {
	if (w<=0 || h<=0 || !img->data) return false;
	unsigned char *tdata = malloc(w*h*img->bytespp);
	int nscanline = 0;
	int oscanline = 0;
	int erry = 0;
	unsigned long nlinebytes = w*img->bytespp;
	unsigned long olinebytes = img->width*img->bytespp;
	for (int j=0; j<img->height; j++) {
		int errx = img->width-w;
		int nx   = -img->bytespp;
		int ox   = -img->bytespp;
		for (int i=0; i<img->width; i++) {
			ox += img->bytespp;
			errx += w;
			while (errx>=(int)img->width) {
				errx -= img->width;
				nx += img->bytespp;
				memcpy(tdata+nscanline+nx, img->data+oscanline+ox, img->bytespp);
			}
		}
		erry += h;
		oscanline += olinebytes;
		while (erry>=(int)img->height) {
			if (erry>=(int)img->height<<1) // it means we jump over a scanline
				memcpy(tdata+nscanline+nlinebytes, tdata+nscanline, nlinebytes);
			erry -= img->height;
			nscanline += nlinebytes;
		}
	}
	free(img->data);
	img->data = tdata;
	img->width = w;
	img->height = h;
	return true;
}
#endif

TGAColor tgai_get(TGAImage *img, int x, int y) {
	TGAColor color = tgacolor_val(0, 1);
	if (!img->data || x<0 || y<0 || x>=img->width || y>=img->height) {
		return color;
	}
	return tgacolor_raw(img->data+(x+y*img->width)*img->bytespp, img->bytespp);
}

bool tgai_set(TGAImage *img, int x, int y, TGAColor c) {
	if (!img->data || x<0 || y<0 || x>=img->width || y>=img->height) {
		return false;
	}
	memcpy(img->data+(x+y*img->width)*img->bytespp, c.raw, img->bytespp);
	return true;
}

bool tgai_flip_vertically(TGAImage *img) 
{
	if (!img->data) return false;
	unsigned long bytes_per_line = img->width*img->bytespp;
	unsigned char *line = malloc(bytes_per_line);
	int half = img->height>>1;
	for (int j=0; j<half; j++) {
		unsigned long l1 = j*bytes_per_line;
		unsigned long l2 = (img->height-1-j)*bytes_per_line;
		memmove((void *)line,      (void *)(img->data+l1), bytes_per_line);
		memmove((void *)(img->data+l1), (void *)(img->data+l2), bytes_per_line);
		memmove((void *)(img->data+l2), (void *)line,      bytes_per_line);
	}
	free(line);
	return true;
}

bool tgai_flip_horizontally(TGAImage *img) {
	if (!img->data) return false;
	int half = img->width>>1;
	for (int i=0; i<half; i++) {
		for (int j=0; j<img->height; j++) {
			TGAColor c1 = tgai_get(img, i, j);
			TGAColor c2 = tgai_get(img, img->width-1-i, j);
			tgai_set(img, i, j, c2);
			tgai_set(img, img->width-1-i, j, c1);
		}
	}
	return true;
}

bool load_rle_data(TGAImage *img, FILE *fp) {
	unsigned long pixelcount = img->width*img->height;
	unsigned long currentpixel = 0;
	unsigned long currentbyte  = 0;
	TGAColor colorbuffer;
	int ret;
	do {
		unsigned char chunkheader = 0;
		chunkheader = fgetc(fp);
		if ((char)chunkheader == EOF) {
			printf("an error occured while reading the data\n");
			return false;
		}
		if (chunkheader<128) {
			chunkheader++;
			for (int i=0; i<chunkheader; i++) {
				ret = fread((char *)colorbuffer.raw, img->bytespp, 1, fp);
				if (ret != 1) {
					printf("an error occured while reading the header\n");
					return false;
				}
				for (int t=0; t<img->bytespp; t++)
					img->data[currentbyte++] = colorbuffer.raw[t];
				currentpixel++;
				if (currentpixel>pixelcount) {
					printf("Too many pixels read\n");
					return false;
				}
			}
		} else {
			chunkheader -= 127;
			ret = fread((char *)colorbuffer.raw, img->bytespp, 1, fp);
			if (ret != 1) {
				printf("an error occured while reading the header\n");
				return false;
			}
			for (int i=0; i<chunkheader; i++) {
				for (int t=0; t<img->bytespp; t++)
					img->data[currentbyte++] = colorbuffer.raw[t];
				currentpixel++;
				if (currentpixel>pixelcount) {
					printf("Too many pixels read\n");
					return false;
				}
			}
		}
	} while (currentpixel < pixelcount);
	return true;
}

// TODO: it is not necessary to break a raw chunk for two equal pixels (for the matter of the resulting size)
bool unload_rle_data(TGAImage *img, FILE *fp) {
	const unsigned char max_chunk_length = 128;
	unsigned long npixels = img->width*img->height;
	unsigned long curpix = 0;
	int ret;
	while (curpix<npixels) {
		unsigned long chunkstart = curpix*img->bytespp;
		unsigned long curbyte = curpix*img->bytespp;
		unsigned char run_length = 1;
		bool raw = true;
		while (curpix+run_length<npixels && run_length<max_chunk_length) {
			bool succ_eq = true;
			for (int t=0; succ_eq && t<img->bytespp; t++) {
				succ_eq = (img->data[curbyte+t]==img->data[curbyte+t+img->bytespp]);
			}
			curbyte += img->bytespp;
			if (1==run_length) {
				raw = !succ_eq;
			}
			if (raw && succ_eq) {
				run_length--;
				break;
			}
			if (!raw && !succ_eq) {
				break;
			}
			run_length++;
		}
		curpix += run_length;
		ret = fputc(raw?run_length-1:run_length+127, fp);
		if (ret == EOF) {
			printf("can't dump the tga file\n");
			return false;
		}
		ret = fwrite((char *)(img->data+chunkstart), (raw?run_length*img->bytespp:img->bytespp), 1, fp);
		if (ret != 1) {
			printf("can't dump the tga file\n");
			return false;
		}
	}
	return true;
}

bool tgai_read_file(TGAImage *img, const char *filename)
{
	if(img->data != NULL)
		free(img->data);
	img->data = NULL;

	FILE *fp = fopen(filename, "rb");
	if(fp == NULL){
		printf("can't open file %s\n", filename);
		return false;
	}
	
	TGA_Header header;
	int ret;
	ret = fread(&header, sizeof(TGA_Header), 1, fp);
	if(ret != 1){
		printf("an error occured while reading the header\n");
		fclose(fp);
		return false;
	}
	img->width = header.width;
	img->height = header.height;
	img->bytespp = header.bitsperpixel>>3;
	if (img->width<=0 || img->height<=0 || (img->bytespp!=GRAYSCALE && img->bytespp!=RGB && img->bytespp!=RGBA)){
		fclose(fp);
		printf("bad bpp (or width/height) value\n");
		return false;
	}

	unsigned long nbytes = img->bytespp*img->width*img->height;
	img->data = malloc(nbytes);
	//2:uncompressed true-color image , 3:uncompressed grayscale image
	if (3==header.datatypecode || 2==header.datatypecode) {
		ret = fread(img->data, nbytes, 1, fp);
		if (ret != 1) {
			fclose(fp);
			printf("an error occured while reading the data\n");
			return false;
		}
	//10:run-length encoded true-color image, 11 run-length encoded grayscale image
	}else if (10==header.datatypecode||11==header.datatypecode) {
	if (!load_rle_data(img, fp)) {
		fclose(fp);
		printf("an error occured while reading the data\n");
		return false;
	}
	} else {
		fclose(fp);
		printf("unknown file format %d\n", (int)header.datatypecode);
		return false;
	}
	//Image descriptor (1 byte): bits 3–0 give the alpha channel depth, bits 5–4 give pixel ordering
	if (!(header.imagedescriptor & 0x20)) {
		tgai_flip_vertically(img);
	}
	if (header.imagedescriptor & 0x10) {
		tgai_flip_horizontally(img);
	}

	printf("%d X %d / %d\n", img->width, img->height, img->bytespp*8);
	fclose(fp);
	return true;
}

bool tgai_write_tga_file(TGAImage *img, const char *filename, bool rle) {
	unsigned char developer_area_ref[4] = {0, 0, 0, 0};
	unsigned char extension_area_ref[4] = {0, 0, 0, 0};
	unsigned char footer[18] = {'T','R','U','E','V','I','S','I','O','N','-','X','F','I','L','E','.','\0'};
	FILE *fp = fopen(filename, "w");
	if(fp == NULL){
		printf("can't open file %s\n", filename);
		return false;
	} 
	
	TGA_Header header;
	memset((void *)&header, 0, sizeof(header));
	header.bitsperpixel = img->bytespp<<3;
	header.width  = img->width;
	header.height = img->height;
	header.datatypecode = (img->bytespp==GRAYSCALE?(rle?11:3):(rle?10:2));
	header.imagedescriptor = 0x20; // top-left origin
	int ret;
	ret = fwrite((char *)&header, sizeof(header), 1, fp);
	if (ret != 1) {
		fclose(fp);
		printf("can't dump the tga file\n");
		return false;
	}
	if (!rle) {
		ret = fwrite((char *)img->data, img->width*img->height*img->bytespp, 1, fp);
		if (ret != 1) {
			printf("can't unload raw data\n");
			fclose(fp);
			return false;
		}
	} else {
		if (!unload_rle_data(img, fp)) {
			fclose(fp);
			printf("can't unload rle data\n");
			return false;
		}
	}
	ret = fwrite((char *)developer_area_ref, sizeof(developer_area_ref), 1, fp);
	if (ret != 1) {
		printf("can't dump the tga file\n");
		fclose(fp);
		return false;
	}
	ret = fwrite((char *)extension_area_ref, sizeof(extension_area_ref), 1, fp);
	if (ret != 1) {
		printf("can't dump the tga file\n");
		fclose(fp);
		return false;
	}
	ret = fwrite((char *)footer, sizeof(footer), 1, fp);
	if (ret != 1) {
		printf("can't dump the tga file\n");
		fclose(fp);
		return false;
	}
	fclose(fp);
	return true;
}

void tgai_test(void)
{
	TGAColor color1 = tgacolor_rgba(255, 0,   0,   255);
	printf("color1: r:%02x, g:%02x, b:%02x, a:%02x\n", color1.r, color1.g, color1.b, color1.a);
	TGAColor color2 = tgacolor_val(0x103f55ff, RGBA);
	printf("color2: r:%02x, g:%02x, b:%02x, a:%02x\n", color2.r, color2.g, color2.b, color2.a);
	unsigned char raw[4] = {0x01, 0x20, 0x44, 0x5f};
	TGAColor color3 = tgacolor_raw(raw, RGBA);
	printf("color3: r:%02x, g:%02x, b:%02x, a:%02x\n", color3.r, color3.g, color3.b, color3.a);

	TGAImage image = tgai_create(40, 40, RGB);
	int bytespp = tgai_get_bytespp(&image);
	int width = tgai_get_width(&image);
	int height = tgai_get_height(&image);
	printf("image width: %d, height: %d, bytespp: %d\n", width, height, bytespp);
	unsigned char *buffer = tgai_buffer(&image);
	memset(buffer, 0x55, image.width*image.height*image.bytespp);
	tgai_write_tga_file(&image, "output1.tga", true);
	//tgai_clear(&image);
	//tgai_write_tga_file(&image, "output.tga", true);
	//memset(buffer, 0x55, image.width*image.height*image.bytespp);
	//tgai_scale(&image, 80, 80);
	//tgai_write_tga_file(&image, "output2.tga", true);

    tgai_clear(&image);
	for(int i = 0; i < 10; i++){
		for(int j = 0; j < 20; j++){
			if(j > i){
				tgai_set(&image, i, j, color1);
			}
		}
	}
	tgai_write_tga_file(&image, "output3.tga", true);

	TGAColor color4 = tgacolor_rgba(0, 255,   0,   255);
    TGAImage image2 = tgai_create(40, 40, RGB);
	tgai_read_file(&image2, "output3.tga");
	for(int i = 20; i < 30; i++){
		for(int j = 20; j < 30; j++){
			tgai_set(&image, i, j, color4);
		}
	}
	tgai_write_tga_file(&image, "output4.tga", true);
}