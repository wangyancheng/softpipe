#ifndef __GEOMETRY_H__
#define __GEOMETRY_H__

#include <math.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef struct{
	union{
		struct {float u, v;};
		struct {float x, y;};
		float raw[2];
	};
}Vec2f;
static inline Vec2f Vec2f_init(const float u, const float v)
{
	Vec2f result;
	result.u = u;
	result.v = v;
	return result;
}
static inline Vec2f Vec2f_add(const Vec2f v0, const Vec2f v1){
	Vec2f result;
	result.u = v0.u + v1.u;
	result.v = v0.v + v1.v;
	return result;
}
static inline Vec2f Vec2f_sub(const Vec2f v0, const Vec2f v1){
	Vec2f result;
	result.u = v0.u - v1.u;
	result.v = v0.v - v1.v;
	return result;
}
static inline Vec2f Vec2f_mul(const Vec2f v0, float f){
	Vec2f result;
	result.u = v0.u * f;
	result.v = v0.v * f;
	return result;
}

typedef struct{
	union{
		struct {int u, v;};
		struct {int x, y;};
		int raw[2];
	};
}Vec2i;
static inline Vec2i Vec2i_init(const int u, const int v)
{
	Vec2i result;
	result.u = u;
	result.v = v;
	return result;
}
static inline Vec2i Vec2i_add(const Vec2i v0, const Vec2i v1){
	Vec2i result;
	result.u = v0.u + v1.u;
	result.v = v0.v + v1.v;
	return result;
}
static inline Vec2i Vec2i_sub(const Vec2i v0, const Vec2i v1){
	Vec2i result;
	result.u = v0.u - v1.u;
	result.v = v0.v - v1.v;
	return result;
}
static inline Vec2i Vec2i_mul(const Vec2i v0, int f){
	Vec2i result;
	result.u = v0.u * f;
	result.v = v0.v * f;
	return result;
}

typedef struct{
	union{
		struct {float x, y, z;};
		struct {float ivert, iuv, inorm;};
		float raw[3];
	};
}Vec3f;

static inline Vec3f Vec3f_init(float x, float y, float z)
{
	Vec3f result;
	result.x = x;
	result.y = y;
	result.z = z;
	return result;
}
static inline Vec3f Vec3f_add(const Vec3f v0, const Vec3f v1){
	Vec3f result;
	result.x = v0.x + v1.x;
	result.y = v0.y + v1.y;
	result.z = v0.z + v1.z;
	return result;
}
static inline Vec3f Vec3f_sub(const Vec3f v0, const Vec3f v1){
	Vec3f result;
	result.x = v0.x - v1.x;
	result.y = v0.y - v1.y;
	result.z = v0.z - v1.z;
	return result;
}
static inline Vec3f Vec3f_mul(const Vec3f v0, float f){
	Vec3f result;
	result.x = v0.x * f;
	result.y = v0.y * f;
	result.z = v0.z *f;
	return result;
}

static inline float Vec3f_dot(const Vec3f v0, const Vec3f v1){
	return v0.x*v1.x+v0.y*v1.y+v0.z+v1.z;
}

static inline Vec3f Vec3f_cross(const Vec3f v0, const Vec3f v1){
	Vec3f result;
	result.x = v0.y*v1.z-v0.z*v1.y;
	result.y = v0.z*v1.x-v0.x*v1.z;
	result.z = v0.x*v1.y-v0.y*v1.x;
	return result;
}

static inline float Vec3f_norm(const Vec3f v0){
	return sqrt(v0.x*v0.x+v0.y*v0.y+v0.z*v0.z);
}

static inline Vec3f Vec3f_normalize(const Vec3f v0, float l){
	return Vec3f_mul(v0, l/Vec3f_norm(v0));
}

typedef struct{
	union{
		struct {int x, y, z;};
		struct {int ivert, iuv, inorm;};
		int raw[3];
	};
}Vec3i;

static inline Vec3i Vec3i_init(int x, int y, int z)
{
	Vec3i result;
	result.x = x;
	result.y = y;
	result.z = z;
	return result;
}
static inline Vec3i Vec3i_add(const Vec3i v0, const Vec3i v1){
	Vec3i result;
	result.x = v0.x + v1.x;
	result.y = v0.y + v1.y;
	result.z = v0.z + v1.z;
	return result;
}
static inline Vec3i Vec3i_sub(const Vec3i v0, const Vec3i v1){
	Vec3i result;
	result.x = v0.x - v1.x;
	result.y = v0.y - v1.y;
	result.z = v0.z - v1.z;
	return result;
}
static inline Vec3i Vec3i_mul(const Vec3i v0, int f){
	Vec3i result;
	result.x = v0.x * f;
	result.y = v0.y * f;
	result.z = v0.z *f;
	return result;
}

static inline int Vec3i_dot(const Vec3i v0, const Vec3i v1){
	return v0.x*v1.x+v0.y*v1.y+v0.z+v1.z;
}

static inline Vec3i Vec3i_cross(const Vec3i v0, const Vec3i v1){
	Vec3i result;
	result.x = v0.y*v1.z-v0.z*v1.y;
	result.y = v0.z*v1.x-v0.x*v1.z;
	result.z = v0.x*v1.y-v0.y*v1.x;
	return result;
}

static inline float Vec3i_norm(const Vec3i v0){
	return sqrt(v0.x*v0.x+v0.y*v0.y+v0.z*v0.z);
}

static inline Vec3i Vec3i_normalize(const Vec3i v0, int l){
	return Vec3i_mul(v0, l/Vec3i_norm(v0));
}

#endif //__GEOMETRY_H__