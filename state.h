#ifndef __STATE_H__
#define __STATE_H__

#include <stdint.h>

typedef unsigned char      ubyte;

#define PIPE_MAX_CLIP_PLANES       8
#define PIPE_MAX_COLOR_BUFS        1   //应该为8,暂时改为1
#define PIPE_MAX_SHADER_INPUTS    80 /* 32 GENERIC + 32 PATCH + 16 others */
#define PIPE_MAX_SHADER_OUTPUTS   80 /* 32 GENERIC + 32 PATCH + 16 others */
#define PIPE_MAX_VIEWPORTS        16

/**
 * Inequality functions.  Used for depth test, stencil compare, alpha
 * test, shadow compare, etc.
 */
enum pipe_compare_func {
   PIPE_FUNC_NEVER,
   PIPE_FUNC_LESS,
   PIPE_FUNC_EQUAL,
   PIPE_FUNC_LEQUAL,
   PIPE_FUNC_GREATER,
   PIPE_FUNC_NOTEQUAL,
   PIPE_FUNC_GEQUAL,
   PIPE_FUNC_ALWAYS,
};

/**
 * 光栅化状态信息
 */
struct pipe_rasterizer_state
{
   unsigned flatshade:1;				//渲染模式，ShadeModel SMOOTH, FLAT
   unsigned light_twoside:1;			//双面光照，glLightModeli(LIGHT_TWO_SIDE, GL_TRUE)
   unsigned clamp_vertex_color:1;		//顶点颜色截断，ClampColor(CLAMP_VERTEX_COLOR)
   unsigned clamp_fragment_color:1;		//片段颜色截断，ClampColor(CLAMP_FRAGMENT_COLOR)
   unsigned front_ccw:1;				//多边形的正向，FrontFace CCW, CW
   unsigned cull_face:2;      			//多边形剔除，glCullFace(GL_FRONT)， PIPE_FACE_x
   unsigned fill_front:2;     			//多边形正面渲染模式, PolygonMode(face, mode), PIPE_POLYGON_MODE_x 
   unsigned fill_back:2;      			//多边形背面渲染模式，PolygonMode(face, mode), PIPE_POLYGON_MODE_x
   unsigned offset_point:1;				//多边形点深度偏移使能位 glEnable(POLYGON_OFFSET_POINT)
   unsigned offset_line:1;				//多边形线深度偏移使能位 glEnable(POLYGON_OFFSET_LINE)
   unsigned offset_tri:1;				//多边形填充深度偏移使能位 glEnable(POLYGON_OFFSET_FILL)
   unsigned scissor:1;					//裁剪测试使能位, glEnable(SCISSOR_TEST)
   unsigned poly_smooth:1;				//多边形抗锯齿使能位，glEnable(GL_POLYGON_SMOOTH)
   unsigned poly_stipple_enable:1;		//多边形点画使能位，glEnable(GL_POLYGON_STIPPLE)
   unsigned point_smooth:1;			 	//点抗锯齿使能位，glEnable(GL_POINT_SMOOTH)
   unsigned sprite_coord_mode:1;    	//点精灵纹理坐标原点 glPointParameter(GL_POINT_SPRITE_COORD_ORIGIN,GL_LOWER_LEFT/GL_UPPER_LEFT) PIPE_SPRITE_COORD_ 
   unsigned point_quad_rasterization:1; //点精灵使能位 glEnable(GL_POINT_SPRITE) points rasterized as quads or points
   unsigned point_tri_clip:1; /** large points clipped as tris or points */
   unsigned point_size_per_vertex:1; 	//启动顶点程序大小位， glEnable(GL_PROGRAM_POINT_SIZE) size computed in vertex shader
   unsigned multisample:1;         		//MSAA启动位， glEnable(GL_MULTISAMPLE) XXX maybe more ms state in future 
   unsigned force_persample_interp:1;	//超级采样抗锯齿SSAA启动位, glEnable(SAMPLE_SHADING_ARB) 
   unsigned line_smooth:1;				//线抗锯齿使能位，glEnable(GL_LINE_SMOOT)
   unsigned line_stipple_enable:1;		//线点画使能位, glEnable(GL_LINE_STIPPLE)
   unsigned line_last_pixel:1;			//
   unsigned conservative_raster_mode:2; //保守光栅化模式，ConservativeRasterParameteriNV (CONSERVATIVE_RASTER_MODE_NV, param)  CONSERVATIVE_RASTER_MODE_POST_SNAP_NV/CONSERVATIVE_RASTER_MODE_PRE_SNAP_TRIANGLES_NV
   unsigned flatshade_first:1;		//挑衅顶点， ProvokingVertex(LAST_VERTEX_CONVENTION/FIRST_VERTEX_CONVENTION)
   unsigned half_pixel_center:1;	//光栅化像素偏移位（没有通过API设置
   unsigned bottom_edge_rule:1;		//屏幕是否下边开始 glClipControl(origin, depth) GL_LOWER_LEFT/GL_UPPER_LEF
   unsigned subpixel_precision_x:4;	//保守光栅化子像素精度偏差 SubpixelPrecisionBiasNV(uint xbits, uint ybits)
   unsigned subpixel_precision_y:4;	//保守光栅化子像素精度偏差 SubpixelPrecisionBiasNV(uint xbits, uint ybits
   unsigned rasterizer_discard:1;   //光栅化禁止位，为 true 时，光栅化将被禁用，不会写入任何像素。这只对 Stream Out 功能有意义。
   /**
    * 由 PIPE_CAP_TILE_RASTER_ORDER 公开。
    * 当为 true 时，tile_raster_order_increasing_* 指示光栅化器应渲染图块的顺序，
    * 以满足 GL_MESA_tile_raster_order 的要求。
    */
   unsigned tile_raster_order_fixed:1;
   unsigned tile_raster_order_increasing_x:1;
   unsigned tile_raster_order_increasing_y:1;

   /**
    * 如果为 false，则深度裁剪将被禁用，并且深度值将在深度测试之前在每个像素级别被限制。
    * 这取决于 PIPE_CAP_DEPTH_CLIP_DISABLE。
    *
    * 如果不支持 PIPE_CAP_DEPTH_CLIP_DISABLE_SEPARATE，则depth_clip_near 等于depth_clip_far。
    */
   unsigned depth_clip_near:1;			//glEnable(GL_DEPTH_CLAMP) 深度截取
   unsigned depth_clip_far:1;			//glEnable(GL_DEPTH_CLAMP)
   /**
    * 当 true 时，z 轴上的裁剪空间从 [0..1] 开始 (D3D)。当 false 时，从 [-1, 1] 开始 (GL)。
    * 
    * 注意：D3D 将始终使用深度限制。
    */
   unsigned clip_halfz:1;			//glClipControl:depth深度控制的深度模式

   /**
    * 为 true 时不缩放 offset_units，并对 unorm 和浮点深度缓冲区 (D3D9) 使用相同的规则。
    * 为 false 时使用 GL/D3D1X 行为。这取决于 PIPE_CAP_POLYGON_OFFSET_UNITS_UNSCALED。
    */
   unsigned offset_units_unscaled:1;

   /**
    * 启用用于剪切半空间的位。这适用于用户剪切平面和着色器剪切距离。
    * 请注意，如果绑定的着色器导出任何剪切距离，这些距离将替换所有用户剪切平面，
    * 并且此处启用但着色器未写入的剪切半空间计数为禁用。
    */
   unsigned clip_plane_enable:PIPE_MAX_CLIP_PLANES;		//glEnable(GL_CLIP_DISTANCEi)用户自定义裁剪平面使能

   unsigned line_stipple_factor:8;  /**< [1..256] actually */
   unsigned line_stipple_pattern:16;

   /**
    * 用点坐标替换给定的 TEXCOORD 输入，最多 8 个输入。
    * 如果不支持 TEXCOORD（包括 PCOORD），则用 GENERIC 输入替换。
    * 最多 9 个输入：8 个 GENERIC 模拟 TEXCOORD，1 个 GENERIC 模拟 PCOORD。
    */
   uint16_t sprite_coord_enable; /* 0-7: TEXCOORD/GENERIC, 8: PCOORD */

   float line_width;
   float point_size;           /**< used when no per-vertex size */
   float offset_units;			//glPolygonOffset()
   float offset_scale;			//glPolygonOffset()
   float offset_clamp;			//EXT_polygon_offset_clamp()
   float conservative_raster_dilate;	//保守光栅化附加的配盖章值（像素为单位）ConservativeRasterParameterfNV(CONSERVATIVE_RASTER_DILATE_NV, float value)
};

/****************************************************************
 * 纹理、表面和顶点数据的格式
 *****************************************************************/
enum pipe_format {
   PIPE_FORMAT_NONE                    = 0,
   PIPE_FORMAT_B8G8R8A8_UNORM          = 1,
   PIPE_FORMAT_B8G8R8X8_UNORM          = 2,
   PIPE_FORMAT_A8R8G8B8_UNORM          = 3,
   PIPE_FORMAT_X8R8G8B8_UNORM          = 4,
   PIPE_FORMAT_B5G5R5A1_UNORM          = 5,
   PIPE_FORMAT_B4G4R4A4_UNORM          = 6,
   PIPE_FORMAT_B5G6R5_UNORM            = 7,
   PIPE_FORMAT_R10G10B10A2_UNORM       = 8,
   PIPE_FORMAT_L8_UNORM                = 9,    /**< ubyte luminance */
   PIPE_FORMAT_A8_UNORM                = 10,   /**< ubyte alpha */
   PIPE_FORMAT_I8_UNORM                = 11,   /**< ubyte intensity */
   PIPE_FORMAT_L8A8_UNORM              = 12,   /**< ubyte alpha, luminance */
   PIPE_FORMAT_L16_UNORM               = 13,   /**< ushort luminance */
   PIPE_FORMAT_UYVY                    = 14,
   PIPE_FORMAT_YUYV                    = 15,
   PIPE_FORMAT_Z16_UNORM               = 16,
   PIPE_FORMAT_Z32_UNORM               = 17,
   PIPE_FORMAT_Z32_FLOAT               = 18,
   PIPE_FORMAT_Z24_UNORM_S8_UINT       = 19,
   PIPE_FORMAT_S8_UINT_Z24_UNORM       = 20,
   PIPE_FORMAT_Z24X8_UNORM             = 21,
   PIPE_FORMAT_X8Z24_UNORM             = 22,
   PIPE_FORMAT_S8_UINT                 = 23,   /**< ubyte stencil */
   PIPE_FORMAT_R64_FLOAT               = 24,
   PIPE_FORMAT_R64G64_FLOAT            = 25,
   PIPE_FORMAT_R64G64B64_FLOAT         = 26,
   PIPE_FORMAT_R64G64B64A64_FLOAT      = 27,
   PIPE_FORMAT_R32_FLOAT               = 28,
   PIPE_FORMAT_R32G32_FLOAT            = 29,
   PIPE_FORMAT_R32G32B32_FLOAT         = 30,
   PIPE_FORMAT_R32G32B32A32_FLOAT      = 31,
   PIPE_FORMAT_R32_UNORM               = 32,
   PIPE_FORMAT_R32G32_UNORM            = 33,
   PIPE_FORMAT_R32G32B32_UNORM         = 34,
   PIPE_FORMAT_R32G32B32A32_UNORM      = 35,
   PIPE_FORMAT_R32_USCALED             = 36,
   PIPE_FORMAT_R32G32_USCALED          = 37,
   PIPE_FORMAT_R32G32B32_USCALED       = 38,
   PIPE_FORMAT_R32G32B32A32_USCALED    = 39,
   PIPE_FORMAT_R32_SNORM               = 40,
   PIPE_FORMAT_R32G32_SNORM            = 41,
   PIPE_FORMAT_R32G32B32_SNORM         = 42,
   PIPE_FORMAT_R32G32B32A32_SNORM      = 43,
   PIPE_FORMAT_R32_SSCALED             = 44,
   PIPE_FORMAT_R32G32_SSCALED          = 45,
   PIPE_FORMAT_R32G32B32_SSCALED       = 46,
   PIPE_FORMAT_R32G32B32A32_SSCALED    = 47,
   PIPE_FORMAT_R16_UNORM               = 48,
   PIPE_FORMAT_R16G16_UNORM            = 49,
   PIPE_FORMAT_R16G16B16_UNORM         = 50,
   PIPE_FORMAT_R16G16B16A16_UNORM      = 51,
   PIPE_FORMAT_R16_USCALED             = 52,
   PIPE_FORMAT_R16G16_USCALED          = 53,
   PIPE_FORMAT_R16G16B16_USCALED       = 54,
   PIPE_FORMAT_R16G16B16A16_USCALED    = 55,
   PIPE_FORMAT_R16_SNORM               = 56,
   PIPE_FORMAT_R16G16_SNORM            = 57,
   PIPE_FORMAT_R16G16B16_SNORM         = 58,
   PIPE_FORMAT_R16G16B16A16_SNORM      = 59,
   PIPE_FORMAT_R16_SSCALED             = 60,
   PIPE_FORMAT_R16G16_SSCALED          = 61,
   PIPE_FORMAT_R16G16B16_SSCALED       = 62,
   PIPE_FORMAT_R16G16B16A16_SSCALED    = 63,
   PIPE_FORMAT_R8_UNORM                = 64,
   PIPE_FORMAT_R8G8_UNORM              = 65,
   PIPE_FORMAT_R8G8B8_UNORM            = 66,
   PIPE_FORMAT_R8G8B8A8_UNORM          = 67,
   PIPE_FORMAT_X8B8G8R8_UNORM          = 68,
   PIPE_FORMAT_R8_USCALED              = 69,
   PIPE_FORMAT_R8G8_USCALED            = 70,
   PIPE_FORMAT_R8G8B8_USCALED          = 71,
   PIPE_FORMAT_R8G8B8A8_USCALED        = 72,
   PIPE_FORMAT_R8_SNORM                = 74,
   PIPE_FORMAT_R8G8_SNORM              = 75,
   PIPE_FORMAT_R8G8B8_SNORM            = 76,
   PIPE_FORMAT_R8G8B8A8_SNORM          = 77,
   PIPE_FORMAT_R8_SSCALED              = 82,
   PIPE_FORMAT_R8G8_SSCALED            = 83,
   PIPE_FORMAT_R8G8B8_SSCALED          = 84,
   PIPE_FORMAT_R8G8B8A8_SSCALED        = 85,
   PIPE_FORMAT_R32_FIXED               = 87,
   PIPE_FORMAT_R32G32_FIXED            = 88,
   PIPE_FORMAT_R32G32B32_FIXED         = 89,
   PIPE_FORMAT_R32G32B32A32_FIXED      = 90,
   PIPE_FORMAT_R16_FLOAT               = 91,
   PIPE_FORMAT_R16G16_FLOAT            = 92,
   PIPE_FORMAT_R16G16B16_FLOAT         = 93,
   PIPE_FORMAT_R16G16B16A16_FLOAT      = 94,

   /* sRGB formats */
   PIPE_FORMAT_L8_SRGB                 = 95,
   PIPE_FORMAT_L8A8_SRGB               = 96,
   PIPE_FORMAT_R8G8B8_SRGB             = 97,
   PIPE_FORMAT_A8B8G8R8_SRGB           = 98,
   PIPE_FORMAT_X8B8G8R8_SRGB           = 99,
   PIPE_FORMAT_B8G8R8A8_SRGB           = 100,
   PIPE_FORMAT_B8G8R8X8_SRGB           = 101,
   PIPE_FORMAT_A8R8G8B8_SRGB           = 102,
   PIPE_FORMAT_X8R8G8B8_SRGB           = 103,
   PIPE_FORMAT_R8G8B8A8_SRGB           = 104,

   /* compressed formats */
   PIPE_FORMAT_DXT1_RGB                = 105,
   PIPE_FORMAT_DXT1_RGBA               = 106,
   PIPE_FORMAT_DXT3_RGBA               = 107,
   PIPE_FORMAT_DXT5_RGBA               = 108,

   /* sRGB, compressed */
   PIPE_FORMAT_DXT1_SRGB               = 109,
   PIPE_FORMAT_DXT1_SRGBA              = 110,
   PIPE_FORMAT_DXT3_SRGBA              = 111,
   PIPE_FORMAT_DXT5_SRGBA              = 112,

   /* rgtc compressed */
   PIPE_FORMAT_RGTC1_UNORM             = 113,
   PIPE_FORMAT_RGTC1_SNORM             = 114,
   PIPE_FORMAT_RGTC2_UNORM             = 115,
   PIPE_FORMAT_RGTC2_SNORM             = 116,

   PIPE_FORMAT_R8G8_B8G8_UNORM         = 117,
   PIPE_FORMAT_G8R8_G8B8_UNORM         = 118,

   /* mixed formats */
   PIPE_FORMAT_R8SG8SB8UX8U_NORM       = 119,
   PIPE_FORMAT_R5SG5SB6U_NORM          = 120,

   /* TODO: re-order these */
   PIPE_FORMAT_A8B8G8R8_UNORM          = 121,
   PIPE_FORMAT_B5G5R5X1_UNORM          = 122,
   PIPE_FORMAT_R10G10B10A2_USCALED     = 123,
   PIPE_FORMAT_R11G11B10_FLOAT         = 124,
   PIPE_FORMAT_R9G9B9E5_FLOAT          = 125,
   PIPE_FORMAT_Z32_FLOAT_S8X24_UINT    = 126,
   PIPE_FORMAT_R1_UNORM                = 127,
   PIPE_FORMAT_R10G10B10X2_USCALED     = 128,
   PIPE_FORMAT_R10G10B10X2_SNORM       = 129,
   PIPE_FORMAT_L4A4_UNORM              = 130,
   PIPE_FORMAT_B10G10R10A2_UNORM       = 131,
   PIPE_FORMAT_R10SG10SB10SA2U_NORM    = 132,
   PIPE_FORMAT_R8G8Bx_SNORM            = 133,
   PIPE_FORMAT_R8G8B8X8_UNORM          = 134,
   PIPE_FORMAT_B4G4R4X4_UNORM          = 135,

   /* some stencil samplers formats */
   PIPE_FORMAT_X24S8_UINT              = 136,
   PIPE_FORMAT_S8X24_UINT              = 137,
   PIPE_FORMAT_X32_S8X24_UINT          = 138,

   PIPE_FORMAT_B2G3R3_UNORM            = 139,
   PIPE_FORMAT_L16A16_UNORM            = 140,
   PIPE_FORMAT_A16_UNORM               = 141,
   PIPE_FORMAT_I16_UNORM               = 142,

   PIPE_FORMAT_LATC1_UNORM             = 143,
   PIPE_FORMAT_LATC1_SNORM             = 144,
   PIPE_FORMAT_LATC2_UNORM             = 145,
   PIPE_FORMAT_LATC2_SNORM             = 146,

   PIPE_FORMAT_A8_SNORM                = 147,
   PIPE_FORMAT_L8_SNORM                = 148,
   PIPE_FORMAT_L8A8_SNORM              = 149,
   PIPE_FORMAT_I8_SNORM                = 150,
   PIPE_FORMAT_A16_SNORM               = 151,
   PIPE_FORMAT_L16_SNORM               = 152,
   PIPE_FORMAT_L16A16_SNORM            = 153,
   PIPE_FORMAT_I16_SNORM               = 154,

   PIPE_FORMAT_A16_FLOAT               = 155,
   PIPE_FORMAT_L16_FLOAT               = 156,
   PIPE_FORMAT_L16A16_FLOAT            = 157,
   PIPE_FORMAT_I16_FLOAT               = 158,
   PIPE_FORMAT_A32_FLOAT               = 159,
   PIPE_FORMAT_L32_FLOAT               = 160,
   PIPE_FORMAT_L32A32_FLOAT            = 161,
   PIPE_FORMAT_I32_FLOAT               = 162,

   PIPE_FORMAT_YV12                    = 163,
   PIPE_FORMAT_YV16                    = 164,
   PIPE_FORMAT_IYUV                    = 165,  /**< aka I420 */
   PIPE_FORMAT_NV12                    = 166,
   PIPE_FORMAT_NV21                    = 167,

   PIPE_FORMAT_A4R4_UNORM              = 168,
   PIPE_FORMAT_R4A4_UNORM              = 169,
   PIPE_FORMAT_R8A8_UNORM              = 170,
   PIPE_FORMAT_A8R8_UNORM              = 171,

   PIPE_FORMAT_R10G10B10A2_SSCALED     = 172,
   PIPE_FORMAT_R10G10B10A2_SNORM       = 173,

   PIPE_FORMAT_B10G10R10A2_USCALED     = 174,
   PIPE_FORMAT_B10G10R10A2_SSCALED     = 175,
   PIPE_FORMAT_B10G10R10A2_SNORM       = 176,

   PIPE_FORMAT_R8_UINT                 = 177,
   PIPE_FORMAT_R8G8_UINT               = 178,
   PIPE_FORMAT_R8G8B8_UINT             = 179,
   PIPE_FORMAT_R8G8B8A8_UINT           = 180,

   PIPE_FORMAT_R8_SINT                 = 181,
   PIPE_FORMAT_R8G8_SINT               = 182,
   PIPE_FORMAT_R8G8B8_SINT             = 183,
   PIPE_FORMAT_R8G8B8A8_SINT           = 184,

   PIPE_FORMAT_R16_UINT                = 185,
   PIPE_FORMAT_R16G16_UINT             = 186,
   PIPE_FORMAT_R16G16B16_UINT          = 187,
   PIPE_FORMAT_R16G16B16A16_UINT       = 188,

   PIPE_FORMAT_R16_SINT                = 189,
   PIPE_FORMAT_R16G16_SINT             = 190,
   PIPE_FORMAT_R16G16B16_SINT          = 191,
   PIPE_FORMAT_R16G16B16A16_SINT       = 192,

   PIPE_FORMAT_R32_UINT                = 193,
   PIPE_FORMAT_R32G32_UINT             = 194,
   PIPE_FORMAT_R32G32B32_UINT          = 195,
   PIPE_FORMAT_R32G32B32A32_UINT       = 196,

   PIPE_FORMAT_R32_SINT                = 197,
   PIPE_FORMAT_R32G32_SINT             = 198,
   PIPE_FORMAT_R32G32B32_SINT          = 199,
   PIPE_FORMAT_R32G32B32A32_SINT       = 200,

   PIPE_FORMAT_A8_UINT                 = 201,
   PIPE_FORMAT_I8_UINT                 = 202,
   PIPE_FORMAT_L8_UINT                 = 203,
   PIPE_FORMAT_L8A8_UINT               = 204,

   PIPE_FORMAT_A8_SINT                 = 205,
   PIPE_FORMAT_I8_SINT                 = 206,
   PIPE_FORMAT_L8_SINT                 = 207,
   PIPE_FORMAT_L8A8_SINT               = 208,

   PIPE_FORMAT_A16_UINT                = 209,
   PIPE_FORMAT_I16_UINT                = 210,
   PIPE_FORMAT_L16_UINT                = 211,
   PIPE_FORMAT_L16A16_UINT             = 212,

   PIPE_FORMAT_A16_SINT                = 213,
   PIPE_FORMAT_I16_SINT                = 214,
   PIPE_FORMAT_L16_SINT                = 215,
   PIPE_FORMAT_L16A16_SINT             = 216,

   PIPE_FORMAT_A32_UINT                = 217,
   PIPE_FORMAT_I32_UINT                = 218,
   PIPE_FORMAT_L32_UINT                = 219,
   PIPE_FORMAT_L32A32_UINT             = 220,

   PIPE_FORMAT_A32_SINT                = 221,
   PIPE_FORMAT_I32_SINT                = 222,
   PIPE_FORMAT_L32_SINT                = 223,
   PIPE_FORMAT_L32A32_SINT             = 224,

   PIPE_FORMAT_B10G10R10A2_UINT        = 225, 

   PIPE_FORMAT_ETC1_RGB8               = 226,

   PIPE_FORMAT_R8G8_R8B8_UNORM         = 227,
   PIPE_FORMAT_G8R8_B8R8_UNORM         = 228,

   PIPE_FORMAT_R8G8B8X8_SNORM          = 229,
   PIPE_FORMAT_R8G8B8X8_SRGB           = 230,
   PIPE_FORMAT_R8G8B8X8_UINT           = 231,
   PIPE_FORMAT_R8G8B8X8_SINT           = 232,
   PIPE_FORMAT_B10G10R10X2_UNORM       = 233,
   PIPE_FORMAT_R16G16B16X16_UNORM      = 234,
   PIPE_FORMAT_R16G16B16X16_SNORM      = 235,
   PIPE_FORMAT_R16G16B16X16_FLOAT      = 236,
   PIPE_FORMAT_R16G16B16X16_UINT       = 237,
   PIPE_FORMAT_R16G16B16X16_SINT       = 238,
   PIPE_FORMAT_R32G32B32X32_FLOAT      = 239,
   PIPE_FORMAT_R32G32B32X32_UINT       = 240,
   PIPE_FORMAT_R32G32B32X32_SINT       = 241,

   PIPE_FORMAT_R8A8_SNORM              = 242,
   PIPE_FORMAT_R16A16_UNORM            = 243,
   PIPE_FORMAT_R16A16_SNORM            = 244,
   PIPE_FORMAT_R16A16_FLOAT            = 245,
   PIPE_FORMAT_R32A32_FLOAT            = 246,
   PIPE_FORMAT_R8A8_UINT               = 247,
   PIPE_FORMAT_R8A8_SINT               = 248,
   PIPE_FORMAT_R16A16_UINT             = 249,
   PIPE_FORMAT_R16A16_SINT             = 250,
   PIPE_FORMAT_R32A32_UINT             = 251,
   PIPE_FORMAT_R32A32_SINT             = 252,
   PIPE_FORMAT_R10G10B10A2_UINT        = 253,

   PIPE_FORMAT_B5G6R5_SRGB             = 254,

   PIPE_FORMAT_BPTC_RGBA_UNORM         = 255,
   PIPE_FORMAT_BPTC_SRGBA              = 256,
   PIPE_FORMAT_BPTC_RGB_FLOAT          = 257,
   PIPE_FORMAT_BPTC_RGB_UFLOAT         = 258,

   PIPE_FORMAT_A8L8_UNORM              = 259,
   PIPE_FORMAT_A8L8_SNORM              = 260,
   PIPE_FORMAT_A8L8_SRGB               = 261,
   PIPE_FORMAT_A16L16_UNORM            = 262,

   PIPE_FORMAT_G8R8_UNORM              = 263,
   PIPE_FORMAT_G8R8_SNORM              = 264,
   PIPE_FORMAT_G16R16_UNORM            = 265,
   PIPE_FORMAT_G16R16_SNORM            = 266,

   PIPE_FORMAT_A8B8G8R8_SNORM          = 267,
   PIPE_FORMAT_X8B8G8R8_SNORM          = 268,

   PIPE_FORMAT_ETC2_RGB8               = 269,
   PIPE_FORMAT_ETC2_SRGB8              = 270,
   PIPE_FORMAT_ETC2_RGB8A1             = 271,
   PIPE_FORMAT_ETC2_SRGB8A1            = 272,
   PIPE_FORMAT_ETC2_RGBA8              = 273,
   PIPE_FORMAT_ETC2_SRGBA8             = 274,
   PIPE_FORMAT_ETC2_R11_UNORM          = 275,
   PIPE_FORMAT_ETC2_R11_SNORM          = 276,
   PIPE_FORMAT_ETC2_RG11_UNORM         = 277,
   PIPE_FORMAT_ETC2_RG11_SNORM         = 278,

   PIPE_FORMAT_ASTC_4x4                = 279,
   PIPE_FORMAT_ASTC_5x4                = 280,
   PIPE_FORMAT_ASTC_5x5                = 281,
   PIPE_FORMAT_ASTC_6x5                = 282,
   PIPE_FORMAT_ASTC_6x6                = 283,
   PIPE_FORMAT_ASTC_8x5                = 284,
   PIPE_FORMAT_ASTC_8x6                = 285,
   PIPE_FORMAT_ASTC_8x8                = 286,
   PIPE_FORMAT_ASTC_10x5               = 287,
   PIPE_FORMAT_ASTC_10x6               = 288,
   PIPE_FORMAT_ASTC_10x8               = 289,
   PIPE_FORMAT_ASTC_10x10              = 290,
   PIPE_FORMAT_ASTC_12x10              = 291,
   PIPE_FORMAT_ASTC_12x12              = 292,

   PIPE_FORMAT_ASTC_4x4_SRGB           = 293,
   PIPE_FORMAT_ASTC_5x4_SRGB           = 294,
   PIPE_FORMAT_ASTC_5x5_SRGB           = 295,
   PIPE_FORMAT_ASTC_6x5_SRGB           = 296,
   PIPE_FORMAT_ASTC_6x6_SRGB           = 297,
   PIPE_FORMAT_ASTC_8x5_SRGB           = 298,
   PIPE_FORMAT_ASTC_8x6_SRGB           = 299,
   PIPE_FORMAT_ASTC_8x8_SRGB           = 300,
   PIPE_FORMAT_ASTC_10x5_SRGB          = 301,
   PIPE_FORMAT_ASTC_10x6_SRGB          = 302,
   PIPE_FORMAT_ASTC_10x8_SRGB          = 303,
   PIPE_FORMAT_ASTC_10x10_SRGB         = 304,
   PIPE_FORMAT_ASTC_12x10_SRGB         = 305,
   PIPE_FORMAT_ASTC_12x12_SRGB         = 306,

   PIPE_FORMAT_P016                    = 307,

   PIPE_FORMAT_R10G10B10X2_UNORM       = 308,
   PIPE_FORMAT_A1B5G5R5_UNORM          = 309,
   PIPE_FORMAT_X1B5G5R5_UNORM          = 310,
   PIPE_FORMAT_A4B4G4R4_UNORM          = 311,

   PIPE_FORMAT_COUNT
};

//表面描述
union pipe_surface_desc {
   struct {
      unsigned level;
      unsigned first_layer:16;
      unsigned last_layer:16;
   } tex;
   struct {
      unsigned first_element;
      unsigned last_element;
   } buf;
};
/**
 * 纹理类型
*/
enum pipe_texture_target
{
   PIPE_BUFFER,
   PIPE_TEXTURE_1D,
   PIPE_TEXTURE_2D,
   PIPE_TEXTURE_3D,
   PIPE_TEXTURE_CUBE,
   PIPE_TEXTURE_RECT,
   PIPE_TEXTURE_1D_ARRAY,
   PIPE_TEXTURE_2D_ARRAY,
   PIPE_TEXTURE_CUBE_ARRAY,
   PIPE_MAX_TEXTURE_TYPES,
};

#define PIPE_TEXTURE_BUFFER_SIZE 1920*1080
#define PIPE_TEXTURE_BUFFER_NUMBER 1
 /**************************************************************************
 * 内存对象/资源，例如顶点缓冲区或纹理
 ***************************************************************************/
struct pipe_resource
{
   unsigned width0; 		//buffer的大小（48）/**< Used by both buffers and textures. */
   uint16_t height0; 		//buffer固定为1 /* Textures: The maximum height/depth/array_size is 16k. */
   uint16_t depth0;			//深度（1）
   uint16_t array_size;		//矩阵大小（1）层数量

   enum pipe_format format:16;      //纹理、表面和顶点数据的格式（64）     /**< PIPE_FORMAT_x */
   enum pipe_texture_target target:8; 		//纹理类型（0） /**< PIPE_TEXTURE_x */
   unsigned last_level:8;    //当前/定义的最后mipmap的索引（多级渐远纹理）（0）	/**< Index of last mipmap level present/defined */

   /** Number of samples determining quality, driving rasterizer, shading,
    *  and framebuffer.
    */
   unsigned nr_samples:8;	//决定质量、驱动光栅器、着色和帧缓冲区的样本数（0）

    /***
    *一个像素内的多个样本可以具有相同的值。 nr_storage_samples 确定每个像素有多少个不同值的插槽。 
    *只有颜色缓冲区可以将其设置为低于 nr_samples
    **/
   unsigned nr_storage_samples:8;		//（0）

   unsigned usage:8;         /**< PIPE_USAGE_x (not a bitmask) */
   unsigned bind;            /**< bitmask of PIPE_BIND_x */
   unsigned flags;           /**< bitmask of PIPE_RESOURCE_FLAG_x */

   void *userdata; 
};

 /*
 * 可以绑定到颜色渲染目标/深度模板附着点的纹理视图
 */
struct pipe_surface
{
   enum pipe_format format:16;				//数据格式
   unsigned writable:1;          //可写着色器资源		/**< writable shader resource */
   struct pipe_resource *texture; 	//这是一个视图的资源				/**< resource into which this is a view  */

   /* XXX width/height should be removed */
   //应删除XXX宽度/高度
   uint16_t width;         //以像素为单位的逻辑宽度		      /**< logical width in pixels */
   uint16_t height;        //以像素为单位的逻辑高度		      /**< logical height in pixels */

   union pipe_surface_desc u;		//表面描述
};

 /*
 * 传输对象使用标志
 */
enum pipe_transfer_usage
{
   //资源内容在传输创建时回读（或直接访问）
   PIPE_TRANSFER_READ = (1 << 0),
   //资源内容在transfer_unmap时写回（或直接访问而修改）
   PIPE_TRANSFER_WRITE = (1 << 1),
   //读修改写入
   PIPE_TRANSFER_READ_WRITE = PIPE_TRANSFER_READ | PIPE_TRANSFER_WRITE,

   /** 
    * 传输应直接映射纹理存储。驾驶员可以
    * 如果不可能，返回null，并且状态跟踪器需要应对
    * 这样并使用没有此标志的替代路径。
    *
    * 例如。状态跟踪器可以具有更简单的路径，该路径可以映射纹理和
    * 会直接在它们上读取/修改/写周期，并且更复杂
    * 使用最小读写转移的路径。
    */
   PIPE_TRANSFER_MAP_DIRECTLY = (1 << 2),

   /**
    * 丢弃映射区域内的内存。
    *
    * 它不应与pipe_transfer_read一起使用。
    *
    * 参见：
    * -  OpenGL的ARB_MAP_BUFFER_RANGE扩展名，map_invalidate_range_bit标志。
    */
   PIPE_TRANSFER_DISCARD_RANGE = (1 << 8),

   /**
    * 如果无法立即映射资源，则失败。
    *
    * See also:
    * - Direct3D's D3DLOCK_DONOTWAIT flag.
    * - Mesa's MESA_MAP_NOWAIT_BIT flag.
    * - WDDM's D3DDDICB_LOCKFLAGS.DonotWait flag.
    */
   PIPE_TRANSFER_DONTBLOCK = (1 << 9),

   /**
    * 在映射时，请勿尝试同步在资源上的待处理操作。
    *
    * 它不应与pipe_transfer_read一起使用。
    *
    * See also:
    * - OpenGL's ARB_map_buffer_range extension, MAP_UNSYNCHRONIZED_BIT flag.
    * - Direct3D's D3DLOCK_NOOVERWRITE flag.
    * - WDDM's D3DDDICB_LOCKFLAGS.IgnoreSync flag.
    */
   PIPE_TRANSFER_UNSYNCHRONIZED = (1 << 10),

   /**
    * 书面范围将在以后通知
    * pipe_context :: Transfer_flush_region。
    *
    * 它不应与pipe_transfer_read一起使用。
    *
    * See also:
    * - pipe_context::transfer_flush_region
    * - OpenGL's ARB_map_buffer_range extension, MAP_FLUSH_EXPLICIT_BIT flag.
    */
   PIPE_TRANSFER_FLUSH_EXPLICIT = (1 << 11),

   /**
    * 丢弃所有备份资源的内存。
    *
    * 它不应与pipe_transfer_read一起使用。
    *
    * This is equivalent to:
    * - OpenGL's ARB_map_buffer_range extension, MAP_INVALIDATE_BUFFER_BIT
    * - BufferData(NULL) on a GL buffer
    * - Direct3D's D3DLOCK_DISCARD flag.
    * - WDDM's D3DDDICB_LOCKFLAGS.Discard flag.
    * - D3D10 DDI's D3D10_DDI_MAP_WRITE_DISCARD flag
    * - D3D10's D3D10_MAP_WRITE_DISCARD flag.
    */
   PIPE_TRANSFER_DISCARD_WHOLE_RESOURCE = (1 << 12),

   /**
    * 允许在映射时将资源用于渲染。
    *
    * PIPE_RESOURCE_FLAG_MAP_PERSISTENT must be set when creating
    * the resource.
    *
    * If COHERENT is not set, memory_barrier(PIPE_BARRIER_MAPPED_BUFFER)
    * must be called to ensure the device can see what the CPU has written.
    */
   PIPE_TRANSFER_PERSISTENT = (1 << 13),

   /**
    * 如果设置了持久性，则可以确保设备完成的任何写作
    * 立即可见CPU，反之亦然。
    *
    * PIPE_RESOURCE_FLAG_MAP_COHERENT must be set when creating
    * the resource.
    */
   PIPE_TRANSFER_COHERENT = (1 << 14)
};

/**
 * 1D/2D/3D 图像资源的子区域。
*/
struct pipe_box
{
   /* 仅由纹理使用的字段使用int16_t而不是int。
    * X和宽度由缓冲区使用，因此他们需要完整的32位范围。
    */
   int x;
   int16_t y;
   int16_t z;
   int width;
   int16_t height;
   int16_t depth;
};

#define TRANSFER_MAX 16
/**
 * 转移对象，用于向/从资源传输数据
*/
struct pipe_transfer
{
   struct pipe_resource *resource; 		//管道资源
   unsigned level;                     //纹理mipmap level
   enum pipe_transfer_usage usage;	   //传输使用标志
   struct pipe_box box;                //要访问的资源区域
   unsigned stride;                    //步幅
   unsigned layer_stride;              //图像/层跨度
};

/**
 * 描述如何从规定的格式中打包/将像素打包到/将像素打包。
 *
 * xxx：可以将其重命名为util_format_pack之类的东西，也可以分解
 * 在util_format_block中的标志中，说出了我们想要的。
 */
enum util_format_layout {
   /**
    * Formats with util_format_block::width == util_format_block::height == 1
    * that can be described as an ordinary data structure.
    */
    /*
    * util_format_block::width=util_format_block::height==1的格式可以描述为普通数据格式。
    */
   UTIL_FORMAT_LAYOUT_PLAIN = 0,

   /**
    * Formats with sub-sampled channels.
    *
    * This is for formats like YVYU where there is less than one sample per
    * pixel.
    */
    /*
    * 具有子采样通道的格式。
    * 这适用于像YVYU这样每个像素样本少于一个的格式
    */
   UTIL_FORMAT_LAYOUT_SUBSAMPLED = 3,

   /**
    * S3 Texture Compression formats.
    */
    /*
    * S3 纹理压缩格式
    */
   UTIL_FORMAT_LAYOUT_S3TC = 4,

   /**
    * Red-Green Texture Compression formats.
    */
    /*
    * 红-绿纹理压缩格式
    */
   UTIL_FORMAT_LAYOUT_RGTC = 5,

   /**
    * Ericsson Texture Compression
    */
    /*
    * 爱立信纹理压缩
    */
   UTIL_FORMAT_LAYOUT_ETC = 6,

   /**
    * BC6/7 Texture Compression
    */
    /*
    * BC6/7 纹理压缩
    */
   UTIL_FORMAT_LAYOUT_BPTC = 7,

   /**
    * ASTC
    */
   UTIL_FORMAT_LAYOUT_ASTC = 8,

   /**
    * Everything else that doesn't fit in any of the above layouts.
    */
    /*
    * 不适合上述任何布局的所有其他内容
    */
   UTIL_FORMAT_LAYOUT_OTHER = 9
};

/**
 * 像素块尺寸
 */
struct util_format_block
{
   /** Block width in pixels */
   unsigned width;
   
   /** Block height in pixels */
   unsigned height;

   /** Block size in bits */
   unsigned bits;
};

enum util_format_colorspace {
   UTIL_FORMAT_COLORSPACE_RGB = 0,
   UTIL_FORMAT_COLORSPACE_SRGB = 1,
   UTIL_FORMAT_COLORSPACE_YUV = 2,
   UTIL_FORMAT_COLORSPACE_ZS = 3
};

/**
 * 通道描述
 */
struct util_format_channel_description
{
   unsigned type:5;        //通道类型/**< UTIL_FORMAT_TYPE_x */
   unsigned normalized:1;	//归一化
   unsigned pure_integer:1;	//纯整数
   unsigned size:9;        	//每个通道的bit数/**< bits per channel */
   unsigned shift:16;      //偏移/** number of bits from lsb */
};

/*
* 通用描述符形式
*/
struct util_format_description
{
   enum pipe_format format;		//数据格式

   const char *name;		         //名字

   /**
    * 简短名称，前缀条纹，较低案例。
    */
   const char *short_name;

   /**
    * Pixel block dimensions.
    */
    //像素块尺寸
   struct util_format_block block;

   enum util_format_layout layout;		//像素的打包/解包格式

   /**
    * The number of channels.
    */
   unsigned nr_channels:3;		//通道数量(应该是图像通道)

   /**
    * 所有通道是否具有相同数量的（整个）字节和类型。
    */
   unsigned is_array:1;

   /**
    * 像素格式是否可以描述为比特场结构。
    *
    * In particular:
    * - pixel depth must be 8, 16, or 32 bits;
    * - all channels must be unsigned, signed, or void
    */
   unsigned is_bitmask:1;

   /**
    * Whether channels have mixed types (ignoring UTIL_FORMAT_TYPE_VOID).
    */
   unsigned is_mixed:1;

   /**
    * 输入通道描述，按XYZW顺序。
    *
    * Only valid for UTIL_FORMAT_LAYOUT_PLAIN formats.
    *
    * If each channel is accessed as an individual N-byte value, X is always
    * at the lowest address in memory, Y is always next, and so on.  For all
    * currently-defined formats, the N-byte value has native endianness.
    *
    * If instead a group of channels is accessed as a single N-byte value,
    * the order of the channels within that value depends on endianness.
    * For big-endian targets, X is the most significant subvalue,
    * otherwise it is the least significant one.
    *
    * For example, if X is 8 bits and Y is 24 bits, the memory order is:
    *
    *                 0  1  2  3
    *  little-endian: X  Yl Ym Yu    (l = lower, m = middle, u = upper)
    *  big-endian:    X  Yu Ym Yl
    *
    * If X is 5 bits, Y is 5 bits, Z is 5 bits and W is 1 bit, the layout is:
    *
    *                        0        1
    *                 msb  lsb msb  lsb
    *  little-endian: YYYXXXXX WZZZZZYY
    *  big-endian:    XXXXXYYY YYZZZZZW
    */
   struct util_format_channel_description channel[4];	//通道描述
   /**
    * 输出通道swizzle。
    *
    * The order is either:
    * - RGBA
    * - YUV(A)
    * - ZS
    * depending on the colorspace.
    */
   unsigned char swizzle[4];

   /**
    * Colorspace转换。
    */
   enum util_format_colorspace colorspace;

};

/***********************************************************************
 * framebuffer状态
 * 请注意，pipe_surface是"用于渲染的纹理视图",因此在ARB_framebuffer_no_attachment
 * 的情况下，没有可用的pipe_surface状态，因此我们可以提取样本和层的数量。
 ************************************************************************/
struct pipe_framebuffer_state
{
   uint16_t width, height;		//宽，高
   uint16_t layers;  	// 无附件帧缓冲区中的层数	
   ubyte samples; 	// 无附件帧缓冲区中的样本数	

   ubyte nr_cbufs;  //多个渲染目标的多个颜色缓冲区
   struct pipe_surface *cbufs[PIPE_MAX_COLOR_BUFS];

   struct pipe_surface *zsbuf;     //Z/模板缓冲区 /**< Z/stencil buffer */
};

struct pipe_scissor_state
{
   unsigned minx:16;
   unsigned miny:16;
   unsigned maxx:16;
   unsigned maxy:16;
};

//深度状态
struct pipe_depth_state
{
   unsigned enabled:1;         //深度测试使能位（glEnable(GL_DEPTH_TEST) /**< depth test enabled? */
   unsigned writemask:1;       //允许深度buffer写位(glDepthMask(GL_FALSE)) /**< allow depth buffer writes? */
   unsigned func:3;            //深度测试功能 (glDepthFunc(GL_LESS)) /**< depth test func (PIPE_FUNC_x) */
   unsigned bounds_test:1;     //深度界限测试使能(不知道) /**< depth bounds test enabled? */
   float bounds_min;           //深度的最小值（glDepthRange(GLclampd near, GLclampd far)） /**< minimum depth bound */
   float bounds_max;           //深度的最大值（glDepthRange(GLclampd near, GLclampd far)） /**< maximum depth bound */
};

//模板状态
struct pipe_stencil_state
{
   unsigned enabled:1;  	//模板测试使能（glEnable(GL_STENCIL_TEST)/**< stencil[0]: stencil enabled, stencil[1]: two-side enabled */
   unsigned func:3;     	//模板测试功能 (void glStencilFunc(GLenum func, GLint ref, GLuint mask)) /**< PIPE_FUNC_x */
   unsigned fail_op:3;  	//如果模板测试失败将采取的动作(void glStencilOp(GLenum sfail, GLenum dpfail, GLenum dppass)) /**< PIPE_STENCIL_OP_x */
   unsigned zpass_op:3; 	//如果深度测试和模板测试都通过，将采取的动作 /**< PIPE_STENCIL_OP_x */
   unsigned zfail_op:3; 	//如果模板测试通过，但深度测试失败时将采取的动作/**< PIPE_STENCIL_OP_x */
   unsigned valuemask:8;	//模板值(void glStencilFunc(GLenum func, GLint ref, GLuint mask))
   unsigned writemask:8;	//模板可写入状态(void glStencilFunc(GLenum func, GLint ref, GLuint mask))
};

//阿尔法状态
struct pipe_alpha_state
{
   unsigned enabled:1;	//阿尔法使能glEnable（GL_ALPHA_TEST）
   unsigned func:3;     //阿尔法功能glAlphaFunc（GLenum_func，GLclampf ref） /**< PIPE_FUNC_x */
   float ref_value;     //阿尔法参考值glAlphaFunc（GLenum_func，GLclampf ref） /**< reference value */
};

//深度、模板、阿尔法状态
struct pipe_depth_stencil_alpha_state
{
   struct pipe_depth_state depth;		//深度状态
   struct pipe_stencil_state stencil[2]; //模板状态	/**< [0] = front, [1] = back */
   struct pipe_alpha_state alpha;		//阿尔法状态
};

enum pipe_blendfactor {
   PIPE_BLENDFACTOR_ONE = 1,
   PIPE_BLENDFACTOR_SRC_COLOR,
   PIPE_BLENDFACTOR_SRC_ALPHA,
   PIPE_BLENDFACTOR_DST_ALPHA,
   PIPE_BLENDFACTOR_DST_COLOR,
   PIPE_BLENDFACTOR_SRC_ALPHA_SATURATE,
   PIPE_BLENDFACTOR_CONST_COLOR,
   PIPE_BLENDFACTOR_CONST_ALPHA,
   PIPE_BLENDFACTOR_SRC1_COLOR,
   PIPE_BLENDFACTOR_SRC1_ALPHA,

   PIPE_BLENDFACTOR_ZERO = 0x11,
   PIPE_BLENDFACTOR_INV_SRC_COLOR,
   PIPE_BLENDFACTOR_INV_SRC_ALPHA,
   PIPE_BLENDFACTOR_INV_DST_ALPHA,
   PIPE_BLENDFACTOR_INV_DST_COLOR,

   PIPE_BLENDFACTOR_INV_CONST_COLOR = 0x17,
   PIPE_BLENDFACTOR_INV_CONST_ALPHA,
   PIPE_BLENDFACTOR_INV_SRC1_COLOR,
   PIPE_BLENDFACTOR_INV_SRC1_ALPHA,
};

enum pipe_blend_func {
   PIPE_BLEND_ADD,
   PIPE_BLEND_SUBTRACT,
   PIPE_BLEND_REVERSE_SUBTRACT,
   PIPE_BLEND_MIN,
   PIPE_BLEND_MAX,
};

enum pipe_logicop {
   PIPE_LOGICOP_CLEAR,
   PIPE_LOGICOP_NOR,
   PIPE_LOGICOP_AND_INVERTED,
   PIPE_LOGICOP_COPY_INVERTED,
   PIPE_LOGICOP_AND_REVERSE,
   PIPE_LOGICOP_INVERT,
   PIPE_LOGICOP_XOR,
   PIPE_LOGICOP_NAND,
   PIPE_LOGICOP_AND,
   PIPE_LOGICOP_EQUIV,
   PIPE_LOGICOP_NOOP,
   PIPE_LOGICOP_OR_INVERTED,
   PIPE_LOGICOP_COPY,
   PIPE_LOGICOP_OR_REVERSE,
   PIPE_LOGICOP_OR,
   PIPE_LOGICOP_SET,
};

struct pipe_rt_blend_state
{
   unsigned blend_enable:1;      //混合使能

   unsigned rgb_func:3;          // < PIPE_BLEND_x 
   unsigned rgb_src_factor:5;    //< PIPE_BLENDFACTOR_x 
   unsigned rgb_dst_factor:5;    //< PIPE_BLENDFACTOR_x 

   unsigned alpha_func:3;        //< PIPE_BLEND_x 
   unsigned alpha_src_factor:5;  //< PIPE_BLENDFACTOR_x 
   unsigned alpha_dst_factor:5;  //< PIPE_BLENDFACTOR_x 

   unsigned colormask:4;         //< bitmask of PIPE_MASK_R/G/B/A 
};

//混合状态
struct pipe_blend_state
{
   unsigned independent_blend_enable:1;      //独立混合使能
   unsigned logicop_enable:1;                //逻辑运算使能
   unsigned logicop_func:4;                  //逻辑运算 < PIPE_LOGICOP_x >
   unsigned dither:1;                        //抖动使能
   unsigned alpha_to_coverage:1;             //SAMPLE_ALPHA_TO_COVERAGE使能
   unsigned alpha_to_one:1;                  //GL_SAMPLE_ALPHA_TO_ONE使能
   struct pipe_rt_blend_state rt[PIPE_MAX_COLOR_BUFS];   //混合状态
};

#endif