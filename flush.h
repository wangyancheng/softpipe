#ifndef __FLUSH_H__
#define __FLUSH_H__

#include <stdio.h>

/**
 * Flags for the flush function.
 */
enum pipe_flush_flags
{
   PIPE_FLUSH_END_OF_FRAME = (1 << 0),
   PIPE_FLUSH_DEFERRED = (1 << 1),
   PIPE_FLUSH_FENCE_FD = (1 << 2),
   PIPE_FLUSH_ASYNC = (1 << 3),
   PIPE_FLUSH_HINT_FINISH = (1 << 4),
   PIPE_FLUSH_TOP_OF_PIPE = (1 << 5),
   PIPE_FLUSH_BOTTOM_OF_PIPE = (1 << 6),
};


#define SP_FLUSH_TEXTURE_CACHE  0x2


#define util_snprintf snprintf

void softpipe_flush( struct softpipe_context *softpipe,
                    unsigned flags);
                    
#endif