/* This file is autogenerated by tracetool, do not edit. */

#ifndef TRACE_HW_BLOCK_DATAPLANE_GENERATED_TRACERS_H
#define TRACE_HW_BLOCK_DATAPLANE_GENERATED_TRACERS_H

#include "trace/control.h"

extern TraceEvent _TRACE_VIRTIO_BLK_DATA_PLANE_START_EVENT;
extern TraceEvent _TRACE_VIRTIO_BLK_DATA_PLANE_STOP_EVENT;
extern uint16_t _TRACE_VIRTIO_BLK_DATA_PLANE_START_DSTATE;
extern uint16_t _TRACE_VIRTIO_BLK_DATA_PLANE_STOP_DSTATE;
#define TRACE_VIRTIO_BLK_DATA_PLANE_START_ENABLED 1
#define TRACE_VIRTIO_BLK_DATA_PLANE_STOP_ENABLED 1
#include "qemu/log-for-trace.h"
#include "qemu/error-report.h"


#define TRACE_VIRTIO_BLK_DATA_PLANE_START_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_VIRTIO_BLK_DATA_PLANE_START) || \
    false)

static inline void _nocheck__trace_virtio_blk_data_plane_start(void * s)
{
    if (trace_event_get_state(TRACE_VIRTIO_BLK_DATA_PLANE_START) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 4 "/FEMU/hw/block/dataplane/trace-events"
            qemu_log("%d@%zu.%06zu:virtio_blk_data_plane_start " "dataplane %p" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , s);
#line 34 "trace/trace-hw_block_dataplane.h"
        } else {
#line 4 "/FEMU/hw/block/dataplane/trace-events"
            qemu_log("virtio_blk_data_plane_start " "dataplane %p" "\n", s);
#line 38 "trace/trace-hw_block_dataplane.h"
        }
    }
}

static inline void trace_virtio_blk_data_plane_start(void * s)
{
    if (true) {
        _nocheck__trace_virtio_blk_data_plane_start(s);
    }
}

#define TRACE_VIRTIO_BLK_DATA_PLANE_STOP_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_VIRTIO_BLK_DATA_PLANE_STOP) || \
    false)

static inline void _nocheck__trace_virtio_blk_data_plane_stop(void * s)
{
    if (trace_event_get_state(TRACE_VIRTIO_BLK_DATA_PLANE_STOP) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 5 "/FEMU/hw/block/dataplane/trace-events"
            qemu_log("%d@%zu.%06zu:virtio_blk_data_plane_stop " "dataplane %p" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , s);
#line 65 "trace/trace-hw_block_dataplane.h"
        } else {
#line 5 "/FEMU/hw/block/dataplane/trace-events"
            qemu_log("virtio_blk_data_plane_stop " "dataplane %p" "\n", s);
#line 69 "trace/trace-hw_block_dataplane.h"
        }
    }
}

static inline void trace_virtio_blk_data_plane_stop(void * s)
{
    if (true) {
        _nocheck__trace_virtio_blk_data_plane_stop(s);
    }
}
#endif /* TRACE_HW_BLOCK_DATAPLANE_GENERATED_TRACERS_H */
