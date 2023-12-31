/* This file is autogenerated by tracetool, do not edit. */

#ifndef TRACE_HW_WATCHDOG_GENERATED_TRACERS_H
#define TRACE_HW_WATCHDOG_GENERATED_TRACERS_H

#include "trace/control.h"

extern TraceEvent _TRACE_CMSDK_APB_WATCHDOG_READ_EVENT;
extern TraceEvent _TRACE_CMSDK_APB_WATCHDOG_WRITE_EVENT;
extern TraceEvent _TRACE_CMSDK_APB_WATCHDOG_RESET_EVENT;
extern TraceEvent _TRACE_CMSDK_APB_WATCHDOG_LOCK_EVENT;
extern TraceEvent _TRACE_ASPEED_WDT_READ_EVENT;
extern TraceEvent _TRACE_ASPEED_WDT_WRITE_EVENT;
extern TraceEvent _TRACE_SPAPR_WATCHDOG_START_EVENT;
extern TraceEvent _TRACE_SPAPR_WATCHDOG_STOP_EVENT;
extern TraceEvent _TRACE_SPAPR_WATCHDOG_QUERY_EVENT;
extern TraceEvent _TRACE_SPAPR_WATCHDOG_QUERY_LPM_EVENT;
extern TraceEvent _TRACE_SPAPR_WATCHDOG_EXPIRED_EVENT;
extern TraceEvent _TRACE_WATCHDOG_PERFORM_ACTION_EVENT;
extern TraceEvent _TRACE_WATCHDOG_SET_ACTION_EVENT;
extern uint16_t _TRACE_CMSDK_APB_WATCHDOG_READ_DSTATE;
extern uint16_t _TRACE_CMSDK_APB_WATCHDOG_WRITE_DSTATE;
extern uint16_t _TRACE_CMSDK_APB_WATCHDOG_RESET_DSTATE;
extern uint16_t _TRACE_CMSDK_APB_WATCHDOG_LOCK_DSTATE;
extern uint16_t _TRACE_ASPEED_WDT_READ_DSTATE;
extern uint16_t _TRACE_ASPEED_WDT_WRITE_DSTATE;
extern uint16_t _TRACE_SPAPR_WATCHDOG_START_DSTATE;
extern uint16_t _TRACE_SPAPR_WATCHDOG_STOP_DSTATE;
extern uint16_t _TRACE_SPAPR_WATCHDOG_QUERY_DSTATE;
extern uint16_t _TRACE_SPAPR_WATCHDOG_QUERY_LPM_DSTATE;
extern uint16_t _TRACE_SPAPR_WATCHDOG_EXPIRED_DSTATE;
extern uint16_t _TRACE_WATCHDOG_PERFORM_ACTION_DSTATE;
extern uint16_t _TRACE_WATCHDOG_SET_ACTION_DSTATE;
#define TRACE_CMSDK_APB_WATCHDOG_READ_ENABLED 1
#define TRACE_CMSDK_APB_WATCHDOG_WRITE_ENABLED 1
#define TRACE_CMSDK_APB_WATCHDOG_RESET_ENABLED 1
#define TRACE_CMSDK_APB_WATCHDOG_LOCK_ENABLED 1
#define TRACE_ASPEED_WDT_READ_ENABLED 1
#define TRACE_ASPEED_WDT_WRITE_ENABLED 1
#define TRACE_SPAPR_WATCHDOG_START_ENABLED 1
#define TRACE_SPAPR_WATCHDOG_STOP_ENABLED 1
#define TRACE_SPAPR_WATCHDOG_QUERY_ENABLED 1
#define TRACE_SPAPR_WATCHDOG_QUERY_LPM_ENABLED 1
#define TRACE_SPAPR_WATCHDOG_EXPIRED_ENABLED 1
#define TRACE_WATCHDOG_PERFORM_ACTION_ENABLED 1
#define TRACE_WATCHDOG_SET_ACTION_ENABLED 1
#include "qemu/log-for-trace.h"
#include "qemu/error-report.h"


#define TRACE_CMSDK_APB_WATCHDOG_READ_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_CMSDK_APB_WATCHDOG_READ) || \
    false)

static inline void _nocheck__trace_cmsdk_apb_watchdog_read(uint64_t offset, uint64_t data, unsigned size)
{
    if (trace_event_get_state(TRACE_CMSDK_APB_WATCHDOG_READ) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 4 "/FEMU/hw/watchdog/trace-events"
            qemu_log("%d@%zu.%06zu:cmsdk_apb_watchdog_read " "CMSDK APB watchdog read: offset 0x%" PRIx64 " data 0x%" PRIx64 " size %u" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , offset, data, size);
#line 67 "trace/trace-hw_watchdog.h"
        } else {
#line 4 "/FEMU/hw/watchdog/trace-events"
            qemu_log("cmsdk_apb_watchdog_read " "CMSDK APB watchdog read: offset 0x%" PRIx64 " data 0x%" PRIx64 " size %u" "\n", offset, data, size);
#line 71 "trace/trace-hw_watchdog.h"
        }
    }
}

static inline void trace_cmsdk_apb_watchdog_read(uint64_t offset, uint64_t data, unsigned size)
{
    if (true) {
        _nocheck__trace_cmsdk_apb_watchdog_read(offset, data, size);
    }
}

#define TRACE_CMSDK_APB_WATCHDOG_WRITE_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_CMSDK_APB_WATCHDOG_WRITE) || \
    false)

static inline void _nocheck__trace_cmsdk_apb_watchdog_write(uint64_t offset, uint64_t data, unsigned size)
{
    if (trace_event_get_state(TRACE_CMSDK_APB_WATCHDOG_WRITE) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 5 "/FEMU/hw/watchdog/trace-events"
            qemu_log("%d@%zu.%06zu:cmsdk_apb_watchdog_write " "CMSDK APB watchdog write: offset 0x%" PRIx64 " data 0x%" PRIx64 " size %u" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , offset, data, size);
#line 98 "trace/trace-hw_watchdog.h"
        } else {
#line 5 "/FEMU/hw/watchdog/trace-events"
            qemu_log("cmsdk_apb_watchdog_write " "CMSDK APB watchdog write: offset 0x%" PRIx64 " data 0x%" PRIx64 " size %u" "\n", offset, data, size);
#line 102 "trace/trace-hw_watchdog.h"
        }
    }
}

static inline void trace_cmsdk_apb_watchdog_write(uint64_t offset, uint64_t data, unsigned size)
{
    if (true) {
        _nocheck__trace_cmsdk_apb_watchdog_write(offset, data, size);
    }
}

#define TRACE_CMSDK_APB_WATCHDOG_RESET_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_CMSDK_APB_WATCHDOG_RESET) || \
    false)

static inline void _nocheck__trace_cmsdk_apb_watchdog_reset(void)
{
    if (trace_event_get_state(TRACE_CMSDK_APB_WATCHDOG_RESET) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 6 "/FEMU/hw/watchdog/trace-events"
            qemu_log("%d@%zu.%06zu:cmsdk_apb_watchdog_reset " "CMSDK APB watchdog: reset" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     );
#line 129 "trace/trace-hw_watchdog.h"
        } else {
#line 6 "/FEMU/hw/watchdog/trace-events"
            qemu_log("cmsdk_apb_watchdog_reset " "CMSDK APB watchdog: reset" "\n");
#line 133 "trace/trace-hw_watchdog.h"
        }
    }
}

static inline void trace_cmsdk_apb_watchdog_reset(void)
{
    if (true) {
        _nocheck__trace_cmsdk_apb_watchdog_reset();
    }
}

#define TRACE_CMSDK_APB_WATCHDOG_LOCK_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_CMSDK_APB_WATCHDOG_LOCK) || \
    false)

static inline void _nocheck__trace_cmsdk_apb_watchdog_lock(uint32_t lock)
{
    if (trace_event_get_state(TRACE_CMSDK_APB_WATCHDOG_LOCK) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 7 "/FEMU/hw/watchdog/trace-events"
            qemu_log("%d@%zu.%06zu:cmsdk_apb_watchdog_lock " "CMSDK APB watchdog: lock %" PRIu32 "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , lock);
#line 160 "trace/trace-hw_watchdog.h"
        } else {
#line 7 "/FEMU/hw/watchdog/trace-events"
            qemu_log("cmsdk_apb_watchdog_lock " "CMSDK APB watchdog: lock %" PRIu32 "\n", lock);
#line 164 "trace/trace-hw_watchdog.h"
        }
    }
}

static inline void trace_cmsdk_apb_watchdog_lock(uint32_t lock)
{
    if (true) {
        _nocheck__trace_cmsdk_apb_watchdog_lock(lock);
    }
}

#define TRACE_ASPEED_WDT_READ_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_ASPEED_WDT_READ) || \
    false)

static inline void _nocheck__trace_aspeed_wdt_read(uint64_t addr, uint32_t size)
{
    if (trace_event_get_state(TRACE_ASPEED_WDT_READ) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 10 "/FEMU/hw/watchdog/trace-events"
            qemu_log("%d@%zu.%06zu:aspeed_wdt_read " "@0x%" PRIx64 " size=%d" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , addr, size);
#line 191 "trace/trace-hw_watchdog.h"
        } else {
#line 10 "/FEMU/hw/watchdog/trace-events"
            qemu_log("aspeed_wdt_read " "@0x%" PRIx64 " size=%d" "\n", addr, size);
#line 195 "trace/trace-hw_watchdog.h"
        }
    }
}

static inline void trace_aspeed_wdt_read(uint64_t addr, uint32_t size)
{
    if (true) {
        _nocheck__trace_aspeed_wdt_read(addr, size);
    }
}

#define TRACE_ASPEED_WDT_WRITE_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_ASPEED_WDT_WRITE) || \
    false)

static inline void _nocheck__trace_aspeed_wdt_write(uint64_t addr, uint32_t size, uint64_t data)
{
    if (trace_event_get_state(TRACE_ASPEED_WDT_WRITE) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 11 "/FEMU/hw/watchdog/trace-events"
            qemu_log("%d@%zu.%06zu:aspeed_wdt_write " "@0x%" PRIx64 " size=%d value=0x%"PRIx64 "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , addr, size, data);
#line 222 "trace/trace-hw_watchdog.h"
        } else {
#line 11 "/FEMU/hw/watchdog/trace-events"
            qemu_log("aspeed_wdt_write " "@0x%" PRIx64 " size=%d value=0x%"PRIx64 "\n", addr, size, data);
#line 226 "trace/trace-hw_watchdog.h"
        }
    }
}

static inline void trace_aspeed_wdt_write(uint64_t addr, uint32_t size, uint64_t data)
{
    if (true) {
        _nocheck__trace_aspeed_wdt_write(addr, size, data);
    }
}

#define TRACE_SPAPR_WATCHDOG_START_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_SPAPR_WATCHDOG_START) || \
    false)

static inline void _nocheck__trace_spapr_watchdog_start(uint64_t flags, uint64_t num, uint64_t timeout)
{
    if (trace_event_get_state(TRACE_SPAPR_WATCHDOG_START) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 14 "/FEMU/hw/watchdog/trace-events"
            qemu_log("%d@%zu.%06zu:spapr_watchdog_start " "Flags 0x%" PRIx64 " num=%" PRId64 " %" PRIu64 "ms" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , flags, num, timeout);
#line 253 "trace/trace-hw_watchdog.h"
        } else {
#line 14 "/FEMU/hw/watchdog/trace-events"
            qemu_log("spapr_watchdog_start " "Flags 0x%" PRIx64 " num=%" PRId64 " %" PRIu64 "ms" "\n", flags, num, timeout);
#line 257 "trace/trace-hw_watchdog.h"
        }
    }
}

static inline void trace_spapr_watchdog_start(uint64_t flags, uint64_t num, uint64_t timeout)
{
    if (true) {
        _nocheck__trace_spapr_watchdog_start(flags, num, timeout);
    }
}

#define TRACE_SPAPR_WATCHDOG_STOP_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_SPAPR_WATCHDOG_STOP) || \
    false)

static inline void _nocheck__trace_spapr_watchdog_stop(uint64_t num, uint64_t ret)
{
    if (trace_event_get_state(TRACE_SPAPR_WATCHDOG_STOP) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 15 "/FEMU/hw/watchdog/trace-events"
            qemu_log("%d@%zu.%06zu:spapr_watchdog_stop " "num=%" PRIu64 " ret=%" PRId64 "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , num, ret);
#line 284 "trace/trace-hw_watchdog.h"
        } else {
#line 15 "/FEMU/hw/watchdog/trace-events"
            qemu_log("spapr_watchdog_stop " "num=%" PRIu64 " ret=%" PRId64 "\n", num, ret);
#line 288 "trace/trace-hw_watchdog.h"
        }
    }
}

static inline void trace_spapr_watchdog_stop(uint64_t num, uint64_t ret)
{
    if (true) {
        _nocheck__trace_spapr_watchdog_stop(num, ret);
    }
}

#define TRACE_SPAPR_WATCHDOG_QUERY_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_SPAPR_WATCHDOG_QUERY) || \
    false)

static inline void _nocheck__trace_spapr_watchdog_query(uint64_t caps)
{
    if (trace_event_get_state(TRACE_SPAPR_WATCHDOG_QUERY) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 16 "/FEMU/hw/watchdog/trace-events"
            qemu_log("%d@%zu.%06zu:spapr_watchdog_query " "caps=0x%" PRIx64 "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , caps);
#line 315 "trace/trace-hw_watchdog.h"
        } else {
#line 16 "/FEMU/hw/watchdog/trace-events"
            qemu_log("spapr_watchdog_query " "caps=0x%" PRIx64 "\n", caps);
#line 319 "trace/trace-hw_watchdog.h"
        }
    }
}

static inline void trace_spapr_watchdog_query(uint64_t caps)
{
    if (true) {
        _nocheck__trace_spapr_watchdog_query(caps);
    }
}

#define TRACE_SPAPR_WATCHDOG_QUERY_LPM_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_SPAPR_WATCHDOG_QUERY_LPM) || \
    false)

static inline void _nocheck__trace_spapr_watchdog_query_lpm(uint64_t caps)
{
    if (trace_event_get_state(TRACE_SPAPR_WATCHDOG_QUERY_LPM) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 17 "/FEMU/hw/watchdog/trace-events"
            qemu_log("%d@%zu.%06zu:spapr_watchdog_query_lpm " "caps=0x%" PRIx64 "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , caps);
#line 346 "trace/trace-hw_watchdog.h"
        } else {
#line 17 "/FEMU/hw/watchdog/trace-events"
            qemu_log("spapr_watchdog_query_lpm " "caps=0x%" PRIx64 "\n", caps);
#line 350 "trace/trace-hw_watchdog.h"
        }
    }
}

static inline void trace_spapr_watchdog_query_lpm(uint64_t caps)
{
    if (true) {
        _nocheck__trace_spapr_watchdog_query_lpm(caps);
    }
}

#define TRACE_SPAPR_WATCHDOG_EXPIRED_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_SPAPR_WATCHDOG_EXPIRED) || \
    false)

static inline void _nocheck__trace_spapr_watchdog_expired(uint64_t num, unsigned action)
{
    if (trace_event_get_state(TRACE_SPAPR_WATCHDOG_EXPIRED) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 18 "/FEMU/hw/watchdog/trace-events"
            qemu_log("%d@%zu.%06zu:spapr_watchdog_expired " "num=%" PRIu64 " action=%u" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , num, action);
#line 377 "trace/trace-hw_watchdog.h"
        } else {
#line 18 "/FEMU/hw/watchdog/trace-events"
            qemu_log("spapr_watchdog_expired " "num=%" PRIu64 " action=%u" "\n", num, action);
#line 381 "trace/trace-hw_watchdog.h"
        }
    }
}

static inline void trace_spapr_watchdog_expired(uint64_t num, unsigned action)
{
    if (true) {
        _nocheck__trace_spapr_watchdog_expired(num, action);
    }
}

#define TRACE_WATCHDOG_PERFORM_ACTION_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_WATCHDOG_PERFORM_ACTION) || \
    false)

static inline void _nocheck__trace_watchdog_perform_action(unsigned int action)
{
    if (trace_event_get_state(TRACE_WATCHDOG_PERFORM_ACTION) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 21 "/FEMU/hw/watchdog/trace-events"
            qemu_log("%d@%zu.%06zu:watchdog_perform_action " "action=%u" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , action);
#line 408 "trace/trace-hw_watchdog.h"
        } else {
#line 21 "/FEMU/hw/watchdog/trace-events"
            qemu_log("watchdog_perform_action " "action=%u" "\n", action);
#line 412 "trace/trace-hw_watchdog.h"
        }
    }
}

static inline void trace_watchdog_perform_action(unsigned int action)
{
    if (true) {
        _nocheck__trace_watchdog_perform_action(action);
    }
}

#define TRACE_WATCHDOG_SET_ACTION_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_WATCHDOG_SET_ACTION) || \
    false)

static inline void _nocheck__trace_watchdog_set_action(unsigned int action)
{
    if (trace_event_get_state(TRACE_WATCHDOG_SET_ACTION) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 22 "/FEMU/hw/watchdog/trace-events"
            qemu_log("%d@%zu.%06zu:watchdog_set_action " "action=%u" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , action);
#line 439 "trace/trace-hw_watchdog.h"
        } else {
#line 22 "/FEMU/hw/watchdog/trace-events"
            qemu_log("watchdog_set_action " "action=%u" "\n", action);
#line 443 "trace/trace-hw_watchdog.h"
        }
    }
}

static inline void trace_watchdog_set_action(unsigned int action)
{
    if (true) {
        _nocheck__trace_watchdog_set_action(action);
    }
}
#endif /* TRACE_HW_WATCHDOG_GENERATED_TRACERS_H */
