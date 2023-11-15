/* This file is autogenerated by tracetool, do not edit. */

#ifndef TRACE_HW_S390X_GENERATED_TRACERS_H
#define TRACE_HW_S390X_GENERATED_TRACERS_H

#include "trace/control.h"

extern TraceEvent _TRACE_CSS_ENABLE_FACILITY_EVENT;
extern TraceEvent _TRACE_CSS_CRW_EVENT;
extern TraceEvent _TRACE_CSS_CHPID_ADD_EVENT;
extern TraceEvent _TRACE_CSS_NEW_IMAGE_EVENT;
extern TraceEvent _TRACE_CSS_ASSIGN_SUBCH_EVENT;
extern TraceEvent _TRACE_CSS_IO_INTERRUPT_EVENT;
extern TraceEvent _TRACE_CSS_ADAPTER_INTERRUPT_EVENT;
extern TraceEvent _TRACE_CSS_DO_SIC_EVENT;
extern TraceEvent _TRACE_VIRTIO_CCW_INTERPRET_CCW_EVENT;
extern TraceEvent _TRACE_VIRTIO_CCW_NEW_DEVICE_EVENT;
extern TraceEvent _TRACE_VIRTIO_CCW_SET_IND_EVENT;
extern TraceEvent _TRACE_S390_PCI_CLP_CAP_EVENT;
extern TraceEvent _TRACE_S390_PCI_CLP_CAP_SIZE_EVENT;
extern TraceEvent _TRACE_S390_PCI_CLP_DEV_INFO_EVENT;
extern uint16_t _TRACE_CSS_ENABLE_FACILITY_DSTATE;
extern uint16_t _TRACE_CSS_CRW_DSTATE;
extern uint16_t _TRACE_CSS_CHPID_ADD_DSTATE;
extern uint16_t _TRACE_CSS_NEW_IMAGE_DSTATE;
extern uint16_t _TRACE_CSS_ASSIGN_SUBCH_DSTATE;
extern uint16_t _TRACE_CSS_IO_INTERRUPT_DSTATE;
extern uint16_t _TRACE_CSS_ADAPTER_INTERRUPT_DSTATE;
extern uint16_t _TRACE_CSS_DO_SIC_DSTATE;
extern uint16_t _TRACE_VIRTIO_CCW_INTERPRET_CCW_DSTATE;
extern uint16_t _TRACE_VIRTIO_CCW_NEW_DEVICE_DSTATE;
extern uint16_t _TRACE_VIRTIO_CCW_SET_IND_DSTATE;
extern uint16_t _TRACE_S390_PCI_CLP_CAP_DSTATE;
extern uint16_t _TRACE_S390_PCI_CLP_CAP_SIZE_DSTATE;
extern uint16_t _TRACE_S390_PCI_CLP_DEV_INFO_DSTATE;
#define TRACE_CSS_ENABLE_FACILITY_ENABLED 1
#define TRACE_CSS_CRW_ENABLED 1
#define TRACE_CSS_CHPID_ADD_ENABLED 1
#define TRACE_CSS_NEW_IMAGE_ENABLED 1
#define TRACE_CSS_ASSIGN_SUBCH_ENABLED 1
#define TRACE_CSS_IO_INTERRUPT_ENABLED 1
#define TRACE_CSS_ADAPTER_INTERRUPT_ENABLED 1
#define TRACE_CSS_DO_SIC_ENABLED 1
#define TRACE_VIRTIO_CCW_INTERPRET_CCW_ENABLED 1
#define TRACE_VIRTIO_CCW_NEW_DEVICE_ENABLED 1
#define TRACE_VIRTIO_CCW_SET_IND_ENABLED 1
#define TRACE_S390_PCI_CLP_CAP_ENABLED 1
#define TRACE_S390_PCI_CLP_CAP_SIZE_ENABLED 1
#define TRACE_S390_PCI_CLP_DEV_INFO_ENABLED 1
#include "qemu/log-for-trace.h"
#include "qemu/error-report.h"


#define TRACE_CSS_ENABLE_FACILITY_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_CSS_ENABLE_FACILITY) || \
    false)

static inline void _nocheck__trace_css_enable_facility(const char * facility)
{
    if (trace_event_get_state(TRACE_CSS_ENABLE_FACILITY) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 4 "/FEMU/hw/s390x/trace-events"
            qemu_log("%d@%zu.%06zu:css_enable_facility " "CSS: enable %s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , facility);
#line 70 "trace/trace-hw_s390x.h"
        } else {
#line 4 "/FEMU/hw/s390x/trace-events"
            qemu_log("css_enable_facility " "CSS: enable %s" "\n", facility);
#line 74 "trace/trace-hw_s390x.h"
        }
    }
}

static inline void trace_css_enable_facility(const char * facility)
{
    if (true) {
        _nocheck__trace_css_enable_facility(facility);
    }
}

#define TRACE_CSS_CRW_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_CSS_CRW) || \
    false)

static inline void _nocheck__trace_css_crw(uint8_t rsc, uint8_t erc, uint16_t rsid, const char * chained)
{
    if (trace_event_get_state(TRACE_CSS_CRW) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 5 "/FEMU/hw/s390x/trace-events"
            qemu_log("%d@%zu.%06zu:css_crw " "CSS: queueing crw: rsc=0x%x, erc=0x%x, rsid=0x%x %s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , rsc, erc, rsid, chained);
#line 101 "trace/trace-hw_s390x.h"
        } else {
#line 5 "/FEMU/hw/s390x/trace-events"
            qemu_log("css_crw " "CSS: queueing crw: rsc=0x%x, erc=0x%x, rsid=0x%x %s" "\n", rsc, erc, rsid, chained);
#line 105 "trace/trace-hw_s390x.h"
        }
    }
}

static inline void trace_css_crw(uint8_t rsc, uint8_t erc, uint16_t rsid, const char * chained)
{
    if (true) {
        _nocheck__trace_css_crw(rsc, erc, rsid, chained);
    }
}

#define TRACE_CSS_CHPID_ADD_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_CSS_CHPID_ADD) || \
    false)

static inline void _nocheck__trace_css_chpid_add(uint8_t cssid, uint8_t chpid, uint8_t type)
{
    if (trace_event_get_state(TRACE_CSS_CHPID_ADD) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 6 "/FEMU/hw/s390x/trace-events"
            qemu_log("%d@%zu.%06zu:css_chpid_add " "CSS: add chpid %x.%02x (type 0x%02x)" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , cssid, chpid, type);
#line 132 "trace/trace-hw_s390x.h"
        } else {
#line 6 "/FEMU/hw/s390x/trace-events"
            qemu_log("css_chpid_add " "CSS: add chpid %x.%02x (type 0x%02x)" "\n", cssid, chpid, type);
#line 136 "trace/trace-hw_s390x.h"
        }
    }
}

static inline void trace_css_chpid_add(uint8_t cssid, uint8_t chpid, uint8_t type)
{
    if (true) {
        _nocheck__trace_css_chpid_add(cssid, chpid, type);
    }
}

#define TRACE_CSS_NEW_IMAGE_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_CSS_NEW_IMAGE) || \
    false)

static inline void _nocheck__trace_css_new_image(uint8_t cssid, const char * default_cssid)
{
    if (trace_event_get_state(TRACE_CSS_NEW_IMAGE) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 7 "/FEMU/hw/s390x/trace-events"
            qemu_log("%d@%zu.%06zu:css_new_image " "CSS: add css image 0x%02x %s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , cssid, default_cssid);
#line 163 "trace/trace-hw_s390x.h"
        } else {
#line 7 "/FEMU/hw/s390x/trace-events"
            qemu_log("css_new_image " "CSS: add css image 0x%02x %s" "\n", cssid, default_cssid);
#line 167 "trace/trace-hw_s390x.h"
        }
    }
}

static inline void trace_css_new_image(uint8_t cssid, const char * default_cssid)
{
    if (true) {
        _nocheck__trace_css_new_image(cssid, default_cssid);
    }
}

#define TRACE_CSS_ASSIGN_SUBCH_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_CSS_ASSIGN_SUBCH) || \
    false)

static inline void _nocheck__trace_css_assign_subch(const char * do_assign, uint8_t cssid, uint8_t ssid, uint16_t schid, uint16_t devno)
{
    if (trace_event_get_state(TRACE_CSS_ASSIGN_SUBCH) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 8 "/FEMU/hw/s390x/trace-events"
            qemu_log("%d@%zu.%06zu:css_assign_subch " "CSS: %s %x.%x.%04x (devno 0x%04x)" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , do_assign, cssid, ssid, schid, devno);
#line 194 "trace/trace-hw_s390x.h"
        } else {
#line 8 "/FEMU/hw/s390x/trace-events"
            qemu_log("css_assign_subch " "CSS: %s %x.%x.%04x (devno 0x%04x)" "\n", do_assign, cssid, ssid, schid, devno);
#line 198 "trace/trace-hw_s390x.h"
        }
    }
}

static inline void trace_css_assign_subch(const char * do_assign, uint8_t cssid, uint8_t ssid, uint16_t schid, uint16_t devno)
{
    if (true) {
        _nocheck__trace_css_assign_subch(do_assign, cssid, ssid, schid, devno);
    }
}

#define TRACE_CSS_IO_INTERRUPT_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_CSS_IO_INTERRUPT) || \
    false)

static inline void _nocheck__trace_css_io_interrupt(int cssid, int ssid, int schid, uint32_t intparm, uint8_t isc, const char * conditional)
{
    if (trace_event_get_state(TRACE_CSS_IO_INTERRUPT) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 9 "/FEMU/hw/s390x/trace-events"
            qemu_log("%d@%zu.%06zu:css_io_interrupt " "CSS: I/O interrupt on sch %x.%x.%04x (intparm 0x%08x, isc 0x%x) %s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , cssid, ssid, schid, intparm, isc, conditional);
#line 225 "trace/trace-hw_s390x.h"
        } else {
#line 9 "/FEMU/hw/s390x/trace-events"
            qemu_log("css_io_interrupt " "CSS: I/O interrupt on sch %x.%x.%04x (intparm 0x%08x, isc 0x%x) %s" "\n", cssid, ssid, schid, intparm, isc, conditional);
#line 229 "trace/trace-hw_s390x.h"
        }
    }
}

static inline void trace_css_io_interrupt(int cssid, int ssid, int schid, uint32_t intparm, uint8_t isc, const char * conditional)
{
    if (true) {
        _nocheck__trace_css_io_interrupt(cssid, ssid, schid, intparm, isc, conditional);
    }
}

#define TRACE_CSS_ADAPTER_INTERRUPT_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_CSS_ADAPTER_INTERRUPT) || \
    false)

static inline void _nocheck__trace_css_adapter_interrupt(uint8_t isc)
{
    if (trace_event_get_state(TRACE_CSS_ADAPTER_INTERRUPT) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 10 "/FEMU/hw/s390x/trace-events"
            qemu_log("%d@%zu.%06zu:css_adapter_interrupt " "CSS: adapter I/O interrupt (isc 0x%x)" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , isc);
#line 256 "trace/trace-hw_s390x.h"
        } else {
#line 10 "/FEMU/hw/s390x/trace-events"
            qemu_log("css_adapter_interrupt " "CSS: adapter I/O interrupt (isc 0x%x)" "\n", isc);
#line 260 "trace/trace-hw_s390x.h"
        }
    }
}

static inline void trace_css_adapter_interrupt(uint8_t isc)
{
    if (true) {
        _nocheck__trace_css_adapter_interrupt(isc);
    }
}

#define TRACE_CSS_DO_SIC_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_CSS_DO_SIC) || \
    false)

static inline void _nocheck__trace_css_do_sic(uint16_t mode, uint8_t isc)
{
    if (trace_event_get_state(TRACE_CSS_DO_SIC) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 11 "/FEMU/hw/s390x/trace-events"
            qemu_log("%d@%zu.%06zu:css_do_sic " "CSS: set interruption mode 0x%x on isc 0x%x" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , mode, isc);
#line 287 "trace/trace-hw_s390x.h"
        } else {
#line 11 "/FEMU/hw/s390x/trace-events"
            qemu_log("css_do_sic " "CSS: set interruption mode 0x%x on isc 0x%x" "\n", mode, isc);
#line 291 "trace/trace-hw_s390x.h"
        }
    }
}

static inline void trace_css_do_sic(uint16_t mode, uint8_t isc)
{
    if (true) {
        _nocheck__trace_css_do_sic(mode, isc);
    }
}

#define TRACE_VIRTIO_CCW_INTERPRET_CCW_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_VIRTIO_CCW_INTERPRET_CCW) || \
    false)

static inline void _nocheck__trace_virtio_ccw_interpret_ccw(int cssid, int ssid, int schid, int cmd_code)
{
    if (trace_event_get_state(TRACE_VIRTIO_CCW_INTERPRET_CCW) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 14 "/FEMU/hw/s390x/trace-events"
            qemu_log("%d@%zu.%06zu:virtio_ccw_interpret_ccw " "VIRTIO-CCW: %x.%x.%04x: interpret command 0x%x" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , cssid, ssid, schid, cmd_code);
#line 318 "trace/trace-hw_s390x.h"
        } else {
#line 14 "/FEMU/hw/s390x/trace-events"
            qemu_log("virtio_ccw_interpret_ccw " "VIRTIO-CCW: %x.%x.%04x: interpret command 0x%x" "\n", cssid, ssid, schid, cmd_code);
#line 322 "trace/trace-hw_s390x.h"
        }
    }
}

static inline void trace_virtio_ccw_interpret_ccw(int cssid, int ssid, int schid, int cmd_code)
{
    if (true) {
        _nocheck__trace_virtio_ccw_interpret_ccw(cssid, ssid, schid, cmd_code);
    }
}

#define TRACE_VIRTIO_CCW_NEW_DEVICE_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_VIRTIO_CCW_NEW_DEVICE) || \
    false)

static inline void _nocheck__trace_virtio_ccw_new_device(int cssid, int ssid, int schid, int devno, const char * devno_mode)
{
    if (trace_event_get_state(TRACE_VIRTIO_CCW_NEW_DEVICE) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 15 "/FEMU/hw/s390x/trace-events"
            qemu_log("%d@%zu.%06zu:virtio_ccw_new_device " "VIRTIO-CCW: add subchannel %x.%x.%04x, devno 0x%04x (%s)" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , cssid, ssid, schid, devno, devno_mode);
#line 349 "trace/trace-hw_s390x.h"
        } else {
#line 15 "/FEMU/hw/s390x/trace-events"
            qemu_log("virtio_ccw_new_device " "VIRTIO-CCW: add subchannel %x.%x.%04x, devno 0x%04x (%s)" "\n", cssid, ssid, schid, devno, devno_mode);
#line 353 "trace/trace-hw_s390x.h"
        }
    }
}

static inline void trace_virtio_ccw_new_device(int cssid, int ssid, int schid, int devno, const char * devno_mode)
{
    if (true) {
        _nocheck__trace_virtio_ccw_new_device(cssid, ssid, schid, devno, devno_mode);
    }
}

#define TRACE_VIRTIO_CCW_SET_IND_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_VIRTIO_CCW_SET_IND) || \
    false)

static inline void _nocheck__trace_virtio_ccw_set_ind(uint64_t ind_loc, uint8_t ind_old, uint8_t ind_new)
{
    if (trace_event_get_state(TRACE_VIRTIO_CCW_SET_IND) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 16 "/FEMU/hw/s390x/trace-events"
            qemu_log("%d@%zu.%06zu:virtio_ccw_set_ind " "VIRTIO-CCW: indicator at %" PRIu64 ": 0x%x->0x%x" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , ind_loc, ind_old, ind_new);
#line 380 "trace/trace-hw_s390x.h"
        } else {
#line 16 "/FEMU/hw/s390x/trace-events"
            qemu_log("virtio_ccw_set_ind " "VIRTIO-CCW: indicator at %" PRIu64 ": 0x%x->0x%x" "\n", ind_loc, ind_old, ind_new);
#line 384 "trace/trace-hw_s390x.h"
        }
    }
}

static inline void trace_virtio_ccw_set_ind(uint64_t ind_loc, uint8_t ind_old, uint8_t ind_new)
{
    if (true) {
        _nocheck__trace_virtio_ccw_set_ind(ind_loc, ind_old, ind_new);
    }
}

#define TRACE_S390_PCI_CLP_CAP_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_S390_PCI_CLP_CAP) || \
    false)

static inline void _nocheck__trace_s390_pci_clp_cap(const char * id, uint32_t cap)
{
    if (trace_event_get_state(TRACE_S390_PCI_CLP_CAP) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 19 "/FEMU/hw/s390x/trace-events"
            qemu_log("%d@%zu.%06zu:s390_pci_clp_cap " "PCI: %s: missing expected CLP capability %u" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , id, cap);
#line 411 "trace/trace-hw_s390x.h"
        } else {
#line 19 "/FEMU/hw/s390x/trace-events"
            qemu_log("s390_pci_clp_cap " "PCI: %s: missing expected CLP capability %u" "\n", id, cap);
#line 415 "trace/trace-hw_s390x.h"
        }
    }
}

static inline void trace_s390_pci_clp_cap(const char * id, uint32_t cap)
{
    if (true) {
        _nocheck__trace_s390_pci_clp_cap(id, cap);
    }
}

#define TRACE_S390_PCI_CLP_CAP_SIZE_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_S390_PCI_CLP_CAP_SIZE) || \
    false)

static inline void _nocheck__trace_s390_pci_clp_cap_size(const char * id, uint32_t size, uint32_t cap)
{
    if (trace_event_get_state(TRACE_S390_PCI_CLP_CAP_SIZE) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 20 "/FEMU/hw/s390x/trace-events"
            qemu_log("%d@%zu.%06zu:s390_pci_clp_cap_size " "PCI: %s: bad size (%u) for CLP capability %u" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , id, size, cap);
#line 442 "trace/trace-hw_s390x.h"
        } else {
#line 20 "/FEMU/hw/s390x/trace-events"
            qemu_log("s390_pci_clp_cap_size " "PCI: %s: bad size (%u) for CLP capability %u" "\n", id, size, cap);
#line 446 "trace/trace-hw_s390x.h"
        }
    }
}

static inline void trace_s390_pci_clp_cap_size(const char * id, uint32_t size, uint32_t cap)
{
    if (true) {
        _nocheck__trace_s390_pci_clp_cap_size(id, size, cap);
    }
}

#define TRACE_S390_PCI_CLP_DEV_INFO_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_S390_PCI_CLP_DEV_INFO) || \
    false)

static inline void _nocheck__trace_s390_pci_clp_dev_info(const char * id)
{
    if (trace_event_get_state(TRACE_S390_PCI_CLP_DEV_INFO) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 21 "/FEMU/hw/s390x/trace-events"
            qemu_log("%d@%zu.%06zu:s390_pci_clp_dev_info " "PCI: %s: cannot read vfio device info" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , id);
#line 473 "trace/trace-hw_s390x.h"
        } else {
#line 21 "/FEMU/hw/s390x/trace-events"
            qemu_log("s390_pci_clp_dev_info " "PCI: %s: cannot read vfio device info" "\n", id);
#line 477 "trace/trace-hw_s390x.h"
        }
    }
}

static inline void trace_s390_pci_clp_dev_info(const char * id)
{
    if (true) {
        _nocheck__trace_s390_pci_clp_dev_info(id);
    }
}
#endif /* TRACE_HW_S390X_GENERATED_TRACERS_H */
