/* This file is autogenerated by tracetool, do not edit. */

#ifndef TRACE_HW_PCI_GENERATED_TRACERS_H
#define TRACE_HW_PCI_GENERATED_TRACERS_H

#include "trace/control.h"

extern TraceEvent _TRACE_PCI_UPDATE_MAPPINGS_DEL_EVENT;
extern TraceEvent _TRACE_PCI_UPDATE_MAPPINGS_ADD_EVENT;
extern TraceEvent _TRACE_PCI_ROUTE_IRQ_EVENT;
extern TraceEvent _TRACE_PCI_CFG_READ_EVENT;
extern TraceEvent _TRACE_PCI_CFG_WRITE_EVENT;
extern TraceEvent _TRACE_MSIX_WRITE_CONFIG_EVENT;
extern TraceEvent _TRACE_SRIOV_REGISTER_VFS_EVENT;
extern TraceEvent _TRACE_SRIOV_UNREGISTER_VFS_EVENT;
extern TraceEvent _TRACE_SRIOV_CONFIG_WRITE_EVENT;
extern uint16_t _TRACE_PCI_UPDATE_MAPPINGS_DEL_DSTATE;
extern uint16_t _TRACE_PCI_UPDATE_MAPPINGS_ADD_DSTATE;
extern uint16_t _TRACE_PCI_ROUTE_IRQ_DSTATE;
extern uint16_t _TRACE_PCI_CFG_READ_DSTATE;
extern uint16_t _TRACE_PCI_CFG_WRITE_DSTATE;
extern uint16_t _TRACE_MSIX_WRITE_CONFIG_DSTATE;
extern uint16_t _TRACE_SRIOV_REGISTER_VFS_DSTATE;
extern uint16_t _TRACE_SRIOV_UNREGISTER_VFS_DSTATE;
extern uint16_t _TRACE_SRIOV_CONFIG_WRITE_DSTATE;
#define TRACE_PCI_UPDATE_MAPPINGS_DEL_ENABLED 1
#define TRACE_PCI_UPDATE_MAPPINGS_ADD_ENABLED 1
#define TRACE_PCI_ROUTE_IRQ_ENABLED 1
#define TRACE_PCI_CFG_READ_ENABLED 1
#define TRACE_PCI_CFG_WRITE_ENABLED 1
#define TRACE_MSIX_WRITE_CONFIG_ENABLED 1
#define TRACE_SRIOV_REGISTER_VFS_ENABLED 1
#define TRACE_SRIOV_UNREGISTER_VFS_ENABLED 1
#define TRACE_SRIOV_CONFIG_WRITE_ENABLED 1
#include "qemu/log-for-trace.h"
#include "qemu/error-report.h"


#define TRACE_PCI_UPDATE_MAPPINGS_DEL_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_PCI_UPDATE_MAPPINGS_DEL) || \
    false)

static inline void _nocheck__trace_pci_update_mappings_del(const char * dev, uint32_t bus, uint32_t slot, uint32_t func, int bar, uint64_t addr, uint64_t size)
{
    if (trace_event_get_state(TRACE_PCI_UPDATE_MAPPINGS_DEL) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 4 "/FEMU/hw/pci/trace-events"
            qemu_log("%d@%zu.%06zu:pci_update_mappings_del " "%s %02x:%02x.%x %d,0x%"PRIx64"+0x%"PRIx64 "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , dev, bus, slot, func, bar, addr, size);
#line 55 "trace/trace-hw_pci.h"
        } else {
#line 4 "/FEMU/hw/pci/trace-events"
            qemu_log("pci_update_mappings_del " "%s %02x:%02x.%x %d,0x%"PRIx64"+0x%"PRIx64 "\n", dev, bus, slot, func, bar, addr, size);
#line 59 "trace/trace-hw_pci.h"
        }
    }
}

static inline void trace_pci_update_mappings_del(const char * dev, uint32_t bus, uint32_t slot, uint32_t func, int bar, uint64_t addr, uint64_t size)
{
    if (true) {
        _nocheck__trace_pci_update_mappings_del(dev, bus, slot, func, bar, addr, size);
    }
}

#define TRACE_PCI_UPDATE_MAPPINGS_ADD_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_PCI_UPDATE_MAPPINGS_ADD) || \
    false)

static inline void _nocheck__trace_pci_update_mappings_add(const char * dev, uint32_t bus, uint32_t slot, uint32_t func, int bar, uint64_t addr, uint64_t size)
{
    if (trace_event_get_state(TRACE_PCI_UPDATE_MAPPINGS_ADD) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 5 "/FEMU/hw/pci/trace-events"
            qemu_log("%d@%zu.%06zu:pci_update_mappings_add " "%s %02x:%02x.%x %d,0x%"PRIx64"+0x%"PRIx64 "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , dev, bus, slot, func, bar, addr, size);
#line 86 "trace/trace-hw_pci.h"
        } else {
#line 5 "/FEMU/hw/pci/trace-events"
            qemu_log("pci_update_mappings_add " "%s %02x:%02x.%x %d,0x%"PRIx64"+0x%"PRIx64 "\n", dev, bus, slot, func, bar, addr, size);
#line 90 "trace/trace-hw_pci.h"
        }
    }
}

static inline void trace_pci_update_mappings_add(const char * dev, uint32_t bus, uint32_t slot, uint32_t func, int bar, uint64_t addr, uint64_t size)
{
    if (true) {
        _nocheck__trace_pci_update_mappings_add(dev, bus, slot, func, bar, addr, size);
    }
}

#define TRACE_PCI_ROUTE_IRQ_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_PCI_ROUTE_IRQ) || \
    false)

static inline void _nocheck__trace_pci_route_irq(int dev_irq, const char * dev_path, int parent_irq, const char * parent_path)
{
    if (trace_event_get_state(TRACE_PCI_ROUTE_IRQ) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 6 "/FEMU/hw/pci/trace-events"
            qemu_log("%d@%zu.%06zu:pci_route_irq " "IRQ %d @%s -> IRQ %d @%s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , dev_irq, dev_path, parent_irq, parent_path);
#line 117 "trace/trace-hw_pci.h"
        } else {
#line 6 "/FEMU/hw/pci/trace-events"
            qemu_log("pci_route_irq " "IRQ %d @%s -> IRQ %d @%s" "\n", dev_irq, dev_path, parent_irq, parent_path);
#line 121 "trace/trace-hw_pci.h"
        }
    }
}

static inline void trace_pci_route_irq(int dev_irq, const char * dev_path, int parent_irq, const char * parent_path)
{
    if (true) {
        _nocheck__trace_pci_route_irq(dev_irq, dev_path, parent_irq, parent_path);
    }
}

#define TRACE_PCI_CFG_READ_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_PCI_CFG_READ) || \
    false)

static inline void _nocheck__trace_pci_cfg_read(const char * dev, uint32_t bus, uint32_t slot, uint32_t func, unsigned offs, unsigned val)
{
    if (trace_event_get_state(TRACE_PCI_CFG_READ) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 9 "/FEMU/hw/pci/trace-events"
            qemu_log("%d@%zu.%06zu:pci_cfg_read " "%s %02x:%02x.%x @0x%x -> 0x%x" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , dev, bus, slot, func, offs, val);
#line 148 "trace/trace-hw_pci.h"
        } else {
#line 9 "/FEMU/hw/pci/trace-events"
            qemu_log("pci_cfg_read " "%s %02x:%02x.%x @0x%x -> 0x%x" "\n", dev, bus, slot, func, offs, val);
#line 152 "trace/trace-hw_pci.h"
        }
    }
}

static inline void trace_pci_cfg_read(const char * dev, uint32_t bus, uint32_t slot, uint32_t func, unsigned offs, unsigned val)
{
    if (true) {
        _nocheck__trace_pci_cfg_read(dev, bus, slot, func, offs, val);
    }
}

#define TRACE_PCI_CFG_WRITE_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_PCI_CFG_WRITE) || \
    false)

static inline void _nocheck__trace_pci_cfg_write(const char * dev, uint32_t bus, uint32_t slot, uint32_t func, unsigned offs, unsigned val)
{
    if (trace_event_get_state(TRACE_PCI_CFG_WRITE) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 10 "/FEMU/hw/pci/trace-events"
            qemu_log("%d@%zu.%06zu:pci_cfg_write " "%s %02x:%02x.%x @0x%x <- 0x%x" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , dev, bus, slot, func, offs, val);
#line 179 "trace/trace-hw_pci.h"
        } else {
#line 10 "/FEMU/hw/pci/trace-events"
            qemu_log("pci_cfg_write " "%s %02x:%02x.%x @0x%x <- 0x%x" "\n", dev, bus, slot, func, offs, val);
#line 183 "trace/trace-hw_pci.h"
        }
    }
}

static inline void trace_pci_cfg_write(const char * dev, uint32_t bus, uint32_t slot, uint32_t func, unsigned offs, unsigned val)
{
    if (true) {
        _nocheck__trace_pci_cfg_write(dev, bus, slot, func, offs, val);
    }
}

#define TRACE_MSIX_WRITE_CONFIG_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_MSIX_WRITE_CONFIG) || \
    false)

static inline void _nocheck__trace_msix_write_config(char * name, bool enabled, bool masked)
{
    if (trace_event_get_state(TRACE_MSIX_WRITE_CONFIG) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 13 "/FEMU/hw/pci/trace-events"
            qemu_log("%d@%zu.%06zu:msix_write_config " "dev %s enabled %d masked %d" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , name, enabled, masked);
#line 210 "trace/trace-hw_pci.h"
        } else {
#line 13 "/FEMU/hw/pci/trace-events"
            qemu_log("msix_write_config " "dev %s enabled %d masked %d" "\n", name, enabled, masked);
#line 214 "trace/trace-hw_pci.h"
        }
    }
}

static inline void trace_msix_write_config(char * name, bool enabled, bool masked)
{
    if (true) {
        _nocheck__trace_msix_write_config(name, enabled, masked);
    }
}

#define TRACE_SRIOV_REGISTER_VFS_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_SRIOV_REGISTER_VFS) || \
    false)

static inline void _nocheck__trace_sriov_register_vfs(const char * name, int slot, int function, int num_vfs)
{
    if (trace_event_get_state(TRACE_SRIOV_REGISTER_VFS) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 16 "/FEMU/hw/pci/trace-events"
            qemu_log("%d@%zu.%06zu:sriov_register_vfs " "%s %02x:%x: creating %d vf devs" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , name, slot, function, num_vfs);
#line 241 "trace/trace-hw_pci.h"
        } else {
#line 16 "/FEMU/hw/pci/trace-events"
            qemu_log("sriov_register_vfs " "%s %02x:%x: creating %d vf devs" "\n", name, slot, function, num_vfs);
#line 245 "trace/trace-hw_pci.h"
        }
    }
}

static inline void trace_sriov_register_vfs(const char * name, int slot, int function, int num_vfs)
{
    if (true) {
        _nocheck__trace_sriov_register_vfs(name, slot, function, num_vfs);
    }
}

#define TRACE_SRIOV_UNREGISTER_VFS_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_SRIOV_UNREGISTER_VFS) || \
    false)

static inline void _nocheck__trace_sriov_unregister_vfs(const char * name, int slot, int function, int num_vfs)
{
    if (trace_event_get_state(TRACE_SRIOV_UNREGISTER_VFS) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 17 "/FEMU/hw/pci/trace-events"
            qemu_log("%d@%zu.%06zu:sriov_unregister_vfs " "%s %02x:%x: Unregistering %d vf devs" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , name, slot, function, num_vfs);
#line 272 "trace/trace-hw_pci.h"
        } else {
#line 17 "/FEMU/hw/pci/trace-events"
            qemu_log("sriov_unregister_vfs " "%s %02x:%x: Unregistering %d vf devs" "\n", name, slot, function, num_vfs);
#line 276 "trace/trace-hw_pci.h"
        }
    }
}

static inline void trace_sriov_unregister_vfs(const char * name, int slot, int function, int num_vfs)
{
    if (true) {
        _nocheck__trace_sriov_unregister_vfs(name, slot, function, num_vfs);
    }
}

#define TRACE_SRIOV_CONFIG_WRITE_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_SRIOV_CONFIG_WRITE) || \
    false)

static inline void _nocheck__trace_sriov_config_write(const char * name, int slot, int fun, uint32_t offset, uint32_t val, uint32_t len)
{
    if (trace_event_get_state(TRACE_SRIOV_CONFIG_WRITE) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 18 "/FEMU/hw/pci/trace-events"
            qemu_log("%d@%zu.%06zu:sriov_config_write " "%s %02x:%x: sriov offset 0x%x val 0x%x len %d" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , name, slot, fun, offset, val, len);
#line 303 "trace/trace-hw_pci.h"
        } else {
#line 18 "/FEMU/hw/pci/trace-events"
            qemu_log("sriov_config_write " "%s %02x:%x: sriov offset 0x%x val 0x%x len %d" "\n", name, slot, fun, offset, val, len);
#line 307 "trace/trace-hw_pci.h"
        }
    }
}

static inline void trace_sriov_config_write(const char * name, int slot, int fun, uint32_t offset, uint32_t val, uint32_t len)
{
    if (true) {
        _nocheck__trace_sriov_config_write(name, slot, fun, offset, val, len);
    }
}
#endif /* TRACE_HW_PCI_GENERATED_TRACERS_H */
