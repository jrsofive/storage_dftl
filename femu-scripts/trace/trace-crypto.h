/* This file is autogenerated by tracetool, do not edit. */

#ifndef TRACE_CRYPTO_GENERATED_TRACERS_H
#define TRACE_CRYPTO_GENERATED_TRACERS_H

#include "trace/control.h"

extern TraceEvent _TRACE_QCRYPTO_TLS_CREDS_LOAD_DH_EVENT;
extern TraceEvent _TRACE_QCRYPTO_TLS_CREDS_GET_PATH_EVENT;
extern TraceEvent _TRACE_QCRYPTO_TLS_CREDS_ANON_LOAD_EVENT;
extern TraceEvent _TRACE_QCRYPTO_TLS_CREDS_PSK_LOAD_EVENT;
extern TraceEvent _TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_EVENT;
extern TraceEvent _TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_BASIC_CONSTRAINTS_EVENT;
extern TraceEvent _TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_KEY_USAGE_EVENT;
extern TraceEvent _TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_KEY_PURPOSE_EVENT;
extern TraceEvent _TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_CERT_EVENT;
extern TraceEvent _TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_CERT_LIST_EVENT;
extern TraceEvent _TRACE_QCRYPTO_TLS_SESSION_NEW_EVENT;
extern TraceEvent _TRACE_QCRYPTO_TLS_SESSION_CHECK_CREDS_EVENT;
extern TraceEvent _TRACE_QCRYPTO_TLS_CIPHER_SUITE_PRIORITY_EVENT;
extern TraceEvent _TRACE_QCRYPTO_TLS_CIPHER_SUITE_INFO_EVENT;
extern TraceEvent _TRACE_QCRYPTO_TLS_CIPHER_SUITE_COUNT_EVENT;
extern uint16_t _TRACE_QCRYPTO_TLS_CREDS_LOAD_DH_DSTATE;
extern uint16_t _TRACE_QCRYPTO_TLS_CREDS_GET_PATH_DSTATE;
extern uint16_t _TRACE_QCRYPTO_TLS_CREDS_ANON_LOAD_DSTATE;
extern uint16_t _TRACE_QCRYPTO_TLS_CREDS_PSK_LOAD_DSTATE;
extern uint16_t _TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_DSTATE;
extern uint16_t _TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_BASIC_CONSTRAINTS_DSTATE;
extern uint16_t _TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_KEY_USAGE_DSTATE;
extern uint16_t _TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_KEY_PURPOSE_DSTATE;
extern uint16_t _TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_CERT_DSTATE;
extern uint16_t _TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_CERT_LIST_DSTATE;
extern uint16_t _TRACE_QCRYPTO_TLS_SESSION_NEW_DSTATE;
extern uint16_t _TRACE_QCRYPTO_TLS_SESSION_CHECK_CREDS_DSTATE;
extern uint16_t _TRACE_QCRYPTO_TLS_CIPHER_SUITE_PRIORITY_DSTATE;
extern uint16_t _TRACE_QCRYPTO_TLS_CIPHER_SUITE_INFO_DSTATE;
extern uint16_t _TRACE_QCRYPTO_TLS_CIPHER_SUITE_COUNT_DSTATE;
#define TRACE_QCRYPTO_TLS_CREDS_LOAD_DH_ENABLED 1
#define TRACE_QCRYPTO_TLS_CREDS_GET_PATH_ENABLED 1
#define TRACE_QCRYPTO_TLS_CREDS_ANON_LOAD_ENABLED 1
#define TRACE_QCRYPTO_TLS_CREDS_PSK_LOAD_ENABLED 1
#define TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_ENABLED 1
#define TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_BASIC_CONSTRAINTS_ENABLED 1
#define TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_KEY_USAGE_ENABLED 1
#define TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_KEY_PURPOSE_ENABLED 1
#define TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_CERT_ENABLED 1
#define TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_CERT_LIST_ENABLED 1
#define TRACE_QCRYPTO_TLS_SESSION_NEW_ENABLED 1
#define TRACE_QCRYPTO_TLS_SESSION_CHECK_CREDS_ENABLED 1
#define TRACE_QCRYPTO_TLS_CIPHER_SUITE_PRIORITY_ENABLED 1
#define TRACE_QCRYPTO_TLS_CIPHER_SUITE_INFO_ENABLED 1
#define TRACE_QCRYPTO_TLS_CIPHER_SUITE_COUNT_ENABLED 1
#include "qemu/log-for-trace.h"
#include "qemu/error-report.h"


#define TRACE_QCRYPTO_TLS_CREDS_LOAD_DH_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_CREDS_LOAD_DH) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_creds_load_dh(void * creds, const char * filename)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_CREDS_LOAD_DH) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 4 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_creds_load_dh " "TLS creds load DH creds=%p filename=%s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , creds, filename);
#line 73 "trace/trace-crypto.h"
        } else {
#line 4 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_creds_load_dh " "TLS creds load DH creds=%p filename=%s" "\n", creds, filename);
#line 77 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_creds_load_dh(void * creds, const char * filename)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_creds_load_dh(creds, filename);
    }
}

#define TRACE_QCRYPTO_TLS_CREDS_GET_PATH_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_CREDS_GET_PATH) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_creds_get_path(void * creds, const char * filename, const char * path)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_CREDS_GET_PATH) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 5 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_creds_get_path " "TLS creds path creds=%p filename=%s path=%s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , creds, filename, path);
#line 104 "trace/trace-crypto.h"
        } else {
#line 5 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_creds_get_path " "TLS creds path creds=%p filename=%s path=%s" "\n", creds, filename, path);
#line 108 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_creds_get_path(void * creds, const char * filename, const char * path)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_creds_get_path(creds, filename, path);
    }
}

#define TRACE_QCRYPTO_TLS_CREDS_ANON_LOAD_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_CREDS_ANON_LOAD) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_creds_anon_load(void * creds, const char * dir)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_CREDS_ANON_LOAD) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 8 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_creds_anon_load " "TLS creds anon load creds=%p dir=%s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , creds, dir);
#line 135 "trace/trace-crypto.h"
        } else {
#line 8 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_creds_anon_load " "TLS creds anon load creds=%p dir=%s" "\n", creds, dir);
#line 139 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_creds_anon_load(void * creds, const char * dir)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_creds_anon_load(creds, dir);
    }
}

#define TRACE_QCRYPTO_TLS_CREDS_PSK_LOAD_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_CREDS_PSK_LOAD) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_creds_psk_load(void * creds, const char * dir)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_CREDS_PSK_LOAD) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 11 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_creds_psk_load " "TLS creds psk load creds=%p dir=%s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , creds, dir);
#line 166 "trace/trace-crypto.h"
        } else {
#line 11 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_creds_psk_load " "TLS creds psk load creds=%p dir=%s" "\n", creds, dir);
#line 170 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_creds_psk_load(void * creds, const char * dir)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_creds_psk_load(creds, dir);
    }
}

#define TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_CREDS_X509_LOAD) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_creds_x509_load(void * creds, const char * dir)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_CREDS_X509_LOAD) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 14 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_creds_x509_load " "TLS creds x509 load creds=%p dir=%s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , creds, dir);
#line 197 "trace/trace-crypto.h"
        } else {
#line 14 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_creds_x509_load " "TLS creds x509 load creds=%p dir=%s" "\n", creds, dir);
#line 201 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_creds_x509_load(void * creds, const char * dir)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_creds_x509_load(creds, dir);
    }
}

#define TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_BASIC_CONSTRAINTS_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_BASIC_CONSTRAINTS) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_creds_x509_check_basic_constraints(void * creds, const char * file, int status)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_BASIC_CONSTRAINTS) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 15 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_creds_x509_check_basic_constraints " "TLS creds x509 check basic constraints creds=%p file=%s status=%d" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , creds, file, status);
#line 228 "trace/trace-crypto.h"
        } else {
#line 15 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_creds_x509_check_basic_constraints " "TLS creds x509 check basic constraints creds=%p file=%s status=%d" "\n", creds, file, status);
#line 232 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_creds_x509_check_basic_constraints(void * creds, const char * file, int status)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_creds_x509_check_basic_constraints(creds, file, status);
    }
}

#define TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_KEY_USAGE_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_KEY_USAGE) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_creds_x509_check_key_usage(void * creds, const char * file, int status, int usage, int critical)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_KEY_USAGE) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 16 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_creds_x509_check_key_usage " "TLS creds x509 check key usage creds=%p file=%s status=%d usage=%d critical=%d" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , creds, file, status, usage, critical);
#line 259 "trace/trace-crypto.h"
        } else {
#line 16 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_creds_x509_check_key_usage " "TLS creds x509 check key usage creds=%p file=%s status=%d usage=%d critical=%d" "\n", creds, file, status, usage, critical);
#line 263 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_creds_x509_check_key_usage(void * creds, const char * file, int status, int usage, int critical)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_creds_x509_check_key_usage(creds, file, status, usage, critical);
    }
}

#define TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_KEY_PURPOSE_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_KEY_PURPOSE) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_creds_x509_check_key_purpose(void * creds, const char * file, int status, const char * usage, int critical)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_CREDS_X509_CHECK_KEY_PURPOSE) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 17 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_creds_x509_check_key_purpose " "TLS creds x509 check key usage creds=%p file=%s status=%d usage=%s critical=%d" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , creds, file, status, usage, critical);
#line 290 "trace/trace-crypto.h"
        } else {
#line 17 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_creds_x509_check_key_purpose " "TLS creds x509 check key usage creds=%p file=%s status=%d usage=%s critical=%d" "\n", creds, file, status, usage, critical);
#line 294 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_creds_x509_check_key_purpose(void * creds, const char * file, int status, const char * usage, int critical)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_creds_x509_check_key_purpose(creds, file, status, usage, critical);
    }
}

#define TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_CERT_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_CERT) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_creds_x509_load_cert(void * creds, int isServer, const char * file)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_CERT) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 18 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_creds_x509_load_cert " "TLS creds x509 load cert creds=%p isServer=%d file=%s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , creds, isServer, file);
#line 321 "trace/trace-crypto.h"
        } else {
#line 18 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_creds_x509_load_cert " "TLS creds x509 load cert creds=%p isServer=%d file=%s" "\n", creds, isServer, file);
#line 325 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_creds_x509_load_cert(void * creds, int isServer, const char * file)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_creds_x509_load_cert(creds, isServer, file);
    }
}

#define TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_CERT_LIST_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_CERT_LIST) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_creds_x509_load_cert_list(void * creds, const char * file)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_CREDS_X509_LOAD_CERT_LIST) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 19 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_creds_x509_load_cert_list " "TLS creds x509 load cert list creds=%p file=%s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , creds, file);
#line 352 "trace/trace-crypto.h"
        } else {
#line 19 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_creds_x509_load_cert_list " "TLS creds x509 load cert list creds=%p file=%s" "\n", creds, file);
#line 356 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_creds_x509_load_cert_list(void * creds, const char * file)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_creds_x509_load_cert_list(creds, file);
    }
}

#define TRACE_QCRYPTO_TLS_SESSION_NEW_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_SESSION_NEW) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_session_new(void * session, void * creds, const char * hostname, const char * authzid, int endpoint)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_SESSION_NEW) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 22 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_session_new " "TLS session new session=%p creds=%p hostname=%s authzid=%s endpoint=%d" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , session, creds, hostname, authzid, endpoint);
#line 383 "trace/trace-crypto.h"
        } else {
#line 22 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_session_new " "TLS session new session=%p creds=%p hostname=%s authzid=%s endpoint=%d" "\n", session, creds, hostname, authzid, endpoint);
#line 387 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_session_new(void * session, void * creds, const char * hostname, const char * authzid, int endpoint)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_session_new(session, creds, hostname, authzid, endpoint);
    }
}

#define TRACE_QCRYPTO_TLS_SESSION_CHECK_CREDS_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_SESSION_CHECK_CREDS) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_session_check_creds(void * session, const char * status)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_SESSION_CHECK_CREDS) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 23 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_session_check_creds " "TLS session check creds session=%p status=%s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , session, status);
#line 414 "trace/trace-crypto.h"
        } else {
#line 23 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_session_check_creds " "TLS session check creds session=%p status=%s" "\n", session, status);
#line 418 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_session_check_creds(void * session, const char * status)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_session_check_creds(session, status);
    }
}

#define TRACE_QCRYPTO_TLS_CIPHER_SUITE_PRIORITY_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_CIPHER_SUITE_PRIORITY) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_cipher_suite_priority(const char * name)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_CIPHER_SUITE_PRIORITY) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 26 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_cipher_suite_priority " "priority: %s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , name);
#line 445 "trace/trace-crypto.h"
        } else {
#line 26 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_cipher_suite_priority " "priority: %s" "\n", name);
#line 449 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_cipher_suite_priority(const char * name)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_cipher_suite_priority(name);
    }
}

#define TRACE_QCRYPTO_TLS_CIPHER_SUITE_INFO_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_CIPHER_SUITE_INFO) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_cipher_suite_info(uint8_t data0, uint8_t data1, const char * version, const char * name)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_CIPHER_SUITE_INFO) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 27 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_cipher_suite_info " "data=[0x%02x,0x%02x] version=%s name=%s" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , data0, data1, version, name);
#line 476 "trace/trace-crypto.h"
        } else {
#line 27 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_cipher_suite_info " "data=[0x%02x,0x%02x] version=%s name=%s" "\n", data0, data1, version, name);
#line 480 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_cipher_suite_info(uint8_t data0, uint8_t data1, const char * version, const char * name)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_cipher_suite_info(data0, data1, version, name);
    }
}

#define TRACE_QCRYPTO_TLS_CIPHER_SUITE_COUNT_BACKEND_DSTATE() ( \
    trace_event_get_state_dynamic_by_id(TRACE_QCRYPTO_TLS_CIPHER_SUITE_COUNT) || \
    false)

static inline void _nocheck__trace_qcrypto_tls_cipher_suite_count(unsigned count)
{
    if (trace_event_get_state(TRACE_QCRYPTO_TLS_CIPHER_SUITE_COUNT) && qemu_loglevel_mask(LOG_TRACE)) {
        if (message_with_timestamp) {
            struct timeval _now;
            gettimeofday(&_now, NULL);
#line 28 "/FEMU/crypto/trace-events"
            qemu_log("%d@%zu.%06zu:qcrypto_tls_cipher_suite_count " "count: %u" "\n",
                     qemu_get_thread_id(),
                     (size_t)_now.tv_sec, (size_t)_now.tv_usec
                     , count);
#line 507 "trace/trace-crypto.h"
        } else {
#line 28 "/FEMU/crypto/trace-events"
            qemu_log("qcrypto_tls_cipher_suite_count " "count: %u" "\n", count);
#line 511 "trace/trace-crypto.h"
        }
    }
}

static inline void trace_qcrypto_tls_cipher_suite_count(unsigned count)
{
    if (true) {
        _nocheck__trace_qcrypto_tls_cipher_suite_count(count);
    }
}
#endif /* TRACE_CRYPTO_GENERATED_TRACERS_H */
