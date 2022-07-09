#ifndef UTILS_MARK_HPP
#define UTILS_MARK_HPP 

#include "config.hpp"
#include "utils/trace.hpp"

#ifdef ENABLE_TRACE

#define MARK_FUNCTION() prl::scoped_trace_t _trace_(__func__)
#define MARK_REGION(name) prl::scoped_trace_t _trace_(name)
#define MARK_REGION_BEGIN(name) prl::tracer_t::instance().start(name)
#define MARK_REGION_END(name) prl::tracer_t::instance().stop(name)
#define DUMP_TRACE(...) prl::tracer_t::instance().dump(__VA_ARGS__)

#else

#define MARK_FUNCTION()
#define MARK_REGION(name)
#define MARK_REGION_BEGIN(name)
#define MARK_REGION_END(name)
#define DUMP_TRACE(...) do { } while(false)

#endif

#endif
