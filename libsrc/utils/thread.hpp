#ifndef THREAD_HPP
#define THREAD_HPP

#include "config.hpp"
#include "errors.hpp"

#ifdef HAVE_THREADS
#include "pthread.h"
#endif

#include<functional>

namespace prl {

void * thread_driver(void * data);

////////////////////////////////////////////////////////////////////////////////
/// Comm queue thread tasks
////////////////////////////////////////////////////////////////////////////////
enum class thread_status_t {
  idle,
  busy,
  stop
};

////////////////////////////////////////////////////////////////////////////////
/// Thread helper
////////////////////////////////////////////////////////////////////////////////
class thread_t {
  
  std::function<void(void)> task_;
  std::function<void(void)> wait_;

#ifdef HAVE_THREADS
  pthread_attr_t attr_;
  pthread_t thread_;
  pthread_mutex_t mutex_ = PTHREAD_MUTEX_INITIALIZER;
  pthread_cond_t cond_ = PTHREAD_COND_INITIALIZER;
  thread_status_t status_ = thread_status_t::stop;
  
public:

  bool has_work() const
  { return status_ == thread_status_t::busy; }
  bool has_stop() const
  { return status_ == thread_status_t::stop; }
  bool has_idle() const
  { return status_ == thread_status_t::idle; }
  
  void set_idle()
  {	status_ = thread_status_t::idle; }
	
  void mutex_lock() { pthread_mutex_lock(&mutex_);	}
  void mutex_unlock() { pthread_mutex_unlock(&mutex_);	}
  void cond_wait() { pthread_cond_wait(&cond_, &mutex_); }
  void cond_signal() { pthread_cond_signal(&cond_); }

  void execute_task() {
    task_();
    wait_();
  }

#endif

private:

  /// start a thread
  void start();

public:
  
  ~thread_t();

  void launch();
  void wait();

  /// Constructor
  template<typename Fun>
  thread_t(Fun f) : task_(f)
  { start(); }   
  
  template<typename TaskFun, typename WaitFun>
  thread_t(TaskFun t, WaitFun w) : task_(t), wait_(w)
  { start(); }   

};
  

} // namespace

#endif
