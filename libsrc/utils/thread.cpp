
#include "thread.hpp"
#include "errors.hpp"

namespace prl {

#ifdef HAVE_THREADS

////////////////////////////////////////////////////////////////////////////////
/// Thread driver
////////////////////////////////////////////////////////////////////////////////
void * thread_driver(void * data) {
	auto & thread = *static_cast<thread_t*>(data);

  //signal start
  thread.mutex_lock();
  thread.set_idle();
  thread.cond_signal();
  thread.mutex_unlock();

  while(1) 
  {
    // go to sleep
    thread.mutex_lock();
    while (!thread.has_work() && !thread.has_stop())
      thread.cond_wait();

    if (thread.has_stop()) {
      thread.mutex_unlock();
      break;
    }
      
    thread.mutex_unlock();
    
    thread.execute_task();

    // signal done
    thread.mutex_lock();
    thread.set_idle();
    thread.cond_signal();
    thread.mutex_unlock();
  }

  // signal stopping
  thread.mutex_lock();
  thread.set_idle();
  thread.cond_signal();
  thread.mutex_unlock();

	pthread_exit(nullptr);
	return nullptr;
}

#endif


////////////////////////////////////////////////////////////////////////////////
/// Thread destructor
////////////////////////////////////////////////////////////////////////////////
thread_t::~thread_t()
{ 
#ifdef HAVE_THREADS

  // wait for thread to finish work
  mutex_lock();
  while (has_work()) cond_wait();

  // signal stop
  status_ = thread_status_t::stop;
  cond_signal();
  mutex_unlock();

  // wait for thread to stop
  mutex_lock();
  while (has_work()) cond_wait();
  mutex_unlock();

  pthread_join(thread_, nullptr);

  pthread_attr_destroy(&attr_);
  pthread_mutex_destroy(&mutex_);
  pthread_cond_destroy(&cond_);

#endif 
}

////////////////////////////////////////////////////////////////////////////////
/// process 
////////////////////////////////////////////////////////////////////////////////
void thread_t::launch() {
#ifdef HAVE_THREADS

  // wait for thread to be free
  mutex_lock();
  while (has_work())
    cond_wait();
  
  // signal thread there is work
  status_ = thread_status_t::busy;
  cond_signal();
  mutex_unlock();
  
#else

  task_();

#endif
}

////////////////////////////////////////////////////////////////////////////////
/// Wait for a thread to finish
////////////////////////////////////////////////////////////////////////////////
void thread_t::wait() {
#ifdef HAVE_THREADS

  // wait for threead to be finished
  mutex_lock();
  while (has_work())
    cond_wait();
  mutex_unlock();
  
#else

  wait_();

#endif
}

////////////////////////////////////////////////////////////////////////////////
/// start a thread
////////////////////////////////////////////////////////////////////////////////
void thread_t::start() {
#ifdef HAVE_THREADS

  if (pthread_mutex_init(&mutex_, NULL) != 0)
    THROW_ERROR("pthread_mutex_init() error");
  
  if (pthread_cond_init(&cond_, NULL) != 0)
  	THROW_ERROR("pthread_cond_init() error");

  if (pthread_attr_init(&attr_) != 0)
    THROW_ERROR("pthread_attr_init() error");

  if (pthread_attr_setdetachstate(&attr_, PTHREAD_CREATE_JOINABLE) != 0)
    THROW_ERROR("pthread_attr_setdetachstate() error");

  auto t = pthread_create(&thread_, &attr_, thread_driver, (void*)this);
  if (t != 0)
  	THROW_ERROR("pthread_create() error");

  // wait for thread to signal start
  mutex_lock();
  while (!has_idle())
    cond_wait();
  mutex_unlock();
#endif
}

} // namespace
