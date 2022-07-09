#ifndef MPI_MUTEX_HPP
#define MPI_MUTEX_HPP

#include "config.hpp"

#include <stdio.h>

#ifdef HAVE_THREADS
#include "pthread.h"
#endif

namespace prl {

class comm_mutex_t {
  
#ifdef HAVE_THREADS
  pthread_mutex_t mutex_ = PTHREAD_MUTEX_INITIALIZER;
  
  comm_mutex_t() {
    if (pthread_mutex_init(&mutex_, NULL) != 0)
      perror("pthread_mutex_init() error");
  }
  
public:
  
  ~comm_mutex_t() { pthread_mutex_destroy(&mutex_); }


  void lock() { /*pthread_mutex_lock(&mutex_);*/ }
  void unlock() { /*pthread_mutex_unlock(&mutex_);*/ }

#else

public:

  void lock() const {}
  void unlock() const {}

#endif
  
  static comm_mutex_t& instance()
  {
    static comm_mutex_t s;
    return s;
  }

};
  
} // namespace

#endif
