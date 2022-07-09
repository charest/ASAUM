//! system includes
#include <time.h>
#include <sys/time.h>

namespace prl {

////////////////////////////////////////////////////////////////////////////////
//! \brief Return the wall time counter.
//! \return the wall time
////////////////////////////////////////////////////////////////////////////////
double wall_time(void)
{
  struct timeval time;
  if (gettimeofday(&time,NULL)){
      //  Handle error
      return 0;
  }
  return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

}
