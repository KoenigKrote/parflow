#ifndef _TIMER_H
#define _TIMER_H

#include <time.h>

#define SETTIME(tp)                             \
  gettimeofday(&tp, NULL);

#define GETTIME(tval)                           \
  (double)tval.tv_sec + (double)tval.tv_usec * 1.e-6;

#define TIMER_LOOP(itr, fp, body)               \
  struct timeval tp;                            \
  gettimeofday(&tp, NULL);                      \
  double start_time = GETTIME(tp);              \
  for (int ii = 0; ii < itr; ii++)              \
    body;                                       \
                                                \
  gettimeofday(&tp, NULL);                          \
  double end_time = GETTIME(tp);                    \
  fprintf(fp, "%d,%lf,%lf\n",                       \
          itr, end_time - start_time,               \
          (end_time - start_time) / (double)(itr));

#define TIMER(fp, tp, start, end, body)         \
  SETTIME(tp);                                  \
  start = GETTIME(tp);                          \
  body;                                         \
  SETTIME(tp);                                  \
  end = GETTIME(tp);                            \
  fprintf(fp, "%lf\n", end - start);

#endif // _TIMER_H
