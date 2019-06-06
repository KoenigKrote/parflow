#ifndef _FLUX_RICHARDS_H_
#define _FLUX_RICHARDS_H_

#define Richards_Flux_Contrib \
  RJE_Flux_Prologue,                            \
    RJE_Flux_Epilogue,                          \
    RJE_Flux_Contrib_Left,                      \
    RJE_Flux_Contrib_Right,                     \
    RJE_Flux_Contrib_Up,                        \
    RJE_Flux_Contrib_Down,                      \
    RJE_Flux_Contrib_Front,                     \
    RJE_Flux_Contrib_Back

#define RJE_Flux_Prologue                       \
  {                                             \
    im = SubmatrixEltIndex(J_sub, i, j, k);     \
  }

#define RJE_Flux_Epilogue                       \
  {                                             \
    cp[im] += op[im];                           \
    op[im] = 0.0;                               \
  }

#define RJE_Flux_Contrib_Left                   \
  FACE(Left, { op = wp; })
#define RJE_Flux_Contrib_Right                  \
  FACE(Right, { op = ep; })
#define RJE_Flux_Contrib_Up                     \
  FACE(Up, { op = sop; })
#define RJE_Flux_Contrib_Down                   \
  FACE(Down, { op = np; })
#define RJE_Flux_Contrib_Front                  \
  FACE(Front, { op = lp; })
#define RJE_Flux_Contrib_Back                   \
  FACE(Back, { op = up; })

#endif // _FLUX_RICHARDS_H_
