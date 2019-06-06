#ifndef _OVERLAND_RICHARDS_H_
#define _OVERLAND_RICHARDS_H_

#define Richards_Overland_Contrib_Default       \
  Richards_Overland_Contrib_Prologue,           \
    Richards_Overland_Contrib_Epilogue,         \
    RJE_Overland_Contrib_Left,                  \
    RJE_Overland_Contrib_Right,                 \
    RJE_Overland_Contrib_Up,                    \
    RJE_Overland_Contrib_Down,                  \
    RJE_Overland_Contrib_Front,                 \
    FACE(Back, RJE_Overland_Contrib_Back)


#define Richards_Overland_Contrib_Simple \
  Richards_Overland_Contrib_Prologue,    \
    Richards_Overland_Contrib_Epilogue,  \
    RJE_Overland_Contrib_Left,           \
    RJE_Overland_Contrib_Right,          \
    RJE_Overland_Contrib_Up,             \
    RJE_Overland_Contrib_Down,           \
    RJE_Overland_Contrib_Front,          \
    FACE(Back, RJE_Overland_Contrib_Back_Simple)

#define Richards_Overland_Contrib_Spinup        \
  Richards_Overland_Contrib_Prologue,           \
    Richards_Overland_Contrib_Epilogue,         \
    RJE_Overland_Contrib_Left,                  \
    RJE_Overland_Contrib_Right,                 \
    RJE_Overland_Contrib_Up,                    \
    RJE_Overland_Contrib_Down,                  \
    RJE_Overland_Contrib_Front,                 \
    FACE(Back, RJE_Overland_Contrib_Back_Spinup)

#define Richards_Overland_Contrib_Prologue      \
  {                                             \
    im = SubmatrixEltIndex(J_sub, i, j, k);     \
  }

#define Richards_Overland_Contrib_Epilogue      \
  {                                             \
    cp[im] += op[im];                           \
    op[im] = 0.0;                               \
  }

#define RJE_Overland_Contrib_Left                       \
  FACE(Left, { op = wp; })
#define RJE_Overland_Contrib_Right                      \
  FACE(Right, { op = ep; })
#define RJE_Overland_Contrib_Up                         \
  FACE(Up, { op = sop; })
#define RJE_Overland_Contrib_Down                       \
  FACE(Down, { op = np; })
#define RJE_Overland_Contrib_Front                      \
  FACE(Front, { op = lp; })
#define RJE_Overland_Contrib_Back                       \
  { op = wp; }

#define RJE_Overland_Contrib_Back_Simple                  \
  RJE_Overland_Contrib_Back                               \
  {                                                       \
    ip = SubvectorEltIndex(p_sub, i, j, k);               \
    if (pp[ip] > 0.0) {                                   \
      cp[im] += (vol * z_mult_dat[ip]) /                  \
        (dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v]))  \
        * (dt + 1);                                         \
    }                                                       \
  }

#define RJE_Overland_Contrib_Back_Spinup        \
  RJE_Overland_Contrib_Back                     \
  {                                             \
    ip = SubvectorEltIndex(p_sub, i, j, k);     \
    vol = dx * dy * dz;                         \
    if (pp[ip] >= 0.0) {                        \
      cp[im] += (vol / dz) * dt * (1.0 + 0.0);  \
    }                                           \
  }

/*
{                                               \
    if (diffusive == 0) {                       \
      PFModuleInvokeType(OverlandFlowEvalInvoke, overlandflow_module,   \
                         (grid, is, bc_struct, ipatch, problem_data, pressure, \
                          ke_der, kw_der, kn_der, ks_der, NULL, NULL, CALCDER)); \
    } else {                                                            \
      PFModuleInvokeType(OverlandFlowEvalDiffInvoke, overlandflow_module_diff, \
                         (grid, is, bc_struct, ipatch, problem_data, pressure, \
                          ke_der, kw_der, kn_der, ks_der,               \
                          kens_der, kwns_der, knns_der, ksns_der, NULL, NULL, CALCDER)); \
    }                                                                   \
  }
*/
#endif // _OVERLAND_RICHARDS_H_
