#ifndef _DIRICHLET_RICHARDS_H_
#define _DIRICHLET_RICHARDS_H_

#define Richards_DirichletBC_Contrib  \
  RJE_Dirichlet_Contrib_Prologue,     \
    RJE_Dirichlet_Contrib_Epilogue,   \
    RJE_Dirichlet_Contrib_Left,       \
    RJE_Dirichlet_Contrib_Right,      \
    RJE_Dirichlet_Contrib_Up,         \
    RJE_Dirichlet_Contrib_Down,       \
    RJE_Dirichlet_Contrib_Front,      \
    RJE_Dirichlet_Contrib_Back

#define RJE_Dirichlet_Contrib_Prologue                              \
  {                                                                 \
    value = bc_patch_values[ival];                                  \
    ip = SubvectorEltIndex(p_sub, i, j, k);                         \
    im = SubmatrixEltIndex(J_sub, i, j, k);                         \
    prod = rpp[ip] * dp[ip];                                        \
    prod_der = rpdp[ip] * dp[ip] + rpp[ip] * ddp[ip];               \
    PFModuleInvokeType(PhaseDensityInvoke, density_module,          \
                       (0, NULL, NULL, &value, &den_d, CALCFCN));   \
    PFModuleInvokeType(PhaseDensityInvoke, density_module,          \
                       (0, NULL, NULL, &value, &dend_d, CALCDER));  \
  }

#define RJE_Dirichlet_Contrib_Epilogue          \
  {                                             \
    cp[im] += op[im];                           \
    cp[im] -= o_temp;                           \
    op[im] = 0.0;                               \
  }

#define RJE_Dirichlet_Contrib_Left                          \
  FACE(Left,                                                \
       {                                                    \
         op = wp;                                           \
         coeff = dt * ffx * z_mult_dat[ip]                  \
           * (2.0 / dx) * permxp[ip] / viscosity;           \
         prod_val = rpp[ip - 1] * den_d;                    \
         diff = value - pp[ip];                             \
         o_temp = coeff                                     \
           * (diff * RPMean(value, pp[ip], 0.0, prod_der)   \
              - RPMean(value, pp[ip], prod_val, prod));     \
       })

#define RJE_Dirichlet_Contrib_Right                           \
  FACE(Right,                                                 \
       {                                                      \
         op = ep;                                             \
         coeff = dt * ffx * z_mult_dat[ip]                    \
           * (2.0 / dx) * permxp[ip] / viscosity;             \
         prod_val = rpp[ip + 1] * den_d;                      \
         diff = pp[ip] - value;                               \
         o_temp = -coeff                                      \
           * (diff * RPMean(value, pp[ip], prod_der, 0.0)     \
              + RPMean(value, pp[ip], prod, prod_val));       \
       })

#define RJE_Dirichlet_Contrib_Up                             \
  FACE(Up,                                                   \
       {                                                     \
         op = sop;                                           \
         coeff = dt * ffy * z_mult_dat[ip]                   \
           * (2.0 / dy) * permyp[ip] / viscosity;            \
         prod_val = rpp[ip - sy_v] * den_d;                  \
         diff = value - pp[ip];                              \
         o_temp = coeff                                      \
           * (diff * RPMean(value, pp[ip], 0.0, prod_der)    \
              - RPMean(value, pp[ip], prod_val, prod));      \
       })

#define RJE_Dirichlet_Contrib_Down                             \
  FACE(Down,                                                   \
       {                                                       \
       op = np;                                                \
       coeff = dt * ffx * z_mult_dat[ip]                       \
         * (2.0 / dx) * permxp[ip] / viscosity;                \
       prod_val = rpp[ip + sy_v] * den_d;                      \
       diff = pp[ip] - value;                                  \
       o_temp = -coeff                                         \
         * (diff * RPMean(value, pp[ip], prod_der, 0.0)        \
            + RPMean(value, pp[ip], prod, prod_val));          \
     })

#define RJE_Dirichlet_Contrib_Front                                     \
  FACE(Front,                                                           \
  {                                                                     \
    op = lp;                                                            \
    coeff = dt * ffz *                                                  \
      (2.0 / (dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])))        \
      * permzp[ip] / viscosity;                                         \
    prod_val = rpp[ip - sz_v] * den_d;                                  \
    lower_cond = value - 0.5 * dz * z_mult_dat[ip] * den_d * gravity;   \
    upper_cond = pp[ip] + 0.5 * dz * z_mult_dat[ip] * dp[ip] * gravity; \
    diff = lower_cond - upper_cond;                                     \
    o_temp = coeff * (diff * RPMean(lower_cond, upper_cond, 0.0, prod_der) \
                      + ((-1.0 - gravity * 0.5 * dz                     \
                          * Mean(z_mult_dat[ip], z_mult_dat[ip - sz_v]) * ddp[ip]) \
                         * RPMean(lower_cond, upper_cond, prod_val, prod))); \
  })

#define RJE_Dirichlet_Contrib_Back                                      \
  FACE(Back,                                                            \
  {                                                                     \
    op = up;                                                            \
    coeff = dt * ffz *                                                  \
      (2.0 / (dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])))        \
      * permzp[ip] / viscosity;                                         \
    prod_val = rpp[ip + sz_v] * den_d;                                  \
    lower_cond = pp[ip] - 0.5 * dz                                      \
      * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])                     \
      * dp[ip] * gravity;                                               \
    upper_cond = value + 0.5 * dz                                       \
      * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])                     \
      * den_d * gravity;                                                \
    diff = lower_cond - upper_cond;                                     \
    o_temp = -coeff * (diff * RPMean(lower_cond, upper_cond, prod_der, 0.0) \
                       + ((1.0 - gravity * 0.5 * dz                     \
                           * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v]) * ddp[ip]) \
                          * RPMean(lower_cond, upper_cond, prod, prod_val))); \
  })


#endif // _DIRICHLET_RICHARDS_H_
