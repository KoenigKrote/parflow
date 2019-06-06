#ifndef _RICHARDS_SYMM_CORRECTION_H_
#define _RICHARDS_SYMM_CORRECTION_H_

#define Richards_SymmCorrection(J_sub, p_sub, rpp, dp, i, j, k)     \
  Richards_SymmCorrection_Prologue(J_sub, p_sub, rpp, dp, i, j, k), \
    Richards_SymmCorrection_Epilogue,                               \
    Richards_SymmCorrection_Left,                                   \
    Richards_SymmCorrection_Right,                                  \
    Richards_SymmCorrection_Up,                                     \
    Richards_SymmCorrection_Down,                                   \
    Richards_SymmCorrection_Front,                                  \
    Richards_SymmCorrection_Back                                    \

#define Richards_SymmCorrection_Prologue(J_sub, p_sub, rpp, dp, i, j, k) \
  {                                                                   \
    ip = SubvectorEltIndex((p_sub), (i), (j), (k));                   \
    im = SubmatrixEltIndex((J_sub), (i), (j), (k));                   \
    prod = (rpp[ip]) * (dp[ip]);                                      \
  }

#define Richards_SymmCorrection_Epilogue {}

#define Richards_SymmCorrection_Left                          \
  FACE(Left,                                                  \
       {                                                      \
         diff = pp[ip - 1] - pp[ip];                          \
         prod_der = rpdp[ip - 1] * dp[ip - 1]                 \
           + rpp[ip - 1] * ddp[ip - 1];                       \
         coeff = dt * z_mult_dat[ip] * ffx * (1.0 / dx)       \
           * PMean(pp[ip - 1], pp[ip],                        \
                   permxp[ip - 1], permxp[ip])                \
           / viscosity;                                       \
         wp[im] = -coeff * diff                               \
           * RPMean(pp[ip - 1], pp[ip], prod_der, 0.0);       \
       })

#define Richards_SymmCorrection_Right                         \
  FACE(Right,                                                 \
       {                                                      \
         diff = pp[ip] - pp[ip + 1];                          \
         prod_der = rpdp[ip + 1] * dp[ip + 1]                 \
           + rpp[ip + 1] * ddp[ip + 1];                       \
         coeff = dt * z_mult_dat[ip] * ffx * (1.0 / dx)       \
           * PMean(pp[ip], pp[ip + 1],                        \
                   permxp[ip], permxp[ip + 1])                \
           / viscosity;                                       \
         ep[im] = coeff * diff                                \
           * RPMean(pp[ip], pp[ip + 1], 0.0, prod_der);       \
       })

#define Richards_SymmCorrection_Up                                 \
  FACE(Up,                                                         \
       {                                                           \
         diff = pp[ip - sy_v] - pp[ip];                            \
         prod_der = rpdp[ip - sy_v] * dp[ip - sy_v]                \
           + rpp[ip - sy_v] * ddp[ip - sy_v];                      \
         coeff = dt * z_mult_dat[ip] * ffy * (1.0 / dy)            \
           * PMean(pp[ip - sy_v], pp[ip],                          \
                   permxp[ip - sy_v], permxp[ip])                  \
           / viscosity;                                            \
         sop[im] = -coeff * diff                                   \
           * RPMean(pp[ip - sy_v], pp[ip], prod_der, 0.0);         \
       })

#define Richards_SymmCorrection_Down                               \
  FACE(Down,                                                       \
       {                                                           \
         diff = pp[ip] - pp[ip + sy_v];                            \
         prod_der = rpdp[ip + sy_v] * dp[ip + sy_v]                \
           + rpp[ip + sy_v] * ddp[ip + sy_v];                      \
         coeff = dt * z_mult_dat[ip] * ffy * (1.0 / dy)            \
           * PMean(pp[ip], pp[ip + sy_v],                          \
                   permxp[ip], permxp[ip + sy_v])                  \
           / viscosity;                                            \
         np[im] = -coeff * diff                                    \
           * RPMean(pp[ip], pp[ip + sy_v], 0.0, prod_der);         \
       })

#define Richards_SymmCorrection_Front                                   \
  FACE(Front,                                                           \
       {                                                                \
         mult_dat_mean = Mean(z_mult_dat[ip], z_mult_dat[ip - sz_v]);   \
         lower_cond = pp[ip - sz_v]                                     \
           - 0.5 * dz * mult_dat_mean                                   \
           * dp[ip - sz_v] * gravity;                                   \
         upper_cond = pp[ip] + 0.5 * dz                                 \
           * mult_dat_mean                                              \
           * dp[ip] * gravity;                                          \
         diff = lower_cond - upper_cond;                                \
         prod_der = rpdp[ip - sz_v] * dp[ip - sz_v]                     \
           + rpp[ip - sz_v] * ddp[ip - sz_v];                           \
         prod_lo = rpp[ip - sz_v] * dp[ip - sz_v];                      \
         coeff = dt * ffz * (1.0 / (dz * mult_dat_mean))                \
           * PMeanDZ(permzp[ip - sz_v], permzp[ip],                     \
                     z_mult_dat[ip - sz_v], z_mult_dat[ip])             \
           / viscosity;                                                 \
         lp[im] = -coeff *                                              \
           (diff * RPMean(lower_cond, upper_cond, prod_der, 0.0)        \
            - gravity * 0.5  * dz * mult_dat_mean * ddp[ip]             \
            * RPMean(lower_cond, upper_cond, prod_lo, prod));           \
       })

#define Richards_SymmCorrection_Back                                    \
  FACE(Back,                                                            \
       {                                                                \
         mult_dat_mean = Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v]);   \
         lower_cond = pp[ip]                                            \
           - 0.5 * dz * mult_dat_mean                                   \
           * dp[ip] * gravity;                                          \
         upper_cond = pp[ip + sz_v]                                     \
           + 0.5 * dz * mult_dat_mean                                   \
           * dp[ip + sz_v] * gravity;                                   \
         diff = lower_cond - upper_cond;                                \
         prod_der = rpdp[ip + sz_v] * dp[ip + sz_v]                     \
           + rpp[ip + sz_v] * ddp[ip + sz_v];                           \
         prod_lo = rpp[ip + sz_v] * dp[ip + sz_v];                      \
         coeff = dt * ffz * (1.0 / (dz * mult_dat_mean))                \
           * PMeanDZ(permzp[ip], permzp[ip + sz_v],                     \
                     z_mult_dat[ip], z_mult_dat[ip + sz_v])             \
           / viscosity;                                                 \
         up[im] = -coeff *                                              \
           (diff * RPMean(lower_cond, upper_cond, 0.0, prod_der)        \
            - gravity * 0.5  * dz * mult_dat_mean * ddp[ip]             \
            * RPMean(lower_cond, upper_cond, prod, prod_up));           \
       })

/*
#define Richards_SymmCorrection_Left                          \
  {                                                           \
    diff = pp[ip - 1] - pp[ip];                               \
    prod_der = rpdp[ip - 1] * dp[ip - 1]                      \
      + rpp[ip - 1] * ddp[ip - 1];                            \
    coeff = dt * z_mult_dat[ip] * ffx * (1.0 / dx)            \
      * PMean(pp[ip - 1], pp[ip],                             \
              permxp[ip - 1], permxp[ip])                     \
      / viscosity;                                            \
    wp[im] = -coeff * diff                                    \
      * RPMean(pp[ip - 1], pp[ip], prod_der, 0.0);            \
  }

#define Richards_SymmCorrection_Right                         \
  {                                                           \
    diff = pp[ip] - pp[ip + 1];                               \
    prod_der = rpdp[ip + 1] * dp[ip + 1]                      \
      + rpp[ip + 1] * ddp[ip + 1];                            \
    coeff = dt * z_mult_dat[ip] * ffx * (1.0 / dx)            \
      * PMean(pp[ip], pp[ip + 1],                             \
              permxp[ip], permxp[ip + 1])                     \
      / viscosity;                                            \
    ep[im] = coeff * diff                                     \
      * RPMean(pp[ip], pp[ip + 1], 0.0, prod_der);            \
  }

#define Richards_SymmCorrection_Up                                 \
  {                                                                \
    diff = pp[ip - sy_v] - pp[ip];                                 \
    prod_der = rpdp[ip - sy_v] * dp[ip - sy_v]                     \
      + rpp[ip - sy_v] * ddp[ip - sy_v];                           \
    coeff = dt * z_mult_dat[ip] * ffy * (1.0 / dy)                 \
      * PMean(pp[ip - sy_v], pp[ip],                               \
              permxp[ip - sy_v], permxp[ip])                       \
      / viscosity;                                                 \
    sop[im] = -coeff * diff                                        \
      * RPMean(pp[ip - sy_v], pp[ip], prod_der, 0.0);              \
  }

#define Richards_SymmCorrection_Down                             \
  {                                                              \
    diff = pp[ip] - pp[ip + sy_v];                               \
    prod_der = rpdp[ip + sy_v] * dp[ip + sy_v]                   \
      + rpp[ip + sy_v] * ddp[ip + sy_v];                         \
    coeff = dt * z_mult_dat[ip] * ffy * (1.0 / dy)               \
      * PMean(pp[ip], pp[ip + sy_v],                             \
              permxp[ip], permxp[ip + sy_v])                     \
      / viscosity;                                               \
    np[im] = -coeff * diff                                       \
      * RPMean(pp[ip], pp[ip + sy_v], 0.0, prod_der);            \
  }

#define Richards_SymmCorrection_Front                                   \
  {                                                                     \
    mult_dat_mean = Mean(z_mult_dat[ip], z_mult_dat[ip - sz_v]);        \
    lower_cond = pp[ip - sz_v]                                          \
      - 0.5 * dz * mult_dat_mean                                        \
      * dp[ip - sz_v] * gravity;                                        \
    upper_cond = pp[ip] + 0.5 * dz                                      \
      * mult_dat_mean                                                   \
      * dp[ip] * gravity;                                               \
    diff = lower_cond - upper_cond;                                     \
    prod_der = rpdp[ip - sz_v] * dp[ip - sz_v]                          \
      + rpp[ip - sz_v] * ddp[ip - sz_v];                                \
    prod_lo = rpp[ip - sz_v] * dp[ip - sz_v];                           \
    coeff = dt * ffz * (1.0 / (dz * mult_dat_mean))                     \
      * PMeanDZ(permzp[ip - sz_v], permzp[ip],                          \
                z_mult_dat[ip - sz_v], z_mult_dat[ip])                  \
      / viscosity;                                                      \
    lp[im] = -coeff *                                                   \
      (diff * RPMean(lower_cond, upper_cond, prod_der, 0.0)             \
       - gravity * 0.5  * dz * mult_dat_mean * ddp[ip]                  \
       * RPMean(lower_cond, upper_cond, prod_lo, prod));                \
  }

#define Richards_SymmCorrection_Back                                    \
  {                                                                     \
    mult_dat_mean = Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v]);        \
    lower_cond = pp[ip]                                                 \
      - 0.5 * dz * mult_dat_mean                                        \
      * dp[ip] * gravity;                                               \
    upper_cond = pp[ip + sz_v]                                          \
      + 0.5 * dz * mult_dat_mean                                        \
      * dp[ip + sz_v] * gravity;                                        \
    diff = lower_cond - upper_cond;                                     \
    prod_der = rpdp[ip + sz_v] * dp[ip + sz_v]                          \
      + rpp[ip + sz_v] * ddp[ip + sz_v];                                \
    prod_lo = rpp[ip + sz_v] * dp[ip + sz_v];                           \
    coeff = dt * ffz * (1.0 / (dz * mult_dat_mean))                     \
      * PMeanDZ(permzp[ip], permzp[ip + sz_v],                          \
                z_mult_dat[ip], z_mult_dat[ip + sz_v])                  \
      / viscosity;                                                      \
    up[im] = -coeff *                                                   \
      (diff * RPMean(lower_cond, upper_cond, 0.0, prod_der)             \
       - gravity * 0.5  * dz * mult_dat_mean * ddp[ip]                  \
       * RPMean(lower_cond, upper_cond, prod, prod_up));                \
  }
*/

#endif // _RICHARDS_SYMM_CORRECTION_H_
