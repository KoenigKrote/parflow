#ifndef _BC_BRANCHING_H
#define _BC_BRANCHING_H

#include "problem_bc.h"

/*****************************************************************************
 * Header file for defining equation macros for usage in Boundary Condition Loops
 *
 *****************************************************************************/

#define DirichletBC_Pressure(bc_struct, ipatch, is, p_sub, pp, sy_v, sz_v) \
  bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is);         \
                                                                        \
  switch (BCStructBCType(bc_struct, ipatch))                            \
  {                                                                     \
    case DirichletBC:                                                   \
    {                                                                   \
      BCStructPatchLoop(i, j, k, fdir, ival, bc_struct, ipatch, is,     \
      {                                               \
        ip = SubvectorEltIndex(p_sub, i, j, k);                         \
        value = bc_patch_values[ival];                                  \
        pp[ip + fdir[0] * 1 + fdir[1] * sy_v + fdir[2] * sz_v] = value; \
      });                                                               \
      break;                                                            \
    }                                                                   \
  }                                                                     \


#define TemperatureContrib(cp, wp, ep, sop, np, lp, up, im, sy_m, sz_m, \
                           west_temp, east_temp, north_temp, south_temp, \
                           upper_temp, lower_temp, symm_part,           \
                           sym_east_temp, sym_north_temp, sym_upper_temp) \
  cp[im] -= west_temp + south_temp + lower_temp;                        \
  cp[im + 1] -= east_temp;                                              \
  cp[im + sy_m] -= north_temp;                                          \
  cp[im + sz_m] -= upper_temp;                                          \
  if (!symm_part)                                                       \
  {                                                                     \
    ep[im] += east_temp;                                                \
    np[im] += north_temp;                                               \
    up[im] += upper_temp;                                               \
    wp[im + 1] += west_temp;                                            \
    sop[im + sy_m] += south_temp;                                       \
    lp[im + sz_m] += lower_temp;                                        \
  }                                                                     \
  else                                                                  \
  {                                                                     \
    ep[im] += sym_east_temp;                                            \
    np[im] += sym_north_temp;                                           \
    up[im] += sym_upper_temp;                                           \
  }

#define L830_lower(idx, pos, neg)                                       \
  {                                                                     \
    lower_cond = pp[idx + neg] - 0.5 * dz *                             \
      Mean(z_mult_dat[idx], z_mult_dat[idx + pos + neg])                \
      * dp[idx + neg] * gravity;                                        \
  }

#define L830_upper(idx, pos, neg)                         \
  {                                                       \
    upper_cond = pp[idx + pos] + 0.5 * dz *               \
      Mean(z_mult_dat[idx], z_mult_dat[idx + pos + neg])  \
      * dp[idx + pos] * gravity;                          \
  }

#define L830_prod_der(idx, offset)                    \
  {                                                   \
    prod_der = rpdp[idx + offset] * dp[idx + offset]  \
      + rpp[ip + offset] * ddp[ip + offset];          \
  }

#define L830_prod_xtra(idx, offset)                   \
  {                                                   \
    prod_xtra = rpp[idx + offset] * dp[idx + offset]; \
  }

#define L830_XY_coeff(idx, pos, neg, ff, der, perm)          \
  {                                                         \
  coeff = dt * z_mult_dat[idx] * ff * (1.0 / der)           \
    * PMean(pp[idx + neg], pp[idx + pos],                   \
            perm[idx + neg], perm[idx + pos]) / viscosity;  \
  }

#define L830_XY_calc(op, idx, jdx, pos, neg, meanA, meanB, coeff_sign)  \
  {                                                                    \
    op[jdx] = (coeff_sign * coeff) * diff                              \
      * RPMean(pp[idx + neg], pp[idx + pos], meanA, meanB);            \
  }

#define L830_Z_coeff(idx, ff, pos, neg)                                 \
  {                                                                     \
    coeff = dt * ff * (1.0 / (dz * Mean(z_mult_dat[idx],                \
                                        z_mult_dat[idx + pos + neg])))  \
      * PMeanDZ(permzp[idx + neg], permzp[idx + pos],                   \
                z_mult_dat[idx + neg], z_mult_dat[idx + pos])           \
      / viscosity;                                                      \
  }

#define L830_Z_calc(op, idx, offset, prod_der_mean, prod_xtra_mean)     \
  {                                                                     \
    op[idx] = -coeff * (diff * prod_der_mean - gravity * 0.5 * dz *     \
                        Mean(z_mult_dat[idx], z_mult_dat[idx + offset]) * \
                        ddp[idx] * prod_xtra_mean);                     \
  }

#define L830_XY(op, idx, jdx, ff, pos, neg, der, perm, meanA, meanB, coeff_sign)  \
  {                                                                     \
    diff = pp[idx + neg] - pp[idx + pos];                               \
    L830_prod_der(idx, pos + neg);                                      \
    L830_XY_coeff(idx, pos, neg, ff, der, perm);                         \
    L830_XY_calc(op, idx, jdx, pos, neg, meanA, meanB, coeff_sign);      \
}

#define L830_Z(op, idx, jdx, ff, pos, neg, prod_der_mean, prod_xtra_mean) \
  {                                                                     \
    L830_lower(idx, pos, neg);                                          \
    L830_upper(idx, pos, neg);                                          \
    diff = lower_cond - upper_cond;                                     \
    L830_prod_der(idx, pos + neg);                                      \
    L830_prod_xtra(idx, pos + neg);                                     \
    L830_Z_coeff(idx, ff, pos, neg);                                    \
    L830_Z_calc(op, jdx, pos + neg, prod_der_mean, prod_xtra_mean);     \
  }


#endif // _BC_BRANCHING_H
