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
      BCStructPatchLoopX(i, j, k, fdir2, ival, bc_struct, ipatch, is,   \
      {                                              \
        ip = SubvectorEltIndex(p_sub, i, j, k);                         \
        value = bc_patch_values[ival];                                  \
        PatchFaceSwitch(fdir2,                                          \
                        {pp[ip - 1] = value;},                          \
                        {pp[ip + 1] = value;},                          \
                        {pp[ip - sy_v] = value;},                       \
                        {pp[ip + sy_v] = value;},                       \
                        {pp[ip - sz_v] = value;},                       \
                        {pp[ip + sz_v] = value;},                       \
                        {});                                            \
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

/**
 *  @brief Switch macro for iterating over faces, to be used in conjunction with BCStructPatchLoopX
 *
 *  @note Need to make this more meaningful and useful
 *
 *  @param[in] fdir Face direction to switch on
 *  @param[in] L Expression to execute on Left face
 *  @param[in] R Expression to execute on Right face
 *  @param[in] D Expression to execute on Down face
 *  @param[in] U Expression to execute on Upper face
 *  @param[in] B Expression to execute on Back face
 *  @param[in] F Expression to execute on Front face
 *  @param[in] DEFAULT Expression to execute for default case
 **/
#define PatchFaceSwitch(fdir, L, R, D, U, B, F, DEFAULT)  \
  {                                                       \
    switch (fdir)                                         \
    {                                                     \
      case GrGeomOctreeFaceL:                             \
      {                                                   \
        L;                                                \
        break;                                            \
      }                                                   \
      case GrGeomOctreeFaceR:                             \
      {                                                   \
        R;                                                \
        break;                                            \
      }                                                   \
      case GrGeomOctreeFaceD:                             \
      {                                                   \
        D;                                                \
        break;                                            \
      }                                                   \
      case GrGeomOctreeFaceU:                             \
      {                                                   \
        U;                                                \
        break;                                            \
      }                                                   \
      case GrGeomOctreeFaceB:                             \
      {                                                   \
        B;                                                \
        break;                                            \
      }                                                   \
      case GrGeomOctreeFaceF:                             \
      {                                                   \
        F;                                                \
        break;                                            \
      }                                                   \
      default:                                            \
      {                                                   \
        DEFAULT;                                          \
        break;                                            \
      }                                                   \
    }                                                     \
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

/**
 *  @brief Symmetric boundary condition corrections in richards jacobian in XY plane
 *
 *  @param[in] op Submatrix Stencil Data to write into
 *  @param[in] idx Base index for reading out of SubvectorData
 *  @param[in] jdx Index to write into op
 *  @param[in] ff ffx, ffy, or ffz to use
 *  @param[in] pos Integer offset in positive direction
 *  @param[in] neg Integer offset in negative direction
 *  @param[in] der Derivative to use (dx, dy, dz)
 *  @param[in] perm Permability data to use (permxp, permzp, etc)
 *  @param[in] meanA Value to use for when RPMean is true
 *  @param[in] meanB Value to use for when RPMean is false
 *  @param[in] coeff_sign Sign to use on coefficient value in final calculation
 **/
#define L830_XY(op, idx, jdx, ff, pos, neg, der, perm, meanA, meanB, coeff_sign)  \
  {                                                                     \
    diff = pp[idx + neg] - pp[idx + pos];                               \
    L830_prod_der(idx, pos + neg);                                      \
    L830_XY_coeff(idx, pos, neg, ff, der, perm);                         \
    L830_XY_calc(op, idx, jdx, pos, neg, meanA, meanB, coeff_sign);      \
}

/**
 *  @brief Symmetric boundary condition corrections in richards jacobian in Z plane
 *
 *  @param[in] op Submatrix Stencil Data to write into
 *  @param[in] idx Base index for reading out of SubvectorData
 *  @param[in] jdx Index to write into op
 *  @param[in] ff ffx, ffy, or ffz to use
 *  @param[in] pos Integer offset in positive direction
 *  @param[in] neg Integer offset in negative direction
 *  @param[in] prod_der_mean RPMean with prod_der value used in final calculation
 *  @param[in] prod_xtra_mean RPMean with prod_xtra value used in final calculation
 **/
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






  L830_Z(up, ip, im, ffz, sz_v, 0,
           RPMean(lower_cond, upper_cond, 0.0, prod_der),
           RPMean(lower_cond, upper_cond, prod, prod_xtra)),


#define L830_CALC_LEFT L830_XY(wp, ip, im, ffx, 0, -1, dx, permxp, prod_der, 0.0, 1.0)
#define L830_CALC_RIGHT L830_XY(ep, ip, im, ffx, 1,  0, dx, permxp, 0.0, prod_der, 1.0)
#define L830_CALC_UP L830_XY(sop, ip, im, ffy, 0, -sy_v, dy, permyp, prod_der, 0.0, -1.0)
#define L830_CALC_DOWN L830_XY(np, ip, im, ffy, sy_v, 0, dy, permyp, 0.0, prod_der, -1.0)
#define L830_CALC_FRONT L830_Z(lp, ip, im, ffz, 0, -sz_v, \
                               RPMean(lower_cond, upper_cond, prod_der, 0.0), \
                               RPMean(lower_cond, upper_cond, prod_xtra, prod))
#define L830_CALC_BACK L830_Z(up, ip, im, ffz, sz_v, 0, \
                              RPMean(lower_cond, upper_cond, 0.0, prod_der), \
                              RPMean(lower_cond, upper_cond, prod, prod_xtra))

#define L830_CALC(face) L830_CALC_##face


/**
 *  @brief Calculation used in Dirichlet Boundary Condition case
 *
 *  @param[in] zp Submatrix Stencil Data to assign op to
 *  @param[in] idx Base index for reading out of SubvectorData
 *  @param[in] jdx Index to write into op
 *  @param[in] ff ffx, ffy, or ffz to use
 *  @param[in] der Derivative to use (dx, dy, dz)
 *  @param[in] perm Permability data to use (permxp, permzp, etc)
 *  @param[in] offset Offset to use when reading out of rpp data
 *  @param[in] lval Value to use on the left side of the rhs of diff calculation
 *  @param[in] rval Value to use on the right side of the rhs of diff calculation
 *  @param[in] coeff_sign Sign to use on coefficient value in final calculation
 *  @param[in] A Value to use for when first RPMean is true
 *  @param[in] B Value to use for when first RPMean is false
 *  @param[in] C Value to use for when second RPMean is true
 *  @param[in] D Value to use for when second RPMean is false
 **/
#define L970_Dirichlet_calc(zp, idx, ff, der, perm, offset, \
                       lval, rval, coeff_sign, A, B, C, D)  \
  {                                                    \
    coeff = dt * ff * z_mult_dat[idx] * (2.0 / der)    \
      * perm[idx] / viscosity;                         \
    op = zp;                                           \
    prod_val = rpp[idx + offset] * den_d;              \
    diff = lval + pp[idx] + rval;                      \
    o_temp = (coeff_sign * coeff)                      \
      * (diff * RPMean(value, pp[idx], A, B)           \
         - RPMean(value, pp[idx], C, D));              \
  }

#define L970_Dirichlet_Left \
  L970_Dirichlet_calc(wp, ip, ffx, dx, permxp, -1, -value,  \
                      0, 1, 0.0, prod_der, prod, prod_val)
#define L970_Dirichlet_Right \
  L970_Dirichlet_calc(ep, ip, ffx, dx, permxp, 1, 0, -value,  \
                      -1, prod_der, 0.0, prod, prod_val)
#define L970_Dirichlet_Up \
  L970_Dirichlet_calc(sop, ip, ffy, dy, permyp, -sy_v, -value, 0, \
                      1, 0.0, prod_der, prod_val, prod)
#define L970_Dirichlet_Down                                             \
  L970_Dirichlet_calc (np, ip, ffy, dy, permyp, sy_v, 0, -value,        \
                       -1, prod_der, 0.0, prod, prod_val)

#define L970_Dirichlet_Front \
  {                          \
  coeff = dt * ffz *                                          \
    (2.0 / (dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])))  \
    * permzp[ip] / viscosity;                                   \
  op = lp;                                                      \
  prod_val = rpp[ip - sz_v] * den_d;                                \
  lower_cond = (value) - 0.5 * dz * z_mult_dat[ip] * den_d * gravity; \
  upper_cond = (pp[ip]) + 0.5 * dz * z_mult_dat[ip] * dp[ip] * gravity; \
  diff = lower_cond - upper_cond;                                       \
  o_temp = coeff                                                        \
    * (diff * RPMean(lower_cond, upper_cond, 0.0, prod_der)             \
       + ((-1.0 - gravity * 0.5 * dz                                    \
           * Mean(z_mult_dat[ip], z_mult_dat[ip - sz_v]) * ddp[ip])     \
          * RPMean(lower_cond, upper_cond, prod_val, prod)));           \
  }

#define L970_Dirichlet_Back \
  {                                             \
  coeff = dt * ffz *                                                    \
    (2.0 / (dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])))          \
    * permzp[ip] / viscosity;                                           \
  op = up;                                                              \
  prod_val = rpp[ip + sz_v] * den_d;                                    \
  lower_cond = (pp[ip]) - 0.5 * dz *                                    \
    Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v]) * dp[ip] * gravity;     \
  upper_cond = (value) + 0.5 * dz *                                     \
    Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v]) * den_d * gravity;      \
  diff = lower_cond - upper_cond;                                       \
  o_temp = -coeff                                                       \
    * (diff * RPMean(lower_cond, upper_cond, prod_der, 0.0)             \
       + ((1.0 - gravity * 0.5 * dz                                     \
           * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v]) * ddp[ip])     \
          * RPMean(lower_cond, upper_cond, prod, prod_val)));           \
  }

#define L970_Dirichlet(face) L970_Dirichlet_##face


/*
  Found while reading about x-macros:
  By merely writing this macro, I am a horrible person who deserves terrible things.
*/

#define WITH_BODY(_case, body)                  \
  case _case:                                   \
  {                                             \
    body;                                       \
    break;                                      \
  }

#define ONLY_CASE(_case) \
  case _case:

#define CASE_BODY_SELECTION(arg1, arg2, arg3, ...) arg3
#define XCASE(...) CASE_BODY_SELECTION(__VA_ARGS__, WITH_BODY, ONLY_CASE)(__VA_ARGS__)

#define EXPAND_CASES_1(a, ...) a
#define EXPAND_CASES_2(a, ...) a EXPAND_CASES_1(__VA_ARGS__)
#define EXPAND_CASES_3(a, ...) a EXPAND_CASES_2(__VA_ARGS__)
#define EXPAND_CASES_4(a, ...) a EXPAND_CASES_3(__VA_ARGS__)
#define EXPAND_CASES_5(a, ...) a EXPAND_CASES_4(__VA_ARGS__)
#define EXPAND_CASES_6(a, ...) a EXPAND_CASES_5(__VA_ARGS__)

#define XSWITCH(key, numcase, ...)														\
  switch (key) {																							\
		EXPAND_CASES_ ##numcase(__VA_ARGS__)											\
    default:                                                  \
    break;                                                    \
  }


#endif // _BC_BRANCHING_H
