#ifndef _BC_BRANCHING_H
#define _BC_BRANCHING_H

#include "overland_richards.h"
#include "dirichlet_richards.h"
#include "flux_richards.h"
#include "richards_symm_correction.h"
#include "problem_bc.h"

/*****************************************************************************
 * Header file for defining equation macros for usage in Boundary Condition Loops
 *
 *****************************************************************************/

// NOTE: This macro will cause breakage on the switch if too many arguments are passed in!
// Solution is to just increase the total number of args it can take
#define EVAL_COUNT_VARARGS(_1, _2, _3, _4, _5, _6, _7, _8, _9, _10, N, ...) N
#define COUNT_VARARGS(...) EVAL_COUNT_VARARGS("n", ##__VA_ARGS__, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0)

// This is silly but the preprocessor requires indirection to evaluate this correctly
#define EXPAND_INDIRECTION(n) EXPAND_CASES_ ## n
#define EXPAND_CASES(n) EXPAND_INDIRECTION(n)

/* TODO: Replace PatchFaceSwitch
#define DirichletBC_Pressure(bc_struct, ipatch, is, p_sub, pp, sy_v, sz_v) \
  bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is);         \
                                                                        \
  switch (BCStructBCType(bc_struct, ipatch))                            \
  {                                                                     \
    case DirichletBC:                                                   \
    {                                                                   \
      BCStructPatchLoopXX(i, j, k, ival, bc_struct, ipatch, is,         \
      {                                             \
        ip = SubvectorEltIndex(p_sub, i, j, k);                         \
        value = bc_patch_values[ival];                                  \
      },                                                                \
                          {},                                           \
                          PatchFaceSwitch(fdir2,                      \
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
*/

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

#define L970_Dirichlet_Left                                 \
  L970_Dirichlet_calc(wp, ip, ffx, dx, permxp, -1, -value,  \
                      0, 1, 0.0, prod_der, prod, prod_val)
#define L970_Dirichlet_Right                                  \
  L970_Dirichlet_calc(ep, ip, ffx, dx, permxp, 1, 0, -value,  \
                      -1, prod_der, 0.0, prod, prod_val)
#define L970_Dirichlet_Up                                         \
  L970_Dirichlet_calc(sop, ip, ffy, dy, permyp, -sy_v, -value, 0, \
                      1, 0.0, prod_der, prod_val, prod)
#define L970_Dirichlet_Down                                             \
  L970_Dirichlet_calc (np, ip, ffy, dy, permyp, sy_v, 0, -value,        \
                       -1, prod_der, 0.0, prod, prod_val)

#define L970_Dirichlet_Front                                            \
  {                                                                     \
    coeff = dt * ffz *                                                  \
      (2.0 / (dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])))        \
      * permzp[ip] / viscosity;                                         \
    op = lp;                                                            \
    prod_val = rpp[ip - sz_v] * den_d;                                  \
    lower_cond = (value) - 0.5 * dz * z_mult_dat[ip] * den_d * gravity; \
    upper_cond = (pp[ip]) + 0.5 * dz * z_mult_dat[ip] * dp[ip] * gravity; \
    diff = lower_cond - upper_cond;                                     \
    o_temp = coeff                                                      \
      * (diff * RPMean(lower_cond, upper_cond, 0.0, prod_der)           \
         + ((-1.0 - gravity * 0.5 * dz                                  \
             * Mean(z_mult_dat[ip], z_mult_dat[ip - sz_v]) * ddp[ip])   \
            * RPMean(lower_cond, upper_cond, prod_val, prod)));         \
  }

#define L970_Dirichlet_Back                                             \
  {                                                                     \
    coeff = dt * ffz *                                                  \
      (2.0 / (dz * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v])))        \
      * permzp[ip] / viscosity;                                         \
    op = up;                                                            \
    prod_val = rpp[ip + sz_v] * den_d;                                  \
    lower_cond = (pp[ip]) - 0.5 * dz *                                  \
      Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v]) * dp[ip] * gravity;   \
    upper_cond = (value) + 0.5 * dz *                                   \
      Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v]) * den_d * gravity;    \
    diff = lower_cond - upper_cond;                                     \
    o_temp = -coeff                                                     \
      * (diff * RPMean(lower_cond, upper_cond, prod_der, 0.0)           \
         + ((1.0 - gravity * 0.5 * dz                                   \
             * Mean(z_mult_dat[ip], z_mult_dat[ip + sz_v]) * ddp[ip])   \
            * RPMean(lower_cond, upper_cond, prod, prod_val)));         \
  }

#define L970_Dirichlet(face) L970_Dirichlet_##face


/*
  Found while reading about x-macros:
  By merely writing this macro, I am a horrible person who deserves terrible things.
*/


#define Left 0
#define Right 1
#define Up 2
#define Down 3
#define Front 4
#define Back 5
#define FACE(a, b) XCASE(a, b)

#define WITH_BODY(_case, body)                  \
  case _case:                                   \
  {                                             \
    body;                                       \
    break;                                      \
  }

#define ONLY_CASE(_case) \
  case _case:

#define CASE_BODY_SELECTION(arg1, arg2, arg3, ...) arg3
#define XCASE(...) \
  CASE_BODY_SELECTION(__VA_ARGS__, WITH_BODY, ONLY_CASE)(__VA_ARGS__)

#define EXPAND_CASES_1(a, ...) a
#define EXPAND_CASES_2(a, ...) a EXPAND_CASES_1(__VA_ARGS__)
#define EXPAND_CASES_3(a, ...) a EXPAND_CASES_2(__VA_ARGS__)
#define EXPAND_CASES_4(a, ...) a EXPAND_CASES_3(__VA_ARGS__)
#define EXPAND_CASES_5(a, ...) a EXPAND_CASES_4(__VA_ARGS__)
#define EXPAND_CASES_6(a, ...) a EXPAND_CASES_5(__VA_ARGS__)

#define XSWITCH(key, ...)                                     \
  switch (key) {																							\
		EXPAND_CASES(COUNT_VARARGS(__VA_ARGS__))(__VA_ARGS__)     \
  }


// More silly indirection to expand the actual equation bodies
#define DoRichards_SymmCorrection(J_sub, p_sub, rpp, dp, i, j, k)       \
  DoRichards_SymmCorrection_(Richards_SymmCorrection(p_sub, J_sub, rpp, dp, i, j, k))

#define DoRichards_SymmCorrection_(equations)  DoRichards_SymmCorrection__(equations)

#define DoRichards_SymmCorrection__(prologue, epilogue, ...)    \
  ForBCStructNumPatches(ipatch, bc_struct)                     \
  {                                                            \
    BCStructPatchLoopXX(i, j, k, ival, bc_struct, ipatch, is,  \
                        prologue, epilogue, __VA_ARGS__);      \
  }

#define DoRichards_BC_Contrib(equations)                            \
  ForBCStructNumPatches(ipatch, bc_struct)                          \
  {                                                                 \
    bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is);   \
    switch (BCStructBCType(bc_struct, ipatch))                      \
    {                                                               \
      equations;                                                    \
    }                                                               \
  }

#define EXPAND_EQUATIONS(prologue, epilogue, ...)           \
  BCStructPatchLoopXX(i, j, k, ival, bc_struct, ipatch, is, \
                      prologue, epilogue, __VA_ARGS__);     \

#define EXPAND_EQUATIONS_(equations) EXPAND_EQUATIONS(equations)

#define ApplyPatch(_case, equations)                          \
  case _case:                                                 \
  {                                                           \
    EXPAND_EQUATIONS(equations);                              \
    break;                                                    \
  }

#define ApplyPatchSubtypes(_case, _subcase, equations) \
  case _case:                                    \
  {                                              \
    switch(_subcase)                             \
    {                                            \
      equations;                               \
    }                                            \
  }

#define PatchSubtype(_case, equations) \
  case _case:                          \
  {                                    \
    EXPAND_EQUATIONS(equations);       \
    break;                             \
  }

#define DoRichardsDirichlet_Contrib(equations)  \
  case DirichletBC:                             \
  {                                             \
    DoRichardsDirichlet_Contrib_(equations);    \
    break;                                      \
  }
#define DoRichardsDirichlet_Contrib_(prologue, epilogue, ...) \
  BCStructPatchLoopXX(i, j, k, ival, bc_struct, ipatch, is,   \
                      prologue, epilogue, __VA_ARGS__);       \

#define DoRichardsFlux_Contrib(equations)  \
  case FluxBC:                             \
  {                                        \
    DoRichardsFlux_Contrib_(equations);    \
    break;                                 \
  }

#define DoRichardsFlux_Contrib_(prologue, epilogue, ...) \
  BCStructPatchLoopXX(i, j, k, ival, bc_struct, ipatch, is, \
                      prologue, epilogue, __VA_ARGS__);

#endif // _BC_BRANCHING_H
