#ifndef _BC_RICHARDS_H
#define _BC_RICHARDS_H

/**
 * @brief Indirection Macro to expand Patch equations
 *
 * Expects that the provided equations should be applied over the entire
 * patch, hence the BCStructPatchLoop call.
 *
 * @param[in] prologue To be done before every face equation
 * @param[in] epilogue To be done after every face equation
 * @param[in] ... Variadic set of possible face equations
 **/
#define EXPAND_PATCH_EQUATIONS(_list, prologue, epilogue, ...)           \
  BCStructPatchLoopXX(i, j, k, ival, bc_struct, ipatch, is, prologue, epilogue, __VA_ARGS__);

#define EXPAND_PATCH_EQUATIONS2(_list, prologue, epilogue, ...) \
  BCStructPatchLoop_Collected(ipatch, _list, i, j, k,           \
                              prologue, epilogue, __VA_ARGS__);

// Generates case statement for patch type
#define ApplyPatch(_case, ...)                                      \
  case _case:                                                       \
  {                                                                 \
    EXPAND_PATCH_EQUATIONS2(__VA_ARGS__);                            \
    break;                                                          \
  }                                                                 \

// For nested Overland BC types
// @TODO: This is being replaced by giving each subtype an actual enum
#define ApplyPatchSubtypes(_case, _subcase, equations)  \
  case _case:                                           \
  {                                                     \
    switch ((_subcase))                                 \
    {                                                   \
      equations;                                        \
    }                                                   \
  }

// For clarity, no different than the ApplyPatch macro
#define PatchSubtype(_case, equations) ApplyPatch(_case, equations)

// Add DirichletBC pressures
#define Do_DirichletBCPressureContrib(patch_idx, bc_struct, list, subgrid_idx, \
                                      pressure_subvector, pressure_data, \
                                      y_offset, z_offset)               \
  {                                                                     \
    ForBCStructNumPatches(patch_idx, bc_struct)                         \
    {                                                                   \
      bc_patch_values = BCStructPatchValues(bc_struct, patch_idx, subgrid_idx); \
      switch (BCStructBCType(bc_struct, patch_idx))                     \
      {                                                                 \
        case DirichletBC:                                               \
        {                                                               \
          BCStructPatchLoop_Collected(patch_idx, list, i, j, k, \
          {                                                             \
            ip = SubvectorEltIndex(pressure_subvector, i, j, k);        \
            value = bc_patch_values[ival];                              \
            pressure_data[ip + fdir[0] * 1                              \
            + fdir[1] * y_offset                                        \
            + fdir[2] * z_offset] = value;                              \
          },{});                                                        \
        }                                                               \
      }                                                                 \
    }                                                                   \
  }


#define Do_RichardsGravityX(prologue, flow_pos, flow_neg) \
  RichardsGravity(prologue, flow_pos, flow_neg);
#define Do_RichardsGravityY(prologue, flow_pos, flow_neg) \
  RichardsGravity(prologue, flow_pos, flow_neg);
#define Do_RichardsGravityZ(prologue, flow_pos, flow_neg) \
  RichardsGravity(prologue, flow_pos, flow_neg);

// Add all boundary condition contributions in Richards Jacobian
#define Do_RichardsBCContrib(equations)                           \
  ForBCStructNumPatches(ipatch, bc_struct)                        \
  {                                                               \
    bc_patch_values = BCStructPatchValues(bc_struct, ipatch, is); \
    switch (BCStructBCType(bc_struct, ipatch))                    \
    {                                                             \
      equations;                                                  \
    }                                                             \
  }

#define Do_RichardsBuildJC(equations)           \
  ForBCStructNumPatches(ipatch, bc_struct)      \
  {                                             \
    switch (BCStructBCType(bc_struct, ipatch))  \
    {                                           \
      equations;                                \
    }                                           \
  }

#define FLOW(dir, body) body

// a - b >= 0 ? c : d
#define RichardsGravity(prologue, flow_pos, flow_neg)       \
  {                                                         \
    prologue;                                               \
    if (updir >= 0)                                         \
    {                                                       \
      flow_pos;                                             \
    } else {                                                \
      flow_neg;                                             \
    }                                                       \
  }


// Add symmetric contributions for Richards Jacobian
#define Do_RichardsSymmCorrection(prologue, epilogue, ...) \
  ForBCStructNumPatches(ipatch, bc_struct)                 \
  {                                                        \
    BCStructPatchLoopXX(i, j, k, ival, bc_struct, ipatch, is, prologue, epilogue, __VA_ARGS__); \
  }


//#define Do_DerivativeAndGravityContrib()

#endif // _BC_RICHARDS_H
