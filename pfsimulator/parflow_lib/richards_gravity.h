#ifndef _RICHARDS_GRAVITY_H
#define _RICHARDS_GRAVITY_H

/*********************************************
 *
 * Header file for RichardsJacobianEval gravity
 * and 2nd order derivative contributions
 *
 *********************************************/

#define Do_RichardsGravityX(prologue, flow_positive, flow_negative) \
  RichardsGravity(prologue, flow_positive, flow_negative);
#define Do_RichardsGravityY(prologue, flow_positive, flow_negative) \
  RichardsGravity(prologue, flow_positive, flow_negative);
#define Do_RichardsGravityZ(prologue, flow_positive, flow_negative) \
  RichardsGravity(prologue, flow_positive, flow_negative);

#define PROLOGUE(x) x
#define FLOW(dir, body) body

// a - b >= 0 ? c : d
#define RichardsGravity(prologue, flow_positive, flow_negative) \
  {                                                             \
    prologue;                                                   \
    if (updir >= 0)                                             \
    {                                                           \
      flow_positive;                                            \
    } else {                                                    \
      flow_negiatve;                                            \
    }                                                           \
  }


#endif // _RICHARDS_GRAVITY_H
