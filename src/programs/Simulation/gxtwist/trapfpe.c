#ifdef TRAPFPE
#include <fpu_control.h>
#endif

void trapfpe_ ()
{

#ifdef TRAPFPE

    fpu_control_t cw = _FPU_MASK_PM   |      // bypass PrecisionLoss traps
                       _FPU_MASK_UM   |      // bypass Underflow traps
                 //    _FPU_MASK_OM   |      // bypass Overflow traps
                 //    _FPU_MASK_DM   |      // bypass Denormalized traps
                 //    _FPU_MASK_IM   |      // bypass Invalid traps
                 //    _FPU_MASK_ZM   |      // bypass ZeroDivide traps
                       _FPU_EXTENDED;        // enable extended precision

    _FPU_SETCW(cw);

#endif

}
