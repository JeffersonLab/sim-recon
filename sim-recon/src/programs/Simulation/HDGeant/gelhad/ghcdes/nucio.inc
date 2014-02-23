**** Sequence /nucio from PAM geanh321 ****
      COMMON/NUCIN /TECM,AMASS(18),NPG,KGENEV
      COMMON/NUCOUT/PCM(5,18),WGT
#if  ! _SINGLE_
#include "nuciod.inc"
#endif
C
