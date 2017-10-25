#include <stdio.h>
#include <iostream>

#include "wave.h"

wave::wave(const wave& o)
{
  l = o.l;
  m = o.m;
  name = o.name;
  idx = o.idx;
  phaseLocked = o.phaseLocked;
}


waveset::waveset() { return; }

