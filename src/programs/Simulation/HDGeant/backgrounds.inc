c
c  These parameters describe the background particle flux that is
c  superimposed on top of the standard generation.
c
c  author: Richard Jones
c  date: March 22, 2005
c
c  Notes:
c  -------
c  bgrate: the rate (1/ns) of beam photons above the threshold BCUTE
c          for discrete bremsstrahlung to generate during the ADC gate;
c          default value of 0 represents a background-free simulation.
c  bggate: time (ns) of the start (1) and end (2) of the gate during
c          which random background is simulated, relative to the time
c          of the photon that originated the event.
c  bgtagonly: flag saying whether to tag and track all of the background
c          beam photons generated within bggate (bgtagonly=0), or to
c          only record their tags but not track them (bgtagonly=1).

      real bgrate
      real bggate(2)
      integer bgtagonly
      common /backgrounds/bgrate,bggate,bgtagonly
