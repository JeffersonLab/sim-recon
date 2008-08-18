*
* Goptimize - do any Geant3 geometry/tracking optimizations
*
* The actual definitions of the materials, tracking media and volume
* tree are found in the file hddsGeant3.f which is generated automatically
* from the HDDS xml geometry database by the translator hdds-geant.
*
* NOTE: It is tempting to put Geant geometry and tracking optimization
*	commands into hddsGeant3.f at the point where the geometry is
*	being defined.  DO NOT DO THAT.  Put them here in this file.
*
* This routine is part of the HDGeant simulation package
*
* Author: Richard Jones
*	  University of Connecticut
*	  July 5, 2001
*------------------------------------
 
      subroutine Goptimize
      
c     User optimizations go here
c     such as: gsnext, gunear, gsord, 
c              tracking medium parameter modifications,
c              graphical attributes of volumes,
c              etc...

      end
