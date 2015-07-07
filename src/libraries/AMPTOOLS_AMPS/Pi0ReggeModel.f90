subroutine DiffCross(Eg,t,DSG)
! return the differential cross section
! Input:
!   Eg      beam energy in the lab frame
!   t       the madelstam variable
! Output:
!   DSG     differential cross section
    double precision PI, xmg, xmpi, xmp, s, som, sint2
    double precision g1p, g4p, g1c, g4c, g2p, g2c, alp1, alp2, alp3
    double precision, dimension(4)  :: mass
    double complex tc,kt,pt,cost, f01p,f11p,f01m,f11m, res1, Fvec
    double complex   :: II = (0.d0,1.0d0)
    double complex F1,F2,F3,F4, xi1,xi2,xi3
    double precision, dimension(6)  :: paralp, parcou

    ! constant parameters
    PI=3.141592653589793238462643383279502884197
    xmg   = 0.0     ! photon mass
    xmpi = 0.138    ! pion mass
    xmp  = 0.938    ! proton mass
    mass = (/ xmg,xmp,xmpi,xmp /)

    ! kinematic quantities
    s = xmp**2 + 2*xmp*Eg
    som = 2*xmp**2 + xmpi**2
    tc = complex(t,0)
    kt = (t - xmpi**2)/(2*sqrt(tc))
    pt = sqrt(tc-4*xmp**2)/2.0
    cost = (2*s + t - som)/(4*kt*pt)
    sint2 = 1- abs(cost)**2

    ! insert Ampl_F function in single subroutine for AmpTools model

    ! parameters: trajectories and couplings
    ! the (linear) trajectories are
    !   vector pole : alp1 = paralp(1) + t*paralp(2)
    !   vector cut  : alp2 = paralp(3) + t*paralp(4)
    !   axial  pole : alp3 = paralp(5) + t*paralp(6)
    ! the couplings are (p=pole and c=cut)
    !   parcoup = (/ g1p, g4p, g1c, g4c, g2p, g2c /)
    paralp = (/ 0.442457d0, 1.09991d0, 0.46181d0, 0.16654d0, -0.193332d0, 1.0213d0  /)
    parcou = (/ 3.8873d0, -10.1403d0, -1.76187d0, -3.58089d0, -8.887d0, 0.0d0/)

    g1p = parcou(1)
    g4p = parcou(2)
    g1c = parcou(3)
    g4c = parcou(4)
    g2p = parcou(5)
    g2c = parcou(6)

    alp1 = paralp(1) + paralp(2)*t      ! vector pole trajectory
    alp2 = paralp(3) + paralp(4)*t      ! vector cut trajectory
    alp3 = paralp(5) + paralp(6)*t      ! axial pole trajectory

    xi1 = gamma(1 - alp1)*( 1 - exp(-II*PI*alp1) )/2.d0*s**(alp1-1)
    xi2 = gamma(1 - alp2)*( 1 - exp(-II*PI*alp2) )/2.d0*s**(alp2-1)/log(s)
    xi3 = gamma(1 - alp3)*( 1 - exp(-II*PI*alp3) )/2.d0*s**(alp3-1)

    F1 = -(t*g1p - 2*xmp*g4p)*xi1 - (t*g1c - 2*xmp*g4c)*xi2
    F2 = t*g2p*xi3
    F3 = t*(2*xmp*g1p - g4p)*xi1 + t*(2*xmp*g1c - g4c)*xi2
    F4 = 0

    ! use F1-F4 from old Ampl_F subroutine
    f01p = kt*sqrt(tc)/2.0*F1
    f11p = kt/2.0*F3;
    f01m = pt*kt*F2;
    f11m = -pt*kt*sqrt(tc)*F4

    res1 = 2*(abs(f01p)**2  + abs(f01m)**2)*abs(sint2)
    res1 = res1  + abs(f11p - f11m)**2*(1 - cost)**2 + abs(f11p + f11m)**2*(1 + cost)**2
    DSG = real(res1)/(64*PI*xmp**2*Eg**2)*389.4

end subroutine
