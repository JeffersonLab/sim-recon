      program mcfast
      call makeGeometry
      end

      subroutine detectorDef
      character*40 name
      character*40 geom_id
      common /detector/name,geom_id
      name = "AMEZON"
      geom_id = "CENTRAL"
      call makedetector
      end

      subroutine MaterialDef
      character*40 name
      real a
      real z
      real density
      real radlen
      real abslen
      real collen
      real dedx
      common /Material/name,a,z,density,radlen,abslen,collen,dedx
      name = "VACU"
      a = 0.1000E-15
      z = 0.1000E-15
      density = 0.1000E-15
      radlen = 0.1000E+17
      abslen = 0.1000E+17
      collen = 0.1000E+17
      dedx = 0.1000E-15
      call makeMaterial
      name = "AIR"
      a = 0.1461E+02
      z = 0.7300E+01
      density = 0.1205E-02
      radlen = 0.3042E+05
      abslen = 0.7469E+05
      collen = 0.4960E+05
      dedx = 0.2269E-05
      call makeMaterial
      name = "ALUM"
      a = 0.2698E+02
      z = 0.1300E+02
      density = 0.2700E+01
      radlen = 0.8900E+01
      abslen = 0.3941E+02
      collen = 0.2615E+02
      dedx = 0.4360E-02
      call makeMaterial
      name = "IRON"
      a = 0.5585E+02
      z = 0.2600E+02
      density = 0.7870E+01
      radlen = 0.1760E+01
      abslen = 0.1676E+02
      collen = 0.1052E+02
      dedx = 0.1142E-01
      call makeMaterial
      name = "SILI"
      a = 0.2809E+02
      z = 0.1400E+02
      density = 0.2330E+01
      radlen = 0.9360E+01
      abslen = 0.4549E+02
      collen = 0.3030E+02
      dedx = 0.3877E-02
      call makeMaterial
      name = "SCIN"
      a = 0.1200E+02
      z = 0.6000E+01
      density = 0.1032E+01
      radlen = 0.4240E+02
      abslen = 0.7646E+02
      collen = 0.5659E+02
      dedx = 0.2016E-02
      call makeMaterial
      name = "CU"
      a = 0.6355E+02
      z = 0.2900E+02
      density = 0.8960E+01
      radlen = 0.1430E+01
      abslen = 0.1506E+02
      collen = 0.9554E+01
      dedx = 0.1257E-01
      call makeMaterial
      name = "TITA"
      a = 0.4788E+02
      z = 0.2200E+02
      density = 0.4540E+01
      radlen = 0.3560E+01
      abslen = 0.2751E+02
      collen = 0.1760E+02
      dedx = 0.6701E-02
      call makeMaterial
      name = "BE"
      a = 0.9012E+01
      z = 0.4000E+01
      density = 0.1848E+01
      radlen = 0.3530E+02
      abslen = 0.4069E+02
      collen = 0.3019E+02
      dedx = 0.2946E-02
      call makeMaterial
      name = "CARB"
      a = 0.1201E+02
      z = 0.6000E+01
      density = 0.2265E+01
      radlen = 0.1880E+02
      abslen = 0.3810E+02
      collen = 0.2658E+02
      dedx = 0.3952E-02
      call makeMaterial
      name = "CFIB"
      a = 0.1201E+02
      z = 0.6000E+01
      density = 0.2265E+01
      radlen = 0.2140E+02
      abslen = 0.4990E+02
      collen = 0.2660E+02
      dedx = 0.4030E-02
      call makeMaterial
      name = "MXSO"
      a = 0.5572E+02
      z = 0.2550E+02
      density = 0.6026E+01
      radlen = 0.2372E+01
      abslen = 0.2481E+02
      collen = 0.0000E+00
      dedx = 0.9875E-02
      call makeMaterial
      name = "ARGON"
      a = 39.95
      z = 18.0
      density = 1.78E-3
      radlen = 11706.
      abslen = 5.E+4
      collen = 1.E+16
      dedx = 2.70E-6
      call makeMaterial
      name = "C2H6"
      a = 30.8
      z = 18.0
      density = 1.36E-3
      radlen = 34035.
      abslen = 5.E+4
      collen = 1.E+16
      dedx = 3.05E-6
      call makeMaterial
      name = "MIX1"
      a = 0.5077E+02
      z = 0.2368E+02
      density = 0.4451E+01
      radlen = 0.3380E+01
      abslen = 0.3181E+02
      collen = 0.0000E+00
      dedx = 0.1053E-01
      call makeMaterial
      name = "TUNGS"
      a = 0.1839E+03
      z = 0.7400E+02
      density = 0.1930E+02
      radlen = 0.3500E+00
      abslen = 0.9585E+01
      collen = 0.5715E+01
      dedx = 0.2210E-01
      call makeMaterial
      name = "KAPTON"
      a = 0.1200E+02
      z = 0.6000E+01
      density = 0.1420E+01
      radlen = 0.1000E+17
      abslen = 0.1000E+17
      collen = 0.1000E+17
      dedx = 0.1000E-15
      call makeMaterial
      name = "LEAD"
      a = 0.2072E+03
      z = 0.8200E+02
      density = 0.1135E+02
      radlen = 0.5600E+00
      abslen = 0.1709E+02
      collen = 0.1024E+02
      dedx = 0.1275E-01
      call makeMaterial
      name = "PETHYL"
      a = 0.1200E+02
      z = 0.6000E+01
      density = 0.9200E+00
      radlen = 0.4790E+02
      abslen = 0.8565E+02
      collen = 0.6185E+02
      dedx = 0.1910E-02
      call makeMaterial
      name = "GLASS"
      a = 0.5953E+02
      z = 0.2958E+02
      density = 0.2230E+01
      radlen = 0.1270E+02
      abslen = 0.4377E+02
      collen = 0.2969E+02
      dedx = 0.3780E-02
      call makeMaterial
      name = "CF4"
      a = 88.
      z = 42.
      density = 3.9286E-3
      radlen = 64000.
      abslen = 5.0E+04
      collen = 1.0E+16
      dedx = 6.75E-06
      call makeMATERIAL
      end

      subroutine MixtureDef
      character*40 name
      integer nmat
      character*40 matnames(5)
      real prop(5)
      common /Mixture/name,nmat,matnames,prop
      name = "ARGON-ETHAN"
      nmat = 2
      matnames(1) = "ARGON"
      matnames(2) = "C2H6"
      matnames(3) = "-"
      matnames(4) = "-"
      matnames(5) = "-"
      prop(1) = 0.65
      prop(2) = 0.35
      prop(3) = 0.00
      prop(4) = 0.00
      prop(5) = 0.00
      call makeMixture
      name = "AR-ETH-CF4"
      nmat = 3
      matnames(1) = "ARGON"
      matnames(2) = "C2H6"
      matnames(3) = "CF4"
      matnames(4) = "-"
      matnames(5) = "-"
      prop(1) = 0.5
      prop(2) = 0.35
      prop(3) = 0.15
      call makeMIXTURE
      name = "GAS-MX"
      nmat = 3
      matnames(1) = "ARGON"
      matnames(2) = "C2H6"
      matnames(3) = "CF4"
      matnames(4) = "ALUM"
      matnames(5) = "-"
      prop(1) = 0.500
      prop(2) = 0.35
      prop(3) = 0.15
      call makeMIXTURE
      name = "PBSC"
      nmat = 2
      matnames(1) = "LEAD"
      matnames(2) = "SCIN"
      matnames(3) = "-"
      matnames(4) = "-"
      matnames(5) = "-"
      prop(1) = 0.4
      prop(2) = 0.6
      call makeMIXTURE
      end

      subroutine BPipeDef
      character*40 name
      real rmin
      real rmax
      real z0
      real zlen
      character*40 mat_fill
      real bndrthk(4)
      character*40 matrbnd(4)
      common /BPipe/name,rmin,rmax,z0,zlen,mat_fill,bndrthk,matrbnd
      name = "beam"
      rmin = 0.00
      rmax = 2.50
      z0 = 225.0
      zlen = 450.0
      mat_fill = "VACU"
      bndrthk(1) = 0.00
      bndrthk(2) = 0.05
      bndrthk(3) = 0.00
      bndrthk(4) = 0.00
      matrbnd(1) = "VACU"
      matrbnd(2) = "BE"
      matrbnd(3) = "VACU"
      matrbnd(4) = "VACU"
      call makeBPipe
      end

      subroutine SolenoidDef
      character*40 name
      real bfield
      real rmin
      real rmax
      real z0
      real zlen
      character*40 mat_fill
      real thick_boun(4)
      character*40 mat_boun(4)
      common /Solenoid/name,bfield,rmin,rmax,z0,zlen,mat_fill,
     + thick_boun,mat_boun
      name = "LASS"
      bfield = 2.24000
      rmin = 0.0000
      rmax = 95.0000
      z0 = 210.0
      zlen = 420.0
      mat_fill = "AIR"
      thick_boun(1) = 0.0000
      thick_boun(2) = 2.0000
      thick_boun(3) = 0.0000
      thick_boun(4) = 0.0000
      mat_boun(1) = "----"
      mat_boun(2) = "MXSO"
      mat_boun(3) = "----"
      mat_boun(4) = "----"
      call makeSolenoid
      end

      subroutine EMCalDef
      character*40 name
      character*40 shape
      integer type
      real rmin(2)
      real rmax(2)
      real z0
      real zlen
      character*40 material
      character*40 active
      integer nphi
      integer neta
      real siga_em
      real sigb_em
      real siga_had
      real sigb_had
      real em_had_ratio
      integer nlayers
      common /EMCal/name,shape,type,rmin,rmax,z0,zlen,material,active,
     + nphi,neta,siga_em,sigb_em,siga_had,sigb_had,em_had_ratio,nlayers
      name = "BCAL"
      shape = "TUBE"
      type = 1
      rmin(1) = 65
      rmin(2) = 65
      rmax(1) = 90
      rmax(2) = 90
      z0 = 217.5
      zlen = 405
      material = "PBSC"
      active = "PBSC"
      nphi = 48
      neta = 1
      siga_em = 0.06
      sigb_em = 0.01
      siga_had = 0.3
      sigb_had = 0.3
      em_had_ratio = 4
      nlayers = 1
      call makeEMCal
      end

      subroutine DriftDef
      integer num
      character*40 name
      integer num_anode
      integer num_cathode
      real rmin
      real rmax
      real z0
      real zlen
      character*40 material
      real thick_boun(4)
      character*40 mat_boun(4)
      common /Drift/num,name,num_anode,num_cathode,rmin,rmax,z0,zlen,
     + material,thick_boun,mat_boun
      num = 1
      name = "VTX"
      num_anode = 3
      num_cathode = 0
      rmin = 4.95
      rmax = 5.65
      z0 = 80.0
      zlen = 90.000
      material = "SCIN"
      thick_boun(1) = 0.2000
      thick_boun(2) = 0.0500
      thick_boun(3) = 0.2000
      thick_boun(4) = 0.2000
      mat_boun(1) = "CFIB"
      mat_boun(2) = "CFIB"
      mat_boun(3) = "CFIB"
      mat_boun(4) = "CFIB"
      call makeDrift
      num = 2
      name = "CDC"
      num_anode = 22
      num_cathode = 1
      rmin = 15.0
      rmax = 60.0
      z0 = 117.0
      zlen = 200.0
      material = "GAS-MX"
      thick_boun(1) = 0.2000
      thick_boun(2) = 0.4000
      thick_boun(3) = 2.0000
      thick_boun(4) = 2.0000
      mat_boun(1) = "CFIB"
      mat_boun(2) = "CFIB"
      mat_boun(3) = "ALUM"
      mat_boun(4) = "ALUM"
      call makeDrift
      end

      subroutine LayerDRFAnoDef
      integer det
      integer lyr
      real radius
      real zlen
      real cell_height
      integer nwires
      integer ID_readout
      integer ID_cathode
      real phi0
      real stereo_tau
      real stereo_offset
      real eff_hit
      real eff_dedx
      real siga
      real sigb
      real sigc
      common /LayerDRFAno/det,lyr,radius,zlen,cell_height,nwires,
     + ID_readout,ID_cathode,phi0,stereo_tau,stereo_offset,eff_hit,
     + eff_dedx,siga,sigb,sigc
      det = 1
      lyr = 1
      radius = 5.0
      zlen = 90.0
      cell_height = 1.0
      nwires = 25
      ID_readout = -1
      ID_cathode = 1
      phi0 = 0.
      stereo_tau = -0.10
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = .04
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 1
      lyr = 2
      radius = 5.1
      zlen = 90.0
      cell_height = 1.0
      nwires = 25
      ID_readout = -1
      ID_cathode = 1
      phi0 = 0.
      stereo_tau = 0.0
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = .04
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 1
      lyr = 3
      radius = 5.2
      zlen = 90.0
      cell_height = 1.0
      nwires = 25
      ID_readout = -1
      ID_cathode = 1
      phi0 = 0.
      stereo_tau = 0.10
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = .04
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 1
      radius = 16.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 100
      ID_readout = -1
      ID_cathode = 1
      phi0 = 0.
      stereo_tau = 0.
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 2
      radius = 18.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 113
      ID_readout = -1
      ID_cathode = 0
      phi0 = 0.
      stereo_tau = 0.
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 3
      radius = 20.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 126
      ID_readout = -1
      ID_cathode = 0
      phi0 = 0.
      stereo_tau = 0.104
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 4
      radius = 22.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 139
      ID_readout = -1
      ID_cathode = 0
      phi0 = 0.
      stereo_tau = 0.104
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 5
      radius = 24.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 152
      ID_readout = -1
      ID_cathode = 0
      phi0 = 0.
      stereo_tau = 0.
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 6
      radius = 26.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 165
      ID_readout = -1
      ID_cathode = 0
      phi0 = 0.
      stereo_tau = 0.
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 7
      radius = 28.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 178
      ID_readout = -1
      ID_cathode = 0
      phi0 = 0.
      stereo_tau = 0.
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 8
      radius = 30.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 191
      ID_readout = -1
      ID_cathode = 0
      phi0 = 0.
      stereo_tau = -.104
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 9
      radius = 32.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 204
      ID_readout = -1
      ID_cathode = 0
      phi0 = 0.
      stereo_tau = -.104
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 10
      radius = 34.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 217
      ID_readout = -1
      ID_cathode = 0
      phi0 = 0.
      stereo_tau = 0.
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 11
      radius = 36.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 230
      ID_readout = -1
      ID_cathode = 0
      phi0 = 0.
      stereo_tau = 0.
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 12
      radius = 38.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 243
      ID_readout = -1
      ID_cathode = 0
      phi0 = 0.
      stereo_tau = 0.
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 13
      radius = 40.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 256
      ID_readout = -1
      ID_cathode = 0
      phi0 = 0.
      stereo_tau = 0.104
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 14
      radius = 42.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 269
      ID_readout = -1
      ID_cathode = 0
      phi0 = 0.
      stereo_tau = 0.104
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 15
      radius = 44.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 282
      ID_readout = -1
      ID_cathode = 0
      phi0 = 0.
      stereo_tau = 0.
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 16
      radius = 46.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 295
      ID_readout = -1
      ID_cathode = 2
      phi0 = 0.
      stereo_tau = 0.
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 17
      radius = 48.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 308
      ID_readout = -1
      ID_cathode = 2
      phi0 = 0.
      stereo_tau = 0.
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 18
      radius = 50.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 321
      ID_readout = -1
      ID_cathode = 2
      phi0 = 0.
      stereo_tau = -.104
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 19
      radius = 52.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 334
      ID_readout = -1
      ID_cathode = 2
      phi0 = 0.
      stereo_tau = -.104
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 20
      radius = 54.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 347
      ID_readout = -1
      ID_cathode = 2
      phi0 = 0.
      stereo_tau = 0.
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 21
      radius = 56.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 360
      ID_readout = -1
      ID_cathode = 2
      phi0 = 0.
      stereo_tau = 0.
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      det = 2
      lyr = 22
      radius = 58.0
      zlen = 200.0
      cell_height = 1.0
      nwires = 373
      ID_readout = -1
      ID_cathode = 2
      phi0 = 0.
      stereo_tau = 0.
      stereo_offset = 0.
      eff_hit = .96
      eff_dedx = 0.96
      siga = 0.020
      sigb = 0.0
      sigc = 0.0
      call makeLayerDRFAno
      end

      subroutine OffsetDRFAnoDef
      integer det
      integer lyr
      real cell_offset
      real sag
      real offset(3)
      real dircos(3)
      common /OffsetDRFAno/det,lyr,cell_offset,sag,offset,dircos
      det = 1
      lyr = 1
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 1
      lyr = 2
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 1
      lyr = 3
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 1
      lyr = 4
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 1
      lyr = 5
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 1
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 2
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 3
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 4
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 5
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 6
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 7
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 8
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 9
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 10
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 11
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 12
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 13
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 14
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 15
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 16
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 17
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 18
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 19
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 20
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 21
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      det = 2
      lyr = 22
      cell_offset = 0.0
      sag = 0.0
      offset(1) = 0.0
      offset(2) = 0.0
      offset(3) = 0.0
      dircos(1) = 0.0
      dircos(2) = 0.0
      dircos(3) = 0.0
      call makeOffsetDRFAno
      end

      subroutine AbsorberDef
      character*40 name
      character*40 shape
      integer type
      real rmin(2)
      real rmax(2)
      real z0
      real zlen
      character*40 material
      common /Absorber/name,shape,type,rmin,rmax,z0,zlen,material
      name = "SHIT01"
      shape = "TUBE"
      type = 41
      rmin(1) = 17.00
      rmin(2) = 17.00
      rmax(1) = 17.04
      rmax(2) = 17.04
      z0 = 117.0
      zlen = 200.0
      material = "ALUM"
      call makeAbsorber
      name = "SHIT03"
      shape = "TUBE"
      type = 41
      rmin(1) = 21.00
      rmin(2) = 21.00
      rmax(1) = 21.04
      rmax(2) = 21.04
      z0 = 117.0
      zlen = 200.0
      material = "ALUM"
      call makeAbsorber
      name = "SHIT05"
      shape = "TUBE"
      type = 41
      rmin(1) = 25.00
      rmin(2) = 25.00
      rmax(1) = 25.04
      rmax(2) = 25.04
      z0 = 117.0
      zlen = 200.0
      material = "ALUM"
      call makeAbsorber
      name = "SHIT07"
      shape = "TUBE"
      type = 41
      rmin(1) = 29.00
      rmin(2) = 29.00
      rmax(1) = 29.04
      rmax(2) = 29.04
      z0 = 117.0
      zlen = 200.0
      material = "ALUM"
      call makeAbsorber
      name = "SHIT08"
      shape = "TUBE"
      type = 41
      rmin(1) = 31.00
      rmin(2) = 31.00
      rmax(1) = 31.04
      rmax(2) = 31.04
      z0 = 117.0
      zlen = 200.0
      material = "ALUM"
      call makeAbsorber
      name = "SHIT10"
      shape = "TUBE"
      type = 41
      rmin(1) = 35.00
      rmin(2) = 35.00
      rmax(1) = 35.04
      rmax(2) = 35.04
      z0 = 117.0
      zlen = 200.0
      material = "ALUM"
      call makeAbsorber
      name = "SHIT11"
      shape = "TUBE"
      type = 41
      rmin(1) = 37.00
      rmin(2) = 37.00
      rmax(1) = 37.04
      rmax(2) = 37.04
      z0 = 117.0
      zlen = 200.0
      material = "ALUM"
      call makeAbsorber
      name = "SHIT13"
      shape = "TUBE"
      type = 41
      rmin(1) = 41.00
      rmin(2) = 41.00
      rmax(1) = 41.04
      rmax(2) = 41.04
      z0 = 117.0
      zlen = 200.0
      material = "ALUM"
      call makeAbsorber
      name = "SHIT15"
      shape = "TUBE"
      type = 41
      rmin(1) = 45.00
      rmin(2) = 45.00
      rmax(1) = 45.04
      rmax(2) = 45.04
      z0 = 117.0
      zlen = 200.0
      material = "ALUM"
      call makeAbsorber
      name = "SHIT17"
      shape = "TUBE"
      type = 41
      rmin(1) = 49.00
      rmin(2) = 49.00
      rmax(1) = 49.04
      rmax(2) = 49.04
      z0 = 117.0
      zlen = 200.0
      material = "ALUM"
      call makeAbsorber
      name = "SHIT19"
      shape = "TUBE"
      type = 41
      rmin(1) = 53.00
      rmin(2) = 53.00
      rmax(1) = 53.04
      rmax(2) = 53.04
      z0 = 117.0
      zlen = 200.0
      material = "ALUM"
      call makeAbsorber
      name = "SHIT21"
      shape = "TUBE"
      type = 41
      rmin(1) = 57.00
      rmin(2) = 57.00
      rmax(1) = 57.04
      rmax(2) = 57.04
      z0 = 117.0
      zlen = 200.0
      material = "ALUM"
      call makeAbsorber
      name = "SHL7"
      shape = "TUBE"
      type = 41
      rmin(1) = 60.0
      rmin(2) = 60.0
      rmax(1) = 61.0
      rmax(2) = 61.0
      z0 = 230.0
      zlen = 12.0
      material = "ALUM"
      call makeAbsorber
      name = "SHL8"
      shape = "TUBE"
      type = 41
      rmin(1) = 60.0
      rmin(2) = 60.0
      rmax(1) = 61.0
      rmax(2) = 61.0
      z0 = 282.0
      zlen = 12.0
      material = "ALUM"
      call makeAbsorber
      name = "SHL9"
      shape = "TUBE"
      type = 41
      rmin(1) = 60.0
      rmin(2) = 60.0
      rmax(1) = 61.0
      rmax(2) = 61.0
      z0 = 338.0
      zlen = 12.0
      material = "ALUM"
      call makeAbsorber
      name = "SHL10"
      shape = "TUBE"
      type = 41
      rmin(1) = 60.0
      rmin(2) = 60.0
      rmax(1) = 61.0
      rmax(2) = 61.0
      z0 = 394.0
      zlen = 12.0
      material = "ALUM"
      call makeAbsorber
      name = "CTOF"
      shape = "TUBE"
      type = 41
      rmin(1) = 65.0
      rmin(2) = 65.0
      rmax(1) = 66.0
      rmax(2) = 66.0
      z0 = 217.5
      zlen = 405.0
      material = "SCIN"
      call makeAbsorber
      end

      subroutine AbsorberBoxDef
      character*40 name
      character*40 shape
      integer type
      real xlimit(2)
      real ylimit(2)
      real xlimit_gap(2)
      real ylimit_gap(2)
      real z0
      real zlen
      character*40 material
      common /AbsorberBox/name,shape,type,xlimit,ylimit,xlimit_gap,
     + ylimit_gap,z0,zlen,material
      name = "FTOF"
      shape = "BOX"
      type = 42
      xlimit(1) = -125.0
      xlimit(2) = 125.0
      ylimit(1) = -125.0
      ylimit(2) = 125.0
      xlimit_gap(1) = 0.0
      xlimit_gap(2) = 0.0
      ylimit_gap(1) = 0.0
      ylimit_gap(2) = 0.0
      z0 = 566.27
      zlen = 2.54
      material = "SCIN"
      call makeAbsorberBox
      name = "CERENKOV"
      shape = "BOX"
      type = 42
      xlimit(1) = -95.0
      xlimit(2) = 95.0
      ylimit(1) = -95.0
      ylimit(2) = 95.0
      xlimit_gap(1) = 0.0
      xlimit_gap(2) = 0.0
      ylimit_gap(1) = 0.0
      ylimit_gap(2) = 0.0
      z0 = 420.5
      zlen = 1.0
      material = "SCIN"
      call makeAbsorberBox
      end

      subroutine SiDiskDef
      integer num
      character*40 name
      integer nlyr
      real zpos
      common /SiDisk/num,name,nlyr,zpos
      num = 1
      name = "FDC1"
      nlyr = 9
      zpos = 227.0
      call makeSiDisk
      num = 2
      name = "FDC2"
      nlyr = 9
      zpos = 233.0
      call makeSiDisk
      num = 3
      name = "FDC3"
      nlyr = 9
      zpos = 279.0
      call makeSiDisk
      num = 4
      name = "FDC4"
      nlyr = 9
      zpos = 285.0
      call makeSiDisk
      num = 5
      name = "FDC5"
      nlyr = 9
      zpos = 335.0
      call makeSiDisk
      num = 6
      name = "FDC6"
      nlyr = 9
      zpos = 341.0
      call makeSiDisk
      num = 7
      name = "FDC7"
      nlyr = 9
      zpos = 391.0
      call makeSiDisk
      num = 8
      name = "FDC8"
      nlyr = 9
      zpos = 397.0
      call makeSiDisk
      end

      subroutine LayerSiDiDef
      integer det
      integer lyr
      character*40 mat
      integer nwed
      real z_local
      real thick
      real rmin
      real rmax
      real phi(2)
      real dphi
      integer type
      common /LayerSiDi/det,lyr,mat,nwed,z_local,thick,rmin,rmax,phi,
     + dphi,type
      det = 1
      lyr = 1
      mat = "GAS-MX"
      nwed = 1
      z_local = -2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 1
      lyr = 2
      mat = "GAS-MX"
      nwed = 1
      z_local = 0.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 1
      lyr = 3
      mat = "GAS-MX"
      nwed = 1
      z_local = 2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 1
      lyr = 4
      mat = "PETHYL"
      nwed = 1
      z_local = -2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 1
      lyr = 5
      mat = "PETHYL"
      nwed = 1
      z_local = -1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 1
      lyr = 6
      mat = "PETHYL"
      nwed = 1
      z_local = -0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 1
      lyr = 7
      mat = "PETHYL"
      nwed = 1
      z_local = 0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 1
      lyr = 8
      mat = "PETHYL"
      nwed = 1
      z_local = 1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 1
      lyr = 9
      mat = "PETHYL"
      nwed = 1
      z_local = 2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 2
      lyr = 1
      mat = "GAS-MX"
      nwed = 1
      z_local = -2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 2
      lyr = 2
      mat = "GAS-MX"
      nwed = 1
      z_local = 0.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 2
      lyr = 3
      mat = "GAS-MX"
      nwed = 1
      z_local = 2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 2
      lyr = 4
      mat = "PETHYL"
      nwed = 1
      z_local = -2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 2
      lyr = 5
      mat = "PETHYL"
      nwed = 1
      z_local = -1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 2
      lyr = 6
      mat = "PETHYL"
      nwed = 1
      z_local = -0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 2
      lyr = 7
      mat = "PETHYL"
      nwed = 1
      z_local = 0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 2
      lyr = 8
      mat = "PETHYL"
      nwed = 1
      z_local = 1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 2
      lyr = 9
      mat = "PETHYL"
      nwed = 1
      z_local = 2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 3
      lyr = 1
      mat = "GAS-MX"
      nwed = 1
      z_local = -2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 3
      lyr = 2
      mat = "GAS-MX"
      nwed = 1
      z_local = 0.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 3
      lyr = 3
      mat = "GAS-MX"
      nwed = 1
      z_local = 2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 3
      lyr = 4
      mat = "PETHYL"
      nwed = 1
      z_local = -2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 3
      lyr = 5
      mat = "PETHYL"
      nwed = 1
      z_local = -1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 3
      lyr = 6
      mat = "PETHYL"
      nwed = 1
      z_local = -0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 3
      lyr = 7
      mat = "PETHYL"
      nwed = 1
      z_local = 0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 3
      lyr = 8
      mat = "PETHYL"
      nwed = 1
      z_local = 1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 3
      lyr = 9
      mat = "PETHYL"
      nwed = 1
      z_local = 2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 4
      lyr = 1
      mat = "GAS-MX"
      nwed = 1
      z_local = -2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 4
      lyr = 2
      mat = "GAS-MX"
      nwed = 1
      z_local = 0.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 4
      lyr = 3
      mat = "GAS-MX"
      nwed = 1
      z_local = 2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 4
      lyr = 4
      mat = "PETHYL"
      nwed = 1
      z_local = -2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 4
      lyr = 5
      mat = "PETHYL"
      nwed = 1
      z_local = -1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 4
      lyr = 6
      mat = "PETHYL"
      nwed = 1
      z_local = -0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 4
      lyr = 7
      mat = "PETHYL"
      nwed = 1
      z_local = 0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 4
      lyr = 8
      mat = "PETHYL"
      nwed = 1
      z_local = 1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 4
      lyr = 9
      mat = "PETHYL"
      nwed = 1
      z_local = 2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 5
      lyr = 1
      mat = "GAS-MX"
      nwed = 1
      z_local = -2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 5
      lyr = 2
      mat = "GAS-MX"
      nwed = 1
      z_local = 0.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 5
      lyr = 3
      mat = "GAS-MX"
      nwed = 1
      z_local = 2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 5
      lyr = 4
      mat = "PETHYL"
      nwed = 1
      z_local = -2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 5
      lyr = 5
      mat = "PETHYL"
      nwed = 1
      z_local = -1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 5
      lyr = 6
      mat = "PETHYL"
      nwed = 1
      z_local = -0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 5
      lyr = 7
      mat = "PETHYL"
      nwed = 1
      z_local = 0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 5
      lyr = 8
      mat = "PETHYL"
      nwed = 1
      z_local = 1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 5
      lyr = 9
      mat = "PETHYL"
      nwed = 1
      z_local = 2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 6
      lyr = 1
      mat = "GAS-MX"
      nwed = 1
      z_local = -2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 6
      lyr = 2
      mat = "GAS-MX"
      nwed = 1
      z_local = 0.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 6
      lyr = 3
      mat = "GAS-MX"
      nwed = 1
      z_local = 2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 6
      lyr = 4
      mat = "PETHYL"
      nwed = 1
      z_local = -2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 6
      lyr = 5
      mat = "PETHYL"
      nwed = 1
      z_local = -1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 6
      lyr = 6
      mat = "PETHYL"
      nwed = 1
      z_local = -0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 6
      lyr = 7
      mat = "PETHYL"
      nwed = 1
      z_local = 0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 6
      lyr = 8
      mat = "PETHYL"
      nwed = 1
      z_local = 1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 6
      lyr = 9
      mat = "PETHYL"
      nwed = 1
      z_local = 2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 7
      lyr = 1
      mat = "GAS-MX"
      nwed = 1
      z_local = -2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 7
      lyr = 2
      mat = "GAS-MX"
      nwed = 1
      z_local = 0.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 7
      lyr = 3
      mat = "GAS-MX"
      nwed = 1
      z_local = 2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 7
      lyr = 4
      mat = "PETHYL"
      nwed = 1
      z_local = -2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 7
      lyr = 5
      mat = "PETHYL"
      nwed = 1
      z_local = -1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 7
      lyr = 6
      mat = "PETHYL"
      nwed = 1
      z_local = -0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 7
      lyr = 7
      mat = "PETHYL"
      nwed = 1
      z_local = 0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 7
      lyr = 8
      mat = "PETHYL"
      nwed = 1
      z_local = 1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 7
      lyr = 9
      mat = "PETHYL"
      nwed = 1
      z_local = 2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 8
      lyr = 1
      mat = "GAS-MX"
      nwed = 1
      z_local = -2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 8
      lyr = 2
      mat = "GAS-MX"
      nwed = 1
      z_local = 0.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 8
      lyr = 3
      mat = "GAS-MX"
      nwed = 1
      z_local = 2.0
      thick = 2.00
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 8
      lyr = 4
      mat = "PETHYL"
      nwed = 1
      z_local = -2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 8
      lyr = 5
      mat = "PETHYL"
      nwed = 1
      z_local = -1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 8
      lyr = 6
      mat = "PETHYL"
      nwed = 1
      z_local = -0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 8
      lyr = 7
      mat = "PETHYL"
      nwed = 1
      z_local = 0.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 8
      lyr = 8
      mat = "PETHYL"
      nwed = 1
      z_local = 1.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      det = 8
      lyr = 9
      mat = "PETHYL"
      nwed = 1
      z_local = 2.5
      thick = 0.008
      rmin = 3.50
      rmax = 60.00
      phi(1) = 0.
      phi(2) = 6.283
      dphi = 6.283
      type = 3
      call makeLayerSiDi
      end

      subroutine WedgeDef
      character*40 speci
      integer det
      integer lyr
      integer nwed
      integer nstrip
      real c0_r
      real c0_phi
      real pitch
      real stereo
      real eff_hit
      real siga
      real sigb
      real sigc
      common /Wedge/speci,det,lyr,nwed,nstrip,c0_r,c0_phi,pitch,stereo,
     + eff_hit,siga,sigb,sigc
      speci = "ALL"
      det = 1
      lyr = 1
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.000
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 1
      lyr = 2
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 1
      lyr = 3
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 1
      lyr = 4
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 1
      lyr = 5
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 1
      lyr = 6
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 1
      lyr = 7
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 1
      lyr = 8
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 1
      lyr = 9
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 2
      lyr = 1
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.000
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 2
      lyr = 2
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 2
      lyr = 3
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 2
      lyr = 4
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 2
      lyr = 5
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 2
      lyr = 6
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 2
      lyr = 7
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 2
      lyr = 8
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 2
      lyr = 9
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 3
      lyr = 1
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.000
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 3
      lyr = 2
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 3
      lyr = 3
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 3
      lyr = 4
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 3
      lyr = 5
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 3
      lyr = 6
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 3
      lyr = 7
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 3
      lyr = 8
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 3
      lyr = 9
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 4
      lyr = 1
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.000
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 4
      lyr = 2
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 4
      lyr = 3
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 4
      lyr = 4
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 4
      lyr = 5
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 4
      lyr = 6
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 4
      lyr = 7
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 4
      lyr = 8
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 4
      lyr = 9
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 5
      lyr = 1
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.000
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 5
      lyr = 2
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 5
      lyr = 3
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 5
      lyr = 4
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 5
      lyr = 5
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 5
      lyr = 6
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 5
      lyr = 7
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 5
      lyr = 8
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 5
      lyr = 9
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 6
      lyr = 1
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.000
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 6
      lyr = 2
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 6
      lyr = 3
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 6
      lyr = 4
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 6
      lyr = 5
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 6
      lyr = 6
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 6
      lyr = 7
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 6
      lyr = 8
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 6
      lyr = 9
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 7
      lyr = 1
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.000
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 7
      lyr = 2
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 7
      lyr = 3
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 7
      lyr = 4
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 7
      lyr = 5
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 7
      lyr = 6
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 7
      lyr = 7
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 7
      lyr = 8
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 7
      lyr = 9
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 8
      lyr = 1
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.000
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 8
      lyr = 2
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 8
      lyr = 3
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.96
      siga = 0.015
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 8
      lyr = 4
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 8
      lyr = 5
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -1.0472
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 8
      lyr = 6
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 8
      lyr = 7
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 8
      lyr = 8
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = 0.0000
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      speci = "ALL"
      det = 8
      lyr = 9
      nwed = 1
      nstrip = 119
      c0_r = 0.000
      c0_phi = 0.00
      pitch = 1.0
      stereo = -2.0944
      eff_hit = 0.90
      siga = 1.000
      sigb = 0.
      sigc = 0.
      call makeWedge
      end

      subroutine HitsOnTrackDef
      integer all
      integer z
      integer svx
      common /HitsOnTrack/all,z,svx
      all = 4
      z = 0
      svx = 0
      call makeHitsOnTrack
      end

      subroutine BeamVrtxDef
      real xyz(3)
      real sig(3)
      common /BeamVrtx/xyz,sig
      xyz(1) = 0.
      xyz(2) = 0.
      xyz(3) = 65
      sig(1) = 0.30
      sig(2) = 0.30
      sig(3) = 15
      call makeBeamVrtx
      end

      subroutine makeEMCal
      character*40 name
      character*40 shape
      integer type
      real rmin(2)
      real rmax(2)
      real z0
      real zlen
      character*40 material
      character*40 active
      integer nphi
      integer neta
      real siga_em
      real sigb_em
      real siga_had
      real sigb_had
      real em_had_ratio
      integer nlayers
      common /EMCal/name,shape,type,rmin,rmax,z0,zlen,material,active,
     + nphi,neta,siga_em,sigb_em,siga_had,sigb_had,em_had_ratio,nlayers

CC---> add the appropriate GEANT calls here

      end

      subroutine makeAbsorberBox
      character*40 name
      character*40 shape
      integer type
      real xlimit(2)
      real ylimit(2)
      real xlimit_gap(2)
      real ylimit_gap(2)
      real z0
      real zlen
      character*40 material
      common /AbsorberBox/name,shape,type,xlimit,ylimit,xlimit_gap,
     + ylimit_gap,z0,zlen,material

CC---> add the appropriate GEANT calls here

      end

      subroutine makeGeometry
      call EMCalDef
      call AbsorberBoxDef
      end
