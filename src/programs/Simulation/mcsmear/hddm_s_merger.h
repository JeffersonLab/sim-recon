//
// hddm_s_merger.h - Utility class for merging hits from two hddm_s element lists
//
// author: richard.t.jones at uconn.edu
// version: march 20, 2017

#ifndef _HDDM_S_MERGER_H_
#define _HDDM_S_MERGER_H_

#include <HDDM/hddm_s.hpp>

namespace hddm_s_merger {
   double get_t_shift_ns();
   void set_t_shift_ns(double dt_ns);

   // hits merging / truncation parameters for the CDC
   int get_cdc_max_hits();
   void set_cdc_max_hits(int maxhits);
   double get_cdc_integration_window_ns();
   void set_cdc_integration_window_ns(double dt_ns);

   // hits merging / truncation parameters for the FDC
   int get_fdc_wires_max_hits();
   void set_fdc_wires_max_hits(int maxhits);
   double get_fdc_wires_min_delta_t_ns();
   void set_fdc_wires_min_delta_t_ns(double dt_ns);
   int get_fdc_strips_max_hits();
   void set_fdc_strips_max_hits(int maxhits);
   double get_fdc_strips_integraton_window_ns();
   void set_fdc_strips_integration_window_ns(double dt_ns);

   // hits merging / truncation parameters for the STC
   int get_stc_adc_max_hits();
   void set_stc_adc_max_hits(int maxhits);
   int get_stc_tdc_max_hits();
   void set_stc_tdc_max_hits(int maxhits);
   double get_stc_min_delta_t_ns();
   void set_stc_min_delta_t_ns(double dt_ns);
   double get_stc_integration_window_ns();
   void set_stc_integration_window_ns(double dt_ns);

   // hits merging / truncation parameters for the BCAL
   int get_bcal_adc_max_hits();
   void set_bcal_adc_max_hits(int maxhits);
   int get_bcal_tdc_max_hits();
   void set_bcal_tdc_max_hits(int maxhits);
   double get_bcal_min_delta_t_ns();
   void set_bcal_min_delta_t_ns(double dt_ns);
   double get_bcal_integration_window_ns();
   void set_bcal_integration_window_ns(double dt_ns);
   double get_bcal_fadc_counts_per_ns();
   void set_bcal_fadc_counts_per_ns(double slope);
   double get_bcal_tdc_counts_per_ns();
   void set_bcal_tdc_counts_per_ns(double slope);

   // hits merging / truncation parameters for the TOF
   int get_ftof_adc_max_hits();
   void set_ftof_adc_max_hits(int maxhits);
   int get_ftof_tdc_max_hits();
   void set_ftof_tdc_max_hits(int maxhits);
   double get_ftof_min_delta_t_ns();
   void set_ftof_min_delta_t_ns(double dt_ns);
   double get_ftof_integration_window_ns();
   void set_ftof_integration_window_ns(double dt_ns);

   // hits merging / truncation parameters for the FCAL
   int get_fcal_max_hits();
   void set_fcal_max_hits(int maxhits);
   double get_fcal_min_delta_t_ns();
   void set_fcal_min_delta_t_ns(double dt_ns);
   double get_fcal_integration_window_ns();
   void set_fcal_integration_window_ns(double dt_ns);

   // hits merging / truncation parameters for the CCAL
   int get_ccal_max_hits();
   void set_ccal_max_hits(int maxhits);
   double get_ccal_min_delta_t_ns();
   void set_ccal_min_delta_t_ns(double dt_ns);
   double get_ccal_integration_window_ns();
   void set_ccal_integration_window_ns(double dt_ns);

   // hits merging / truncation parameters for the PS
   int get_ps_max_hits();
   void set_ps_max_hits(int maxhits);
   double get_ps_integration_window_ns();
   void set_ps_integration_window_ns(double dt_ns);
   int get_psc_adc_max_hits();
   void set_psc_adc_max_hits(int maxhits);
   int get_psc_tdc_max_hits();
   void set_psc_tdc_max_hits(int maxhits);
   double get_psc_min_delta_t_ns();
   void set_psc_min_delta_t_ns(double dt_ns);
   double get_psc_integration_window_ns();
   void set_psc_integration_window_ns(double dt_ns);

   // hits merging / truncation parameters for the TAGM/TAGH
   int get_tag_adc_max_hits();
   void set_tag_adc_max_hits(int maxhits);
   int get_tag_tdc_max_hits();
   void set_tag_tdc_max_hits(int maxhits);
   double get_tag_min_delta_t_ns();
   void set_tag_min_delta_t_ns(double dt_ns);
   double get_tag_integration_window_ns();
   void set_tag_integration_window_ns(double dt_ns);

   // hits merging / truncation parameters for the TPOL
   int get_tpol_max_hits();
   void set_tpol_max_hits(int maxhits);
   double get_tpol_integration_window_ns();
   void set_tpol_integration_window_ns(double dt_ns);
 
   // hits merging / truncation parameters for the FWMPC
   int get_fmwpc_max_hits();
   void set_fmwpc_max_hits(int maxhits);
   double get_fmwpc_min_delta_t_ns();
   void set_fmwpc_min_delta_t_ns(double dt_ns);

   void truncate_hits(hddm_s::HDDM &record);
   void truncate_cdc_hits(hddm_s::CdcStrawHitList &hits);
   void truncate_fdc_wire_hits(hddm_s::FdcAnodeHitList &hits);
   void truncate_fdc_strip_hits(hddm_s::FdcCathodeHitList &hits);
   void truncate_stc_hits(hddm_s::StcHitList &hits);
   void truncate_bcal_adc_hits(hddm_s::BcalfADCHitList &hits);
   void truncate_bcal_tdc_hits(hddm_s::BcalTDCHitList &hits);
   void truncate_bcal_adc_digihits(hddm_s::BcalfADCDigiHitList &hits);
   void truncate_bcal_tdc_digihits(hddm_s::BcalTDCDigiHitList &hits);
   void truncate_ftof_hits(hddm_s::FtofHitList &hits);
   void truncate_fcal_hits(hddm_s::FcalHitList &hits);
   void truncate_ccal_hits(hddm_s::CcalHitList &hits);
   void truncate_tag_hits(hddm_s::TaggerHitList &hits);
   void truncate_ps_hits(hddm_s::PsHitList &hits);
   void truncate_psc_hits(hddm_s::PscHitList &hits);
   void truncate_tpol_hits(hddm_s::TpolHitList &hits);
   void truncate_fmwpc_hits(hddm_s::FmwpcHitList &hits);
}

hddm_s::HDDM &operator+=(hddm_s::HDDM &dst, hddm_s::HDDM &src);
hddm_s::PhysicsEventList &operator+=(hddm_s::PhysicsEventList &dst,
                                     hddm_s::PhysicsEventList &src);
hddm_s::HitViewList &operator+=(hddm_s::HitViewList &dst,
                                hddm_s::HitViewList &src);
hddm_s::CentralDCList &operator+=(hddm_s::CentralDCList &dst,
                                  hddm_s::CentralDCList &src);
hddm_s::CdcStrawList &operator+=(hddm_s::CdcStrawList &dst,
                                 hddm_s::CdcStrawList &src);
hddm_s::CdcStrawHitList &operator+=(hddm_s::CdcStrawHitList &dst,
                                    hddm_s::CdcStrawHitList &src);
hddm_s::ForwardDCList &operator+=(hddm_s::ForwardDCList &dst,
                                  hddm_s::ForwardDCList &src);
hddm_s::FdcChamberList &operator+=(hddm_s::FdcChamberList &dst,
                                   hddm_s::FdcChamberList &src);
hddm_s::FdcAnodeWireList &operator+=(hddm_s::FdcAnodeWireList &dst,
                                     hddm_s::FdcAnodeWireList &src);
hddm_s::FdcAnodeHitList &operator+=(hddm_s::FdcAnodeHitList &dst,
                                    hddm_s::FdcAnodeHitList &src);
hddm_s::FdcCathodeStripList &operator+=(hddm_s::FdcCathodeStripList &dst,
                                        hddm_s::FdcCathodeStripList &src);
hddm_s::FdcCathodeHitList &operator+=(hddm_s::FdcCathodeHitList &dst,
                                      hddm_s::FdcCathodeHitList &src);
hddm_s::StartCntrList &operator+=(hddm_s::StartCntrList &dst,
                                  hddm_s::StartCntrList &src);
hddm_s::StcPaddleList &operator+=(hddm_s::StcPaddleList &dst,
                                  hddm_s::StcPaddleList &src);
hddm_s::StcHitList &operator+=(hddm_s::StcHitList &dst,
                               hddm_s::StcHitList &src);
hddm_s::BarrelEMcalList &operator+=(hddm_s::BarrelEMcalList &dst,
                                    hddm_s::BarrelEMcalList &src);
hddm_s::BcalCellList &operator+=(hddm_s::BcalCellList &dst,
                                 hddm_s::BcalCellList &src);
hddm_s::BcalfADCHitList &operator+=(hddm_s::BcalfADCHitList &dst,
                                    hddm_s::BcalfADCHitList &src);
hddm_s::BcalTDCHitList &operator+=(hddm_s::BcalTDCHitList &dst,
                                   hddm_s::BcalTDCHitList &src);
hddm_s::BcalfADCDigiHitList &operator+=(hddm_s::BcalfADCDigiHitList &dst,
                                        hddm_s::BcalfADCDigiHitList &src);
hddm_s::BcalTDCDigiHitList &operator+=(hddm_s::BcalTDCDigiHitList &dst,
                                        hddm_s::BcalTDCDigiHitList &src);
hddm_s::ForwardTOFList &operator+=(hddm_s::ForwardTOFList &dst,
                                   hddm_s::ForwardTOFList &src);
hddm_s::FtofCounterList &operator+=(hddm_s::FtofCounterList &dst,
                                    hddm_s::FtofCounterList &src);
hddm_s::FtofHitList &operator+=(hddm_s::FtofHitList &dst,
                                hddm_s::FtofHitList &src);
hddm_s::ForwardEMcalList &operator+=(hddm_s::ForwardEMcalList &dst,
                                     hddm_s::ForwardEMcalList &src);
hddm_s::FcalBlockList &operator+=(hddm_s::FcalBlockList &dst,
                                  hddm_s::FcalBlockList &src);
hddm_s::FcalHitList &operator+=(hddm_s::FcalHitList &dst,
                                hddm_s::FcalHitList &src);
hddm_s::ComptonEMcalList &operator+=(hddm_s::ComptonEMcalList &dst,
                                     hddm_s::ComptonEMcalList &src);
hddm_s::CcalBlockList &operator+=(hddm_s::CcalBlockList &dst,
                                  hddm_s::CcalBlockList &src);
hddm_s::CcalHitList &operator+=(hddm_s::CcalHitList &dst,
                                hddm_s::CcalHitList &src);
hddm_s::TaggerList &operator+=(hddm_s::TaggerList &dst,
                               hddm_s::TaggerList &src);
hddm_s::MicroChannelList &operator+=(hddm_s::MicroChannelList &dst,
                                     hddm_s::MicroChannelList &src);
hddm_s::HodoChannelList &operator+=(hddm_s::HodoChannelList &dst,
                                    hddm_s::HodoChannelList &src);
hddm_s::TaggerHitList &operator+=(hddm_s::TaggerHitList &dst,
                                  hddm_s::TaggerHitList &src);
hddm_s::PairSpectrometerFineList &operator+=(
                                 hddm_s::PairSpectrometerFineList &dst,
                                 hddm_s::PairSpectrometerFineList &src);
hddm_s::PsTileList &operator+=(hddm_s::PsTileList &dst,
                               hddm_s::PsTileList &src);
hddm_s::PsHitList &operator+=(hddm_s::PsHitList &dst,
                              hddm_s::PsHitList &src);
hddm_s::PairSpectrometerCoarseList &operator+=(
                              hddm_s::PairSpectrometerCoarseList &dst,
                              hddm_s::PairSpectrometerCoarseList &src);
hddm_s::PscPaddleList &operator+=(hddm_s::PscPaddleList &dst,
                                  hddm_s::PscPaddleList &src);
hddm_s::PscHitList &operator+=(hddm_s::PscHitList &dst,
                               hddm_s::PscHitList &src);
hddm_s::TripletPolarimeterList &operator+=(
                                  hddm_s::TripletPolarimeterList &dst,
                                  hddm_s::TripletPolarimeterList &src);
hddm_s::TpolSectorList &operator+=(hddm_s::TpolSectorList &dst,
                                   hddm_s::TpolSectorList &src);
hddm_s::TpolHitList &operator+=(hddm_s::TpolHitList &dst,
                                hddm_s::TpolHitList &src);
hddm_s::ForwardMWPCList &operator+=(hddm_s::ForwardMWPCList &dst,
                                    hddm_s::ForwardMWPCList &src);
hddm_s::FmwpcChamberList &operator+=(hddm_s::FmwpcChamberList &dst,
                                     hddm_s::FmwpcChamberList &src);
hddm_s::FmwpcHitList &operator+=(hddm_s::FmwpcHitList &dst,
                                 hddm_s::FmwpcHitList &src);

#endif

