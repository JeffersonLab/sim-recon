# PS/PSC Timing Offsets
Follow these steps after running `hd_root` with the **PS_timing** plugin; The JANA configuration parameter `PSCHit:DELTA_T_ADC_TDC_MAX` should initially be set to a large value, such as 500 ns.

1. `root -b -q 'fits.C("hd_root.root",true)'`
2. `root -b -q 'offsets.C("fits-csv")'`
3. `ccdb add /PHOTON_BEAM/pair_spectrometer/base_time_offset -v default -r $RunNo-$RunNo offsets/base_time_offset.txt`
4. `ccdb add /PHOTON_BEAM/pair_spectrometer/coarse/tdc_timing_offsets -v default -r $RunNo-$RunNo offsets/tdc_timing_offsets_psc.txt`
5. `ccdb add /PHOTON_BEAM/pair_spectrometer/coarse/adc_timing_offsets -v default -r $RunNo-$RunNo offsets/adc_timing_offsets_psc.txt`
6. `ccdb add /PHOTON_BEAM/pair_spectrometer/fine/adc_timing_offsets -v default -r $RunNo-$RunNo offsets/adc_timing_offsets_ps.txt`
