
import sbms

Import('*')

subdirs = ['genr8', 'genr8_2_hddm', 'HDGeant', 'mcsmear', 'gamp2hddm', 'bggen', 'gen_2pi', 'gen_2pi_amp', 'gen_3pi', 'gen_pi0', 'gen_omega_3pi', 'gen_2k']

SConscript(dirs=subdirs, exports='env osname', duplicate=0)

sbms.OptionallyBuild(env, ['genphoton', 'genpi', 'gen_2mu', 'genEtaRegge'])
