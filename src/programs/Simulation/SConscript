
import sbms

Import('*')

<<<<<<< HEAD
subdirs = ['genr8', 'genr8_2_hddm', 'HDGeant', 'mcsmear', 'gamp2hddm', 'bggen', 'gen_2pi', 'gen_2pi_amp', 'gen_3pi', 'gen_pi0', 'gen_omega_3pi', 'gen_2k']
=======
subdirs = ['genr8', 'genr8_2_hddm', 'HDGeant', 'mcsmear', 'bggen', 'gen_2pi', 'gen_2pi_amp', 'gen_2pi_primakoff','gen_3pi', 'gen_pi0', 'gen_omega_3pi', 'nullgen']

>>>>>>> master


sbms.OptionallyBuild(env, ['genphoton', 'genpi', 'gen_2mu', 'genEtaRegge'])
SConscript(dirs=subdirs, exports='env osname', duplicate=0)