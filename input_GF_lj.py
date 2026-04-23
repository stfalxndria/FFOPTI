#THIS ARE ALL INFORMATIONS THAT WILL BE MOVED AND SHOULD BELONG IN A STRUCTURE_INFORMATION.PY FILE
Structure_name = 'GF334'
structure_carfile = './data/GF334_v2/GF334.car'
structure_mdffile = './data/GF334_v2/GF334.mdf'
DFT_ref_path = './data/GF334_v2/DFT_reference'
GMX_ref_path = '/home/uccaset/Scratch/FFOPTI/data/GF334_v2/classical'
charge_file = '/home/uccaset/Scratch/FFOPTI/data/GF334_v2/REPEAT_charges.resp'
groupid='MOF'
atom_types = ['Zn3', 'C_2', 'C_1', 'H_', 'N_R']  # fixed order is important
train = ['lj']
max_iter = 1000
size = [26.494, 26.494, 29.203]
reference_force_field = 'D2'

