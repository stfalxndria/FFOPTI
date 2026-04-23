# collection of all the classes needed for FFOPTI
from all_needed_modules import *

#GETTING INFORMATIONS ABOUT UFF TO USE AS INITIAL INFORMATION
_uff='./data/uff_tu'
with open(_uff,'r') as f:
    lines = f.readlines()
uff = []
for line in lines:
    words = [x.strip() for x in line.split('\t')]
    uff.append(words)
uff_table = pd.DataFrame(uff, columns=['atoms', 'ri', 'phi', 'xi', 'di', 'psi', 'zmm', 'vsp3', 'vsp2', 'chi', 'nc','mass'])
#grouping systems needed for the separation of atoms
#checking whether the atom is a group_6 or not
g6 = ['O','S','Se','Te','Po']


_angleff='data/angle.ff'
with open(_angleff,'r') as f:
    lines = f.readlines()
angle_ff = []
for line in lines:
    words = [x.strip() for x in line.split('\t')]
    angle_ff.append(words)
    
_torsionff='data/torsion.ff'
with open(_torsionff,'r') as f:
    lines = f.readlines()
torsion_ff = []
for line in lines:
    words = [x.strip() for x in line.split('\t')]
    torsion_ff.append(words)
    
_inversionff='data/inversion.ff'
with open(_inversionff,'r') as f:
    lines = f.readlines()
inversion_ff = []
for line in lines:
    words = [x.strip() for x in line.split('\t')]
    inversion_ff.append(words)



class About_Structure:
    '''This file is only able to gain informations through .car and .mdf file, with charge file being optional'''
    class about_atoms:
        def __init__(self, uff_table):
            self.uff_table = uff_table
    
        def get_hybridisation(self, atom):
            """Return hybridisation from uff-style atom name. Safe checks."""
            _hyb = 'sp3'
            if not atom:
                return _hyb
    
            if len(atom) >= 3:
                c = atom[2]
                # c is a char; try to interpret as digit
                if isinstance(c, str) and c.isdigit():
                    val = int(c)
                else:
                    val = None
    
                if val == 1:
                    _hyb = 'sp'
                elif val == 2:
                    _hyb = 'sp2'
                elif atom.lower().find('ar') != -1:
                    # fallback heuristic for aromatic names (if present)
                    _hyb = 'aromatic'
                else:
                    _hyb = 'sp3'
            else:
                _hyb = 'sp3'
            return _hyb
    
        def get_BO(self, atom_A, atom_B):
            """Bond order heuristic using hybridisation and special cases."""
            A_hyb = self.get_hybridisation(atom_A)
            B_hyb = self.get_hybridisation(atom_B)
    
            # special-cases preserved from your original code
            if (atom_A == 'C_R' and atom_B == 'N_R') or (atom_A == 'N_R' and atom_B == 'C_R'):
                return 1.41
            if (atom_A == 'O_3z' and atom_B == 'Si3') or (atom_A == 'Si3' and atom_B == 'O_3z'):
                return 1.44
    
            if A_hyb == 'sp' and B_hyb == 'sp':
                return 3.0
            if A_hyb == 'sp2' and B_hyb == 'sp2':
                return 2.0
            if A_hyb == 'aromatic' and B_hyb == 'aromatic':
                return 1.5
    
            return 1.0
    
        def get_natural_bond(self, atom_A, atom_B):
            """Compute natural bond length using uff parameters table.
            Returns None if required parameters are missing."""
            r_i = None
            x_i = None
            r_j = None
            x_j = None
    
            # find parameters for atom_A
            for n, name in enumerate(self.uff_table['atoms']):
                if atom_A == name:
                    try:
                        r_i = float(self.uff_table['ri'][n])
                        x_i = float(self.uff_table['xi'][n])
                    except Exception:
                        r_i = None
                        x_i = None
                    break
    
            # find parameters for atom_B
            for a, name in enumerate(self.uff_table['atoms']):
                if atom_B == name:
                    try:
                        r_j = float(self.uff_table['ri'][a])
                        x_j = float(self.uff_table['xi'][a])
                    except Exception:
                        r_j = None
                        x_j = None
                    break
    
            # If any parameter missing, return None so caller can handle it
            if r_i is None or r_j is None or x_i is None or x_j is None:
                return None
    
            # bond order
            n_bo = self.get_BO(atom_A, atom_B)
            # bond order correction term
            r_bo_ij = -0.1332 * (r_i + r_j) * math.log(n_bo)
    
            # equilibrium bond distance correction term (avoid division by zero)
            denom = (x_i * r_i + x_j * r_j)
            if denom == 0:
                ren_ij = 0.0
            else:
                ren_ij = r_i * r_j * ((math.sqrt(x_i) - math.sqrt(x_j)) ** 2) / denom
    
            # final natural bond distance
            r_ab = r_i + r_j - ren_ij + r_bo_ij
            return r_ab


    def __init__ (self, name, carfile, mdffile, chargefile, uff_table=uff_table):
        self.name = name
        self.uff_table = uff_table
        self.carfile = carfile
        self.mdffile = mdffile
        self.chargefile = chargefile
        self.atom_tools = self.about_atoms(self.uff_table)

    def read_car(self): #just processes the .car file to make it readable 
        with open(self.carfile, 'r') as f:
            lines = f.readlines()
        car_ = []
        for line in lines:
            word = line.split()
            car_.append(word)
        return car_
    
    def get_charges(self): #processing the charge files
        if (self.chargefile.endswith("pacman.cif")) or (self.chargefile.endswith("pacmof.cif")): #gets us charges that comes from PACMAN and PACMOF
            with open(self.chargefile,'r') as f:
                lines = f.readlines()
            pac = []
            for line in lines:
                word = line.split()
                pac.append(word)

            pac = [x for x in pac if x != []]
            
            charge_list = []
            for i, row in enumerate(pac):
                if len(row)==9:
                    charge_list.append([pac[i][0],pac[i][1]])
            return charge_list

        
        elif self.chargefile.endswith(".resp") : #get charges from CP2K
            with open(self.chargefile,'r') as f:
                lines = f.readlines()
            charge_list = []
            
            for line in lines:
                if len(line.split()) == 4:
                    charge_list.append(line.split()[3])
            return charge_list
        
        else:
            charge_list = None
            print('no charge file found/ readble!')
            return charge_list 
        
    def get_atom_list(self): #obtaining the full list of atoms
        car_ = self.read_car()
        car_table = pd.DataFrame()
        atom_list_unmod = []
        atom_names = []
        ff_names = []
        atom_elements = []
        target1 = '(P1)'
        target2 = 'end'
        inside_target = False  
        for i, line in enumerate(car_):
            if target1 in line:
                line = car_[i+1]
                inside_target = True
            elif target2 in line:
                inside_target = False

            if (inside_target == True) and (line not in atom_list_unmod):
                atom_list_unmod.append(line)
                atom_names.append(line[0])
                ff_names.append(line[6])
                atom_elements.append(line[7])
                
                
        car_table['index'] = np.arange(1,len(atom_names)+1,1)
        car_table['labels']= atom_names
        car_table['elements']=atom_elements
        car_table['uff names']=ff_names         
        return(car_table)
    
    def read_mdf(self): #reading the mdf file
        with open(self.mdffile, 'r') as f:
            lines = f.readlines()
        mdf_unmod = []
        mdf_ = []
        for line in lines:
            word = line.split()
            mdf_unmod.append(word)
        
        for i, line in enumerate(mdf_unmod):
            start_ = 21
            target = '!'
            if target in line:
                end_ = i-1
                mdf_ = mdf_unmod[start_:end_]
        mdf_ = [[s.replace("XXXX_1:", "") for s in inner_list] for inner_list in mdf_]
        return(mdf_)
    
    def raw_full_bond_list(self): #obtaining the full but raw list of connectivitiy
        mdf_ = self.read_mdf()
        all_bond_list = pd.DataFrame()
        atoms_t = []
        bonded = []
        for row in mdf_:
            x = 12
            y = len(row)+1
            atoms_t.append(row[0])
            bonded.append(row[x:y])
        all_bond_list['atom']=atoms_t
        all_bond_list['bonded'] = bonded

        all_bond_list['bonded'] = [[atom.split('%')[0] for atom in row] for row in all_bond_list['bonded']]
        all_bond_list['bonded'] = [[atom.split('/')[0] for atom in row] for row in all_bond_list['bonded']]

        return all_bond_list
    
    def unique_atoms(self):
        u_atoms = []
        car_tab = self.get_atom_list()
        for i, line in enumerate(car_tab['uff names']):
            atom = car_tab['uff names'][i]
            elements = car_tab['elements'][i]
            if atom not in [item[0] for item in u_atoms]:
                u_atoms.append([atom, elements])
        
        return u_atoms


    def atom_types(self):
        return [line[0] for line in self.unique_atoms()]
    
    def get_full_bonds(self):
        mdf_ = self.read_mdf()
        bond_connectivity = []
        for i, line in enumerate(mdf_):
            atom_A = line[0].split(':')[-1]
            x = 12
            y = len(line)
            for b in range(x,y):
                if ('%' in line[b]) and ('/' in line[b]):
                    temp = line[b].split('%')[0]
                    atom_B = temp.split('/')[0]
                if '%' in line[b]:
                    atom_B = line[b].split('%')[0]
                elif '/' in line[b]:
                    atom_B = line[b].split('/')[0]
                else:
                    atom_B = line[b]
                if ([atom_A,atom_B] not in bond_connectivity) and ([atom_B,atom_A] not in bond_connectivity):
                    bond_connectivity.append([atom_A,atom_B])
        return bond_connectivity
    
    def unique_bonds(self):
        u_bonds = []
        bond_con = self.get_full_bonds()
        car_table = self.get_atom_list()
        for b, ll in enumerate(bond_con):
            uff_A = None
            uff_B = None
            atom_A = ll[0]
            atom_B = ll[1]
            for a, line in enumerate(car_table['labels']):
                if atom_A == car_table['labels'][a]:
                    uff_A = car_table['uff names'][a]
            for c, l3 in enumerate(car_table['labels']):
                if atom_B == car_table['labels'][c]:
                    uff_B = car_table['uff names'][c]
            
            if ([uff_A,uff_B] not in u_bonds) and ([uff_B, uff_A] not in u_bonds):
                u_bonds.append([uff_A,uff_B])
        return u_bonds

    def get_full_angles(self):
        angle_connectivity=[]
        bond_connectivity = self.get_full_bonds()
        car_table = self.get_atom_list()
        for i, line in enumerate(bond_connectivity):
            atom_A, atom_B = line[0], line[1]
            for a, line in enumerate(car_table['labels']):
                atom_C = car_table['labels'][a]
                angle1 = None
                angle2 = None
                if atom_B != atom_C: #checking for bond A-C, C must not be B
                    if ([atom_A, atom_C] in bond_connectivity) or ([atom_C,atom_A] in bond_connectivity):
                        angle1 = [atom_C, atom_A, atom_B]
                if atom_A != atom_C:#checking for bond B-C, C must not be A
                    if ([atom_B, atom_C] in bond_connectivity) or ([atom_C, atom_B] in bond_connectivity):
                        angle2 = [atom_A, atom_B, atom_C]
                if ([atom_C, atom_A, atom_B] not in angle_connectivity) and ([atom_B, atom_A, atom_C] not in angle_connectivity) and (angle1 != None):
                    angle_connectivity.append(angle1)
                if ([atom_A, atom_B, atom_C] not in angle_connectivity) and ([atom_C, atom_B, atom_A] not in angle_connectivity) and (angle2 != None):
                    angle_connectivity.append(angle2)

        return angle_connectivity
    
    def unique_angles(self):
        u_angles=[]
        angle_connectivity = self.get_full_angles()
        car_table = self.get_atom_list()

        for i, line in enumerate(angle_connectivity):
            uff_A = None
            uff_B = None
            uff_C = None
            atom_A, atom_B, atom_C = line[0], line[1], line[2]
            for a, l1 in enumerate(car_table['labels']):
                if atom_A == car_table['labels'][a]:
                    uff_A = car_table['uff names'][a]
            for b, l2 in enumerate(car_table['labels']):
                if atom_B == car_table['labels'][b]:
                    uff_B = car_table['uff names'][b]
            for c, l3 in enumerate(car_table['labels']):
                if atom_C == car_table['labels'][c]:
                    uff_C = car_table['uff names'][c]
            if ([uff_A,uff_B,uff_C] not in u_angles) and ([uff_C,uff_B, uff_A] not in u_angles):
                u_angles.append([uff_A,uff_B,uff_C])
        
        return u_angles
    
    def get_full_torsions(self):
        torsion_connectivity=[]
        angle_connectivity = self.get_full_angles()
        bond_connectivity = self.get_full_bonds()
        car_table = self.get_atom_list()
        for i, line in enumerate(angle_connectivity):
            atom_A, atom_B, atom_C = line[0], line[1], line[2]
            for a, line in enumerate(car_table['labels']):
                atom_D = car_table['labels'][a]
                
                tors1 = None
                tors2 = None
                
                if atom_B != atom_D:
                    if ([atom_A, atom_D] in bond_connectivity) or ([atom_D, atom_A] in bond_connectivity):
                        tors1 = [atom_D, atom_A, atom_B, atom_C]
                if atom_B != atom_D:
                    if ([atom_C, atom_D] in bond_connectivity) or ([atom_D, atom_C] in bond_connectivity):
                        tors2 = [atom_A, atom_B, atom_C, atom_D]
                        
                if ([atom_D, atom_A, atom_B, atom_C] not in torsion_connectivity) and ([atom_C, atom_B, atom_A, atom_D] not in torsion_connectivity) and (tors1 != None):
                    torsion_connectivity.append(tors1)
                    
                if ([atom_A, atom_B, atom_C, atom_D] not in torsion_connectivity) and ([atom_D, atom_C, atom_B, atom_A] not in torsion_connectivity) and (tors2 != None):
                    torsion_connectivity.append(tors2)
        
        return torsion_connectivity
    
    def unique_torsions(self):
        u_torsions = []
        car_table = self.get_atom_list()
        torsion_connectivity = self.get_full_torsions()
        for i, line in enumerate(torsion_connectivity):
            uff_A = None
            uff_B = None
            uff_C = None
            uff_D = None
            atom_A, atom_B, atom_C, atom_D = line[0], line[1], line[2], line[3]
            for a, l1 in enumerate(car_table['labels']):
                if atom_A == car_table['labels'][a]:
                    uff_A = car_table['uff names'][a]
            for b, l2 in enumerate(car_table['labels']):
                if atom_B == car_table['labels'][b]:
                    uff_B = car_table['uff names'][b]
            for c, l3 in enumerate(car_table['labels']):
                if atom_C == car_table['labels'][c]:
                    uff_C = car_table['uff names'][c]
            for d, l3 in enumerate(car_table['labels']):
                if atom_D == car_table['labels'][d]:
                    uff_D = car_table['uff names'][d]
            
            if ([uff_A,uff_B,uff_C, uff_D] not in u_torsions) and ([uff_D, uff_C,uff_B, uff_A] not in u_torsions):
                u_torsions.append([uff_A,uff_B,uff_C,uff_D])    
                    
        return u_torsions
    
    def get_full_inversions(self):
        from itertools import combinations as combi
        inversion_connectivity = []
        all_bond_list = self.raw_full_bond_list()
        for i, atom in enumerate(all_bond_list['atom']):
            atom_j = atom
            bonded_ = all_bond_list['bonded'][i]
            inversion_ = None
            if len(bonded_) >= 3:
                y = list(combi(bonded_, 3))
                inversion_ = [(atom_j,) + combo for combo in y]
                
                for x in inversion_:
                    if x not in inversion_connectivity:
                        inversion_connectivity.append(x)
        return inversion_connectivity
    
    def unique_inversions(self):
        inversion_connectivity = self.get_full_inversions()
        car_table = self.get_atom_list()
        u_inversions = []

        for i, line in enumerate(inversion_connectivity):
            uff_A = None
            uff_B = None
            uff_C = None
            uff_D = None
            atom_A, atom_B, atom_C, atom_D = line[0], line[1], line[2], line[3]
            for a, l1 in enumerate(car_table['labels']):
                if atom_A == car_table['labels'][a]:
                    uff_A = car_table['uff names'][a]
            for b, l2 in enumerate(car_table['labels']):
                if atom_B == car_table['labels'][b]:
                    uff_B = car_table['uff names'][b]
            for c, l3 in enumerate(car_table['labels']):
                if atom_C == car_table['labels'][c]:
                    uff_C = car_table['uff names'][c]
            for d, l3 in enumerate(car_table['labels']):
                if atom_D == car_table['labels'][d]:
                    uff_D = car_table['uff names'][d]
            
            if ([uff_A,uff_B,uff_C, uff_D] not in u_inversions) and ([uff_D, uff_C,uff_B, uff_A] not in u_inversions):
                u_inversions.append([uff_A,uff_B,uff_C,uff_D]) 
        return u_inversions

    def get_lj(self):
        masses = []
        for i in range(len(self.unique_atoms())):
            atom_uff = self.unique_atoms()[i][0]
            _atom = self.unique_atoms()[i][1]
            atom_mass = Formula(_atom).mass
            masses.append([atom_uff,atom_mass])

        atom_types = []
        #Lennard Jones Potential in GROMACS UNITS
        for n in masses:
            _atom=n[0]
            _mass=round(n[1],3)
            for i, line in enumerate(self.uff_table['atoms']):
                if _atom == self.uff_table['atoms'][i]:
                    _energy=round(float(self.uff_table['di'][i])*4.184,4)
                    _radius=round((float(self.uff_table['xi'][i]))*0.89089871814*0.1,6)
                    atom_types.append([_atom,'1',_mass,'0.000','A',_radius,_energy])
            
        return atom_types
    
    def get_full_hybridisation(self):
        label_list = []
        for line in self.unique_atoms():
            label_list.append(line[0])

        atom_hybridisation=pd.DataFrame()
        hyb_list = []
        for i in range(len(label_list)):
            hyb_list.append(self.atom_tools.get_hybridisation(label_list[i]))

        atom_hybridisation['atoms']=label_list
        atom_hybridisation['hybridisation']=hyb_list

        return atom_hybridisation
    
    def get_bond_BO(self):

        #now assigning bond orders for each pairs
        BO_list=[]
        for i in range(len(self.unique_bonds())):
            atom_A = self.unique_bonds()[i][0]
            atom_B = self.unique_bonds()[i][1]
            BO_list.append([atom_A,atom_B,self.atom_tools.get_BO(atom_A, atom_B)])
        return BO_list
    
    def get_bond_distance(self):
        #DISTANCE: calculating the r_ij value 
            #r_ij = r_i + r_j + r_BO - r_EN
        r_ij=[] 
        for i in range(len(self.unique_bonds())):
            atom_A = self.unique_bonds()[i][0]
            atom_B = self.unique_bonds()[i][1]
            _r = self.atom_tools.get_natural_bond(atom_A, atom_B)
            
            r_ij.append([atom_A, atom_B, round(_r,4)])

        return r_ij    #returning the results in GROMACS units
    
    def get_bond_energy(self):
        k_ij=[] #kcal/mol/A^2
        for i in range(len(self.unique_bonds())):
            atom_A, atom_B = self.unique_bonds()[i]
            z_i = None
            z_j = None
            rab=[]
            for n, line in enumerate(self.uff_table['atoms']):
                if atom_A == self.uff_table['atoms'][n]:
                    z_i = float(self.uff_table['zmm'][n])
            for a, ll in enumerate(self.uff_table['atoms']):
                if atom_B == self.uff_table['atoms'][a]:
                    z_j = float(self.uff_table['zmm'][a]) 

            for b, lll in enumerate(self.get_bond_distance()):
                if atom_A in lll[0] and atom_B in lll[1]:
                    rab=lll[2]
                    k_ab=664.12*((z_i*z_j)/(rab))
                    k_ij.append([atom_A, atom_B,k_ab])
        return k_ij
    
    def get_bond_ff(self):
        #This file converts the informations to gromacs units for the bond section
        k_ij_G = []
        bond_types = []
        k_ij = self.get_bond_energy()
        r_ij = self.get_bond_distance()
        for i in range(len(k_ij)):
            kab=k_ij[i][2]*418.4    #kJ/mol/nm^2
            kab=round(kab,2)
            k_ij_G.append(kab)
            
        for i, line in enumerate(k_ij):
            atom_A = line[0]
            atom_B = line[1]
            _energy = k_ij_G[i]
            _radius = r_ij[i][2]*0.1
            
            bond_types.append([atom_A, atom_B, '1', _radius, _energy])
        
        return bond_types
    
    def get_natural_angles(self):
        #obtaining the middle atom natural angle
        theta=[] #in degrees
        for i in range(len(self.unique_angles())):
            _atom=self.unique_angles()[i][1]
            for n, row in enumerate(self.uff_table['atoms']):
                if _atom == self.uff_table['atoms'][n]:
                    theta.append([_atom,self.uff_table['phi'][n]])

        return theta
    
    def get_angle_force_constant(self):
        #calculate the force constant K_ijk 
        K_ijk=[] #kcal/mol/rad*2
        theta = self.get_natural_angles()

        for i, line in enumerate(self.unique_angles()):
            atom_A, atom_B, atom_C =line[0:3]  #i, j, k
            z_i = None
            z_k = None
            angle = None
            
            #calculating the natural bonds

            _r_ij=self.atom_tools.get_natural_bond(atom_A, atom_B)  #atom 1 and 2
            _r_jk=self.atom_tools.get_natural_bond(atom_B, atom_C)  #atom 2 and 3
            for n, row in enumerate(theta):
                if atom_B == theta[n][0]:
                    angle=float(theta[n][1])
                    theta_rad = math.radians(angle)
                    
            for b, ll in enumerate(self.uff_table['atoms']): 
                if atom_A == self.uff_table['atoms'][b]:
                    z_i = float(self.uff_table['zmm'][b])
                    
            for d, llll in enumerate(self.uff_table['atoms']):
                if atom_C == self.uff_table['atoms'][d]:
                    z_k = float(self.uff_table['zmm'][d])
                    _r_ik=math.sqrt((_r_ij**2)*(_r_jk**2)-2*(_r_ij*_r_jk*math.cos(theta_rad)))  #atom 1 and 3
                    
                    
                    #K_ijk equation
                    _beta = 664.12/(_r_ij*_r_jk)
                    _back = 3*_r_ij*_r_jk*(1-((math.cos(theta_rad))**2))-(((_r_ik)**2)*math.cos(theta_rad))
                    _kijk= _beta*((z_i*z_k)/(_r_ik)**5)*(_r_ij*_r_jk)*_back
                    K_ijk.append([atom_A, atom_B, atom_C,_kijk])
        return K_ijk
    
    def get_angle_ff(self, force_field_mixing=False, angle_ff='./data/angle.ff'):
        angle_types = []
        theta = self.get_natural_angles()
        K_ijk = self.get_angle_force_constant()
    
        # If force_field_mixing=True, parse angle_ff file, else empty list
        if force_field_mixing:
            with open(angle_ff, 'r') as f:
                angle_ff_lines = [line.strip().split() for line in f if line.strip()]
        else:
            angle_ff_lines = []
    
        for line in self.unique_angles():
            atom_A, atom_B, atom_C = line[0:3]
    
            K = None
            angle = None
    
            # Find angle from theta
            for row in theta:
                if atom_B == row[0]:
                    angle = float(row[1])
                    break  # found, stop
    
            # Find force constant from K_ijk
            for sub in K_ijk:
                if atom_B == sub[1]:
                    K = round(4.184 * float(sub[3]), 2)  # kcal/mol to kJ/mol
                    break
    
            # Handle force_field_mixing == True case
            if force_field_mixing:
                for l2 in angle_ff_lines:
                    if [atom_A, atom_B, atom_C] == l2[:3]:
                        deg_harm = float(l2[3])
                        k_harm = float(l2[4])
                        angle_types.append([atom_A, atom_B, atom_C, '1', deg_harm, k_harm])
                        break
                # Move to next angle after mixing handled
                continue
    
            # For force_field_mixing == False, do your angle logic:
            if (len(atom_B) >= 3) and (atom_B[2] in ['1', '2', '4', '6']):
                if atom_B[2] == '1':
                    k_harm = K
                    deg_harm = 180
                elif atom_B[2] == '2':
                    k_harm = 4 * K / 3
                    deg_harm = 120
                else:  # '4' or '6'
                    k_harm = 4 * K / 3
                    deg_harm = 120
            else:
                if angle == 180.0:
                    k_harm = K
                    deg_harm = 180
                elif angle == 120.0:
                    k_harm = 4 * K / 3
                    deg_harm = 120
                else:
                    k_harm = K
                    deg_harm = angle
    
            angle_types.append([atom_A, atom_B, atom_C, '1', deg_harm, k_harm])
    
        return angle_types

    
    def get_torsion_ff(self, force_field_mixing=False, torsion_ff=None):
        dihedral_types = []
        g6 = ['S', 'Se', 'Te']  # define group 6 atoms if not defined elsewhere
    
        for line in self.unique_torsions():
            atom_A, atom_B, atom_C, atom_D = line[0:4]
            B_hyb = self.atom_tools.get_hybridisation(atom_B)
            C_hyb = self.atom_tools.get_hybridisation(atom_C)
    
            V_tot = None
            nval = None
            _phi = None
            k_gro = None
            deg_gro = None
    
            if force_field_mixing and torsion_ff is not None:
                # Loop over torsion force field entries
                found_match = False
                for gamma in torsion_ff:
                    # Assuming gamma is like [atom1, atom2, atom3, atom4, deg, k, n]
                    if ([atom_A, atom_B, atom_C, atom_D] == gamma[:4]) or ([atom_D, atom_C, atom_B, atom_A] == gamma[:4]):
                        deg_gro = gamma[4]
                        k_gro = gamma[5]
                        nval = gamma[6]
                        dihedral_types.append([atom_A, atom_B, atom_C, atom_D, '1', deg_gro, k_gro, nval])
                        found_match = True
                        break
                    elif (['*', atom_B, atom_C, '*'] == gamma[:4]) or (['*', atom_C, atom_B, '*'] == gamma[:4]):
                        deg_gro = gamma[4]
                        k_gro = gamma[5]
                        nval = gamma[6]
                        dihedral_types.append([atom_A, atom_B, atom_C, atom_D, '1', deg_gro, k_gro, nval])
                        found_match = True
                        break
                if found_match:
                    continue  # skip further calculation since FF data found
    
            # Calculate torsion parameters when no force field mixing or no match found
            if B_hyb in ['sp3', 'sp2'] and C_hyb in ['sp3', 'sp2']:
    
                if (B_hyb == 'sp3') and (C_hyb == 'sp3'):
                    # sp3-sp3 pairs
                    if any(elem in atom_B for elem in g6) and any(elem in atom_C for elem in g6):
                        if atom_B.startswith('O') or atom_C.startswith('O'):
                            V_tot = math.sqrt(2.0 * 6.8)  # kcal/mol
                            _phi = 0
                            nval = 2
                        else:
                            V_tot = 2  # kcal/mol
                            _phi = 0
                            nval = 2
                    else:
                        # Other sp3-sp3 bonds
                        nval = 3
                        _phi = 180
                        vj = None
                        vk = None
                        for a, atom_name in enumerate(self.uff_table['atoms']):
                            if atom_B == atom_name:
                                vj = float(self.uff_table['vsp3'][a])
                            if atom_C == atom_name:
                                vk = float(self.uff_table['vsp3'][a])
                        if vj is not None and vk is not None:
                            V_tot = math.sqrt(vj * vk)  # kcal/mol
    
                elif (B_hyb == 'sp2') and (C_hyb == 'sp2'):
                    # sp2-sp2 pairs
                    vj = None
                    vk = None
                    nval = 2
                    _phi = 180
                    for c, atom_name in enumerate(self.uff_table['atoms']):
                        if atom_B == atom_name:
                            vj = float(self.uff_table['vsp2'][c])
                        if atom_C == atom_name:
                            vk = float(self.uff_table['vsp2'][c])
                    if vj is not None and vk is not None:
                        rbo_jk = self.atom_tools.get_BO(atom_B, atom_C)
                        V_tot = 5 * math.sqrt(vj * vk) * (1 + 4.18 * math.log(rbo_jk))  # kcal/mol
    
                elif (B_hyb == 'sp2' and C_hyb == 'sp3') or (B_hyb == 'sp3' and C_hyb == 'sp2'):
                    # sp2-sp3 pairs
                    A_hyb = self.atom_tools.get_hybridisation(atom_A)
                    D_hyb = self.atom_tools.get_hybridisation(atom_D)
    
                    if (A_hyb == 'sp2' and B_hyb == 'sp2') or (C_hyb == 'sp2' and D_hyb == 'sp2'):
                        V_tot = 2  # kcal/mol
                        nval = 3
                        _phi = 180
                    else:
                        nval = 6
                        _phi = 0
                        V_tot = 1  # kcal/mol
    
                if V_tot is not None and nval is not None and _phi is not None:
                    k_gro = round(0.5 * V_tot * 4.184, 4)  # convert kcal/mol to kJ/mol
                    deg_gro = nval * _phi - 180
    
                    dihedral_types.append([atom_A, atom_B, atom_C, atom_D, '1', deg_gro, k_gro, nval])
    
        return dihedral_types
         

    def get_inversion_ff(self, force_field_mixing=False, inversion_ff='./data/inversion.ff'):
        inversion_types = []
        inv_target = ['O_2','O_R', 'N_R', 'N_2']
        for a, line in enumerate(self.unique_inversions()):
            _c0 = 0
            _c1 = 0
            _c2 = 0
            _k = 0  #kcal/mol
            k_gro = None #kj/mol
            A_gro = None
            
            atom_A=line[0]
            atom_B=line[1]
            atom_C=line[2]
            atom_D=line[3]
            
            if force_field_mixing==True:
                for b, l2 in enumerate(inversion_ff):
                    if [atom_A,atom_B,atom_C,atom_D] in l2:
                        A_gro = l2[4]
                        k_gro = l2[5]
                        inversion_types.append([atom_A,atom_B, atom_C, atom_D, '4', k_gro, A_gro, '#inversion'])
                        break
            else:
                continue
            
            if (atom_A == 'C_2') or (atom_A == 'C_R'):
                _c0 = 1
                _c1 = -1
                _c2 = 0
                if 'O_2' in line:
                    _k = 50 
                else: 
                    _k = 6
                k_gro = 4*_c2*(1-(_c1/(4*_c2))**2)*k*4.184 #kj/mol
                A_gro = math.pi- math.acos(_c1/(4*c_c2))
                    
            elif atom_A in inv_target:
                _c0 = 1
                _c1 = -1
                _c2 = 0
                _k = 6
                
                k_gro = 4*_c2*(1-(_c1/(4*_c2))**2)*k*4.184 #kj/mol
                A_gro = math.pi- math.acos(_c1/(4*c_c2))
                
                
            else:
                k_gro = 0 #kj/mol
                A_gro = 0
                
            inversion_types.append([atom_A,atom_B, atom_C, atom_D, '4', k_gro, A_gro, '#inversion'])
    
    def write_ff_itp(self,print_directory='./'):
        ff_list = []
        ff_list.append(['[ defaults ]'])
        ff_list.append(['; nbfunc','comb-rule','gen-pairs','fudgeLJ','fudgeQQ'])
        ff_list.append(['1','2','No','0.5','0.5'])
        ff_list.append(['[ atomtypes ]'])
        ff_list.append([';name','bond_type','mass','charge','ptype','sigma','epsilon'])
        print('starting LJ')
        ff_list += self.get_lj()
        ff_list.append([''])
        ff_list.append(['[ bondtypes ]'])
        print('starting bond')
        ff_list += self.get_bond_ff()
        ff_list.append([''])
        ff_list.append(['[ angletypes ]'])
        print('starting angle')
        ff_list += self.get_angle_ff()
        ff_list.append([''])
        
        ff_list.append(['[ dihedraltypes ]'])
        print('starting dihedral')
        ff_list += self.get_torsion_ff()
        #print('starting inversion')
        #ff_list += self.get_inversion_ff()    
        ff_name=f'{print_directory}/force_field.itp'
        with open(ff_name, 'w') as file:
            for row in ff_list:
                row_str = '\t'.join(map(str, row)) 
                file.write(row_str + '\n')
                
        print("force_field.itp file printed")

    def write_top_file(self,groupid='MOF',print_directory='./'):
        top_list = []
        top_list.append(['#include "force_field.itp"'])
        top_list.append([f'#include "{groupid}.itp" '])
        top_list.append([''])
        top_list.append(['[ system ]'])
        top_list.append([f'{groupid}','',''])
        top_list.append([''])
        top_list.append(['[ molecules ]'])
        top_list.append([groupid,1])


        top_name=f'{print_directory}/topol.top'
        with open(top_name, 'w') as file:
            for row in top_list:
                row_str = '\t'.join(map(str, row)) 
                file.write(row_str + '\n')
                
        print("topol.top file printed")        

    def write_gmx_connectivity_itp(self, groupid='MOF',print_directory='./'):
        atom_con = []
        car_table = self.get_atom_list()
        charge_list=self.get_charges()
        if charge_list is None:
            charge_list = [(0, 0)] * len(car_table)

        
        for i, row in enumerate(car_table['index']):
            atom_charge = charge_list[i][1]
            atom_con.append([row,car_table['uff names'][i],'1',groupid,car_table['uff names'][i],' ',atom_charge])

        bond_con = []
        bond_connectivity = self.get_full_bonds()
        for i, row in enumerate(bond_connectivity):
            atom_A = row[0]
            atom_B = row[1]
            index_A = None
            index_B = None
            for a, l1 in enumerate(car_table['labels']):
                if atom_A == l1:
                    index_A = car_table['index'][a]
            for b, l2 in enumerate(car_table['labels']):
                if atom_B == l2:
                    index_B = car_table['index'][b]
            
            bond_con.append([index_A, index_B, '1'])
            
        angle_con = []
        angle_connectivity = self.get_full_angles()
        for i, row in enumerate(angle_connectivity):
            atom_A = row[0]
            atom_B = row[1]
            atom_C = row[2]
            
            index_A = None
            index_B = None
            index_C = None
            for a, l1 in enumerate(car_table['labels']):
                if atom_A == l1:
                    index_A = car_table['index'][a]
            for b, l2 in enumerate(car_table['labels']):
                if atom_B == l2:
                    index_B = car_table['index'][b]
            for c, l3 in enumerate(car_table['labels']):
                if atom_C == l3:
                    index_C = car_table['index'][c]
            
            angle_con.append([index_A, index_B,index_C, '1'])
            
        torsion_con = []
        torsion_connectivity = self.get_full_torsions()
        for i, row in enumerate(torsion_connectivity):
            atom_A = row[0]
            atom_B = row[1]
            atom_C = row[2]
            atom_D = row[3]
            
            index_A = None
            index_B = None
            index_C = None
            index_D = None
            for a, l1 in enumerate(car_table['labels']):
                if atom_A == l1:
                    index_A = car_table['index'][a]
            for b, l2 in enumerate(car_table['labels']):
                if atom_B == l2:
                    index_B = car_table['index'][b]
            for c, l3 in enumerate(car_table['labels']):
                if atom_C == l3:
                    index_C = car_table['index'][c]
            for d, l3 in enumerate(car_table['labels']):
                if atom_D == l3:
                    index_D = car_table['index'][d]
            
            torsion_con.append([index_A, index_B,index_C,index_D,'1'])
            
        #inversion_con = []
        #inversion_connectivity = self.get_full_inversions()
        #for i, row in enumerate(inversion_connectivity):
            #atom_A = row[0]
            #atom_B = row[1]
            #atom_C = row[2]
            #atom_D = row[3]
            
            #index_A = None
            #index_B = None
            #index_C = None
            #index_D = None
            #for a, l1 in enumerate(car_table['labels']):
                #if atom_A == l1:
                    #index_A = car_table['index'][a]
            #for b, l2 in enumerate(car_table['labels']):
                #if atom_B == l2:
                    #index_B = car_table['index'][b]
            #for c, l3 in enumerate(car_table['labels']):
                #if atom_C == l3:
                    #index_C = car_table['index'][c]
            #for d, l3 in enumerate(car_table['labels']):
                #if atom_D == l3:
                    #index_D = car_table['index'][d]
            
            #inversion_con.append([index_A, index_B,index_C,index_D,'4'])
                    
        con_list = []
        con_list.append(['[ moleculetype ]'])
        con_list.append([';','','Name','nrexcl'])
        con_list.append([f'{groupid}','','3'])
        con_list.append([])
        con_list.append(['[ atoms ]'])
        con_list += atom_con
        con_list.append([])
        con_list.append(['[ bonds ]'])
        con_list += bond_con
        con_list.append([])
        con_list.append(['[ angles ]'])
        con_list += angle_con
        con_list.append([])
        con_list.append(['[ dihedrals ]'])
        con_list += torsion_con
        #con_list.append(['[ dihedrals ]'])
        #con_list += inversion_con

        itpname = f'{print_directory}/{groupid}.itp'
        with open(itpname, 'w') as file:
            for row in con_list:
                row_str = '\t'.join(map(str, row)) 
                file.write(row_str + '\n')

        print(f"Data written to {itpname}")
            

   

class DFT_reference:
    def __init__(self, Structure_obj, name, DFT_path, structure_file):
        self.Structure = Structure_obj
        self.name = name 
        self.DFT_path = DFT_path
        self.force_path = f"{DFT_path}/md-frc-1.xyz"
        self.pdb_path = f"{DFT_path}/md-pos-1.pdb"
        self.xyz_path = f"{DFT_path}/md-pos-1.xyz"
        self.energy_path = f"{DFT_path}/md-1.ener"
        self.outfile_path =f"{DFT_path}/output.out"
        self.cif_file = f"{DFT_path}/{structure_file}"
        self.aseobj = read(self.cif_file)
        self.atom_nums = len(self.aseobj)
    
#READING ENERGIES
    def get_energy(self):
        _name = self.energy_path
        with open(_name,'r') as f:
            lines = f.readlines()
        _ener = []
        for line in lines:
            words = [x.strip() for x in line.split()]
            _ener.append(words)
        
        ener = pd.DataFrame()
        step_list = []
        time_fs = []
        energy_list = []
        temperature = []
        
        for a, l1 in enumerate(_ener):
            a = a+1
            if a < len(_ener):
                step_list.append(float(_ener[a][0]))
                time_fs.append(float(_ener[a][1]))
                temperature.append(float(_ener[a][3]))
                energy_list.append(float(_ener[a][4])) 
    
        ener['steps'] = step_list
        ener['time (fs)'] = time_fs
        ener['temp (k)'] = temperature
        ener['energy (HF)'] = energy_list
            
        return ener
    
    def get_number_of_frames(self):
        return(len(self.get_energy()))
    
    def get_average_charges(self):
        charges_per_atom = []
        n_atoms = None
        n_frames = 0

        with open(self.outfile_path) as f:   # change filename if needed
            inside_block = False
            frame_charges = []
            for line in f:
                if "Hirshfeld Charges" in line:
                    inside_block = True
                    frame_charges = []
                    n_frames+=1
                    continue

                if inside_block and "Total Charge" in line:
                    inside_block = False
                
                    # Initialize storage on first frame
                    if n_atoms is None:
                        n_atoms = len(frame_charges)
                        charges_per_atom = [[] for _ in range(n_atoms)]

                    # Store charges
                    for i, q in enumerate(frame_charges):
                        charges_per_atom[i].append(q)

                    continue
                # Parse data lines
                if inside_block:
                    parts = line.split()
                    if len(parts) >= 8 and parts[0].isdigit():
                        net_charge = float(parts[-1])
                        frame_charges.append(net_charge)

        # Convert to numpy for statistics
        charges_per_atom = np.array(charges_per_atom)
        avg_charges = charges_per_atom.mean(axis=1)
        if n_atoms is None:
            raise ValueError("No Hirshfeld charges found in the file.")
        return avg_charges
    
    def get_type_charges(self):
        ab = self.DFT_trajectories_ASE_obj()[0]
        symbols = self.Structure.get_atom_list()['uff names']

        if self.Structure.get_charges() == None:
            charges = self.get_average_charges()
    
            # group charges by element
            grouped = defaultdict(list)
    
            for elem, q in zip(symbols, charges):
                grouped[elem].append(q)
    
            # compute averages
            avg_charges = {elem: np.mean(qs) for elem, qs in grouped.items()}
    
            avg_charge_list = [float(f"{avg_charges[elem]:}") for elem in symbols]
    
            return avg_charge_list
        else:
            charges = np.array(self.Structure.get_charges())
            charges = [float(c) for c in charges]
                        # group charges by element
            grouped = defaultdict(list)
    
            for elem, q in zip(symbols, charges):
                grouped[elem].append(q)
    
            # compute averages
            avg_charges = {elem: np.mean(qs) for elem, qs in grouped.items()}
    
            avg_charge_list = [float(f"{avg_charges[elem]:}") for elem in symbols]
    
            return avg_charge_list

    #READING FORCES OF EACH FRAMES AND ORGANISED INTO COLUMNS
    
    def get_force_xyz(self):
        _name = self.force_path
        atom_number = self.atom_nums
        with open(_name,'r') as f:
            lines = f.readlines()
        _pos = []
        for line in lines:
            words = [x.strip() for x in line.split()]
            _pos.append(words)
        
        pos = pd.DataFrame()
        step_list = []
        force_list = []
        target1 = 'i'
        for a, l1 in enumerate(_pos):
            if target1 in l1:
                num = a+1
                _forces = _pos[num:num+atom_number]
                _step = l1[2]
                step_list.append(_step)
                force_list.append(_forces)
                
        step_list = [s.replace(",", "") for s in step_list]   
        
        pos['steps']=step_list
        pos['force']=force_list    
        

        return pos

    #READING the PDB trajectory 
    def get_pos_pdb(self):
        _name = self.pdb_path
        atom_number = self.atom_nums
        with open(_name,'r') as f:
            lines = f.readlines()
        _pos = []
        for line in lines:
            words = [x.strip() for x in line.split()]
            _pos.append(words)
        
        pos = pd.DataFrame()
        step_list = []
        coord_list = []
        dimensions = []
        angles = []
        target1 = 'Step'
        for a, l1 in enumerate(_pos):
            if target1 in l1:
                start = a+2
                num_dim = a+1
                _coords = _pos[start:start+atom_number]
                #removing ATOMS, INDEX and VECTORS 
                _coords = [s[2:-3] for s in _coords]
                _step = l1[2].replace(',','')
                _dim = _pos[num_dim][1:4]
                _ang = _pos[num_dim][4:]
                
                step_list.append(_step)
                #removing the commas from the number
                #step_list = [step.replace(',','') for step in step_list]
                coord_list.append(_coords)
                dimensions.append(_dim)
                angles.append(_ang)
        
        
        pos['steps']=step_list
        pos['coordinates']=coord_list   
        pos['dimensions']=dimensions
        pos['angles']=angles
        
        return pos

    #READING XYZ TRAJECTORY
    def get_pos_xyz(self):
        _name = self.xyz_path
        atom_number = self.atom_nums
        with open(_name,'r') as f:
            lines = f.readlines()
        _pos = []
        for line in lines:
            words = [x.strip() for x in line.split()]
            _pos.append(words)
        
        pos = pd.DataFrame()
        step_list = []
        coord_list = []
        target1 = 'i'
        for a, l1 in enumerate(_pos):
            if target1 in l1:
                num = a+1
                _coords = _pos[num:num+atom_number]
                _step = l1[2]
                step_list.append(_step)
                coord_list.append(_coords)
                
        step_list = [s.replace(",", "") for s in step_list]   
        
        pos['steps']=step_list
        pos['coordinates']=coord_list   
        
        return pos
    
    def get_positions(self):
        try:
            self.get_pos_pdb()
        except FileNotFoundError:
            self.get_pos_xyz()
            
    def DFT_trajectories_ASE_obj(self):
        obj_trj_ = read(self.pdb_path, index=':')
        return obj_trj_
    
#    def all_DFT_data(self):
#        dft_data = pd.DataFrame()
#        dft_data['Energy'] = self.get_energy()
#        dft_data['Force'] = self.get_force_xyz()
        
#        dft_data['coordinates'] = self.get_positions()
#        dft_data['ASE trajectories'] = self.DFT_trajectories_ASE_obj()
        
#        return dft_data
    def force_unit_convert(self, force_list, CONVERSION):
        """
        force_list: list of [element, fx, fy, fz]
        returns: same structure with scaled forces
        """
        converted = []
        for atom in force_list:
            element = atom[0]
            forces = np.array(atom[1:], dtype=float) * CONVERSION
            converted.append([element, *forces])
        return converted

    def force_to_kjmol(self):
        DFT_forces=self.get_force_xyz()
        #CONVERTING THE UNITS OF CP2K TO kJmol-1nm-1:
        DFT_forces['force'] = DFT_forces['force'].apply(lambda x: self.force_unit_convert(x, 49614.626)) #kJmol-1nm-1
        
        return DFT_forces

    #def normalised_and_rms_ref_forces(self, eps = 1e-4):
        #DFT_forces = self.force_to_kjmol() #load DFT forces
        #NORMALISATION:
        #F_dft_norm = []
        #F_dft_rms = []
        #for a in range(len(DFT_forces['force'])):
            #F_dft = np.array([atom[1:] for atom in DFT_forces['force'][a]], dtype=float) # Convert DFT list to numpy array)
            #rms_per_atom = np.sqrt(np.mean(F_dft**2, axis=1)) # calculating the RMS of each
            #to prevent dividing by a value too small, use eps
            #rms_per_atom = np.maximum(rms_per_atom, eps)
            #F_dft_norm.append(F_dft / rms_per_atom[:, None])  # broadcasting divides each atom's 3 components
           # F_dft_rms.append(rms_per_atom)
        #DFT_df = DFT_forces
        #DFT_df['normalised_forces']= F_dft_norm

        #return DFT_df


class ParOpti_preparation:
    from ase.io import read
    
    def __init__(self, name, struct_info, reference_path_DFT, reference_path_GMX):
        self.name = name 
        self.pdb_file = f"{reference_path_DFT}/{name}.pdb"
        self.reference_path_DFT = reference_path_DFT
        self.reference_path_GMX = reference_path_GMX
        self.about_structure = struct_info
        self.aseobj = read(self.pdb_file)
        atom_number = len(self.aseobj)
        self.DFT_data = DFT_reference(self.about_structure, self.name, self.reference_path_DFT,f'{self.name}.pdb')
        self.GMX_data = GMX_read_reference(self.name, self.reference_path_GMX)

        
        
        
    #def default_CoulombMatrix(self):
        #from dscribe.descriptors import CoulombMatrix
        #cm_ds = CoulombMatrix(n_atoms_max=self.atom_number,permutation='eigenspectrum')
        #return cm_ds
    
    def DFT_trajectories_to_descriptors(self, descriptor=None):
        DFT_trajectories = self.DFT_data.DFT_trajectories_ASE_obj()
        descriptor_frames_ = []
        descriptor = None
        if descriptor is None:
            descriptor = cm
        else:
            descriptor = descriptor
            
        #for i,ase_struct in enumerate(DFT_trajectories):
            #dscribe_matrix = descriptor.create(ase_struct)
        #    #dscribe_matrix = np.real(dscribe_matrix)
            #descriptor_frames_.append(dscribe_matrix)
            
        return descriptor_frames_


    def normalised_gmx_forces(self, eps=1e-4):
        gmx_forces = self.GMX_data.extract_forces()
        ## Normalize all gmx forces in the DataFrame
        #for a in range(len(gmx_forces)):
            #frame_num = int(gmx_forces[a][0])
            #print(frame_num)
            #rms_per_atom = DFT_df['rms_per_atoms'][frame_num]
            #_F_gmx = np.array([atom for atom in gmx_forces[a][2]], dtype=float)
            #_rms_per_atom = np.sqrt(np.mean(_F_gmx**2, axis=1)) # calculating the RMS of each
            #rms_per_atom = np.maximum(rms_per_atom, eps)
            #F_gmx_norm.append(_F_gmx / rms_per_atom[:, None])  # broadcasting divides each atom's 3 components
            #F_gmx_rms.append(rms_per_atom)

        Classical_df = pd.DataFrame()
        Classical_df['frame'] = [row[0] for row in gmx_forces]
        Classical_df['parameter'] = [row[1] for row in gmx_forces]
        Classical_df['forces'] = [row[2] for row in gmx_forces]
        
        return Classical_df
    
    from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, mean_squared_log_error, median_absolute_error, mean_absolute_percentage_error

    METRICS = {
    "MSE": mean_squared_error,
    "MAE": mean_absolute_error,
    "R2": r2_score,
    "MSLE": mean_squared_log_error,
    "MEDIANAE": median_absolute_error,
    "MAPE": mean_absolute_percentage_error}
    
    def compute_loss_per_frame(self, metrics="MSE", sample_weight=None, eps=1e-4):
        '''
        All loss for each frame, each parameter each atom {frames x parameters, atom_number, 3 dimensions}
        '''
        
        metricf = self.METRICS.get(metrics.upper()) 
        if metricf is None:
            raise ValueError(f"Unknown metric '{metrics}'")
        
        loss_all = []
        classical_data = self.GMX_data.extract_forces()
        all_dft = self.DFT_data.force_to_kjmol()
        
        #normalising the DFT forces
        all_dft = np.array([np.array(f)[:, 1:].astype(float)for f in all_dft['force']]) #extracting just the forces
        all_dft_flat = all_dft.reshape(-1, 3) #flattening it, this is to obtain the mean and the STD
        mean = all_dft_flat.mean(axis=0) #mean in X, Y, Z
        std = all_dft_flat.std(axis=0) #std in X, Y, Z
        std[std == 0] = 1.0
        all_dft_normalised = (all_dft - mean) / std #normalising the original data

        for frame, param_idx, forces in classical_data:
                dft_forces = all_dft_normalised[frame]
                classic_forces = (np.array(forces) - mean)/std
                per_atom_loss = np.array([metricf(classic_forces[i], dft_forces[i], sample_weight=sample_weight)for i in range(len(classic_forces))])
                loss_all.append(per_atom_loss)
        
        return loss_all  # list of arrays, one per frame
    
    def full_loss_per_D(self, metrics="MSE", sample_weight=None, eps=1e-4):
        """
        Compute average loss per atom for each D value.
        Returns: dictionary {D_value: array of per-atom average losses}
        """
        frame_losses = self.compute_loss_per_frame(metrics=metrics, sample_weight=sample_weight, eps=eps)
        Classical_df = self.normalised_gmx_forces(eps=eps)
    
        # Collect losses per D
        loss_dict = {}
        for loss_per_frame, D_val in zip(frame_losses, Classical_df['parameter']):
            loss_dict.setdefault(D_val, []).append(loss_per_frame)  # append array per frame
        
        return loss_dict


    
    def average_loss_per_D(self, metrics="MSE", sample_weight=None, eps=1e-4):
        """
        Compute average loss per atom for each D value.
        Returns: dictionary {D_value: array of per-atom average losses}
        """
        frame_losses = self.compute_loss_per_frame(metrics=metrics, sample_weight=sample_weight, eps=eps)
        Classical_df = self.normalised_gmx_forces(eps=eps)
    
        # Collect losses per D
        loss_dict = {}
        for loss_per_frame, D_val in zip(frame_losses, Classical_df['parameter']):
            loss_dict.setdefault(D_val, []).append(loss_per_frame)  # append array per frame
        
        # Average per atom for each D
        avg_loss_dict = {D_val: np.mean(np.vstack(losses), axis=0)  # mean over frames, per atom
                         for D_val, losses in loss_dict.items()}
        
        return avg_loss_dict

    def average_type_loss_per_D(self, metrics='MSE', sample_weight=None, eps=1e-4):
        """This takes the averaged per atom(index) and average it further into their atom types"""
        symbols = self.about_structure.get_atom_list()['uff names']
        ave_full_loss = self.average_loss_per_D(metrics=metrics, sample_weight=sample_weight, eps=eps)
        
        D_type_loss = []
        D_vals = []
        for D_val, losses in ave_full_loss.items():
            losses = losses.tolist()
            grouped = defaultdict(list)
            for elem, loss in zip(symbols, losses):
                    grouped[elem].append(loss)
            avg_losses = {elem: np.mean(loss) for elem, loss in grouped.items()}
            #avg_loss_list = [float(f"{avg_losses[elem]:.5g}") for elem in symbols]
            D_type_loss.append(avg_losses)
            D_vals.append(D_val)

        all_ave_loss_df = pd.DataFrame()
        all_ave_loss_df[D_vals] = D_type_loss

        return all_ave_loss_df
            
    
    def gmx_param_shape_extract(self):
        '''NEED TO REDO THIS ONE'''
        param_forcefields = self.GMX_data.ref_FF_parameters()
        param = 'D1'
        #EXTRACTING SHAPE
        lj_shape = param_forcefields[param][['sigma','epsilon']].shape
        bonds_shape = param_forcefields[param]['bonds'][['length', 'force_const']].shape
        angles_shape = param_forcefields[param]['angles'][['angle', 'force_const']].shape
        dihedrals_shape = param_forcefields[param]['dihedrals'][['periodicity', 'phase', 'force_const']].shape

        #adding to the data list
        flat = result.x
        idx = 0
        nb = bonds_shape[0] * bonds_shape[1]
        bonds_new = flat[idx:idx+nb].reshape(bonds_shape)
        idx += nb
        na = angles_shape[0] * angles_shape[1]
        angles_new = flat[idx:idx+na].reshape(angles_shape)
        idx += na
        nd = dihedrals_shape[0] * dihedrals_shape[1]
        dihedrals_new = flat[idx:idx+nd].reshape(dihedrals_shape)

        # Bonds
        bonds_df_new = param_forcefields['D1']['bonds'].copy()
        bonds_df_new['length'] = bonds_new[:, 0]
        bonds_df_new['force_const'] = bonds_new[:, 1]
        bonds_df = bonds_df_new.values.tolist()

        # Angles
        angles_df_new = param_forcefields['D1']['angles'].copy()
        angles_df_new['angle'] = angles_new[:, 0]
        angles_df_new['force_const'] = angles_new[:, 1]
        angles_df = angles_df_new.values.tolist()

        # Dihedrals
        dihedrals_df_new = param_forcefields['D1']['dihedrals'].copy()
        dihedrals_df_new['periodicity'] = dihedrals_new[:, 0]
        dihedrals_df_new['phase'] = dihedrals_new[:, 1]
        dihedrals_df_new['force_const'] = dihedrals_new[:, 2]
        dihedrals_df = dihedrals_df_new.values.tolist()


        
    def write_mdp_file(self, groupid='system', print_directory='./'):
        file_insert = []
        file_insert.append(['integrator               = md ; leap-frog'])
        file_insert.append(['tinit     = 0'])
        file_insert.append(['dt	  = 0.001		; 1 fs'])
        file_insert.append(['nsteps	  = 0   		; 1 μs'])
        file_insert.append([''])
        file_insert.append(['nstxout		= 1		; save coordinates'])
        file_insert.append(['nstvout		= 1			; save velocities'])
        file_insert.append(['nstenergy	= 1		; save energies'])
        file_insert.append(['nstlog		= 1		; update log'])
        file_insert.append(['nstfout         = 1'])
        file_insert.append([''])
        file_insert.append(['cutoff-scheme = Verlet       ; pair list with buffering'])
        file_insert.append(['nstlist       = 10            ; Frequency to update the neighbor list'])
        file_insert.append(['ns-type       = grid         ; Make a grid in the box and only check atoms in neighboring grid cells'])
        file_insert.append(['pbc           = xyz          ; Periodic boundary conditions in all directions'])
        file_insert.append(['rlist		= 1.2		; Neighbor list search cut-off distance (nm)'])
        file_insert.append(['rcoulomb	= 1.2		; Short-range Coulombic interactions cut-off distance (nm)'])
        file_insert.append(['rvdw		= 1.2		; Short-range van der Waals cutoff distance (nm)'])
        file_insert.append(['DispCorr	= EnerPres	; Dispersion correction accounts for cut-off vdW scheme'])
        file_insert.append(['periodic-molecules       = yes'])
        file_insert.append([''])
        file_insert.append(['coulombtype	    = PME	; Particle Mesh Ewald for long-range electrostatics'])
        file_insert.append(['pme_order	    = 4		; cubic interpolation'])
        file_insert.append(['fourierspacing	    = 0.12	; grid spacing for FFT'])
        file_insert.append([''])
        file_insert.append(['tcoupl = nose-hoover'])
        file_insert.append([f'tc-grps                  = {groupid}   ; mention the residue names if wanted multiple thermostats'])
        file_insert.append(['tau_t                    = 0.2 ;mention tau_t as many time as number of residues in tc-grps'])
        file_insert.append(['ref_t                    = 298 ; reference temperature'])
        file_insert.append([''])
        file_insert.append(['pcoupl		= no'])
        file_insert.append([''])
        file_insert.append(['continuation	        = no'])
        file_insert.append([''])
        file_insert.append(['gen_vel		= yes ;'])
        file_insert.append(['gen-temp                 = 298'])
        file_insert.append(['gen-seed                 = -1'])
        file_insert.append([''])
        file_insert.append(['constraints = none'])
        
        file_name = f'{print_directory}/SP.mdp'
        with open(file_name, 'w') as file:
            for row in file_insert:
                row_str = '\t'.join(map(str, row)) 
                file.write(row_str + '\n')

        print(f"{file_name} written")
        
    
    def write_job_script(self,print_directory,budget='Free',time_hours_min_sec='24:00:00', memory_gb='1', core='80', job_name='SP run'):
        file_insert = []
        file_insert.append(['#!/bin/bash -l'])
        file_insert.append(['#$ -S /bin/bash'])
        file_insert.append([f'#$ -l h_rt={time_hours_min_sec}'])
        file_insert.append([f'#$ -l mem={memory_gb}G'])
        file_insert.append([f'#$ -pe mpi {core}'])
        file_insert.append([f'#$ -P {budget}'])
        file_insert.append(['#$ -A UCL_chemE_Yazaydin'])
        file_insert.append([f'#$ -N {job_name}'])
        file_insert.append(['#$ -cwd'])
        file_insert.append(['export OMP_NUM_THREADS=64'])
        file_insert.append(['module unload -f compilers mpi gcc-libs'])
        file_insert.append([''])
        file_insert.append(['module unload -f compilers mpi'])
        file_insert.append(['module load compilers/intel/2018/update3'])
        file_insert.append(['module load mpi/intel/2018/update3/intel'])
        file_insert.append(['module load libmatheval'])
        file_insert.append(['module load flex'])
        file_insert.append(['module load plumed/2.5.2/intel-2018'])
        file_insert.append(['module load gromacs/2019.3/plumed/intel-2018'])
        file_insert.append(['TF_ENABLE_ONEDNN_OPTS=0'])
        file_insert.append([f'gmx_mpi grompp -f SP.mdp -c {self.name}.gro -p topol.top -o SP.tpr -maxwarn 3'])
        file_insert.append(['gmx_mpi mdrun -deffnm SP'])
        file_insert.append(['gmx_mpi dump -f SP.trr | grep \'f\[\' > forces.txt'])
        file_insert.append([''])

        file_name = f'{print_directory}/gmx_sp.sh'
        with open(file_name, 'w') as file:
            for row in file_insert:
                row_str = '\t'.join(map(str, row)) 
                file.write(row_str + '\n')

        print(f"{file_name} written")


        

    def write_gmx_structure_itp_simple(self, groupid='MOF',print_directory='./'):
        atom_con = []
        car_table = self.about_structure.get_atom_list()
        charge_list=self.DFT_data.get_type_charges()
        
        for i, row in enumerate(car_table['index']):
            atom_charge = charge_list[i]
            atom_con.append([row,car_table['uff names'][i],'1',groupid,car_table['uff names'][i],' ',atom_charge])

        bond_con = []
        bond_connectivity = self.about_structure.get_full_bonds()
        for i, row in enumerate(bond_connectivity):
            atom_A = row[0]
            atom_B = row[1]
            index_A = None
            index_B = None
            for a, l1 in enumerate(car_table['labels']):
                if atom_A == l1:
                    index_A = car_table['index'][a]
            for b, l2 in enumerate(car_table['labels']):
                if atom_B == l2:
                    index_B = car_table['index'][b]
            
            bond_con.append([index_A, index_B, '1'])
            
        angle_con = []
        angle_connectivity = self.about_structure.get_full_angles()
        for i, row in enumerate(angle_connectivity):
            atom_A = row[0]
            atom_B = row[1]
            atom_C = row[2]
            
            index_A = None
            index_B = None
            index_C = None
            for a, l1 in enumerate(car_table['labels']):
                if atom_A == l1:
                    index_A = car_table['index'][a]
            for b, l2 in enumerate(car_table['labels']):
                if atom_B == l2:
                    index_B = car_table['index'][b]
            for c, l3 in enumerate(car_table['labels']):
                if atom_C == l3:
                    index_C = car_table['index'][c]
            
            angle_con.append([index_A, index_B,index_C, '1'])
            
        torsion_con = []
        torsion_connectivity = self.about_structure.get_full_torsions()
        for i, row in enumerate(torsion_connectivity):
            atom_A = row[0]
            atom_B = row[1]
            atom_C = row[2]
            atom_D = row[3]
            
            index_A = None
            index_B = None
            index_C = None
            index_D = None
            for a, l1 in enumerate(car_table['labels']):
                if atom_A == l1:
                    index_A = car_table['index'][a]
            for b, l2 in enumerate(car_table['labels']):
                if atom_B == l2:
                    index_B = car_table['index'][b]
            for c, l3 in enumerate(car_table['labels']):
                if atom_C == l3:
                    index_C = car_table['index'][c]
            for d, l3 in enumerate(car_table['labels']):
                if atom_D == l3:
                    index_D = car_table['index'][d]
            
            torsion_con.append([index_A, index_B,index_C,index_D,'3'])
                    
        con_list = []
        con_list.append(['[ moleculetype ]'])
        con_list.append([';','','Name','nrexcl'])
        con_list.append([f'{groupid}','','3'])
        con_list.append([])
        con_list.append(['[ atoms ]'])
        con_list += atom_con
        con_list.append([])
        con_list.append(['[ bonds ]'])
        con_list += bond_con
        con_list.append([])
        con_list.append(['[ angles ]'])
        con_list += angle_con
        con_list.append([])
        con_list.append(['[ dihedrals ]'])
        con_list += torsion_con

        itpname = f'{print_directory}/{groupid}.itp'
        with open(itpname, 'w') as file:
            for row in con_list:
                row_str = '\t'.join(map(str, row)) 
                file.write(row_str + '\n')

        print(f"Data written to {itpname}")



class GMX_read_reference:
    def __init__(self, name, GMX_ref_path):
        self.name = name 
        self.GMX_ref_path = GMX_ref_path
        #self.cif_file = f"{GMX_ref_path}/{name}.cif"
        #self.aseobj = read(self.cif_file)
        #self.atom_nums = len(self.aseobj)
    
    def read_gmx_forces_dump(self, file):
        forces = []
        pattern = re.compile(r'\{([^}]+)\}')
        with open(file, 'r') as f:
            for line in f:
                match = pattern.search(line)
                if match:
                    nums_str = match.group(1)
                    force_vec = [float(x.strip()) for x in nums_str.split(',')]
                    forces.append(force_vec)
        return forces
        
    def extract_forces(self):
        """
        Expects structure: GMX_ref_path/D*/F*D*/forces.txt
        Returns a list: [frame_number, param_number, forces_list]
        """
        gmx_forces = []
    
        # Find all D* folders
        d_folders = [d for d in glob.glob(os.path.join(self.GMX_ref_path, 'D*')) if os.path.isdir(d)]
        d_folders.sort(key=lambda x: int(re.search(r'D(\d+)', x).group(1)))  # sort numerically by D number
    
        for d_folder in d_folders:
            # Look inside D*/ for F*D* folders
            f_folders = [f for f in glob.glob(os.path.join(d_folder, 'F*D*')) if os.path.isdir(f)]
            f_folders.sort(key=lambda x: int(re.search(r'F(\d+)D(\d+)', os.path.basename(x)).group(1)))  # sort by frame
    
            for f_folder in f_folders:
                folder_name = os.path.basename(f_folder)
                match = re.search(r'F(\d+)D(\d+)', folder_name)
                if not match:
                    print(f"Skipping folder {folder_name} (doesn't match pattern F#D#)")
                    continue
                frame = int(match.group(1))
                param_idx = int(match.group(2))
    
                force_file = os.path.join(f_folder, 'forces.txt')
                if not os.path.exists(force_file):
                    print(f"forces.txt not found in {f_folder}, skipping")
                    continue
    
                forces = self.read_gmx_forces_dump(force_file)
                gmx_forces.append([frame, param_idx, forces])
    
        print(f"Extracted forces from {len(gmx_forces)} frames.")
        return gmx_forces

    
    def extract_gmx_force_fields(self, filename):
        import pandas as pd

        sections = {
            'atomtypes': [],
            'bondtypes': [],
            'angletypes': [],
            'dihedraltypes': []
        }

        current_section = None
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith(';'):
                    continue

                # Detect section start
                if line.startswith('[') and line.endswith(']'):
                    sec_name = line[1:-1].strip().lower()
                    current_section = sec_name if sec_name in sections else None
                    continue  # skip the section header line

                if not current_section:
                    continue

                parts = line.split()

                # Parse atomtypes
                if current_section == 'atomtypes':
                    if len(parts) >= 6:
                        atom = parts[0]
                        func = '1'
                        mass = float(parts[2])
                        charge = float(parts[3])
                        ptype = parts[4]
                        sigma = float(parts[5])
                        epsilon = float(parts[6])
                        sections['atomtypes'].append([atom,func, mass, charge, ptype, sigma, epsilon])

                # Parse bondtypes
                elif current_section == 'bondtypes':
                    if len(parts) >= 5:
                        atom1, atom2 = parts[0], parts[1]
                        func = '1'
                        length = float(parts[-2])
                        force_const = float(parts[-1])
                        sections['bondtypes'].append([atom1, atom2,func, length, force_const])

                # Parse angletypes
                elif current_section == 'angletypes':
                    if len(parts) >= 6:
                        atom1, atom2, atom3 = parts[0], parts[1], parts[2]
                        angle = float(parts[-2])
                        func = '1'
                        force_const = float(parts[-1])
                        sections['angletypes'].append([atom1, atom2, atom3,func, angle, force_const])

                # Parse dihedraltypes
                elif current_section == 'dihedraltypes':
                    if len(parts) == 8:
                        atom1, atom2, atom3, atom4 = parts[0], parts[1], parts[2], parts[3]
                        periodicity = float(parts[-3])
                        phase = float(parts[-2])
                        func = '1'
                        force_const = float(parts[-1])
                        sections['dihedraltypes'].append([atom1, atom2, atom3, atom4,func, periodicity, phase, force_const])
                    if len(parts) == 11:
                        atom1, atom2, atom3, atom4 = parts[0], parts[1], parts[2], parts[3]
                        T1, T2, T3, T4, T5, T6 = parts[5],parts[6],parts[7],parts[8],parts[9],parts[10]
                        func = parts[4]
                        sections['dihedraltypes'].append([atom1, atom2, atom3, atom4,func, T1, T2, T3, T4, T5, T6])

        # Convert to DataFrames
        lj_df = pd.DataFrame(
            sections['atomtypes'],
            columns=['atom', 'func','mass', 'charge', 'ptype', 'sigma', 'epsilon']
        )

        bonds_df = pd.DataFrame(
            sections['bondtypes'],
            columns=['atom1', 'atom2','func', 'length', 'force_const']
        )

        angles_df = pd.DataFrame(
            sections['angletypes'],
            columns=['atom1', 'atom2','func', 'atom3', 'angle', 'force_const']
        )

        if len(sections['dihedraltypes'][0]) == 8:
            dihedrals_df = pd.DataFrame(
                sections['dihedraltypes'],
                columns=['atom1', 'atom2', 'atom3', 'atom4','func', 'periodicity', 'phase', 'force_const']
            )
        if len(sections['dihedraltypes'][0]) == 11:
            dihedrals_df = pd.DataFrame(
                sections['dihedraltypes'],
                columns=['atom1', 'atom2', 'atom3', 'atom4','func','T1','T2','T3','T4','T5','T6']
            )

        return lj_df, bonds_df, angles_df, dihedrals_df

    def ref_FF_parameters(self):
        import os, glob

        # Get all D* folders
        d_folders = [d for d in glob.glob(os.path.join(self.GMX_ref_path, 'D*')) if os.path.isdir(d)]
        d_folders.sort(key=lambda x: int(x.split('D')[-1]))  # sort by D number

        param_forcefields = {}

        for d_folder in d_folders:
            # Look at all folders inside D* (no assumptions about name)
            subfolders = [f for f in glob.glob(os.path.join(d_folder, '*')) if os.path.isdir(f)]

            found = False
            for sub in subfolders:
                force_field_path = os.path.join(sub, 'force_field.itp')
                if os.path.exists(force_field_path):
                    lj_df, bonds_df, angles_df, dihedrals_df = self.extract_gmx_force_fields(force_field_path)
                    D_number = os.path.basename(d_folder)  # e.g., "D0", "D1"
                    param_forcefields[D_number] = {
                        'lj': lj_df,
                        'bonds': bonds_df,
                        'angles': angles_df,
                        'dihedrals': dihedrals_df
                    }
                    found = True
                    break  # stop after first valid folder
            if not found:
                print(f"No force_field.itp found in any subfolder of {d_folder}, skipping this D")

        print(f"Obtained force field parameters for {len(param_forcefields)} unique parameter sets.")
        return param_forcefields

    
    


class FFOPTI_tools:

    @staticmethod
    def ref_ffitp_remake(print_directory, param_list):
        ff_list = []
        ff_list.append(['[ defaults ]'])
        ff_list.append(['; nbfunc','comb-rule','gen-pairs','fudgeLJ','fudgeQQ'])
        ff_list.append(['1','2','No','0.5','0.5'])
        ff_list.append(['[ atomtypes ]'])
        ff_list.append([';name','bond_type','mass','charge','ptype','sigma','epsilon'])
        ff_list += param_list['lj'].values.tolist()
        ff_list.append([''])
        ff_list.append(['[ bondtypes ]'])
        ff_list += param_list['bonds'].values.tolist()
        ff_list.append([''])
        ff_list.append(['[ angletypes ]'])
        ff_list += param_list['angles'].values.tolist()
        ff_list.append([''])
        ff_list.append(['[ dihedraltypes ]'])
        ff_list += param_list['dihedrals'].values.tolist()
        #print('starting inversion')
        #ff_list += param_list['inversions']  
        ff_name=f'{print_directory}/force_field.itp'
        with open(ff_name, 'w') as file:
            for row in ff_list:
                row_str = '\t'.join(map(str, row)) 
                file.write(row_str + '\n')
                    
        print("force_field.itp file printed")

    def gmx_SP_job_submission(self, workdir, job_script="gmx_sp.sh"):
        """
        Submits a single job in workdir and returns the numeric SGE job ID.
        """
    
        script_path = os.path.join(workdir, job_script)
        if not os.path.isfile(script_path):
            raise RuntimeError(f"Missing {job_script} in {workdir}")
    
        result = subprocess.run(
            ["qsub", job_script],
            cwd=workdir,
            capture_output=True,
            text=True,
            check=True
        )
    
        import re
        match = re.search(r"Your job (\d+)", result.stdout)
        if not match:
            raise RuntimeError(f"Could not parse job ID from qsub output:\n{result.stdout}")
    
        job_id = match.group(1)
        print(f"Submitted job {job_id} in {workdir}")
    
        return job_id



    def wait_for_sp_jobs(self, job_id_list, poll_interval=120):
        """
        Args:
            job_id_list (list of str): SGE job IDs
            poll_interval (int): Seconds between qstat checks
        """
    
        job_ids = [str(jid) for jid in job_id_list if jid]
    
        if not job_ids:
            print("No jobs to wait for.")
            return
    
        print(f"Waiting for {len(job_ids)} jobs to finish...")
    
        while True:
            running = []
    
            for job_id in job_ids:
                result = subprocess.run(
                    ["qstat","-j", job_id],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL
                )
    
                if result.returncode == 0:
                    running.append(job_id)
    
            if not running:
                break
    
            print(f"Still running: {running}")
            time.sleep(poll_interval)
    
        print("All jobs in batch completed.")



    def submit_param_batches(self, GMX_path, param, poll_interval=120):
    
        for a, _ in enumerate(param):
            D = f"D{a+1}"
            print(f"\nSubmitting batch for {D}")
    
            D_path = os.path.join(GMX_path, D)
            if not os.path.isdir(D_path):
                print(f"Warning: {D_path} does not exist, skipping")
                continue
    
            job_ids = []
    
            for subdir in sorted(os.listdir(D_path)):
                if subdir.startswith("F"):
                    workdir = os.path.join(D_path, subdir)
                    script_path = os.path.join(workdir, "gmx_sp.sh")
    
                    if os.path.isfile(script_path):
                        job_id = self.gmx_SP_job_submission(workdir)
                        job_ids.append(job_id)
                    else:
                        print(f"Warning: Job script not found in {workdir}")
    
            if job_ids:
                self.wait_for_sp_jobs(job_ids, poll_interval=poll_interval)
            else:
                print(f"No jobs submitted for {D}, skipping wait.")


    @staticmethod 
    def job_script_all(name, GMX_path, new_D, frame):
        file_insert = []
        file_insert.append(['#!/bin/bash -l'])
        #file_insert.append(['module unload -f compilers mpi gcc-libs'])
        file_insert.append([''])
        file_insert.append(['module unload -f compilers mpi'])
        file_insert.append(['module load compilers/intel/2018/update3'])
        file_insert.append(['module load mpi/intel/2018/update3/intel'])
        file_insert.append(['module load gromacs/2019.3/intel-2018'])
        #file_insert.append(['module load flex'])
        #file_insert.append(['module load plumed/2.5.2/intel-2018'])
        #file_insert.append(['module load gromacs/2019.3/plumed/intel-2018'])
        #file_insert.append(['TF_ENABLE_ONEDNN_OPTS=0'])

        for i in range(len(frame)):
            file_insert.append([f'cd /{GMX_path}/{new_D}/F{i}{new_D}'])
            file_insert.append([f'gmx grompp -f SP.mdp -c {name}.gro -p topol.top -o SP.tpr -maxwarn 3'])
            file_insert.append(['gmx mdrun -deffnm SP'])
            file_insert.append(['gmx dump -f SP.trr | grep \'f\[\' > forces.txt'])
            file_insert.append([''])

        file_name = f'{GMX_path}/{new_D}/run_all'
        with open(file_name, 'w') as file:
            for row in file_insert:
                row_str = '\t'.join(map(str, row)) 
                file.write(row_str + '\n')

        print(f"{file_name} written")
        return(f"job script for {new_D} created in {GMX_path}/{new_D}")

    @staticmethod 
    def job_script_mpi(name, GMX_path, new_D, idx):
        file_insert = []
        file_insert.append(['#!/bin/bash -l'])
        file_insert.append([''])
        file_insert.append(['module unload -f compilers mpi'])
        file_insert.append(['module load compilers/intel/2018/update3'])
        file_insert.append(['module load mpi/intel/2018/update3/intel'])
        file_insert.append(['module load libmatheval '])
        file_insert.append(['module load flex'])
        file_insert.append(['module load plumed/2.5.2/intel-2018'])
        file_insert.append(['module load gromacs/2019.3/plumed/intel-2018'])

        for i in (idx):
            file_insert.append([f'cd /{GMX_path}/{new_D}/F{i}{new_D}'])
            file_insert.append([f'gmx_mpi grompp -f SP.mdp -c {name}.gro -p topol.top -o SP.tpr -maxwarn 3'])
            file_insert.append(['gmx_mpi mdrun -deffnm SP'])
            file_insert.append(['gmx_mpi dump -f SP.trr | grep \'f\[\' > forces.txt'])
            file_insert.append([''])

        file_name = f'{GMX_path}/{new_D}/run_all'
        with open(file_name, 'w') as file:
            for row in file_insert:
                row_str = '\t'.join(map(str, row)) 
                file.write(row_str + '\n')

        print(f"{file_name} written")
        return(f"job script for {new_D} created in {GMX_path}/{new_D}")


    def submit_all_jobs_in_D(self, GMX_path, D_number):
        D_path = os.path.join(GMX_path, D_number)    
        if not os.path.isdir(D_path):
            raise RuntimeError(f"No such D folder: {D_path}")
    
        job_ids = []
    
        print(f"\nSubmitting all jobs in {D_number}")
    
        for subdir in sorted(os.listdir(D_path)):
            if not subdir.startswith("F"):
                continue
    
            workdir = os.path.join(D_path, subdir)
            script_path = os.path.join(workdir, "gmx_sp.sh")
    
            if not os.path.isfile(script_path):
                print(f"Warning: gmx_sp.sh not found in {workdir}, skipping.")
                continue
    
            job_id = self.gmx_SP_job_submission(workdir)
            job_ids.append(job_id)
    
        if not job_ids:
            print(f"No jobs submitted in {D_number}")
    
        return job_ids


    @staticmethod       
    def write_gro(name, ase_obj, write_path, Structure):
        atom_UFF_names = Structure.get_atom_list()['uff names']
        io.write(f'{name}.gro',ase_obj)
        file = f'{name}.gro'

        char_list = []
        with open(file, 'r') as file:
        # Read each line
            for line in file:
                # Strip the newline character at the end of each line
                line = line.rstrip('\n')
                # Convert the line into a list of characters
                char_list.append(list(line))
                
        residue_number= [sub[0:5] for sub in char_list[2:-1]]
        residue_name= [sub[5:10] for sub in char_list[2:-1]]
        atom_name= [sub[10:15] for sub in char_list[2:-1]]
        atom_number= [sub[15:20] for sub in char_list[2:-1]]
        x_coord= [sub[20:28] for sub in char_list[2:-1]]
        y_coord= [sub[28:36] for sub in char_list[2:-1]]
        z_coord= [sub[36:44] for sub in char_list[2:-1]]
        final = [char_list[-1]]
    
        atom_name = [list(f"{name:<5}"[:5]) for name in atom_UFF_names]
   
        #COMBINING THEM BACK
        combine = list(zip(residue_number, residue_name, atom_name, atom_number, x_coord, y_coord, z_coord))
        flat_comb = [
            [item for sublist in tup for item in sublist]
            for tup in combine
        ]      
        file_content = []
        file_content.append(char_list[0])
        file_content.append(str(len(residue_number)))
        file_content += flat_comb
        file_content.append(char_list[-1])
        
        with open(f'{write_path}/{name}.gro', 'w') as file:
            for row in file_content:
                row_str = ''.join(map(str, row)) 
                file.write(row_str + '\n')
            
        print(f'grofile created')


class single_loss_calc:
    '''
    for obtaining loss value for a single parameter in the same format as the ones on top
    '''
    def __init__(self, name, struct_info, reference_path_DFT, d_folder):
        self.name = name
        self.about_structure = struct_info
        self.reference_path_DFT = reference_path_DFT
        self.d_folder = d_folder
        self.DFT_data = DFT_reference(self.about_structure, self.name, self.reference_path_DFT,f'{self.name}.pdb')
    
    def read_gmx_forces_dump(self, file):
        forces = []
        pattern = re.compile(r'\{([^}]+)\}')
        with open(file, 'r') as f:
            for line in f:
                match = pattern.search(line)
                if match:
                    nums_str = match.group(1)
                    force_vec = [float(x.strip()) for x in nums_str.split(',')]
                    forces.append(force_vec)
        return forces
    
    def extract_forces_singles(self):
        """
        Expects structure: GMX_ref_path/D*/F*D*/forces.txt
        Returns a list: [frame_number, param_number,forces_list]
        """
        gmx_forces = []
        
        # Find all D* folders
        d_folder = self.d_folder
        
        # Look inside D_folder/ for F*D* folders
        f_folders = [f for f in glob.glob(os.path.join(d_folder, 'F*D*')) if os.path.isdir(f)]
        f_folders.sort(key=lambda x: int(re.search(r'F(\d+)D(\d+)', os.path.basename(x)).group(1)))  # sort by frame
        
        for f_folder in f_folders:
            folder_name = os.path.basename(f_folder)
            match = re.search(r'F(\d+)D(\d+)', folder_name)
            if not match:
                print(f"Skipping folder {folder_name} (doesn't match pattern F#D#)")
                continue
            frame = int(match.group(1))
            param_idx = int(match.group(2))
        
            force_file = os.path.join(f_folder, 'forces.txt')
            if not os.path.exists(force_file):
                print(f"forces.txt not found in {f_folder}, skipping")
                continue
        
            forces = self.read_gmx_forces_dump(force_file)
            gmx_forces.append([frame, param_idx, forces])
        
        print(f"Extracted forces from {len(gmx_forces)} frames.")
        return gmx_forces

    def normalised_gmx_forces(self):
        gmx_forces = self.extract_forces_singles()

        Classical_df = pd.DataFrame()
        Classical_df['frame'] = [row[0] for row in gmx_forces]
        Classical_df['parameter'] = [row[1] for row in gmx_forces]
        Classical_df['forces'] = [row[2] for row in gmx_forces]
        
        return Classical_df

    from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score, mean_squared_log_error, median_absolute_error, mean_absolute_percentage_error

    METRICS = {
    "MSE": mean_squared_error,
    "MAE": mean_absolute_error,
    "R2": r2_score,
    "MSLE": mean_squared_log_error,
    "MEDIANAE": median_absolute_error,
    "MAPE": mean_absolute_percentage_error}
    

    def compute_loss_per_frame(self, metrics="MSE", sample_weight=None, eps=1e-4):
        '''
        All loss for each frame, each parameter each atom {frames x parameters, atom_number, 3 dimensions}
        '''
        
        metricf = self.METRICS.get(metrics.upper()) 
        if metricf is None:
            raise ValueError(f"Unknown metric '{metrics}'")
        
        loss_all = []
        classical_data = self.extract_forces_singles()
        all_dft = self.DFT_data.force_to_kjmol()
        
        #normalising the DFT forces
        all_dft = np.array([np.array(f)[:, 1:].astype(float)for f in all_dft['force']]) #extracting just the forces
        all_dft_flat = all_dft.reshape(-1, 3) #flattening it, this is to obtain the mean and the STD
        mean = all_dft_flat.mean(axis=0) #mean in X, Y, Z
        std = all_dft_flat.std(axis=0) #std in X, Y, Z
        std[std == 0] = 1.0
        all_dft_normalised = (all_dft - mean) / std #normalising the original data

        for frame, param_idx, forces in classical_data:
                dft_forces = all_dft_normalised[frame]
                classic_forces = (np.array(forces) - mean)/std
                per_atom_loss = np.array([metricf(classic_forces[i], dft_forces[i], sample_weight=sample_weight)for i in range(len(classic_forces))])
                loss_all.append(per_atom_loss)
        
        return loss_all  # list of arrays, one per frame
    
    def full_loss_per_D(self, metrics="MSE", sample_weight=None, eps=1e-4):
        """
        Compute average loss per atom for each D value.
        Returns: dictionary {D_value: array of per-atom average losses}
        """
        frame_losses = self.compute_loss_per_frame(metrics=metrics, sample_weight=sample_weight, eps=eps)
        Classical_df = self.normalised_gmx_forces()
    
        # Collect losses per D
        loss_dict = {}
        for loss_per_frame, D_val in zip(frame_losses, Classical_df['parameter']):
            loss_dict.setdefault(D_val, []).append(loss_per_frame)  # append array per frame
        
        return loss_dict


    
    def average_loss_per_D(self, metrics="MSE", sample_weight=None, eps=1e-4):
        """
        Compute average loss per atom for each D value.
        Returns: dictionary {D_value: array of per-atom average losses}
        """
        frame_losses = self.compute_loss_per_frame(metrics=metrics, sample_weight=sample_weight, eps=eps)
        Classical_df = self.normalised_gmx_forces()
    
        # Collect losses per D
        loss_dict = {}
        for loss_per_frame, D_val in zip(frame_losses, Classical_df['parameter']):
            loss_dict.setdefault(D_val, []).append(loss_per_frame)  # append array per frame
        
        # Average per atom for each D
        avg_loss_dict = {D_val: np.mean(np.vstack(losses), axis=0)  # mean over frames, per atom
                         for D_val, losses in loss_dict.items()}
        
        return avg_loss_dict

    def average_type_loss_per_D(self, metrics='MSE', sample_weight=None, eps=1e-4):
        """This takes the averaged per atom(index) and average it further into their atom types"""
        symbols = self.about_structure.get_atom_list()['uff names']
        ave_full_loss = self.average_loss_per_D(metrics=metrics, sample_weight=sample_weight, eps=eps)
        
        D_type_loss = []
        D_vals = []
        for D_val, losses in ave_full_loss.items():
            losses = losses.tolist()
            grouped = defaultdict(list)
            for elem, loss in zip(symbols, losses):
                    grouped[elem].append(loss)
            avg_losses = {elem: np.mean(loss) for elem, loss in grouped.items()}
            #avg_loss_list = [float(f"{avg_losses[elem]:.5g}") for elem in symbols]
            D_type_loss.append(avg_losses)
            D_vals.append(D_val)

        all_ave_loss_df = pd.DataFrame()
        all_ave_loss_df[D_vals] = D_type_loss

        return all_ave_loss_df

        

###FOR BAYESIAN OPTIMISATION MODEL
    #AQUISITION FUNCTIONS
def expected_improvement(X, gp, y_best, xi=0.01):
    """
    Standard EI. Balances exploration/exploitation.
    """
    mu, sigma = gp.predict(X, return_std=True)
    with np.errstate(divide='warn'):
        imp = y_best - mu - xi
        Z = imp / sigma
        ei = imp * norm.cdf(Z) + sigma * norm.pdf(Z)
        ei[sigma == 0.0] = 0.0
    return ei

def probability_of_improvement(X, gp, y_best, xi=0.01):
    """
    Good for strictly finding a better point, regardless of magnitude.
    """
    mu, sigma = gp.predict(X, return_std=True)
    with np.errstate(divide='warn'):
        imp = y_best - mu - xi
        Z = imp / sigma
        pi = norm.cdf(Z)
        pi[sigma == 0.0] = 0.0
    return pi

def lower_confidence_bound(X, gp, y_best=None, kappa=1.96):
    """
    LCB = Mean - Kappa * Std.
    Kappa ~ 1.96 corresponds to 95% confidence exploration.
    Higher Kappa = More exploration.
    """
    mu, sigma = gp.predict(X, return_std=True)
    return -(mu - kappa * sigma) # We return negative because we maximize acquisition

def log_expected_improvement(X, gp, y_best, xi=0.01):
    """
    Log Expected Improvement (LogEI).
    Mathematically more robust for surfaces with large dynamic ranges.
    """
    mu, sigma = gp.predict(X, return_std=True)

    # Avoid division by zero
    sigma = np.maximum(sigma, 1e-9)

    imp = y_best - mu - xi
    z = imp / sigma

    # LogEI formulation
    # We use the log of the standard EI formula components
    # log(sigma * (z * Phi(z) + phi(z)))
    phi_z = norm.pdf(z)
    Phi_z = norm.cdf(z)

    # To maintain numerical stability for very small values:
    ei = sigma * (imp/sigma * Phi_z + phi_z)

    # We return the log, but since we usually maximize acquisition,
    # we take log of a very small positive number.
    # For the optimizer, we handle it carefully:
    return np.log(ei + 1e-300)