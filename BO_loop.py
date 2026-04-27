from all_needed_modules import *
from all_FFOPT_classes import *
from input_GF_lj import *
import warnings
from sklearn.preprocessing import StandardScaler
from sklearn.gaussian_process.kernels import DotProduct, WhiteKernel
warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)


#information extraction
Structure = About_Structure(name=Structure_name, carfile=structure_carfile, mdffile=structure_mdffile, chargefile=charge_file)
DFT_ref = DFT_reference(Structure, Structure_name, DFT_ref_path, structure_file='GF334.pdb')
GMX_ref = GMX_read_reference(Structure_name, GMX_ref_path)
initial_opti_info = ParOpti_preparation(Structure_name,Structure, DFT_ref_path, GMX_ref_path)
#Obtaining the initial GMX parameters
param = GMX_ref.ref_FF_parameters()
starting_step = len(param)
print('information obtained')


##########################################################################
#GETTING TRAJECTORIES
##########################################################################
traj_dcd_path = f'{DFT_ref_path}/md-pos-1.dcd'
frame_list = read(traj_dcd_path, index=':') #turning each frames into ase objects
# reading indx list
with open(f"{GMX_ref_path}/ffopti_backup/numbers.json", "r") as f:
    indx = json.load(f)
frames = [frame_list[i] for i in indx]

ref = frames[0] #choosing a reference frame for NVT to set PBC and cell dimensions

for atoms in frames:
    atoms.set_cell(size)
    atoms.set_pbc(True)
    atoms.wrap()  # optional, wraps positions inside box
n_of_types = Structure.atom_types()

##########################################################################
#Extracting parameters being optimised
##########################################################################
order = ['lj', 'bond', 'angle', 'dihedral']

flattened_parameter = []
all_indices = []
for a, l1 in enumerate(list(param)):
    # enforce ordering regardless of input
    active = [x for x in order if x in train]

    parts = []
    indices = {}
    offset = 0

    for key in active:

        if key == 'lj':
            arr = param[f'{l1}']['lj'][['sigma', 'epsilon']].to_numpy().ravel()
            parts.append(arr)

            n = len(arr)
            if 'lj_sigma' not in indices:
                indices['lj_sigma']   = list(range(offset, offset+n, 2))
                indices['lj_epsilon'] = list(range(offset+1, offset+n, 2))
            offset += n

        elif key == 'bond':
            arr = param[f'{l1}']['bonds'][['length', 'force_const']].to_numpy().ravel()
            parts.append(arr)

            n = len(arr)
            if 'bond_rij' not in indices:
                indices['bond_rij'] = list(range(offset, offset+n, 2))
                indices['bond_kij'] = list(range(offset+1, offset+n, 2))
            offset += n

        elif key == 'angle':
            arr = param[f'{l1}']['angles'][['angle', 'force_const']].to_numpy().ravel()
            parts.append(arr)

            n = len(arr)
            if 'angle_angle' not in indices:
                indices['angle_angle'] = list(range(offset, offset+n, 2))
                indices['angle_kijk'] = list(range(offset+1, offset+n, 2))
            offset += n

        elif key == 'dihedral':
            arr = param[f'{l1}']['dihedrals'][['T1', 'T2', 'T3', 'T4', 'T5', 'T6']].to_numpy().ravel()
            parts.append(arr)

            n = len(arr)
            #if 'dih_T1' not in indices:
            #indices['dih_T1'] = list(range(offset, offset+n, 2))
            #indices['dih_T2'] = list(range(offset+1, offset+n, 2))
            #offset += n

    # safely concatenate
    flatie = np.concatenate(parts) if parts else np.array([])

    flattened_parameter.append(flatie)
    all_indices.append(indices)

feature_length = len(flattened_parameter[0])
print(feature_length)

#EXTRACTING SHAPE
lj_shape = param[reference_force_field]['lj'][['sigma','epsilon']].shape
bonds_shape = param[reference_force_field]['bonds'][['length', 'force_const']].shape
angles_shape = param[reference_force_field]['angles'][['angle', 'force_const']].shape
#dihedrals_shape = param['D3']['dihedrals'][['periodicity', 'phase', 'force_const']].shape
dihedrals_shape = param[reference_force_field]['dihedrals'][['T1', 'T2', 'T3','T4', 'T5', 'T6']].shape

##########################################################################
#EXTRACTING THE LOSSES
##########################################################################
file_path = f"{GMX_ref_path}/ffopti_backup/param_type_loss.pkl"
if os.path.exists(file_path):
    average_loss = pd.read_pickle(file_path)
else:
    average_loss = initial_opti_info.average_type_loss_per_D(metrics="MAPE")
    average_loss.to_pickle(f"{GMX_ref_path}/ffopti_backup/param_type_loss.pkl")

#If we did read from the pkl file, we now check if there are any losses that were not updated
D_vals_files = list(param.keys())
pkl_list = average_loss.columns

missing_param_loss = []
for a in D_vals_files:
    num=a[1:]
    if int(num) not in pkl_list:
        missing_param_loss.append(num)
print(f'{len(missing_param_loss)} missing from backup file')

for a in missing_param_loss:
    d_folder = f'{GMX_ref_path}/D{a}'
    single_calc_tools = single_loss_calc(Structure_name,Structure, DFT_ref_path,d_folder)
    new_loss = single_calc_tools.average_type_loss_per_D(metrics="MAPE")
    average_loss = average_loss.join(new_loss)

print('loss obtained')
##########################################################################
#Setting bounds
##########################################################################
refere_indx = int(reference_force_field[1:])-1
initial_guess = np.array(flattened_parameter)[refere_indx]
bounds = []

lj_sigma = set(indices.get('lj_sigma', []))
lj_epsilon = set(indices.get('lj_epsilon', []))

bond_rij = set(indices.get('bond_rij', []))
bond_kij = set(indices.get('bond_kij', []))


angle_angle = set(indices.get('angle_angle', []))
angle_kijk = set(indices.get('angle_kijk', []))

#dih_t1 = set(indices.get('dih_T1', []))
#dih_t2 = set(indices.get('dih_T2', []))

for i, x in enumerate(initial_guess):

    # LJ rules
    if i in lj_sigma or i in bond_rij:
        bounds.append((0.90 * x, 1.10 * x))
    if i in lj_epsilon:
        bounds.append((0.75 * x, 1.25 * x))

    # bond rules
    if 'bond' in train and i in bond_rij:
        bounds.append((0.90 * x, 1.10 * x))
    if 'bond' in train and i in bond_kij:
        bounds.append((0.70 * x, 1.30 * x))

        
    # angle rules
    if 'angle' in train and i in angle_angle:
        bounds.append((0.90 * x, 1.10 * x))
    if 'angle' in train and i in angle_kijk:
        bounds.append((0.70 * x, 1.30 * x))

    # dihedral rules
    #if 'dihedral' in train and i in dih_T1:
        #bounds.append((0.70 * x, 1.30 * x))
    #if 'dihedral' in train and i in dih_T2:
        #bounds.append((0.70 * x, 1.30 * x))
        

    # fallback
    #if len(bounds) < i + 1:
        #bounds.append((0.95 * x, 1.05 * x))
print('bound set')
##########################################################################
#Setting up loop
##########################################################################
# Find iteration index with lowest total loss
total_loss_per_iter = average_loss.sum(axis=0)
best_iter = total_loss_per_iter.idxmin()

for a in range(max_iter - starting_step):
    param = GMX_ref.ref_FF_parameters()
    total_loss_per_iter = average_loss.sum(axis=0)
    best_iter = total_loss_per_iter.idxmin()
    current_iter = len(average_loss.columns) - 1
    print(current_iter, best_iter)
    if current_iter - int(best_iter) <= 200:
        print('starting epoch:',len(param))
        ########################
        #Setting up data for BO
        ########################
        flattened_parameter = []
        for a, l1 in enumerate(list(param)):
            # enforce ordering regardless of input
            active = [x for x in order if x in train]
        
            parts = []
            indices = {}
            offset = 0
        
            for key in active:
        
                if key == 'lj':
                    arr = param[f'{l1}']['lj'][['sigma', 'epsilon']].to_numpy().ravel()
                    parts.append(arr)
        
                    n = len(arr)
                    offset += n
        
                elif key == 'bond':
                    arr = param[f'{l1}']['bonds'][['length', 'force_const']].to_numpy().ravel()
                    parts.append(arr)
        
                    n = len(arr)
                    offset += n
        
                elif key == 'angle':
                    arr = param[f'{l1}']['angles'][['angle', 'force_const']].to_numpy().ravel()
                    parts.append(arr)
        
                    n = len(arr)
                    offset += n
        
                elif key == 'dihedral':
                    arr = param[f'{l1}']['dihedrals'][['T1', 'T2', 'T3', 'T4', 'T5', 'T6']].to_numpy().ravel()
                    parts.append(arr)        
            # safely concatenate
            flatie = np.concatenate(parts) if parts else np.array([])
        
            flattened_parameter.append(flatie)

        
        loss_scalar = []
        for keys in param.keys():
            param_code = str(keys)
            D_val = int((param_code)[1:]) #just the number
            loss = average_loss[D_val]
            loss_scalar.append(loss)
        
        
        y_train = np.array([
            [d[atom] for atom in atom_types]
            for d in loss_scalar
        ])
        summed_y_train = [sum(l1) for l1 in y_train]
        y_train = summed_y_train
        
        #print("y_train shape:", y_train.shape)
        
        X = np.array(flattened_parameter)
        X_train = np.array(flattened_parameter)
        
        
        print("X_train shape:", X_train.shape)
        #print("y_train shape:", y_train.shape)
        
    
    
        ###########################################
        #BUILDING THE GAUSSIAN PROCESS OF BO METHOD
        ###########################################
        # 2. Initialize the Scalers
        scaler_X = StandardScaler()
        scaler_y = StandardScaler()
         
        # 3. Scale the training data
        # It is vital to fit the scaler only on existing data and transform it
        X_scaled = scaler_X.fit_transform(X_train)
        y_train = np.array(y_train)          # convert list → numpy array
        y_train = y_train.reshape(-1, 1) 
        y_scaled = scaler_y.fit_transform(y_train)
         
        # 4. Define and fit the GP model with an improved kernel
        # Added a WhiteKernel to account for noise in the GROMACS simulations
        #kernel = C(1.0) * Matern(length_scale=1.0, nu=1.5) + WhiteKernel(noise_level=1e-5)
        #kernel = C(1.0, (1e-3, 1e4)) * Matern(
            #length_scale=np.ones(feature_length),          # 8 features in X_train
            #length_scale_bounds=(1e-2, 1e2),  # optimizer can tune each scale
            #nu=2.5                            # smoothness of function
        #) + WhiteKernel(noise_level=1e-5)
        kernel = C(1.0, (1e-3, 1e4)) * Matern(
            length_scale=np.ones(feature_length),          # 8 features in X_train
            length_scale_bounds=(1e-2, 1e2),  # optimizer can tune each scale
            nu=2.5                            # smoothness of function
        )          
        #gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10,random_state=42)
        gp = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=10, alpha=1e-6, normalize_y=True, random_state=42)
        gp.fit(X_scaled, y_scaled)
         
        # 5. When predicting the next sample, remember to scale the input
        # and inverse-transform the output to get physical units back
        def predict_loss(new_params):
            # Ensure input is a 2D array for the scaler
            new_params_scaled = scaler_X.transform(np.atleast_2d(new_params))
            y_pred_scaled, sigma_scaled = gp.predict(new_params_scaled, return_std=True)
            # Transform back to original units
            y_pred = scaler_y.inverse_transform(y_pred_scaled.reshape(-1, 1))
            # Sigma needs scaling by the standard deviation of y
            sigma = sigma_scaled * np.sqrt(scaler_y.var_)
            print(f'predicted_loss: {y_pred.flatten()}, sigma: {sigma.flatten()}')
            return y_pred.flatten(), sigma.flatten()
        
        # NOTE: normalize_y=True fixes the ConvergenceWarning by scaling target data
        
        print("Learned kernel:", gp.kernel_)
    
        ACQ_FUNC = log_expected_improvement
        # Direction of optimisation
        def min_obj(X):
            return -ACQ_FUNC(X.reshape(1, -1), gp, np.min(y_train)) #use this if just the average
    
        # Optimization Loop: Select the number of BO rounds  ---
        n_iterations = 100
        print("\n--- Starting Optimization ---")
        
        for i in range(n_iterations):
            y_best = np.min(y_train)
        
            best_x_next = None
            max_acq = -np.inf
        
            # Restart optimizer to avoid local traps in acquisition surface
            for _ in range(10):
                bounds = np.array(bounds)  # shape should be (n_features, 2)
                x0 = np.random.uniform(bounds[:, 0], bounds[:, 1])
                res = minimize(min_obj, x0=x0, bounds=bounds, method='L-BFGS-B')
                if -res.fun > max_acq:
                    max_acq = -res.fun
                    best_x_next = res.x
        
        print(f"Iteration {i+1}: Suggested new parameters: {best_x_next}")
        loss_best_x = predict_loss(best_x_next)
        print(f"loss_best_x")
        idx = 0
        x = best_x_next
        #LJ parameter
        if 'lj' in train:
            nlj = lj_shape[0] * lj_shape[1]
            lj_new = x[idx:idx + nlj].reshape(lj_shape)
            idx += nlj
        
            lj_df_new = param[reference_force_field]['lj'].copy()
            lj_df_new['sigma'] = lj_new[:, 0]
            lj_df_new['epsilon'] = lj_new[:, 1]
        
        else:
            lj_df_new = param[reference_force_field]['lj'].copy()
        
        #Bond Parameter
        if 'bond' in train:
            nbond = bonds_shape[0] * bonds_shape[1]
            bonds_new = x[idx:idx + nbond].reshape(bonds_shape)
            idx += nbond
            
            bonds_df_new = param[reference_force_field]['bonds'].copy()
            bonds_df_new['length'] = bonds_new[:, 0]
            bonds_df_new['force_const'] = bonds_new[:, 1]
        else:
            bonds_df_new = param[reference_force_field]['bonds'].copy()
        
        #angle Parameter
        if 'angle' in train:
            nangle = angles_shape[0] * angles_shape[1]
            angles_new = x[idx:idx + nangle].reshape(angles_shape)
            idx += nangle
            
            angles_df_new = param[reference_force_field]['angles'].copy()
            angles_df_new['angle'] = angles_new[:, 0]
            angles_df_new['force_const'] = angles_new[:, 1]
        else:
            angles_df_new = param[reference_force_field]['angles'].copy()
        
        
        #angle Parameter
        if 'dihedral' in train and dihedrals_shape[1] == 3:
            ndihedral = dihedrals_shape[0] * dihedrals_shape[1]
            dihedrals_new = x[idx:idx + ndihedral].reshape(dihedrals_shape)
            idx += ndihedral
            
            dihedrals_df_new = param[reference_force_field]['dihedrals'].copy()
            dihedrals_df_new['periodicity'] = dihedrals_new[:, 0]
            dihedrals_df_new['phase'] = dihedrals_new[:, 1]
            dihedrals_df_new['force_const'] = dihedrals_new[:, 2]
        elif 'dihedral' in train and dihedrals_shape[1] == 6:
            ndihedral = dihedrals_shape[0] * dihedrals_shape[1]
            dihedrals_new = x[idx:idx + ndihedral].reshape(dihedrals_shape)
            idx += ndihedral
            
            dihedrals_df_new = param[reference_force_field]['dihedrals'].copy()
            dihedrals_df_new['T1'] = dihedrals_new[:, 0]
            dihedrals_df_new['T2'] = dihedrals_new[:, 1]
            dihedrals_df_new['T3'] = dihedrals_new[:, 2]
            dihedrals_df_new['T4'] = dihedrals_new[:, 3]
            dihedrals_df_new['T5'] = dihedrals_new[:, 4]
            dihedrals_df_new['T6'] = dihedrals_new[:, 5]
        else:
            dihedrals_df_new = param[reference_force_field]['dihedrals'].copy()
        
        new_D = f'D{len(list(param))+1}'
        param[new_D] = {
            'lj': lj_df_new,
            'bonds': bonds_df_new,
            'angles': angles_df_new,
            'dihedrals': dihedrals_df_new
        }
        
        print(param[new_D])
    
        #Setting up loop for to create initial GMX reference data
        for a, l1 in enumerate(frames):
            frame_number = f'F{indx[a]}'
            D_number = new_D
        
            # Create the directory D*/F*D* if it doesn't exist
            folder_path = os.path.join(GMX_ref_path, D_number, f'{frame_number}{D_number}')
            os.makedirs(folder_path, exist_ok=True)
        
                #.gro file
            FFOPTI_tools.write_gro(Structure_name, l1, f'{folder_path}', Structure)
            param_list = param[f'{D_number}']
                #MOF.itp file
            subprocess.run(["cp", f"{GMX_ref_path}/MOF.itp", f"{folder_path}"],check=True)
                #Force_field.itp file
            FFOPTI_tools.ref_ffitp_remake(f'{folder_path}', param_list)
                #topol.top
            Structure.write_top_file(print_directory=f'{folder_path}')
                #SP.mdp
            initial_opti_info.write_mdp_file(print_directory=f'{folder_path}')
            print(f'files created for {D_number}{frame_number}')
        print('initial files created')
    
        tools = FFOPTI_tools()
        tools.job_script_mpi(Structure_name, GMX_ref_path, new_D, indx)
        print('submitted')
        subprocess.run(["bash", "-l", f'{GMX_ref_path}/{new_D}/run_all'],check=True)
    
        #calculating the loss for the new output
        d_folder = f'{GMX_ref_path}/{new_D}'
        single_calc_tools = single_loss_calc(Structure_name,Structure, DFT_ref_path,d_folder)
        new_loss = single_calc_tools.average_type_loss_per_D(metrics="MAPE")
        average_loss = average_loss.join(new_loss)
        average_loss.to_pickle(f"{GMX_ref_path}/ffopti_backup/param_type_loss.pkl")

    else:
        print('no further improvement')
        break
        
