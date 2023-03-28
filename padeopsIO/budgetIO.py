import numpy as np
import os
import re
import warnings
import glob
from scipy.io import savemat
import matplotlib.pyplot as plt

try: 
    import f90nml
except ImportError: 
    warnings.warn("Could not find package f90nml in system path. ")
    # TODO - handle this somehow
    f90nml = None

import padeopsIO.budgetkey as budgetkey  # defines key pairing
import padeopsIO.inflow as inflow  # interface to retrieve inflow profiles
import padeopsIO.turbineArray as turbineArray  # reads in a turbine array similar to turbineMod.F90

class BudgetIO(): 
    """
    # TODO fix class docstring
    """

    key = budgetkey.get_key()
    key_xy = budgetkey.get_key_xy()

    def __init__(self, dir_name, **kwargs): 
        """
        Class initiator. Creates different instance variables depending on the keyword arguments given. 
        
        Every instance needs a directory name. If this object is reading information from output files dumped
        by PadeOps, then this is the directory where those files are stored. This object may also read information
        from a local subset of saved data.  
        
        
        The BudgetIO class will try to initialize from source files if kwarg `padeops=True` is given. 
        
        Alternatively, BudgetIO will read source files from PadeOps if the following keyword arguments are given: 
            runid (int)
            Lx (int)
            Ly (int)
            Lz (int)
            tidx (int) - optional, default is zero. 
            
        If not all of those keyword arguments are present, then the directory name will (attempt to) read from 
        budgets of saved .npz files. 

        Regardless of the method of reading budget files, __init__ will initialize the following fields: 
        RUN INFORMATION: 
            filename, dir_name, 
        DOMAIN VARIABLES: 
            Lx, Ly, Lz, nx, ny, nz, dx, dy, dz, xLine, yLine, zLine, 
        TURBINE VARIABLES: 
            nTurb, 
        PHYSICS: 
            Re, Ro, 
        BUDGET VARIABLES: 
            last_tidx, last_n, 
        
        """
        
        # print statements? default False
        if 'verbose' in kwargs and kwargs['verbose']: 
            self.verbose = True 
        else: 
            self.verbose = False
        
        self.dir_name = dir_name
        
        # all files associated with this case will begin with <filename>
        if not 'filename' in kwargs: 
            # defaults to the directory name, split on non-word characters
            dir_list = re.split('\W+', dir_name)
            # pick the last non-empty string
            self.filename = next(s for s in reversed(dir_list) if s)
        else: 
            self.filename = kwargs['filename']

        self.filename_budgets = self.filename + '_budgets.npz'  # standardize this
        
        # ========== Associate files ==========

        # if we are given the required keywords, try to initialize from PadeOps source files
        self.associate_padeops = False
        self.associate_npz = False
        self.associate_nml = False
        self.associate_fields = False 
        self.associate_budgets = False
        self.associate_grid = False
        self.associate_field = False
        self.associate_turbines = False
        self.normalized_xyz = False

        if all(x in kwargs for x in ['runid', 'Lx', 'Ly', 'Lz']) or ('padeops' in kwargs): 
            try: 
                self._init_padeops(**kwargs)
                self.associate_padeops = True

                if self.verbose: 
                    print('Initialized BudgetIO at ' + dir_name + ' from PadeOps source files. ')

            except OSError as err: 
                print('Attempted to read PadeOps output files, but at least one was missing.')
                print(err)
                raise
        
        # if padeops wasn't specifically requested, try with npz files: 
        else: 
            self._init_npz(**kwargs)
            self.associate_npz = True

            if self.verbose: 
                print('Initialized BudgetIO at ' + dir_name + ' from .npz files. ')

        self.budget = {}  # empty dictionary
        self.budget_xy = {}  # empty dictionary

        if 'read_budgets' in kwargs: 
            # if read_budgets passed in as keyword argument, read budgets on initialization
            self.read_budgets(budget_terms=kwargs['read_budgets'])
    

    def _init_padeops(self, **kwargs): 
        """
        Initializes source files to be read from output files in PadeOps. 
        
        Raises OSError if source files cannot be read
        """

        # Begin with keyword arguments - initialize if these were passed in
        for kwarg in kwargs: 
            if kwarg in ['runid', 'Lx', 'Ly', 'Lz']: 
                self.__dict__[kwarg] = kwargs[kwarg]
                
        # do namelist stuff - this is necessary if Lx, Ly, ... etc were not passed in. 
        try: 
            self._read_inputfile(**kwargs)
            
        except FileNotFoundError as err:  # TODO fix this
            warnings.warn("_init_padeops(): Could not find input file. ")
            print(err)
            raise err
            
        except IndexError as err: 
            warnings.warn("_init_padeops(): Could not find input file. ")
            print(err)
            raise err
        
        # default time ID for initialization
        if not 'tidx' in kwargs: 
            self.tidx = 0  # tidx defualts to zero
        else: 
            self.tidx = kwargs['tidx']

        # READ TURBINES, only do this if usewindturbines = True
        if self.associate_nml and self.input_nml['windturbines']['usewindturbines']: 
            
            if self.verbose: 
                print('_init_padeops(): Initializing wind turbine array object')
                
            turb_dir = self.input_nml['windturbines']['turbinfodir']
            ADM_type = self.input_nml['windturbines']['adm_type']
            num_turbines = self.input_nml['windturbines']['num_turbines']
            self.turbineArray = turbineArray.TurbineArray(turb_dir, 
                                                          ADM_type=ADM_type, 
                                                          num_turbines=num_turbines, 
                                                          verbose=self.verbose)
            self.associate_turbines = True
            
            if self.verbose: 
                print('_init_padeops(): Finished initializing wind turbine array with {:d} turbine(s)'.format(num_turbines))

        # These lines are taken almost verbatim from PadeOpsViz.py
        # read accompanying info file
        # may throw OSError
        
        # TODO this may all be in the input file. 
        
        if 'runid' not in self.__dict__:
            raise AttributeError("No RunID found. To explicitly pass one in, use kwarg: runid=")
        
        info_fname = self.dir_name + '/Run{:02d}_info_t{:06d}.out'.format(self.runid, self.tidx)
        # self.info = np.genfromtxt(info_fname, dtype=None)  
        # TODO: fix reading .info files
        self.time = 0  # self.info[0]
        
        if self.associate_nml: 
            self.nx = self.input_nml['input']['nx']
            self.ny = self.input_nml['input']['ny']
            self.nz = self.input_nml['input']['nz']
        else: 
            self.nx = int(self.info[1])
            self.ny = int(self.info[2])
            self.nz = int(self.info[3])
    
        # loads the grid, normalizes if associate_turbines = True (Should be done AFTER loading turbines)
        if not self.associate_grid: 
            self._load_grid(**kwargs)
            
        # object is reading from PadeOps output files directly
        if self.verbose: 
            print('BudgetIO initialized using info files at time:' + '{:.06f}'.format(self.time))
        
        self.field = {}
        
        try: 
            self.last_tidx = self.unique_tidx(return_last=True)  # last tidx in the run with fields
            
            # The following are initialized as the final saved instanteous field and budget: 
            self.field_tidx = self.last_tidx
            self.associate_fields = True

        except ValueError as e: 
            warnings.warn("_init_padeops(): No field files found!")
            if self.verbose: 
                print("\tNo field files found!")
            
        try: 
            self.all_budget_tidx = self.unique_budget_tidx(return_last=False)
            self.associate_budgets = True
        except ValueError as e: 
            warnings.warn("_init_padeops(): No budget files found!")
            if self.verbose: 
                print("\tNo budget files found.")
        
        if self.associate_budgets: 
            self.last_n = self.last_budget_n()  # last tidx with an associated budget
            self.budget_tidx = self.unique_budget_tidx()  # but may be changed by the user
            self.budget_n = self.last_n
            
    
    def _read_inputfile(self, **kwargs): 
        """
        Reads the input file (Fortran 90 namelist) associated with the CFD simulation. 

        Dependencies: f90nml, see https://github.com/marshallward/f90nml 
        """
        
        if f90nml is None: 
            warnings.warn('_read_inputfile(): No namelist reader loaded. ')
            return

        inputfile_ls = []  
                
        if 'runid' in kwargs.keys(): 
            # if given a runid, search for an inputfile with Run{runid} in the name
            inputfile_ls = glob.glob(self.dir_name + os.sep + 'Run{:02d}*.dat'.format(kwargs['runid']))
                        
            if len(inputfile_ls) == 0: 
                pass
#                 if self.verbose: 
#                     warnings.warn('_read_inputfile(): runid {:d} requested, \
#                         but no matching inputfile found was found.'.format(kwargs['runid']))

                # try to search all input files '*.dat' for the proper run and match it
                for inputfile in glob.glob(self.dir_name + os.sep + '*.dat'): 
                    input_nml = f90nml.read(inputfile) 
                    if self.verbose: 
                        print('_read_inputfile(): trying inputfile', inputfile)
                        
                    try: 
                        tmp_runid = input_nml['IO']['runid']
                    except KeyError as e: 
                        if self.verbose: 
                            print('_read_inputfile(): no runid for', inputfile)
                        tmp_runid = None

                    if tmp_runid == kwargs['runid']: 
                        self.input_nml = input_nml
                        self._convenience_variables(**kwargs)  # make some variables in the metadata more accessible, also loads grid
                        self.associate_nml = True  # successfully loaded input file
                        
                        if self.verbose: 
                            print('_read_inputfile(): matched RunID with', inputfile)
                        return
                        
        # if no runid given, or if the previous search failed, search all files ending in '*.dat' 
        if len(inputfile_ls) == 0: 
            inputfile_ls = glob.glob(self.dir_name + os.sep + '*.dat')  # for now, just search this. # TODO improve later? 

        # there should only be one input file for each run 
        if len(inputfile_ls) > 1: 
            warnings.warn('_read_inputfile(): Multiple files ending in *.dat found')
            if self.verbose: 
                print("    Found the following files:", inputfile_ls)
        
        # if there are still no input files found, we've got a problem 
        # TODO: trim dir_name() to remove trailing spaces
        
        if self.verbose: 
            print("Reading namelist file from {}".format(inputfile_ls[0]))
            
        self.input_nml = f90nml.read(inputfile_ls[0])
        self._convenience_variables(**kwargs)  # make some variables in the metadata more accessible
        self.associate_nml = True  # successfully loaded input file
        
        # BUDGET VARIABLES: 
        
        # TODO - fix this in the metadata
        
#         self.last_tidx = self.unique_budget_tidx()
#         self.last_n = self.last_budget_n() 

        
    def _convenience_variables(self, **kwargs): 
        """
        Aside from reading in the Namelist, which has all of the metadata, also make some
        parameters more accessible. 
        
        Called by _read_inputfile() and by _init_npz()
        
        Special note: these are all lower case when reading from the dictionary or namelist! 
        """
                
        # RUN VARIABLES: 
        self.runid = self.input_nml['io']['runid']
        
        # DOMAIN VARIABLES: 
        term_name = ['nx', 'ny', 'nz', 'Lx', 'Ly', 'Lz'] 
        term_search = ['nx', 'ny', 'nz', 'lx', 'ly', 'lz']  # same ordering
        
        for tname, tsearch in zip(term_name, term_search): 
            self.__dict__[tname] = budgetkey.key_search_r(self.input_nml, tsearch)
        
        # TURBINE VARIABLES: 
        self.nTurb = self.input_nml['windturbines']['num_turbines']

        # PHYSICS: 
        if self.input_nml['physics']['isinviscid']:  # boolean
            self.Re = np.Inf
        else: 
            self.Re = self.input_nml['physics']['re']

        if self.input_nml['physics']['usecoriolis']: 
            self.Ro = self.input_nml['physics']['ro']
        else: 
            self.Ro = np.Inf
        if self.input_nml['physics']['usecoriolis']: 
            self.latitude = self.input_nml['physics']['latitude']

        if self.input_nml['physics']['isstratified']: 
            self.Fr = self.input_nml['physics']['fr']
        else: 
            self.Fr = np.Inf

        if self.input_nml['problem_input']['Tref']:
            self.Tref = self.input_nml['problem_input']['Tref']
        
        

    
    def _load_grid(self, **kwargs): 
        """
        Creates dx, dy, dz, and xLine, yLine, zLine variables. 
        
        Expects (self.)Lx, Ly, Lz, nx, ny, nz to exist. 
        """

        if self.associate_grid and self.verbose: 
            print("_load_grid(): Grid already exists. ")
            return
        
        self.dx = self.Lx/self.nx
        self.dy = self.Ly/self.ny
        self.dz = self.Lz/self.nz

        # initialize grid
        self.xLine = np.linspace(0,self.Lx-self.dx,self.nx)
        self.yLine = np.linspace(0,self.Ly-self.dy,self.ny)
        self.zLine = np.linspace(self.dz/2,self.Lz-(self.dz/2),self.nz)

        self.associate_grid = True

        if self.associate_turbines: 
            if self.turbineArray.num_turbines == 1: 
                if self.verbose: 
                    print('_load_grid: attempting to normalize the origin to the turbine')
                    
                if 'normalize_origin' in kwargs and not kwargs['normalize_origin']: 
                    print("One turbine found, but keeping domain coordinates")
                    
                else: 
                    if self.verbose: 
                        print("Reading 1 turbine, normalizing origin. To turn off, initialize with `normalize_origin=False`")
                    self.xLine -= self.turbineArray.xloc
                    self.yLine -= self.turbineArray.yloc
                    self.zLine -= self.turbineArray.zloc
                    
                    self.normalized_xyz = True


    def _init_npz(self, **kwargs): 
        """
        Initializes the BudgetIO object by attempting to read .npz files saved from a previous BudgetIO object 
        from write_npz(). 

        Expects target files: 
        One filename including "{filename}_budgets.npz"
        One filename including "_metadata.npz"
        """
        
        # check budget files
        
        budget_files = glob.glob(self.dir_name + os.sep + '*_budgets.npz')
        if len(budget_files) == 0: 
            warnings.warn("No associated budget files found")
        else: 
            self.associate_budgets = True
        
        # load metadata: expects a file named <filename>_metadata.npy
        
        filepath = self.dir_name + os.sep + self.filename + '_metadata.npy'
        try: 
            self.input_nml = np.load(filepath, allow_pickle=True).item()
        except FileNotFoundError as e: 
            print(e)
            return
        
        # attempt to load turbine file - need this before loading grid
        if 'turbineArray' in self.input_nml['auxiliary']: 
            self.turbineArray = turbineArray.TurbineArray(init_dict=self.input_nml['auxiliary']['turbineArray'])
            self.associate_turbines = True
        
        self._convenience_variables(**kwargs)
        self.associate_nml = True
        
        if not self.associate_grid: 
            self._load_grid(**kwargs)

        if self.verbose: 
            print('_init_npz(): BudgetIO initialized using .npz files.')


    def set_filename(self, filename): 
        """
        Changes the filename associated with this object. 

        Make sure filename_budgets is consistent with __init__()
        """
        self.filename = filename
        self.filename_budgets = filename + "_budgets.npz"
        
        
    def write_npz(self, write_dir=None, budget_terms='default', filename=None, overwrite=False, 
                  xlim=None, ylim=None, zlim=None): 
        """
        Saves budgets as .npz files. Each budget receives its own .npz file, with the fourth dimension representing
        the budget term number (minus one, because python indexing starts at zero). 
        
        Budgets are defined in e.g. PadeOps/src/incompressible/budget_time_avg.F90. See budget_key.py
        From a high level: 
            Budget 0: Mean quantities (1st and 2nd order)
            Budget 1: Momentum budget terms
            Budget 2: MKE budget terms
            Budget 3: TKE budget terms
            
        parameters 
        ----------
        write_dir (str) : location to write .npz files. Default: same directory as self.outputdir_name
        budget_terms : list of budget terms to be saved (see ._parse_budget_terms()). Alternatively, 
            use 'current' to save the budget terms that are currently loaded. 
        filename (str) : calls self.set_filename()
        overwrite (bool) : if true, will overwrite existing .npz files. 
        xlim, ylim, zlim : slice bounds, see BudgetIO.slice()  # TODO: SAVE X,Y,Z information of slices
        """
        
        if not self.associate_budgets: 
            warnings.warn('write_npz(): No budgets associated! ') 
            return 
        
        # declare directory to write to, default to the working directory
        if write_dir is None: 
            write_dir = self.dir_name
        
        if budget_terms=='current': 
            # these are the currently loaded budgets
            key_subset = self.budget.keys()
            
        else: 
            # need to parse budget_terms with the key
            key_subset = self._parse_budget_terms(budget_terms)

            # load budgets
            sl = self.slice(budget_terms=key_subset, xlim=xlim, ylim=ylim, zlim=zlim)

        # if `filename` is provided, change this in the object
        # importantly, this needs to be done AFTER reading budgets! 
        if filename is not None: 
            self.set_filename(filename)

        filepath = write_dir + os.sep + self.filename_budgets
        
        # don't unintentionally overwrite files... 
        write_arrs = False  # this variable doesn't actually do anything
        if not os.path.exists(filepath): 
            write_arrs = True

        elif overwrite: 
            warnings.warn("File already exists, overwriting... ")
            write_arrs = True

        else: 
            warnings.warn("Existing files found. Failed to write; try passing overwrite=True to override.")
            return

        save_arrs = {}
        for key in key_subset: 
#             save_arrs[key] = self.budget[key]  # only save the requested budgets

            # TODO: We might want to crop the domain of the budgets here to reduce filesize... do that here: 
            save_arrs[key] = sl[key]

        # write npz files! 
        if write_arrs: 
            np.savez(filepath, **save_arrs)
            
            self.write_metadata(write_dir)  # TODO, need to save slice axes
            
            if self.verbose: 
                print("write_npz: Successfully saved the following budgets: ", list(key_subset))
                print("at " + filepath)
        
        
    def write_metadata(self, write_dir, xlim=None, ylim=None, zlim=None): 
        """
        The saved budgets aren't useful on their own unless we also save some information like the mesh
        used in the simulation and some other information like the physical setup. That goes here. 
        
        #TODO - switch this to an npz file. 
        """
        
        if self.associate_padeops: 
            # create a copy of the input namelist
            
            meta = self.input_nml.todict().copy() 
            meta['auxiliary'] = {  # this is new stuff not in the namelist
                                 'last_n': self.last_n, 
                                 'last_tidx': self.last_tidx
                                 # add more things here
                                }
            
            if self.associate_turbines: 
                meta['auxiliary'].update({
                    'turbineArray': self.turbineArray.todict()
                })
        else: 
            meta = self.input_nml.todict().copy()  # copy is probably unnecessary
        
        filename = write_dir + os.sep + self.filename + '_metadata.npy'
        np.save(filename, meta)
        
        if self.verbose: 
            print('write_metadata(): metadata written to {}'.format(filename))
            
            
    def write_mat(self, write_dir=None, budget_terms='default', filename=None, overwrite=False, 
                  xlim=None, ylim=None, zlim=None): 
        """
        Saves budgets as .mat (MATLAB) files. This is lazy code copied from write_npz(). 
        
        Budgets are defined in e.g. PadeOps/src/incompressible/budget_time_avg.F90. See budget_key.py
        From a high level: 
            Budget 0: Mean quantities (1st and 2nd order)
            Budget 1: Momentum budget terms
            Budget 2: MKE budget terms
            Budget 3: TKE budget terms
            
        parameters 
        ----------
        write_dir (str) : location to write .npz files. Default: same directory as self.outputdir_name
        budget_terms : list of budget terms to be saved (see ._parse_budget_terms()). Alternatively, 
            use 'current' to save the budget terms that are currently loaded. 
        filename (str) : calls self.set_filename()
        overwrite (bool) : if true, will overwrite existing .npz files. 
        xlim, ylim, zlim : slice bounds, see BudgetIO.slice()  # TODO: SAVE X,Y,Z information of slices
        """
        
        if not self.associate_budgets: 
            warnings.warn('write_mat(): No budgets associated! ') 
            return 
        
        # declare directory to write to, default to the working directory
        if write_dir is None: 
            write_dir = self.dir_name
        
        if budget_terms=='current': 
            key_subset = self.budget.keys()

        else: 
            key_subset = self._parse_budget_terms(budget_terms)

        # load budgets
        sl = self.slice(budget_terms=key_subset, xlim=xlim, ylim=ylim, zlim=zlim)

        if filename is not None: 
            self.set_filename(filename)
        
        filepath = write_dir + os.sep + self.filename_budgets[:-4] + '.mat'
        
        # don't unintentionally overwrite files... 
        write_arrs = False  
        if not os.path.exists(filepath): 
            write_arrs = True

        elif overwrite: 
            warnings.warn("File already exists, overwriting... ")
            write_arrs = True

        else: 
            warnings.warn("Existing files found. Failed to write; try passing overwrite=True to override.")
            return

        save_arrs = {}
        for key in key_subset: 
            save_arrs[key] = sl[key]

        # write mat files! 
        if write_arrs: 
            if self.verbose: 
                print('write_mat(): attempting to save budgets to', filepath)
                
            savemat(filepath, save_arrs)
            
            # SAVE METADATA HERE: 
            save_vars = ['input_nml'] #, 'xLine', 'yLine', 'zLine']
            save_dict = {key: self.__dict__[key] for key in save_vars}
            for key in ['x', 'y', 'z']: 
                save_dict[key] = sl[key]  # copy over x,y,z from slice
                
            filepath_meta = os.path.join(write_dir, self.filename + '_metadata.mat')
            savemat(filepath_meta, save_dict)
            
            if self.verbose: 
                print("write_mat(): Successfully saved the following budgets: ", list(key_subset))
                print("at" + filepath)
                print("write_mat(): Successfully saved metadata at" + filepath_meta)
                print('with fields', save_dict.keys())
        

    def read_metadata(self): 
        """
        Reads the saved .npy written in write_metadata(). 
        """
        # a rather boring function...
        self._init_npz()
     
            
    def read_fields(self, field_terms=None, tidx=None): 
        """
        Reads fields from PadeOps output files into the self.field dictionary. 
        
        Parameters
        ----------
        field_terms (list) : list of field terms to read, must be be limited to: 
            'u', 'v', 'w', 'p', 'T'
        tidx (int) : reads fields from the specified time ID. Default: self.last_tidx 
        
        Returns
        -------
        None
        
        """
        
        if not self.associate_fields: 
           raise AttributeError("read_fields(): No fields linked. ")
        
        dict_match = {'u':'uVel', 
                      'v':'vVel', 
                      'w':'wVel', 
                      'p':'prss', 
                      'T':'potT',
                      'kSGS':'kSGS',
                      'nSGS':'nSGS'}  # add more? 
        
        # parse terms: 
        if field_terms is None:
            terms = dict_match.keys()
            
        else:
            terms = [t for t in field_terms if t in dict_match.keys()]
        
        # parse tidx
        if tidx is None: 
            tidx = self.last_tidx
            
        # update self.time and self.tidx: 
        self.tidx = tidx
        
        info_fname = self.dir_name + '/Run{:02d}_info_t{:06d}.out'.format(self.runid, self.tidx)
        self.info = np.genfromtxt(info_fname, dtype=None)
        self.time = self.info[0]
        
        # the following is very similar to PadeOpsViz.ReadVelocities()
        
        for term in terms:             
            fname = self.dir_name + '/Run{:02d}_{:s}_t{:06d}.out'.format(self.runid, dict_match[term], tidx)
            tmp = np.fromfile(fname, dtype=np.dtype(np.float64), count=-1)
            self.field[term] = tmp.reshape((self.nx,self.ny,self.nz), order='F')  # reshape into a 3D array
            
        self.associate_field = True
            
        print('BudgetIO loaded fields {:s} at time: {:.06f}'.format(str(list(terms)), self.time))


    def time_slices(self, field_terms=None, tidx=None, xlim=None, ylim=None, zlim=None): 
        """
        Reads instantaneous slices from fields from PadeOps output files into the self.field dictionary.
        Each dictionary item is an array where the three dimensions are 2 spatial dimensions and time.
        
        Parameters
        ----------
        field_terms (list or list of strings) : fields to use
        tidx (list) : reads fields from the specified list of time ID. Default: self.unique_tidx
        keys (fields in slice `sl`) : keys to slice into from the input slice `sl`
        tidx (int) : time ID to read budgets from, see read_budgets(). Default None
        xlim, ylim, zlim (tuple) : in physical domain coordinates, the slice limits. If an integer is given, then the 
            dimension of the slice will be reduced by one. If None is given (default), then the entire domain extent is sliced. 
        
        Returns
        -------
        slices (dict) : dictionary organized with all of the sliced fields, keyed by the budget name, and additional keys for
            the slice domain 'x', 'y', 'z', 't'
                
        
        """

        if tidx is None:
            tidx = self.unique_tidx()

        slices = {}  # build from empty dict

        # find slice indices
        xid, yid, zid = self.get_xids(x=xlim, y=ylim, z=zlim, return_none=True, return_slice=True)
        xLine = self.xLine
        yLine = self.yLine
        zLine = self.zLine    


        dict_match = {'u':'uVel', 
                      'v':'vVel', 
                      'w':'wVel', 
                      'p':'prss', 
                      'T':'potT',
                      'kSGS':'kSGS',
                      'nSGS':'nSGS'}

        
        # parse terms: 
        if field_terms is None:
            terms = dict_match.keys()
        else:
            terms = [t for t in field_terms if t in dict_match.keys()]

        # initialize temporary array 
        for key in terms:
            slices[key] = np.squeeze(np.zeros([xid.stop-xid.start,yid.stop-yid.start,zid.stop-zid.start,len(tidx)]))
                
        # iterate through list of time ids
        for i in range(len(tidx)):
            self.read_fields(tidx=tidx[i])

            # iterate through field terms
            for key in terms:
                slices[key][:,:,i] = np.squeeze(self.field[key][xid, yid, zid])

        # also save domain information
        slices['x'] = xLine[xid]
        slices['y'] = yLine[yid]
        slices['z'] = zLine[zid]
        slices['tidx'] = tidx
        
        return slices

    def clear_budgets(self): 
        """
        Clears any loaded budgets. 

        Returns
        -------
        keys (list) : list of cleared budgets. 
        """
        if not self.associate_budgets: 
            if self.verbose: 
                print('clear_budgets(): no budgets to clear. ')
            return
        
        loaded_keys = self.budget.keys()
        self.budget = {}  # empty dictionary
#         self.associate_budgets = False

        if self.verbose: 
            print('clear_budgets(): Cleared loaded budgets: {}'.format(loaded_keys))
        
        return loaded_keys

    
    def read_budgets(self, budget_terms='default', mmap=None, overwrite=False, tidx=None): 
        """
        Accompanying method to write_budgets. Reads budgets saved as .npz files 
        
        Arguments 
        ---------
        budget_terms : list of terms (see ._parse_budget_terms())
        mmap : default None. Sets the memory-map settings in numpy.load(). Expects None, 'r+', 'r', 'w+', 'c'
        overwrite (bool) : if True, re-loads budgets that have already been loaded. Default False; checks  
            existing budgets before loading new ones. 
        tidx (int) : if given, requests budget dumps at a specific time ID. Default None. This only affects
            reading from PadeOps output files; .npz are limited to one saved tidx. 
        """
        
        if not self.associate_budgets: 
           raise AttributeError("read_budgets(): No budgets linked. ")
        
        # we need to handle computed quantities differently... 
        if any(t in ['uwake', 'vwake', 'wwake'] for t in budget_terms): 
            self.calc_wake()
            
            if self.verbose: 
                print("read_budgets: Successfully loaded wake budgets. ")


        # parse budget_terms with the key
        key_subset = self._parse_budget_terms(budget_terms, include_wakes=False)
        
        if not overwrite:
            remove_keys = [key for key in key_subset if key in self.budget.keys()]
            if len(remove_keys) > 0 and self.verbose: 
                print("read_budgets(): requested budgets that have already been loaded. \
                    \n  Removed the following: {}. Pass overwrite=True to read budgets anyway.".format(remove_keys))
            
            # remove items that have already been loaded in  
            key_subset = {key:key_subset[key] for key in key_subset if key not in self.budget.keys()}

        if self.associate_padeops: 
            self._read_budgets_padeops(key_subset, tidx=tidx)  # this will not include wake budgets
        
        elif self.associate_npz: 
            self._read_budgets_npz(key_subset, mmap=mmap)
        
        if self.verbose and len(key_subset) > 0: 
            print("read_budgets: Successfully loaded budgets. ")

    def read_budgets_xy(self, budget_terms='default', mmap=None, overwrite=False, tidx=None): 
        """
        Accompanying method to write_budgets. Reads budgets saved as .npz files 
        
        Arguments 
        ---------
        budget_terms : list of terms (see ._parse_budget_terms())
        mmap : default None. Sets the memory-map settings in numpy.load(). Expects None, 'r+', 'r', 'w+', 'c'
        overwrite (bool) : if True, re-loads budgets that have already been loaded. Default False; checks  
            existing budgets before loading new ones. 
        tidx (int) : if given, requests budget dumps at a specific time ID. Default None. This only affects
            reading from PadeOps output files; .npz are limited to one saved tidx. 
        """

        if not self.associate_budgets: 
           raise AttributeError("read_budgets(): No budgets linked. ")
        
        # we need to handle computed quantities differently... 
        if any(t in ['uwake', 'vwake', 'wwake'] for t in budget_terms): 
            self.calc_wake()
            
            if self.verbose: 
                print("read_budgets: Successfully loaded wake budgets. ")


        # parse budget_terms with the key
        key_subset = self._parse_budget_xy_terms(budget_terms, include_wakes=False)
        
        if not overwrite:
            remove_keys = [key for key in key_subset if key in self.budget.keys()]
            if len(remove_keys) > 0 and self.verbose: 
                print("read_budgets(): requested budgets that have already been loaded. \
                    \n  Removed the following: {}. Pass overwrite=True to read budgets anyway.".format(remove_keys))
            
            # remove items that have already been loaded in  
            key_subset = {key:key_subset[key] for key in key_subset if key not in self.budget.keys()}

        if self.associate_padeops: 
            self._read_budgets_xy_padeops(key_subset, tidx=tidx)  # this will not include wake budgets
        
        elif self.associate_npz: 
            self._read_budgets_npz(key_subset, mmap=mmap)
        
        if self.verbose and len(key_subset) > 0: 
            print("read_budgets: Successfully loaded budgets. ")
        

    def _read_budgets_padeops(self, key_subset, tidx): 
        """
        Uses a method similar to ReadVelocities_Budget() in PadeOpsViz to read and store full-field budget terms. 
        """

        budget4_components = [11, 22, 33, 13, 23]
        
        if tidx is None: 
            tidx = self.budget_tidx
            
        elif tidx not in self.all_budget_tidx: 
            # find the nearest that actually exists
            tidx_arr = np.array(self.all_budget_tidx)
            closest_tidx = tidx_arr[np.argmin(np.abs(tidx_arr-tidx))]
            
            print("Requested budget tidx={:d} could not be found. Using tidx={:d} instead.".format(tidx, closest_tidx))
            tidx = closest_tidx 
            
        print(tidx)
        # update self.time and self.tidx: 
#         self.tidx = tidx
        
#         info_fname = self.dir_name + '/Run{:02d}_info_t{:06d}.out'.format(self.runid, self.tidx)
#         self.info = np.genfromtxt(info_fname, dtype=None)
#         self.time = self.info[0]

        # these lines are almost verbatim from PadeOpsViz.py
        for key in key_subset:
            budget, term = BudgetIO.key[key]
            if budget==4:
                component = budget4_components[int(np.floor((term-1)/10))]

                if term > 10:
                    term = term % 10
                    if term == 0:
                        term = 10
                searchstr =  self.dir_name + '/Run{:02d}_budget{:01d}_{:02d}_term{:02d}_t{:06d}_*.s3D'.format(self.runid, budget, component, term, tidx)
                u_fname = glob.glob(searchstr)[0]  

            else:
                searchstr =  self.dir_name + '/Run{:02d}_budget{:01d}_term{:02d}_t{:06d}_*.s3D'.format(self.runid, budget, term, tidx)
                u_fname = glob.glob(searchstr)[0]  
            
            self.budget_n = int(re.findall('.*_t\d+_n(\d+)', u_fname)[0])  # extract n from string
            self.budget_tidx = tidx
            
            temp = np.fromfile(u_fname, dtype=np.dtype(np.float64), count=-1)
            self.budget[key] = temp.reshape((self.nx,self.ny,self.nz), order='F')  # reshape into a 3D array

        if self.verbose and len(key_subset) > 0: 
            print('PadeOpsViz loaded the budget fields at time:' + '{:.06f}'.format(tidx))

    def _read_budgets_xy_padeops(self, key_subset, tidx): 
        """
        Uses a method similar to ReadVelocities_Budget() in PadeOpsViz to read and store full-field budget terms. 
        """
        
        if tidx is None: 
            tidx = self.budget_tidx
            
        elif tidx not in self.all_budget_tidx: 
            # find the nearest that actually exists
            tidx_arr = np.array(self.all_budget_tidx)
            closest_tidx = tidx_arr[np.argmin(np.abs(tidx_arr-tidx))]
            
            print("Requested budget tidx={:d} could not be found. Using tidx={:d} instead.".format(tidx, closest_tidx))
            tidx = closest_tidx 
            
        # update self.time and self.tidx: 
        #         self.tidx = tidx
        
        #         info_fname = self.dir_name + '/Run{:02d}_info_t{:06d}.out'.format(self.runid, self.tidx)
        #         self.info = np.genfromtxt(info_fname, dtype=None)
        #         self.time = self.info[0]

        for budget in key_subset:
            if budget != 4:
                searchstr =  self.dir_name + '/Run{:02d}_budget{:01d}_t{:06d}_*.stt'.format(self.runid, budget, tidx)
                u_fname = glob.glob(searchstr)[0]  
            
                self.budget_n = int(re.findall('.*_t\d+_n(\d+)', u_fname)[0])  # extract n from string
                self.budget_tidx = tidx
            
                temp = np.genfromtxt(u_fname, dtype=np.dtype(np.float64))
                for key in key_subset[budget]:
                    self.budget_xy[key] = temp[:,BudgetIO.key_xy[key][1]-1]  # reshape into a 3D array
            elif budget == 4:
                for component in key_subset[budget]:
                    searchstr =  self.dir_name + '/Run{:02d}_budget{:01d}_{:02d}_t{:06d}*.stt'.format(self.runid, budget, component, tidx)
                    u_fname = glob.glob(searchstr)[0]  
            
                    self.budget_n = int(re.findall('.*_t\d+_n(\d+)', u_fname)[0])  # extract n from string
                    self.budget_tidx = tidx
            
                    temp = np.genfromtxt(u_fname, dtype=np.dtype(np.float64))
                    
                    for key in key_subset[budget][component]:
                        self.budget_xy[key] = temp[:,key_subset[budget][component][key][1]-1]  
                    
        if self.verbose and len(key_subset) > 0: 
            print('PadeOpsViz loaded the budget fields at time:' + '{:.06f}'.format(tidx))


    def _read_budgets_npz(self, key_subset, mmap=None): 
        """
        Reads budgets written by .write_npz() and loads them into memory
        """

        # load the npz file and keep the requested budget keys
        for key in key_subset: 
            npz = np.load(self.dir_name + os.sep + self.filename_budgets)
            self.budget[key] = npz[key]  

        if self.verbose: 
            print('PadeOpsViz loaded the following budgets from .npz: ', list(key_subset.keys()))


    def _parse_budget_terms(self, budget_terms, include_wakes=False): 
        """
        Takes a list of budget terms, either keyed in index form (budget #, term #) or in common form (e.g. ['u_bar', 'v_bar'])
        and returns a subset of the `keys` dictionary that matches two together. `keys` dictionary is always keyed in common form. 

        budget_terms can also be a string: 'all', or 'default'. 

        'default' tries to load the following: 
            Budget 0 terms: ubar, vbar, wbar, all Reynolds stresses, and p_bar
            Budget 1 terms: all momentum terms
        'all' checks what budgets exist and tries to load them all. 

        For more information on the bi-directional keys, see budget_key.py
        
        Arguments
        ---------
        budget_terms : list of strings or string, see above
        include_wakes (bool) : optional, includes wake budgets if True, default False. 
        """

        # add string shortcuts here... # TODO move shortcuts to budgetkey.py? 
        if budget_terms=='default': 
            budget_terms = ['ubar', 'vbar', 'wbar', 
                            'tau11', 'tau12', 'tau13', 'tau22', 'tau23', 'tau33', 
                            'pbar']

        elif budget_terms=='all': 
            budget_terms = self.existing_terms(include_wakes=include_wakes)
            
        elif budget_terms=='RANS': 
            budget_terms = ['ubar', 'vbar', 'wbar', 'pbar', 'uu', 'uv', 'uw', 'vv', 'vw', 'ww', 'dpdx', 'dpdy', 'dpdz']

        elif type(budget_terms)==str: 
            warnings.warn("keyword argument budget_terms must be either 'default', 'all', 'RANS' or a list.")
            return {}  # empty dictionary

        # parse through terms: they are either 1) valid, 2) missing (but valid keys), or 3) invalid (not in BudgetIO.key)

        existing_keys = self.existing_terms(include_wakes=include_wakes)
        existing_tup = [BudgetIO.key[key] for key in existing_keys]  # corresponding associated tuples (#, #)

        valid_keys = [t for t in budget_terms if t in existing_keys]
        missing_keys = [t for t in budget_terms if t not in existing_keys and t in BudgetIO.key]
        invalid_terms = [t for t in budget_terms if t not in BudgetIO.key and t not in BudgetIO.key.inverse]
        
        valid_tup = [tup for tup in budget_terms if tup in existing_tup]  # existing tuples
        missing_tup = [tup for tup in budget_terms if tup not in existing_tup and tup in BudgetIO.key.inverse]
        
        # now combine existing valid keys and valid tuples, removing any duplicates

        valid_terms = set(valid_keys + [BudgetIO.key.inverse[tup][0] for tup in valid_tup])  # combine and remove duplicates
        missing_terms = set(missing_keys + [BudgetIO.key.inverse[tup][0] for tup in missing_tup])

        # generate the key
        key_subset = {key: BudgetIO.key[key] for key in valid_terms}
        
        # warn the user if some requested terms did not exist
        if len(key_subset) == 0: 
            warnings.warn('_parse_budget_terms(): No keys being returned; none matched.')

        if len(missing_terms) > 0: 
            warnings.warn('_parse_budget_terms(): Several terms were requested but the following could not be found: \
                {}'.format(missing_terms))

        if len(invalid_terms) > 0: 
            warnings.warn('_parse_budget_terms(): The following budget terms were requested but the following do not exist: \
                {}'.format(invalid_terms))

        # TODO - fix warning messages for the wakes

        return key_subset

    def _parse_budget_xy_terms(self, budget_terms, include_wakes=False): 
        """
        Takes a list of budget terms, either keyed in index form (budget #, term #) or in common form (e.g. ['u_bar', 'v_bar'])
        and returns a subset of the `keys` dictionary that matches two together. `keys` dictionary is always keyed in common form. 

        budget_terms can also be a string: 'all', or 'default'. 

        'default' tries to load the following: 
            Budget 0 terms: ubar, vbar, wbar, all Reynolds stresses, and p_bar
            Budget 1 terms: all momentum terms
        'all' checks what budgets exist and tries to load them all. 

        For more information on the bi-directional keys, see budget_key.py
        
        Arguments
        ---------
        budget_terms : list of strings or string, see above
        include_wakes (bool) : optional, includes wake budgets if True, default False. 
        """

        # add string shortcuts here... # TODO move shortcuts to budgetkey.py? 
        if budget_terms=='default': 
            budget_terms = ['ubar', 'vbar', 'wbar', 
                            'tau11', 'tau12', 'tau13', 'tau22', 'tau23', 'tau33', 
                            'pbar']

        elif budget_terms=='all': 
            budget_terms = self.existing_terms_xy(include_wakes=include_wakes)
            
        elif budget_terms=='RANS': 
            budget_terms = ['ubar', 'vbar', 'wbar', 'pbar', 'uu', 'uv', 'uw', 'vv', 'vw', 'ww', 'dpdx', 'dpdy', 'dpdz']

        elif type(budget_terms)==str: 
            warnings.warn("keyword argument budget_terms must be either 'default', 'all', 'RANS' or a list.")
            return {}  # empty dictionary

        # parse through terms: they are either 1) valid, 2) missing (but valid keys), or 3) invalid (not in BudgetIO.key)

        existing_keys = self.existing_terms_xy(include_wakes=include_wakes)
        existing_tup = [BudgetIO.key_xy[key] for key in existing_keys]  # corresponding associated tuples (#, #)
        
        valid_keys = [t for t in budget_terms if t in existing_keys]
        missing_keys = [t for t in budget_terms if t not in existing_keys and t in BudgetIO.key_xy]
        invalid_terms = [t for t in budget_terms if t not in BudgetIO.key_xy and t not in BudgetIO.key_xy.inverse]

        valid_tup = [BudgetIO.key_xy[key] for key in budget_terms if BudgetIO.key_xy[key] in existing_tup]  # existing tuples
        missing_tup = [BudgetIO.key_xy[key] for key in budget_terms if BudgetIO.key_xy[key] not in existing_tup and BudgetIO.key_xy[key] in BudgetIO.key_xy.inverse]

        # now combine existing valid keys and valid tuples, removing any duplicates

        valid_terms = set(valid_keys + [BudgetIO.key_xy.inverse[tup][0] for tup in valid_tup])  # combine and remove duplicates
        missing_terms = set(missing_keys + [BudgetIO.key_xy.inverse[tup][0] for tup in missing_tup])

        # find the budgets
        budgets = []
        [budgets.append(tup[0]) for tup in valid_tup if tup not in budgets]

        # initialize key_subset nested dictionary
        # for the xy budgets the key_subset dict is of the form:
        # key_subset = {budget: {term: tuple}}
        # if budget 4 is included then there is further nesting for the component eg:
        # key_subset = {4: {component : {term: tuple}}}
        if 4 in budgets:
            key_subset = {4:{}}
        else:
            key_subset = {}
        
        for budget in budgets:
            tmp_dict = {}
            if budget != 4:
                for key in valid_terms:
                    if BudgetIO.key_xy[key][0] == budget and budget != 4:
                        tmp_dict[key] = BudgetIO.key_xy[key]
                key_subset[budget] = tmp_dict
        
            elif budget == 4:
                tmp_uu_dict = {}
                tmp_uw_dict = {}
                tmp_vw_dict = {}
                tmp_ww_dict = {}
        
                for key in valid_terms:
                    if BudgetIO.key_xy[key][0] == 4 and BudgetIO.key_xy[key][1] >= 1 and BudgetIO.key_xy[key][1] <= 9:
                        tmp_uu_dict[key] = BudgetIO.key_xy[key]
                    elif BudgetIO.key_xy[key][0] == 4 and BudgetIO.key_xy[key][1] >= 10 and BudgetIO.key_xy[key][1] <= 18:
                        tmp_uw_dict[key] = [4, BudgetIO.key_xy[key][1] - 9]
                    elif BudgetIO.key_xy[key][0] == 4 and BudgetIO.key_xy[key][1] >= 19 and BudgetIO.key_xy[key][1] <= 27:
                        tmp_vw_dict[key] = [4, BudgetIO.key_xy[key][1] - 18]
                    elif BudgetIO.key_xy[key][0] == 4 and BudgetIO.key_xy[key][1] >= 28 and BudgetIO.key_xy[key][1] <= 36:
                        tmp_ww_dict[key] = [4, BudgetIO.key_xy[key][1] - 27]
                
                key_subset[budget][11] = tmp_uu_dict
                key_subset[budget][13] = tmp_uw_dict
                key_subset[budget][23] = tmp_vw_dict
                key_subset[budget][33] = tmp_ww_dict

        '''
        # warn the user if some requested terms did not exist
        if len(key_subset) == 0: 
            warnings.warn('_parse_budget_terms(): No keys being returned; none matched.')

        if len(missing_terms) > 0: 
            warnings.warn('_parse_budget_terms(): Several terms were requested but the following could not be found: \
                {}'.format(missing_terms))

        if len(invalid_terms) > 0: 
            warnings.warn('_parse_budget_terms(): The following budget terms were requested but the following do not exist: \
                {}'.format(invalid_terms))
        '''
        # TODO - fix warning messages for the wakes

        return key_subset 


    def _get_inflow(self, offline=False, wInflow=False): 
        """
        Calls the appropriate functions in inflow.py to retrieve the inflow profile for the corresponding flow. 

        Arguments
        ---------
        offline (bool) : if True, uses the target inflow profile prescribed by initialize.F90. Default (False) reads
            the inflow profile from the first x-index of the domain and average over y. 
        wInflow (bool) : if True, returns an array of w-inflow velocities. Default (False) only returns u, v. 
        
        Returns
        -------
        u (array) : [nz x 1] array of u-velocities as a function of height
        v (array) : [nz x 1] array of v-velocities as a function of height
        w (array) : [nz x 1] array of w-velocities as a function of height. Nominally this is all zero. 

        """
        
        # load using InflowParser: 

        if offline: 
            if self.associate_nml: 
                u, v = inflow.InflowParser.inflow_offline(**dict(self.input_nml['AD_coriolisinput']), zLine=self.zLine)
            
            # reading from the budgets
            else: 
                warnings.warn('_get_inflow: Requested offline inflow, but namelist not associated. Trying online.')
                u, v = inflow.InflowParser.inflow_budgets(self)
                
        else: 
            u, v = inflow.InflowParser.inflow_budgets(self) 

        # return requested information
 
        if wInflow: 
            # If this is not nominally zero, then this will need to be fixed. 
            w = np.zeros(self.zLine.shape)
            return np.array([u, v, w])

        else: 
            return np.array([u, v])

    
    def calc_wake(self, offline=False, wInflow=False, overwrite=False):
        """
        Computes the wake deficit by subtracting the target inflow from the flow field. 
        # TODO - right now this must compute at minimum uwake and vwake. Fix?  

        Arguments
        ---------
        see _get_inflow()
        overwrite (bool) : re-computes wakes if they are already read in. 

        Returns
        -------
        None (updates self.budget[] with 'uwake' and 'vwake' keys)

        """ 

        target_terms = ['uwake', 'vwake']
        req_terms = ['ubar', 'vbar']  # required budget terms
        if wInflow: 
            target_terms.append('wwake')
            req_terms.append('wbar') # might also need w

        # check to see if the terms exist already
        if all(t in self.budget.keys() for t in target_terms): 
            if not overwrite: 
                warnings.warn('Wake terms already computed, returning. To compute anyway, use keyword overwrite=True')
                return
        
        # Need mean velocity fields to be loaded
        if not all(t in self.budget.keys() for t in req_terms): 
            self.read_budgets(budget_terms=req_terms)

        # retrieve inflow
        if wInflow: 
            u, v, w = self._get_inflow(offline=offline, wInflow=wInflow)
            self.budget['wwake'] = self.budget['wbar'] - w[np.newaxis, np.newaxis, :]
        
        else: 
            u, v = self._get_inflow(offline=offline, wInflow=wInflow)

        self.budget['uwake'] = u[np.newaxis, np.newaxis, :] - self.budget['ubar']
        self.budget['vwake'] = v[np.newaxis, np.newaxis, :] - self.budget['vbar']

        if self.verbose: 
            print("calc_wake(): Computed wake velocities. ")

    
    def slice(self, budget_terms=None, field=None, sl=None, keys=None, tidx=None, 
              xlim=None, ylim=None, zlim=None, 
              overwrite=False, round_extent=True): 
        """
        Returns a slice of the requested budget term(s) as a dictionary. 

        Arguments
        ---------
        budget_terms (list or string) : budget term or terms to slice from. If None, expects a value for `field`
        field (arraylike or dict of arraylike) : fields similar to self.budget[]
        sl (slice from self.slice()) : dictionary of fields to be sliced into again. 
        keys (fields in slice `sl`) : keys to slice into from the input slice `sl`
        tidx (int) : time ID to read budgets from, see read_budgets(). Default None
        xlim, ylim, zlim (tuple) : in physical domain coordinates, the slice limits. If an integer is given, then the 
            dimension of the slice will be reduced by one. If None is given (default), then the entire domain extent is sliced. 
        overwrite (bool) : Overwrites loaded budgets, see read_budgets(). Default True
        round_extent (bool) : Rounds extents to the nearest integer. Default True
        
        Returns
        -------
        slices (dict) : dictionary organized with all of the sliced fields, keyed by the budget name, and additional keys for
            the slice domain 'x', 'y', and 'z'
        
        """

        if sl is None: 
            xid, yid, zid = self.get_xids(x=xlim, y=ylim, z=zlim, return_none=True, return_slice=True)
            xLine = self.xLine
            yLine = self.yLine
            zLine = self.zLine
        else: 
            xid, yid, zid = self.get_xids(x=xlim, y=ylim, z=zlim, 
                                          x_ax=sl['x'], y_ax=sl['y'], z_ax=sl['z'], 
                                          return_none=True, return_slice=True)
            xLine = sl['x']
            yLine = sl['y']
            zLine = sl['z']

        slices = {}  # build from empty dict

        if field is not None: 
            # slice the given field
            
            if type(field) == dict: 
                # iterate through dictionary of fields
                for key in field.keys(): 
                    slices[key] = np.squeeze(field[key][xid, yid, zid])
            else: 
                slices['field'] = np.squeeze(field[xid, yid, zid])

        elif budget_terms is not None: 
            # read budgets
            key_subset = self._parse_budget_terms(budget_terms)
            self.read_budgets(budget_terms=key_subset, tidx=tidx, overwrite=overwrite)

            for term in key_subset: 
                slices[term] = np.squeeze(self.budget[term][xid, yid, zid])  

        elif keys is not None and sl is not None: 
            # slice into slices

            if type(keys) == list: 
                for key in keys: 
                    slices[key] = np.squeeze(sl[key][xid, yid, zid])

            else: 
                slices[key] = np.squeeze(sl[xid, yid, zid])
                
        else: 
            warnings.warn("BudgetIO.slice(): either budget_terms= or field= must be initialized.")
            return None
        
        # also save domain information
        slices['x'] = xLine[xid]
        slices['y'] = yLine[yid]
        slices['z'] = zLine[zid]
        
        # build and save the extents, either in 1D, 2D, or 3D
        ext = []
        for term in ['x', 'y', 'z']: 
            if len(slices[term]) > 1:  # if this is actually a slice (not a number), then add it to the extents
                ext += [np.min(slices[term]), np.max(slices[term])]
        
        if round_extent: 
            slices['extent'] = np.array(ext).round()
        else: 
            slices['extent'] = np.array(ext)

        return slices


    def get_xids(self, x=None, y=None, z=None, 
                 x_ax=None, y_ax=None, z_ax=None, 
                 return_none=False, return_slice=False): 
        """
        Translates x, y, and z limits in the physical domain to indices based on self.xLine, self.yLine, and self.zLine

        Arguments
        ---------
        x, y, z : float or iterable (tuple, list, etc.) of physical locations to return the nearest index for
        return_none : if True, populates output tuple with None if input is None. Default False. 
        return_slice : if True, returns a tuple of slices instead a tuple of lists. 

        Returns
        -------
        xid, yid, zid : list or tuple of lists with indices for the requested x, y, z, args in the order: x, y, z. 
            If, for example, y and z are requested, then the returned tuple will have (yid, zid) lists. 
            If only one value (float or int) is passed in for e.g. x, then an integer will be passed back in xid. 
        """

        if not self.associate_grid: 
            raise(AttributeError('No grid associated. '))

        # set up this way in case we want to introduce an offset later on (i.e. turbine-centered coordinates)
        if x_ax is None: 
            x_ax = self.xLine  # - offset_x  # TODO? 
        if y_ax is None: 
            y_ax = self.yLine
        if z_ax is None:        
            z_ax = self.zLine

        ret = ()

        # iterate through x, y, z, index matching for each term
        for s, s_ax in zip([x, y, z], [x_ax, y_ax, z_ax]): 
            if s is not None: 
                if hasattr(s, '__iter__'): 
                    xids = [np.argmin(np.abs(s_ax-xval)) for xval in s]
                else: 
                    xids = np.argmin(np.abs(s_ax-s))
                
                if return_slice:  # append slices to the return tuple
                    ret = ret + (slice(np.min(xids), np.max(xids)+1), )
                    
                else:  # append index list to the return tuple
                    ret = ret + (xids, )
            
            elif return_none:  # fill with None or slice(None)
                if return_slice: 
                    ret = ret + (slice(None), )
                    
                else: 
                    ret = ret + (None, )
        
        if len(ret)==1: 
            return ret[0]  # don't return a length one tuple 
        else: 
            return ret
    

    def unique_tidx(self, return_last=False): 
        """
        Pulls all the unique tidx values from a directory. 
        
        Arguments
        ---------
        return_last (bool) : If True, returns only the largest value of TIDX. Default False. 
        """

        if not self.associate_padeops:
            pass  # TODO - is this lost information? is it useful information? 
            
        # retrieves filenames and parses unique integers, returns an array of unique integers
        filenames = os.listdir(self.dir_name)
        runid = self.runid
        
        # searches for the formatting *_t(\d+)* in all filenames
        t_list = [int(re.findall('Run{:02d}.*_t(\d+).*.out'.format(runid), name)[0]) 
                  for name in filenames 
                  if re.findall('Run{:02d}.*_t(\d+).*.out'.format(runid), name)]
        
        if return_last: 
            return np.max(t_list)
        else: 
            return np.unique(t_list)

    
    def unique_budget_tidx(self, return_last=True): 
        """
        Pulls all the unique tidx values from a directory. 
        
        Parameters
        ----------
        return_last (bool) : If False, returns only the largest TIDX associated with budgets. 
            Else, returns an entire list of unique tidx associated with budgets. Default True
        """

        # TODO: fix for .npz

        # retrieves filenames and parses unique integers, returns an array of unique integers
        filenames = os.listdir(self.dir_name)
        runid = self.runid
        
        # searches for the formatting *_t(\d+)* in budget filenames
        t_list = [int(re.findall('Run{:02d}.*budget.*_t(\d+).*'.format(runid), name)[0]) 
                  for name in filenames 
                  if re.findall('Run{:02d}.*budget.*_t(\d+).*'.format(runid), name)]
        
        if len(t_list) == 0: 
            raise ValueError("No budget times found. ")
        
        if return_last: 
            return np.max(t_list)
        else: 
            return np.unique(t_list)
        
        
    def unique_times(self, return_last=False): 
        """
        Reads the .out file of each unique time and returns an array of [physical] times corresponding
        to the time IDs from unique_tidx(). 
        
        Returns
        -------
        times (arr) : list of times associated with each time ID in unique_tidx()
        return_false (bool) : if True, only returns the final element in the array. Default: False
        """
        
        times = []; 
        
        if return_last:  # save time by only reading the final TIDX
            tidx = self.unique_tidx(return_last=return_last)
            fname = os.path.join(self.dir_name, "Run{:02d}_info_t{:06d}.out".format(self.runid, tidx))
            t = np.genfromtxt(fname, dtype=None)[0]
            return t
        
        for tidx in self.unique_tidx(): 
            fname = os.path.join(self.dir_name, "Run{:02d}_info_t{:06d}.out".format(self.runid, tidx))
            t = np.genfromtxt(fname, dtype=None)[0]
            times.append(t)

        return np.array(times)

    
    def last_budget_n(self): 
        """
        Pulls all unique n from budget terms in a directory and returns the largest value. 
        """

        # TODO: fix for .npz

        filenames = os.listdir(self.dir_name)
        runid = self.runid

        # capturing *_n(\d+)* in filenames
        t_list = [int(re.findall('Run{:02d}.*_n(\d+).*'.format(runid), name)[0]) 
                  for name in filenames 
                  if re.findall('Run{:02d}.*_n(\d+).*'.format(runid), name)]
        return np.max(t_list)
    
    
    def existing_budgets(self): 
        """
        Checks file names for which budgets were output.  
        """
        filenames = os.listdir(self.dir_name)

        if self.associate_npz: 
            filename = self.dir_name + os.sep + self.filename_budgets
            with np.load(filename) as npz: 
                t_list = npz.files  # load all the budget filenames
            
            budget_list = [BudgetIO.key[t][0] for t in t_list]

        elif self.associate_padeops: 
            runid = self.runid
            # capturing *_budget(\d+)* in filenames
            budget_list = [int(re.findall('Run{:02d}.*_budget(\d+).*'.format(runid), name)[0]) 
                           for name in filenames 
                           if re.findall('Run{:02d}.*_budget(\d+).*'.format(runid), name)]

        else: 
            warnings.warn('existing_budgets(): No associated budget files found. ')
            return []
        
        if 0 in budget_list: 
            budget_list.append(5)  # wake budgets can be recovered from mean budgets

        return list(np.unique(budget_list))
    
    
    def existing_terms(self, budget=None, include_wakes=True): 
        """
        Checks file names for a particular budget and returns a list of all the existing terms.  

        Arguments 
        ---------
        budget (integer) : optional, default None. If provided, searches a particular budget for existing terms. 
            Otherwise, will search for all existing terms. `budget` can also be a list of integers. 
            Budget 0: mean statistics
            Budget 1: momentum
            Budget 2: MKE
            Budget 3: TKE
            Budget 5: Wake deficit
        include_wakes (bool) : Includes wakes in the returned budget terms if True, default True. 

        Returns
        -------
        t_list (list) : list of tuples of budgets found

        """

        t_list = []

        budget4_comp_dict = {11 : 0,
                             22 : 10,
                             33 : 20,
                             13 : 30, 
                             23 : 40}
        
        
        # if no budget is given, look through all saved budgets
        if budget is None: 
            budget_list = self.existing_budgets()
        
        else: 
            # convert to list if integer is given
            if type(budget) != list: 
                budget_list = [budget]
            else: 
                budget_list = budget

        # find budgets matching .npz convention in write_npz()
        if self.associate_npz: 
            filename = self.dir_name + os.sep + self.filename_budgets
            with np.load(filename) as npz: 
                all_terms = npz.files

            if budget is None:  # i.e. requesting all budgets
                return all_terms  # we can stop here without sorting through each budget
            
            tup_list = [BudgetIO.key[t] for t in all_terms]  # list of associated tuples
            t_list = []  # this is the list to be built and returned

            for b in budget_list: 
                t_list += [tup for tup in tup_list if tup[0] == b]

        # find budgets by name matching with PadeOps output conventions
        elif self.associate_padeops: 

            filenames = os.listdir(self.dir_name)
            runid = self.runid
            
            tup_list = []
            # loop through budgets
            for b in budget_list: 
                # capturing *_term(\d+)* in filenames
                terms = [int(re.findall('Run{:02d}_budget{:01d}_term(\d+).*'.format(runid, b), name)[0]) 
                        for name in filenames if 
                        re.findall('Run{:02d}_budget{:01d}_term(\d+).*'.format(runid, b), name)]
                tup_list += [((b, term)) for term in set(terms)]  # these are all tuples

                # reynolds stress budgets
                if b == 4:
                    for component in budget4_comp_dict:
                        terms = [int(re.findall('Run{:02d}_budget{:01d}_{:01d}_term(\d+).*'.format(runid, b, component), name)[0]) + budget4_comp_dict[component]
                                 for name in filenames if 
                                 re.findall('Run{:02d}_budget{:01d}_{:01d}_term(\d+).*'.format(runid, b, component), name)]
                        tup_list += [((b, term)) for term in set(terms)]  # these are all tuples
                    
                
                # wake budgets: 
                wake_budgets = (1, 2, 3)
                if include_wakes and b == 5:  
                    terms = [int(re.findall('Run{:02d}_budget{:01d}_term(\d+).*'.format(runid, 0), name)[0]) 
                            for name in filenames if 
                            re.findall('Run{:02d}_budget{:01d}_term(\d+).*'.format(runid, 0), name)]  # read from mean budgets

                    tup_list += [((b, term)) for term in wake_budgets if term in terms]
            
            # convert tuples to keys
            t_list = [BudgetIO.key.inverse[key][0] for key in tup_list]
        
        else: 
            warnings.warn('existing_terms(): No terms found for budget ' + str(budget))
            return []
        
        return t_list

    def existing_terms_xy(self, budget=None, include_wakes=True): 
        """
        Checks file names for a particular budget and returns a list of all the existing terms.  

        Arguments 
        ---------
        budget (integer) : optional, default None. If provided, searches a particular budget for existing terms. 
            Otherwise, will search for all existing terms. `budget` can also be a list of integers. 
            Budget 0: mean statistics
            Budget 1: momentum
            Budget 2: MKE
            Budget 3: TKE
            Budget 4: uiuj 
            Budget 5: Wake deficit
        include_wakes (bool) : Includes wakes in the returned budget terms if True, default True. 

        Returns
        -------
        t_list (list) : list of tuples of budgets found

        """

        t_list = []
        
        # if no budget is given, look through all saved budgets
        if budget is None: 
            budget_list = self.existing_budgets()
        
        else: 
            # convert to list if integer is given
            if type(budget) != list: 
                budget_list = [budget]
            else: 
                budget_list = budget

        # find budgets matching .npz convention in write_npz()
        if self.associate_npz: 
            filename = self.dir_name + os.sep + self.filename_budgets
            with np.load(filename) as npz: 
                all_terms = npz.files

            if budget is None:  # i.e. requesting all budgets
                return all_terms  # we can stop here without sorting through each budget
            
            tup_list = [BudgetIO.key[t] for t in all_terms]  # list of associated tuples
            t_list = []  # this is the list to be built and returned

            for b in budget_list: 
                t_list += [tup for tup in tup_list if tup[0] == b]

        # find budgets by name matching with PadeOps output conventions
        elif self.associate_padeops: 

            filenames = os.listdir(self.dir_name)
            runid = self.runid
            
            tup_list = []
            
            for b in budget_list:
                if b == 0:
                    terms = np.arange(1,22,1)
                elif b == 1:
                    terms = np.arange(1,15,1)
                elif b == 2:
                    terms = np.arange(1,8,1)
                elif b == 3:
                    terms = np.arange(1,9,1)
                elif b == 4:
                    terms = np.arange(1,37,1) # TODO this assumes that all 4 components of the Rij tensor are present
                else:
                    continue
                
                tup_list += [((b, term)) for term in set(terms)]  # these are all tuples
                
            # convert tuples to keys
            t_list = [BudgetIO.key_xy.inverse[key][0] for key in tup_list]
        
        else: 
            warnings.warn('existing_terms(): No terms found for budget ' + str(budget))
            return []
        
        return t_list

    def Read_x_slice(self, xid, label_list=['u'], tidx_list=[]):
        """
        Reads slices of dumped quantities at a time ID or time IDs. 
        
        Arguments
        ---------
        xid (int) : integer of xid dumped by initialize.F90. NOTE: Fortran indexing starts at 1. 
        label_list (list) : list of terms to read in. Available is typically: "u", "v", "w", and "P" (case-sensitive)
        tidx_list (list) : list of time IDs. 
        
        Returns
        -------
        sl (dict) : formatted dictionary similar to BudgetIO.slice()
        """
        
        sl = {}
        if type(label_list)==str: 
            label_list = [label_list]
        
        for tidx in tidx_list:  
            for lab in label_list: 
                fname = "{:s}/Run{:02d}_t{:06d}_{:s}{:05d}.pl{:s}".format(self.dir_name, self.runid, tidx, 'x', xid, lab)

                key_name = "{:s}_{:d}".format(lab, tidx)
                sl[key_name] = np.fromfile(
                    fname, dtype=np.dtype(np.float64), count=-1).reshape((self.ny,self.nz), order='F')
            
        sl['x'] = self.xLine[[xid-1]]
        sl['y'] = self.yLine
        sl['z'] = self.zLine

        # build and save the extents, either in 1D, 2D, or 3D
        ext = []
        for term in ['x', 'y', 'z']: 
            if len(sl[term]) > 1:  # if this is actually a slice (not a number), then add it to the extents
                ext += [np.min(sl[term]), np.max(sl[term])]

        sl['extent'] = ext

        return sl

    
    def Read_y_slice(self, yid, label_list=['u'], tidx_list=[]):
        """
        Reads slices of dumped quantities at a time ID or time IDs. 
        
        Arguments
        ---------
        yid (int) : integer of yid dumped by initialize.F90
        label_list (list) : list of terms to read in. Available is typically: "u", "v", "w", and "P" (case-sensitive)
        tidx_list (list) : list of time IDs. 
        
        Returns
        -------
        sl (dict) : formatted dictionary similar to BudgetIO.slice()
        """
        
        sl = {}
        if type(label_list)==str: 
            label_list = [label_list]
        
        for tidx in tidx_list:  
            for lab in label_list: 
                fname = "{:s}/Run{:02d}_t{:06d}_{:s}{:05d}.pl{:s}".format(self.dir_name, self.runid, tidx, 'y', yid, lab)

                key_name = "{:s}_{:d}".format(lab, tidx)
                sl[key_name] = np.fromfile(
                    fname, dtype=np.dtype(np.float64), count=-1).reshape((self.nx,self.nz), order='F')
            
        sl['x'] = self.xLine
        sl['y'] = self.yLine[[yid-1]]
        sl['z'] = self.zLine

        # build and save the extents, either in 1D, 2D, or 3D
        ext = []
        for term in ['x', 'y', 'z']: 
            if len(sl[term]) > 1:  # if this is actually a slice (not a number), then add it to the extents
                ext += [np.min(sl[term]), np.max(sl[term])]

        sl['extent'] = ext

        return sl
    
    
    def Read_z_slice(self, zid, label_list=['u'], tidx_list=[]):
        """
        Reads slices of dumped quantities at a time ID or time IDs. 
        
        Arguments
        ---------
        zid (int) : integer of zid dumped by initialize.F90
        label_list (list) : list of terms to read in. Available is typically: "u", "v", "w", and "P" (case-sensitive)
        tidx_list (list) : list of time IDs. 
        
        Returns
        -------
        sl (dict) : formatted dictionary similar to BudgetIO.slice()
        """
        
        sl = {}
        if type(label_list)==str: 
            label_list = [label_list]
        
        for tidx in tidx_list:  
            for lab in label_list: 
                fname = "{:s}/Run{:02d}_t{:06d}_{:s}{:05d}.pl{:s}".format(self.dir_name, self.runid, tidx, 'z', zid, lab)

                key_name = "{:s}_{:d}".format(lab, tidx)
                sl[key_name] = np.fromfile(
                    fname, dtype=np.dtype(np.float64), count=-1).reshape((self.nx,self.ny), order='F')
            
        sl['x'] = self.xLine
        sl['y'] = self.yLine
        sl['z'] = self.zLine[[zid-1]]

        # build and save the extents, either in 1D, 2D, or 3D
        ext = []
        for term in ['x', 'y', 'z']: 
            if len(sl[term]) > 1:  # if this is actually a slice (not a number), then add it to the extents
                ext += [np.min(sl[term]), np.max(sl[term])]

        sl['extent'] = ext

        return sl
    
    
    def read_turb_power(self, tidx=None, turb=1, steady=True, nondim=False): 
        """
        Reads the turbine power from the output file *.pow. 

        tidx (int) : time ID to read turbine power from. Default: calls self.unique_budget_tidx()
        turb (int) : Turbine number. Default 1
        steady (bool) : Averages results if True. If False, returns an array containing the contents of `*.pow`. 
        nondim (bool) : Non-dimensionalizes the output power. This applies to ADM Type 5, and to any other
            turbine model that outputs dimensional power (0.5*rho*u^3*A*C_P)
        """
        if tidx is None: 
            try: 
                tidx = self.last_tidx
            except ValueError as e:   # TODO - Fix this!! 
                tidx = self.unique_tidx(return_last=True)
        
        fname = self.dir_name + '/Run{:02d}_t{:06d}_turbP{:02}.pow'.format(self.runid, tidx, turb)
        if self.verbose: 
            print("  Reading", fname)
            
        power = np.genfromtxt(fname, dtype=float)
                
        pnormfact = 1.
        if nondim: 
            try: 
                # TODO this will break for multiple turbines, probably
                pnormfact = 0.5*np.pi/4*self.turbineArray.diam**2
            except AttributeError as e: 
                print(e)
                print("Tried to normalize power, but failed. Perhaps there is more than one turbine in the simulation?")

        power /= pnormfact

        if steady: 
            return np.mean(power)

        return power/pnormfact  # this is an array

    #### Functions for plotting budgets ####

    # production      - blues
    # transport       - reds
    # redistributions - yellows
    # dissipation     - purples
    # coriolis        - greens
    
    budget_colors = {"shear_production": "tab:blue",
                     "buoyancy" : "tab:gray",
                     "AD" : "tab:cyan",
                     "adv" : "tab:orange", #
                     "p_transport" : "tab:purple",  #
                     "SGS_transport" : "tab:brown",   #
                     "turb_transport" : "tab:green", #
                     "p_strain" : "tab:red", #
                     "dissipation" : "tab:pink", #
                     "coriolis" : 'tab:olive'} #
    '''
    {"shear_production": "#7e82ed", #
                     "buoyancy" : "#0b5394",
                     "AD" : "#1e90ff",
                     "adv" : "gray", #
                     "p_transport" : "#439a1d",  #
                     "SGS_transport" : "#7926c7",   #
                     "turb_transport" : "#f01094", #
                     "p_strain" : "#ffae1d", #
                     "dissipation" : "red", #
                     "coriolis" : '#42b79b'} #
    '''
    def plot_budget_momentum(self, component=None, coords=None):
        '''
        Plots the mean momentum budget for a given component
        
        Arguments
        ---------
        component (int) : vector component of the mean momentum budget to plot
                          1 - streamwise (x)
                          2 - lateral    (y)
                          3 - vertical   (w)

        Returns
        -------
        fig : figure object 
        ax : axes object
        '''

        if component is None:
            # default to the streamwise mean momentum budget
            component = 1
        if component == 1:
            comp_str = ['u', 'x']
        elif component == 2:
            comp_str = ['v', 'y']
        elif component == 3:
            comp_str = ['w', 'z']
        else:
            print("Please enter a valid component number. Valid options are 1, 2, 3, or None (defaults to 1).")
            return None

        if coords is None:
            coords = [(0,len(self.xLine)-1), (0, len(self.yLine)-1), (0, len(self.zLine)-1)]
        
        keys = [key for key in budgetkey.get_key() if budgetkey.get_key()[key][0] == 1]

        keys = [key for key in keys if comp_str[0] in key or comp_str[1] in key]

        key_labels = budgetkey.key_labels()

        fig, ax = plt.subplots()

        residual = 0
        
        for key in keys:
            ax.plot(np.mean(np.mean(self.budget[key][slice(*coords[0]), slice(*coords[1]), slice(*coords[2])], axis=1), axis=0), 
                self.zLine, label=key_labels[key])
            residual += self.budget[key]

        ax.plot(np.mean(np.mean(residual[slice(*coords[0]), slice(*coords[1]), slice(*coords[2])], axis=1), axis=0), 
                self.zLine, label='Residual', linestyle='--', color='black')

        ax.set_ylabel('$z/L$')
            
        return fig, ax

    def plot_budget_tke(self, fig=None, ax=None, linestyle=None, alpha=None, coords=None):
        '''
        Plots the tke budget
        
        Arguments
        ---------
        
        Returns
        -------
        fig : figure object 
        ax : axes object
        '''
        
        keys = [key for key in budgetkey.get_key() if budgetkey.get_key()[key][0] == 3]
        key_labels = budgetkey.key_labels()

        if coords is None:
            coords = [(0,len(self.xLine)-1), (0, len(self.yLine)-1), (0, len(self.zLine)-1)]

        if not linestyle:
            linestyle = '-'
        if not alpha:
            alpha = 1

        if not ax and not fig:
            fig, ax = plt.subplots()

        residual = 0
        
        for key in keys:
            ax.plot(np.mean(np.mean(self.budget[key][slice(*coords[0]), slice(*coords[1]), slice(*coords[2])], axis=1), axis=0), 
                self.zLine, label=key,  color = self.budget_colors[key.replace("TKE_", "")], 
                linestyle=linestyle, alpha=alpha) 
            residual += self.budget[key]

        ax.plot(np.mean(np.mean(residual[slice(*coords[0]), slice(*coords[1]), slice(*coords[2])], axis=1), axis=0), 
            self.zLine, label='Residual', linestyle='--', color='black')

        ax.set_ylabel('$z/L$')
            
        return fig, ax

    def plot_budget_mke(self, coords=None):
        '''
        Plots the tke budget
        
        Arguments
        ---------
        
        Returns
        -------
        fig : figure object 
        ax : axes object
        '''
        
        if coords is None:
            coords = [(0,len(self.xLine)-1), (0, len(self.yLine)-1), (0, len(self.zLine)-1)]


        keys = [key for key in budgetkey.get_key() if budgetkey.get_key()[key][0] == 2]
        key_labels = budgetkey.key_labels()

        fig, ax = plt.subplots()

        residual = 0
        
        for key in keys:
            ax.plot(np.mean(np.mean(self.budget[key][slice(*coords[0]), slice(*coords[1]), slice(*coords[2])], axis=1), axis=0), self.zLine, label=key)
            residual += self.budget[key]

        ax.plot(np.mean(np.mean(residual[slice(*coords[0]), slice(*coords[1]), slice(*coords[2])], axis=1), axis=0), 
            self.zLine, label='Residual', linestyle='--', color='black')

        ax.set_ylabel('$z/L$')
            
        return fig, ax

    def plot_budget_uiuj(self, component, fig=None, ax=None, linestyle=None, alpha=None, coords=None):
        '''
        Plots the tke budget
        
        Arguments
        ---------
        
        Returns
        -------
        fig : figure object 
        ax : axes object
        '''

        if coords is None:
            coords = [(0,len(self.xLine)-1), (0, len(self.yLine)-1), (0, len(self.zLine)-1)]


        comp_dict = {11 : [1, 10],
                     22 : [11, 20],
                     33 : [21, 30],
                     13 : [31, 40],
                     23 : [41, 50]}

        comp_str_dict = {11 : 'uu',
                         22 : 'vv',
                         33 : 'ww',
                         13 : 'uw',
                         23 : 'vw'}

        
        if not linestyle:
            linestyle = '-'

        if not alpha:
            alpha = 1

        keys = [key for key in budgetkey.get_key() if budgetkey.get_key()[key][0] == 4 and budgetkey.get_key()[key][1] in range(comp_dict[component][0], comp_dict[component][1])]
        
        if not ax and not fig:
            fig, ax = plt.subplots()

        residual = 0
        
        for key in keys:
            ax.plot(np.mean(np.mean(self.budget[key][slice(*coords[0]), slice(*coords[1]), slice(*coords[2])], axis=1), axis=0), 
                self.zLine, label=key, color = self.budget_colors[key.replace(comp_str_dict[component] + "_", "")], 
                linestyle=linestyle, alpha=alpha)
            residual += self.budget[key]

        ax.plot(np.mean(np.mean(residual[slice(*coords[0]), slice(*coords[1]), slice(*coords[2])], axis=1), axis=0), 
            self.zLine, label='Residual', linestyle=linestyle, color='black')

        ax.set_ylabel('$z/L$')
            
        return fig, ax

    
    def plot_budget_xy_momentum(self, component=None):
        '''
        Plots the mean momentum budget for a given component
        
        Arguments
        ---------
        component (int) : vector component of the mean momentum budget to plot
                          1 - streamwise (x)
                          2 - lateral    (y)
                          3 - vertical   (w)

        Returns
        -------
        fig : figure object 
        ax : axes object
        '''

        if component is None:
            # default to the streamwise mean momentum budget
            component = 1

        if component == 1:
            comp_str = ['u', 'x']
        elif component == 2:
            comp_str = ['v', 'y']
        elif component == 3:
            comp_str = ['w', 'z']
        else:
            print("Please enter a valid component number. Valid options are 1, 2, 3, or None (defaults to 1).")
            return None
        
        keys = [key for key in budgetkey.get_key_xy() if budgetkey.get_key_xy()[key][0] == 1]

        keys = [key for key in keys if comp_str[0] in key or comp_str[1] in key]
        
        fig, ax = plt.subplots()

        residual = 0
        
        for key in keys:
            ax.plot(self.budget_xy[key], self.zLine, label = key)
            residual += self.budget_xy[key]

        ax.plot(residual, self.zLine, label='Residual', linestyle='--', color='black')

        ax.set_ylabel('$z/L$')
            
        return fig, ax
    

    def plot_budget_xy_tke(self, fig=None, ax=None, linestyle=None, alpha=None):
        '''
        Plots the tke budget
        
        Arguments
        ---------
        
        Returns
        -------
        fig : figure object 
        ax : axes object
        '''

        if not linestyle:
            linestyle = '-'
        if not alpha:
            alpha = 1
        
        keys = [key for key in budgetkey.get_key_xy() if budgetkey.get_key_xy()[key][0] == 3]
        key_labels = budgetkey.key_labels()

        if not ax and not fig:
            fig, ax = plt.subplots()

        residual = 0
        
        for key in keys:
            ax.plot(self.budget_xy[key], self.zLine, label=key,color = self.budget_colors[key.replace("TKE_", "")], 
                    linestyle=linestyle, alpha=alpha)
            residual += self.budget_xy[key]

        ax.plot(residual, self.zLine, label='Residual', linestyle='--', color='black')

        ax.set_ylabel('$z/L$')
            
        return fig, ax


    def plot_budget_xy_mke(self):
        '''
        Plots the tke budget
        
        Arguments
        ---------
        
        Returns
        -------
        fig : figure object 
        ax : axes object
        '''
        
        keys = [key for key in budgetkey.get_key_xy() if budgetkey.get_key_xy()[key][0] == 2]
        key_labels = budgetkey.key_labels()

        fig, ax = plt.subplots()

        residual = 0
        
        for key in keys:
            ax.plot(self.budget_xy[key], self.zLine, label=key)
            residual += self.budget_xy[key]

        ax.plot(residual, self.zLine, label='Residual', linestyle='--', color='black')

        ax.set_ylabel('$z/L$')
            
        return fig, ax

    
    def plot_budget_xy_uu(self, fig=None, ax=None, linestyle=None, alpha=None):
        '''
        Plots the streamwise variance <uu> budget
        
        Arguments
        ---------
        
        Returns
        -------
        fig : figure object 
        ax : axes object
        '''

        if not linestyle:
            linestyle = '-'

        if not alpha:
            alpha = 1
        
        keys = [key for key in budgetkey.get_key_xy() if budgetkey.get_key_xy()[key][0] == 4 and budgetkey.get_key_xy()[key][1] >= 1 and budgetkey.get_key_xy()[key][1] <= 9]
        key_labels = budgetkey.key_labels()

        if not ax and not fig:
            fig, ax = plt.subplots()

        residual = 0
        
        for key in keys:
            ax.plot(self.budget_xy[key], self.zLine, label = key, color = self.budget_colors[key.replace("uu_", "")], linestyle=linestyle, alpha=alpha)
            residual += self.budget_xy[key]

        ax.plot(residual, self.zLine, label='Residual', linestyle='--', color='black')

        ax.set_ylabel('$z/L$')
            
        return fig, ax


    def plot_budget_xy_uw(self, fig=None, ax=None, linestyle=None, alpha=None):
        '''
        Plots the streamwise variance <uu> budget
        
        Arguments
        ---------
        
        Returns
        -------
        fig : figure object 
        ax : axes object
        '''
        
        if not linestyle:
            linestyle = '-'

        if not alpha:
            alpha = 1

        keys = [key for key in budgetkey.get_key_xy() if budgetkey.get_key_xy()[key][0] == 4 and budgetkey.get_key_xy()[key][1] >= 10 and budgetkey.get_key_xy()[key][1] <= 18]
        key_labels = budgetkey.key_labels()

        if not ax and not fig:
            fig, ax = plt.subplots()

        residual = 0
        
        for key in keys:
            ax.plot(self.budget_xy[key], self.zLine, label = key, color = self.budget_colors[key.replace("uw_", "")], linestyle=linestyle, alpha=alpha)
            residual += self.budget_xy[key]

        ax.plot(residual, self.zLine, label='Residual', linestyle='--', color='black')

        ax.set_ylabel('$z/L$')
            
        return fig, ax


    def plot_budget_xy_vw(self, fig=None, ax=None, linestyle=None, alpha=None):
        '''
        Plots the streamwise variance <uu> budget
        
        Arguments
        ---------
        
        Returns
        -------
        fig : figure object 
        ax : axes object
        '''

        
        if not linestyle:
            linestyle = '-'

        if not alpha:
            alpha = 1
        
        keys = [key for key in budgetkey.get_key_xy() if budgetkey.get_key_xy()[key][0] == 4 and budgetkey.get_key_xy()[key][1] >= 19 and budgetkey.get_key_xy()[key][1] <= 27]
        key_labels = budgetkey.key_labels()

        if not ax and not fig:
            fig, ax = plt.subplots()

        residual = 0
        
        for key in keys:
            ax.plot(self.budget_xy[key], self.zLine, label = key, color = self.budget_colors[key.replace("vw_", "")], linestyle=linestyle, alpha=alpha)
            residual += self.budget_xy[key]

        ax.plot(residual, self.zLine, label='Residual', linestyle='--', color='black')

        ax.set_ylabel('$z/L$')
            
        return fig, ax


    def plot_budget_xy_ww(self, fig=None, ax=None, linestyle=None, alpha=None):
        '''
        Plots the streamwise variance <uu> budget
        
        Arguments
        ---------
        
        Returns
        -------
        fig : figure object 
        ax : axes object
        '''
        
        
        if not linestyle:
            linestyle = '-'

        if not alpha:
            alpha = 1

        keys = [key for key in budgetkey.get_key_xy() if budgetkey.get_key_xy()[key][0] == 4 and budgetkey.get_key_xy()[key][1] >= 28 and budgetkey.get_key_xy()[key][1] <= 36]
        key_labels = budgetkey.key_labels()

        if not ax and not fig:
            fig, ax = plt.subplots()

        residual = 0
        
        for key in keys:
            ax.plot(self.budget_xy[key], self.zLine, label = key, color = self.budget_colors[key.replace("ww_", "")], linestyle=linestyle, alpha=alpha)
            residual += self.budget_xy[key]

        ax.plot(residual, self.zLine, label='Residual', linestyle='--', color='black')

        ax.set_ylabel('$z/L$')
            
        return fig, ax

if __name__ == "__main__": 
    """
    TODO - add unit tests to class
    """
    print("padeopsIO: No unit tests included yet. ")
