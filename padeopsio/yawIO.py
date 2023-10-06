import numpy as np
import os
import re
import warnings
import glob

import padeopsio.budgetIO as pio

class YawIO(pio.BudgetIO): 
    """
    Class that extends BudgetIO which adds helper functions for reading turbine power, yaw, etc. 
    """
    
    def __init__(self, dir_name, **kwargs): 
        """
        Calls the constructor of BudgetIO
        """
        
        super().__init__(dir_name, **kwargs)
        
        if self.associate_nml: 
            try: 
                self.yaw = self.input_nml['ad_coriolisinput']['yaw']
                self.uInflow = self.input_nml['ad_coriolisinput']['uinflow'] * np.cos(self.yaw*np.pi/180.)
                self.vInflow = self.input_nml['ad_coriolisinput']['uinflow'] * -np.sin(self.yaw*np.pi/180.)
            except KeyError as e: 
                if self.verbose: 
                    print(repr(e)) 
                    print('Initializing YawIO, could not find `ad_coriolisinput` in the Namelist input.')
                warnings.warn('Could not find `ad_coriolisinput` in the Namelist input.')
                
        if self.verbose: 
            print("Initialized YawIO object")


    def read_turb_vel(self, tidx=None, turb=1, steady=True, u=True, v=True, rotate=False): 
        """
        Reads the turbine power from the output file *.pow. 

        tidx (int) : time ID to read turbine power from. Default: calls self.unique_budget_tidx()
        turb (int) : Turbine number. Default 1
        steady (bool) : Averages results if True. If False, returns an array containing the contents of `*.pow`. 
        u, v (bool) : dictates whether to return u, v, or both. Default: u=True, v=True
        rotate (bool) : rotates uTurb, vTurb to be aligned with the freestream flow. Will return both
            u and v velocities if rotate=True. Default False. 
        """
        
        if tidx is None: 
            try: 
                tidx = self.last_tidx
            except ValueError as e:   # TODO - Fix this!! 
                tidx = self.unique_tidx(return_last=True)
        
        ret = ()
        
        if rotate: 
            # if rotating the domain, need to read U and V regardless of what is output
            u_read, v_read = True, True
        else: 
            u_read, v_read = u, v
        
        # read in velocity files
        for i, ui in enumerate((u_read, v_read)): 
            if ui: 
                if i == 0: 
                    u_string = "U"
                else: 
                    u_string = "V"
                    
                fname = self.dir_name + '/Run{:02d}_t{:06d}_turb{:s}{:02}.vel'.format(self.runid, tidx, u_string, turb)
                uturb = np.genfromtxt(fname, dtype=float)

                if steady: 
                    ret += (np.mean(uturb), )

                else: 
                    ret += (uturb, )  # this is an array
                    
        if rotate: 
            # Do the euler rotations now
            yaw = self.yaw * np.pi / 180  # yaw in radians
            new_ret = ()  # this will be returned
            
            if u: 
                new_ret += (ret[0] * np.cos(yaw) - ret[1] * np.sin(yaw), )
            if v: 
                new_ret += (ret[0] * np.sin(yaw) + ret[1] * np.cos(yaw), )
            
            ret = new_ret  # there are an infinite number of better ways to have done this
                    
        if len(ret) == 0: 
            raise ValueError("u or v must be True, function cannot return nothing")
        if len(ret) == 1: 
            return ret[0]
        else: 
            return ret
        
        
    def rotate_uv(self, overwrite=True, load_fields=False, u=None, v=None): 
        """
        Rotate the ubar and vbar budget field. Loads 'ubar', 'vbar' in the YawIO object. 
        
        Parameters 
        ----------
        overwrite (bool) : Overwrites 'ubar' and 'vbar' in self.budget if True; default True. If False, then 
            rotate_uv() returns the rotated u_bar, v_bar arrays. 
        load_fields (bool) : uses instanteous field variables instead of budgets. 
        """

        # always reload u, v, to make sure this isn't a doubly-rotated field
        if u is None and v is None: 
            if load_fields: 
                self.read_fields(field_terms=['u', 'v']) # TODO: there is no overwrite check on this! 
                u = self.field['u']
                v = self.field['v']
            else: 
                self.read_budgets(budget_terms=['ubar', 'vbar'], overwrite=True)  
                u = self.budget['ubar']
                v = self.budget['vbar']
        else: 
            overwrite=False
            if self.verbose: 
                print('\trotate_uv(): Using input u, v fields. ')

        yaw = self.yaw * np.pi / 180  # yaw in radians
        new_ubar = u * np.cos(yaw) - v * np.sin(yaw)
        new_vbar = u * np.sin(yaw) + v * np.cos(yaw)

        # overwrite

        if self.verbose: 
            print("Rotated ubar and vbar fields to align with the freestream. ")

        if overwrite: 
            if load_fields: 
                self.field['u'] = new_ubar
                self.field['v'] = new_vbar
            else: 
                self.budget['ubar'] = new_ubar
                self.budget['vbar'] = new_vbar

        else: 
            return (new_ubar, new_vbar)
        

