import numpy as np

import padeopsIO.budgetIO as pio
import padeopsIO.deficitkey as deficitkey

class DeficitIO(): 
    """
    Interface class for computing the base and deficit budgets using a precursor BudgetIO object
    and a primary BudgetIO object. 
    """

    key = deficitkey.get_key()

    def __init__(self, pre, prim, dir_name): 
        """
        Initaliazes the instance of DeficitIO object
        ----------------------------------
        parameters (required):
        pre  - BudgetIO object for precursor
        prim - BudgetIO object for primary
        ----------------------------------

        """

        self.pre = pre
        self.prim = prim

        self.dir_name = dir_name

        self.budget = {}

    
    def read_fields(self, tidx):
        self.pre.read_fields(tidx=tidx)
        self.prim.read_fields(tidx=tidx)

    def do_budgets(self, tidxs=None, save=True, budget=0):
        """
        Compute and save the momentum budget
        """

        if tidxs is None:
            tidxs = self.unique_tidx()
            
        self.counter = 0

        self.last_tidx = tidxs[-1]

        if budget==0:
            self.initialize_budget_0()
            for tidx in tidxs:
                self.read_fields(tidx)

                self.compute_budget_0()
            
                # increment counter
                self.counter += 1
        
            self.assemble_budget_0(save=save)

        elif budget==1:
            self.initialize_budget_0()
            self.initialize_budget_1()

            for tidx in tidxs:
                self.read_fields(tidx)

                self.compute_budget_0()
            
                # increment counter
                self.counter += 1
        
            self.assemble_budget_0(save=save)
            self.assemble_budget_1(save=save)

    def initialize_budget_0(self):
        """
        Initialize the terms in the base and deficit moment budgets
        """
        self.budget = {}

        # util vars
        self.budget['ui_pre']            = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz, 3])
        self.budget['ui_prim']            = np.zeros([self.prim.nx, self.prim.ny, self.prim.nz, 3])

        # base flow variables
        self.budget['ui_base']       = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz, 3])               
        self.budget['p_base']        = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz])                  
        self.budget['T_base']        = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz])                  
        self.budget['kSGS_base']     = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz])                  
        self.budget['nSGS_base']     = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz])                  
    
        self.budget['dpdxi_base']        = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz, 3])           
        self.budget['duidxj_base']       = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz, 3, 3])        
        self.budget['uiuj_base']         = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz, 3, 3])        
        self.budget['duiujdxk_base']     = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz, 3, 3, 3])     
        self.budget['tauij_base']        = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz, 3, 3])        
        self.budget['dtauijdxk_base']    = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz, 3, 3, 3])     
        self.budget['pui_base']          = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz, 3])  

        # deficit variables
        self.budget['delta_ui']      = np.zeros([self.prim.nx, self.prim.ny, self.prim.nz,3])            
        self.budget['delta_p']       = np.zeros([self.prim.nx, self.prim.ny, self.prim.nz])               
        self.budget['delta_T']       = np.zeros([self.prim.nx, self.prim.ny, self.prim.nz])               
        self.budget['delta_kSGS']    = np.zeros([self.prim.nx, self.prim.ny, self.prim.nz])               
        self.budget['delta_nSGS']    = np.zeros([self.prim.nx, self.prim.ny, self.prim.nz])               
        
        self.budget['delta_dpdxi']       = np.zeros([self.prim.nx, self.prim.ny, self.prim.nz, 3])        
        self.budget['delta_duidxj']      = np.zeros([self.prim.nx, self.prim.ny, self.prim.nz, 3, 3])     
        self.budget['delta_uiuj']        = np.zeros([self.prim.nx, self.prim.ny, self.prim.nz, 3, 3])     
        self.budget['delta_duiujdxk']    = np.zeros([self.prim.nx, self.prim.ny, self.prim.nz, 3, 3, 3])   
        self.budget['delta_tauij']       = np.zeros([self.prim.nx, self.prim.ny, self.prim.nz, 3, 3])         
        self.budget['delta_dtauijdxk']   = np.zeros([self.prim.nx, self.prim.ny, self.prim.nz, 3, 3, 3]) 
        self.budget['delta_pui']         = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz, 3])     

        # mixed variables
        self.budget['ui_base_delta_uj']      = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz, 3, 3])
        self.budget['ddxk_ui_base_delta_uj'] = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz, 3, 3, 3])
        self.budget['delta_p_ui_base']       = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz, 3])         
        self.budget['p_base_delta_ui']       = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz, 3])         

    def initialize_budget_1(self):
        """
        Initialize the terms in the base and deficit moment budgets
        """
    def initialize_budget_2(self):
        """
        Initialize the terms in the base and deficit mean kinetic energy budgets
        """

    def compute_budget_0(self):
        """
        Computes the base and deficit momentum budget terms
        """
        
        # these are used to more efficiently compute Reynolds stresses
        self.budget['ui_pre'][:,:,:,0] = self.pre.field['u']
        self.budget['ui_pre'][:,:,:,1] = self.pre.field['v']
        self.budget['ui_pre'][:,:,:,2] = self.pre.field['w']

        self.budget['ui_prim'][:,:,:,0] = self.prim.field['u']
        self.budget['ui_prim'][:,:,:,1] = self.prim.field['v']
        self.budget['ui_prim'][:,:,:,2] = self.prim.field['w']

        self.budget['p_base'] = np.add(self.budget['p_base'], self.pre.field['p'])
        self.budget['T_base'] = np.add(self.budget['T_base'], self.pre.field['T'])
        self.budget['kSGS_base'] = np.add(self.budget['kSGS_base'], self.pre.field['kSGS'])
        self.budget['nSGS_base'] = np.add(self.budget['nSGS_base'], self.pre.field['nSGS'])
    
        self.budget['delta_p'] = np.add(self.budget['delta_p'], np.subtract(self.prim.field['p'], self.pre.field['p']))
        self.budget['delta_T'] = np.add(self.budget['delta_T'], np.subtract(self.prim.field['T'], self.pre.field['T']))
        self.budget['delta_kSGS'] = np.add(self.budget['delta_kSGS'], np.subtract(self.prim.field['kSGS'], self.pre.field['kSGS']))
        self.budget['delta_nSGS'] = np.add(self.budget['delta_nSGS'], np.subtract(self.prim.field['nSGS'], self.pre.field['nSGS']))
        
        # Compute base flow velocity
        self.budget['ui_base'] = np.add(self.budget['ui_base'], self.budget['ui_pre'])
            
        # Compute wake flow velocity
        self.budget['delta_ui'] = np.add(self.budget['delta_ui'], np.subtract(self.budget['ui_prim'], self.budget['ui_pre']))

        # Compute base flow Reynolds stresses
        self.budget['uiuj_base'] = np.add(self.budget['uiuj_base'], self.budget['ui_pre'][:,:,:,:,None]*self.budget['ui_pre'][:,:,:,None,:])
            
        # Compute wake flow Reynolds stresses
        self.budget['delta_uiuj'] = np.add(self.budget['delta_uiuj'], np.subtract(self.budget['ui_prim'], self.budget['ui_pre'])[:,:,:,:,None]
                                *np.subtract(self.budget['ui_prim'], self.budget['ui_pre'])[:,:,:,None,:])

        # Compute mixed Reynolds stresses
        self.budget['ui_base_delta_uj'] = np.add(self.budget['ui_base_delta_uj'], self.budget['ui_pre'][:,:,:,:,None]
                                                 *np.subtract(self.budget['ui_prim'],self.budget['ui_pre'])[:,:,:,None,:]) 

        # Compute velocity-pressure correlation
        self.budget['pui_base'] = np.add(self.budget['pui_base'], self.budget['ui_pre']*self.pre.field['p'][:,:,:,None])
        self.budget['delta_pui'] = np.add(self.budget['delta_pui'], np.subtract(self.budget['ui_prim'],self.budget['ui_pre'])\
                                    *np.subtract(self.prim.field['p'],self.pre.field['p'])[:,:,:,None])
        self.budget['delta_p_ui_base'] = np.add(self.budget['delta_p_ui_base'], self.budget['ui_pre']*np.subtract(self.prim.field['p'],self.pre.field['p'])[:,:,:,None])
        self.budget['p_base_delta_ui'] = np.add(self.budget['p_base_delta_ui'], np.subtract(self.budget['ui_prim'],self.budget['ui_pre'])*self.pre.field['p'][:,:,:,None])

    def assemble_budget_0(self, save=False):
        
        # initalize final budget variables (that will be saved if files are written)
        keys = [key for key in deficitkey.get_key() if deficitkey.get_key()[key][1] == 0]
        for key in keys:
            self.budget['key'] = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz])               

        # variables to help build final variables

        counter = self.counter
        
        x_pre = self.pre.xLine
        y_pre = self.pre.yLine
        z_pre = self.pre.zLine
        
        x_prim = self.prim.xLine
        y_prim = self.prim.yLine
        z_prim = self.prim.zLine

        # average
        self.budget['ui_base']       = self.budget['ui_base']/counter
        self.budget['p_base']        = self.budget['p_base']/counter
        self.budget['T_base']        = self.budget['T_base']/counter
        self.budget['kSGS_base']     = self.budget['kSGS_base']/counter
        self.budget['nSGS_base']     = self.budget['nSGS_base']/counter
        self.budget['uiuj_base']     = self.budget['uiuj_base']/counter
        
        self.budget['delta_ui']      = self.budget['delta_ui']/counter
        self.budget['delta_p']       = self.budget['delta_p']/counter
        self.budget['delta_T']       = self.budget['delta_T']/counter
        self.budget['delta_kSGS']    = self.budget['delta_kSGS']/counter
        self.budget['delta_nSGS']    = self.budget['delta_nSGS']/counter
        self.budget['delta_uiuj']    = self.budget['delta_uiuj']/counter
    
        self.budget['ui_base_delta_uj']    = self.budget['ui_base_delta_uj']/counter
    
        # calculate derivatives and other statistics
        for ii in range(3):
            self.budget['duidxj_base'][:,:,:,ii,:] = np.transpose(np.gradient(self.budget['ui_base'][:,:,:,ii], x_pre, y_pre, z_pre), [1,2,3,0])
        
            self.budget['delta_duidxj'][:,:,:,ii,:] = np.transpose(np.gradient(self.budget['delta_ui'][:,:,:,ii], x_prim, y_prim, z_prim), [1,2,3,0])
        
        
        self.budget['uiuj_base'] = np.subtract(self.budget['uiuj_base'], self.budget['ui_base'][:,:,:,:,None]*self.budget['ui_base'][:,:,:,None,:])
            
        self.budget['delta_uiuj'] = np.subtract(self.budget['delta_uiuj'], self.budget['delta_ui'][:,:,:,:,None]*self.budget['delta_ui'][:,:,:,None,:])
            
        self.budget['ui_base_delta_uj'] = np.subtract(self.budget['ui_base_delta_uj'], 
                                                      self.budget['ui_base'][:,:,:,:,None]*self.budget['delta_ui'][:,:,:,None,:])
        
        self.budget['sij_base'] = 0.5*np.add(self.budget['duidxj_base'], np.transpose(self.budget['duidxj_base'], [0,1,2,4,3]))
        self.budget['tauij_base'] = 2.0*self.budget['nSGS_base'][:,:,:,None,None]*self.budget['sij_base']

        self.budget['delta_sij'] = 0.5*np.add(self.budget['delta_duidxj'], np.transpose(self.budget['delta_duidxj'], [0,1,2,4,3]))
        self.budget['delta_tauij'] = 2.0*self.budget['delta_nSGS'][:,:,:,None,None]*self.budget['delta_sij']
      
        #####################
        ##### BASE FLOW #####
        #####################
        
        # Mean velocity
        self.budget['u_base'] = self.budget['ui_base'][:,:,:,0]
        self.budget['v_base'] = self.budget['ui_base'][:,:,:,1]
        self.budget['w_base'] = self.budget['ui_base'][:,:,:,2]

        # Reynolds stresses
        self.budget['uu_base'] = self.budget['uiuj_base'][:,:,:,0,0] 
        self.budget['vv_base'] = self.budget['uiuj_base'][:,:,:,1,1]
        self.budget['ww_base'] = self.budget['uiuj_base'][:,:,:,2,2] 
        self.budget['uw_base'] = self.budget['uiuj_base'][:,:,:,0,2]         
        self.budget['vw_base'] = self.budget['uiuj_base'][:,:,:,1,2] 
        self.budget['uv_base'] = self.budget['uiuj_base'][:,:,:,0,1] 

        # SGS
        self.budget['tau11_base'] = self.budget['tauij_base'][:,:,:,0,0]
        self.budget['tau22_base'] = self.budget['tauij_base'][:,:,:,1,1]
        self.budget['tau33_base'] = self.budget['tauij_base'][:,:,:,2,2]
        self.budget['tau13_base'] = self.budget['tauij_base'][:,:,:,0,2]
        self.budget['tau23_base'] = self.budget['tauij_base'][:,:,:,1,2]
        self.budget['tau12_base'] = self.budget['tauij_base'][:,:,:,0,1]

        # PRESSURE
        self.budget['pu_base'] = self.budget['pui_base'][:,:,:,0]
        self.budget['pv_base'] = self.budget['pui_base'][:,:,:,1]
        self.budget['pw_base'] = self.budget['pui_base'][:,:,:,2]

        #####################
        ##### WAKE FLOW #####
        #####################
        
        # Mean velocity    
        self.budget['delta_u'] = self.budget['delta_ui'][:,:,:,0]
        self.budget['delta_v'] = self.budget['delta_ui'][:,:,:,1]
        self.budget['delta_w'] = self.budget['delta_ui'][:,:,:,2]
        
        # Reynolds stresses
        self.budget['delta_uu'] = self.budget['delta_uiuj'][:,:,:,0,0]
        self.budget['delta_vv'] = self.budget['delta_uiuj'][:,:,:,1,1]
        self.budget['delta_ww'] = self.budget['delta_uiuj'][:,:,:,2,2]
        self.budget['delta_uw'] = self.budget['delta_uiuj'][:,:,:,0,2]
        self.budget['delta_vw'] = self.budget['delta_uiuj'][:,:,:,1,2]
        self.budget['delta_uv'] = self.budget['delta_uiuj'][:,:,:,0,1]

        # Mixed Reynolds stresses
        self.budget['u_base_delta_u'] = self.budget['ui_base_delta_uj'][:,:,:,0,0]
        self.budget['v_base_delta_v'] = self.budget['ui_base_delta_uj'][:,:,:,1,1]
        self.budget['w_base_delta_w'] = self.budget['ui_base_delta_uj'][:,:,:,2,2]
        self.budget['u_base_delta_w'] = self.budget['ui_base_delta_uj'][:,:,:,0,2]
        self.budget['w_base_delta_u'] = self.budget['ui_base_delta_uj'][:,:,:,2,0]
        self.budget['v_base_delta_w'] = self.budget['ui_base_delta_uj'][:,:,:,1,2]
        self.budget['w_base_delta_v'] = self.budget['ui_base_delta_uj'][:,:,:,2,1]
        
        # SGS
        self.budget['delta_tau11'] = self.budget['delta_tauij'][:,:,:,0,0]
        self.budget['delta_tau22'] = self.budget['delta_tauij'][:,:,:,1,1]
        self.budget['delta_tau33'] = self.budget['delta_tauij'][:,:,:,2,2]
        self.budget['delta_tau13'] = self.budget['delta_tauij'][:,:,:,0,2]
        self.budget['delta_tau23'] = self.budget['delta_tauij'][:,:,:,1,2]
        self.budget['delta_tau12'] = self.budget['delta_tauij'][:,:,:,0,1]

        # PRESSURE
        self.budget['delta_pu'] = self.budget['delta_pui'][:,:,:,0]
        self.budget['delta_pv'] = self.budget['delta_pui'][:,:,:,1]
        self.budget['delta_pw'] = self.budget['delta_pui'][:,:,:,2]

        # mixed
        self.budget['delta_p_u_base'] = self.budget['delta_p_ui_base'][:,:,:,0]
        self.budget['delta_p_v_base'] = self.budget['delta_p_ui_base'][:,:,:,1]
        self.budget['delta_p_w_base'] = self.budget['delta_p_ui_base'][:,:,:,2]
        self.budget['p_base_delta_u'] = self.budget['p_base_delta_ui'][:,:,:,0]
        self.budget['p_base_delta_v'] = self.budget['p_base_delta_ui'][:,:,:,1]
        self.budget['p_base_delta_w'] = self.budget['p_base_delta_ui'][:,:,:,2]

        if save:
            self.write_budget(budget=0)

    def assemble_budget_1(self, save=False):

        x_pre = self.pre.xLine
        y_pre = self.pre.yLine
        z_pre = self.pre.zLine

        x_prim = self.prim.xLine
        y_prim = self.prim.yLine
        z_prim = self.prim.zLine

        latitude = self.prim.latitude
        g_alpha = -16.072 # this is hardcoded for now

        # initialize final budget variables (that will be saved if files are written)
        keys = [key for key in deficitkey.get_key() if deficitkey.get_key()[key][1] == 1]
        for key in keys:
            self.budget[key] = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz])               

        # calculate derivatives and other statistics
        self.budget['dpdxi_base'] = np.transpose(np.gradient(self.budget['p_base'], x_pre, y_pre, z_pre), [1,2,3,0])
        
        self.budget['delta_dpdxi'] = np.transpose(np.gradient(self.budget['delta_p'], x_prim, y_prim, z_prim), [1,2,3,0])
        
        for ii in range(3):
            for jj in range(3):
                self.budget['duiujdxk_base'][:,:,:,ii,jj,:] = np.transpose(np.gradient(self.budget['uiuj_base'][:,:,:,ii,jj], 
                                                                                       x_pre, y_pre, z_pre),[1,2,3,0])            
                self.budget['delta_duiujdxk'][:,:,:,ii,jj,:] = np.transpose(np.gradient(self.budget['delta_uiuj'][:,:,:,ii,jj], 
                                                                                        x_prim, y_prim, z_prim), [1,2,3,0])
                self.budget['ddxk_ui_base_delta_uj'][:,:,:,ii,jj,:] = np.transpose(np.gradient(self.budget['ui_base_delta_uj'][:,:,:,ii,jj], 
                                                                                        x_prim, y_prim, z_prim), [1,2,3,0])
                self.budget['dtauijdxk_base'][:,:,:,ii,jj,:] = np.transpose(np.gradient(self.budget['tauij_base'][:,:,:,ii,jj], 
                                                                                        x_pre, y_pre, z_pre), [1,2,3,0])        
                self.budget['delta_dtauijdxk'][:,:,:,ii,jj,:] = np.transpose(np.gradient(self.budget['delta_tauij'][:,:,:,ii,jj], 
                                                                                         x_prim, y_prim, z_prim), [1,2,3,0])
      
        #####################
        ##### BASE FLOW #####
        #####################

        # momentum budget terms
        # Advection
        self.budget['xadv_base'] = - np.sum(np.multiply(self.budget['ui_base'], self.budget['duidxj_base'][:,:,:,0,:]),3)
        self.budget['yadv_base'] = - np.sum(np.multiply(self.budget['ui_base'], self.budget['duidxj_base'][:,:,:,1,:]),3)
        self.budget['zadv_base'] = - np.sum(np.multiply(self.budget['ui_base'], self.budget['duidxj_base'][:,:,:,2,:]),3)
        
        # Pressure gradient
        self.budget['dpdx_base'] = - self.budget['dpdxi_base'][:,:,:,0]
        self.budget['dpdy_base'] = - self.budget['dpdxi_base'][:,:,:,1]
        self.budget['dpdz_base'] = - self.budget['dpdxi_base'][:,:,:,2]

        # SGS and Reynolds stress divergence
        for jj in range(3):
            self.budget['xSGS_base'] += self.budget['dtauijdxk_base'][:,:,:,0,jj,jj]
            self.budget['ySGS_base'] += self.budget['dtauijdxk_base'][:,:,:,1,jj,jj]
            self.budget['zSGS_base'] += self.budget['dtauijdxk_base'][:,:,:,2,jj,jj]
        
            self.budget['xturb_base'] -= self.budget['duiujdxk_base'][:,:,:,0,jj,jj]
            self.budget['yturb_base'] -= self.budget['duiujdxk_base'][:,:,:,1,jj,jj]
            self.budget['zturb_base'] -= self.budget['duiujdxk_base'][:,:,:,2,jj,jj]
        
        # Coriolis
        self.budget['xCor_base'] = -np.sin(np.deg2rad(latitude))*(2/self.pre.Ro)\
                            *(np.sin(np.deg2rad(g_alpha)) - self.budget['ui_base'][:,:,:,1])
        self.budget['yCor_base'] = np.sin(np.deg2rad(latitude))*(2/self.pre.Ro)\
                            *(np.cos(np.deg2rad(g_alpha)) - self.budget['ui_base'][:,:,:,0])
        
        # Buoyancy
        self.budget['zbuoy_base'] =  -self.budget['T_base']*(1/(self.prim.Fr*self.prim.Fr*self.prim.Tref))


        #####################
        ##### WAKE FLOW #####
        #####################

        # momentum budget terms
        # Advection
        self.budget['xadv_delta_base'] = - np.sum(np.multiply(self.budget['ui_base'], self.budget['delta_duidxj'][:,:,:,0,:]),3)
        self.budget['yadv_delta_base'] = - np.sum(np.multiply(self.budget['ui_base'], self.budget['delta_duidxj'][:,:,:,1,:]),3)
        self.budget['zadv_delta_base'] = - np.sum(np.multiply(self.budget['ui_base'], self.budget['delta_duidxj'][:,:,:,2,:]),3)
        self.budget['xadv_delta_delta'] = - np.sum(np.multiply(self.budget['delta_ui'], self.budget['delta_duidxj'][:,:,:,0,:]),3)
        self.budget['yadv_delta_delta'] = - np.sum(np.multiply(self.budget['delta_ui'], self.budget['delta_duidxj'][:,:,:,1,:]),3)
        self.budget['zadv_delta_delta'] = - np.sum(np.multiply(self.budget['delta_ui'], self.budget['delta_duidxj'][:,:,:,2,:]),3)
        
        # Pressure gradient
        self.budget['dpdx_delta'] = - self.budget['delta_dpdxi'][:,:,:,0]
        self.budget['dpdy_delta'] = - self.budget['delta_dpdxi'][:,:,:,1]
        self.budget['dpdx_delta'] = - self.budget['delta_dpdxi'][:,:,:,2]

        # SGS and Reynolds stress divergence
        for jj in range(3):
            self.budget['xSGS_delta'] += self.budget['delta_dtauijdxk'][:,:,:,0,jj,jj]
            self.budget['ySGS_delta'] += self.budget['delta_dtauijdxk'][:,:,:,1,jj,jj]
            self.budget['zSGS_delta'] += self.budget['delta_dtauijdxk'][:,:,:,2,jj,jj]
        
            self.budget['xturb_delta'] -= self.budget['delta_duiujdxk'][:,:,:,0,jj,jj]
            self.budget['yturb_delta'] -= self.budget['delta_duiujdxk'][:,:,:,1,jj,jj]
            self.budget['zturb_delta'] -= self.budget['delta_duiujdxk'][:,:,:,2,jj,jj]

            self.budget['xturb_mixed'] -= self.budget['ddxk_ui_base_delta_uj'][:,:,:,0,jj,jj] \
                                          + self.budget['ddxk_ui_base_delta_uj'][:,:,:,jj,0,jj]
            self.budget['yturb_mixed'] -= self.budget['ddxk_ui_base_delta_uj'][:,:,:,1,jj,jj] \
                                          + self.budget['ddxk_ui_base_delta_uj'][:,:,:,jj,1,jj]
            self.budget['zturb_mixed'] -= self.budget['ddxk_ui_base_delta_uj'][:,:,:,2,jj,jj] \
                                          + self.budget['ddxk_ui_base_delta_uj'][:,:,:,jj,2,jj]


        
        # Coriolis
        self.budget['xCor_delta'] = -np.sin(np.deg2rad(latitude))*(2/self.pre.Ro)\
                            *(np.sin(np.deg2rad(g_alpha)) - self.budget['delta_ui'][:,:,:,1])
        self.budget['yCor_delta'] = np.sin(np.deg2rad(latitude))*(2/self.pre.Ro)\
                            *(np.cos(np.deg2rad(g_alpha)) - self.budget['delta_ui'][:,:,:,0])
        
        # Buoyancy
        self.budget['zbuoy_delta'] =  -self.budget['delta_T']*(1/(self.prim.Fr*self.prim.Fr*self.prim.Tref))


        if save:
            self.write_budget(budget=1)

    def assemble_budget_2(self, save=False):

        x_pre = self.pre.xLine
        y_pre = self.pre.yLine
        z_pre = self.pre.zLine

        x_prim = self.prim.xLine
        y_prim = self.prim.yLine
        z_prim = self.prim.zLine

        latitude = self.prim.latitude
        g_alpha = -16.072 # this is hardcoded for now

        # initialize final budget variables (that will be saved if files are written)
        keys = [key for key in deficitkey.get_key() if deficitkey.get_key()[key][1] == 1]
        for key in keys:
            self.budget[key] = np.zeros([self.pre.nx, self.pre.ny, self.pre.nz])               


        #####################
        ##### BASE FLOW #####
        #####################

        # mean kinetic energy budget terms
        # TKE loss
        self.budget['MKE_TKE_loss_base'] = np.sum(np.multiply(self.budget['uiuj_base'], self.budget['duidxj_base']), axis=[3,4])

        # Advection
        self.budget['MKE_adv_base'] = np.multiply(self.budget['xadv_base'],self.budget['u_base']) + np.multiply(self.budget['yadv_base'],self.budget['v_base']) \
                                        + np.multiply(self.budget['zadv_base'],self.budget['w_base'])

    def unique_tidx(self):
        """
        Returns the list of tidxs that are the same in the precursor and primary
        """

        pre_tidxs = self.pre.unique_tidx()
        prim_tidxs = self.prim.unique_tidx()
    
        tidxs = [tidx for tidx in prim_tidxs if tidx in pre_tidxs]

        return tidxs

    def write_budget(self, budget):
        print("Writing budget {} files".format(budget))
        base_keys = [key for key in deficitkey.get_key() if deficitkey.get_key()[key][1]==budget and deficitkey.get_key()[key][0]==0]
        
        for key in base_keys:
            term = DeficitIO.key[key][2]
            filename = self.dir_name + '/base_budget{:01d}_term{:02d}_t{:06d}_n{:06d}.s3D'.format(budget, term, self.last_tidx, self.counter)
            self.budget[key].astype('float64').tofile(filename)

        wake_keys = [key for key in deficitkey.get_key() if deficitkey.get_key()[key][1]==budget and deficitkey.get_key()[key][0]==1]
        
        for key in wake_keys:
            term = DeficitIO.key[key][2]
            filename = self.dir_name + '/wake_budget{:01d}_term{:02d}_t{:06d}_n{:06d}.s3D'.format(budget, term, self.last_tidx, self.counter)
            self.budget[key].astype('float64').tofile(filename)


    def read_budgets(self, budget):
        base_keys = [key for key in deficitkey.get_key() if deficitkey.get_key()[key][1]==budget and deficitkey.get_key()[key][0]==0]
        
        for key in base_keys:
            term = DeficitIO.key[key][2]
            filename = self.dir_name + '/base_budget{:01d}_term{:02d}_t{:06d}_n{:06d}.s3D'.format(budget, term, self.last_tidx, self.counter)
            temp = np.fromfile(filename, dtype=np.dtype(np.float64), count=-1) 
            self.budget[key] = temp.reshape((self.pre.nx,self.pre.ny,self.pre.nz))  # reshape into a 3D array

        wake_keys = [key for key in deficitkey.get_key() if deficitkey.get_key()[key][1]==budget and deficitkey.get_key()[key][0]==1]
        
        for key in wake_keys:
            term = DeficitIO.key[key][2]
            filename = self.dir_name + '/wake_budget{:01d}_term{:02d}_t{:06d}_n{:06d}.s3D'.format(budget, term, self.last_tidx, self.counter)
            temp = np.fromfile(filename, dtype=np.dtype(np.float64), count=-1) 
            self.budget[key] = temp.reshape((self.prim.nx,self.prim.ny,self.prim.nz))  # reshape into a 3D array

    #def write_input(self):
    #    """
    #    Writes an input file for the deficit budgets
    #    Called when write_budgets is called
    #    """

        
