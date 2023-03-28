
"""
Module containing a bi-directional has thable for the mapping and inverse-mapping
for budget terms in index-tuple form (flow #, budget #, term #) based off of the 
PadeOps budget definitions and common form as defined in get_key()

Kerry Klemmer
March 2023
(borrowed from budgetkey.py)
"""

class bidict(dict):
    """
    Bi-directional dictionary, courtesy of: 
    https://stackoverflow.com/questions/3318625/how-to-implement-an-efficient-bidirectional-hash-table

    Used to look up a mapping and inverse mapping (not necessarily unique) between key parings. 
    """
    def __init__(self, *args, **kwargs):
        super(bidict, self).__init__(*args, **kwargs)
        self.inverse = {}
        for key, value in self.items():
            self.inverse.setdefault(value,[]).append(key) 

    def __setitem__(self, key, value):
        if key in self:
            self.inverse[self[key]].remove(key) 
        super(bidict, self).__setitem__(key, value)
        self.inverse.setdefault(value,[]).append(key)        

    def __delitem__(self, key):
        self.inverse.setdefault(self[key],[]).remove(key)
        if self[key] in self.inverse and not self.inverse[self[key]]: 
            del self.inverse[self[key]]
        super(bidict, self).__delitem__(key)


def key_search_r(nested_dict, key):
    """
    Performs a recursive search for the dictionary key `key` in any of the dictionaries contained 
    inside `nested_dict`. Returns the value of nested_dict[key] at the first match. 
    
    Parameters
    ----------
    nested_dict (dict-like) : dictionary [possibly of dictionaries]
    key (str) : dictionary key to match
    
    Returns
    -------
    nested_dict[key] if successful, None otherwise. 
    """
    
    try: 
        for k in nested_dict.keys(): 
            if k == key: 
                return nested_dict[k]
            else: 
                res = key_search_r(nested_dict[k], key)
                if res is not None: 
                    return res
        
    except AttributeError as e: 
        return
    
    
def get_key(): 
    """
    Returns a bidirectional hash table between colloquial string values for budget terms and 
    a tuple-look up ordered (flow #, budget #, term #), where flow # is 0 for base flow and
    1 in the wake deficit, and budget and terms are as defined in e.g. budget_time_avg.F90. 
    """
    key = {  # FLOW 0: BASE FLOW
        # BUDGET 0 TERMS: (1st and second order averages, scalars excluded)
        'u_base': (0, 0, 1), 
        'v_base': (0, 0, 2), 
        'w_base': (0, 0, 3), 
        'uu_base': (0, 0, 4), 
        'uv_base': (0, 0, 5), 
        'uw_base': (0, 0, 6), 
        'vv_base': (0, 0, 7), 
        'vw_base': (0, 0, 8), 
        'ww_base': (0, 0, 9), 
        'p_base': (0, 0, 10), 
        'tau11_base': (0, 0, 11), 
        'tau12_base': (0, 0, 12), 
        'tau13_base': (0, 0, 13), 
        'tau22_base': (0, 0, 14), 
        'tau23_base': (0, 0, 15), 
        'tau33_base': (0, 0, 16), 
        'pu_base': (0, 0, 17), 
        'pv_base': (0, 0, 18), 
        'pw_base': (0, 0, 19), 
        #'T_base': (0, 20), 
        #'uT_base': (0, 21), 
        #'vT_base': (0, 22), 
        #'wT_base': (0, 23), 
        #'TT_base': (0, 24), 
        #'theta': (0, 31), 
        # BUDGET 1 TERMS: (momentum)
        'xadv_base': (0, 1, 1),  # x-advection
        'dpdx_base': (0, 1, 2),  # x-pressure gradient
        'xSGS_base': (0, 1, 3),  # x-sub grid stresses
        'xturb_base': (0, 1, 4),   # x-Actuator disk
        'xCor_base': (0, 1, 5), # x-coriolis
        'yadv_base': (0, 1, 6),  
        'dpdy_base': (0, 1, 7), 
        'ySGS_base': (0, 1, 9), 
        'yturb_base': (0, 1, 10), 
        'yCor_base': (0, 1, 11), 
        'zadv_base': (0, 1, 12),  
        'dpdz_base': (0, 1, 13), 
        'zSGS_base': (0, 1, 14), 
        'zturb_base': (0, 1, 15), 
        'zbuoy_base': (0, 1, 16),
        # BUDGET 2 TERMS: (MKE) 
        'MKE_TKE_loss_base': (0, 2, 1), 
        'MKE_adv_base': (0, 2, 2), 
        'MKE_turb_transport_base': (0, 2, 3), 
        'MKE_p_transport_base': (0, 2, 4), 
        'MKE_SGS_transport_base': (0, 2, 5), 
        'MKE_dissipation_base': (0, 2, 6), 
        'MKE_AD_base': (0, 2, 7), 
        'MKE_geostrophic_base': (0, 2, 8), 
        'MKE_coriolis_base': (0, 2, 9),
        # FLOW 1: WAKE DEFICIT
        # BUDGET 0 TERMS: (1st and second order averages, scalars excluded)
        'delta_u': (1, 0, 1), 
        'delta_v': (1, 0, 2), 
        'delta_w': (1, 0, 3), 
        'delta_uu': (1, 0, 4), 
        'delta_uv': (1, 0, 5), 
        'delta_uw': (1, 0, 6), 
        'delta_vv': (1, 0, 7), 
        'delta_vw': (1, 0, 8), 
        'delta_ww': (1, 0, 9), 
        'delta_p': (1, 0, 10), 
        'delta_tau11': (1, 0, 11), 
        'delta_tau12': (1, 0, 12), 
        'delta_tau13': (1, 0, 13), 
        'delta_tau22': (1, 0, 14), 
        'delta_tau23': (1, 0, 15), 
        'delta_tau33': (1, 0, 16),
        'u_base_delta_u': (1, 0, 17), 
        'u_base_delta_w': (1, 0, 18), 
        'w_base_delta_u': (1, 0, 19), 
        'v_base_delta_v': (1, 0, 20), 
        'v_base_delta_w': (1, 0, 21), 
        'w_base_delta_v': (1, 0, 22), 
        'w_base_delta_w': (1, 0, 23), 
        'delta_pu': (1, 0, 24), 
        'delta_pv': (1, 0, 25), 
        'delta_pw': (1, 0, 26), 
        'delta_p_u_base': (1, 0, 27), 
        'delta_p_v_base': (1, 0, 28), 
        'delta_p_w_base': (1, 0, 29),  
        'p_base_delta_u': (1, 0, 30), 
        'p_base_delta_v': (1, 0, 31), 
        'p_base_delta_w': (1, 0, 32),        
         # BUDGET 1 TERMS: (momentum)
        'xadv_delta_base': (1, 1, 1),  # x-advection
        'xadv_delta_delta': (1, 1, 2),  # x-advection
        'dpdx_delta': (1, 1, 3),  # x-pressure gradient
        'xSGS_delta': (1, 1, 4),  # x-sub grid stresses
        'xturb_delta': (1, 1, 5),   # x-Actuator disk
        'xturb_mixed': (1, 1, 7),   # x-Actuator disk
        'xCor_delta': (1, 1, 8), # x-coriolis
        'yadv_delta_base': (1, 1, 9),  # x-advection
        'yadv_delta_delta': (1, 1, 10),  # x-advection
        'dpdy_delta': (1, 1, 11),  # x-pressure gradient
        'ySGS_delta': (1, 1, 12),  # x-sub grid stresses
        'yturb_delta': (1, 1, 13),   # x-Actuator disk
        'yturb_mixed': (1, 1, 14),   # x-Actuator disk
        'yCor_delta': (1, 1, 15), # x-coriolis
        'zadv_delta_base': (1, 1, 16),  # x-advection
        'zadv_delta_delta': (1, 1, 17),  # x-advection
        'dpdz_delta': (1, 1, 18),  # x-pressure gradient
        'zSGS_delta': (1, 1, 19),  # x-sub grid stresses
        'zturb_delta': (1, 1, 20),   # x-Actuator disk
        'zturb_mixed': (1, 1, 21),   # x-Actuator disk
        'zbuoy_delta': (1, 1, 22)
        }

    
    return bidict(key)

