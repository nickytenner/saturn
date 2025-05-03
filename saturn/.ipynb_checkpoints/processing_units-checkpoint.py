import numpy as np
import periodictable as pt


class unit(object):
    """
    This is to declare a process unit
    options:
        inlet: specify all inlet lines. eg: inlet = [line1,line2]
        outlet: specify all outlet lines. eg: outlet = [line3,line4]
    """
    def __init__(self,*args,**kwargs):
        kw = {}
        kw.update(kwargs)
        self.inlet = list(kw.get('inlets',[]))  #inlet lines
        self.outlet = list(kw.get('outlets',[]))  #outlet lines
    def mass_balance_check (self):
        mass_in = np.sum(np.array([getattr(line,"mdot") for line in self.inlet]))
        mass_out = np.sum(np.array([getattr(line,"mdot") for line in self.outlet]))
        if abs(mass_in - mass_out) > 1E-10:
            print("Mass balance around the unit failed by {}".format(mass_in-mass_out))
            print("mass in: {} mass out: {}".format(mass_in,mass_out))
        else:
            print("Mass balance check around the unit: success")
    def update_outlines(self):
        for line in self.outlet:
            getattr(line,"update_line")()
            getattr(line,"get_Hline")()

        
        
