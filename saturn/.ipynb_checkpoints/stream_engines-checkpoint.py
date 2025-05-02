import numpy as np
import periodictable as pt

Tref = 298.15 #K
phases = ['solid', 'liquid','gas']
R = 8.314 #kJ/kmolK

def get_molar_mass(mol="CH4"):
    molar_mass=0.0
    elem = []
    for i, v in enumerate(mol):
        if v.isupper():
            if len(elem)!=0:
                molar_mass += getattr(pt,''.join(elem)).mass
                elem = []
            elem.append(v)
        elif v.islower():
            elem.append(v)
        else:
            count=int(v)
            molar_mass += getattr(pt,''.join(elem)).mass*count
            elem = []
    if len(elem)!=0:
        molar_mass += getattr(pt,''.join(elem)).mass
        elem = []
    return(molar_mass)

class line(object):
    def __init__(self,*args,**kwargs):
        kw = {}
        kw.update(kwargs)
        self.mdot = float(kw.get('mdot',0.0))  #total mass flow rate
        self.ndot = float(kw.get('ndot',0.0))  #total molar flow rate
        self.x = np.array(kw.get('x',[]))   #mass fractions
        self.y = np.array(kw.get('y',[]))   #mole fractions
        self.mdoti = np.array(kw.get("mdoti",[]))  #mass flow rates of individual species
        self.ndoti = np.array(kw.get("ndoti",[]))  #molar flow rates of individual species
        self.T = float(kw.get('T',298))   #line temperature
        self.P = float(kw.get('P',101325)) #line pressure
        #self.cp = cp_coeff               # heat capacity
        self.phase = kw.get('phase','gas') 
        self.update_phaseID() # phase ID
        self.rho = kw.get('density',None)  # density
        self.species_list = kw.get('species_list',None)
        assert self.species_list is not None
        self.cp = kw.get('cp_coeff',None)
        self.Mw = np.array([get_molar_mass(elem) for elem in self.species_list])    #molar mass

    def update_mdoti (self):    #update mass flow rates of species using mi = mdot * xi
        self.mdoti = self.mdot*np.array(self.x)
    def update_x(self):
        assert self.mdot != 0;
        self.x = np.array(self.mdoti)/self.mdot
    def update_mdot(self):
        self.mdot = np.sum(np.array(self.mdoti))
    def update_ndot(self):
        self.ndot = np.sum(np.array(self.ndoti))
    def update_ndoti(self):
        self.ndoti = np.array(self.mdoti)/np.array(self.Mw)
    def update_ndoti_from_ndot(self):
        self.ndoti = np.array(np.array(self.ndot)*np.array(self.y))
    def update_mdoti_from_ndoti (self):
        self.mdoti = np.array(self.ndoti)*np.array(self.Mw)
    def update_y (self):
        assert self.ndot != 0.0
        self.y = np.array(self.ndoti)/self.ndot
    def update_phaseID (self):
        self.phaseID=phases.index(self.phase)
    def get_Hline(self):
        self.update_phaseID()
        self.Hi = []
        print("Enthalphy calculated in {} phase".format(phases[self.phaseID]))
        for n,cp in zip(self.ndoti,self.cp):
            buf = R*n*(cp['A'][self.phaseID]*(self.T-Tref)+
                          cp['B'][self.phaseID]*(self.T**2-Tref**2)/2+
                          cp['C'][self.phaseID]*(self.T**3-Tref**3)/3-
                          cp['D'][self.phaseID]*(1/self.T-1/Tref)+
                          cp['E'][self.phaseID]*(self.T**4-Tref**4)/4
                         )
            #print (buf,Tref)
            self.Hi.append(buf)
        self.Hi = np.array(self.Hi)
        self.Hline = np.sum (self.Hi)
        if self.phaseID == 2:   #gas
            self.Mw_avg = np.sum(self.y*self.Mw)
            self.rho_avg = self.P*self.Mw_avg/(R*self.T)
        else:
            pass
            #self.rho_avg = np.sum(self.x*self.rho.T[self.phaseID])
        #self.PrE = self.mdot*self.P/self.rho_avg
    def check(self,verbose=True):
        self.missing_list= {"Total mass flow rate":False,
                       "Total molar flow rate":False,
                       "Species mass flow rate":False,
                       "Species molar flow rate":False,
                       "Mass fractions":False,
                       "Mole fractions":False
                      }
        if self.mdot != 0.0:
            self.missing_list["Total mass flow rate"] = True
        if self.ndot != 0.0:
            self.missing_list["Total molar flow rate"] = True
        if list(self.x):
            self.missing_list["Mass fractions"] = True
        if list(self.y):
            self.missing_list["Mole fractions"] = True
        if list(self.mdoti):
            self.missing_list["Species mass flow rate"] = True
        if list(self.ndoti):
            self.missing_list["Species molar flow rate"] = True
        if verbose:
            for key,val in self.missing_list.items():
                #val = missing_list[key]
                if val:
                    buf = "Updated"
                else:
                    buf = "Missing"
                print (key, ":", buf)
        try:
            if min(self.mdoti) < 0.0 or min(self.ndoti) < 0.0 or min(self.x) < 0.0 or min(self.y) < 0.0:
                print("WARNING: Negative values found in the stream. Please recheck the stream.")
        except:
            pass
        if verbose:
            if abs(np.sum(self.x)-1)>1E-1: 
                print("WARNING: sum of mass fractions is {}".format(np.sum(self.x)))
            if abs(np.sum(self.y)-1)>1E-3: 
                print("WARNING: sum of mole fractions is {}".format(np.sum(self.y)))
            if abs(np.sum(self.mdoti)-self.mdot)>1E-3:
                print("WARNING: Sum of species mass flow rates {} is different from total mass flow rate {}".\
                format(np.sum(self.mdoti),self.mdot))
            if abs(np.sum(self.ndoti)-self.ndot)>1E-1:
                print("WARNING: Sum of species molar flow rates {} is different from total molar flow rate {}".\
                format(np.sum(self.ndoti),self.ndot))
                  
        if not False in list(self.missing_list.values()):
            print ("Line is fully updated.")
    def update_line(self,verbose=False):
        self.check(verbose=False)
        missing_list_old = list(self.missing_list.values())
        #print(self.missing_list)
        if self.missing_list["Total mass flow rate"] and self.missing_list["Mass fractions"]:
            if verbose: print("Total mass flow rate and mass fractions exists")
            self.update_mdoti()
            self.update_ndoti()
            self.update_ndot()
            self.update_y()
            self.check(verbose=True)
        elif self.missing_list["Total molar flow rate"] and self.missing_list["Mole fractions"]:
            if verbose: print("Total molar flow rate and mole fractions exists")
            self.update_ndoti_from_ndot()
            self.update_mdoti_from_ndoti()
            self.update_mdot()
            self.update_x()
            self.check(verbose=True)
        elif self.missing_list["Species molar flow rate"]:
            if verbose: print("Species molar flow rates exists")
            self.update_ndot()
            self.update_y()
            self.update_mdoti_from_ndoti()
            self.update_mdot()
            self.update_x()
            self.check(verbose=True)
        elif self.missing_list["Species mass flow rate"]:
            if verbose: print("Species mass flow rates exists")
            self.update_mdot()
            self.update_x()
            self.update_ndoti()
            self.update_ndot()
            self.update_y()
            self.check(verbose=True)
        else:
            print("ERROR: Line is not defined properly. Start the line by specifying any of the following\
            \n-Total mass flow rate (mdot) and mass fractions(x)\n-Total molar flow rate(ndot) and mole fractions (y)\
            \n-All species mass flow rates (mdoti)\n-All species molar flow rates (ndoti)")
            self.check(verbose=True)
        missing_list_new = list(self.missing_list.values())
        if missing_list_old == missing_list_new and "No" in missing_list_new and verbose:
            print("Unable to update further. Please check if the line is specified correctly.")
        elif False in list(self.missing_list.values()) and verbose:
            print("Line is not fully updated yet. Try again.")
    def update_all(self,verbose = False):
        print("Depracted method. Do not use. Use update_line instead")
        self.check(verbose=False)
        missing_list_old = list(self.missing_list.values())
        #print(missing_list_old)
        try:
            if not self.missing_list["Total mass flow rate"]:
                if verbose: print("Trying to update mdot")
                self.update_mdot()
        except: 
            if verbose: print("Failed")
            else:pass
        try:
            if not self.missing_list["Species mass flow rate"]:
                if verbose: print("Trying to update mdoti")
                self.update_mdoti()
        except:
            if verbose: print("Failed")
            else:pass
        try:
            if not self.missing_list["Mass fractions"]:
                if verbose: print("Trying to update mass fraction")
                self.update_x()
        except: 
            if verbose: print("Failed")
            else:pass
        try:
            if not self.missing_list["Total molar flow rate"]:
                if verbose: print("Trying to update ndot")
                self.update_ndot()
        except: 
            if verbose: print("Failed")
            else:pass
        try:
            if not self.missing_list["Species molar flow rate"]:
                if verbose: print("Trying to update ndoti")
                self.update_ndoti()
        except: 
            if verbose: print("Failed")
            else:pass
        try:
            if not self.missing_list["Species molar flow rate"]:
                if verbose: print("Trying to update ndoti from ndot")
                self.update_ndoti_from_ndot()
        except: 
            if verbose: print("Failed")
            else:pass
        try:
            if not self.missing_list["Species mass flow rate"]:
                if verbose: print("Trying to update mdoti from ndoti")
                self.update_mdoti_from_ndoti()
        except: 
            if verbose: print("Failed")
            else:pass
        try:
            if not self.missing_list["Mole fractions"]:
                if verbose: print("Trying to update y")
                self.update_y()
        except: 
            if verbose: print("Failed")
            else:pass
        self.check(verbose=True)
        missing_list_new = list(self.missing_list.values())
        if missing_list_old == missing_list_new and "No" in missing_list_new:
            print("Unable to update further. Please check if the line is specified correctly.")
        elif False in list(self.missing_list.values()):
            print("Line is not fully updated yet. Try again.")
    def react(self, stoich=None,conversion=None):
        assert stoich is not None
        assert conversion is not None
        stoich = np.array(stoich)
        limiting_react_ndot = self.ndoti[conversion[0]]
        normal_stoich = stoich/abs(stoich[conversion[0]])
        delta_ndoti = normal_stoich*conversion[1]*limiting_react_ndot
        self.ndoti = self.ndoti + delta_ndoti
            