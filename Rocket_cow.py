import cantera as ct
import numpy as np
import scipy as sp
import pandas as pd
from pint import UnitRegistry
import plotly.express as px
import plotly.graph_objects as go
# import yaml
# import fluids
import streamlit as st
# import time
ureg = UnitRegistry()
state = st.session_state

# initialize buttons
if 'run_button' not in state:
    state.run_button = False

if 'plot_button' not in state:
    state.plot_button = False

if 'contour_button' not in state:
    state.contour_button = False

def run_button():
    state.run_button = True  

def plot_button():
    state.plot_button = True

def contour_button():
    state.contour_button = True  

class EngineState:
    '''
    This class will define the properties of an rocket engine at the chamber, throat and exit 
    given the propellant, and conditions that engine is in.

    Parameters
    --------------------------------------------------
    oxidizer : string
        A string of a species name which is used to define a cantera solution representing the oxidizer.
        Define based off the list of phase names from the propellants.yaml file TODO: kill propellants.yaml and create a list users can choose from 
    fuel : string
        Similar to oxidizer, this is a string of a species name which is used to define a cantera solution representing the fuel.
        Define based off the list of phase names from the propellants.yaml file TODO: kill propellants.yaml and create a list users can choose from 
    of_ratio : float
        the mass ratio of oxidizer to fuel defined as mass_oxidizer / mass_fuel.
    pressure : float
        The pressure of the combustion chamber in units of Pascals (Pa).
    exit_value : float
        The exit conditions of the engine. Depending on what 'exit parameter' is
        defined to be the exit condition can be defined in the following ways:
            'pressure' (default): Exit pressure of the gases in Pa should be greater than 0 but lower than chamber pressure
            'area ratio': Defined as area_exit / area_throat should be value greater than 1 TODO: add functionality
    exit_parameter : string
        A string which determines how you are defining the exit condition of your engine. Can either be set to 'pressure' (default)
        or 'area ratio'. 
    temp_oxidizer : float
        A positive float value defining the temperature of the oxidizer in Kelvin (K). If the oxidizer is liquid the temperature 
        is assumed to be saturated temperatures at standard pressure, and gaseous oxidizer will be set to a room tempearture of 295.15K. TODO: add functionality
    temp_fuel : float
        Similar to above. positive float value defining the temperature of the fuel in Kelvin (K). If the fuel is liquid the temperature 
        is assumed to be saturated temperatures at standard pressure, and gaseous oxidizer will be set to a room tempearture of 295.15K. TODO: add functionality
    transport_species Optional: [list | None]
        A list of species 
    assumption : str
        The assumption determines the state of the combustion products and whether their composition is in equilibrium, or fixed/frozen:
            'equilibrium' (default): Composition is always in equilibrium
            'frozen': Composition is fixed after initial combustion in the chamber. 

        

    TODO board: 
    --------------------------------------------------
        TODO: apply efficiency.
        TODO: apply flow seperation
        TODO: improve chemical database
        TODO: allow fluid temp setting
        TODO: condensed species capabilities
        TODO: Improve convergence *especially for locations along the length of the chamber*
        TODO: Ambient pressure
        TODO: Ambient pressure Array
        TODO: add Pint https://pint.readthedocs.io/en/stable/advanced/wrapping.html#wrapping
        TODO: add FAC capabilities
        TODO: add frozen

    Attributes:
    --------------------------------------------------
    Properties : Pandas Dataframe
    
    Methods
    --------------------------------------------------

    A class for the state of a rocket engine. Given the propellants and conditions of a rocket engine this 
    class will define the properties of an engine at
    '''

    def __init__(self, oxidizer, fuel, of_ratio, pressure, exit_value, gas, condensate = None, size_value=None, assumption = 'equilibrium', size_parameter ="thrust",  
                 exit_parameter="pressure", transport=None, throat_inlet_radius_ratio=1.5, temperature_wall=500):
        '''
        Initializes instance of class. See class description for details. 
        '''
        # initializes propellants
        self.assumption = assumption
        self._of_ratio = of_ratio
        self._pressure = pressure
        oxidizer.TP = oxidizer.T, pressure
        oxidizer.equilibrate('TP')
        self._oxidizer = oxidizer
        # self._oxidizer = self.__clone_solution(oxidizer)
        fuel.TP = fuel.T, pressure
        fuel.equilibrate('TP')
        self._fuel = fuel
        self._gas = gas
        self._condensate = condensate
        self._transport = transport
        # self._fuel = self.__clone_solution(fuel)
        # self.chamber_gas = self.__clone_solution(gas)
        # self.throat_gas = self.__clone_solution(gas)
        # self.exit_gas = self.__clone_solution(gas)
        # self._gas = self.__clone_solution(gas)

        # if condensate:
        #     self._condensate = self.__clone_solution(condensate)
        #     self.exit_condensate = self.__clone_solution(condensate)
        #     self.throat_condensate = self.__clone_solution(condensate)
        #     self.chamber_condensate = self.__clone_solution(condensate)
        # else:
        #     self._condensate = None
        #     self.exit_condensate = None
        #     self.throat_condensate = None
        #     self.chamber_condensate = None
        # if transport: 
        #     self._transport = self.__clone_solution(transport)
        #     self.exit_transport = self.__clone_solution(transport)
        #     self.throat_transport = self.__clone_solution(transport)
        #     self.chamber_transport = self.__clone_solution(transport)
        # else: 
        #     self._transport = None
        #     self.exit_transport = None
        #     self.throat_transport = None
        #     self.chamber_transport = None

    
        # finds chamber, throat, and exit properties
        self.__chamber_properties(self._gas, self._condensate, self._transport)
        self.__throat_properties(self._gas, self._condensate, self._transport)
        self.__exit_properties(self._gas, self._condensate, self._transport, exit_value, exit_parameter=exit_parameter)

        self.throat_inlet_radius_ratio = throat_inlet_radius_ratio # this value is not needed for sizing but is necessary to find heat tranfer coeff. via bartz
        self.temperature_wall = temperature_wall
        
        # sizes engine
        if  size_value:
            self.size_engine(size_value, size_parameter=size_parameter)


        # self.engine_state_dict = {'chamber' : self.chamber, 'throat' : self.throat, 'exit' : self.exit }
        self.engine_state = pd.DataFrame([self.chamber, self.throat, self.exit], index=["chamber", "throat", "exit"])

    def __call__(self):
        return self.engine_state
    
    def __str__(self):
        return self.engine_state.to_string()
    
    # def __clone_solution(self, sol: ct.Solution) -> ct.Solution:
    #     """Return a deep, independant copy of a Cantera Solution"""
    #     new = ct.Solution(thermo=sol.thermo_model, species=sol.species())
    #     if sol.transport_model != 'none': 
    #         new.transport_model = sol.transport_model
    #     new.TPX = sol.T, sol.P, sol.X
    #     return new

    def __get_thermo_derivatives(self, gas, condensate):
        '''
        This is an internal method not meant for use outside the class.
        given a mixture of gases & condensed species this method will find thermodynamic derivatives.
        These derivatives can then be used to calculate thermodynamic properties, or create interpolation/extrapolation of engine properties 

        The theory used for these calculations can be found in section 2.5 and 2.6 of RP-1311.
        The equations used are 2.50, 2.51, 2.56, 2.57, 2.58, 2.64, 2.65, 2.66 

        TODO: add condensed species to implementation
        TODO: explore interpolation of engine properties using derivatives
        
        Parameters
        --------------------------------------------------
        mixture : Cantera Mixture
            A solution or mixture of primarily gases (ideal gas assumptions will hold up to several percent condensed species by mass)
        
        Returns
        --------------------------------------------------
        Derivatives : dictionary:
            A dictionary consisting of the derivatives this function calculates. The keys and descritions of said dictionaries are as follows
        
            dpi_dlnT_P : list
                A list of the derivative of pi with respect to ln(T) at constant pressure for each element, where pi is -lambda/RT, 
                and lambda is the langrangian multiplier
            dlnn_dlnT_P : float
                The derivative of ln(n) with respect to ln(T) at constant pressure. Where n is the number of moles of gas
            dpi_dlnP_T : list
                A list of the derivative of pi with respect to ln(P) at constant temperature for each element.
            dlnn_dlnP_T : float
                The derivative of ln(n) with respect to ln(P) at constant tempearture.
            dlnV_dlnT_P : float
                The derivative of ln(V) with respect to ln(T) at constant pressure.
            dlnV_dlnP_T : float
                The derivative of ln(V) with respect to ln(P) at constant tempearture.
        '''

        if self.assumption == 'equilibrium':             
            # Defines the number of moles of each species in the mixture
            moles = gas.X * (1/ gas.mean_molecular_weight)
            num_variables = 2 * gas.n_elements + 2

            # Initializing Solution Matrices table 2.3 and 2.4 in RP-1311
            coeff_matrix = np.zeros((num_variables, num_variables))
            right_hand_side = np.zeros(num_variables)

            # Initializes a_ij
            a_ij = np.zeros((gas.n_elements, gas.n_species))
            for i, element in enumerate(gas.element_names): 
                for j, species in enumerate(gas.species_names):
                    a_ij[i,j] = gas.n_atoms(species, element)


            # Coefficients for equation 2.56 TODO: add terms for condensed species
            for k in range(gas.n_elements):
                for i in range(gas.n_elements):
                    coeff_matrix[k,i] = np.sum(a_ij[k,:] * a_ij[i,:] * moles)
                coeff_matrix[k, gas.n_elements] = np.sum(a_ij[k,:] * moles)
                right_hand_side[k] = -np.sum(a_ij[k,:] * moles * gas.standard_enthalpies_RT)

            # TODO add equation 2.57 (it is for condensed species)
            # for i in range():
            #     None


            # Coefficients for equation 2.58 
            for i in range(gas.n_elements):
                coeff_matrix[gas.n_elements, i] = np.sum(a_ij[i, :] * moles)
            right_hand_side[gas.n_elements] = -np.sum(moles * gas.standard_enthalpies_RT)

            # Coefficients for equation 2.64 TODO: add terms for condensed species
            for k in range(gas.n_elements):
                for i in range(gas.n_elements):
                    coeff_matrix[gas.n_elements+1+k,gas.n_elements+1+i] = np.sum(a_ij[k,:] * a_ij[i,:] * moles)
                coeff_matrix[gas.n_elements+1+k, 2*gas.n_elements+1] = np.sum(a_ij[k,:] * moles)
                right_hand_side[gas.n_elements+1+k] = np.sum(a_ij[k,:] * moles)

            # TODO add equation 2.65 (it is for condensed species)
            

            # Coefficeints for equation 2.66
            for i in range(gas.n_elements):
                coeff_matrix[2*gas.n_elements+1, gas.n_elements+1+i] = np.sum(a_ij[i, :] * moles)
            right_hand_side[2*gas.n_elements+1] = np.sum(moles)

            # Solve for the derivatives define them based off table 2.3, 2.4 and equation 2.50 and 2.51 
            derivs = np.linalg.solve(coeff_matrix, right_hand_side)
            derivatives = { "dpi_dlnT_P"    : derivs[0 : gas.n_elements], 
                            "dlnn_dlnT_P"   : derivs[gas.n_elements], 
                            "dpi_dlnP_T"    : derivs[gas.n_elements + 1: 2 * gas.n_elements + 1],
                            "dlnn_dlnP_T"   : derivs[2 * gas.n_elements + 1], 
                            "dlnV_dlnT_P"   : 1 + derivs[gas.n_elements], 
                            "dlnV_dlnP_T"   : -1 + derivs[2 * gas.n_elements + 1]}

        elif self.assumption == 'frozen': 
            derivatives = { "dpi_dlnT_P"    : [np.nan]*gas.n_elements, 
                            "dlnn_dlnT_P"   : 0,
                            "dpi_dlnP_T"    : [np.nan]*gas.n_elements,
                            "dlnn_dlnP_T"   : 0, 
                            "dlnV_dlnT_P"   : 1, 
                            "dlnV_dlnP_T"   : -1}
        else: 
            raise ValueError("Assumption must be 'frozen' or 'equilibrium'")

        return derivatives

    def __get_thermo_properties(self, gas, condensate, dpi_dlnT_P, dlnn_dlnT_P, dlnV_dlnT_P, dlnV_dlnP_T):
        '''
        This is an internal method not meant for use outside the class.
        given a mixture of gases & condensed species, as well as certain thermdynamic derivatives of said mixture, this function will find 
        the thermodynamic properties of the mixture

        The theory used for these calculations can be found in section 2.5 and 2.6 of RP-1311.

        TODO: add condensed species to implementation
        TODO: add capability to solve for thermal conductivity, viscocity, and prandtl number  

        Parameters
        --------------------------------------------------
        mixture : Cantera Mixture
            A solution or mixture of primarily gases (ideal gas assumptions will hold up to several percent condensed species by mass)
        dpi_dlnT_P : list
            A list of the derivative of pi with respect to ln(T) at constant pressure for each element, where pi is -lambda/RT, 
            and lambda is the langrangian multiplier
        dlnn_dlnT_P : float
            The derivative of ln(n) with respect to ln(T) at constant pressure. Where n is the number of moles of gas.
        dlnV_dlnT_P : float
            The derivative of ln(V) with respect to ln(T) at constant pressure.
        dlnV_dlnP_T : float
            The derivative of ln(V) with respect to ln(P) at constant tempearture.
            
        Returns
        --------------------------------------------------
        properties : dictionaries
            A dictionary consisting of the properties this function calculates. The keys and descritions of said dictionaries are as follows

            Pressure : float
            Temperature : float
            density : float
            specific_volume : float
            enthalpy : float
            internal_energy : float
            gibbs : float
            entropy : float
            molar_mass : float
            c_p : float
                the specific heat at constant pressure
            c_v : float
                the specific heat at constant volume
            gamma : float
                the specific heat ratio
            gamma_s : float
                defined as derivative ln(P) with respect to ln(rho) at constant entropy per equation 2.71 in RP-1311
            speed_sound : float
                the speed of sound in the mixture
        '''

        if self.assumption == 'equilibrium':
            # Defines the number of moles of each species in the micture
            moles = gas.X * (1/ gas.mean_molecular_weight)

            # Initializes a_ij
            a_ij = np.zeros((gas.n_elements, gas.n_species))
            for i, element in enumerate(gas.element_names): 
                for j, species in enumerate(gas.species_names):
                    a_ij[i,j] = gas.n_atoms(species, element)

            # Finds specific heat at constant pressure based on equation 2.59 TODO: add terms for condensed species
            c_p =  ct.gas_constant * (
                np.sum([dpi_dlnT_P[i] * np.sum(a_ij[i,:] * moles * gas.standard_enthalpies_RT) for i in range(gas.n_elements)]) +
                np.sum(moles * gas.standard_enthalpies_RT) * dlnn_dlnT_P +
                np.sum(moles * gas.standard_cp_R) +
                np.sum(moles * gas.standard_enthalpies_RT**2)
            )

            # Finds specifc heat at constant volume, based on equation 2.70, specific heat ratio, The isentropic exponent based on equation 2.73, 
            # and speed of sound based on equation 2.74
            c_v = c_p + gas.P * gas.v / gas.T * dlnV_dlnT_P**2 / dlnV_dlnP_T
            gamma = c_p / c_v
            gamma_s = -gamma/dlnV_dlnP_T
            speed_sound = np.sqrt(ct.gas_constant * gas.T * gamma_s/gas.mean_molecular_weight)
        
        elif self.assumption == 'frozen':
            c_p = gas.cp
            c_v = gas.cv
            gamma = c_p / c_v
            gamma_s = gamma
            speed_sound = np.sqrt(ct.gas_constant * gas.T * gamma_s/gas.mean_molecular_weight)
        
        else: 
            raise ValueError("Assumption must be 'frozen' or 'equilibrium'")

        properties = {  "of ratio"          : self._of_ratio,
                        "pressure"          : gas.P,
                        "temperature"       : gas.T,
                        "density"           : gas.density_mass,
                        "specific volume"   : gas.volume_mass,
                        "enthalpy"          : gas.enthalpy_mass,
                        "internal energy"   : gas.int_energy_mass,
                        "gibbs"             : gas.gibbs_mass,
                        "entropy"           : gas.entropy_mass,
                        "molar mass"        : gas.mean_molecular_weight,
                        "c_p"               : c_p,
                        "c_v"               : c_v,
                        "gamma"             : gamma,                                                                            
                        "gamma_s"           : gamma_s,    
                        "speed sound"       : speed_sound          
                        }

        return properties
    
    def __get_transport_properties(self, gas, transport):
        '''
        TODO: Update  comment
        TODO: find equilibrium thermal conductivity
        '''
        if transport:
            transport.TPX = gas.T, gas.P, {k: v for k, v in zip(gas.species_names, gas.X) if k in transport.species_names}        
            # frozen transport properties
            properties = {  "viscosity"             : transport.viscosity,                                                  # frozen viscosity in Pa*s
                            "thermal conductivity"  : transport.thermal_conductivity,                                       # thermal conductivity in W/m/K
                            "prandtl number"        : transport.cp * transport.viscosity / transport.thermal_conductivity   # frozen prandtl number
                            }
        else: 
            properties = {  "viscosity"             : np.nan,
                            "thermal conductivity"  : np.nan,
                            "prandtl number"        : np.nan
                            }
        # TODO: implement reacting transport properties
        

        return properties, transport.mole_fraction_dict()

    def _heat_flux(self, temperature, pressure, transport_mole_fractions, area, mach, viscosity_exponent = 0.6):

        if not self.temperature_wall:
            # print(  "temperature_wall is not defined. If you want heat transfer coeff. try initiating class with temperature wall defined, or " \
                    # "setting [object_name].temperature_wall to chosen value.")
            heat_prop = {'heat transfer coefficient': np.nan,
                         'heat flux': np.nan} 
        else: 
            T_ref = temperature * (1+0.032* mach**2 + 0.58*(self.temperature_wall/temperature - 1))
            self._transport.TPX = T_ref, pressure, transport_mole_fractions

            prandtl_number = self._transport.viscosity * self._transport.cp / self._transport.thermal_conductivity
            d_t = self.throat['diameter']
            r_t_inlet = self.throat_inlet_radius_ratio*d_t/2
            m = viscosity_exponent #exponent of temperature such that mu proportional to T^m

            sigma = 1/ ( 1/2*self.temperature_wall/self.chamber['temperature'] * (1+ (self.chamber['gamma']-1)/2 * mach**2)+ 1/2)**(0.8-m/5) / (1+(self.chamber['gamma']-1)/2 * mach**2)**(m/5)
            
            heat_transfer_coeff =   (0.026/d_t**0.2 * (self._transport.viscosity**0.2*self.chamber['c_p']/prandtl_number**0.6) 
                                    * (self.chamber['pressure']/self.chamber['c*'])**0.8* (d_t/r_t_inlet)**0.1 * ((np.pi*d_t**2/4)/area)**0.9* sigma)

            T_aw = self.chamber['temperature'] * (1+ prandtl_number**0.33*(self.chamber['gamma']-1)/2*mach**2)/(1+(self.chamber['gamma']-1)/2*mach**2)

            heat_prop = {'heat transfer coefficient': heat_transfer_coeff,
                         'heat flux': heat_transfer_coeff * (T_aw-self.temperature_wall)} 

        return heat_prop

    def __chamber_properties(self, gas, condensate, transport):
        '''
        This is an internal method not meant for use outside the class.
        given a mixture of gases & condensed species, as well as certain thermdynamic properties of said mixture, this function will find 
        properties of the engine

        The theory used for these calculations can be found in section 2.5 and 2.6 of RP-1311.

        Parameters
        --------------------------------------------------
        products : Cantera Mixture 
            A solution of primarily gases (ideal gas assumptions will hold up to several percent condensed species by mass)
            
        Returns
        --------------------------------------------------
        properties : dictionary
        '''

        # Equilibriates chamber returning a cantera mixture with the properties & composition at the chamber
        molar_ratio = self._of_ratio / (self._oxidizer.mean_molecular_weight / self._fuel.mean_molecular_weight)
        moles_ox = molar_ratio / (1 + molar_ratio)
        moles_f = 1 - moles_ox

        chamber_mixture = ct.Mixture([(self._fuel, moles_f), (self._oxidizer, moles_ox), (gas, 0)])
        chamber_mixture.equilibrate('HP', estimate_equil=-1)
        # Finds thermodynamic derivatives 
        derivatives = self.__get_thermo_derivatives(gas, condensate)

        # Finds thermodynamic properties
        therm_prop = self.__get_thermo_properties(gas, condensate, derivatives["dpi_dlnT_P"], derivatives["dlnn_dlnT_P"], derivatives["dlnV_dlnT_P"], derivatives["dlnV_dlnP_T"])
        
        # Finds tranpsort properties
        transport_prop, transport_moles = self.__get_transport_properties(gas, transport)

        # Calculate c* per
        char_velocity = (np.sqrt(ct.gas_constant * therm_prop["temperature"] / (therm_prop["molar mass"] * therm_prop["gamma"])) * 
                            np.power(2 / (therm_prop["gamma"] + 1), -(therm_prop["gamma"] + 1) / (2*(therm_prop["gamma"] - 1))))

        # velocity and Mach are 0 in a FAC combustor, area ratio, Isp, Ivac, and Cf are not defined at chamber
        chamber_prop = {"velocity"                  : 0, 
                        "mach"                      : 0, 
                        "area ratio"                : np.nan,
                        "I_sp"                      : np.nan,
                        "I_vac"                     : np.nan,
                        "c*"                        : char_velocity,
                        "C_f"                       : np.nan,
                        "mole fraction"             : gas.mole_fraction_dict(),
                        "mole fraction transport"   : transport_moles}
        
        self.chamber = therm_prop | transport_prop | chamber_prop | derivatives

    def __throat_properties(self, gas, condensate, transport):
        '''
        Description
        --------------------------------------------------
        This is an internal method not meant for use outside the class.
        given a solution of gases & condensed species, this function will find 
        properties of the engine

        The theory used for these calculations can be found in section 2.5 and 2.6 of RP-1311.

        Parameters
        --------------------------------------------------
        products : Cantera Solution 
            A solution of primarily gases (ideal gas assumptions will hold up to several percent condensed species by mass)
            
        Returns
        --------------------------------------------------
        None
        '''
        chamber = self.chamber
        # initial guess at throat pressure using specific heat ratio gamma
        pressure_throat = chamber["pressure"] / np.power((chamber["gamma_s"] + 1) / 2., chamber["gamma_s"] / (chamber["gamma_s"] - 1))

        # Setting up for iteration
        max_iter_throat = 10            # NOTE exceeds value of 4 from RP-1311
        tolerance_throat = 4e-5
        mach = 1.0
        num_iter = 0
        residual = 1

        gas.SPX = chamber["entropy"], pressure_throat, self.chamber["mole fraction"]
 
        while (residual > tolerance_throat and num_iter < max_iter_throat ) :
            num_iter += 1
            gas.SP = chamber["entropy"], pressure_throat

            if self.assumption == 'equilibrium': 
                gas.equilibrate('SP')   
            elif self.assumption == 'frozen':
                pass
            else:
                raise ValueError("Assumption must be 'frozen' or 'equilibrium'")
            
            throat_derivatives = self.__get_thermo_derivatives(gas, condensate)
            throat_properties = self.__get_thermo_properties(gas, condensate, throat_derivatives["dpi_dlnT_P"], throat_derivatives["dlnn_dlnT_P"], throat_derivatives["dlnV_dlnT_P"], throat_derivatives["dlnV_dlnP_T"])
            
            velocity = np.sqrt(2 * (chamber["enthalpy"] - throat_properties["enthalpy"]))
            speed_sound = np.sqrt(ct.gas_constant * throat_properties["temperature"] * throat_properties["gamma_s"]  / throat_properties["molar mass"])
            mach = velocity / speed_sound
            
            pressure_throat = pressure_throat * (1 + throat_properties["gamma_s"] * mach**2) / (1 + throat_properties["gamma_s"] )

            residual = np.abs((velocity**2 - speed_sound**2)/velocity**2)

        if num_iter >= max_iter_throat:
            print(f'Warning: Convergance took {num_iter} iterations which exceeds the limit of {max_iter_throat} max iterations. residual is {residual} which exceeds tolerance of {tolerance_throat}.')
        
        transport_prop, transport_moles = self.__get_transport_properties(gas, transport)
        
        # velocity and Mach are the same at throat, area ratio, Isp, Ivac, and Cf are not defined at chamber
        throat_prop = { "velocity"                  : speed_sound, 
                        "mach"                      : 1, 
                        "area ratio"                : np.nan,
                        "I_sp"                      : np.nan,
                        "I_vac"                     : np.nan,
                        "c*"                        : np.nan,
                        "C_f"                       : np.nan,
                        "mole fraction"             : gas.mole_fraction_dict(),
                        "mole fraction transport"   : transport_moles}

        self.throat = throat_properties | transport_prop | throat_prop | throat_derivatives

    def __exit_properties(self, gas, condensate, transport, exit_value, exit_parameter='pressure'):
        '''
        This is an internal method not meant for use outside the class.
        equilbriates solution and finds thermal derivatives and properties at the exit of the engine.
      
        Parameters
        --------------------------------------------------
        products : Cantera Solution 
            A solution of primarily gases (ideal gas assumptions will hold up to several percent condensed species by mass)
        exit_value : float
            The exit conditions of the engine. Depending on what 'exit parameter' is
            defined to be the exit condition can be defined in the following ways:
                'pressure' (default): Exit pressure of the gases in Pa should be greater than 0 but lower than chamber pressure
                'area ratio': Definfed as area_exit / area_throat should be value greater than 1 TODO: add functionality
        exit_parameter : string
            A string which determines how you are defining the exit condition of your engine. Can either be set to 'pressure' (default)
            or 'area ratio'. 
            
        Returns
        --------------------------------------------------
        None
        '''

        if exit_parameter == "area ratio":
            exit_properties = self.state_at_area(gas, condensate, transport, exit_value, speed = "supersonic")
            ae_mdot = 1/(exit_properties["density"]*exit_properties["velocity"])
            exit_properties["I_sp"] = exit_properties["velocity"] / sp.constants.g
            exit_properties["I_vac"] = exit_properties["velocity"] / sp.constants.g + exit_properties["pressure"] * ae_mdot / sp.constants.g
            exit_properties["C_f"] = exit_properties["velocity"] / self.chamber["c*"]

            self.exit =  exit_properties
        
        elif exit_parameter == "pressure":
            pressure = exit_value
            gas.SPX = self.chamber["entropy"], pressure, self.throat['mole fraction']
            if self.assumption == 'equilibrium': 
                gas.equilibrate('SP')   
            elif self.assumption == 'frozen':
                pass
            else:
                raise ValueError("Assumption must be 'frozen' or 'equilibrium'")

            exit_derivatives = self.__get_thermo_derivatives(gas, condensate)
            exit_properties = self.__get_thermo_properties(gas, condensate, exit_derivatives["dpi_dlnT_P"], exit_derivatives["dlnn_dlnT_P"], exit_derivatives["dlnV_dlnT_P"], exit_derivatives["dlnV_dlnP_T"])

            velocity = np.sqrt(2* (self.chamber["enthalpy"] - exit_properties["enthalpy"]))
            speed_sound = np.sqrt(ct.gas_constant * exit_properties["temperature"] * exit_properties["gamma_s"]  / exit_properties["molar mass"])

            at_mdot = 1 / (self.throat["density"]*self.throat["velocity"])
            ae_mdot = 1 / (exit_properties["density"]*velocity)
            ae_at = ae_mdot/at_mdot
            Isp = velocity / sp.constants.g
            Ivac = Isp + pressure * ae_mdot / sp.constants.g
            Cf = velocity / self.chamber["c*"]

            transport_prop, transport_moles = self.__get_transport_properties(gas, transport)

            exit_prop = {   "velocity"                  : velocity, 
                            "mach"                      : velocity/speed_sound, 
                            "area ratio"                : ae_at,
                            "I_sp"                      : Isp,
                            "I_vac"                     : Ivac,
                            "c*"                        : np.nan,
                            "C_f"                       : Cf,
                            "mole fraction"             : gas.mole_fraction_dict(),
                            "mole fraction transport"   : transport_moles}

            exit_properties = exit_properties | transport_prop | exit_prop | exit_derivatives

            self.exit = exit_properties
            
        else: 
            raise ValueError("Invalid input. exit_parameter was not defined as 'pressure' or 'area ratio'")

    def size_engine(self, size_value, size_parameter = 'thrust'): 
        '''
        Description
        --------------------------------------------------
        Given thrust, throat diameter, or massflow will define the other three and update engine_state accordingly

        Parameters
        --------------------------------------------------
        size_value : float
            Depending on what 'size_parameter' is defined in the following ways
                'thrust' : the target thrust of the engine in units newtons
                'mass flow' : the target mass flow of the engine in units kilogram/second
                'throat diameter' : the diameter of the throat in units meter
        size_parameter : string
            a string with three possible values:
                'thrust' 
                'mass flow'
                'throat diameter' 

        Returns
        --------------------------------------------------
        None
        '''

        if size_parameter == 'thrust': 
            thrust = size_value
            mass_flow = thrust / self.exit['velocity']
            throat_area = mass_flow / (self.throat['velocity'] * self.throat['density'])
            # throat_area_2 = eneter ideal compressible choked equation
            # throat_area_3 = thrust / (self.exit['C_f'] * self.chamber['pressure'])
            throat_diameter = np.sqrt(throat_area/np.pi) * 2
            # throat_diameter_2 = np.sqrt(throat_area_2/np.pi) * 2
            # throat_diameter_3 = np.sqrt(throat_area_3/np.pi) * 2


        elif size_parameter == 'mass flow':
            mass_flow = size_value
            thrust = mass_flow * self.exit['velocity']
            throat_area = mass_flow / (self.throat['velocity'] * self.throat['density'])
            throat_diameter = np.sqrt(throat_area/np.pi) * 2
        elif size_parameter == 'throat diameter':
            throat_diameter = size_value
            throat_area = np.pi * (throat_diameter/2)**2
            mass_flow = throat_area * self.throat['velocity'] * self.throat['density']
            thrust = mass_flow * self.exit['velocity']
        else:
            raise ValueError("Invalid input. size_parameter was not defined as 'thrust', 'mass flow' or 'throat diameter'")

        chamber_size = {'thrust' :      np.nan,
                        'mass flow' :   mass_flow,
                        'area' :        np.inf,
                        'diameter' :    np.inf}
        
        throat_size =  {'thrust' :      np.nan,
                        'mass flow' :   mass_flow,
                        'area' :        throat_area,
                        # 'area 2' :      throat_area_2,
                        # 'area 3' :      throat_area_3,
                        'diameter' :    throat_diameter}
                        # 'diameter 2' :  throat_diameter_2,
                        # 'diameter 3' :  throat_diameter_3

        exit_size =    {'thrust' :      thrust,
                        'mass flow' :   mass_flow,
                        'area' :        throat_area * self.exit['area ratio'],
                        'diameter' :    np.sqrt(throat_area * self.exit['area ratio']/np.pi) * 2}
        
        self.chamber.update(chamber_size)
        self.throat.update(throat_size)
        self.exit.update(exit_size)

        self.chamber.update({'heat transfer coefficient': np.nan,'heat flux': np.nan})
        self.throat.update(self._heat_flux(self.throat["temperature"],self.throat["pressure"], self.throat["mole fraction transport"], self.throat['area'], self.throat['mach']))
        self.exit.update(self._heat_flux(self.exit["temperature"],self.exit["pressure"], self.exit["mole fraction transport"], self.exit['area'], self.exit['mach']))

        self.engine_state = pd.DataFrame([self.chamber, self.throat, self.exit], index=["chamber", "throat", "exit"])
            
    def state_at_area(self, gas, condensate, transport, area_ratio, speed = "supersonic"):
        '''
        Description
        --------------------------------------------------
        Iteratively find the properties at a given location in along the nozzle given the area ratio and whether it is in the subsonic or supersonic section of the nozzle.

        Parameters
        --------------------------------------------------
        area_ratio : float
            The ratio at of the area at a given location over the area of the throat.
        speed : string
            a string with two possible values:
                'subsonic' : the converging section of the nozzle
                'supersonic' : the diverging section of the nozzle        

        Returns
        --------------------------------------------------
        local_properties : Dict
            A dictionary of the thermodynamic properties and derivatives 
        
        '''

        # checking for valid input
        if area_ratio <= 1:
            raise ValueError("Area ratio was less than or equal to 1")
        
        # getting initial guess of ln(pc/pe).
        if speed == "subsonic": 
            if area_ratio > 1.000 and area_ratio < 1.09:
                lnpc_p = 0.9*np.log(self.chamber["pressure"]/self.throat["pressure"])/(area_ratio+10.587*np.log(area_ratio)**3+9.454*np.log(area_ratio))
            elif area_ratio >= 1.09:
                lnpc_p = np.log(self.chamber["pressure"]/self.throat["pressure"])/(area_ratio+10.587*np.log(area_ratio)**3+9.454*np.log(area_ratio)) 
                
        if speed == "supersonic": 
            if area_ratio > 1.000 and area_ratio < 2:
                lnpc_p = np.log(self.chamber["pressure"]/self.throat["pressure"]) +np.sqrt(3.294*np.log(area_ratio)**2+1.534*np.log(area_ratio))
            elif area_ratio >= 2:
                lnpc_p = self.chamber["gamma_s"] + 1.4 * np.log(area_ratio)
        
        # initial guess at equilibrium
        pressure = self.chamber["pressure"]/np.exp(lnpc_p)
        gas.SPX = self.throat["entropy"], pressure, self.throat["mole fraction"]
        if self.assumption == 'equilibrium': 
                gas.equilibrate('SP')   
        elif self.assumption == 'frozen':
            pass
        else:
            raise ValueError("Assumption must be 'frozen' or 'equilibrium'")            

        # defining iteration limits and convergence
        num_iter = 0
        max_iter = 10
        tolerance = 4e-5
        residual = 1

        # defines throat area / throat massflow
        at_mdot = 1 / (self.throat["density"]*self.throat["velocity"])

        # iterative solver for state at position _blank_
        while residual > tolerance:
            num_iter += 1

            if num_iter == max_iter:
                print(f"exceeded {max_iter} iterations, residual is {residual} which is not below tolerance of {tolerance}")
                break

            derivatives = self.__get_thermo_derivatives(gas, condensate)
            properties = self.__get_thermo_properties(gas, condensate, derivatives["dpi_dlnT_P"], derivatives["dlnn_dlnT_P"], derivatives["dlnV_dlnT_P"], derivatives["dlnV_dlnP_T"])

            velocity = np.sqrt(2 * (self.chamber["enthalpy"] - properties["enthalpy"]))
            speed_sound = np.sqrt(ct.gas_constant * properties["temperature"] * properties["gamma_s"] / properties["molar mass"])
            a_mdot = 1/(properties["density"]*velocity)
            a_at = a_mdot/at_mdot

            dlnpc_p_dlna_at = properties["gamma_s"] * velocity**2 / (velocity**2 - speed_sound**2)
            lnpc_p = lnpc_p + dlnpc_p_dlna_at * (np.log(area_ratio) - np.log(a_at))
            residual = abs(dlnpc_p_dlna_at * (np.log(area_ratio) - np.log(a_at)))

            pressure = self.chamber["pressure"]/np.exp(lnpc_p)

            gas.SP = self.throat["entropy"], pressure
            if self.assumption == 'equilibrium': 
                gas.equilibrate('SP')   
            elif self.assumption == 'frozen':
                pass
            else:
                raise ValueError("Assumption must be 'frozen' or 'equilibrium'")

        transport_prop, transport_moles = self.__get_transport_properties(gas, transport)

        local_prop  = { "velocity"        : velocity, 
                        "mach"            : velocity/speed_sound, 
                        "area ratio"      : area_ratio,
                        "I_sp"            : np.nan,
                        "I_vac"           : np.nan,
                        "c*"              : np.nan,
                        "C_f"             : np.nan,
                        "mole fraction"   : gas.mole_fraction_dict()}

        local_proprties = properties | transport_prop |local_prop | derivatives

        return local_proprties
    
    def property(self, location, variable):
        if location in self.engine_state.index and variable in self.engine_state.columns:
            return self.engine_state[variable][location]

class Engine(EngineState):     
    def __init__(self,  oxidizer, fuel, of_ratio, pressure, exit_value, gas,  size_value, condensate = None, transport=None, assumption='equilibrium', size_parameter ="thrust",
                 temp_oxidizer=None, temp_fuel=None, combustion_products=None, exit_parameter="pressure", throat_inlet_radius_ratio=1.5, temperature_wall=None):
        '''
        This class will define the properties of an rocket engine at the chamber, throat and exit 
        given the propellant, and conditions that engine is in.

        Parameters
        --------------------------------------------------
        oxidizer : string
            A string of a species name which is used to define a cantera solution representing the oxidizer.
            Define based off the list of phase names from the propellants.yaml file TODO: kill propellants.yaml and create a list users can choose from 
        fuel : string
            Similar to oxidizer, this is a string of a species name which is used to define a cantera solution representing the fuel.
            Define based off the list of phase names from the propellants.yaml file TODO: kill propellants.yaml and create a list users can choose from 
        of_ratio : float
            the mass ratio of oxidizer to fuel defined as mass_oxidizer / mass_fuel.
        pressure : float
            The pressure of the combustion chamber in units of Pascals (Pa).
        exit_value : float
            The exit conditions of the engine. Depending on what 'exit parameter' is
            defined to be the exit condition can be defined in the following ways:
                'pressure' (default): Exit pressure of the gases in Pa should be greater than 0 but lower than chamber pressure
                'area ratio': Definfed as area_exit / area_throat should be value greater than 1 TODO: add functionality
        size_value : float 
            Determines the size of the engine. depending on what 'size parameter' is defined to be. 
            exit condition can be defined in the following way
        exit_parameter : string
            A string which determines how you are defining the exit condition of your engine. Can either be set to 'pressure' (default)
            or 'area ratio'. 
        temp_oxidizer : float
            A positive float value defining the temperature of the oxidizer in Kelvin (K). If the oxidizer is liquid the temperature 
            is assumed to be saturated temperatures at standard pressure, and gaseous oxidizer will be set to a room tempearture of 295.15K. TODO: add functionality
        temp_fuel : float
            Similar to above. positive float value defining the temperature of the fuel in Kelvin (K). If the fuel is liquid the temperature 
            is assumed to be saturated temperatures at standard pressure, and gaseous oxidizer will be set to a room tempearture of 295.15K. TODO: add functionality
        

        Attributes:
        --------------------------------------------------
        Properties : Pandas Dataframe
        
        Methods
        --------------------------------------------------

        A class for the state of a rocket engine. Given the propellants and conditions of a rocket engine this 
        class will define the properties of an engine at
        '''
        super().__init__(oxidizer, fuel, of_ratio, pressure, exit_value, gas, size_value = size_value, condensate=condensate, transport=transport, size_parameter= size_parameter,
                        exit_parameter=exit_parameter, throat_inlet_radius_ratio=throat_inlet_radius_ratio, temperature_wall=temperature_wall, assumption=assumption)
        
    def __chamber_contour(self, length_value, contraction_ratio, contraction_angle=30, nozzle_inlet_radius_ratio=0.5, 
                         throat_inlet_radius_ratio = None, length_parameter = 'characteristic length'):
        '''
        Description
        --------------------------------------------------
        will create

        Parameters
        --------------------------------------------------
        leng_value : float
            Depending on what 'size_parameter' is defined in the following ways
                'thrust' : the target thrust of the engine in units newtons
                'mass flow' : the target mass flow of the engine in units kilogram/second
                'throat diameter' : the diameter of the throat in units meter
        length_parameter : string
            a string with three possible values:
                'characteristic length' 
                'chamber length'
                'throat diameter' 

        Returns
        --------------------------------------------------
        None
        '''

        theta_c = np.radians(contraction_angle)

        r_t = self.throat["diameter"]/2
        r_c = np.sqrt(contraction_ratio * r_t**2)
        r_tin = self.throat_inlet_radius_ratio * r_t
        r_ninmax = (r_c-r_t)/(1-np.cos(theta_c)) - r_tin
        r_nin = r_ninmax * nozzle_inlet_radius_ratio

        self._r_nin = r_nin
        self._r_tin = r_tin

        if self.throat_inlet_radius_ratio != throat_inlet_radius_ratio and throat_inlet_radius_ratio:
            self.throat_inlet_radius_ratio = throat_inlet_radius_ratio
            self.throat.update(self._heat_flux(self.throat["temperature"],self.throat["pressure"], self.throat_transport.X, self.throat['area'], self.throat['mach']))
            self.exit.update(self._heat_flux(self.exit["temperature"],self.exit["pressure"], self.exit_transport.X, self.exit['area'], self.exit['mach']))
        
        # checks for valid inputs: 
        if not nozzle_inlet_radius_ratio <=1:
            raise ValueError("the nozzle_inlet_radius_ratio must not exceed 1. This is the ratio of the radius/ maximum possible radius")
        if r_tin *(1- np.cos(theta_c)) > r_c-r_t:
            raise ValueError("throat_inlet_radius_ratio is too high. The nozzle geometry is not possible with the given contraction ratio and contraction angle")
        

        # defines bounds for x.        
        x_tin = -r_tin * np.sin(theta_c)
        x_nin = -1/np.tan(theta_c)*(r_c-r_nin*(1-np.cos(theta_c))-
                                                    (r_t+r_tin*(1-np.cos(theta_c))))-r_tin*np.sin(theta_c)
        x_c = x_nin - r_nin * np.sin(theta_c)

        # find volume of converging section
        volume_throat_inlet = np.pi*(-(x_nin-x_c)**3/3 + (r_c-r_nin)*r_nin**2*(np.asin((x_nin-x_c)/r_nin)+1/2*np.sin(2*np.asin((x_nin-x_c)/r_nin)))+(r_nin**2+(r_c-r_nin)**2)*(x_nin-x_c))
        sub = r_t+r_tin*(1-np.cos(theta_c)-np.tan(theta_c)*np.sin(theta_c))
        volume_nozzle_line = np.pi*(np.tan(theta_c)**2*x_tin**3/3-np.tan(theta_c)*sub*x_tin**2+sub**2*x_tin) - np.pi*(np.tan(theta_c)**2*x_nin**3/3-np.tan(theta_c)*sub*x_nin**2+sub**2*x_nin) 
        volume_nozzle_inlet =  np.pi*(-x_nin**3/3+x_c*x_nin**2+(r_c-r_nin)*r_nin**2*(np.asin((x_nin-x_c)/r_nin)+np.sin(2*np.asin((x_nin-x_c)/r_nin))/2)+((r_c-r_nin)**2+r_nin**2-x_c**2)*x_nin) -np.pi*(2*x_c**3/3+((r_c-r_nin)**2+r_nin**2-x_c**2)*x_c) 
        volume_converging = volume_throat_inlet + volume_nozzle_line + volume_nozzle_inlet
        

        # find length of chamber
        if length_parameter == 'characteristic length':
            l_star = length_value            
            volume_total = self.throat["area"] * length_value

            if volume_converging > volume_total:
                raise ValueError("the converging section volume is larger than the total volume of your chamber. try increasing your L* or reducing your contraction ratio.")
            volume_chamber = volume_total - volume_converging

            l_c = volume_chamber/(np.pi*r_c**2)
        elif length_parameter == 'chamber length':
            l_c = length_value
            volume_chamber = np.pi * r_c **2 * l_c
            volume_total = volume_chamber + volume_converging
            l_star = volume_total / self.throat["area"]

        x_inj = x_c - l_c

        # print(f"r_t: {r_t}, r_c: {r_c}, r_tin: {r_tin}, r_nin = {r_nin}")
        # print(f"throat_inlet: {x_tin}, nozzle_inlet: {x_nin}, chamber_end: {x_c}, x_inj: {x_inj}")
        # print(f"volume_throat_inlet: {volume_throat_inlet}, volume_nozzle_line: {volume_nozzle_line}, volume_nozzle_inlet: {volume_nozzle_inlet}, volume_converging: {volume_converging}, volume_chamber: {volume_chamber}, volume_total: {volume_total}")

        x_dict = {"x_inj" : 0,  "x_c" : x_c-x_inj, "x_nin" : x_nin-x_inj, "x_tin" : x_tin-x_inj, "x_t" : -x_inj}

        def chamber_contour(x): 
            x = x+x_inj
            # Defines throat
            if x == 0:
                return r_t

            # defines radius before throat
            elif x_tin <= x and x < 0:
                return -np.sqrt(r_tin**2 - x**2) + (r_tin + r_t)
            
            # defines contraction line
            elif x_nin <= x and x < x_tin:
                return -np.tan(theta_c) * (x + r_tin * np.sin(theta_c))+r_t+r_tin*(1-np.cos(theta_c))
            
            # defines radius befor nozzle
            elif x_c <= x and x < x_nin: 
                return np.sqrt(r_nin**2 - (x-x_c)**2) + (r_c - r_nin)

            # defines chamber section
            elif x_inj <= x and x < x_c:
                return r_c
            
            else: 
                raise ValueError(f"x coordinate exceeds bounds of chamber. please make sure 0 <= x <= {-x_inj}")
            
        return chamber_contour, x_dict

    def conical_contour(self, length_value, contraction_ratio, contraction_angle=30, nozzle_inlet_radius_ratio=0.5, throat_inlet_radius_ratio = None, 
                                length_parameter = 'characteristic length', throat_outlet_radius_ratio=0.382, expansion_angle=15, fidelity = 500):
        
        chamber_contour, chamber_x_dict = self.__chamber_contour(length_value, contraction_ratio, contraction_angle, nozzle_inlet_radius_ratio, throat_inlet_radius_ratio, length_parameter = length_parameter)

        theta_e = np.radians(expansion_angle)
        r_t = self.throat["diameter"]/2
        r_tout = throat_outlet_radius_ratio * r_t
        self._r_tout = r_tout

        x_tout = r_tout * np.sin(theta_e) + chamber_x_dict["x_t"]
        x_e = (self.exit["diameter"]/2-r_t-r_tout*(1-np.cos(theta_e)))/np.tan(theta_e)+r_tout*np.sin(theta_e) + chamber_x_dict["x_t"]
        
        nozzle_x_dict = {"x_tout" : x_tout, "x_e" : x_e}

        x_dict = chamber_x_dict | nozzle_x_dict
        
        def contour(x):
            if 0 <= x and x <= x_dict["x_t"]: 
                return chamber_contour(x)
            elif x_dict["x_t"] < x <= x_dict["x_tout"]:
                x = x-x_dict["x_t"]
                return -np.sqrt(-x**2+r_tout**2) + r_t +r_tout
            elif x_dict["x_tout"] < x and x <= x_dict["x_e"]:
                x = x-x_dict["x_t"]
                return np.tan(theta_e)*(x-r_tout*np.sin(theta_e))+r_t+r_tout*(1-np.cos(theta_e))
            else:
                raise ValueError(f"x coordinate exceeds bounds of nozzle. please make sure {x_dict['x_inj']} < x <= {x_dict['x_e']}")

        x_coords = np.linspace(x_dict["x_inj"], x_dict["x_e"], fidelity)
        x_coords = np.sort(np.append(np.linspace(x_dict["x_inj"], x_dict["x_e"], fidelity), list(x_dict.values())))
        contour_coords = []
        # for x in x_coords: 
        #     contour_coords.append([x, contour(x)])

        self.contour = contour
        self.x_dict = x_dict
        self.contour_coords = contour_coords

        # return contour_coords

    def properties_along_contour(self, fidelity):
        stations = np.hstack((0, np.linspace(self.x_dict["x_c"], self.x_dict["x_e"], fidelity, endpoint = True)))
        
        if self.x_dict["x_t"] not in  stations: 
            idx = np.searchsorted(stations, self.x_dict["x_t"])
            stations = np.insert(stations, idx, self.x_dict["x_t"])
        
        self._gas.TPX = self.chamber['temperature'], self.chamber['pressure'], self.chamber['mole fraction'] 
        list_station_properties = []

        for station in stations:
            area_ratio = np.pi*self.contour(station)**2 / self.throat["area"]
            if 0 <= station and station < self.x_dict['x_t']:   
                station_properties = self.state_at_area(self._gas, self._condensate, self._transport, area_ratio, speed = "subsonic")
                station_size = {'thrust' :      np.nan,
                                'mass flow' :   self.throat['mass flow'],
                                'area' :        self.throat['area'] * area_ratio,
                                'diameter' :    self.contour(station)*2}
                
                station_heat = self._heat_flux(station_properties['temperature'], station_properties['pressure'], self._transport.X, 
                                                          station_size['area'], station_properties['mach'], viscosity_exponent = 0.6)

                station_properties.update(station_size | station_heat)
            
            elif station == self.x_dict['x_t']:
                station_properties = self.throat

            elif self.x_dict['x_t'] < station and station < self.x_dict['x_e']:
                station_properties = self.state_at_area(self._gas, self._condensate, self._transport, area_ratio, speed = "supersonic")
                station_size = {'thrust' :      np.nan,
                                'mass flow' :   self.throat['mass flow'],
                                'area' :        self.throat['area'] * area_ratio,
                                'diameter' :    self.contour(station)*2}
                
                station_heat = self._heat_flux(station_properties['temperature'], station_properties['pressure'], self._transport.X, 
                                                          station_size['area'], station_properties['mach'])

                station_properties.update(station_size | station_heat)
            
            elif station == self.x_dict['x_e']:
                station_properties = self.exit
        
            else:
                raise ValueError(f"Invalid station: station must be greater than 0 (injector face) and less that {self.x_dict['x_e']} (nozzle exit) in order to be valid.")
            
            # channel_geom = self.channel_contour(station)

            list_station_properties.append(station_properties)  # | channel_geom)

        self.property_contour = pd.DataFrame(list_station_properties, stations)

                

        return self.property_contour 

    def parabolic_contour():
        return
    
    def bell_contour():
        return
    
    def channel_design(self, n_channels, stations, wall_thickness,  channel_height, channel_width):
        if len(stations) != len(wall_thickness) or len(stations) != len(channel_height) or len(stations) != len(channel_width):
           raise ValueError("length of inputs inconsistent, insure that the number of stations, thicknesses, heights, and widths defined are equal.")
       
        if len(stations) < 1:
           raise ValueError("1 or more stations need to be defined")

        for i, station in enumerate(stations):
            if  n_channels * 2 * np.atan(channel_width[i]/2/(wall_thickness[i]+self.contour(station))) >= 2*np.pi:
                raise ValueError(  f"Specified channel geometry results in no wall between channels at station: {station}. please reduce channel width or number of channels")
            
            if not 0 <= station and stations <= self.x_dict["x_e"]:
                raise ValueError(f"Invalid station value {station} for the channel contour, all values must be between 0 and {self.x_dict['x_e']}")
        
        def channel_countour(x):
            r = self.contour(x)
            if 0 <= x and x <= stations[0]:
                t = wall_thickness[0]
                h = channel_height[0]
                w = channel_width[0]
                return {'chamber wall'    : t,
                        'channel height'    : h, 
                        'channel width'     : w,
                        'channel wall'      : 2*np.sqrt(w**2+(t+r)**2) * np.sin((2* np.pi - n_channels * np.atan(w/2/(t+r)))/n_channels)}
            
            for i in range(1, len(stations)):
                if stations[i-1] < x and x <= stations[i]:
                    t = wall_thickness[i-1] + (x-stations[i-1]) * (wall_thickness[i]-wall_thickness[i-1]) / (stations[i]-stations[i-1])
                    h = channel_height[i-1] + (x-stations[i-1]) * (channel_height[i]-channel_height[i-1]) / (stations[i]-stations[i-1])
                    w = channel_width[i-1] + (x-stations[i-1]) * (channel_width[i]-channel_width[i-1]) / (stations[i]-stations[i-1])
                    return {'chamber wall'    : t,
                            'channel height'    : h, 
                            'channel width'     : w, 
                            'channel wall'      : 2*np.sqrt(w**2+(t+r)**2) * np.sin((2* np.pi - n_channels * np.atan(w/2/(t+r)))/n_channels)}

            if stations[-1] <= x and x <= self.x_dict["x_e"]:
                t = wall_thickness[-1]
                h = channel_height[-1]
                w = channel_width[-1]
                return {'chamber wall'      : t,
                        'channel height'    : h, 
                        'channel width'     : w,
                        'channel wall'      : 2*np.sqrt(w**2+(t+r)**2) * np.sin((2* np.pi - n_channels * np.atan(w/2/(t+r)))/n_channels)}

            if not 0 <= x and x <= self.x_dict["x_e"]:
                raise ValueError(f"Invalid x value {x} for the channel contour, all values must be between 0 and {self.x_dict['x_e']}")    

        self.n_channels = n_channels
        # self.channel_contour = channel_countour

def chemistry_initializer(oxidizer, fuel, temp_oxidizer=None, temp_fuel=None, combustion_products=None):
    propellant_array = [['oxidizer', oxidizer, temp_oxidizer],['fuel', fuel, temp_fuel]]
    reactants_all   = {S.name: S for S in ct.Species.list_from_file('chem_prop/reactants.yaml')}
    gaseous_all     = {S.name: S for S in ct.Species.list_from_file('chem_prop/gaseous_products.yaml')}
    condensed_all   = {S.name: S for S in ct.Species.list_from_file('chem_prop/condensed_products.yaml')}
    propellant = {}
    propellant_elements = []

    def transport_species(gaseous_species):
        transport_species = [sp for sp in gaseous_species.values() if sp.transport is not None]
        if not transport_species:
            return None
        transport = ct.Solution(thermo="ideal-gas", transport = "mixture-averaged", species = transport_species)
        transport.transport_model = "mixture-averaged"
        return transport
    
    # initialize propellants
    for type, name, temp in propellant_array:
        if name in reactants_all:
            if temp is None:
                if reactants_all[name].thermo.input_data['model'] == 'constant-cp':
                    temp = reactants_all[name].thermo.input_data['T0']
                else:
                    temp= 298.15
            propellant[type] = ct.Solution(thermo='ideal-gas', species=[reactants_all[name]])
            propellant[type].TP = temp, 2e5
            propellant[type].equilibrate('TP')
        
        elif name in gaseous_all:
            if temp is None:
                temp = 298.15
            propellant[type] = ct.Solution(thermo='ideal-gas', species=[gaseous_all[name]])
            propellant[type].TP = temp, 2e5
            propellant[type].equilibrate('TP')

        
        elif name in condensed_all:
            if temp is None:
                temp = 298.15
            propellant[type] = ct.Solution(thermo='ideal-gas', species=[condensed_all[name]])
            propellant[type].TP = temp, 2e5
            propellant[type].equilibrate('TP')

        else:
            raise ValueError(f"{name} does not exist in the thermal database.")    
        
        propellant_elements.extend(propellant[type].element_names)

    # ensures propellant elements are unique
    propellant_elements = list(set(propellant_elements))

    if combustion_products == None: 
        condensed_species = {k: v for k, v in condensed_all.items() if set(v.composition.keys()).issubset(propellant_elements)}
        gaseous_species = {k: v for k, v in gaseous_all.items() if set(v.composition.keys()).issubset(propellant_elements) and ('+' not in k and '-' not in k)}
        transport = transport_species(gaseous_species)
    else: 
        # checks if designated combustion products exist in the thermal database. if it does then filters the data
        if set(combustion_products).issubset(set(condensed_all.keys()) | set(gaseous_all.keys())):
            condensed_species = {k: v for k, v in condensed_all.items() if k in combustion_products and set(v.composition.keys()).issubset(propellant_elements)}
            gaseous_species = {k: v for k, v in gaseous_all.items() if k in combustion_products and set(v.composition.keys()).issubset(propellant_elements) and ('+' not in k and '-' not in k)}
            transport = transport_species(gaseous_species)
        else:
            invalid_species = set(combustion_products) - (set(condensed_all.keys()) | set(gaseous_all.keys()))
            raise ValueError(f"the following species {invalid_species} do not exist in the thermal database.")   
    
    if condensed_species:
        condensate = ct.Solution(thermo='ideal-gas', species=list(condensed_species.values()))
    else: 
        condensate = None
    if gaseous_species:
        gas =  ct.Solution(thermo='ideal-gas', species=list(gaseous_species.values()))
    else: 
        raise ValueError("No valid gas phase combustion products of given propellants. Consider changing propllants or inputted combustion products")        

    return propellant['oxidizer'], propellant['fuel'], gas, condensate, transport

@st.cache_data
def ranged_sim_rocketcow(ox, f, t_ox, t_fuel, of_arr, p_arr, exit_value, size_value=None, size_parameter ="thrust", exit_parameter="pressure", assumption = 'Equilibrium'):
    oxidizer, fuel, gas, condensate, transport = chemistry_initializer(ox, f, temp_oxidizer=t_ox, temp_fuel=t_fuel)
    state_list = []
    for pressure in p_arr:
        for of_ratio in of_arr:
            state = Engine(oxidizer, fuel, of_ratio, pressure, exit_value, gas, size_value, condensate=condensate, transport=transport, size_parameter=size_parameter,
                           exit_parameter=exit_parameter, assumption=assumption, temperature_wall=500)
            state_list.append(state())
    states = pd.concat(state_list, keys = list(range(len(state_list))))
    states.drop(columns=['mole fraction', 'mole fraction transport', 'dpi_dlnT_P', 'dpi_dlnP_T'])
    return states

@st.cache_data
def data_filter(df, pos, var):
    core = df.loc[pd.IndexSlice[:, 'chamber'], ['of ratio', 'pressure']]
    core.index = core.index.droplevel(1)
    param = df.loc[pd.IndexSlice[:, pos], var]
    param.index = param.index.droplevel(1)
    core[var] = param
    core['pressure'] = core['pressure'].round(3)
    output = pd.pivot(core, index= 'pressure',  columns = 'of ratio')
    output = output.iloc[::-1]
    output = output.droplevel(0, 1)
    output.index.name = None
    output.columns.name = None
    output.columns = ['{:.1f}'.format(x) for x in output.columns.round(2)]
    output.attrs = core.attrs
    return output

@st.cache_data
def gen_plot(df, pos, var, _plot_type = 'heatmap'): 
    data_filtered = data_filter(df, pos, var)
    unit_dict = {'temperature': '(K)', 'density': 'kg/m^3', 'specific volume': 'm^3/kg', 'enthalpy': 'J/kg', 'internal energy': 'J/kg', 
                 'gibbs': 'J/kg', 'entropy' : 'J/kg/K', 'molar mass' : 'kg/kmol', 'c_p' : 'J/kg/K','c_v' : 'J/kg/K', 'gamma' : '', 
                 'gamma_s' : '', 'speed sound' : '(m/s)', 'viscosity' : '(Pa*s)','thermal conductivity': 'W/m/K', 'prandtl number' : '', 
                 'velocity' : '(m/s)', 'mach' : '','area ratio' : '', 'I_sp' : '(s)', 'I_vac' : '(s)', 'c*' : '(m/s)', 'C_f' : '', 
                 'thrust' : '(N)', 'mass flow' : 'kg/s','area' : 'm^2', 'diameter' : 'm', 'heat transfer coefficient' : 'W/m^2/K', 
                 'heat flux' : 'W/m^2'}
    if _plot_type == 'heatmap': 
        fig = px.imshow(data_filtered, width=600, height=600, origin='lower', color_continuous_scale='viridis', 
                        labels={"x": "OF Ratio (% weight)", "y": "Pressure (Pa)", "hover": f"{var} {unit_dict[var]}",'color': f"{var} {unit_dict[var]}"}, 
                        title="{} at the {}".format(var, pos), aspect="auto")
    elif _plot_type == 'surface':
        surface = go.Surface(z=data_filtered.values, x = data_filtered.columns.values , y = data_filtered.index.values, colorscale = 'Viridis')
        fig = go.Figure(data = [surface])
        fig.update_layout(xaxis=dict(title = dict(text="OF Ratio (% weight)")),
                          yaxis=dict(title = dict(text="Pressure (Pa)")),
                          zaxis=dict(title = dict(text=f"{var} {unit_dict[var]}")),)
    else: 
        raise Exception('invalid plot type, options are \'heatmap\' or \'surface\'')
    
    return fig, data_filtered

@st.cache_data
def convert_df(df):
    return df.to_csv().encode('utf-8')

@st.cache_data
def contour_design( ox, f, t_ox, t_fuel, of, pressure, pressure_exit, thrust, length_value, contraction_ratio, contraction_angle=30, nozzle_inlet_radius_ratio=0.5, throat_inlet_radius_ratio = None, 
                    length_parameter = 'characteristic length', throat_outlet_radius_ratio=0.382, expansion_angle=15, temperature_wall=500, fidelity = 200):
    oxidizer, fuel, gas, condensate, transport = chemistry_initializer(ox, f, temp_oxidizer=t_ox, temp_fuel=t_fuel)
    da_engine = Engine(oxidizer, fuel, of, pressure, pressure_exit, gas, thrust, condensate=condensate, transport=transport)
    da_engine.conical_contour(length_value, contraction_ratio, contraction_angle=contraction_angle, nozzle_inlet_radius_ratio=nozzle_inlet_radius_ratio, throat_inlet_radius_ratio = throat_inlet_radius_ratio, 
                                length_parameter = length_parameter, throat_outlet_radius_ratio=throat_outlet_radius_ratio, expansion_angle=expansion_angle, fidelity = fidelity)
    properties = da_engine.properties_along_contour(fidelity)
    properties.drop(columns=['mole fraction', 'mole fraction transport', 'dpi_dlnT_P', 'dpi_dlnP_T'])
    return properties

@st.cache_data
def contour_plot(contour, var1, var2, var3, var4):
    fig = go.Figure()
    fig.update_layout(xaxis=dict(title = dict(text="Axial distance (m)"), domain=[0, 0.9]), 
                    title_text="Properties along Engine contour")
    unit_dict = {'temperature': '(K)', 'density': 'kg/m^3', 'specific volume': 'm^3/kg', 'enthalpy': 'J/kg', 'internal energy': 'J/kg', 
                 'gibbs': 'J/kg', 'entropy' : 'J/kg/K', 'molar mass' : 'kg/kmol', 'c_p' : 'J/kg/K','c_v' : 'J/kg/K', 'gamma' : '', 
                 'gamma_s' : '', 'speed sound' : '(m/s)', 'viscosity' : '(Pa*s)','thermal conductivity': 'W/m/K', 'prandtl number' : '', 
                 'velocity' : '(m/s)', 'mach' : '','area ratio' : '', 'I_sp' : '(s)', 'I_vac' : '(s)', 'c*' : '(m/s)', 'C_f' : '', 
                 'thrust' : '(N)', 'mass flow' : 'kg/s','area' : 'm^2', 'diameter' : 'm', 'heat transfer coefficient' : 'W/m^2/K', 
                 'heat flux' : 'W/m^2', 'pressure' : 'Pa'}
    fig.add_trace(go.Scatter(x = contour.index, y = contour['diameter']/2, name= f"radius (m)"))

    for i, prop in enumerate([var1, var2, var3, var4]): 
        fig.add_trace(go.Scatter(x = contour.index, y = contour[prop], name= f"{var1} {unit_dict[var1]}", yaxis=f"y{i+2}"))

    fig.update_layout(yaxis=dict(title=dict(text="Radius (m)"),side="right", rangemode='tozero',position=0.9, automargin = True),
                    yaxis2=dict(title=dict(text=f"{var1} {unit_dict[var1]}"),overlaying="y"),
                    yaxis3=dict(title=dict(text=f"{var2} {unit_dict[var2]}"),anchor="free",overlaying="y",autoshift=True),
                    yaxis4=dict(title=dict(text=f"{var3} {unit_dict[var3]}"),anchor="free",overlaying="y",autoshift=True),
                    yaxis5=dict(title=dict(text=f"{var4} {unit_dict[var4]}"),anchor="free",overlaying="y",autoshift=True))
    return fig

### Creates the UI ###
st.header('Rocket Cow Engine Visualizer')
st.write('''Welcome to V0.1 of the :rainbow[Rocket Cow\u2122]. V0.1 is capable of plotting engine properties over a 
         range of pressures and OF ratios, calculating gas side heat transfer, designing an engine contour (with a 
         conical nozzle) and plotting properties along that contour. Rocket Cow uses Cantera along with the methodology 
         laid out in NASA CEA to find engine properties. Please note some functionality is hidden and only available 
         via script. this UI is just a demo, if there are any bugs or requests please create an issue on github.''')

with st.form('OF ratio & Pressure Trade'):
    # Prop input
    st.subheader("Propellants")
    col1, col2 = st.columns(2)
    col1.selectbox('Fuel:', ['CH4(L)', 'RP-1(L)', 'H2(L)', 'C2H5OH(L)', 'C3H8(L)'], key='fuel')
    col1.number_input('Fuel Temperature (K):', key= 'temp_fuel', min_value= 0.0, value=None, step=5.0, placeholder="Optional")
    col2.selectbox('Oxidizer:',['O2(L)', 'N2O(L)', 'Air(g)'], key='ox')
    col2.number_input('Oxidizer Temperature (K):', key= 'temp_ox', min_value= 0.0, value=None, step=5.0, placeholder="Optional")
    
    # Pressure input
    st.divider()
    st.subheader("Pressure")

    col1, col2, col3 = st.columns([0.38,0.38,0.24   ])
    col1.number_input('Min Chamber Pressure (bar):', key='p_min', min_value=0.0, value=10.0, step=5.0)
    col2.number_input('Max Pressure (bar):', key='p_max', min_value = 0.0, value = 100.0,  step = 5.0)
    col3.number_input('Step Size (bar):', key='p_step', min_value=0.01, value = 10.0, step = 1.0)  

    st.number_input('Ambient Pressure (bar):', key='p_e', min_value=0.0 ,value=1.0 ,step=0.1)

    # OF conditions
    st.divider()
    st.subheader("OF Ratio")
    
    col1, col2, col3 = st.columns(3)
    col1.number_input('Min OF (%wt ratio):', key='of_min', min_value = 0.01, value = 1.5, step=0.1)
    col2.number_input('Max OF (%wt ratio):', key='of_max', min_value = 0.01, value = 4.5,  step = 0.1)
    col3.number_input('Step Size (%wt ratio):', key='of_step', min_value = 0.0001, value= 0.2, step= 0.05)
    
    # Thrust
    st.divider()
    st.subheader("Thrust")
    st.number_input('Thrust (N):', key='thrust', min_value=0.0 ,value=5000.0 ,step=0.1)
    
    # sim type and run
    st.divider()  
    st.subheader("Assumptions")
    st.selectbox('reacting flow condition', ['equilibrium', 'frozen'    ], key='assume')
    st.form_submit_button('Run', use_container_width = True, on_click=run_button)

if state.run_button:
    
    if(state.p_min > state.p_max):
        st.warning('Your minimum pressure is greater than your maximum pressure! Adjust and rerun.', icon="⚠️")
        st.stop()
    if(state.p_min < state.p_e):
        st.warning('Your your exit pressure exceeds your minimum pressure! Your chamber pressure must always exceed ambient pressure. Adjust and rerun.', icon="⚠️")
        st.stop()
    if(state.p_step > (state.p_max-state.p_min)):
        st.warning('Your pressure step size is greater than your pressure range! Adjust and rerun', icon="⚠️")
        st.stop()

    if(state.of_min > state.of_max):
        st.warning('Your minimum OF ratio is greater than your maximum OF ratio! Adjust and rerun', icon="⚠️")
        st.stop()
    if(state.of_step > (state.p_max-state.p_min)):
        st.warning('Your OF ratio step size is greater than your OF ratio range! Adjust and rerun', icon="⚠️")
        st.stop()
    if(state.of_step > (state.p_max-state.p_min)):
        st.warning('Your OF ratio step size is greater than your OF ratio range! Adjust and rerun', icon="⚠️")
        st.stop()

    state.data = ranged_sim_rocketcow(state.ox, state.fuel, state.temp_ox, state.temp_fuel, list(np.arange(state.of_min, state.of_max, state.of_step)), 
                                    list(np.arange(state.p_min*1e5,state.p_max*1e5, state.p_step*1e5)), state.p_e*1e5, size_value=state.thrust, assumption=state.assume)
    
    st.header('OF ratio and Pressure trade Results')
    if st.checkbox('Show raw data'):
        st.subheader('Raw data')
        st.write(state.data)
        csv = convert_df(state.data)
        st.download_button("Press to Download", csv, "file.csv", "text/csv", key='download-csv')

    st.subheader('Visualization')
    col1, col2, col3 = st.columns(3)
    with col1:
        st.selectbox('Position:', ['chamber', 'throat', 'exit'], key='plot_pos')

    with col2:
        st.selectbox('Variable:', [ 'temperature', 'density', 'specific volume',
                                    'enthalpy', 'internal energy', 'gibbs', 'entropy', 'molar mass', 'c_p',
                                    'c_v', 'gamma', 'gamma_s', 'speed sound', 'viscosity',
                                    'thermal conductivity', 'prandtl number', 'velocity', 'mach',
                                    'area ratio', 'I_sp', 'I_vac', 'c*', 'C_f', 'thrust', 'mass flow',
                                    'area', 'diameter', 'heat transfer coefficient', 'heat flux'], key='plot_var')

    with col3:
        st.selectbox('Plot Type:', ['heatmap', 'surface'], key='plot_type')
    
    st.button('Plot', use_container_width=True, on_click=plot_button)

    if state.plot_button:
        state.figure, state.data_filtered = gen_plot(state.data, state.plot_pos, state.plot_var, _plot_type=state.plot_type)
        st.plotly_chart(state.figure)

        if st.checkbox('Show plotted data'):
            st.subheader('Plotted data')
            st.dataframe(state.data_filtered)
            csv = convert_df(state.data_filtered)
            st.download_button("Press to Download", csv, "file.csv", "text/csv")  

        with st.form('Engine Design'):            
            # Pressure input
            st.subheader("Engine Conditions")

            col1, col2= st.columns([0.5,0.5])
            col1.number_input('Chamber Pressure (bar):', key='pressure', min_value=0.0, value=100.0, step=5.0)
            col2.number_input('OF ratio (% weight)', key='of_ratio', min_value = 0.0, value = 3.0,  step = .1)
            st.number_input('Thrust (N):', key='engine_thrust', min_value=0.0 ,value=5000.0 ,step=0.1)
            st.number_input('Ambient Pressure (bar):', key='pressure_exit', min_value=0.0 ,value=1.0 ,step=0.1)

            # OF conditions
            st.divider()
            st.subheader("Engine Geometry")
            col1, col2 = st.columns(2)
            col1.selectbox('Length parameter', [ 'chamber length', 'characteristic length'], key='length_param')
            col2.number_input('length (m):', key='length', min_value=0.0 ,value=0.3 , step=0.1)
            col1, col2, col3 = st.columns(3)
            col1.number_input('contraction ratio:', key='contraction_ratio', min_value = 1.0001, value = 2.0, step=0.1)
            col2.number_input('contraction angle (deg):', key='contraction_angle', min_value = 1.0, value = 30.0, max_value=89.0,  step = 1.0)
            col3.number_input('expansion angle (deg):', key='expansion_angle', min_value = 1.0, value= 15.0, max_value=89.0, step= 1.0)
            col1, col2, col3 = st.columns(3)
            col1.number_input('nozzle inlet radius ratio:', key='nin_ratio', min_value = 0.01, value = 0.5, max_value=1.0, step=0.1)
            col2.number_input('throat inlet radius ratio:', key='tin_ratio', min_value = 0.01, value = 1.5,  step = 0.1)
            col3.number_input('throat outlet radius ratio:', key='tout_ratio', min_value = 0.01, value= 0.382, step= 0.1)


            # Fidelity
            st.divider()
            st.subheader("Fidelity")
            st.number_input('Number of Stations:', key='fidelity', min_value = 5, value=200, step= 1)
            st.form_submit_button('Run', use_container_width = True, on_click=contour_button)

        if state.contour_button:

            state.contour_properties = contour_design(  state.ox, state.fuel, state.temp_ox, state.temp_fuel, state.of_ratio, state.pressure, 
                                                        state.pressure_exit, state.engine_thrust, state.length, state.contraction_ratio, 
                                                        contraction_angle=state.contraction_angle, nozzle_inlet_radius_ratio=state.nin_ratio, 
                                                        throat_inlet_radius_ratio = state.tin_ratio, length_parameter = state.length_param, 
                                                        throat_outlet_radius_ratio=state.tout_ratio, expansion_angle=state.expansion_angle, 
                                                        temperature_wall=500, fidelity = state.fidelity)

            st.header('Properties Along contour')
            if st.checkbox('Show raw data', key="contour_check"):
                st.subheader('Raw data')
                st.write(state.contour_properties)
                csv = convert_df(state.contour_properties)
                st.download_button("Press to Download", csv, "file.csv", "text/csv", key='download-csv-contour')
            
            col1, col2 = st.columns(2)
            col1.selectbox('Variable:', [ 'temperature', 'pressure', 'density', 'specific volume',
                                    'enthalpy', 'internal energy', 'gibbs', 'entropy', 'molar mass', 'c_p',
                                    'c_v', 'gamma', 'gamma_s', 'speed sound', 'viscosity',
                                    'thermal conductivity', 'prandtl number', 'velocity', 'mach',
                                    'area ratio', 'I_sp', 'I_vac', 'c*', 'C_f', 'thrust', 'mass flow',
                                    'area', 'diameter', 'heat transfer coefficient', 'heat flux'], key='first')
            
            col2.selectbox('Variable:', [ 'velocity','temperature', 'pressure', 'density', 'specific volume',
                                    'enthalpy', 'internal energy', 'gibbs', 'entropy', 'molar mass', 'c_p',
                                    'c_v', 'gamma', 'gamma_s', 'speed sound', 'viscosity',
                                    'thermal conductivity', 'prandtl number' , 'mach',
                                    'area ratio', 'I_sp', 'I_vac', 'c*', 'C_f', 'thrust', 'mass flow',
                                    'area', 'diameter', 'heat transfer coefficient', 'heat flux'], key='second')
            col1, col2 = st.columns(2)
            col1.selectbox('Variable:', [ 'pressure','temperature', 'density', 'specific volume',
                                    'enthalpy', 'internal energy', 'gibbs', 'entropy', 'molar mass', 'c_p',
                                    'c_v', 'gamma', 'gamma_s', 'speed sound', 'viscosity',
                                    'thermal conductivity', 'prandtl number', 'velocity', 'mach',
                                    'area ratio', 'I_sp', 'I_vac', 'c*', 'C_f', 'thrust', 'mass flow',
                                    'area', 'diameter', 'heat transfer coefficient', 'heat flux'], key='third')
            
            col2.selectbox('Variable:', ['heat transfer coefficient',  'temperature', 'pressure', 'density', 'specific volume',
                                    'enthalpy', 'internal energy', 'gibbs', 'entropy', 'molar mass', 'c_p',
                                    'c_v', 'gamma', 'gamma_s', 'speed sound', 'viscosity',
                                    'thermal conductivity', 'prandtl number', 'velocity', 'mach',
                                    'area ratio', 'I_sp', 'I_vac', 'c*', 'C_f', 'thrust', 'mass flow',
                                    'area', 'diameter', 'heat flux'], key='fourth')
            
            state.contour_fig = contour_plot(state.contour_properties, state.first, state.second, state.third, state.fourth)
            st.plotly_chart(state.contour_fig)