import numpy as np
from scipy.integrate import odeint
from NetworkModel import NetworkModel

class ABMSimulation:
    def __init__(self, network_model: NetworkModel, ECM_input):
        self.location = [0, 0]
        self.network_model = network_model
        # List with [tgf_beta, il_1_beta, il_6, tnf_alpha]-
        self.ECM_input = ECM_input

    def updateNetworkState(self, results_dict):

        results_dict = self.network_model.get_results_dict()

        # Set constants for latentTGFB_ABM
        kgen, ksec, latentTGFBnet, kdeg, kact = 530000, 23700, results_dict['latentTGFB'], 0.0096, 0.045

        # Set initial condition and time points for latentTGFB_ABM
        initial_latentTGFBABM, t_latentTGFB_ABM = 0.5, np.linspace(0, 10)

        # Solve the ODE for latentTGFB_ABM
        self.latentTGFB_ABM = odeint(self.ode_system_latentTGFB, initial_latentTGFBABM, t_latentTGFB_ABM,
                                     args=(kgen, ksec, latentTGFBnet, kdeg, kact))[-1][0]

        # Solve for TGFB_ABM
        self.TGFB_ABM = 0.045 * self.latentTGFB_ABM

        # Set constants for IL-1b_ABM
        kgen, kdeg = 4847, 0.277

        # Set initial condition and time points for IL-1b_ABM
        initial_IL1B, t_IL1B = 0.5, np.linspace(0, 10)

        # Solve the ODE for IL-1b_ABM
        self.IL1b_ABM = odeint(self.ode_system_IL1B, initial_IL1B, t_IL1B, args=(kgen, kdeg))[-1][0]

        # Set constants for IL6_ABM
        kgen, ksec, kdeg, IL6net = 256000, 79360, 0.277, results_dict['IL6']

        # Set initial condition and time points for IL6_ABM
        initial_IL6, t_IL6 = 0.5, np.linspace(0, 10)

        # Solve the ODE for IL6_ABM
        self.IL6_ABM = odeint(self.ode_system_IL6, initial_IL6, t_IL6, args=(kgen, ksec, IL6net, kdeg))[-1][0]

        # Set constants for TNFA_ABM
        kgen, kdeg = 895.4, 1.386

        # Set initial condition and time points for TNFA_ABM
        initial_TNFA, t_TNFA = 0.25, np.linspace(0, 10)

        # Solve the ODE for TNFA_ABM
        self.TNFA_ABM = odeint(self.ode_system_TNFA, initial_TNFA, t_TNFA, args=(kgen, kdeg))[-1][0]

        # Set constants for Collagen_ABM
        kdep, kdeg, ColIRNAnet, ColIIIRNAnet = 0.0056, 0.0035, results_dict['CImRNA'], results_dict['CIIImRNA']

        # Set initial condition and time points for Collagen_ABM
        initial_Collagen, t_Collagen = 0.25, np.linspace(0, 10)

        # Solve the ODE for Collagen_ABM
        self.Collagen_ABM = odeint(self.ode_system_Collagen, initial_Collagen, t_Collagen,
                                   args=(kdep, kdeg, ColIRNAnet, ColIIIRNAnet))[-1][0]

        # Calculate weights
        self.w_IL6 = self.IL6_ABM / (self.IL6_ABM + 462000)  # kd_IL6 = 462000
        self.w_IL1B = self.IL1b_ABM / (self.IL1b_ABM + 8750)  # kd_IL1B = 8750
        self.w_TNFA = self.TNFA_ABM / (self.TNFA_ABM + 323)  # kd_TNFA = 323
        self.w_TGFB = self.TGFB_ABM / (self.TGFB_ABM + 700)  # kd_TGFB = 700

        return [latentTGFB_ABM, TGFB_ABM, IL1b_ABM, IL6_ABM, TNFA_ABM, w_TGFB, w_TNFA,
                w_IL6, w_IL1B, Collagen_ABM]

    # Add other methods as needed

    def ode_system_latentTGFB(self, latentTGFB_ABM, t, kgen, ksec, latentTGFBnet, kdeg, kact):
        return kgen + ksec * latentTGFBnet - kdeg * latentTGFB_ABM - kact * latentTGFB_ABM

    def ode_system_IL1B(self, IL1B_ABM, t, kgen, kdeg):
        return kgen - kdeg * IL1B_ABM

    def ode_system_IL6(self, IL6_ABM, t, kgen, ksec, IL6net, kdeg):
        return kgen + ksec * IL6net - kdeg * IL6_ABM

    def ode_system_TNFA(self, TNFA_ABM, t, kgen, kdeg):
        return kgen - kdeg * TNFA_ABM

    def ode_system_Collagen(self, Collagen_ABM, t, kdep, kdeg, ColIRNAnet, ColIIIRNAnet):
        return kdep * (ColIRNAnet + ColIIIRNAnet) - kdeg
