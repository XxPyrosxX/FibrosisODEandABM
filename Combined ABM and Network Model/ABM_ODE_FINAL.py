import numpy as np
from scipy.integrate import ode, odeint
import tkinter as tk
import NetfluxODE2


SPECIES_NAMES = ['AngII', 'AT1R', 'AGT', 'ACE', 'NOX', 'ROS', 'ET1', 'ETAR', 'DAG', 'PKC', 'TRPC', 'NE', 'BAR',
                 'Forskolin', 'AC', 'cAMP', 'PKA', 'CREB', 'CBP', 'TGFB', 'TGFB1R', 'smad3', 'smad7', 'latentTGFB',
                 'BAMBI', 'PDGF', 'PDGFR', 'NP', 'NPRA', 'cGMP', 'PKG', 'mechanical', 'B1int', 'Rho', 'ROCK', 'Ca',
                 'calcineurin', 'NFAT', 'IL6', 'gp130', 'STAT', 'IL1', 'IL1RI', 'TNFa', 'TNFaR', 'NFKB', 'PI3K',
                 'Akt', 'p38', 'TRAF', 'ASK1', 'MKK3', 'PP1', 'JNK', 'abl', 'Rac1', 'MEKK1', 'MKK4', 'ERK', 'Ras',
                 'Raf', 'MEK1', 'FAK', 'epac', 'Factin', 'FA', 'migration', 'cmyc', 'CTGF', 'proliferation', 'SRF',
                 'EDAFN', 'aSMA', 'AP1', 'TIMP1', 'TIMP2', 'PAI1', 'proMMP14', 'proMMP1', 'proMMP2', 'proMMP9',
                 'MMP1', 'MMP2', 'MMP9', 'MMP14', 'fibronectin', 'periostin', 'CImRNA', 'CIIImRNA', 'CI', 'CIII']

SPECIFIC_TIME = 1
TIME_STEPS = 10

# copied over from Netflux_ODE
def run_simulation(tspan, y0, params, speciesNames, specific_time):
    t = []
    dt = tspan[1] / ((specific_time * 150) / 10)
    r = ode(NetfluxODE2.ODEfunc).set_integrator('vode', method='adams', order=10, rtol=0, atol=1e-6, with_jacobian=False)
    r.set_initial_value(y0, tspan[0]).set_f_params(*params)
    results = np.empty([0, len(speciesNames)])
    while r.successful() and r.t <= tspan[1]:
        r.integrate(r.t + dt)
        results = np.append(results, [r.y], axis=0)
        t.append(r.t)
    return t, results


# copied over from Netflux_ODE
def create_results_dict(results_at_specific_time, speciesNames):
    results_dict = {}
    for i in range(len(speciesNames)):
        node = speciesNames[i]
        value = results_at_specific_time[i]
        results_dict[node] = value
    return results_dict


def run_simulation_and_return_dict(tspan, y0, params, specific_time, speciesNames):
    t, results = run_simulation(tspan, y0, params, speciesNames, specific_time)
    index = int((specific_time - tspan[0]) / (tspan[1] / ((specific_time * 150) / 10))) - 1
    results_at_specific_time = results[index, :]
    results_dict = create_results_dict(results_at_specific_time, speciesNames)
    return results_dict

GRID_SIZE = 10
value_layers = []

for i in range(GRID_SIZE):
    for j in range(GRID_SIZE):
        tgf_beta = (0.72 - (i * 0.08))
        il_1_beta = il_6 = tnfa = j * 0.08
        value_layers.append([tgf_beta, il_1_beta])

collagen_values = []

for x in range(len(value_layers)):
    results_dict = {}
    value_layer_1 = value_layers[x]
    tau = np.array(
        [1, 1.000000e-01, 10, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01, 1.000000e-01,
         1.000000e-01, 1, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1,
         1.000000e-01, 1.000000e-01, 10, 10, 1.000000e-01, 1, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01, 1.000000e-01,
         1, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 1.000000e-01,
         1.000000e-01, 1, 1.000000e-01, 1, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01,
         1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01,
         1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1.000000e-01, 1, 1, 10,
         1.000000e-01, 1.000000e-01, 10, 1.000000e-01, 10, 10, 1.000000e-01, 10, 10, 10, 1, 1, 1, 1, 10, 10, 10, 10, 10,
         10, 1, 1, 10, 10, ])
    ymax = np.array(
        [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5.000000e-01, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
         1, 1, 1, 1, 5.000000e-01, 5.000000e-01, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ])
    y0 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ]
    n = np.array([1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, 1.400000e+00, ])
    EC50 = np.array([6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, 6.000000e-01, ])
    w = [2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01, 2.500000e-01,
     2.500000e-01, 2.500000e-01, 2.500000e-01, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, ]
    input_parameters = [0, 0, 0, 0]
    speciesNames = ['AngII', 'AT1R', 'AGT', 'ACE', 'NOX', 'ROS', 'ET1', 'ETAR', 'DAG', 'PKC', 'TRPC', 'NE', 'BAR',
                'Forskolin', 'AC', 'cAMP', 'PKA', 'CREB', 'CBP', 'TGFB', 'TGFB1R', 'smad3', 'smad7', 'latentTGFB',
                'BAMBI', 'PDGF', 'PDGFR', 'NP', 'NPRA', 'cGMP', 'PKG', 'mechanical', 'B1int', 'Rho', 'ROCK', 'Ca',
                'calcineurin', 'NFAT', 'IL6', 'gp130', 'STAT', 'IL1', 'IL1RI', 'TNFa', 'TNFaR', 'NFKB', 'PI3K', 'Akt',
                'p38', 'TRAF', 'ASK1', 'MKK3', 'PP1', 'JNK', 'abl', 'Rac1', 'MEKK1', 'MKK4', 'ERK', 'Ras', 'Raf',
                'MEK1', 'FAK', 'epac', 'Factin', 'FA', 'migration', 'cmyc', 'CTGF', 'proliferation', 'SRF', 'EDAFN',
                'aSMA', 'AP1', 'TIMP1', 'TIMP2', 'PAI1', 'proMMP14', 'proMMP1', 'proMMP2', 'proMMP9', 'MMP1', 'MMP2',
                'MMP9', 'MMP14', 'fibronectin', 'periostin', 'CImRNA', 'CIIImRNA', 'CI', 'CIII', 'latentTGFB_ABM',
                'IL1_ABM', 'IL6_ABM', 'TNFa_ABM', 'Collagen']
    w_IL6 = 0
    w_IL1B = 0
    w_TNFA = 0
    w_TGFB = 0
    collagen = 0.25

    for x in range(TIME_STEPS):
        if x == 0:
            w[1] = value_layer_1[0] # TGFB
            w[3] = value_layer_1[1] # IL6
            w[4] = value_layer_1[1] # IL1
            w[5] = value_layer_1[1] # TNFa

        else:
            w[1] = w_IL6  # TGFB
            w[3] = w_IL6  # IL6
            w[4] = w_IL1B  # IL1
            w[5] = w_TNFA  # TNFa
            y0 = list(results_dict.values())


        results_dict = run_simulation_and_return_dict([0, SPECIFIC_TIME], y0,
                                                      [tau, ymax, w, n, EC50],
                                                      SPECIFIC_TIME, SPECIES_NAMES)
        ########################################################################################################################
        # Differential Equation Solution latentTGFB_ABM (5)
        def ode_system_latentTGFB(latentTGFB_ABM, t, kgen, ksec, latentTGFBnet, kdeg, kact):
            return kgen + ksec * latentTGFBnet - kdeg * latentTGFB_ABM - kact * latentTGFB_ABM

        # Set constants
        kgen, ksec, latentTGFBnet, kdeg, kact = 530000, 23700, results_dict['latentTGFB'], 0.0096, 0.045

        # Set initial condition and time points
        initial_latentTGFBABM, t = input_parameters[3], np.linspace(0, 10)

        # Solve the ODE
        latentTGFB_ABM = odeint(ode_system_latentTGFB, initial_latentTGFBABM, t, args=(kgen, ksec, latentTGFBnet, kdeg, kact))[-1][0]
        ########################################################################################################################
        # Equation Solution TGFB_ABM (6)
        TGFB_ABM = 0.045 * latentTGFB_ABM
        ########################################################################################################################
        # Differential Equation Solution IL-1B (7)
        def ode_system_IL1B(IL1B_ABM, t, kgen, kdeg):
            return kgen - kdeg * IL1B_ABM

        # Set constants
        kgen, kdeg = 4847, 0.277

        # Set initial condition and time points
        initial_IL1B, t = input_parameters[1], np.linspace(0, 10)

        # Solve the ODE
        IL1b_ABM = odeint(ode_system_IL1B, initial_IL1B, t, args=(kgen, kdeg))[-1][0]
        ########################################################################################################################
        # Differential Equation Solution IL-6_ABM (8)
        def ode_system_IL6(IL6_ABM, t, kgen, ksec, IL6net, kdeg):
            return kgen + ksec * IL6net - kdeg * IL6_ABM

        # Set constants
        kgen, ksec, kdeg, IL6net = 256000, 79360, 0.277, results_dict['IL6']

        # Set initial condition and time points
        initial_IL6, t = input_parameters[0], np.linspace(0, 10)

        # Solve the ODE
        IL6_ABM = odeint(ode_system_IL6, initial_IL6, t, args=(kgen, ksec, IL6net, kdeg))[-1][0]
        ########################################################################################################################
        # Differential Equation Solution TNFA_ABM (9)
        def ode_system_TNFA(TNFA_ABM, t, kgen, kdeg):
            return kgen - kdeg * TNFA_ABM

        # Set constants
        kgen, kdeg = 895.4, 1.386

        # Set initial condition and time points
        initial_TNFA, t = input_parameters[2], np.linspace(0, 10)

        # Solve the ODE
        TNFA_ABM = odeint(ode_system_TNFA, initial_TNFA, t, args=(kgen, kdeg))[-1][0]
        ########################################################################################################################
        # Differential Equation Solution Collagen (10)
        def ode_system_Collagen(Collagen_ABM, t, kdep, kdeg, ColIRNAnet, ColIIIRNAnet):
            return kdep * (ColIRNAnet + ColIIIRNAnet) - kdeg * Collagen_ABM

        # Set constants
        kdep, kdeg, ColIRNAnet, ColIIIRNAnet = 0.0056, 0.0035, results_dict['CImRNA'], results_dict['CIIImRNA']

        # Set initial condition and time points
        initial_Collagen, t = collagen, np.linspace(0, 10)

        # Solve the ODE
        Collagen_ABM = odeint(ode_system_Collagen, initial_Collagen, t, args=(kdep, kdeg, ColIRNAnet, ColIIIRNAnet))[-1][0]
        collagen = Collagen_ABM
        ########################################################################################################################
        # Weight IL-6 (1)
        w_IL6 = IL6_ABM / (IL6_ABM + 462000) # kd_IL6 = 462000
        ########################################################################################################################
        # Weight IL-1B (2)
        w_IL1B = IL1b_ABM / (IL1b_ABM + 8750) # kd_IL1B = 8750
        ########################################################################################################################
        # Weight TNFA (3)
        w_TNFA = TNFA_ABM / (TNFA_ABM + 323) # kd_TNFA = 323
        ########################################################################################################################
        # Weight TGFB (3)
        w_TGFB = TGFB_ABM / (TGFB_ABM + 700) # kd_TGFB = 700
        # ----------------------------------------------------------------------------------------------------------------------
    collagen_values.append(collagen)

def scale_color(value):
    # Convert value from 0 to 1 to a shade of grey (0:black, 0.6:white)
    shade = int(float(value * 152))  # Scale to range from 0 to 153 (0.6 * 255)
    color = "#{:02x}{:02x}{:02x}".format(shade, shade, shade)  # Convert to hexadecimal color code
    return color

def create_grid(values):
    root = tk.Tk()
    root.title("Grid Display")

    for i in range(GRID_SIZE):
        for j in range(GRID_SIZE):
            if len(values) > 0:
                value = float(values.pop(0))
                color = scale_color(value)
            else:
                value = 0.0
                color = "white"
            label = tk.Label(root, text=str(value), borderwidth=1, relief="solid", width=10, height=5, bg=color)
            label.grid(row=i, column=j)

    root.mainloop()

# Example list
rounded_values = [ '%.4f' % elem for elem in collagen_values]
print(rounded_values)
create_grid(rounded_values)
