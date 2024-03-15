import pandas as pd
import numpy as np
from scipy.integrate import ode
import NetfluxODE


class NetworkModel:
    def __init__(self, specific_time, species_names):
        self.species_names = species_names
        self.specific_time = specific_time
        self.results_dict = None
        self.params = None  # New attribute to store parameters

    def get_results_dict(self):
        return self.results_dict

    def set_params(self, tau_values, y_max_values, w_values, n_values, ec50_values):
        # Validate input lengths here if needed
        self.params = [tau_values, y_max_values, w_values, n_values, ec50_values]

    def run_simulation(self, tspan, y0):
        if self.params is None:
            raise ValueError("Parameters are not set. Call set_params method before running the simulation.")

        t = []
        dt = tspan[1] / ((self.specific_time * 150) / 10)
        r = ode(NetfluxODE.ODEfunc).set_integrator('vode', method='adams', order=10, rtol=0, atol=1e-6,
                                                   with_jacobian=False)
        r.set_initial_value(y0, tspan[0]).set_f_params(*self.params)
        results = np.empty([0, len(self.species_names)])
        while r.successful() and r.t <= tspan[1]:
            r.integrate(r.t + dt)
            results = np.append(results, [r.y], axis=0)
            t.append(r.t)
        return t, results

    def run_simulation_and_return_dict(self, tspan, y0):
        if self.params is None:
            raise ValueError("Parameters are not set. Call set_params method before running the simulation.")

        t, results = self.run_simulation(tspan, y0)
        index = int((self.specific_time - tspan[0]) / (tspan[1] / ((self.specific_time * 150) / 10))) - 1
        results_at_specific_time = results[index, :]
        self.results_dict = self.create_results_dict(results_at_specific_time)
        return self.results_dict

    def create_results_dict(self, results_at_specific_time):
        results_dict = {}
        for i in range(len(self.species_names)):
            node = self.species_names[i]
            value = results_at_specific_time[i]
            results_dict[node] = value
        return results_dict
