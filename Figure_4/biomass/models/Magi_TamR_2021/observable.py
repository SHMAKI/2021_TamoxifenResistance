
import numpy as np
import math
from scipy.integrate import ode

from .name2idx import C, V
from .set_model import DifferentialEquation


observables = [
    "growth_rate",
    "Sensitive_Cells",
    "Pre-resistant_Cells",
    'Resistant_Cells1',
    'Resistant_Cells2',
    # "Total_Cell_number (log10)",
    # "Total_R1_number",
    # "Total_R2_number",
]

class NumericalSimulation(DifferentialEquation):
    """ Simulate a model using scipy.integrate.ode

    Attributes
    ----------
    normalization : bool
        if True, simulation results in each observable are divided by their
        maximum values.

    """
    def __init__(self):
        super().__init__(perturbation={})
        self.normalization = False

    point = 1
    t = np.arange(0.0, 10+10**-point, 10**-point) #0~10 week

    # Experimental conditions
    conditions = ['Control', 'Tamoxifen']

    simulations = np.empty((len(observables), len(t), len(conditions)))

    def simulate(self, x, y0, _perturbation={}):
        if _perturbation:
            self.perturbation = _perturbation
        for i, condition in enumerate(self.conditions):
            if condition == 'Control':
                y0[V.Tam] = 0
                #pass
            elif condition == 'Tamoxifen':
                y0[V.Tam] = 1

            (T, Y) = self._solveode(self.diffeq, y0, self.t, tuple(x)) #tspan = t

#            if T[-1] < self.t[-1]:
            if np.round(T[-1], self.point) != np.round(self.t[-1], self.point) or len(T) != len(self.t):
                return False
            else:
                # calc. growth rate from week1
                span = 10**self.point # delta=0.1, week:1
                delta_t = 10**-self.point #week
                tmp_d = Y[:(Y.shape[0]-span), V.S] + Y[:(Y.shape[0]-span), V.P] + Y[:(Y.shape[0]-span), V.R1] + Y[:(Y.shape[0]-span), V.R2]
                val = (Y[span:, V.S] + Y[span:, V.P] + Y[span:, V.R1] + Y[span:, V.R2]) / tmp_d

                # calc. growth rate from week0 to week1
                # calc. cell number at day -1~-span
                tm1_s = Y[0, V.S] * math.e ** (x[C.v1_k] * np.arange(-span*delta_t, 0, delta_t)) + \
                        Y[0, V.P] * math.e ** (x[C.v5_k] * np.arange(-span*delta_t, 0, delta_t)) + \
                        Y[0, V.R1] * math.e ** (x[C.v9_k] * np.arange(-span*delta_t, 0, delta_t)) + \
                        Y[0, V.R2] * math.e ** (x[C.v12_k] * np.arange(-span*delta_t, 0, delta_t)) #v12_k_tam
                # calc growth rate from week0 to week1
                val1_s = (Y[:span, V.S] + Y[:span, V.P] + Y[:span, V.R1] + Y[:span, V.R2]) / tm1_s

                # all cell pop.
                all_pop = Y[:, V.S] + Y[:, V.P] + Y[:, V.R1] + Y[:, V.R2]

                self.simulations[observables.index('growth_rate'), :, i] = (
                    np.hstack((val1_s, val)) # growth rate of w0tow1, w1~
                )
                self.simulations[observables.index('Sensitive_Cells'), :, i] = (
                    Y[:, V.S] / all_pop
                )
                self.simulations[observables.index('Pre-resistant_Cells'), :, i] = (
                    Y[:, V.P] / all_pop
                )
                self.simulations[observables.index('Resistant_Cells1'), :, i] = (
                    Y[:, V.R1] / all_pop
                )
                self.simulations[observables.index('Resistant_Cells2'), :, i] = (
                    Y[:, V.R2] / all_pop
                )
                # self.simulations[observables.index('Total_Cell_number (log10)'), :, i] = (
                #     np.log10(all_pop)
                # )
                # self.simulations[observables.index('Total_R1_number'), :, i] = (
                #     Y[:, V.R1]
                # )
                # self.simulations[observables.index('Total_R2_number'), :, i] = (
                #     Y[:, V.R2]
                # )

    def _solveode(self, diffeq, y0, tspan, args):
        dt = (self.t[-1] - self.t[0]) / (len(self.t) - 1)
        sol = ode(diffeq)
        sol.set_integrator(
            'vode', with_jacobian=True, method= 'bdf', #adams for non-stiff,
            atol=1e-9, rtol=1e-9, min_step=1e-8
        )
        sol.set_initial_value(y0, tspan[0])
        sol.set_f_params(args)

        T = [tspan[0]]
        Y = [y0]

        # 210727 prevent len(self.t) != len(T)
        rep = 1
        while sol.successful() and sol.t < tspan[-1]  and rep < len(tspan):
            sol.integrate(sol.t+dt)
            T.append(sol.t)
            Y.append(sol.y)
            rep += 1

        return np.array(T), np.array(Y)

    def _get_steady_state(self, diffeq, y0, args, eps=1e-6):
        """
        Run until a time t for which the maximal absolutevalue of the
        regularized relative derivative was smaller than eps.
        """
        while True:
            (T, Y) = self._solveode(diffeq, y0, range(2), args) #tspan = range2
            if T[-1] < 1:
                return None
            elif np.max(np.abs((Y[-1, :] - y0) / (np.array(y0) + eps))) < eps:
                break
            else:
                y0 = Y[-1, :].tolist()

        return y0

class ExperimentalData(object):
    def __init__(self):
        pass

    experiments = [None] * len(observables)
    error_bar = [None] * len(observables)

    t2 = np.arange(11) # (Unit: 1/10 week.)

    experiments[observables.index('growth_rate')] = {
        'Control': [3.27, 3.44, 3.49, 3.3, 3.36, 2.87, 3.2, 3.48, 3.4, 3.3, 3.65], #/Volumes/GoogleDrive/マイドライブ/tamR論文/190710_scRNAseqPaper/fig/Fig1b_gworth/mean.csv
        'Tamoxifen': [3.27, 2.29, 1.56, 0.99, 1.04, 1.04, 1.23, 1.84, 1.83, 1.74, 2.31],
    }
    error_bar[observables.index('growth_rate')] = { #sandard error
        'Control': [0.15, 0.53, 0.3, 0.21, 0.03, 0.24, 0.27, 0.08, 0.27, 0.24, 0.25],
        'Tamoxifen': [0.15, 0.24, 0.1, 0.05, 0.04, 0.04, 0.28, 0.23, 0.05, 0.15, 0.02],
    }

    t3  = np.arange(0,10,3)

    experiments[observables.index('Sensitive_Cells')] = {
        'Control': [0.989247312, 0.989247312, 0.989247312, 0.989247312], #/Volumes/Fast SSD/analysis_190801_mito25/results6_monocle3_n1500/3d/multistate_class_week.csv
        'Tamoxifen': [0.989247312, 0.010582011, 0.025423729, 0.011904762],
    }

    experiments[observables.index('Pre-resistant_Cells')] = {
        #'Control': [0.0, 0.0, 0.0, 0.0],
        'Tamoxifen': [0.0, 0.772486772, 0.610169492, 0.178571429],
    }

    experiments[observables.index('Resistant_Cells1')] = {
        #'Control': [0.0, 0.0, 0.0, 0.0],
        'Tamoxifen': [0.0, 0.042328042, 0.237288136, 0.30952381],
    }

    experiments[observables.index('Resistant_Cells2')] = {
        'Control': [0.010752688, 0.010752688, 0.010752688, 0.010752688],
        'Tamoxifen': [0.010752688, 0.174603175, 0.127118644, 0.5],
    }

    def get_timepoint(self, obs_idx):
        '''
        return list(map(int, exp_t))
        '''
        if obs_idx in [
            observables.index('growth_rate'),
        ]:
            exp_t = self.t2 * 10 #210727
        else:
            exp_t = self.t3 * 10 #210727

        return list(map(int, exp_t))
        #pass

    def get_timepoint_exp(self, obs_idx):
        '''
        return list(map(int, exp_t))
        '''
        if obs_idx in [
            observables.index('growth_rate'),
        ]:
            exp_t = self.t2
        else:
            exp_t = self.t3

        return list(map(int, exp_t))
        #pass
