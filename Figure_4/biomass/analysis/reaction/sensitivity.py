import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from biomass.exec_model import ExecModel
from biomass.analysis import get_signaling_metric, dlnyi_dlnxj

class ReactionSensitivity(ExecModel):
    """ Sensitivity for rate equations
    """
    def __init__(self, model):
        super().__init__(model)

    def _calc_sensitivity_coefficients(self, metric, reaction_indices):
        """ Calculating Sensitivity Coefficients

        Parameters
        ----------
        metric : str
            - 'maximum': The maximum value.
            - 'minimum' : The minimum value.
            - 'duration': The time it takes to decline below 10% of its maximum.
            - 'integral': The integral of concentration over the observation time.

        nonzero_indices : list of int
            for i in nonzero_indices:
                y0[i] != 0.0

        Returns
        -------
        sensitivity_coefficients : numpy array

        """
        rate = 1.01  # 1% change
        n_file = self.get_executable()
        signaling_metric = np.full(
            (
                len(n_file),
                len(reaction_indices)+1,
                len(self.obs),
                len(self.sim.conditions)
            ), np.nan
        )
        for i, nth_paramset in enumerate(n_file): #i: parameter set
            (x, y0) = self.load_param(nth_paramset)
            for j, rxn_idx in enumerate(reaction_indices): #j: reaction_indices
                perturbation = {}
                for idx in reaction_indices:
                    perturbation[idx] = 1
                perturbation[rxn_idx] = rate
                if self.sim.simulate(x, y0, perturbation) is None:
                    for k, _ in enumerate(self.obs):
                        for l, _ in enumerate(self.sim.conditions):
                            signaling_metric[i, j, k, l] = \
                                get_signaling_metric(
                                    metric, self.sim.simulations[k, :, l]
                                )
                sys.stdout.write(
                    '\r{:d} / {:d}'.format(
                        i*len(reaction_indices)+j+1,
                        len(n_file)*len(reaction_indices)
                    )
                )
            if self.sim.simulate(x, y0) is None:
                for k, _ in enumerate(self.obs):
                    for l, _ in enumerate(self.sim.conditions):
                        signaling_metric[i, -1, k, l] = \
                            get_signaling_metric(
                                metric, self.sim.simulations[k, :, l]
                            )
        sensitivity_coefficients = dlnyi_dlnxj(
            signaling_metric, n_file, reaction_indices,
            self.obs, self.sim.conditions, rate
        )

        return sensitivity_coefficients

    # 1106 def
    def _calc_dual_inhibition_scores(self, metric, reaction_indices):
        """ Calculating dual_inhibition_scores on the metric

        Parameters
        ----------
        metric : str
            - 'maximum': The maximum value.
            - 'minimum' : The minimum value.
            - 'duration': The time it takes to decline below 10% of its maximum.
            - 'integral': The integral of concentration over the observation time.
            etc ... now i use only 'mean_integral_after_3w'

        nonzero_indices : list of int
            for i in nonzero_indices:
                y0[i] != 0.0

        reaction_indices??????????????? reaction_index1, reaction_index2 ??????2????????????->?????????

        Returns
        -------
        dual_inhibition_scores : 0~200%???dual moduration 11 x 11 numpy array

        """
        inhibition_rate = np.arange(0,11,1)/10 #0~100%
        n_file = self.get_executable()

        inhibition_metric = np.full(
            (
                len(n_file), #num of paramset
                len(reaction_indices), # reactions 1
                len(reaction_indices), # reactions 2
                len(self.obs), #num of obj.s
                len(self.sim.conditions),  #num of conditions
                len(inhibition_rate), #%inhibition score num1
                len(inhibition_rate), #%inhibition score num2
            ), np.nan
        )

        for i, nth_paramset in enumerate(n_file): #i: parameter set
            (x, y0) = self.load_param(nth_paramset)
            for j1, rxn_idx1 in enumerate(reaction_indices): #j1: reaction_indices
                for j2, rxn_idx2 in enumerate(reaction_indices): #j1: reaction_indices
                    perturbation = {}
                    for idx in reaction_indices:
                        perturbation[idx] = 1
                    for a, inh1 in enumerate(inhibition_rate):
                        for b, inh2 in enumerate(inhibition_rate):
                            perturbation[rxn_idx1] = inh1
                            perturbation[rxn_idx2] = inh2
                            perturbation[6] = 3.0
                            if self.sim.simulate(x, y0, perturbation) is None:
                                for k, _ in enumerate(self.obs):
                                    for l, _ in enumerate(self.sim.conditions):
                                        inhibition_metric[i, j1, j2, k, l, a, b] = \
                                            get_signaling_metric(
                                                metric, self.sim.simulations[k, :, l]
                                            )
                            sys.stdout.write(
                                '\r{:d} / {:d}'.format(
                                    i*len(reaction_indices)**2+j1*len(reaction_indices)+j2+1,
                                    len(n_file)*len(reaction_indices)**2
                                )
                            )

        return inhibition_metric


    def _load_sc(self, metric, reaction_indices):
        os.makedirs(
            self.model_path + '/figure/sensitivity/' \
            'reaction/{}/heatmap'.format(metric), exist_ok=True
        )
        if not os.path.isfile(
                self.model_path + '/sensitivity_coefficients/' \
                'reaction/{}/sc.npy'.format(metric)):
            os.makedirs(
                self.model_path + '/sensitivity_coefficients/' \
                'reaction/{}'.format(metric), exist_ok=True
            )
            sensitivity_coefficients = \
                self._calc_sensitivity_coefficients(metric, reaction_indices)
            np.save(
                self.model_path + '/sensitivity_coefficients/' \
                'reaction/{}/sc'.format(metric), sensitivity_coefficients
            )
        else:
            sensitivity_coefficients = np.load(
                self.model_path + '/sensitivity_coefficients/' \
                'reaction/{}/sc.npy'.format(metric)
            )

        return sensitivity_coefficients

    @staticmethod
    def _draw_vertical_span(biological_processes, width):
        if len(biological_processes) > 1:
            left_end = 0
            for i, proc in enumerate(biological_processes):
                if i % 2 == 0:
                    plt.axvspan(
                        left_end - width,
                        left_end - width + len(proc),
                        facecolor='k', alpha=0.1
                    )
                left_end += len(proc)

    def _write_reaction_indices(self, reaction_indices, average, stdev, width):
        distance = np.max(average) * 0.01 #0.05
        for i, j in enumerate(reaction_indices):
            xp = i + width * 0.5 * (len(self.sim.conditions) - 1)
            yp = average[i, np.argmax(np.abs(average[i, :]))]
            yerr = stdev[i, np.argmax(stdev[i, :])]
            if yp > 0:
                plt.text(
                    xp, yp + yerr + distance, "v" + str(j),
                    #xp, yp + yerr + distance, str(j),
                    ha='center', va='bottom', fontsize=10, rotation=0
                )
            else:
                plt.text(
                    #xp, yp - yerr - distance, str(j),
                    xp, yp - yerr - distance, "v" + str(j),
                    ha='center', va='top', fontsize=10, rotation=0
                )

    def _barplot_sensitivity(
            self,
            metric,
            sensitivity_coefficients,
            biological_processes,
            reaction_indices,
    ):
        options = self.viz.sensitivity_options

        # rcParams
        self.viz.set_sensitivity_rcParams()

        if len(options['cmap']) < len(self.sim.conditions):
            raise ValueError(
                "len(sensitivity_options['cmap']) must be equal to"
                " or greater than len(sim.conditions)."
            )
        for k, obs_name in enumerate(self.obs):
            plt.figure(figsize=options['figsize'])
            self._draw_vertical_span(biological_processes, options['width'])

            sensitivity_array = sensitivity_coefficients[:, :, k, :]
            # Remove NaN
            nan_idx = []
            for i in range(sensitivity_array.shape[0]):
                for j in range(sensitivity_array.shape[1]):
                    if any(np.isnan(sensitivity_array[i, j, :])):
                        nan_idx.append(i)
            #20201008 num to zero
            #sensitivity_array = np.nan_to_num(sensitivity_array)
            #print("sensitivity_array.shape: " + str(sensitivity_array.shape))
            sensitivity_array = np.delete(
                sensitivity_array, nan_idx, axis=0
            )


            if sensitivity_array.size != 0:
                average = np.mean(sensitivity_array, axis=0)
                if sensitivity_array.shape[0] == 1:
                    stdev = np.zeros(
                        (sensitivity_array.shape[1], sensitivity_array.shape[2])
                    )
                else:
                    stdev = np.std(sensitivity_array, axis=0, ddof=1)
                for l, condition in enumerate(self.sim.conditions):
                    # 201008 skip Contorl conditions
                    # if condition == 'Control':
                    #    continue
                    # print("x: ", len(reaction_indices))
                    #print(len(np.arange(len(reaction_indices))))
                    # print("y: ", average.shape)

                    plt.bar(
                        np.arange(len(reaction_indices)) + l * options['width'],
                        average[:, l], yerr=stdev[:, l], ecolor=options['cmap'][l],
                        capsize=2, width=options['width'],
                        color=options['cmap'][l], align='center', label=condition
                    )
                self._write_reaction_indices(
                    reaction_indices, average, stdev, options['width']
                )
                plt.hlines(
                    [0], -options['width'], len(reaction_indices), 'k', lw=1
                )
                plt.xticks([])
                plt.ylabel(
                    'Control coefficients on\n'+metric +
                    ' (' + obs_name.replace('_', ' ') + ')'
                )
                plt.xlim(-options['width'], len(reaction_indices))
                #210603added ylim
                plt.ylim(
                    [np.min([-0.1, np.min(average - stdev)])*1.3,
                     np.max([0.1, np.max(average + stdev)])*1.3,
                     ]
                    )
                plt.legend(loc=options['legend_loc'], frameon=False)
                plt.savefig(
                    self.model_path + '/figure/sensitivity/reaction/'\
                    '{}/{}.pdf'.format(
                        metric, obs_name
                    ), bbox_inches='tight'
                )
                plt.close()

    @staticmethod
    def _remove_nan(sensitivity_matrix, normalize):
        nan_idx = []
        for i in range(sensitivity_matrix.shape[0]):
            if any(np.isnan(sensitivity_matrix[i, :])):
                nan_idx.append(i)
            else:
                pass
            if np.nanmax(np.abs(sensitivity_matrix[i, :])) == 0.0:
                sensitivity_matrix[i, :] = np.zeros(
                    sensitivity_matrix.shape[1]
                )
            else:
                sensitivity_matrix[i, :] = sensitivity_matrix[i, :] / (
                    np.nanmax(
                        np.abs(sensitivity_matrix[i, :])
                    ) if normalize else 1
                )

        return np.delete(sensitivity_matrix, nan_idx, axis=0)


    @staticmethod
    def _remove_nan(sensitivity_matrix, normalize):
        nan_idx = []
        for i in range(sensitivity_matrix.shape[0]):
            if any(np.isnan(sensitivity_matrix[i, :])):
                nan_idx.append(i)
            else:
                pass
            if np.nanmax(np.abs(sensitivity_matrix[i, :])) == 0.0:
                sensitivity_matrix[i, :] = np.zeros(
                    sensitivity_matrix.shape[1]
                )
            else:
                sensitivity_matrix[i, :] = sensitivity_matrix[i, :] / (
                    np.nanmax(
                        np.abs(sensitivity_matrix[i, :])
                    ) if normalize else 1
                )

        return np.delete(sensitivity_matrix, nan_idx, axis=0)

    def _heatmap_sensitivity(
            self,
            metric,
            sensitivity_coefficients,
            biological_processes,
            reaction_indices,
    ):
        options = self.viz.sensitivity_options
        # rcParams
        self.viz.set_sensitivity_rcParams()

        for k, obs_name in enumerate(self.obs):
            for l, condition in enumerate(self.sim.conditions):
                sensitivity_matrix = self._remove_nan(
                    sensitivity_coefficients[:, :, k, l],
                    normalize=False
                )
                if sensitivity_matrix.shape[0] > 1 and \
                        not np.all(sensitivity_matrix == 0.0):
                    sns.clustermap(
                        data=sensitivity_matrix,
                        center=0,
                        robust=True,
                        method='ward',
                        cmap='RdBu_r',
                        linewidth=.5,
                        col_cluster=False,
                        figsize=options['figsize'],
                        xticklabels=[str(j) for j in reaction_indices],
                        yticklabels=[],
                        #cbar_kws={"ticks": [-1, 0, 1]}
                    )
                    plt.savefig(
                        self.model_path + '/figure/sensitivity/reaction/'\
                        '{}/heatmap/{}_{}.pdf'.format(
                            metric, condition, obs_name
                        ), bbox_inches='tight'
                    )
                    plt.close()

    def analyze(self, metric, style):
        if not self.rxn.reactions:
            raise ValueError(
                'Define reaction indices (reactions) in reaction_network.py'
            )
        biological_processes = self.rxn.group()
        reaction_indices = np.sum(biological_processes, axis=0)
        #print("reaction_indices :", reaction_indices)
        sensitivity_coefficients = self._load_sc(metric, reaction_indices)
        #print("sensitivity_coefficients.shape :", sensitivity_coefficients.shape)

        if style == 'barplot':
            self._barplot_sensitivity(
                metric, sensitivity_coefficients, biological_processes,
                reaction_indices
            )
        elif style == 'heatmap':
            self._heatmap_sensitivity(
                metric, sensitivity_coefficients, biological_processes,
                reaction_indices
            )
        else:
            raise ValueError("Available styles are: 'barplot', 'heatmap'")
