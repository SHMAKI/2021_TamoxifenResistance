import numpy as np
from scipy.spatial.distance import cosine

from .observable import observables, ExperimentalData, NumericalSimulation
from .set_search_param import SearchParam


def _compute_objval_rss(sim_data, exp_data):
    """Return Residual Sum of Squares
    """
    #return np.dot((sim_data-exp_data), (sim_data-exp_data)) # Mean squared error
    return np.dot((sim_data-exp_data) / (exp_data + 0.1), (sim_data-exp_data) / (exp_data + 0.1))# Mean squared percemtageerror
    #return np.dot((sim_data-exp_data) / (exp_data), (sim_data-exp_data) / (exp_data))# Mean squared percemtageerror

def _compute_objval_cos(sim_data, exp_data):
    """Return Cosine distance
    """
    return cosine(sim_data, exp_data)


def _diff_sim_and_exp(
        sim_matrix,
        exp_dict,
        exp_timepoint,
        conditions,
        sim_norm_max
):
    sim_val = []
    exp_val = []

    for idx, condition in enumerate(conditions):
        if condition in exp_dict.keys():
            sim_val.extend(sim_matrix[exp_timepoint, idx])
            exp_val.extend(exp_dict[condition])

    return np.array(sim_val) / sim_norm_max, np.array(exp_val)


def objective(indiv_gene, *args):
    """Define an objective function to be minimized
    """
    if len(args) == 0:
        sp = SearchParam()
        indiv = sp.gene2val(indiv_gene)
        (x, y0) = sp.update(indiv)
    elif len(args) == 1:
        raise ValueError('not enough values to unpack (expected 2, got 1)')
    elif len(args) == 2:
        (x, y0) = args
    else:
        raise ValueError('too many values to unpack (expected 2)')

    exp = ExperimentalData()
    sim = NumericalSimulation()

    if sim.simulate(x, y0) is None:
        error = np.zeros(len(observables))
        for i, _ in enumerate(observables):
            if exp.experiments[i] is not None:
                error[i] = _compute_objval_rss(
                    *_diff_sim_and_exp(
                        sim.simulations[i], exp.experiments[i],
                        exp.get_timepoint(i), sim.conditions,
                        sim_norm_max=1 if not sim.normalization \
                            else np.max(sim.simulations[i])
                    )
                )
        return np.sum(error)
    else:
        return np.inf
