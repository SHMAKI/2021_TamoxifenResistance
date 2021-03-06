import numpy as np

from .name2idx import C, V
from .set_model import param_values, initial_values


class SearchParam(object):
    """ Specify model parameters and/or initial values to optimize
    """
    # parameters
    idx_params = [
    C.v1_k,
    C.v1_i,
    C.v2_c,
    C.v2_a,
    C.v2_b,
    C.v3_c,
    C.v4_k,
    C.v5_k,
    C.v6_k,
    C.v7_c,
    C.v7_a,
    C.v7_b,
    C.v8_c,
    C.v9_k,
    C.v10_c,
    C.v10_a,
    C.v10_b,
    C.v11_c,
    C.v12_k,
    ]

    # initial values
    idx_initials = [
        # V.(specie)
        #V.R2
    ]

    def get_region(self):
        x = param_values()
        y0 = initial_values()

        search_param = self._init_search_param(x, y0)

        search_rgn = np.zeros((2, len(x)+len(y0)))
        # Default: 0.1 ~ 10
        for i, j in enumerate(self.idx_params):
            search_rgn[0, j] = search_param[i] * 0.1  # lower bound 0.1
            search_rgn[1, j] = search_param[i] * 10.  # upper bound 10
        # Default: 0.5 ~ 2
        for i, j in enumerate(self.idx_initials):
            search_rgn[0, j+len(x)] = \
                search_param[i+len(self.idx_params)] * 0.5  # lower bound 0.5
            search_rgn[1, j+len(x)] = \
                search_param[i+len(self.idx_params)] * 2.0  # upper bound 2

        # search_rgn[:,C.parameter] = [lower_bound, upper_bound]
        # search_rgn[:,V.specie+len(x)] = [lower_bound, upper_bound]
        search_rgn[:,C.v1_k] = [0.1, 2.0]
        search_rgn[:,C.v5_k] = [0.1, 1.3]
        search_rgn[:,C.v9_k] = [0.1, 1.3]
        search_rgn[:,C.v12_k] = [0.1, 1.3]

        search_rgn[:,C.v2_c] = [0.1, 100]
        search_rgn[:,C.v7_c] = [0.1, 100]
        search_rgn[:,C.v10_c] = [0.1, 100]

        search_rgn[:,C.v1_i] = [0.01, 0.99]

        search_rgn[:,C.v2_b] = [0.1, 10.0]
        search_rgn[:,C.v7_b] = [0.1, 10.0]
        search_rgn[:,C.v10_b] = [0.1, 10.0]

        search_rgn[:,C.v2_a] = [0.01, 5.0]
        search_rgn[:,C.v7_a] = [0.01, 5.0]
        search_rgn[:,C.v10_a] = [0.01, 5.0]

        search_rgn = self._conv_lin2log(search_rgn)

        return search_rgn

    def update(self, indiv):
        x = param_values()
        y0 = initial_values()

        for i, j in enumerate(self.idx_params):
            x[j] = indiv[i]
        for i, j in enumerate(self.idx_initials):
            y0[j] = indiv[i+len(self.idx_params)]

        # constraints --------------------------------------------------------------
        x[C.v4_k] = x[C.v6_k]

        #------------------------------------------
        return x, y0

    def gene2val(self, indiv_gene):
        search_rgn = self.get_region()
        indiv = 10**(
            indiv_gene * (
                search_rgn[1, :] - search_rgn[0, :]
            ) + search_rgn[0, :]
        )

        return indiv

    def val2gene(self, indiv):
        search_rgn = self.get_region()
        indiv_gene = (
            np.log10(indiv) - search_rgn[0, :]
        ) / (
            search_rgn[1, :] - search_rgn[0, :]
        )

        return indiv_gene

    def _init_search_param(self, x, y0):
        """Initialize search_param
        """
        if len(self.idx_params) != len(set(self.idx_params)):
            raise ValueError(
                'Duplicate parameters (C.): {}'.format(
                    [C.NAMES[idx]
                    for idx in [name for name in set(self.idx_params)
                                if self.idx_params.count(name) > 1]]
                )
            )
        elif len(self.idx_initials) != len(set(self.idx_initials)):
            raise ValueError(
                'Duplicate species (V.): {}'.format(
                    [V.NAMES[idx]
                    for idx in [name for name in set(self.idx_initials)
                                if self.idx_initials.count(name) > 1]]
                )
            )
        search_param = np.empty(
            len(self.idx_params) + len(self.idx_initials)
        )
        for i, j in enumerate(self.idx_params):
            search_param[i] = x[j]
        for i, j in enumerate(self.idx_initials):
            search_param[i+len(self.idx_params)] = y0[j]

        if np.any(search_param == 0.):
            message = 'search_param must not contain zero.'
            for idx in self.idx_params:
                if x[int(idx)] == 0.:
                    raise ValueError(
                        '"C.{}" in idx_params: '.format(
                            C.NAMES[int(idx)]
                        ) + message
                    )
            for idx in self.idx_initials:
                if y0[int(idx)] == 0.:
                    raise ValueError(
                        '"V.{}" in idx_initials: '.format(
                            V.NAMES[int(idx)]
                        ) + message
                    )

        return search_param

    def _conv_lin2log(self, search_rgn):
        """Convert Linear scale to Logarithmic scale
        """
        ##commented out for allowing negative
        for i in range(search_rgn.shape[1]):
            if np.min(search_rgn[:, i]) < 0.0:
                msg = 'search_rgn[lower_bound, upper_bound] must be positive.'
                if i <= C.NUM:
                    raise ValueError(
                        '"C.{}": '.format(C.NAMES[i]) + msg
                    )
                else:
                    raise ValueError(
                        '"V.{}": '.format(V.NAMES[i-C.NUM]) + msg
                    )
            elif np.min(search_rgn[:, i]) == 0 and np.max(search_rgn[:, i]) > 0:
                msg = 'lower_bound must be larger than 0.'
                if i <= C.NUM:
                    raise ValueError(
                        '"C.{}" '.format(C.NAMES[i]) + msg
                    )
                else:
                    raise ValueError(
                        '"V.{}" '.format(V.NAMES[i-C.NUM]) + msg
                    )
            elif search_rgn[1, i] - search_rgn[0, i] < 0.0:
                msg = 'lower_bound must be smaller than upper_bound.'
                if i <= C.NUM:
                    raise ValueError(
                        '"C.{}" : '.format(C.NAMES[i]) + msg
                    )
                else:
                    raise ValueError(
                        '"V.{}" : '.format(V.NAMES[i-C.NUM]) + msg
                    )
        difference = list(
            set(
                np.where(
                    np.any(search_rgn != 0., axis=0)
                )[0]
            ) ^ set(
                np.append(
                    self.idx_params, [C.NUM + idx for idx in self.idx_initials]
                )
            )
        )
        if len(difference) > 0:
            msg = 'in both search_idx and search_rgn'
            for idx in difference:
                if idx <= C.NUM:
                    raise ValueError(
                        'Set "C.{}" '.format(C.NAMES[int(idx)]) + msg
                    )
                else:
                    raise ValueError(
                        'Set "V.{}" '.format(V.NAMES[int(idx-C.NUM)]) + msg
                    )
        search_rgn = search_rgn[:, np.any(search_rgn != 0., axis=0)]

        return np.log10(search_rgn)
