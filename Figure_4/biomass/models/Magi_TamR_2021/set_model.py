from .name2idx import C, V
import math

class DifferentialEquation(object):
    """two state model of tamR
    """
    def __init__(self, perturbation):
        self.perturbation = perturbation

    def diffeq(self, t, y, x):
        # Rate equation
        v = {}

        v[1] = x[C.v1_k] * (1 - x[C.v1_i] * y[V.Tam]) * y[V.S] #
        v[2] = x[C.v2_c] * y[V.Tam] * y[V.S] * 1 / (1 + math.e**-(x[C.v2_a] * (y[V.Tam_integ] - x[C.v2_b])))
        v[3] = x[C.v3_c] * y[V.P]
        v[4] = x[C.v4_k] * y[V.Tam] * y[V.S]
        v[5] = x[C.v5_k] * y[V.P]
        v[6] = x[C.v6_k] * y[V.Tam] * y[V.P] #0
        v[7] = x[C.v7_c] * y[V.Tam] * y[V.P] * 1 / (1 + math.e**-(x[C.v7_a] * (y[V.Tam_integ] - x[C.v7_b])))
        v[8] = x[C.v8_c] * y[V.R1]
        v[9] = x[C.v9_k] * y[V.R1]
        v[10] = x[C.v10_c] * y[V.Tam] * y[V.P] * 1 / (1 + math.e**-(x[C.v10_a] * (y[V.Tam_integ] - x[C.v10_b])))
        v[11] = x[C.v11_c] * y[V.R2]
        v[12] = x[C.v12_k] * y[V.R2]

        v[13] = 1 * y[V.Tam]

        if self.perturbation:
            for i, dv in self.perturbation.items():
                v[i] = v[i] * dv

        dydt = [0] * V.NUM
        dydt[V.S] = v[1] - v[2] + v[3] - v[4]
        dydt[V.P] = v[2] + v[5] - v[3] - v[6] - v[7] + v[8] - v[10] + v[11]
        dydt[V.R1] = v[7] - v[8] + v[9]
        dydt[V.R2] = v[10] - v[11] + v[12]
        dydt[V.Tam] = 0
        dydt[V.Tam_integ] = v[13]

        return dydt


def param_values():
    # Parameter values
    x = [0] * C.NUM
    x[C.v1_k] = 1.2
    x[C.v1_i] = 0.1
    x[C.v2_c] = 1.0
    x[C.v2_a] = 1.0
    x[C.v2_b] = 3.0
    x[C.v3_c] = 1.0
    x[C.v4_k] = 0.2
    x[C.v5_k] = 0.7
    x[C.v6_k] = 0.2
    x[C.v7_c] = 1.0
    x[C.v7_a] = 1.0
    x[C.v7_b] = 3.0
    x[C.v8_c] = 1.0
    x[C.v9_k] = 0.5
    x[C.v10_c] = 1.0
    x[C.v10_a] = 1.0
    x[C.v10_b] = 3.0
    x[C.v11_c] = 1.0
    x[C.v12_k] = 0.5

    return x

def initial_values():
    # Value of the initial condition
    y0 = [0] * V.NUM

    y0[V.S] = 100.0
    y0[V.P] = 0.0
    y0[V.R1] = 0.0
    y0[V.R2] = 0.0
    y0[V.Tam] = 0.0
    y0[V.Tam_integ] = 0.0

    return y0
