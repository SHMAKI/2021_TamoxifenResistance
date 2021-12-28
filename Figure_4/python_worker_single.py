import sys
from importlib import import_module
from biomass import optimize_continue
from biomass import optimize

model_module = import_module("biomass.models." + sys.argv[1])

optimize(model_module, sys.argv[2], max_generation=30000)

# # when the calc. was aborted, continue the simulation with the following code:
# optimize_continue(model_module, sys.argv[2])
