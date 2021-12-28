import sys
from importlib import import_module
from biomass import optimize_continue
from biomass import optimize

model_module = import_module("biomass.models." + sys.argv[1])
#optimize_continue(model_module, sys.argv[2])
optimize(model_module, sys.argv[2], max_generation=10000)

#from biomass import run_simulation
#run_simulation(model_module, viz_type='average', show_all=False, stdev=True)
#run_simulation(model_module, viz_type='best', show_all=True, stdev=False)
#run_simulation(model_module, viz_type='original', show_all=True, stdev=True)

#from biomass import run_analysis
#run_analysis(model_module, target='parameter', metric='integral_after_3w', style='barplot')
#run_analysis(model_module, target='parameter', metric='max_after_3w', style='barplot')
#run_analysis(model_module, target='parameter', metric='last', style='barplot')
