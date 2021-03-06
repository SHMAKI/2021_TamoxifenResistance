import matplotlib #matplotlib.colors.TABLEAU_COLORS
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from .observable import observables
from .observable import NumericalSimulation #追加

class Visualization(object):
    """ Plotting parameters for customizing figure properties

    Attributes
    ----------
    timecourse_options : list of dict
        Plotting options for figure/simulation/<viz_type>/<each_observable>.

    multiplot_options : dict
        Plotting options for figure/simulation/<viz_type>/multiplot_observables.

    sensitivity_options : dict
        Plotting options for figure/sensitivity.

    """
    #def __init__(self):
    def __init__(self):
        self.cm = plt.cm.get_cmap('tab20')
        sim_vis = NumericalSimulation() #追加
        cm_range = [range(0, 20, 2), range(20)][len(sim_vis.conditions) > 10] #追加

        self.timecourse_options = [
            {
                'divided_by' : 1,  # to convert time unit. (e.g. sec -> min)
                'xlim' : (),
                'xticks' : None,
                'xlabel': 'Time',
                'ylim' : (),
                'yticks' : None,
                'ylabel': observables[i].replace('__', '\n').replace('_', ' '),
                'exp_data' : True,  # if False, experimental data will not be shown
                'cmap' : [self.cm.colors[j] for j in cm_range], #range(20)
                'shape' : Line2D.filled_markers,
                'dont_show' : [],  # conditions you don't want to plot
            } for i, _ in enumerate(observables)]

        self.multiplot_options = {
            'fig_name' : 'multiplot_observables',
            'observables' : [],
            'condition' : None,
            'xlim' : (),
            'xticks' : None,
            'xlabel': 'Time',
            'ylim' : (),
            'yticks' : None,
            'ylabel': '',
            'cmap' : [self.cm.colors[j] for j in cm_range],
            'shape' : Line2D.filled_markers,
        }

        self.sensitivity_options = {
            'figsize' : (12, 5),
            'width' : 0.3,
            'legend_loc' : 'best', #'upper left',
            'cmap' : [self.cm.colors[j] for j in cm_range],
        }

    def get_timecourse_options(self):
        for i, _ in enumerate(observables):
            #self.timecourse_options[i]['divided_by'] = 60
            #self.timecourse_options[i]['xlim'] = (0, 150)
            #self.timecourse_options[i]['xticks'] = [20*i for i in range(8)]
            self.timecourse_options[i]['xlabel'] = 'TIME (week)'
            #self.timecourse_options[i]['ylim'] = (0, 310)
            #self.timecourse_options[i]['yticks'] = [50*i for i in range(7)]

            self.timecourse_options[
                observables.index('growth_rate')
            ]['ylabel'] = 'growth_rate'
            #self.timecourse_options[
            #    observables.index('unphosphorylated_MAPK')
            #]['ylabel'] = 'ERK'

        return self.timecourse_options

    def multiplot_observables(self):
        self.multiplot_options['fig_name'] = 'subpopulations'
        self.multiplot_options['observables'] = [
        "Sensitive_Cells",
        "Pre-resistant_Cells",
        'Resistant_Cells1',
        'Resistant_Cells2',
        ]
        self.multiplot_options['condition'] = 'Tamoxifen'
        #self.multiplot_options['xlim'] = (0, 150)
        #self.multiplot_options['xticks'] = [20*i for i in range(8)]
        self.multiplot_options['xlabel'] = 'TIME (week)'
        #self.multiplot_options['ylim'] = (0, 310)
        #self.multiplot_options['yticks'] = [50*i for i in range(7)]
        self.multiplot_options['ylabel'] = 'rate of subpopulations'

        return self.multiplot_options

    @staticmethod
    def set_timecourse_rcParams():
        """ figure/simulation
        """
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.linewidth'] = 1.5
        plt.rcParams['xtick.major.width'] = 1.5
        plt.rcParams['ytick.major.width'] = 1.5
        plt.rcParams['lines.linewidth'] = 1.8
        plt.rcParams['lines.markersize'] = 12
        # plt.rcParams['font.family'] = 'Arial'
        # plt.rcParams['mathtext.fontset'] = 'custom'
        # plt.rcParams['mathtext.it'] = 'Arial:italic'

    @staticmethod
    def set_param_range_rcParams():
        """ figure/param_range
        """
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['xtick.major.width'] = 1.2
        plt.rcParams['ytick.major.width'] = 1.2
        # plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def set_sensitivity_rcParams():
        """ figure/sensitivity
        """
        plt.rcParams['font.size'] = 12
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['xtick.major.width'] = 1.2
        plt.rcParams['ytick.major.width'] = 1.2
        # plt.rcParams['font.family'] = 'Arial'

    @staticmethod
    def convert_species_name(name):
        """ figure/sensitivity/initial_condition
        - Sensitivity for species with nonzero initial conditions
        """
        '''
        if name == 'MKKK':
            return 'Raf'
        elif name == 'MKK':
            return 'MEK'
        elif name == 'MAPK':
            return 'ERK'
        '''
        return name
