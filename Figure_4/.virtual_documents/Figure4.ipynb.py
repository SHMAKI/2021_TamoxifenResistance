import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'serif','serif':['Helvetica']})
plt.rcParams['pdf.fonttype'] = 42
from biomass.models import Magi_TamR_2021 as model
from biomass import run_simulation
from biomass import run_analysis

from statannot import add_stat_annotation
import itertools


# draw simulation graph
run_simulation(model, viz_type='average', show_all=False, stdev=True) #Figure 4b,c


# sensitivity analysis
run_analysis(model,target="reaction",metric="integral_after_3w")


from biomass.exec_model import ExecModel
from biomass.analysis import ReactionSensitivity

ExecModel_obj = ExecModel(model)

biological_processes = model.ReactionNetwork().group()
reaction_indices = np.sum(biological_processes, axis=0)
reaction_indices


# load sensicitivity data
sensitivity_coefficients = np.load("./biomass/models/Magi_TamR_2021/sensitivity_coefficients/reaction/integral_after_3w/sc.npy")


# check where is "growth_rate"
# print(ExecModel_obj.obs)
k = ExecModel_obj.obs.index("growth_rate")
# sensitivity_array: sensicitivity data of growth_rate,
# param_set x reaction x condition(ctrl, tam)
sensitivity_array = sensitivity_coefficients[:, :, k, :]


# Remove NaN
nan_idx = []
for i in range(sensitivity_array.shape[0]):
    for j in range(sensitivity_array.shape[1]):
        if any(np.isnan(sensitivity_array[i, j, :])):
            nan_idx.append(i)

sensitivity_array = np.delete(
    sensitivity_array, nan_idx, axis=0
)


# sensitivity_array to tidy pd.DataFrame
sa_df_1 = pd.DataFrame(sensitivity_array[:,:,0]) #ctrl
sa_df_2 = pd.DataFrame(sensitivity_array[:,:,1]) #tam
sa_df_1.columns = ["v" + sub for sub in [str(n) for n in reaction_indices]]
sa_df_2.columns = ["v" + sub for sub in [str(n) for n in reaction_indices]]
sa_df_1["Condition"] = "Control"
sa_df_2["Condition"] = "Tamoxifen"

sa_df = pd.concat([sa_df_1, sa_df_2])

#create individual id, and convert df to long
res=[n + str(s) for (n, s) in zip(sa_df["Condition"].to_list(), sa_df.index.to_list())]
sa_df["id"] = res

sa_df_wide = pd.wide_to_long(sa_df, stubnames='v', i=["id"], j="reaction")
sa_df_wide = sa_df_wide.reset_index()
sa_df_wide["reaction"] = "v" + sa_df_wide["reaction"].astype(str)
#sa_df_wide


# def draw func from biomass
def draw_vertical_span(biological_processes, width):
    if len(biological_processes) > 1:
        left_end = 0
        for i, proc in enumerate(biological_processes):
            if i % 2 == 0:
                plt.axvspan(
                    left_end - width-0.2,
                    left_end - width + len(proc)-0.2,
                    facecolor='k', alpha=0.05
                )
            left_end += len(proc)


# params for drawing graph
options = ExecModel_obj.viz.sensitivity_options

# savefig
fig = plt.figure(figsize=(8,4))

ax = sns.boxplot(x="reaction", y="v", data=sa_df_wide, hue="Condition", width=0.6)
ax.set_ylabel("Local sensitivity to mean growth rate\nfrom W3 to W10")
ax.legend(loc=options['legend_loc'], frameon=False)

add_stat_annotation(ax, data=sa_df_wide, x="reaction", y="v", hue="Condition",
                box_pairs=[(("v9", "Tamoxifen"),("v12", "Tamoxifen")),],
                    test='Wilcoxon', text_format='simple', loc='inside', verbose=2)
                    
draw_vertical_span(biological_processes, options['width'])
plt.savefig("Fig/fig4d.pdf")


# calc dual_inhibition_scores, takes long time because of nested for loops
dual_inhibition_scores = ReactionSensitivity._calc_dual_inhibition_scores(ExecModel_obj, "integral_after_3w", reaction_indices)


# np.save("processed_data/dual_inhibition_scores.npy", dual_inhibition_scores,)


# draw figure

# specify cell transition and growth reaction
x_vec = [7,9]
y_vec = [10,12]

inhibition_rate = np.arange(0,11,1)/10

for x, y in itertools.product(x_vec, y_vec):

    mat = dual_inhibition_scores[:, # a set of param: 20, 
                                 np.where(np.array(reaction_indices) == y)[0][0], # v_x: 12,
                                 np.where(np.array(reaction_indices) == x)[0][0], # v_y: 12,
                                 0, # target (growth_rate, sensitive, preR, R1, R2)
                                 1, # condition (con, tam)
                                 :, # inhibition range_x 0-100get_ipython().run_line_magic(":", " 11,")
                                 :, # inhibition range_y 0-100get_ipython().run_line_magic(":", " 11")
                                ].mean(axis=0) / 70 # (average of 20 param sets) -> average of 70 time points (3-10 week, 0.1 week increments)
    # heatmap
    sns.heatmap(mat, cmap="jet", vmin=0.5, vmax=1.2, annot=True, fmt=".2f",
                xticklabels=inhibition_rate,
               yticklabels=inhibition_rate,)
    plt.xlabel("v" + str(x))
    plt.ylabel("v" + str(y))
    plt.savefig("Fig/fig4e_v" + str(x) + "_v" + str(y) + ".pdf")
    plt.close()
