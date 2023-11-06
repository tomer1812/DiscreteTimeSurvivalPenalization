import pandas as pd
import os
from config import *
import matplotlib.pyplot as plt
import numpy as np
from pydts.data_generation import EventTimesSampler


final_dfs = []
n_patients_list = [250, 500]
for n_patients in n_patients_list:
    case = f'n{n_patients}_beta_results_comparison.csv'

    df = pd.read_csv(os.path.join(OUTPUT_DIR, case))
    df.index = [r'$\beta_{11}$', r'$\beta_{12}$', r'$\beta_{13}$', r'$\beta_{14}$', r'$\beta_{15}$',
                 r'$\beta_{21}$', r'$\beta_{22}$', r'$\beta_{23}$', r'$\beta_{24}$', r'$\beta_{25}$']

    df.columns = pd.MultiIndex.from_tuples([('', 'True'),
                                            ('Lee et al.', 'Estimate'),
                                            ('Lee et al.', 'Estimated SE'),
                                            ('two-step', 'Estimate'),
                                            ('two-step', 'Estimated SE'),
                                            ('two-step', 'Empirical SE'),
                                            ('two-step', 'Coverage Rate')])

    final_dfs.append(df.astype(float).round(3))

final_df = pd.concat(final_dfs, keys=n_patients_list)
print(final_df.to_latex(escape=False))


filename = 'alpha_small_dataset.png'

d_times = 7
n_cov = 5
j_events = 2


real_coef_dict = {
    "alpha": {
        1: lambda t: -1.4 + 0.4 * np.log(t),
        2: lambda t: -1.3 + 0.4 * np.log(t)
    },
    "beta": {
        1: -0.7*np.log([0.8, 3, 3, 2.5, 2]),
        2: -0.6*np.log([1, 3, 4, 3, 2])
    }
}

ets = EventTimesSampler(d_times=d_times, j_event_types=j_events)
covariates = [f'Z{i + 1}' for i in range(n_cov)]

first_model_name = 'Lee et al.'
second_model_name = 'two-step'
times = range(1, d_times + 1)

lee_colors = ['tab:blue', 'tab:green']
two_step_colors = ['navy', 'darkgreen']
true_colors = ['tab:blue', 'tab:green']

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for idn, n_patients in enumerate(n_patients_list):
    case = f'n{n_patients}_results_comparison.csv'

    np.random.seed(idn)
    patients_df = pd.DataFrame(data=pd.DataFrame(data=np.random.uniform(0, 1, size=[n_patients, n_cov]),
                                                 columns=covariates))

    patients_df = ets.sample_event_times(patients_df, hazard_coefs=real_coef_dict, seed=idn)
    patients_df = ets.sample_independent_lof_censoring(patients_df,
                                                       prob_lof_at_t=0.02 * np.ones(d_times))
    patients_df = ets.update_event_or_lof(patients_df)
    patients_df.index.name = 'pid'
    patients_df = patients_df.reset_index()
    counts = patients_df.groupby(['J', 'X'])['pid'].count().unstack('J').fillna(0)

    alpha_res_df = pd.read_csv(os.path.join(OUTPUT_DIR, case))

    ax = axes[int(idn % 2)]
    ax.set_title(f'n={n_patients}', fontsize=15)
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.tick_params(axis='both', which='minor', labelsize=15)

    for j in [1, 2]:
        tmp_lee_alpha = alpha_res_df.iloc[((d_times+n_cov)*(j-1)):((d_times+n_cov)*(j-1)+d_times), 1]
        ax.scatter(times, tmp_lee_alpha.values,
                   label=f'J={j} ({first_model_name})', color=lee_colors[j - 1], marker='o', alpha=0.4, s=40)

        tmp_two_step_alpha = alpha_res_df.iloc[((d_times+n_cov)*(j-1)):((d_times+n_cov)*(j-1)+d_times), 3]
        ax.scatter(times, tmp_two_step_alpha.values,
                   label=f'J={j} ({second_model_name})', color=two_step_colors[j - 1], marker='*', alpha=0.7, s=20)

        true_values = [real_coef_dict['alpha'][j](t) for t in times]
        ax.plot(times, true_values, label=f'J={j} (True)', ls='--', color=true_colors[j - 1])

        ax.set_xlabel(r'Time', fontsize=18)
        ax.set_ylabel(r'$\alpha_{jt}$', fontsize=18)
        ax.legend(loc='upper right', fontsize=12)
        ax.set_ylim([-2.2, 1])
        ax.set_xlim([0, d_times+1])

    ax2 = ax.twinx()
    ax2.bar(counts.index, counts[1].values.squeeze(), label='J=1', color='tab:red', alpha=0.4, width=0.5)
    ax2.bar(counts.index, counts[2].values.squeeze(), label='J=2', color='tab:brown', alpha=0.6, align='edge',
            width=0.5)
    ax2.legend(loc='upper left', fontsize=12)
    ax2.set_ylabel('Number of observed events', fontsize=16, color='red')
    ax2.tick_params(axis='y', colors='red')
    ax2.set_ylim([0, 150])
    ax2.tick_params(axis='both', which='major', labelsize=15)
    ax2.tick_params(axis='both', which='minor', labelsize=15)
fig.tight_layout()

if filename is not None:
    fig.savefig(os.path.join(OUTPUT_DIR, filename), dpi=300)