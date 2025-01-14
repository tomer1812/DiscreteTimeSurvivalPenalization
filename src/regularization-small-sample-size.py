import pandas as pd
import numpy as np
import os
from time import time
from sklearn.model_selection import train_test_split
from pydts.examples_utils.generate_simulations_data import generate_quick_start_df
from pydts.examples_utils.plots import plot_example_pred_output
from pydts.examples_utils.plots import add_panel_text
from pydts.cross_validation import *
from pydts.fitters import TwoStagesFitter, DataExpansionFitter
from pydts.evaluation import *
from pydts.data_generation import EventTimesSampler
from matplotlib import pyplot as plt
import warnings
import pickle
from copy import deepcopy
from sklearn.model_selection import KFold
pd.set_option("display.max_rows", 500)
slicer = pd.IndexSlice


OUTPUT_DIR = ''

n_cov = 35
beta1 = np.zeros(n_cov)
beta1[:5] = np.array([1.2, 1.5, -1, -0.3, -1.2])
beta2 = np.zeros(n_cov)
beta2[:5] = np.array([-1.2, -1, 1.4, 1, 1])


real_coef_dict = {
    "alpha": {
        1: lambda t: -4.4 + 0.3 * t,
        2: lambda t: -4.3 + 0.3 * t
    },
    "beta": {
        1: beta1,
        2: beta2
    }
}

n_patients = 500
d_times = 10
j_events = 2

ets = EventTimesSampler(d_times=d_times, j_event_types=j_events)

seed = 11
means_vector = np.zeros(n_cov)
covariance_matrix = 0.4*np.identity(n_cov)

fig_name = 'regularization_sim_small_sample_corr.png'

covariance_matrix[0, 8] = 0.1
covariance_matrix[8, 0] = 0.1
covariance_matrix[1, 9] = 0.3
covariance_matrix[9, 1] = 0.3
covariance_matrix[3, 7] = -0.3
covariance_matrix[7, 3] = -0.3
covariance_matrix[4, 11] = -0.1
covariance_matrix[11, 4] = -0.1


clip_value = 1.2

covariates = [f'Z{i + 1}' for i in range(n_cov)]

patients_df = pd.DataFrame(data=pd.DataFrame(data=np.random.multivariate_normal(means_vector, covariance_matrix,
                                                                                size=n_patients),
                                             columns=covariates))
patients_df.clip(lower= -1 * clip_value, upper=clip_value, inplace=True)
patients_df = ets.sample_event_times(patients_df, hazard_coefs=real_coef_dict, seed=seed)
patients_df = ets.sample_independent_lof_censoring(patients_df, prob_lof_at_t=0.01 * np.ones_like(ets.times[:-1]),
                                                   seed=seed + 1)
patients_df = ets.update_event_or_lof(patients_df)
patients_df.index.name = 'pid'
patients_df = patients_df.reset_index()


step = 0.25
penalizers = np.exp(np.arange(-6, -3.4, step=step))
n_splits = 3
seed = 2



penalty_cv_search = PenaltyGridSearchCVExact()
gauc_cv_results = penalty_cv_search.cross_validate(full_df=patients_df, l1_ratio=1, penalizers=penalizers,
                                                   n_splits=n_splits, seed=seed)


chosen_eta = np.log(gauc_cv_results['Mean'].idxmax())


chosen_auc_df = pd.DataFrame()
for i_fold in range(n_splits):
    mixed_two_step = penalty_cv_search.folds_grids[i_fold].get_mixed_two_stages_fitter(np.exp(chosen_eta))
    # mixed_two_step = penalty_cv_search.folds_grids[i_fold].get_mixed_two_stages_fitter(chosen_eta)
    test_df = patients_df[patients_df['pid'].isin(penalty_cv_search.test_pids[i_fold])]
    pred_df = mixed_two_step.predict_prob_events(test_df)
    auc_t = events_auc_at_t(pred_df)
    chosen_auc_df = pd.concat([chosen_auc_df, pd.concat([auc_t], keys=[i_fold])])



chosen_gauc = []
chosen_iauc1 = []
chosen_iauc2 = []
chosen_gbs = []
chosen_ibs1 = []
chosen_ibs2 = []

for i_fold in range(n_splits):
    mixed_two_step = penalty_cv_search.folds_grids[i_fold].get_mixed_two_stages_fitter(np.exp(chosen_eta))
    # mixed_two_step = penalty_cv_search.folds_grids[i_fold].get_mixed_two_stages_fitter(chosen_eta)
    test_df = patients_df[patients_df['pid'].isin(penalty_cv_search.test_pids[i_fold])]
    pred_df = mixed_two_step.predict_prob_events(test_df)
    chosen_gauc.append(global_auc(pred_df))
    chosen_gbs.append(global_brier_score(pred_df))
    iauc = events_integrated_auc(pred_df)
    ibs = events_integrated_brier_score(pred_df)
    chosen_iauc1.append(iauc[1])
    chosen_iauc2.append(iauc[2])
    chosen_ibs1.append(ibs[1])
    chosen_ibs2.append(ibs[2])



print(np.mean(chosen_gauc).round(3), np.std(chosen_gauc).round(3))
print(np.mean(chosen_iauc1).round(3), np.std(chosen_iauc1).round(3))
print(np.mean(chosen_iauc2).round(3), np.std(chosen_iauc2).round(3))
print(np.mean(chosen_gbs).round(3), np.std(chosen_gbs).round(3))
print(np.mean(chosen_ibs1).round(3), np.std(chosen_ibs1).round(3))
print(np.mean(chosen_ibs2).round(3), np.std(chosen_ibs2).round(3))


counts = patients_df.groupby(['J', 'X'])['pid'].count().unstack('J').fillna(0)

ticksize = 15
axes_title_fontsize = 17
legend_size = 13

risk_names = []
risk_colors = ['tab:blue', 'tab:green', 'tab:red']
abc_letters = ['a', 'b']
def_letters = ['c', 'd']
ghi_letters = ['e', 'f']

fig, axes = plt.subplots(3, 2, figsize=(13, 15))

for risk in [1, 2]:
    nonzero_count = pd.DataFrame(index=list(range(n_splits)), columns=penalizers)
    for idp, penalizer in enumerate(penalizers):

        tmp_j1_params_df = pd.DataFrame()
        for i_fold in range(n_splits):
            # params_ser = penalty_cv_search.folds_grids[i_fold].meta_models[np.exp(penalizer)].beta_models[risk].params
            params_ser = penalty_cv_search.folds_grids[i_fold].meta_models[penalizer].beta_models[risk].params

            nonzero_count.loc[i_fold, penalizer] = (params_ser.round(3).abs() > 0).sum()
            tmp_j1_params_df = pd.concat([tmp_j1_params_df, params_ser], axis=1)

        ser_1 = tmp_j1_params_df.mean(axis=1)
        ser_1.name = penalizer

        if idp == 0:
            j1_params_df = ser_1.to_frame()
        else:
            j1_params_df = pd.concat([j1_params_df, ser_1], axis=1)

    ax = axes[0, risk - 1]
    add_panel_text(ax, abc_letters[risk - 1])
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.tick_params(axis='both', which='minor', labelsize=ticksize)
    ax.set_xlabel(fr'Log ($\eta_{risk}$)', fontsize=axes_title_fontsize)
    ax.set_ylabel(f'Number of Non-Zero Coefficients', fontsize=axes_title_fontsize)
    ax.set_title(rf'$\beta_{risk}$', fontsize=axes_title_fontsize)
    ax.axhline(5, ls='--', color='k', alpha=0.3, label='True Value')
    ax.axvline(chosen_eta[risk - 1], color=risk_colors[risk - 1], alpha=1, ls='--', lw=1,
               label=rf'Chosen $Log (\eta_{risk})$')
    ax.set_ylim([0, 50])

    for idp, penalizer in enumerate(penalizers):

        count = nonzero_count[penalizer].mean()
        if idp == 0:
            ax.scatter(np.log(penalizer), count, color=risk_colors[risk - 1], alpha=0.8, marker='P',
                       label=f'{n_splits}-Fold mean')
        else:
            ax.scatter(np.log(penalizer), count, color=risk_colors[risk - 1], alpha=0.8, marker='P')
        if penalizer == chosen_eta[risk - 1]:
            print(f"Risk {risk}: {count} non-zero coefficients at chosen eta {chosen_eta[risk - 1]}")

    ax.legend(fontsize=legend_size)

    ax = axes[1, risk - 1]
    add_panel_text(ax, def_letters[risk - 1])
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.tick_params(axis='both', which='minor', labelsize=ticksize)
    for i in range(len(j1_params_df)):
        ax.plot(np.log(penalizers), j1_params_df.iloc[i].values, lw=1)

        if i == 0:
            ax.set_ylabel(f'{n_splits}-Fold Mean Coefficient Value', fontsize=axes_title_fontsize)
            ax.set_xlabel(fr'Log ($\eta_{risk}$)', fontsize=axes_title_fontsize)
            ax.set_title(rf'$\beta_{risk}$', fontsize=axes_title_fontsize)
            ax.axvline(chosen_eta[risk - 1], color=risk_colors[risk - 1], alpha=1, ls='--', lw=1)

    ax = axes[2, risk - 1]
    add_panel_text(ax, ghi_letters[risk - 1])
    ax.tick_params(axis='both', which='major', labelsize=ticksize)
    ax.tick_params(axis='both', which='minor', labelsize=ticksize)
    mean_auc = chosen_auc_df.loc[slicer[:, risk], :].mean(axis=0)
    std_auc = chosen_auc_df.loc[slicer[:, risk], :].std(axis=0)
    ax.errorbar(mean_auc.index, mean_auc.values, yerr=std_auc.values, fmt="o", color=risk_colors[risk - 1], alpha=0.8)
    ax.set_yticks(np.arange(0, 1.1, 0.1))
    ax.set_yticklabels([c.round(1) for c in np.arange(0, 1.1, 0.1)])
    ax.set_xlabel(r'Time', fontsize=axes_title_fontsize)
    ax.set_ylabel(fr'AUC$_{risk} (t)$', fontsize=axes_title_fontsize)
    ax.set_title(fr'Log ($\eta_{risk}$) = {chosen_eta[risk - 1]:.2f}', fontsize=axes_title_fontsize)
    ax.set_ylim([0, 1])
    l = 16
    ax.set_xticks(range(l))
    ax.set_xticklabels([f'{x}' for x in range(l)])

    ax.axhline(0.5, ls='--', color='k', alpha=0.3)
    ax2 = ax.twinx()
    ax2.bar(counts.index, counts[risk].values.squeeze(), color=risk_colors[risk - 1], alpha=0.8, width=0.4)
    ax2.set_ylabel('Number of observed events', fontsize=axes_title_fontsize, color=risk_colors[risk - 1])
    ax2.tick_params(axis='y', colors=risk_colors[risk - 1])
    ax2.set_ylim([0, 250])
    ax2.tick_params(axis='both', which='major', labelsize=ticksize)
    ax2.tick_params(axis='both', which='minor', labelsize=ticksize)

fig.tight_layout()

fig.savefig(os.path.join(OUTPUT_DIR, fig_name), dpi=800)



start = time()
cross_validator_null = TwoStagesCVExact()
cross_validator_null.cross_validate(full_df=patients_df, n_splits=n_splits, seed=seed, nb_workers=1)
end = time()
print(f"Finished {int(end-start)} seconds")


print(np.mean(list(cross_validator_null.global_auc.values())).round(3))

print(np.std(list(cross_validator_null.global_auc.values())).round(3))
