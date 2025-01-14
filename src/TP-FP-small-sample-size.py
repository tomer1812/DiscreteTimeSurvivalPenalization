import pandas as pd
import numpy as np
import os
from time import time
from sklearn.model_selection import train_test_split
from pydts.examples_utils.generate_simulations_data import generate_quick_start_df
from pydts.examples_utils.plots import plot_example_pred_output
from pydts.examples_utils.plots import add_panel_text
from pydts.cross_validation import TwoStagesCV, PenaltyGridSearchCVExact
from pydts.fitters import TwoStagesFitter, DataExpansionFitter

from pydts.data_generation import EventTimesSampler
from matplotlib import pyplot as plt
import warnings
import pickle
from copy import deepcopy
from sklearn.model_selection import KFold
pd.set_option("display.max_rows", 500)
slicer = pd.IndexSlice


OUTPUT_DIR = ''


file_number = 1
runs = 100

n_cov = 35
beta1 = np.zeros(n_cov)
beta1[:5] = [1.2, 1.5, -1, -0.3, -1.2]
beta2 = np.zeros(n_cov)
beta2[:5] = [-1.2, -1, 1.4, 1, 1]

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

step = 0.1
penalizers = np.arange(-6, -3.4, step=step)
n_splits = 3

means_vector = np.zeros(n_cov)
covariance_matrix = 0.4 * np.identity(n_cov)
clip_value = 1.5

for run in range(runs):
    try:
        ets = EventTimesSampler(d_times=d_times, j_event_types=j_events)

        seed = 100 * file_number + run + 1
        print(f"Starting run {run}, seed {seed}")
        np.random.seed(seed)


        covariates = [f'Z{i + 1}' for i in range(n_cov)]

        patients_df = pd.DataFrame(data=pd.DataFrame(data=np.random.multivariate_normal(means_vector, covariance_matrix,
                                                                                        size=n_patients),
                                                     columns=covariates))
        patients_df.clip(lower=-1 * clip_value, upper=clip_value, inplace=True)
        patients_df = ets.sample_event_times(patients_df, hazard_coefs=real_coef_dict, seed=seed)
        patients_df = ets.sample_independent_lof_censoring(patients_df, prob_lof_at_t=0.01 * np.ones_like(ets.times[:-1]),
                                                           seed=seed + 1)
        patients_df = ets.update_event_or_lof(patients_df)
        patients_df.index.name = 'pid'
        patients_df = patients_df.reset_index()

        c = patients_df.groupby(['X', 'J']).size()
        if c.loc[:, [1, 2]].min() < 8:
            print(f'Skipping seed {seed} due to low number of patients at a specific (X, J)')
            continue


        penalty_cv_search = PenaltyGridSearchCVExact()
        gauc_cv_results = penalty_cv_search.cross_validate(full_df=patients_df, l1_ratio=1,
                                                           penalizers=np.exp(penalizers), n_splits=n_splits, seed=seed)

        gauc_cv_results.to_csv(os.path.join(OUTPUT_DIR, f'FP-FN_grid_search_{seed}.csv'))

        tmp_gauc = pd.concat(penalty_cv_search.global_auc).reset_index()
        tmp_gauc.columns = ['fold', 'penalizer 1', 'penalizer 2', 'gauc']
        tmp_gauc['penalizer 1'] = np.log(tmp_gauc['penalizer 1'])
        tmp_gauc['penalizer 2'] = np.log(tmp_gauc['penalizer 2'])
        tmp_gauc.to_csv(os.path.join(OUTPUT_DIR, f'global_auc_{seed}.csv'))

        chosen_eta = np.log(gauc_cv_results['Mean'].idxmax())

        pd.DataFrame(data=chosen_eta).to_csv(os.path.join(OUTPUT_DIR, f'chosen_eta_{seed}.csv'))

        nzc = pd.DataFrame()
        tp_fp_df = pd.DataFrame()
        for risk in [1, 2]:
            nonzero_count = pd.DataFrame(index=list(range(n_splits)), columns=penalizers)
            for idp, penalizer in enumerate(penalizers):
                n_TP, n_FP = [], []
                tmp_j1_params_df = pd.DataFrame()
                for i_fold in range(n_splits):
                    params_ser = penalty_cv_search.folds_grids[i_fold].meta_models[np.exp(penalizer)].beta_models[risk].params
                    # params_ser = penalty_cv_search.folds_grids[i_fold].meta_models[penalizer].beta_models[risk].params

                    nonzero_count.loc[i_fold, penalizer] = (params_ser.round(3).abs() > 0).sum()
                    tmp_j1_params_df = pd.concat([tmp_j1_params_df, params_ser], axis=1)

                    n_TP.append((params_ser.iloc[:5].round(3).abs() > 0).sum())
                    n_FP.append((params_ser.iloc[5:].round(3).abs() > 0).sum())

                tp_fp_df = pd.concat([tp_fp_df, pd.DataFrame(data=[[penalizer, risk, 'TP', np.mean(n_TP)],
                                                                   [penalizer, risk, 'FP', np.mean(n_FP)]],
                                                             columns=['penalizer', 'risk', 'type', '3-fold-n'])])

                ser_1 = tmp_j1_params_df.mean(axis=1)
                ser_1.name = penalizer

                if idp == 0:
                    j1_params_df = ser_1.to_frame()
                else:
                    j1_params_df = pd.concat([j1_params_df, ser_1], axis=1)

            nzc = pd.concat([nzc, pd.concat([nonzero_count], keys=[risk])])

        nzc.to_csv(os.path.join(OUTPUT_DIR, f'nonzero_count_{seed}.csv'))
        tp_fp_df.to_csv(os.path.join(OUTPUT_DIR, f'tp_fp_{seed}.csv'))

    except Exception as e:
       print(f"Run {run} failed: {seed}, {e}")