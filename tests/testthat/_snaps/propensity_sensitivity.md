# 8. print, summary and as.data.frame

    Code
      print(rounded_run)
    Output
      Propensity model sensitivity
      1. Quality gate: 2 of 5 models kept (1 algorithm: random forest); log-loss within 5 % of the best; relative balance rule: max |SMD| main terms at most 0.1 above the best model's and <= 0.25, limit here 0.25
         Excluded:
           logistic_regression_1: log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
           gradient_boosting_1: log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
           gradient_boosting_2: log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
      2. Rank disagreement (study mean over treated): 0.00804; median between/within ratio NA (exploratory; near 1 = seed noise)
         Pair instability across the models in use (study mean over treated; exploratory): 0.368
         Bootstrap rank instability of the logistic model (study mean over treated): 0.0361
      3. Balance spread (range of max |SMD| main terms across the models in use): 0.0183 (0.16 to 0.178)
      Mean scaled control spread over treated (descriptive): 1.02 (1 = controls differ as much as a typical match differs from its treated patient)
      Treated kept: 100 % (lowest model in use)
      Reference model (user): passes the balance gate; mean absolute rank difference from the models in use: 0.101

---

    Code
      print(summary(rounded_run))
    Output
      1. Quality gate: 2 of 5 models kept (1 algorithm: random forest); log-loss within 5 % of the best; relative balance rule: max |SMD| main terms at most 0.1 above the best model's and <= 0.25, limit here 0.25
         Excluded:
           logistic_regression_1: log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
           gradient_boosting_1: log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
           gradient_boosting_2: log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
      2. Rank disagreement (study mean over treated): 0.00804; median between/within ratio NA (exploratory; near 1 = seed noise)
         Pair instability across the models in use (study mean over treated; exploratory): 0.368
         Bootstrap rank instability of the logistic model (study mean over treated): 0.0361
      3. Balance spread (range of max |SMD| main terms across the models in use): 0.0183 (0.16 to 0.178)
      
      Study-level numbers (with the quality gate / without it):
                                              quantity with_gate without_gate
                                 number_of_models_used   2.00000       5.0000
                        mean_rank_disagreement_treated   0.00804       0.0580
                            mean_rank_disagreement_all   0.01730       0.0653
                   median_between_within_ratio_treated        NA       6.8200
                       mean_logit_disagreement_treated   0.11700       0.6740
                    mean_control_spread_scaled_treated   1.02000       1.2400
                     p90_control_spread_scaled_treated   2.04000       1.7400
                                    smd_max_main_range   0.01830       0.2100
                                       smd_max_main_sd   0.01300       0.0847
                                smd_max_extended_range   0.01680       0.2700
                                   smd_max_extended_sd   0.01190       0.1230
                                 energy_distance_range   0.00792       0.1080
                                    energy_distance_sd   0.00560       0.0427
                                          ks_max_range   0.00541       0.2220
                                             ks_max_sd   0.00382       0.0928
                   pair_recurrence_models_mean_treated   0.63200       0.3410
                  pair_instability_models_mean_treated   0.36800       0.6590
                pair_recurrence_bootstrap_mean_treated   0.62000       0.6200
                                    treated_kept_share   1.00000       1.0000
                                treated_kept_share_min   1.00000       1.0000
                    bootstrap_rank_instability_treated   0.03610       0.0361
       reference_mean_absolute_rank_difference_treated   0.10100       0.0857
      
      Models:
                    model_id           algorithm logloss smd_max_main
       logistic_regression_1 logistic_regression   0.412        0.271
             random_forest_1       random_forest   0.363        0.160
             random_forest_2       random_forest   0.360        0.178
         gradient_boosting_1   gradient_boosting   0.378        0.370
         gradient_boosting_2   gradient_boosting   0.379        0.273
                   reference    reference (user)      NA        0.239
       smd_max_extended energy_distance ks_max treated_kept passed_quality_gate
                  0.459          0.0953  0.357            1               FALSE
                  0.244          0.0555  0.141            1                TRUE
                  0.227          0.0634  0.135            1                TRUE
                  0.498          0.1640  0.216            1               FALSE
                  0.333          0.0950  0.157            1               FALSE
                  0.239          0.0840  0.341            1                TRUE
       in_use
        FALSE
         TRUE
         TRUE
        FALSE
        FALSE
        FALSE
                                                                                                                                      gate_note
       log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
                                                                                                                                         passed
                                                                                                                                         passed
       log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
       log-loss above (1 + quality_tolerance) x best; max |SMD| main terms more than balance_margin above the best model's or above balance_cap
                            in-sample score from the user's object, not comparable (log-loss not computed); gate judged on balance only: passed
      
      Quality gate: log-loss within 5 % of the best; relative balance rule: max |SMD| main terms at most 0.1 above the best model's and <= 0.25, limit here 0.25
      Matching: 1:1 nearest neighbour on the logit, ATT, with replacement, order "largest", caliper none

