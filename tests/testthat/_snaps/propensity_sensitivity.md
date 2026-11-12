# 8. print, summary and as.data.frame

    Code
      print(rounded_run)
    Output
      Propensity model sensitivity
      1. Quality gate: 2 of 5 models kept (1 algorithm: random forest); log-loss
         within 5 % of the best; relative balance rule: max |SMD| main terms at most
         0.10 above the best model's and at most 0.25 (limit here 0.25).
         The kept models are random seeds of one algorithm (random forest): rank
            disagreement and pair instability here measure seed noise within that
            algorithm, not disagreement between algorithms.
         Excluded:
           logistic_regression_1: log-loss above (1 + quality_tolerance) x best; max
             |SMD| main terms more than balance_margin above the best model's or
             above balance_cap
           gradient_boosting_1: log-loss above (1 + quality_tolerance) x best; max
             |SMD| main terms more than balance_margin above the best model's or
             above balance_cap
           gradient_boosting_2: log-loss above (1 + quality_tolerance) x best; max
             |SMD| main terms more than balance_margin above the best model's or
             above balance_cap
         Treated patients kept by each kept model: random_forest_1 100 %,
           random_forest_2 100 %
      2. Rank disagreement (primary; study mean over treated): 0.00804
         Equal weight per algorithm (exploratory): NA
         Between/within-algorithm ratio (exploratory; median over treated): NA.
           Disagreement between algorithms divided by disagreement between seeds of
           the same algorithm; with three seeds per algorithm, seed noise alone
           gives about 0.4.
         Pair instability (exploratory; study mean over treated): 0.368 of a
           possible 0.5 with 2 models. It is high whenever many untreated patients
           have similar scores, even when the models agree.
         Bootstrap rank instability of the logistic model (descriptive): 0.0361
      3. Balance spread (descriptive; range of max |SMD| main terms over the kept
         models): 0.0183 (0.16 to 0.178)
      Control spread (descriptive; mean scaled over treated): 1.02; 1 = the partners
         chosen by different models differ as much as a typical matched pair (scale:
         the matches of logistic_regression_1)
      Your own model (user): kept by the gate (its balance is within the limit);
         mean absolute rank difference from the kept models: 0.101
      Matching: 1:1 nearest neighbour on the logit, with replacement, order
         "largest", caliper none; target: ATT.

---

    Code
      print(summary(rounded_run))
    Output
      Propensity model sensitivity: summary
      1. Quality gate: 2 of 5 models kept (1 algorithm: random forest); log-loss
         within 5 % of the best; relative balance rule: max |SMD| main terms at most
         0.10 above the best model's and at most 0.25 (limit here 0.25).
         The kept models are random seeds of one algorithm (random forest): rank
            disagreement and pair instability here measure seed noise within that
            algorithm, not disagreement between algorithms.
         Excluded:
           logistic_regression_1: log-loss above (1 + quality_tolerance) x best; max
             |SMD| main terms more than balance_margin above the best model's or
             above balance_cap
           gradient_boosting_1: log-loss above (1 + quality_tolerance) x best; max
             |SMD| main terms more than balance_margin above the best model's or
             above balance_cap
           gradient_boosting_2: log-loss above (1 + quality_tolerance) x best; max
             |SMD| main terms more than balance_margin above the best model's or
             above balance_cap
         Treated patients kept by each kept model: random_forest_1 100 %,
           random_forest_2 100 %
      2. Rank disagreement (primary; study mean over treated): 0.00804
         Equal weight per algorithm (exploratory): NA
         Between/within-algorithm ratio (exploratory; median over treated): NA.
           Disagreement between algorithms divided by disagreement between seeds of
           the same algorithm; with three seeds per algorithm, seed noise alone
           gives about 0.4.
         Pair instability (exploratory; study mean over treated): 0.368 of a
           possible 0.5 with 2 models. It is high whenever many untreated patients
           have similar scores, even when the models agree.
         Bootstrap rank instability of the logistic model (descriptive): 0.0361
      3. Balance spread (descriptive; range of max |SMD| main terms over the kept
         models): 0.0183 (0.16 to 0.178)
      Control spread (descriptive; mean scaled over treated): 1.02; 1 = the partners
         chosen by different models differ as much as a typical matched pair (scale:
         the matches of logistic_regression_1)
      Your own model (user): kept by the gate (its balance is within the limit);
         mean absolute rank difference from the kept models: 0.101
      Matching: 1:1 nearest neighbour on the logit, with replacement, order
         "largest", caliper none; target: ATT.
      
      Study-level numbers (kept models / every model):
       quantity                                        label       kept_models
       number_of_models_used                           count       2          
       mean_rank_disagreement_treated                  primary     0.00804    
       mean_rank_disagreement_equal_weight_treated     exploratory NA         
       mean_rank_disagreement_all                      descriptive 0.0173     
       median_between_within_ratio_treated             exploratory NA         
       mean_logit_disagreement_treated                 descriptive 0.117      
       mean_control_spread_scaled_treated              descriptive 1.02       
       p90_control_spread_scaled_treated               descriptive 2.04       
       smd_max_main_range                              descriptive 0.0183     
       smd_max_main_sd                                 descriptive 0.013      
       smd_max_extended_range                          descriptive 0.0168     
       smd_max_extended_sd                             descriptive 0.0119     
       energy_distance_range                           descriptive 0.00792    
       energy_distance_sd                              descriptive 0.0056     
       ks_max_range                                    descriptive 0.00541    
       ks_max_sd                                       descriptive 0.00382    
       pair_recurrence_models_mean_treated             descriptive 0.632      
       pair_instability_models_mean_treated            exploratory 0.368      
       pair_recurrence_bootstrap_mean_treated          descriptive 0.62       
       treated_kept_share                              descriptive 1          
       treated_kept_share_min                          descriptive 1          
       bootstrap_rank_instability_treated              descriptive 0.0361     
       reference_mean_absolute_rank_difference_treated descriptive 0.101      
       every_model
       5          
       0.058      
       0.0686     
       0.0653     
       6.82       
       0.674      
       1.24       
       1.74       
       0.21       
       0.0847     
       0.27       
       0.123      
       0.108      
       0.0427     
       0.222      
       0.0928     
       0.341      
       0.659      
       0.62       
       1          
       1          
       0.0361     
       0.0857     
      
      Models (in_use = kept by the quality gate):
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
      
      Gate notes:
        logistic_regression_1: log-loss above (1 + quality_tolerance) x best; max
            |SMD| main terms more than balance_margin above the best model's or
            above balance_cap
        random_forest_1: passed
        random_forest_2: passed
        gradient_boosting_1: log-loss above (1 + quality_tolerance) x best; max
            |SMD| main terms more than balance_margin above the best model's or
            above balance_cap
        gradient_boosting_2: log-loss above (1 + quality_tolerance) x best; max
            |SMD| main terms more than balance_margin above the best model's or
            above balance_cap
        reference: your own model: an in-sample score from the matchit object, so
            its log-loss is not comparable and not computed; judged on balance only:
            kept by the gate (its balance is within the limit)
      
      Quality gate: log-loss within 5 % of the best; relative balance rule: max
         |SMD| main terms at most 0.10 above the best model's and at most 0.25
         (limit here 0.25).
      Matching: 1:1 nearest neighbour on the logit, with replacement, order
         "largest", caliper none; target: ATT.

