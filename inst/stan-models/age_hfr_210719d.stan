functions{
  real dirichlet_multinomial_lpmf(int[] y, int n, vector alpha)
  {

    real l_gamma_alpha = lgamma(sum(alpha));
    real l_fact_n = lgamma(n + 1);
    real l_famm_alpha_n = lgamma(n + sum(alpha));
    real sum_x_alpha = sum( lgamma( to_vector(y) + alpha) - (lgamma(to_vector(y) + 1) + lgamma(alpha)) );

    return l_gamma_alpha + l_fact_n - l_famm_alpha_n + sum_x_alpha;
  }
  
  int[] make_weeks_to_days(int T)
  {
    int idx[T,7];
    for(s in 1:T)
    {
        idx[s,:] = rep_array(s,7);
    }
    return( to_array_1d(idx) );
  }
  
  vector make_logistic_prop_P1(
    // parameters
    real propp1_logistic_growthrate,
    real propp1_logistic_midpoint,
    // data
    real propp1_logistic_max_val,
    vector first_day_of_week_idx
    )
  {
    return( propp1_logistic_max_val * inv(1.0 + exp(-propp1_logistic_growthrate * (first_day_of_week_idx - propp1_logistic_midpoint))) );
  }
  
  matrix make_time_varying_fatality_rate_of_variant(
    // parameters
    vector logit_base_fatality_rate,
    vector fr_regcoeff_overall,
    // data
    int A,
    int Th,
    matrix[] fr_predictors
    )
  {
    matrix[Th,A] logit_fatality_rate_of_variant;
    for(a in 1:A)
    {
        logit_fatality_rate_of_variant[:,a] = fr_predictors[a] * fr_regcoeff_overall;
    }
    logit_fatality_rate_of_variant += rep_matrix(logit_base_fatality_rate',Th);
    return(inv_logit(logit_fatality_rate_of_variant));
  }
  
  matrix make_fr_multiplier(
    // parameters
    vector fr_regcoeff_overall,
    // data
    int A,
    int Th,
    matrix[] fr_predictors
    )
  {
    matrix[Th,A] fr_multiplier_by_week;
    for(a in 1:A)
    {
        fr_multiplier_by_week[:,a] = exp( fr_predictors[a] * fr_regcoeff_overall );
    }
    return(fr_multiplier_by_week);
  }
    
  matrix calculate_age_prop_of_hosps_by_variant(
    // parameters
    vector log_age_rate,
    // data
    int Th,
    int A,
    matrix log_age_pop_mat
    )
  {
    matrix[Th,A] prop_age;
    prop_age = rep_matrix( log_age_rate', Th ) + log_age_pop_mat;
    for(t in 1:Th)
    {
        prop_age[t,] = softmax(prop_age[t,]')';
    }
    return(prop_age);
  }

  matrix calculate_joint_prob_of_deaths_and_hospital_admission_with_variant_conditional_on_age(
    // parameters
    vector prop_1,
    matrix prop_age_hosps_of_1,
    matrix prop_age_hosps_of_2,
    matrix hosp_fatality_rate_1,
    // data
    int Th,
    int A
    )
  {
    matrix[Th,A] ans_1;
    matrix[Th,A] ans_2;
    ans_1 = rep_matrix( prop_1, A) .* prop_age_hosps_of_1;
    ans_2 = rep_matrix( 1.-prop_1, A) .* prop_age_hosps_of_2;
    ans_1 = ans_1 ./ (ans_1+ans_2);
    ans_1 .*= hosp_fatality_rate_1;
    return(ans_1);
  }
  
  matrix calculate_age_prop_of_deaths_in_hosps_by_variant(
    // parameters
    matrix prop_age_hosps_of_variant,
    matrix hosp_fatality_rate_of_variant_by_week,
    // data
    int Th,
    int A
    )
  {
    matrix[Th,A] ans;
    vector[Th] denominator;
    ans = prop_age_hosps_of_variant .* hosp_fatality_rate_of_variant_by_week;
    denominator = rows_dot_product(prop_age_hosps_of_variant, hosp_fatality_rate_of_variant_by_week);
    ans ./= rep_matrix( denominator, A);
    return(ans);
  }

  matrix calculate_means_for_variant(
    // parameters
    vector prop_variant,
    matrix age_composition,
    // data
    int A,
    vector vec_all_hosps
    )
  {
    return( rep_matrix( vec_all_hosps .* prop_variant, A) .* age_composition );
  }
  
  matrix calculate_exp_hospital_admissions(
    // parameters
    matrix exp_deaths_by_week,
    vector hosp_fatality_ratio,
    matrix fr_multiplier_by_week,
    // data
    int Th,
    int Th_days,
    int A,
    int F,
    int offset_in_days_for_deaths,
    int[] weeks_to_days_for_deaths,
    matrix time_hospadmission_to_death_distr
    )
  {
    int F2;
    int F3;
    matrix[Th,A] exp_hosps_by_week;
    matrix[Th_days,A] exp_hosps_by_day;
    matrix[Th_days,A] exp_deaths_by_day;
    matrix[7,A] ones;
    matrix[Th,A] hfr;
    
    exp_deaths_by_day[1:offset_in_days_for_deaths,:] = rep_matrix(0.,offset_in_days_for_deaths,A);
    exp_deaths_by_day[(offset_in_days_for_deaths+1):Th_days,:] = exp_deaths_by_week[weeks_to_days_for_deaths,:] / 7.;
    ones = rep_matrix(1.,7,A);
    
    for(s in 1:Th_days)
    {
        F2 = min(s+F-1, Th_days);
        F3 = min(F, Th_days-s+1);
        exp_hosps_by_day[s,:] =
            columns_dot_product(
                exp_deaths_by_day[s:F2,:],
                time_hospadmission_to_death_distr[1:F3,:]
                );
    }
    
    for(s in 1:Th)
    {
        exp_hosps_by_week[s,:] =
            columns_dot_product(
                exp_hosps_by_day[(7*(s-1)+1):(7*s),:],
                ones
                );
    }
    
    hfr = fr_multiplier_by_week;
    hfr .*= rep_matrix( hosp_fatality_ratio', Th );
    hfr = fmin( 0.999, hfr);
    exp_hosps_by_week ./= hfr;
    
    return( exp_hosps_by_week );
  }
  
  matrix aggregate_latent_parameter_day_to_week(
    // parameters
    matrix exp_latent_by_day,
    // data
    int Th,
    int A
    )
  {
    matrix[Th,A] exp_latent_by_week;
    matrix[7,A] ones = rep_matrix(1.,7,A);
    for(s in 1:Th)
    {
        exp_latent_by_week[s,:] =
            columns_dot_product(
                exp_latent_by_day[(7*(s-1)+1):(7*s),:],
                ones
                );
    }
    return( exp_latent_by_week );
  }
  
  real log_dens_age_composition_of_events_by_variant(
    // parameters
    matrix age_prop_wildtype,
    matrix age_prop_P1,
    vector prop_P1,
    real v_inflation,
    // data
    int T,
    int A,
    int[,] age_events,
    int[] int_overall_events,
    vector vec_overall_events
    )
  {
    real lpmf;
    matrix[T,A] mixed_dirm_pars;
    vector[T] phi;
    phi = ( vec_overall_events - 1 ) / ( 1 + v_inflation );
    mixed_dirm_pars = rep_matrix(prop_P1,A) .* age_prop_P1;
    mixed_dirm_pars += rep_matrix(1-prop_P1,A) .* age_prop_wildtype;
    mixed_dirm_pars .*= rep_matrix( phi, A);
    lpmf = 0.0;
    for(t in 1:T)
    {
      if(int_overall_events[t]>1)
      {
        lpmf += dirichlet_multinomial_lpmf( age_events[t,:] | int_overall_events[t], mixed_dirm_pars[t,]');
      }
    }
    return(lpmf);
  }
  
  real log_dens_age_composition_of_hosps_by_variant(
    // parameters
    real age_hosps_v_inflation,
    matrix exp_hosp_adm_wildtype,
    matrix exp_hosp_adm_P1,
    // data
    int Th,
    int A,
    int[,] age_hosp_admissions
    )
  {
    real logdens;
    vector[Th*A] mean;
    mean = to_vector( (exp_hosp_adm_wildtype + exp_hosp_adm_P1)' );
    logdens = neg_binomial_2_lpmf( to_array_1d(age_hosp_admissions) | mean, mean/age_hosps_v_inflation );
    return( logdens );
  }
  
  real log_dens_p1_sequences(
    // parameters
    real propp1_logistic_growthrate,
    real propp1_logistic_midpoint,
    real propp1_overdisp_inv,
    // data
    real propp1_logistic_max_val,
    int N,
    int[] seq_sampled,
    int[] seq_P1,
    vector seq_timepoints
    )
  {
    real logdens;
    vector[N] prop_P1_fit;
    
    prop_P1_fit = propp1_logistic_max_val * inv(1.0 + exp(-propp1_logistic_growthrate * (seq_timepoints - propp1_logistic_midpoint)));
    prop_P1_fit = fmin( prop_P1_fit, 1.0-1e-10 );
    prop_P1_fit = fmax( prop_P1_fit, 1e-10 );
    
    logdens =
        beta_binomial_lpmf(
            seq_P1 |
            seq_sampled,
            prop_P1_fit / propp1_overdisp_inv,
            (1.-prop_P1_fit) / propp1_overdisp_inv
            );
    return(logdens);
  }
  
  real log_dens_deaths_out_of_hospital_admissions(
    // parameters
    matrix hfr_overall,
    real fr_overdisp_inv,
    // data
    int Th,
    int A,
    int[] fr_wildtype_deaths_row_major_order,
    int[] fr_wildtype_hosps_row_major_order,
    int[] fr_wildtype_obsidx_row_major_order
    )
  {
    real logdens;
    matrix[Th,A] hfr;
    vector[Th*A] hfr_rowmajor_long;
    
    hfr = hfr_overall;
    hfr = fmin( 0.999, hfr);
    hfr = fmax( 1e-10, hfr);
    hfr_rowmajor_long = to_vector(hfr');
    logdens =
        beta_binomial_lpmf(
            fr_wildtype_deaths_row_major_order |
            fr_wildtype_hosps_row_major_order,
            hfr_rowmajor_long[ fr_wildtype_obsidx_row_major_order ] / fr_overdisp_inv,
            (1.0 - hfr_rowmajor_long[ fr_wildtype_obsidx_row_major_order ]) / fr_overdisp_inv
            );
    return(logdens);
  }
}

data{
  //
  // data for age composition of hospitalisations
  //
  int<lower=0> A; // number age bands
  int<lower=1> Th; // number observation weeks for infections, must be larger or equal than T
  int age_hosp_admissions[Th,A]; // hospital admissions by age and week
  matrix[Th,A] age_pop_mat_underlying_hosps; // age composition of population, adjusting for deaths that already occurred, at time of death
  //
  // data for proportion of P1
  //
  int<lower=0> N;
  vector<lower=0>[N] seq_timepoints;
  int<lower=0> seq_sampled[N];
  int<lower=0> seq_P1[N];
  real propp1_logistic_midpoint_priormeansd[2];
  int<lower=1,upper=Th> seq_max_noP1_weekidx;
  //
  // data for fatality rate multiplier
  //
  int P; // number of predictors for fatality rate
  matrix[Th,P] fr_predictors[A]; // matrix of predictors for each age band
  int<lower=0, upper=Th*A> Th_wildtype_fr_obs;   // number of hosp admissions for all age groups for wildtype period, excluding those with zero total
  int fr_wildtype_deaths_row_major_order[Th_wildtype_fr_obs]; // number deaths in row-major order in matrix [weeks,age]
  int fr_wildtype_hosps_row_major_order[Th_wildtype_fr_obs]; // number hosp admissions in row-major order in matrix [weeks,age]
  int fr_wildtype_obsidx_row_major_order[Th_wildtype_fr_obs]; // index which entries in to_array_1d( (matrix[weeks,age])' ) correspond to non-zero totals
  real fr_shrinkage_df; // degrees of freedom for shrinkage prior
}

transformed data{
  int Th_days = 7 * Th;
  matrix[Th,A] log_age_pop_mat_underlying_hosps = log(age_pop_mat_underlying_hosps);
  real propp1_logistic_max_val = 1.0;
  
  vector<lower=0, upper=Th_days>[Th] first_day_of_week_idx; // day index of the beginning of a week
  int seq_noP1_weekidx[seq_max_noP1_weekidx];
  int int_all_hosps[Th];
  vector[Th] vec_all_hosps;
  
  for(t in 1:Th)
  {
    int_all_hosps[t] = sum(age_hosp_admissions[t,:]);
    vec_all_hosps[t] = 1. * int_all_hosps[t];
    first_day_of_week_idx[t] = 7*(t-1)+1;
  }
  for(t in 1:seq_max_noP1_weekidx)
  {
    seq_noP1_weekidx[t] = t;
  }
}

parameters {
  //
  // parameters for age composition of P1 and wildtype hospital admissions
  //
  vector[A] log_age_hosp_rate_wildtype;
  vector[A] log_age_hosp_rate_P1;
  real<lower=0> age_hosps_v_inflation;
  //
  // parameters for hospital fatality rate
  //
  vector[A] logit_hosp_fatality_ratio_wildtype;
  vector[A] logit_hosp_fatality_ratio_P1_rnde;
  real<lower=0> logit_hosp_fatality_ratio_P1_rnde_sd;
  //
  // parameters for proportion of P1
  //
  real propp1_logistic_growthrate;
  real <lower=0,upper=Th*7> propp1_logistic_midpoint;
  real<lower=0> propp1_overdisp_inv;
  //
  // parameters for fatality rate multiplier
  //
  vector<lower=0>[P] fr_regcoeff_overall;
  real<lower=0> fr_overdisp_inv;
  vector<lower=0>[P] fr_shrinkage;
  real<lower=0> fr_scale;
}

transformed parameters {
  matrix[Th,A] age_prop_hosps_wildtype;
  matrix[Th,A] age_prop_hosps_P1;
  matrix[Th,A] hfr_overall;
  vector<lower=0,upper=1>[Th] prop_P1;
  matrix<lower=0,upper=1>[Th,A] hfr_wildtype;
  matrix<lower=0,upper=1>[Th,A] hfr_P1;
    
  {
    prop_P1 =
        make_logistic_prop_P1(
            propp1_logistic_growthrate,
            propp1_logistic_midpoint,
            propp1_logistic_max_val,
            first_day_of_week_idx
            );
    
    hfr_wildtype =
        make_time_varying_fatality_rate_of_variant(
            logit_hosp_fatality_ratio_wildtype,
            fr_regcoeff_overall,
            A,
            Th,
            fr_predictors
            );
    
    hfr_P1 =
        make_time_varying_fatality_rate_of_variant(
            logit_hosp_fatality_ratio_wildtype+logit_hosp_fatality_ratio_P1_rnde,
            fr_regcoeff_overall,
            A,
            Th,
            fr_predictors
            );
            
    age_prop_hosps_wildtype =
        calculate_age_prop_of_hosps_by_variant(
            log_age_hosp_rate_wildtype,
            Th,
            A,
            log_age_pop_mat_underlying_hosps
            );
    
    age_prop_hosps_P1 =
        calculate_age_prop_of_hosps_by_variant(
            log_age_hosp_rate_P1,
            Th,
            A,
            log_age_pop_mat_underlying_hosps
            );
    
    hfr_overall =
        calculate_joint_prob_of_deaths_and_hospital_admission_with_variant_conditional_on_age(
            1.-prop_P1,
            age_prop_hosps_wildtype,
            age_prop_hosps_P1,
            hfr_wildtype,
            Th,
            A
            );
    
    hfr_overall +=
        calculate_joint_prob_of_deaths_and_hospital_admission_with_variant_conditional_on_age(
            prop_P1,
            age_prop_hosps_P1,
            age_prop_hosps_wildtype,
            hfr_P1,
            Th,
            A
            );
  }
  
}

model {
  //
  // priors age composition of hospital admissions
  //
  target += normal_lpdf( log_age_hosp_rate_wildtype | 0, 1 );
  target += normal_lpdf( log_age_hosp_rate_P1 | 0, 1 );
  target += exponential_lpdf( age_hosps_v_inflation | 10 );
    
  //
  // priors hospital fatality ratio
  //
  target += normal_lpdf( logit_hosp_fatality_ratio_wildtype | -.25, 1.5);
  target += normal_lpdf( logit_hosp_fatality_ratio_P1_rnde | 0, logit_hosp_fatality_ratio_P1_rnde_sd );
  target += exponential_lpdf( logit_hosp_fatality_ratio_P1_rnde_sd | 2 );
    
  //
  // priors for P1 proportion
  //
  target += normal_lpdf( propp1_logistic_growthrate | 0, .2);
  target += exponential_lpdf( propp1_overdisp_inv | 20 );
  target += normal_lpdf( propp1_logistic_midpoint | propp1_logistic_midpoint_priormeansd[1], propp1_logistic_midpoint_priormeansd[2]);  
  target += normal_lccdf( prop_P1[seq_noP1_weekidx] | 0, 0.0025 );
  
  //
  // priors for fatality rate multiplier with horseshoe variable selection
  //
  target += cauchy_lpdf( fr_scale | 0, .01); // qcauchy(0.975,0,.01) = 0.12 which is approx the upper 97.5% quantile of the posterior
  target += student_t_lupdf( fr_shrinkage | fr_shrinkage_df, 0, 1);
  target += normal_lpdf( fr_regcoeff_overall | 0, fr_shrinkage*fr_scale);
  target += exponential_lpdf( fr_overdisp_inv | 100);
  
  //
  // likelihood of age composition of hospital admissions
  //
  target += log_dens_age_composition_of_events_by_variant(
    age_prop_hosps_wildtype,
    age_prop_hosps_P1,
    prop_P1,
    age_hosps_v_inflation,
    Th,
    A,
    age_hosp_admissions,
    int_all_hosps,
    vec_all_hosps
    );
    
  //
  // likelihood of P1 genotype data
  //
  target += log_dens_p1_sequences(
    propp1_logistic_growthrate,
    propp1_logistic_midpoint,
    propp1_overdisp_inv,
    propp1_logistic_max_val,
    N,
    seq_sampled,
    seq_P1,
    seq_timepoints
    );
    
  //
  // likelihood of deaths out of hospital admissions
  //
  target += log_dens_deaths_out_of_hospital_admissions(
    hfr_overall,
    fr_overdisp_inv,
    Th,
    A,
    fr_wildtype_deaths_row_major_order,
    fr_wildtype_hosps_row_major_order,
    fr_wildtype_obsidx_row_major_order
    );
}

generated quantities{
    matrix[Th,A] fr_multiplier_by_week;
    matrix[Th,A] exp_hosp_adm_wildtype;
    matrix[Th,A] exp_hosp_adm_P1;
    matrix[Th,A] exp_deaths_wildtype_in_hosp;
    matrix[Th,A] exp_deaths_P1_in_hosp;
    matrix[Th,A] age_prop_deaths_wildtype;
    matrix[Th,A] age_prop_deaths_P1;
    
    fr_multiplier_by_week = hfr_wildtype;
    fr_multiplier_by_week ./= rep_matrix(inv_logit(logit_hosp_fatality_ratio_wildtype)',Th);
        
    exp_hosp_adm_wildtype =
      calculate_means_for_variant(
        1-prop_P1,
        age_prop_hosps_wildtype,
        A,
        vec_all_hosps
        );
    
    exp_hosp_adm_P1 =
      calculate_means_for_variant(
        prop_P1,
        age_prop_hosps_P1,
        A,
        vec_all_hosps
        );
        
    exp_deaths_wildtype_in_hosp =
      calculate_means_for_variant(
        1.-prop_P1,
        age_prop_hosps_wildtype .* hfr_wildtype,
        A,
        vec_all_hosps
        );
    
    exp_deaths_P1_in_hosp =
      calculate_means_for_variant(
        prop_P1,
        age_prop_hosps_P1 .* hfr_P1,
        A,
        vec_all_hosps
        );
        
    age_prop_deaths_wildtype =
      calculate_age_prop_of_deaths_in_hosps_by_variant(
        age_prop_hosps_wildtype,
        hfr_wildtype,
        Th,
        A
        );
    
    age_prop_deaths_P1 =
      calculate_age_prop_of_deaths_in_hosps_by_variant(
        age_prop_hosps_P1,
        hfr_P1,
        Th,
        A
        );
}

