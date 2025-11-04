data{
  // Indices
  int<lower = 1> N_dyad;  
  int<lower = 1> N_ind;  
  int<lower = 1> N_group;  
  
  // Dataset "A": all dyads, and corresponding grp and ind
  array[N_dyad] int A_dyad;        // Dyad id {1, ..., N_dyad}
  array[N_dyad] int A_grp;         // Group id {1, ..., N_group}
  array[N_dyad] int A_ind_a;       // ind a per dyad number
  array[N_dyad] int A_ind_b;       // ind b per dyad number
  
  // Dataset "B": observed states
  int<lower = 1> J;                // Number of observed states
  array[J] int<lower = 1, upper = N_group> B_grp; // Group
  vector<lower = 0>[J] B_x;                       // Holding time
  array[J] int<lower = 1, upper = 4> B_s;  // State {1, 2, 3, 4}
  array[J] int<lower = 0, upper = 1> B_c;  // censoring {0, 1}
  array[J] int<lower = 1, upper = N_dyad> B_dyad; // dyad {1, ..., N_dyad}
  
  // Dataset "C": observed transitions
  int<lower = 1> N_trans;          // Number of observed transitions
  array[N_trans] int<lower = 1, upper = N_group> C_grp;  // Group
  array[N_trans] int<lower = 1, upper = N_dyad> C_dyad;  // dyad {1, ..., N_dyad}
  array[N_trans] int<lower = 1, upper = 4> C_s_from;     // state j
  array[N_trans] int<lower = 1, upper = 4> C_s_to;       // state j + 1
}

transformed data {
  // Count of transitions. One matrix per dyad
  array[N_dyad, 4, 4] int<lower = 0> transition_counts;

  // Initialise the transition matrices with zeros
  for (dyad in 1:N_dyad) {
    for (k in 1:4) {
      for (l in 1:4) {
        transition_counts[dyad, k, l] = 0;
      }
    }
  }

  // Fill the transition matrices with transitions from the input data
  for (tr in 1:N_trans) {
    transition_counts[C_dyad[tr], C_s_from[tr], C_s_to[tr]] += 1;
  }
}

parameters {
  // Intercepts per group
  vector[N_group] alpha_1;
  vector[N_group] alpha_2;
  vector[N_group] alpha_gr;
  vector[N_group] alpha_1_gr;
  vector[N_group] alpha_2_gr;
  vector[N_group] alpha_gr_1;
  vector[N_group] alpha_gr_2;
  
  // Individual varying effects
  matrix[12, N_ind] z_phi; // z-score of phi (wide format: 12 R, N_ind C)
  vector<lower = 0>[12] sigma_phi; // SD of phi
  cholesky_factor_corr[12] L_phi; // Cholesky factor of corr. matrix phi
  
  // Dyadic varying effects
  matrix[12, N_dyad] z_tau; // z-score of tau (wide format: 12 R, N_dyad C)
  vector<lower = 0>[7] sigma_tau; // SD of tau
  cholesky_factor_corr[12] L_tau; // Cholesky factor of corr. matrix for tau
}

transformed parameters{
    // 1. individual varying effects
    matrix[N_ind, 12] phi_raw; // vertical: N_ind R, 12 C)
    phi_raw = (diag_pre_multiply(sigma_phi, L_phi) * z_phi)';
    vector[N_ind] phi_1 = phi_raw[, 1];
    vector[N_ind] phi_2 = phi_raw[, 2];
    vector[N_ind] phi_give_gr = phi_raw[, 3];
    vector[N_ind] phi_rec_gr = phi_raw[, 4];
    vector[N_ind] phi_1_give_gr = phi_raw[, 5];
    vector[N_ind] phi_1_rec_gr = phi_raw[, 6];
    vector[N_ind] phi_2_give_gr = phi_raw[, 7];
    vector[N_ind] phi_2_rec_gr = phi_raw[, 8];
    vector[N_ind] phi_give_gr_1 = phi_raw[, 9];
    vector[N_ind] phi_give_gr_2 = phi_raw[, 10];
    vector[N_ind] phi_rec_gr_1 = phi_raw[, 11];
    vector[N_ind] phi_rec_gr_2 = phi_raw[, 12];
    
  // 2. dyad varying effects
    matrix[N_dyad, 12] tau_raw; // vertical: N_ind R, 12 C
    tau_raw = (diag_pre_multiply([sigma_tau[1], sigma_tau[2], 
                                 sigma_tau[3], sigma_tau[3],
                                 sigma_tau[4], sigma_tau[4],
                                 sigma_tau[5], sigma_tau[5],
                                 sigma_tau[6], sigma_tau[7],
                                 sigma_tau[6], sigma_tau[7]], L_tau) * z_tau)';
    vector[N_dyad] tau_1 = tau_raw[, 1];
    vector[N_dyad] tau_2 = tau_raw[, 2];
    vector[N_dyad] tau_gr_ab = tau_raw[, 3];
    vector[N_dyad] tau_gr_ba = tau_raw[, 4];
    vector[N_dyad] tau_1_gr_ab = tau_raw[, 5];
    vector[N_dyad] tau_1_gr_ba = tau_raw[, 6];
    vector[N_dyad] tau_2_gr_ab = tau_raw[, 7];
    vector[N_dyad] tau_2_gr_ba = tau_raw[, 8];
    vector[N_dyad] tau_gr_1_ab = tau_raw[, 9];
    vector[N_dyad] tau_gr_2_ab = tau_raw[, 10];
    vector[N_dyad] tau_gr_1_ba = tau_raw[, 11];
    vector[N_dyad] tau_gr_2_ba = tau_raw[, 12];
  
  // We declare a bunch of varibales that we define further below  
    array[N_dyad, 4] real theta;
    array[N_dyad, 4, 4] real Psi;
    array[N_dyad, 4, 4] real Gamma;  
    
  // 3. theta
    for (dyad in 1:N_dyad) {
    theta[dyad, 1] = exp(alpha_1[A_grp[dyad]] +
                     phi_1[A_ind_a[dyad]] +
                     phi_1[A_ind_b[dyad]] +
                     tau_1[A_dyad[dyad]]);
    theta[dyad, 2] = exp(alpha_2[A_grp[dyad]] +
                     phi_2[A_ind_a[dyad]] +
                     phi_2[A_ind_b[dyad]] +
                     tau_2[A_dyad[dyad]]);
    theta[dyad, 3] = exp(alpha_gr[A_grp[dyad]] +
                     phi_give_gr[A_ind_a[dyad]] +
                     phi_rec_gr[A_ind_b[dyad]] +
                     tau_gr_ab[A_dyad[dyad]]);
    theta[dyad, 4] = exp(alpha_gr[A_grp[dyad]] +
                     phi_give_gr[A_ind_b[dyad]] +
                     phi_rec_gr[A_ind_a[dyad]] +
                     tau_gr_ba[A_dyad[dyad]]);
     
  // 4. Psi                       
    Psi[dyad, 1, 2] = 0.0;
    Psi[dyad, 1, 3] = alpha_1_gr[A_grp[dyad]] + phi_1_give_gr[A_ind_a[dyad]] +
                             phi_1_rec_gr[A_ind_b[dyad]] + tau_1_gr_ab[A_dyad[dyad]];
    Psi[dyad, 1, 4] = alpha_1_gr[A_grp[dyad]] + phi_1_give_gr[A_ind_b[dyad]] +
                             phi_1_rec_gr[A_ind_a[dyad]] + tau_1_gr_ba[A_dyad[dyad]];
                             
    Psi[dyad, 2, 1] = 0.0;
    Psi[dyad, 2, 3] = alpha_2_gr[A_grp[dyad]] + phi_2_give_gr[A_ind_a[dyad]] +
                             phi_2_rec_gr[A_ind_b[dyad]] + tau_2_gr_ab[A_dyad[dyad]];
    Psi[dyad, 2, 4] = alpha_2_gr[A_grp[dyad]] + phi_2_give_gr[A_ind_b[dyad]] +
                             phi_2_rec_gr[A_ind_a[dyad]] + tau_2_gr_ba[A_dyad[dyad]];
                             
    Psi[dyad, 3, 1] = alpha_gr_1[A_grp[dyad]] + phi_give_gr_1[A_ind_a[dyad]] +
                             phi_rec_gr_1[A_ind_b[dyad]] + tau_gr_1_ab[A_dyad[dyad]];
    Psi[dyad, 3, 2] = alpha_gr_2[A_grp[dyad]] + phi_give_gr_2[A_ind_a[dyad]] +
                             phi_rec_gr_2[A_ind_b[dyad]] + tau_gr_2_ab[A_dyad[dyad]];
    Psi[dyad, 3, 4] = 0.0;
    
    Psi[dyad, 4, 1] = alpha_gr_1[A_grp[dyad]] + phi_give_gr_1[A_ind_b[dyad]] +
                             phi_rec_gr_1[A_ind_a[dyad]] + tau_gr_1_ba[A_dyad[dyad]];
    Psi[dyad, 4, 2] = alpha_gr_2[A_grp[dyad]] + phi_give_gr_2[A_ind_b[dyad]] +
                             phi_rec_gr_2[A_ind_a[dyad]] + tau_gr_2_ba[A_dyad[dyad]];
    Psi[dyad, 4, 3] = 0.0;

  // 5. Gamma
    for (k in 1:4) {
      real row_sum = 0.0;
  
      // Compute the denominator (sum of exponents in row k, excluding diagonal)
      for (m in 1:4) {
        if (m != k) {
          row_sum += exp(Psi[dyad, k, m]);
        } // end if
      } // end for m
  
      // Now compute Gamma for each element in row k
      for (l in 1:4) {
        if (k == l) {
          Gamma[dyad, k, l] = 0.0;  // Diagonal elements are 0
        } else {
          Gamma[dyad, k, l] = exp(Psi[dyad, k, l]) / row_sum;
        } // end if
      } // end for l
    } // end for k
  } // end for dyad
}

model{
  // 1. Likelihood for holding times
  for (j in 1:J){
        // not censored
        if (B_c[j] == 0){
        target += exponential_lpdf(B_x[j] | 1 / theta[B_dyad[j], B_s[j]]);
        } // end non-censored obs
      
        // censored
        if (B_c[j] == 1){
        target += exponential_lccdf(B_x[j] | 1 / theta[B_dyad[j], B_s[j]]);
        } // end censored obs
    } // end for j
    
  // 2. Likelihood for observed transitions
  for (tr in 1:N_trans) {
    target += categorical_lpmf(C_s_to[tr] | 
              to_vector(Gamma[C_dyad[tr], C_s_from[tr]]));
  } // end for tr
  
  // 3. Hyper-priors
    target += normal_lpdf(alpha_1 | 8, 2);
    target += normal_lpdf(alpha_2 | 1.5, 1);
    target += normal_lpdf(alpha_gr | 1.5, 1);
    target += normal_lpdf(alpha_1_gr | 0, 1);
    target += normal_lpdf(alpha_2_gr | 0, 1);
    target += normal_lpdf(alpha_gr_1 | 0, 1);
    target += normal_lpdf(alpha_gr_2 | 0, 1);
    target += normal_lpdf(to_vector(z_phi) | 0, 1);
    target += normal_lpdf(to_vector(z_tau) | 0, 1);
    target += exponential_lpdf(sigma_phi | 1);
    target += exponential_lpdf(sigma_tau | 1);
    target += lkj_corr_cholesky_lpdf(L_phi | 4);
    target += lkj_corr_cholesky_lpdf(L_tau | 4);
} // end model block

generated quantities{
  // Corr. matrix from cholesky factors
    // Varying individual effects
    matrix[12, 12] c_ind;
    c_ind = multiply_lower_tri_self_transpose(L_phi);
  
    // Varying dyadic effects
    matrix[12, 12] c_dyad;
    c_dyad = multiply_lower_tri_self_transpose(L_tau);
    
  // Avg. holding time per state (across dyads)
  vector[3] avg_theta;
  avg_theta[1] = mean(theta[, 1]);
  avg_theta[2] = mean(theta[, 2]);
  avg_theta[3] = mean(append_row(
                      to_vector(theta[, 3]), to_vector(theta[, 4])));
  
  // Avg. transition probabilities (across dyads)
  matrix[3, 3] avg_Gamma;
  avg_Gamma[1, 2] = mean(Gamma[, 1, 2]);
  avg_Gamma[1, 3] = mean(Gamma[, 1, 3]) + mean(Gamma[, 1, 4]);
  
  avg_Gamma[2, 1] = mean(Gamma[, 2, 1]);
  avg_Gamma[2, 3] = mean(Gamma[, 2, 3]) + mean(Gamma[, 2, 4]);
  
  avg_Gamma[3, 1] = mean(append_row(
                      to_vector(Gamma[, 3, 1]), to_vector(Gamma[, 4, 1])));
  avg_Gamma[3, 2] = mean(append_row(
                      to_vector(Gamma[, 3, 2]), to_vector(Gamma[, 4, 2])));
  avg_Gamma[3, 3] = mean(append_row(
                      to_vector(Gamma[, 3, 4]), to_vector(Gamma[, 4, 3])));
}


