// DSE - Full Likelihood (Coverage and Covariate Distribution Parameters)
//
// p(D_obs | N, alpha1, alpha2, beta1, beta2, theta) \propto ... [see Plan C paper]
// 
// prior alpha1 ~ Normal(mu1, sigma1^2)  (non-zero means = X_clus1 %*% gamma1)*
// prior alpha2 ~ Normal(mu2, sigma2^2)  (non-zero means = X_clus2 %*% gamma2)*
// prior beta1  ~ Normal(mb1, sb1)
// prior beta2  ~ Normal(mb2, sb2)
// prior gamma1 ~ Normal(mg1, sg1)
// prior gamma2 ~ Normal(mg2, sg2)
// prior sigma1 ~ uniform[a1, b1] (constant density)
// prior sigma2 ~ uniform[a2, b2] (constant density)
// prior theta  ~ Dirichlet(c)
// prior N ~ Jeffreys

data {
  int<lower = 0> n;        // number of unique covariate combinations
  int<lower = 1> n_clus;   // number of clusters (geographical area)
  int<lower = 1> q1;       // number of individual-level covariate coefficients (list 1)
  int<lower = 1> q2;       // number of individual-level covariate coefficients (list 2)
  int<lower = 1> q_clus1;  // number of cluster-level covariate coefficients (list 1)
  int<lower = 1> q_clus2;  // number of cluster-level covariate coefficients (list 2)
  matrix[n, q1]  X1;       // covariate combinations design matrix (list 1)
  matrix[n, q2]  X2;       // covariate combinations design matrix (list 2
  matrix[n, 3] Y;          // observed counts in the three observed cells
  matrix[n_clus, q_clus1] X_clus1; // cluster covariates design matrix (list 1)
  matrix[n_clus, q_clus2] X_clus2; // cluster covariates design matrix (list 2)
  int<lower = 1, upper = n_clus> cluster[n]; // cluster IDs for individuals
  
  // prior mean and std. dev. of individual-level coefficients vector
  vector[q1] coefs_prior_mean1;
  vector[q2] coefs_prior_mean2;
  vector<lower = 0>[q1] coefs_prior_sd1;
  vector<lower = 0>[q2] coefs_prior_sd2;
  
  // prior mean and std. dev. of cluster-level coefficients vector
  vector[q_clus1] clus_coefs_prior_mean1;
  vector[q_clus2] clus_coefs_prior_mean2;
  vector<lower = 0>[q_clus1] clus_coefs_prior_sd1;
  vector<lower = 0>[q_clus2] clus_coefs_prior_sd2;
  
  // bounds for uniform prior of cluster std. deviations
  real<lower = 0.01> clus_sd_prior_min1;
  real<lower = clus_sd_prior_min1> clus_sd_prior_max1;
  real<lower = 0.01> clus_sd_prior_min2;
  real<lower = clus_sd_prior_min2> clus_sd_prior_max2;
  
  // prior parameters for covariate distribution
  vector<lower = 0>[n] theta_prior;
}

transformed data{
  // cell counts over covariate combinations (check order of Y columns!)
  vector[n] y11 = col(Y, 1);
  vector[n] y10 = col(Y, 2);
  vector[n] y01 = col(Y, 3);
  vector[n] counts = y11 + y10 + y01; // freq. of obs covariate combinations
  real n_obs = sum(counts);
}

parameters {
  // level 1:
  vector[q1] coefs1;               // std.ized covariate coefficients (list 1)
  vector[q2] coefs2;               // std.ized covariate coefficients (list 2)
  vector[n_clus] cluster_effects1; // std.ized cluster effects (list 1)
  vector[n_clus] cluster_effects2; // std.ized cluster effects (list 2)
  
  // level 2 (hyperparameters):
  real<lower = clus_sd_prior_min1, upper = clus_sd_prior_max1> clus_sd1;
  real<lower = clus_sd_prior_min2, upper = clus_sd_prior_max2> clus_sd2;
  
  vector[q_clus1] clus_coefs1;  // std.ized cluster covariate coefficients (list 1)
  vector[q_clus2] clus_coefs2;  // std.ized cluster covariate coefficients (list 2)
  
  // observed covariate distribution
  simplex[n] theta;
  
  // total population count
  real<lower = n_obs> N;
}

model {
  // list indicator probabilities
  vector[n] p1 = inv_logit(X1 * coefs1 + cluster_effects1[cluster]);
  vector[n] p2 = inv_logit(X2 * coefs2 + cluster_effects2[cluster]);
  
  // cell indicator probabilities
  vector[n] p11 = p1 .* p2;
  vector[n] p10 = p1 .* (1 - p2);
  vector[n] p01 = (1 - p1) .* p2;
  vector[n] p00 = (1 - p1) .* (1 - p2);
  
  // expected probability of being missed by both lists (expectation over X)
  real p00_ave = dot_product(p00, theta);
  
  // priors:
  target += normal_lpdf(coefs1 | coefs_prior_mean1, coefs_prior_sd1);
  target += normal_lpdf(coefs2 | coefs_prior_mean1, coefs_prior_sd2);
  target += normal_lpdf(cluster_effects1 | X_clus1 * clus_coefs1, clus_sd1);
  target += normal_lpdf(cluster_effects2 | X_clus2 * clus_coefs2, clus_sd2);
  target += normal_lpdf(clus_coefs1 | clus_coefs_prior_mean1, clus_coefs_prior_sd1);
  target += normal_lpdf(clus_coefs2 | clus_coefs_prior_mean2, clus_coefs_prior_sd2);
  target += dirichlet_lpdf(theta | theta_prior);
  // target += uniform(N, 2 * N);  // Jeffreys prior \propto 1/N

  // likelihood:
  target += sum(y11 .* log(p11) +
                y10 .* log(p10) +
                y01 .* log(p01) +
                counts .* log(theta)) +
            lgamma(N) -   // lgamma(N + 1) if jeffreys prior specified
            lgamma(N - n_obs + 1) +
            (N - n_obs) * log(p00_ave);
}
