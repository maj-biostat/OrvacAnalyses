data {                          
  int<lower=0> N;                // number of observations
  int<lower=1,upper=250> N_pred;

  int y[N];                         // response variable
  int <lower=1>K;
  matrix[N, K] X;                 // design matrix for trt, 0 = ctl 1 = trt
  matrix[N_pred, K] X_new;    // design matrix for future predictions
  
  
}
transformed data {
    //vector[N] age_sq;              // create new variable (4), age 
                                  // squared (no dots in the variable name)
  //age_sq <- age .* age;          // formula for the variable, do not 
                                // forget the . before multiplication
}
parameters {
  real alpha;
  vector[K] b;
}
transformed parameters {
  vector[N] eta;
  vector[N_pred] eta_new;
  eta = alpha +  X * b;
  eta_new = alpha +  X_new * b;

}
model {

  // priors including all constants
  target += student_t_lpdf(alpha | 3, 0, 2.5);// intercept 
  target += student_t_lpdf(b | 3, 0, 2.5);
  
  // likelihood including all constants
  target += bernoulli_logit_lpmf(y | eta);

}
generated quantities { 
  
  /*
  predict the values for those unenrolled
  */
  vector[N_pred] y_rep;

  for (i in 1:N_pred) {
    y_rep[i] = bernoulli_rng(inv_logit(eta_new[i]));
  }

}
