functions {
  real[] mingled(real t, real[] z, real[] par, real[] x_r, int[] x_i) {
    real dzdt[2];
    
    // food intake and growth
    real dFdt = par[1];
    real dWdF = par[2] * (1 - z[2] / par[3]);
    
    dzdt[1] = dFdt;        // food intake (g/day)
    dzdt[2] = dWdF * dFdt; // body weight (g)
    return dzdt;
  }
  
  vector snapshot( vector phi, vector init_par, real[] x_r, int[] x_i ) {
    int N_days = x_i[1];
    int N_metabolites = x_i[2];
    int N_pars = x_i[3];
    real z0[N_metabolites];
    real par[N_pars];
    real z_hat[N_days,N_metabolites];
  
    for ( i in 1:N_metabolites ) {
      z0[i] = init_par[i];
    }
    
    for ( i in 1:N_pars ) {
      par[i] = init_par[N_metabolites + i];
    }
    
    z_hat = integrate_ode_bdf(mingled,  z0, x_r[1], x_r[2:(N_days+1)], par, x_r, x_i);
    return to_vector(to_matrix(z_hat));
  }
}
data {
  int<lower=1> N;                            // number of rows database-like data structure
  real<lower=0> y[N];
  int<lower=0> n_mouse;                      // number of mice
  int<lower=1,upper=n_mouse> mouse_id[N];
  int<lower=1> n_day;                        // number of time points to solve the ODE
  real ts[n_day];                            // days that ODE model solves
  int<lower=1,upper=n_day> day_id[N];
  int<lower=1> n_metabolite;
  int<lower=1,upper=n_metabolite> metabolite_id[N];
  int<lower=1> n_obs_metabolite; // # of metabolites with measurements
  int<lower=1,upper=n_obs_metabolite> obs_metabolite_id[N];
  int<lower=0> n_bw0_obs;
  int<lower=0> n_bw0_mis;
  int<lower=1,upper=n_mouse> which_bw0_obs[n_bw0_obs];
  int<lower=1,upper=n_mouse> which_bw0_mis[n_bw0_mis];
  vector<lower=0>[n_bw0_obs] bw0_obs;
}
transformed data {
  int n_day_pred = 6; // # of days in predictions
  int n_par = 3; // # of parameters used in ODEs
  real t0 = 0; // initial time point
  real F0 = 0;
  vector[0] phi;
  real x_r[n_mouse,1 + n_day];
  int x_i[n_mouse,3];
  real ts_pred[n_day_pred];
  real x_r_pred[0];
  int x_i_pred[0];
  
  for ( i in 1:n_day_pred ) {
    ts_pred[i] = 30 * i;
  }
  
  for ( i in 1:n_mouse ) {
    x_r[i,1] = t0;
    for ( j in 1:n_day ) {
      x_r[i,j+1] = ts[j];
    }
    x_i[i,1] = n_day;
    x_i[i,2] = n_metabolite;
    x_i[i,3] = n_par;
  }
}
parameters {
  vector<lower=0>[n_bw0_mis] bw0_mis;
  real<lower=0> mu_bw0;
  real<lower=0> sigma_bw0;
  vector<lower=0>[n_mouse] fi; // food intake g/day
  real<lower=0> mu_fi; 
  real<lower=0> sigma_fi;
  vector<lower=0,upper=1>[n_mouse] fer_max; // maximal feed efficiency ratio (g BW gain/g food)
  real<lower=0> a_fer_max;
  real<lower=0> b_fer_max;
  vector<lower=0>[n_mouse] bw_max; // maximal body weight
  real<lower=0> mu_bw_max;
  real<lower=0> sigma_bw_max;
  vector<lower=0>[n_metabolite] sigma;
}
transformed parameters {
  real<lower=0> par[n_mouse,n_par]; // each mouse has its own set of parametersW
  
  for ( i in 1:n_mouse ) {
    par[i,1] = fi[i];
    par[i,2] = fer_max[i];
    par[i,3] = bw_max[i];
  }
}
model {
  vector[n_mouse] bw0;
  real z0[n_mouse,n_metabolite];
  vector[n_metabolite + n_par] init_par[n_mouse];
  vector[n_mouse * n_day * n_metabolite] states;
  real z_hat[n_mouse,n_day,n_metabolite]; // pool size
  
  bw0[which_bw0_obs] = bw0_obs;
  bw0[which_bw0_mis] = bw0_mis;
  
  bw0 ~ normal( mu_bw0, sigma_bw0 ); 
  mu_bw0 ~ normal( 32, 1.2 );
  sigma_bw0 ~ lognormal( -0.06, 0.01 ); 
  
  for ( i in 1:n_mouse ) {
    z0[i,1] = F0;
    z0[i,2] = bw0[i];
  }
  
  fi ~ normal( mu_fi, sigma_fi );
  mu_fi ~ normal( 4.4, 0.1 );  // 4.4 g/day N = 28
  sigma_fi ~ lognormal( -1, 0.01 ); 
  
  fer_max ~ beta( a_fer_max, b_fer_max );
  a_fer_max ~ exponential( 1 );
  b_fer_max ~ exponential( 1 );
  
  bw_max ~ normal( mu_bw_max, sigma_bw_max );
  mu_bw_max ~ normal( 49.7, 1.1 ); // N = 51
  sigma_bw_max ~ lognormal( 0, 0.01 ); // 
  
  sigma ~ exponential( 1 );
  
  // assemble initial values and parameters
  for ( i in 1:n_mouse ) {
    for ( m in 1:n_metabolite ) {
      init_par[i,m] = z0[i,m];
    }
    
    for ( p in 1:n_par ) {
      init_par[i,n_metabolite + p] = par[i,p];
    }
  }
  
  // paralleize the model computations
  states = map_rect(snapshot, phi, init_par, x_r, x_i);
  
  // unpack the long vector into the matrix
  for ( i in 1:n_mouse ) {
    int start = n_day * n_metabolite * (i - 1);
    for ( m in 1:n_metabolite ) {
      z_hat[i,,m] = to_array_1d(states[(start + (m - 1) * n_day + 1) : (start + m * n_day)]);
    }
  }
  
  // likelihood
  for ( n in 1:N ) {
    y[n] ~ normal( z_hat[mouse_id[n], day_id[n], metabolite_id[n]], sigma[obs_metabolite_id[n]] );
  }
}
generated quantities {
  real bw0_pred;
  real z0_pred[n_metabolite];
  real par_pred[n_par];
  real z_pred[n_day_pred, n_metabolite];
  real y_pred[1+n_day_pred, n_metabolite];
  
  
  bw0_pred = normal_rng( mu_bw0, sigma_bw0 );
  
  z0_pred[1] = 0;
  z0_pred[2] = bw0_pred;
  
  par_pred[1] = normal_rng( mu_fi, sigma_fi );
  par_pred[2] = beta_rng( a_fer_max, b_fer_max );
  par_pred[3] = normal_rng( mu_bw_max, sigma_bw_max );
  
  z_pred = integrate_ode_bdf(mingled, z0_pred, t0, ts_pred, par_pred, x_r_pred, x_i_pred);
  
  y_pred[1,1] = z0_pred[1];
  y_pred[1,2] = normal_rng(z0_pred[2], sigma[2]);
  for ( i in 1:n_day_pred ) {
    y_pred[i+1,1] = normal_rng(z_pred[i,1], sigma[1]);
    y_pred[i+1,2] = normal_rng(z_pred[i,2], sigma[2]);
  }
}
