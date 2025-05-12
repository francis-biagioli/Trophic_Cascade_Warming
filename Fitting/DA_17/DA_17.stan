///////////////////////////////////////////////////////////////////////////////////////
////////////////////////// STAN CODE FOR TWO SPECIES JARS 17C /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

////////////////////////// DEFINE ODE MODEL /////////////////////////

functions{
  //dR/dt = rR(1-(RQe^Bt))-aRC/1+ahR-w(C-0.0033)  //Alge model
  //dc/dt = eaRC/1+ahR-w(C-0.0033) - dC        //Daphnia model
 
  real[] odemodel(real t, real[] N, real[] p, real[] x_r, int[] x_i) {
    //p[1] = r, p[2] = Q, p[3] = a, p[4] = e, p[5] = d, p[6] = h, p[7] = w, p[8] = B
    real dNdt[2];
    dNdt[1] = (p[1]*N[1]*(1-N[1]*p[2]*exp(p[8]*t)))-(p[3]*N[1]*N[2]/(1+p[7]*(N[2]-0.003333333)+p[3]*p[6]*N[1]));
    dNdt[2] = (p[4]*p[3]*N[1]*N[2]/(1+p[7]*(N[2]-0.003333333)+p[3]*p[6]*N[1]))-p[5]*N[2];
    return dNdt;
  }
}

////////////////////////// FORMAT DATA FOR STAN SOLVER /////////////////////////

data {
  int n; //number of expressions for Daph (most frequent)
  int m; //number of replicates
  
  real t[n]; //time steps
  
  real N1[m,n]; //left is replicate, right is observbations
  real N2[m,n]; 
  
}

////////////////////////// DEFINE PARAMETERS TO BE ESTIMATED BY SOLVER /////////////////////////

parameters {
real<lower = 0, upper = 2> r; //Algae per capita growth rate

 real<lower = 1.0e-8, upper = 1.0e-3>  Q; //Algae inverse carrying capacity
 
 real<lower = 0, upper = 2> a1; //ceriodaphnia space clearance rate
 
 real<lower = 0, upper= 0.01> e1; //ceriodaphnia conversion efficiency
 
 real<lower = 0> d1; //ceriodaphnia background mportality rate
 
 real<lower = 0, upper = 1.0e-3> h1; //cerioaphnia handling time
 
 real<lower = 0> w1; //ceriodaphnia Beddington-DeAngelis interferance
 
 real<lower = 0, upper = 0.4> B; //correctional carrying capacity term
 
 real<lower = 0> sd_obs1; //std dev of algae
 real<lower = 0> sd_obs2; //std dev of daph
 
 real<lower = 0> N1_0[m]; //estimating initial value of Algae per replicate jar
 real<lower = 0> N2_0[m]; //estimating initial value of Ceriodaphnia per replicate jar
}

////////////////////////// DEFINE PPRIORS FOR FIT /////////////////////////

model {
  real p[8]; //number of params
  real Nsim[n-1, 2]; //dimentions of saved Nsim matrix
  
  r ~ normal(9.210132e-02, 1.331021e-02);
  
  Q ~ normal(8.391139e-08, 6.167348e-08);
  
  a1 ~ normal(2.660032e-01, 1.218042e-01);
  
  e1 ~ normal(2.247731e-05, 6.912190e-06); 
  
  d1 ~ normal (1.053474e+00, 2.958796e-01);
  
  h1 ~ normal(4.229568e-06, 3.015631e-06);
 
  w1 ~ normal(1.885828e+00, 3.860036e-01);
  
  B ~ normal(7.973137e-02, 3.010989e-02);

  sd_obs1 ~ exponential(1e-5);
  sd_obs2 ~ exponential(1);
  
  for (i in 1:m) {
  N1_0[m] ~ normal(N1[m,1],1.5e4);
  N2_0[m] ~ normal(N2[m,1],0.2);
  }
  
  ////////////////////////// DEFINE PARS FOR INTEGRATOR ///////////////////////// 
  
  p[1] = r;
  p[2] = Q;
  p[3] = a1;
  p[4] = e1;
  p[5] = d1;
  p[6] = h1;
  p[7] = w1;
  p[8] = B;
  
  ///////////////////////////////////// INTEGRATE ODE ///////////////////////////// 
  
  for (j in 1:m){
  Nsim = integrate_ode_rk45(odemodel, {N1_0[j], N2_0[j]}, t[1], t[2:n], p, rep_array(0.0, 0), rep_array(0,0));
  
////////////////////////// LIKELIEHOOD STATEMENT /////////////////////////

  N1[j,1] ~ normal(N1_0[j], sd_obs1);
  N2[j,1] ~ normal(N2_0[j], sd_obs2);

 for(i in 2:n) {
   
   N2[j,i] ~ normal(Nsim[i-1, 2], sd_obs2); 
   N1[j,i] ~ normal(Nsim[i-1, 1], sd_obs1);
    

  }
}

}
