///////////////////////////////////////////////////////////////////////////////////////
////////////////////////// STAN CODE FOR THREE SPECIES JARS 23C /////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

////////////////////////// DEFINE ODE MODEL /////////////////////////
functions{
  //dR/dt = rR(1-(RQe^(Bt)))-a1RC/1+a1h1R+w1C - Alge model
  //dc/dt = e1a1RC/1+a1h1R+w1C - a2CH/1+a2h2C+w2H - d1C
  //dHdt = e2a2CH/1+a2h2C+w2H - d2H
 
  real[] odemodel(real t, real[] N, real[] p, real[] x_r, int[] x_i) {
   //p[1] = r, p[2] = Q, p[3] = a1, p[4] = e1, p[5] = d1, p[6] = h1, p[7] = w1, p[8] = B,
   //p[9] = a2, p[10] = e2, p[11] = d2, p[12] = h2, p[13] = w2

    real dNdt[3];
    dNdt[1] = (p[1]*N[1]*(1-N[1]*p[2]*exp(p[8]*t)))-(p[3]*N[1]*N[2]/(1+p[7]*(N[2]-0.003333333)+p[3]*p[6]*N[1]));
    
    dNdt[2] = ((p[4]*p[3]*N[1]*N[2])/((1+p[7]*(N[2]-0.003333333))+(p[3]*p[6]*N[1])))-((p[9]*N[2]*N[3])/((1+p[13]*(N[3]-0.003333333))+(p[9]*p[12]*N[2])))-(p[5]*N[2]);
    
    dNdt[3] = ((p[10]*p[9]*N[2]*N[3])/((1+p[13]*(N[3]-0.003333333))+(p[9]*p[12]*N[2])))-(p[11]*N[3]);
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
  real N3[m,n]; 
  
}

////////////////////////// DEFINE PARAMETERS TO BE ESTIMATED BY SOLVER /////////////////////////

 parameters {
 real<lower = 0, upper = 2> r; //Algae per capita growth rate
 
 real<lower = 1.0e-8, upper = 1.0e-3>  Q; //Algae inverse carrying capacity
 
 real<lower = 0> a1; //ceriodaphnia space clearance rate
 
 real<lower = 0, upper = 1.0e-3> h1; //ceriodaphnia handling time
 
 real<lower = 0> w1; //ceriodaphnia Beddington-DeAngelis interferance
 
 real<lower = 0> a2; //Hydra space clearance rate
 
 real<lower = 0, upper =1> e2; //Hydra conversion efficiency
 
 real<lower = 0, upper = 1> d2; //Hydra background mortality rate
 
 real<lower = 0, upper = 1> h2; //Hydra handling time
 
 real<lower = 0> w2; //Hydra Beddington-DeAngelis interferance
 
 real<lower = 0> sd_obs1; //std dev of algae
 real<lower = 0> sd_obs2; //std dev of daph
 real<lower = 0> sd_obs3; //std dev of hydra
 
 real<lower = 0> N1_0[m]; //estimating initial value of Algae
 real<lower = 0> N2_0[m]; //estimating initial value of Algae
 real<lower = 0> N3_0[m]; //estimating initial value of Algae

}

////////////////////////// DEFINE PPRIORS FOR FIT /////////////////////////

model {
  real p[13]; //number of params
  real Nsim[n-1, 3]; //dimentions of saved Nsim matrix
  
  //Priors
  //left num is mean, right numbner is SD
  
  r ~ normal(2.438321e-01, 4.169177e-02); //Prior from DA_23
  
  Q ~ normal(1.290299e-07, 7.812755e-08); //Prior from DA_23

  a1 ~ normal(5, 3);
  
  h1 ~ normal(8.937421e-07, 7.737500e-07); //Prior from DA_23
 
  w1 ~ normal(2.322578e+00, 5.430553e-01); //Prior from DA_23
  
  a2 ~ normal(3.195376e+01, 1.869069e+01);
  
  e2 ~ normal(8.400720e-02, 2.865775e-02);
  
  d2 ~ normal(2.576516e-02, 6.504414e-03);

  h2 ~ normal(1.895108e-01, 7.294823e-02);
  
  w2 ~ normal(1.022434e+01, 1.213240e+01);
  
for (i in 1:m) {
  N1_0[m] ~ normal(N1[m,1],1.5e4);
  N2_0[m] ~ normal(N2[m,1],0.2);
  N3_0[m] ~ normal(N3[m,1],0.1);
  }
  
  sd_obs1 ~ exponential(1e-5);
  sd_obs2 ~ exponential(1);
  sd_obs3 ~ exponential(1);
  
    ////////////////////////// DEFINE PARS FOR INTEGRATOR ///////////////////////// 
  
  ///////////////////////// fixed vals from median of DA_23_fit////////////////////
  p[1] = r;
  p[2] = Q;
  p[3] = a1;
  p[4] = 1.021355e-05; //e1
  p[5] = 1.943547e+00; //d1
  p[6] = h1;
  p[7] = w1;
  p[8] = 7.738450e-02; //B
  p[9] = a2;
  p[10] = e2;
  p[11] = d2;
  p[12] = h2;
  p[13] = w2;
  
  
  ///////////////////////////////////// INTEGRATE ODE ///////////////////////////// 
  
  for (j in 1:m){
  Nsim = integrate_ode_rk45(odemodel, {N1_0[j], N2_0[j], N3_0[j]}, t[1], t[2:n], p, rep_array(0.0, 0), rep_array(0,0));
  
  ////////////////////////// LIKELIHOOD STATEMENT /////////////////////////
  
  N3[j,1] ~ normal(N3_0[j], sd_obs3);
  N2[j,1] ~ normal(N2_0[j], sd_obs2);
  N1[j,1] ~ normal(N1_0[j], sd_obs1);

 for(i in 2:n) {
   
   N3[j,i] ~ normal(Nsim[i-1, 3], sd_obs3);
   N2[j,i] ~ normal(Nsim[i-1, 2], sd_obs2); 
   N1[j,i] ~ normal(Nsim[i-1, 1], sd_obs1);
    

  }
}

}
