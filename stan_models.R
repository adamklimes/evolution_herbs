# Stan_models

#_______________________________________________________________________
# Model with starting values
mod_code_startval <- "
data {
  int Nspec;
  int Ntot;
  vector[Ntot] FinSize;
  vector[Ntot] StaSize;
  vector[Ntot] Treat;
  vector[Nspec] Woody;
  cov_matrix[Nspec] phy;
  int DruhP[Nspec];
}
parameters {
  vector[Nspec] Apom;
  vector[Nspec] B1pom;
  vector[Nspec] B2pom;
  vector<lower=0>[Nspec] sig;
  real muA;
  real muB1;
  real alphaB2;
  real betaB2;
  real<lower=0> CparA;
  real<lower=0> CparB1;
  real<lower=0> CparB2;
  real<lower=0, upper=1> lambdaA;
  real<lower=0, upper=1> lambdaB1;
  real<lower=0, upper=1> lambdaB2;
}
model {
  vector[Nspec] A;
  vector[Nspec] B1;
  vector[Nspec] B2;
  vector[Nspec] vec_muA;
  vector[Nspec] vec_muB1;
  vector[Nspec] vec_muB2;
  matrix[Nspec, Nspec] phyloA;
  matrix[Nspec, Nspec] phyloB1;
  matrix[Nspec, Nspec] phyloB2;
  vector[Ntot] pommu;
  vector[Ntot] pomsig;
  matrix[Nspec, Nspec] AL;
  matrix[Nspec, Nspec] B1L;
  matrix[Nspec, Nspec] B2L;

  vec_muA = rep_vector(muA, Nspec);
  vec_muB1 = rep_vector(muB1, Nspec);
  vec_muB2 = alphaB2 + betaB2 * Woody;
  phyloA = CparA*((phy-diag_matrix(diagonal(phy)))*lambdaA + diag_matrix(diagonal(phy)));
  phyloB1 = CparB1*((phy-diag_matrix(diagonal(phy)))*lambdaB1 + diag_matrix(diagonal(phy)));
  phyloB2 = CparB2*((phy-diag_matrix(diagonal(phy)))*lambdaB2 + diag_matrix(diagonal(phy)));
  AL = cholesky_decompose(phyloA);
  A = vec_muA + AL * Apom;
  B1L = cholesky_decompose(phyloB1);
  B1 = vec_muB1 + B1L * B1pom;
  B2L = cholesky_decompose(phyloB2);
  B2 = vec_muB2 + B2L * B2pom;
  {
  int pom;
  pom = 1;
    for (i in 1:Nspec){
      pommu[pom:(pom+DruhP[i]-1)] = A[i] + B1[i] * StaSize[pom:(pom+DruhP[i]-1)] + B2[i] * Treat[pom:(pom+DruhP[i]-1)];
      pomsig[pom:(pom+DruhP[i]-1)] = rep_vector(sig[i], DruhP[i]);
      pom  += DruhP[i];
    }
  }

  //priors
  B1 ~ cauchy(0,5);
  B2 ~ cauchy(0,5);
  sig ~ cauchy(0,5);
  muB1 ~ cauchy(0,5);
  alphaB2 ~ cauchy(0,5);
  betaB2 ~ cauchy(0,5);
  CparA ~ cauchy(0,5);
  CparB1 ~ cauchy(0,5);
  CparB2 ~ cauchy(0,5);

  FinSize ~ normal(pommu, pomsig);
  Apom ~ std_normal();
  B1pom ~ std_normal();
  B2pom ~ std_normal();
}
"

#_______________________________________________________________________
# Model without staring values
mod_code <- "
data {
  int Nspec;
  int Ntot;
  vector[Ntot] FinSize;
  vector[Ntot] Treat;
  vector[Nspec] Woody;
  cov_matrix[Nspec] phy;
  int DruhP[Nspec];
}
parameters {
  vector[Nspec] Apom;
  vector[Nspec] B1pom;
  vector<lower=0>[Nspec] sig;
  real muA;
  real alphaB1;
  real betaB1;
  real<lower=0> CparA;
  real<lower=0> CparB1;
  real<lower=0, upper=1> lambdaA;
  real<lower=0, upper=1> lambdaB1;
}
model {
  vector[Nspec] A;
  vector[Nspec] B1;
  vector[Nspec] vec_muA;
  vector[Nspec] vec_muB1;
  matrix[Nspec, Nspec] phyloA;
  matrix[Nspec, Nspec] phyloB1;
  vector[Ntot] pommu;
  vector[Ntot] pomsig;
  matrix[Nspec, Nspec] AL;
  matrix[Nspec, Nspec] B1L;

  vec_muA = rep_vector(muA, Nspec);
  vec_muB1 = alphaB1 + betaB1 * Woody;
  phyloA = CparA*((phy-diag_matrix(diagonal(phy)))*lambdaA + diag_matrix(diagonal(phy)));
  phyloB1 = CparB1*((phy-diag_matrix(diagonal(phy)))*lambdaB1 + diag_matrix(diagonal(phy)));
  AL = cholesky_decompose(phyloA);
  A = vec_muA + AL * Apom;
  B1L = cholesky_decompose(phyloB1);
  B1 = vec_muB1 + B1L * B1pom;
  {
  int pom;
  pom = 1;
    for (i in 1:Nspec){
      pommu[pom:(pom+DruhP[i]-1)] = A[i] + B1[i] * Treat[pom:(pom+DruhP[i]-1)];
      pomsig[pom:(pom+DruhP[i]-1)] = rep_vector(sig[i], DruhP[i]);
      pom  += DruhP[i];
    }
  }

  //priors
  B1 ~ cauchy(0,5);
  sig ~ cauchy(0,5);
  alphaB1 ~ cauchy(0,5);
  betaB1 ~ cauchy(0,5);
  CparA ~ cauchy(0,5);
  CparB1 ~ cauchy(0,5);

  FinSize ~ normal(pommu, pomsig);
  Apom ~ std_normal();
  B1pom ~ std_normal();
}
"

#_______________________________________________________________________
# Model of mortality
mod_code_mort <- "
data {
  int Nspec;
  int Ntot;
  int FinSize[Ntot];
  vector[Ntot] Treat;
  vector[Nspec] Woody;
  cov_matrix[Nspec] phy;
  int DruhP[Nspec];
}
parameters {
  vector[Nspec] Apom;
  vector[Nspec] B1pom;
  real muA;
  real alphaB1;
  real betaB1;
  real<lower=0> CparA;
  real<lower=0> CparB1;
  real<lower=0, upper=1> lambdaA;
  real<lower=0, upper=1> lambdaB1;
}
model {
  vector[Nspec] A;
  vector[Nspec] B1;
  vector[Nspec] vec_muA;
  vector[Nspec] vec_muB1;
  matrix[Nspec, Nspec] phyloA;
  matrix[Nspec, Nspec] phyloB1;
  vector[Ntot] pommu;
  matrix[Nspec, Nspec] AL;
  matrix[Nspec, Nspec] B1L;

  vec_muA = rep_vector(muA, Nspec);
  vec_muB1 = alphaB1 + betaB1 * Woody;
  phyloA = CparA*((phy-diag_matrix(diagonal(phy)))*lambdaA + diag_matrix(diagonal(phy)));
  phyloB1 = CparB1*((phy-diag_matrix(diagonal(phy)))*lambdaB1 + diag_matrix(diagonal(phy)));
  AL = cholesky_decompose(phyloA);
  A = vec_muA + AL * Apom;
  B1L = cholesky_decompose(phyloB1);
  B1 = vec_muB1 + B1L * B1pom;
  {
  int pom;
  pom = 1;
    for (i in 1:Nspec){
      pommu[pom:(pom+DruhP[i]-1)] = A[i] + B1[i] * Treat[pom:(pom+DruhP[i]-1)];
      pom  += DruhP[i];
    }
  }

  //priors
  B1 ~ cauchy(0,5);
  alphaB1 ~ cauchy(0,5);
  betaB1 ~ cauchy(0,5);
  CparA ~ cauchy(0,5);
  CparB1 ~ cauchy(0,5);

  FinSize ~ binomial_logit(1, pommu);
  Apom ~ std_normal();
  B1pom ~ std_normal();
}
"

