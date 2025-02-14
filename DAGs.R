install.packages("dagitty")
install.packages("ggdag")

library(dagitty)
library(ggdag)


# OG GDM 
library(DiagrammeR)

# Define the DAG structure
dag <- "
digraph BayesianModel {
  # Nodes
  f_t [label='f(t)', shape=circle]
  g_t_d [label='g(t,d)', shape=circle]
  lambda_t [label='Lambda_t', shape=circle]
  theta [label='Theta', shape=circle]
  y_t [label='Y_t', shape=circle, style=filled, fillcolor=lightgrey]
  z_t_d [label='Z_t_d', shape=circle, style=filled, fillcolor=lightblue]
  z_t_less_d [label='Z_t_<d', shape=circle]
  nu_t_d [label='Nu_t_d', shape=circle]
  phi_t_d [label='Phi_t_d', shape=circle]

  # Relationships
  f_t -> lambda_t
  g_t_d -> nu_t_d
  lambda_t -> y_t
  theta -> y_t
  y_t -> z_t_d [label='via N_t_d']
  z_t_less_d -> z_t_d
  nu_t_d -> z_t_d
  phi_t_d -> z_t_d
}
"

# Render the DAG
grViz(dag)

# NESTED GDM - PARAMS
library(DiagrammeR)

dag <- "
digraph SARI_Model {
  # Nodes
  zeta0_s [label='Zeta0_s', shape=circle]
  zeta_t_s [label='Zeta1_t_s', shape=circle]
  lambda_t_s [label='Lambda_t_s', shape=circle]
  theta_s [label='Theta_s', shape=circle]
  y_t_s [label='Y_t_s', shape=circle, style=filled, fillcolor=lightgrey]
  x_t_s [label='X_t_s', shape=circle, style=filled, fillcolor=lightblue]
  mu_t_s [label='Mu_t_s', shape=circle]
  chi_s [label='Chi_s', shape=circle]
  pi_t_s [label='Pi_t_s', shape=circle]
  omega_s [label='Omega_s', shape=circle]
  c_t_s [label='C_t_s', shape=circle, style=filled, fillcolor=lightblue]
  upsilon_s [label='Upsilon_s', shape=circle]
  nu_t_d_s [label='Nu_t_d_s', shape=circle]
  phi_s_d [label='Phi_s_d', shape=circle]
  psi_d_s [label='Psi_d_s', shape=circle]
  eta_t_s [label='Eta_t_s', shape=circle]
  z_t_d_s [label='Z_t_d_s', shape=circle, style=filled, fillcolor=lightblue]
  probit_p_t_d_s [label='Probit_p_t_d_s', shape=circle]
  beta0_s [label='Beta0_s', shape=circle]
  beta_t_s [label='Beta1_t_s', shape=circle]
  delta_s [label='Delta_s', shape=circle]

  # Relationships
  zeta0_s -> lambda_t_s
  zeta_t_s -> lambda_t_s
  lambda_t_s -> y_t_s
  theta_s -> y_t_s
  y_t_s -> x_t_s
  mu_t_s -> x_t_s
  chi_s -> x_t_s
  y_t_s -> z_t_d_s
  nu_t_d_s -> z_t_d_s
  phi_s_d -> z_t_d_s
  probit_p_t_d_s -> nu_t_d_s
  psi_d_s -> probit_p_t_d_s
  eta_t_s -> probit_p_t_d_s
  x_t_s -> c_t_s
  pi_t_s -> c_t_s
  upsilon_s -> c_t_s
  omega_s -> pi_t_s
  beta0_s -> mu_t_s
  beta_t_s -> mu_t_s
  zeta_t_s -> mu_t_s
  delta_s -> mu_t_s
}
"

# Render the DAG
grViz(dag)

# CASE-LOAD EFFECT DAG

install.packages("DiagrammeR")

library(DiagrammeR)

# Define the DAG
dag <- "
digraph BayesianModel {
  # Nodes
  population_s [label='Population_s', shape=circle]
  iota_s [label='Iota_s', shape=circle]
  zeta_t_s [label='Zeta_t_s', shape=circle]
  xi_weeks_t_s [label='Xi_weeks_t_s', shape=circle]
  theta_s [label='Theta_s', shape=circle]
  lambda_t_s [label='Lambda_t_s', shape=circle]
  y_t_s [label='Y_t_s', shape=circle, style=filled, fillcolor=lightgrey]
  S_t_d_s [label='S_t_d_s', shape=circle]
  nu_t_d_s [label='Nu_t_d_s', shape=circle]
  z_t_d_s [label='Z_t_d_s', shape=circle]
  phi_d_s [label='Phi_d_s', shape=circle]
  psi_d_s [label='Psi_d_s', shape=circle]
  eta_t_s [label='Eta_t_s', shape=circle]
  delta_s [label='Delta_s', shape=circle]

  # Edges
  population_s -> lambda_t_s
  iota_s -> lambda_t_s
  zeta_t_s -> lambda_t_s
  xi_weeks_t_s -> lambda_t_s
  lambda_t_s -> y_t_s
  theta_s -> y_t_s
  y_t_s -> S_t_d_s
  y_t_s -> z_t_d_s
  psi_d_s -> S_t_d_s
  eta_t_s -> S_t_d_s
  delta_s -> S_t_d_s
  S_t_d_s -> nu_t_d_s
  nu_t_d_s -> z_t_d_s
  phi_d_s -> z_t_d_s
}
"

# Render the DAG
grViz(dag)
