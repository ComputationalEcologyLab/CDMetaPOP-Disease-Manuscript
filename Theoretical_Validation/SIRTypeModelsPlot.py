import numpy as np
import sys,pdb
from scipy.integrate import odeint
import matplotlib.pyplot as plt

runthis = sys.argv[1]

def SIR_model(y, t, N, beta, gamma,vital):

    S, I, R = y

    dSdt = vital * N -beta * S * I / N  - vital * S

    dIdt = beta * S * I / N - gamma * I - vital * I

    dRdt = gamma * I - vital * R

    return [dSdt, dIdt, dRdt]

def SIRS_model(y, t, N, beta, gamma,phi,vital):

    S, I, R = y

    dSdt = vital * N -beta * S * I / N + phi * R  - vital * S

    dIdt = beta * S * I / N - gamma * I - vital * I

    dRdt = gamma * I - phi * R - vital * R

    return [dSdt, dIdt, dRdt]
	
def SIRD_model(y, t, N, beta, gamma,mu_d,vital):
	S, I, R, D = y
	dSdt = vital * (S+R) -beta * S * I / N - vital * S
	dIdt =  vital * I + beta * S * I / N - gamma * I - mu_d * I - vital * I
	dRdt = gamma * I - vital * R  
	#dDdt = mu_d * I 
	#dDdt = mu_d * I - mu_d * D
	dDdt = mu_d * I - mu_d * D
	
	return [dSdt, dIdt, dRdt, dDdt]
	
def SEIR_model(y, t, N, beta, epsilon, gamma, vital):

	S, E, I, R = y
	dSdt = vital*N -beta * S * I / N - vital * S

	dEdt = beta * S * I / N - epsilon * E - vital * E
	dIdt = epsilon * E - gamma * I - vital * I

	dRdt = gamma * I - vital * R

	return [dSdt, dEdt, dIdt, dRdt]
	
# SIDP model differential equations
def SIDP_model(y, t, N, beta_d, beta_i, mu_d, alpha, mu,phi,alpha_d):
	S, I, D, P = y

	dSdt = -beta_d * S * I / N - beta_i * S * P + phi *I 
	dIdt = beta_d * S * I / N + beta_i * S * P  - mu_d * I- phi * I
	dDdt = mu_d * I
		
	dPdt = alpha * I - mu * P + alpha_d * D
	return [dSdt, dIdt, dDdt, dPdt]
# SIRP model differential equations
def SIRP_model(y, t, N, beta_d, beta_i, gamma, alpha, mu):
	S, I, R, P = y

	dSdt = -beta_d * S * I / N - beta_i * S * P 
	dIdt = beta_d * S * I / N + beta_i * S * P  - gamma * I
	dRdt = gamma * I
	
	dPdt = alpha * I - mu * P
	return [dSdt, dIdt, dRdt, dPdt]
	
# SIRP model differential equations with P / N
def SIRPdivN_model(y, t, N, beta_d, beta_i, gamma, alpha, mu):
	S, I, R, P = y

	dSdt = -beta_d * S * I / N - beta_i * S * P / N 
	dIdt = beta_d * S * I / N + beta_i * S * P / N  - gamma * I
	dRdt = gamma * I
	
	dPdt = alpha * I - mu * P
	
	return [dSdt, dIdt, dRdt, dPdt]
	


# Time points
t = np.linspace(0, 50, 100)

if runthis == 'SIRS':
	N = 100
	beta = 0.5
	gamma = 0.2
	phi = 0.2
	vital = 0.00
	y0 = [N-10, 10, 0]
	ret = odeint(SIRS_model, y0, t, args=(N, beta, gamma,phi,vital))
	plt.plot(t, ret[:, 0], label='Susceptible')
	plt.plot(t, ret[:, 1], label='Infected')
	plt.plot(t, ret[:, 2], label='Recovered')
	

elif runthis == 'SIRD':
	N = 100
	beta = 0.1
	gamma = .1
	mu_d = 0
	vital = 0
	y0 = [N-10, 10, 0,0]
	ret = odeint(SIRD_model, y0, t, args=(N, beta, gamma,mu_d,vital))
	plt.plot(t, ret[:, 0], label='Susceptible')
	plt.plot(t, ret[:, 1], label='Infected')
	plt.plot(t, ret[:, 2], label='Recovered')
	plt.plot(t, ret[:, 3], label='Death')
	
elif runthis == 'SEIR':
	N = 100
	beta = 0.5
	epsilon = 0.5
	gamma = 0.2
	vital = 0.00
	y0 = [N-10, 0, 10, 0]
	ret = odeint(SEIR_model, y0, t, args=(N, beta, epsilon, gamma,vital))
	plt.plot(t, ret[:, 0], label='Susceptible')
	plt.plot(t, ret[:, 1], label='Exposed')
	plt.plot(t, ret[:, 2], label='Infected')
	plt.plot(t, ret[:, 3], label='Recovered')

elif runthis == 'SIDP':
	beta_d = 0  # Direct transmission rate
	beta_i = 1  # Indirect transmission rate via environment
	gamma = 0  # Host recovery rate
	alpha = 1  # Shedding rate of pathogen by infected hosts
	mu = 0  # Pathogen decay rate
	alpha_d = 1 # Rate death decay
	phi = 0 # I to S rate
	mu_d = 0 # I to D

	# Initial conditions
	S0 = 100  # Initial number of susceptible individuals
	I0 = 0  # Initial number of infected individuals
	R0 = 0  # Initial number of removed individuals
	P0 = 1  # Initial pathogen load in the environment - between 0 and 1
	N = S0 + I0 + P0 # Total population
	y0 = [S0, I0, R0, P0]
	
	ret = odeint(SIDP_model, y0, t, args=(N, beta_d, beta_i, mu_d, alpha, mu,phi,alpha_d))

	# Get solutions
	S, I, D, P = ret.T

	plt.figure(figsize=(12, 8))
	plt.plot(t, S, label='Susceptible (S)')
	plt.plot(t, I, label='Infected (I)')
	plt.plot(t, D, label='Removed (D)')
	plt.plot(t, P, label='Pathogen Load (P)', linestyle='--')
	
elif runthis == 'SIRP':
	beta_d = 0  # Direct transmission rate
	beta_i = 1  # Indirect transmission rate via environment
	gamma = 0  # Host recovery rate
	alpha = 1  # Shedding rate of pathogen by infected hosts
	mu = 0.15  # Pathogen decay rate

	# Initial conditions
	S0 = 100  # Initial number of susceptible individuals
	I0 = 0  # Initial number of infected individuals
	R0 = 0  # Initial number of removed individuals
	P0 = 1  # Initial pathogen load in the environment - between 0 and 1
	N = S0 + I0 + P0 # Total population
	y0 = [S0, I0, R0, P0]
	
	ret = odeint(SIRP_model, y0, t, args=(N, beta_d, beta_i, gamma, alpha, mu))

	# Get solutions
	S, I, R, P = ret.T

	plt.figure(figsize=(12, 8))
	plt.plot(t, S, label='Susceptible (S)')
	plt.plot(t, I, label='Infected (I)')
	plt.plot(t, R, label='Removed (R)')
	plt.plot(t, P, label='Pathogen Load (P)', linestyle='--')

elif runthis == 'SIRP-N':
	beta_d = 0.25  # Direct transmission rate
	beta_i = 0.3  # Indirect transmission rate via environment
	gamma = 0.05  # Host recovery rate
	alpha = 0.1  # Shedding rate of pathogen by infected hosts
	mu = 0.15  # Pathogen decay rate

	# Initial conditions
	S0 = 90  # Initial number of susceptible individuals
	I0 = 10  # Initial number of infected individuals
	R0 = 0  # Initial number of removed individuals
	P0 = 10  # Initial pathogen load in the environment
	N = S0 + I0 + P0 # Total population
	y0 = [S0, I0, R0, P0]
	
	ret = odeint(SIRPdivN_model, y0, t, args=(N, beta_d, beta_i, gamma, alpha, mu))

	# Get solutions
	S, I, R, P = ret.T

	plt.figure(figsize=(12, 8))
	plt.plot(t, S, label='Susceptible (S)')
	plt.plot(t, I, label='Infected (I)')
	plt.plot(t, R, label='Removed (R)')
	plt.plot(t, P, label='Pathogen Load (P)', linestyle='--')

# Solve the differential equations
elif runthis == 'SIR':
	N = 100
	beta = 0.8
	gamma = 0.1
	vital = 0.00
	y0 = [N-10, 10, 0]
	ret = odeint(SIR_model, y0, t, args=(N, beta, gamma,vital))
	plt.plot(t, ret[:, 0], label='Susceptible (S)')
	plt.plot(t, ret[:, 1], label='Infected (I)')
	plt.plot(t, ret[:, 2], label='Recovered (R)')


# Plot the results
#plt.xlabel('Time')
#plt.ylabel('Population Size')
#plt.legend()
plt.show()
