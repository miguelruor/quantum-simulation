import numpy as np
import networkx as nx
from scipy.sparse.linalg import expm_multiply
from scipy.stats import expon
from ibmcloudant.cloudant_v1 import CloudantV1, Document
from ibm_cloud_sdk_core.authenticators import IAMAuthenticator
import os
import time as timing
from datetime import datetime
import uuid

def giveMeHamiltonian(G, gamma, typeMatrix="laplacian"):
  # typeMatrix could be "adjacency" to use the adjacency matrix based Hamiltoninan 
  # or "laplacian" (default) to use the Laplacian matrix based Hamiltonian
  # gamma is the mutation rate 

  if typeMatrix == "adjacency":
    A = nx.adjacency_matrix(G)
  else:
    A = nx.laplacian_matrix(G)

  H = -gamma* A

  return H

def canonical_vector(i, n):
  # i-esimo vector canonico en C^n
  ei = np.zeros(n)
  ei[i] = 1.0
  return ei

def measurement(state, basis):
  # state is a vector in C^N

  qprobs = [abs(state[v])**2 for v in basis]
  collapse = np.random.choice(basis, p=qprobs)

  return collapse

def simulation_qw(gspace, gspace_name, phenotypes, initial_genotype, max_simulation_time, 
  measurement_rate, gamma, max_execution_time, first_checkpoint, waiting_for_checkpoint):

  start = timing.time()
  date = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

  # initial data

  H = giveMeHamiltonian(gspace, gamma)
  M = len(gspace.nodes)

  tau = {} # tau (hitting time) for every phenotype: time it takes quantum walk to find given phenotype
  N = {} # number of measurements the quantum walk takes to find a genotype with a new phenotype
  mutations = {} # mutations the quantum walk took to find a genotype with a given phenotype

  no_measurement = 0
  time = 0
  total_mutations = 0
  #mutation_time = 0

  checkpoint = 0

  _id = uuid.uuid1()

  # database connection
  authenticator = IAMAuthenticator(os.environ['CLOUDANT_APIKEY'])

  cloudant = CloudantV1(authenticator=authenticator)
  cloudant.set_service_url(os.environ['CLOUDANT_URL'])

  try:
    client = cloudant.new_instance()
  except:
    print("Could not create a cloudant client")

  simulation: Document = Document()
  simulation.initial_gen_index = initial_genotype
  simulation._id = _id
  simulation.initial_gen = gspace.nodes[initial_genotype]['sequence']
  simulation.initial_phen = gspace.nodes[initial_genotype]['phenotypeName'][0]
  simulation.measurement_rate = measurement_rate
  simulation.transition_rate = gamma
  simulation.max_simulation_time = max_simulation_time
  simulation.datetime = date

  #initialization
  for phen in phenotypes:
    tau[phen] = -1
    N[phen] = -1
    mutations[phen] = -1
 
  actual_state = initial_genotype 
  phenotypes_actual_state = gspace.nodes[actual_state]['phenotypeName'] # phenotypes of actual state

  # Start of simulation
  while time < max_simulation_time and timing.time()-start <= max_execution_time:
    for phen in phenotypes_actual_state:
      if tau[phen] < 0: # update hitting times of novel phenotypes
        tau[phen] = time
        N[phen] = no_measurement
        mutations[phen] = total_mutations 
      
    T = expon.rvs(scale=measurement_rate, size=1)[0] # time between measurements
    time += T
    
    # evolving quantum walk
    actual_state_vec = canonical_vector(actual_state, M) # vector representing genotype actual_state
    actual_state_vec = expm_multiply(-1j*T*H, actual_state_vec) # evolve quantum walk until time T with actual_state as initial state
    measurement_result = measurement(actual_state_vec, gspace.nodes)
    no_measurement += 1

    if measurement_result != actual_state:
      total_mutations += 1

    actual_state = measurement_result

    # phenotypes of actual state
    phenotypes_actual_state = gspace.nodes[actual_state]['phenotypeName']

    if time >= first_checkpoint + checkpoint*waiting_for_checkpoint:
      checkpoint += 1

      end = timing.time()

      simulation.total_measurements = no_measurement
      simulation.total_mutations = total_mutations
      simulation.computing_time = end-start
      simulation.simulation_time = time

      for phen in phenotypes:
        setattr(simulation, 'tau_'+phen, tau[phen] if tau[phen]>=0.0 else time)
        setattr(simulation, 'N_'+phen, N[phen])
        setattr(simulation, 'mutations_'+phen, mutations[phen])

      try:
        client.post_document(
          db="simulations-"+gspace_name,
          document=simulation
        )
        print(f"Wrote in Cloudant successfully: {checkpoint+1} checkpoint")

      except:
        print("Unexpected error")

  # end of simulation
  
  # writing results of simulation  
  #evolution_paths.to_csv(url_evolution_paths)
  #simulation_results.to_csv(url_simulation_results)