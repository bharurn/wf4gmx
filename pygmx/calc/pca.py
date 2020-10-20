from dask_mpi import initialize
from dask.distributed import Client
import numpy as np

class PCA:
    def __init__(self, uni, sele, explain=0.95, calc=True):
        ### Using Dask with MPI 
        initialize() # init dask mpi
        self.client = Client()  # MPI compatible client
        ###
        
        self.u = uni
        self.sele = sele
        self.explain = explain

        if calc:
            self.calc_mean()
            self.calc_cov()
            self.svd()
            self.transform()
        else:
            self.mean = None
            self.cov = None
            self.variance = None
            self.cumulated_variance = None
            self.dim = None
            self.coords = None

    def calc_mean():
        futures = []

        for i in range(len(self.u.trajectory)):
            self.u.trajectory[i]
            # parallelize ravel operation on positions array for each timestep
            future = client.submit(lambda pos: pos.ravel(), self.u.select_atoms(self.sele).atoms.positions)
            futures.append(future)

        self.mean = sum(client.gather(futures))/len(self.u.trajectory) # add positions for all timesteps and divide by tot. no.

    
    def calc_cov():
        def calc_cov(pos): # function to parallelize
            x = pos.ravel() - mean # diff between positions and mean, at every time step
            return np.dot(x[:, np.newaxis], x[:, np.newaxis].T) # dot product, i.e., (X-M)*(Y-M)
    
        cov_size = len(self.u.select_atoms(receptor).atoms)*3
        self.cov = np.zeros((cov_size, cov_size))

        # cannot parallelize the whole trajectory in one shot
        # insufficient memeory for gathering all numpy arrays 
        # so we split it in chunks of 1000

        interval = 1000
        for k in range(0, len(self.u.trajectory), interval):
            end = k+interval
            if end>len(self.u.trajectory): end = len(self.u.trajectory)
    
            futures = []
    
            print(f"Calculating for range {(k, end)}..")
    
            for i in range(k, end):
                self.u.trajectory[i]
                future = client.submit(calc_cov, self.u.select_atoms(self.sele).atoms.positions)
                futures.append(future)

            for i in client.gather(futures): self.cov += i # add dot products for all timesteps

        self.cov /= len(self.u.trajectory) - 1 # divide by (tot. no. - 1)


    def svd():
        e_vals, e_vects = np.linalg.eig(cov)

        sort_idx = np.argsort(e_vals)[::-1]
        self.variance = e_vals[sort_idx]
        self.p_components = e_vects[:, sort_idx] 
        n = len(variance)
        self.cumulated_variance = (np.cumsum(variance)/np.sum(variance))[:n]

        self.dim = np.where(cumulated_variance > self.explain)[0][0]

    def transform():
        def proj_calc(pos):
            xyz = pos.ravel() - mean
            return np.dot(xyz, p_components[:, :dim])

        dots = []

        for i in range(len(self.u.trajectory)):
            self.u.trajectory[i]
            future = client.submit(proj_calc, self.u.select_atoms(self.sele).atoms.positions)
            dots.append(future)
        
        self.coords = client.gather(dot)