from dask_mpi import initialize
from dask.distributed import Client
import numpy as np

class PCA:
    def __init__(self, analyze, sele, client, frame_start=0, frame_stop=None, explain=0.95, calc=True):
        self.a = analyze
        self.sele = sele
        self.client = client
        self.frame_start = frame_start
        self.frame_stop = frame_stop
        self.explain = explain
        
        self.nframes = self.a.nframes
        self.natoms = len(self.a.mpt.select(sele))

        if calc:
            print("***Starting***")
            self.get_mean()
            print("***Mean calculated***")
            self.get_cov()
            print("***Covariance matrix calculated***")
            self.svd()
            print("***SVD performed***")
            self.transform()
            print("***Coordinates transformed***")
        else:
            self.mean = None
            self.cov = None
            self.variance = None
            self.cumulated_variance = None
            self.dim = None
            self.coords = None

    def get_mean(self):
        calculated = self.a.selectMany(self.sele, self.client, lambda x: x.positions.ravel(), self.frame_start, self.frame_stop)
        self.mean = sum(calculated)/self.nframes # add positions for all timesteps and divide by tot. no.

    
    def get_cov(self):
        mean = self.mean # can't pickle self.vars
        def cov_calc(x): # function to parallelize
            y = x.positions.ravel() - mean # diff between positions and mean, at every time step
            return np.dot(y[:, np.newaxis], y[:, np.newaxis].T) # dot product, i.e., (X-M)*(Y-M)
    
        cov_size = self.natoms*3
        self.cov = np.zeros((cov_size, cov_size))

        calculated = self.a.selectMany(self.sele, self.client, cov_calc, self.frame_start, self.frame_stop)
        for i in calculated: self.cov += i # add dot products for all timesteps
        self.cov /= self.nframes - 1 # divide by (tot. no. - 1)


    def svd(self):
        e_vals, e_vects = np.linalg.eig(self.cov)

        sort_idx = np.argsort(e_vals)[::-1]
        self.variance = e_vals[sort_idx]
        self.p_components = e_vects[:, sort_idx] 
        n = len(self.variance)
        self.cumulated_variance = (np.cumsum(self.variance)/np.sum(self.variance))[:n]

        self.dim = np.where(self.cumulated_variance > self.explain)[0][0]

    def transform(self):
        mean = self.mean
        p_comp = self.p_components
        dim = self.dim
        
        def proj_calc(x):
            xyz = x.positions.ravel() - mean
            return np.dot(xyz, p_comp[:, :dim])

        self.coords = self.a.selectMany(self.sele, self.client, proj_calc, self.frame_start, self.frame_stop)