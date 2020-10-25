import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

class Cluster:
    
    def __init__(self, data, method=None):
        self.data = data
        self.method = method
        
    def traj_cluster(self, plot=True):
        y_km = self.method.fit_predict(self.data)

        if plot:
            fig, ax = plt.subplots()
            ax.scatter(self.data[:,0], self.data[:,1], c= self.method.labels_.astype(float), s=25, alpha=0.5)
            if hasattr(self.method, 'cluster_centers_'):
                ax.scatter(self.method.cluster_centers_[:, 0], self.method.cluster_centers_[:, 1], marker='*',\
                       c='red', edgecolor='black',label='centroids')
    
    def get_centers(self, cutoff=3):
        cen = self.method.cluster_centers_
        pos = (cen[:, np.newaxis] - self.data)
    
        reshaped = pos.reshape(-1, cen.shape[1])
        dist = np.sum(np.abs(reshaped)**2,axis=1)**(1./2)
        dist = dist.reshape(pos.shape[0], pos.shape[1])
        ids = np.where(dist<cutoff)
        vals = dist[ids]

        frames_ = pd.DataFrame([*ids, vals]).T.groupby(0).min()
        return frames_
    
    def heirarchy(self, max_d=2.5, trunc=60):
        from scipy.cluster import hierarchy
        
        linked = hierarchy.linkage(self.data, 'single')

        labelList = range(0, self.data.shape[0])

        plt.figure(figsize=(10, 7))
        hierarchy.dendrogram(linked,
                    orientation='top',
                    truncate_mode='lastp',
                    p=trunc,
                    labels=labelList,
                    distance_sort='descending',
                    show_leaf_counts=True)
        plt.axhline(y=max_d, c='k')
        plt.show()
    
        clus = hierarchy.fcluster(linked, max_d, criterion='distance')
        cluster = pd.DataFrame([clus, labelList]).T
        cluster.columns=['cluster', 'frame']
        return cluster
        
    def optimalK(self, nrefs=3, maxClusters=15):
        """
        Calculates KMeans optimal K using Gap Statistic from Tibshirani, Walther, Hastie
        Params:
            nrefs: number of sample reference datasets to create
            maxClusters: Maximum number of clusters to test for
        Returns: optimalK
        """
        from sklearn.cluster import KMeans
        gaps = np.zeros((len(range(1, maxClusters)),))
        resultsdf = pd.DataFrame({'clusterCount':[], 'gap':[]})
        for gap_index, k in enumerate(range(1, maxClusters)):

            # Holder for reference dispersion results
            refDisps = np.zeros(nrefs)

            # For n references, generate random sample and perform kmeans getting resulting dispersion of each loop
            for i in range(nrefs):

                # Create new random reference set
                randomReference = np.random.random_sample(size=self.data.shape)

                # Fit to it
                km = KMeans(k)
                km.fit(randomReference)

                refDisp = km.inertia_
                refDisps[i] = refDisp

            # Fit cluster to original data and create dispersion
            km = KMeans(k)
            km.fit(self.data)

            origDisp = km.inertia_

            # Calculate gap statistic
            gap = np.log(np.mean(refDisps)) - np.log(origDisp)

            # Assign this loop's gap statistic to gaps
            gaps[gap_index] = gap

            resultsdf = resultsdf.append({'clusterCount':k, 'gap':gap}, ignore_index=True)
        
        self.gapdf = resultsdf
        
        # Plus 1 because index of 0 means 1 cluster is optimal, index 2 = 3 clusters are optimal
        return gaps.argmax() + 1