import numpy as np
import pandas as pd
    
class TrajClust:
    
    def __init__(self):
        self.df = None
    
    def pca(self, u, sele):
        import MDAnalysis.analysis.pca as pca
        uni_pca = pca.PCA(u, select=sele)
        uni_pca.run()

        n_pcs = np.where(uni_pca.cumulated_variance > 0.95)[0][0]
        atomgroup = u.select_atoms(sele)
        pca_space = uni_pca.transform(atomgroup, n_components=n_pcs)

        self.df = pd.DataFrame(pca_space)
        return self.df

    def cluster(self, method, plot=True):
        y_km = method.fit_predict(self.df)

        if plot:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots()
            ax.scatter(self.df[0], self.df[1], c= method.labels_.astype(float), s=25, alpha=0.5)
            if hasattr(method, 'cluster_centers_'):
                ax.scatter(method.cluster_centers_[:, 0], method.cluster_centers_[:, 1], marker='*',\
                       c='red', edgecolor='black',label='centroids')
    
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
                randomReference = np.random.random_sample(size=self.df.shape)

                # Fit to it
                km = KMeans(k)
                km.fit(randomReference)

                refDisp = km.inertia_
                refDisps[i] = refDisp

            # Fit cluster to original data and create dispersion
            km = KMeans(k)
            km.fit(self.df)

            origDisp = km.inertia_

            # Calculate gap statistic
            gap = np.log(np.mean(refDisps)) - np.log(origDisp)

            # Assign this loop's gap statistic to gaps
            gaps[gap_index] = gap

            resultsdf = resultsdf.append({'clusterCount':k, 'gap':gap}, ignore_index=True)
        
        self.gapdf = resultsdf
        
        # Plus 1 because index of 0 means 1 cluster is optimal, index 2 = 3 clusters are optimal
        return gaps.argmax() + 1