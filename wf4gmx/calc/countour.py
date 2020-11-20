from collections import defaultdict
import numpy as np

class CountourMap:
    
    def __init__(self, selectedframe, axis=[0,2]):
        self.data = selectedframe.positions
        self.ax1 = axis[0]
        self.ax2 = axis[1]
        self.data_x, self.data_y = self.data[:, self.ax1], self.data[:, self.ax2]
        self.mesh_x = None
        self.mesh_y = None
        
    def __get_box(self, i):
        xs = [(np.where(self.data_x[i] < self.mesh_x))[0][0], (np.where(self.data_x[i] > self.mesh_x))[0][-1]]
        xs.sort()
        
        ys = [(np.where(self.data_y[i] < self.mesh_y))[0][0], (np.where(self.data_y[i] > self.mesh_y))[0][-1]]
        ys.sort()
        
        tup_str = lambda x: f"{x[0]}-{x[1]}"
        return (tup_str(xs), tup_str(ys))
    
    def __tup2pnt(self, tup, a):
        x1, x2 = tup.split('-')
        return (a[int(x1)]+a[int(x2)])/2
    
    def get_mesh(self, num=10, padding=True):
        max_, min_ = np.max(self.data, axis=0), np.min(self.data, axis=0)
        ax1, ax2 = self.ax1, self.ax2
        if padding:
            pad1 = (max_[ax1]-min_[ax1])/10
            pad2 = (max_[ax2]-min_[ax2])/10
        else:
            pad1 = pad2 = 0
        
        x, y = np.meshgrid(np.linspace(min_[ax1]-pad1, max_[ax1]+pad1, num=10), np.linspace(min_[ax2]-pad2, max_[ax2]+pad2, num=num))
        self.mesh_x = x.flatten()
        self.mesh_y = y.flatten()
    
    def get_countour(self):
        ax1, ax2 = self.ax1, self.ax2
        
        if self.mesh_x is None or self.mesh_y is None:
            self.get_mesh()
        
        dct = defaultdict(int)
        for i in range(self.data.shape[0]):
            box = self.__get_box(i)
            dct[box] += 1

        pointsx, pointsz = list(zip(*dct.keys()))
        feature_x = np.unique(np.array(pointsx))
        feature_y = np.unique(np.array(pointsz))
        
        tup2x = np.vectorize(lambda tup: self.__tup2pnt(tup, self.mesh_x))
        self.X = tup2x(feature_x)
        
        tup2y = np.vectorize(lambda tup: self.__tup2pnt(tup, self.mesh_y))
        self.Y = tup2y(feature_y)
        
        self.X, feature_x = zip(*sorted(zip(self.X, feature_x)))
        self.Y, feature_y = zip(*sorted(zip(self.Y, feature_y)))
        
        X_,Y_ = np.meshgrid(feature_x, feature_y)

        get_Z = np.vectorize(lambda x,y: dct[(x, y)])
        self.Z = get_Z(X_, Y_)
        