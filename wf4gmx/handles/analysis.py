from .base import BaseHandle
import numpy as np
import pandas as pd
from mimicpy import Mpt
from ..io.trr import read_trr, get_trr_frames, TRRReader
from os import path
import pickle


class Analyze(BaseHandle):
    def __init__(self, mpt_file=None, trr_files=None, client=None, status_file=None, status=None):
        super().__init__(status_file, status)
        self.__getMpt(mpt_file)
        
        if trr_files is None:
            trr_files = self.gethistory('trr')[::-1]
        elif not isinstance(trr_files, list):
            trr_files = [trr_files]
            
        each_trr_frames = [get_trr_frames(trr) for trr in trr_files]
        self.trrs = list(zip(trr_files, each_trr_frames))
        
        self.nframes = sum(each_trr_frames)
        
        self.client = client
        
    def __getMpt(self, file):
        if file is None:
            file = self.getcurrent('mpt', exp=False)
            if not file:
                file = self.getcurrent('top')
        
        self.mpt = Mpt.from_file(file)
        
    @staticmethod
    def get_frame(sele, trr_files, frame, around=None):
        if isinstance(sele, list) and len(sele) > 1:
            return [Analyze.get_frame(s, trr_files, frame) for s in sele]
        elif isinstance(sele, list):
            sele = sele[0]
        
        ids = sele.index - 1

        frames_so_far = 0
        trr_file_names, nframes = list(zip(*trr_files))
        trr_data = None
        for i, f in enumerate(nframes):
            if frame < frames_so_far+f:
                trr_data = read_trr(trr_file_names[i], frame-frames_so_far)
                break
            else:
                frames_so_far += f
        
        if trr_data is None:
            raise Exception("Requested frame is greater than total number of frames.")
        
        x = trr_data['x'][ids] if 'x' in trr_data else np.array([])
        v = trr_data['v'][ids] if 'v' in trr_data else np.array([])
        f = trr_data['f'][ids] if 'f' in trr_data else np.array([])
        
        sele = SelectedFrame(trr_data['step'], trr_data['time'], sele, x, v, f)
        
        return sele
    
    @staticmethod
    def log_pickle(x, file_name):
        if file_name is None: return x
        
        xyz = np.array(x)
        
        try:
            if path.isfile(file_name):
                with open(file_name, 'rb') as f:
                    try:
                        old = pickle.load(f)
                        xyz = np.vstack((old, xyz))
                    except EOFError:
                        pass

            with open(file_name, 'wb') as f: pickle.dump(xyz, f)
        except NameError:
            pass
        
        return x
        
    def select(self, *selection, calc=None, frames=0, asFutures=False, **kwargs):
        if isinstance(selection, str):
            sele = self.mpt.select(selection)
        else:
            sele = [self.mpt.select(s) for s in selection]
        
        trr_file = self.trrs
        
        if calc:
            do = lambda frame, **kwargs: calc(Analyze.get_frame(sele, trr_file, frame), **kwargs)
        else:
            do = lambda frame: Analyze.get_frame(sele, trr_file, frame)
        
        if isinstance(frames, int):
            ret = do(frames)
        elif self.client is None:
            ret = [do(i) for i in frames]
        elif asFutures:
            ret = [self.client.submit(do, i, **kwargs) for i in frames]
        else:
            ret = [self.client.submit(do, i, **kwargs).result() for i in frames]
        
        return ret
     
    def write(self, file, *selection, calc=None, frame=0):
        if not isinstance(frame, int):
            raise Exception
            
        sele_frame = self.select(*selection, calc=calc, frames=frame)
        
        if not isinstance(sele_frame, SelectedFrame):
            raise Exception
            
        sele_frame.write(file)
            
    @staticmethod
    def mpi_distribute(u, start, end, pkl_file, loader):
        def wrapper(func):
        
            ###imports
            from mpi4py import MPI
            import pickle

            comm = MPI.COMM_WORLD
            rank = comm.Get_rank()
            size = comm.Get_size()
        
            ###
        
            ###load files on rank 0
            if rank == 0 and loader:
                from tqdm import tqdm
                pbar = tqdm(desc=f"Starting", bar_format='{l_bar}{bar}|')
            else:
                files = None
            ###
        
            ###calculate start and end frames for each rank
            if end == -1: end_ = len(u.trajectory)
            else: end_ = end
            
            perrank = (end_-start)//size
            
            if perrank == 0:
                raise Exception("The number of MPI ranks are too large for the trajectory! Reduce number of ranks or increase number of frames in the trajectory.")
            
            rank_start = start + rank*perrank
            rank_end = start + (rank+1)*perrank
            
            if not loader:    
                if rank == 0: print(f"{size} MPI ranks will each run on {perrank} frames of the trajectory")
                print(f"Rank {rank}: start at frame {rank_start}, end at frame {rank_end}")
            
            ###
        
            ###start calculating
            comm.Barrier()

            x = []
        
            if rank == 0 and loader:
                pbar.total=rank_end-rank_start
                pbar.desc = f"Calculating {rank_end-rank_start} frames on each rank"
        
            for i in range(rank_start, rank_end):
                u.trajectory[i]
                x.append((i,func(u)))
                if rank == 0 and loader: pbar.update(1)
    
            if rank == 0 and loader: pbar.close()
            ###
            ###gather results into rank 0
            
            arr_size = len(pickle.dumps(x, -1))*size
            if arr_size > 2000000000:
                raise Exception(f"Total size of result greater than 2 GB, cannot gather!")
                
            result = comm.gather(x, root=0)
            
            if rank == 0:
                flat = [item for sublist in result for item in sublist]
            
                ##calculate extra frames not equally divided among ranks
                l = len(flat) + start
                if l < end_:
                    if loader:
                        pbar = tqdm(total=end_-l, desc=f"Calculating {end_-l} extra frames on rank 0", bar_format='{l_bar}{bar}|')
                    else:
                        print(f"Calculating {end_-l} extra frames on rank")
 
                    for i in range(l, end_):
                        u.trajectory[i]
                        flat.append((i,func(u)))
                        if loader: pbar.update(1)
                
                ##
                x = [i[1] for i in sorted(flat, key=lambda x: x[0])] # flatten nested list
                
                if pkl_file != None: # pickle results    
                    print("Pickling..")
                    import pickle
                    pickle.dump(x, open(pkl_file, 'wb'))
                    print("Done..")
            ###
                
        return wrapper
    
class SelectedFrame:
    
    def __init__(self, step, time, df, x, v, f):
        self.step = step
        self.time = time
        self.natoms = len(df)
        self.ids = df.index.to_numpy()
        
        self.__repr_list = []
        
        for column in df:
            self.__repr_list.append(column)
            setattr(self, column, df[column].to_numpy())
        
        self.positions = x
        self.velocities = v
        self.forces = f
    
    def __getitem__(self, k):
        if isinstance(k, int): k = [k]
        dct = {}
        for i in self.__repr_list:
            dct[i] = getattr(self, i)[k]
        mpt_df = pd.DataFrame(dct, index=self.ids[k])
        
        try:
            pos = self.positions[k]
        except IndexError:
            pos = []
        
        try:
            v = self.velocities[k]
        except IndexError:
            v = []
        
        try:
            f = self.forces[k]
        except IndexError:
            f = []
            
        return SelectedFrame(self.step, self.time, mpt_df, pos, v, f)
    
    def __repr__(self):
        dct = {}
        for i in self.__repr_list:
            dct[i] = getattr(self, i)
        dct['x'] = self.positions[:, 0]
        dct['y'] = self.positions[:, 1]
        dct['z'] = self.positions[:, 2]
        return repr(pd.DataFrame(dct))
    
    def __add__(self, other):
        dct = {}
        for i in self.__repr_list:
            dct[i] = np.append(getattr(self, i), getattr(other, i))
        mpt_df = pd.DataFrame(dct, index=np.append(self.ids, other.ids))

        pos = np.vstack((self.positions, other.positions))
        v = np.vstack((self.velocities, other.velocities))
        f = np.vstack((self.forces, other.forces))
                              
        return SelectedFrame(self.step, self.time, mpt_df, pos, v, f)
    
    def get_around(self, sele, dist=3):
        cen = np.mean(self.positions, axis=0)
        ids = np.where( np.linalg.norm(sele.positions - cen, axis=1) <= dist)[0].tolist()
        return (self + sele[ids])

    
    def write(self, file):
        dct = {}
        for i in self.__repr_list:
            dct[i] = getattr(self, i)
        mpt_df = pd.DataFrame(dct)
        mpt_df['id'] = self.ids
        
        dct = {}
        dct['x'] = self.positions[:, 0]
        dct['y'] = self.positions[:, 1]
        dct['z'] = self.positions[:, 2]
        coords_df = pd.DataFrame(dct)
        coords_df['id'] = self.ids
        from mimicpy import CoordsIO
        CoordsIO(file, mode='w').write(mpt_df, coords_df)
