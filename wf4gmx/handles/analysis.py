from .base import BaseHandle
from .._global import _Global as _global
from ..utils.logger import LogString
import re
import pandas as pd
import numpy as np
from mimicpy import Mpt
from ..io.trr import read_trr, get_trr_frames, TRRReader

class Analyze(BaseHandle):
    def __init__(self, mpt_file=None, trr_files=None, status_file=None, status=None):
        super().__init__(status_file, status)
        self.__getMpt(mpt_file)
        
        if trr_files is None:
            trr_files = self.gethistory('trr')[::-1]
        elif not isinstance(trr_files, list):
            trr_files = [trr_files]
            
        each_trr_frames = [get_trr_frames(trr) for trr in trr_files]
        self.trrs = list(zip(trr_files, each_trr_frames))
        
        self.nframes = sum(each_trr_frames)
        
    def __getMpt(self, file):
        if file is None:
            file = self.getcurrent('mpt', exp=False)
            if not file:
                file = self.getcurrent('top')
        
        self.mpt = Mpt.from_file(file)
        
    @staticmethod
    def get_frame(sele, trr_files, frame):
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

        return SelectedFrame({k:v for k,v in trr_data.items() if k != 'x'}, sele, 
                                  trr_data['x'][ids])
        
    def select(self, selection, frame=0):
        sele = self.mpt.select(selection)
        return Analyze.get_frame(sele, self.trrs, frame)
    
    def selectMany(self, selection, client, calc=None, frame_start=0, frame_stop=None, asFutures=False):
        
        sele = self.mpt.select(selection)
        
        if frame_stop is None: frame_stop = self.nframes
        
        if calc:
            do = lambda sele, trr_file, frame: calc(Analyze.get_frame(sele, trr_file, frame))
        else:
            do = lambda sele, trr_file, frame: Analyze.get_frame(sele, trr_file, frame)
        
        if asFutures:
            return [client.submit( do, client.scatter(sele), client.scatter(self.trrs), i ) for i in range(frame_start, frame_stop)]
        else:
            return [client.submit( do, client.scatter(sele), client.scatter(self.trrs), i ).result() for i in range(frame_start, frame_stop)]
        
            
    @staticmethod
    def distribute(u, start, end, pkl_file, loader):
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
    
    def __init__(self, header, df, positions):
        for k, v in header.items():
            setattr(self, k, v)
        self.natoms = len(df)
        
        self.__repr_list = []
        
        for column in df:
            self.__repr_list.append(column)
            setattr(self, column, df[column])
        
        self.positions = positions
    
    def __repr__(self):
        dct = {}
        for i in self.__repr_list:
            dct[i] = getattr(self, i)
        dct['x'] = self.positions[:, 0]
        dct['y'] = self.positions[:, 1]
        dct['z'] = self.positions[:, 2]
        return repr(pd.DataFrame(dct))
