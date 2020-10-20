from .base import BaseHandle
from .._global import _Global as _global
from ..utils.logger import LogString
import re
import pandas as pd
import numpy as np

class Analyze(BaseHandle):
    def __init__(self, status_file=None, status=None):
        super().__init__(status_file, status)
        self.xvg = LogString()
        self.logger.add(xvg=self.xvg)
        self.__getUni()
        
    def __getUni(self):
        try:
            files = [self.getcurrent('tpr')]
            files += self.gethistory('trr')[::-1]
        except FileNotFoundError as e:
            self.u = None
        else:
            import MDAnalysis as mda
            self.u = mda.Universe(files[0], files[1:])
        
    def __xvg_df(self, cmd, prev_files):
        new_files = _global.host.ls()
        
        extra_files = [n for n in new_files if n not in prev_files]
        
        if len(extra_files) == 1 and extra_files[0].split('.')[1] == 'xvg':
            o = _global.host.read(extra_files[0])
            self.logger.write('xvg', f"==>XVG raw output from gmx {cmd}\n")
            self.logger.write('xvg', o)
            host.rm(extra_files[0])
            return Analyze.readxvg(o)
        else:
            return extra_files
    
    def gmx(self, cmd, **kwargs):
        
        prev_files = _global.host.ls()
        
        if 's' not in kwargs: kwargs['s'] = self.getcurrent('tpr')
        if 'f' not in kwargs: kwargs['f'] = self.gethistory('trr')
        
        _kwargs = kwargs.copy()
        del _kwargs['f']
        
        df = pd.DataFrame()
        
        for trr in kwargs['f']:
            super().gmx(cmd, f=trr, **_kwargs)
            if df.empty: df = self.__xvg_df(cmd, prev_files)
            else: df = df.append(self.__xvg_df(cmd, prev_files))
        
        return df.sort_index(axis = 0)
    
    @staticmethod
    def readxvg(xvg, readlabel=True):
        title = ''
        x = []
        y = []
    
        if readlabel:
            reg = re.compile(r"@\s*title\s*\"(.*?)\"", re.MULTILINE).findall(xvg)
            if reg == []: title = 'Title'
            else: title = reg[0]
            
            reg = re.compile(r"@\s*xaxis\s*label\s*\"(.*?)\"", re.MULTILINE).findall(xvg)
            if reg == []: xlabel = 'X Axis'
            else: xlabel = reg[0]
    
            reg = re.compile(r"@\s*yaxis\s*label\s*\"(.*?)\"", re.MULTILINE).findall(xvg)
            if reg == []: ylabel = 'Y Axis'
            else: ylabel = reg[0]
    
            reg = re.compile(r"@\s*s\w\s*legend\s*\"(.*?)\"", re.MULTILINE)
        
            cols = [ylabel+': '+m.groups()[0] for m in reg.finditer(xvg)]
            
            if cols:
                cols = [xlabel] + cols
            else: cols = [xlabel, ylabel]
            
            reg = re.compile(r"@.*", re.MULTILINE)
        
            ends = [m.end() for m in reg.finditer(xvg)]
            end = max(ends)
                
        else:
            xlabel = 'X'
            ylabel = 'Y'
        
        data = xvg[end+1:]
        ar = np.array(list(map(float, data.split())))
    
        reshaped = np.reshape(ar, (len(ar)//len(cols),len(cols)))
        df = pd.DataFrame(reshaped, columns=cols).set_index([xlabel])
        df.name = title
        return df
    
    @staticmethod
    def errors(file_eval=None):
    
        def _parse(out):
            output = _global.host.read(out)   
            return BaseHandle._notes(output), BaseHandle._gmxerrhdnl(output, dont_raise=True)
        # include ls() in local
        if file_eval is None:
            files = _global.host.cmd.ls(file_eval=lambda a: True if a.endswith('.log') or a.endswith('.out') else False)
        else:
            files = _global.host.cmd.ls(file_eval=file_eval)
        
        for file in files:
            print(f"\n<========{file}========>")
            
            notes, errors = _parse(file)
            
            if notes: print(f"\nThere were note(s)/warning(s):\n\n{notes}")
            if errors: print(f"\nThere were error(s):\n\n{errors}")
    
    @staticmethod
    def dump(file):
        ext = file.split('.')[-1]
        if ext == 's': cmd = 's'
        elif ext == 'trr' or ext == 'xtc' or ext == 'gro' or ext == 'pdb': cmd = 'f'
        elif ext == 'edr': cmd = 'e'
        elif ext == 'cpt': cmd = 'cp'
        elif ext == 'top': cmd = 'p'
        elif ext == 'mtx': cmd = 'mtx'
    
        kwargs = {cmd:file}
    
        d = BaseHandle()
        out = d.gmx(f'dump', **kwargs)
        return out
    
    @staticmethod
    def log(file):
        f = _global.host.read(file)
        x = re.compile(r"^\s*(Step\s*Time)\n(.+)\n(?:\n.+)+\s*Energies\s*\(kJ/mol\)((?:\n.+)+)\s+\n",\
                   re.MULTILINE)
        res = x.findall(f)
    
        def c(log_txt):
            colms = []
            vals = []
            n = 15
            colms += log_txt[0].split()
            vals += [float(i) for i in log_txt[1].split()]
            for j,l in enumerate(log_txt[2].splitlines()):
                if j%2 != 0:
                    cols = [l[i:i+n].strip() for i in range(0, len(l), n)]
                    colms.extend(cols)
                else:
                    cols = [float(l[i:i+n].strip()) for i in range(0, len(l), n)]
                    vals.extend(cols)
            return (colms, vals)
   
        cols, v1 = c(res[0])
    
        vals = []
    
        for i in res:
            vals.append(c(i)[1])
        
        df = pd.DataFrame(vals, columns=cols)
        df = df.set_index(['Step', 'Time'])
    
        return df
    
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


        
