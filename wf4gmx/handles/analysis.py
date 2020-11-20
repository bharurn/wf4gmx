from os import path
import xmlrpc.client as xmlrpclib
import pickle
import numpy as np
import pandas as pd
from mimicpy import Mpt, CoordsIO
from .base import BaseHandle
from ..io.trr import read_trr, get_trr_frames, TRRReader


class Analyze(BaseHandle):
    def __init__(self, mpt_file=None, trr_files=None, client=None, status_file=None, status=None):
        super().__init__(status_file, status)
        self.__getMpt(mpt_file)
        
        if trr_files is None:
            trr_files = self.gethistory('trr')[::-1]
        elif not isinstance(trr_files, list):
            trr_files = [trr_files]
            
        each_trr_frames = [get_trr_frames(trr) for trr in trr_files]
        self._trrs = list(zip(trr_files, each_trr_frames))
        
        self.nframes = sum(each_trr_frames)
        
        self.client = client
        
    def __getMpt(self, file):
        if file is None:
            file = self.getcurrent('mpt', exp=False)
            if not file:
                file = self.getcurrent('top')
        
        self.mpt = Mpt.from_file(file)
        
    @staticmethod
    def __get_frame(sele, trr_files, frame):
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
        
        return SelectedFrame(trr_data['step'], trr_data['time'], trr_data['box'], sele, x, v, f)
    
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
        
    def select(self, *selection, calc=None, frames=[], frame=0, asFutures=False, **kwargs):
        if isinstance(selection, str):
            sele = self.mpt.select(selection)
        else:
            sele = [self.mpt.select(s) for s in selection]
        
        trr_file = self._trrs
        
        if calc:
            do = lambda frame, **kwargs: calc(Analyze.__get_frame(sele, trr_file, frame), **kwargs)
        else:
            do = lambda frame: Analyze.__get_frame(sele, trr_file, frame)
        
        if frames == []:
            ret = do(frame)
        elif self.client is None:
            ret = [do(i) for i in frames]
        elif asFutures:
            ret = [self.client.submit(do, i, **kwargs) for i in frames]
        else:
            ret = [self.client.submit(do, i, **kwargs).result() for i in frames]
        
        return ret
    
    def select_coords(self, selection, filename):
        sele = self.mpt.select(selection)
        ids = sele.index - 1
        coords_handle = CoordsIO(filename, 'r')
        box = np.diag(coords_handle.box)
        
        x = coords_handle.coords[['x','y','z']].to_numpy()[ids]
        try:
            v = coords_handle.coords[['v_x','v_y','v_z']].to_numpy()[ids]
        except KeyError:
            v = np.array([])
            
        return SelectedFrame(0, 0, box, sele, x, v, np.array([]))

class SelectedFrame:
    
    def __init__(self, step, time, box, df, x, v, f):
        self.step = step
        self.time = time
        self.box = box
        self.natoms = len(df)
        if self.natoms == 1:
            one_only = True
        else:
            one_only = False
        
        if one_only:
            self.ids = df.index.to_numpy()[0]
        else:
            self.ids = df.index.to_numpy()
        
        self.__repr_list = []
        
        for column in df:
            if column == 'resid': dtype = np.int32
            elif column in ['charge', 'mass']: dtype = np.float32
            else: dtype = np.str
            self.__repr_list.append(column)
            
            if one_only:
                setattr(self, column, df[column].to_numpy(dtype=dtype)[0])
            else:
                setattr(self, column, df[column].to_numpy(dtype=dtype))
        
        if one_only and x != []:
            self.positions = x[0]
        else:
            self.positions = x
        
        if one_only and v != []:
            self.velocities = v[0]
        else:
            self.velocities = v
        
        if one_only and f != []:
            self.forces = f[0]
        else:
            self.forces = f
        self.i = 0
        self.__pbc_box = None
    
    def __getitem__(self, k):
        if isinstance(k, int): k = [k]
        dct = {}
        for i in self.__repr_list:
            dct[i] = getattr(self, i)[k]
        mpt_df = pd.DataFrame(dct, index=self.ids[k])
        
        try:
            pos = self.positions[k]
        except (TypeError,IndexError):
            pos = []
        
        try:
            v = self.velocities[k]
        except (TypeError,IndexError):
            v = []
        
        try:
            f = self.forces[k]
        except (TypeError,IndexError):
            f = []
            
        return SelectedFrame(self.step, self.time, self.box, mpt_df, pos, v, f)
    
    def sort(self, by=None, in_place=True):
        if by is None:
            by = self.ids
        sorted_ids = np.argsort(by)
        frame = self.__getitem__(sorted_ids)
        if in_place:
            self.__dict__.update(frame.__dict__)
        else:
            return frame
        
    def __repr__(self):
        s = f"\tStep: {self.step}\tTime: {self.time}\t#Atoms: {self.natoms}\n==========================================================\n\n"

        for k in self._SelectedFrame__repr_list:
            v = getattr(self, k)
            s += f"{k}: {v}\n"

        if self.positions.size != 0:
            s += f"\nPositions:\n==========\n{self.positions}\n"

        if self.velocities.size != 0:
            s += f"\nVelocities:\n==========\n{self.velocities}\n"

        if self.forces.size != 0:
            s += f"\nForces:\n==========\n{self.forces}\n"
        
        return s
    
    def __add__(self, other):
        dct = {}
        for i in self.__repr_list:
            dct[i] = np.append(getattr(self, i), getattr(other, i))
        mpt_df = pd.DataFrame(dct, index=np.append(self.ids, other.ids))

        pos = np.vstack((self.positions, other.positions))
        v = np.vstack((self.velocities, other.velocities))
        f = np.vstack((self.forces, other.forces))
                              
        return SelectedFrame(self.step, self.time, self.box, mpt_df, pos, v, f)
    
    def __iter__(self):
        return self

    def __next__(self):
        num = self.i
        self.i += 1
        if self.i > self.natoms:
            self.i = 0
            raise StopIteration
            
        return self.__getitem__(num)
    
    def get_around(self, sele, dist=3, inplace=True):
        cen = np.mean(self.positions, axis=0)
        ids = np.where( np.linalg.norm(sele.positions - cen, axis=1) <= dist)[0].tolist()
        if inplace:
            return (self + sele[ids])
        else:
            return sele[ids]

    def where(self, bool_id):
        return self.__getitem__(np.where(bool_id)[0].tolist())
    
    def write(self, file, extension='gro'):
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
        
        if file is None:
            return CoordsIO('dummy.'+extension, mode='w').write(mpt_df, coords_df, as_str=True)
        else:
            CoordsIO(file, mode='w').write(mpt_df, coords_df)
        
    def to_pymol(self, url='http://localhost:9123', name='mol1'):
        pymol = xmlrpclib.ServerProxy('http://localhost:9123')
        pymol.zoom()
        
        s = self.write(None).replace('\n','@')
        
        try:
            pymol.do(f'wf4gmx_load("{s}", "{name}")')
        except Fault:
            raise Exception("wf4gmx pymol file not sourced")
    
    def rmsd(self, other):
        pos1 = self.positions
        pos2 = other.positions
        if pos1.shape != pos2.shape:
            raise Exception
        N = pos1.shape[0]
        
        pos1 = pos1 - np.average(pos1, axis=0)
        pos2 = pos2 - np.average(pos2, axis=0)
        
        return np.sqrt(np.sum((pos1 - pos2) ** 2) / N)
    
    def __get_subset(self, s, frame):
            b = [frame.box[0,0], frame.box[1,1], frame.box[2,2]]

            box = [(b[0],0), (b[1],0), (b[2],0)]

            tup = (int(s[0]), int(s[1]), int(s[2]))

            for i, t in enumerate(tup):
                if t == 1:
                    box[i] = (b[i]/2, -b[i]/2) 

            return frame.where(
                (frame.positions[:, 0] <= box[0][0]) & 
                (frame.positions[:, 0] >= box[0][1]) & 
                (frame.positions[:, 1] <= box[1][0]) & 
                (frame.positions[:, 1] >= box[1][1]) & 
                (frame.positions[:, 2] <= box[2][0]) & 
                (frame.positions[:, 2] >= box[2][1])).sort(in_place=False)
        
    def __get_images(self):
        b1, b2, b3 = self.box[0,0], self.box[1,1], self.box[2,2]
        
        import copy
        sele1 = copy.deepcopy(self)
        sele2 = copy.deepcopy(self)
        sele3 = copy.deepcopy(self)
        sele4 = copy.deepcopy(self)
        sele5 = copy.deepcopy(self)
        sele6 = copy.deepcopy(self)
        sele7 = copy.deepcopy(self)
        
        sele1.positions[:,0] -= b1
        sele2.positions[:,1] -= b2
        sele3.positions[:,2] -= b3

        sele4.positions[:,0] -= b1
        sele4.positions[:,1] -= b2

        sele5.positions[:,0] -= b1
        sele5.positions[:,2] -= b3

        sele6.positions[:,0] -= b1
        sele6.positions[:,1] -= b2
        sele6.positions[:,2] -= b3

        sele7.positions[:,1] -= b2
        sele7.positions[:,2] -= b3
        
        return self+sele1+sele2+sele3+sele4+sele5+sele6+sele7
        
    def fix_pbc(self, ref, in_place=True):
        sel = self.__get_images()
        
        if ref.__pbc_box is None:
            rmsd = {}
            for i in range(7):
                s = "{0:03b}".format(i)
                e = self.__get_subset(s, sel)
                rmsd[s] = [e, e.rmsd(ref)]

            min_rmsd = min(rmsd, key=lambda x: rmsd[x][1])
            frame = rmsd[min_rmsd][0]
            frame.__pbc_box = min_rmsd
        else:
            frame = self.__get_subset(ref.__pbc_box, sel)
            frame.__pbc_box = ref.__pbc_box
        
        if in_place:
            self.__dict__.update(frame.__dict__)
        else:
            return frame