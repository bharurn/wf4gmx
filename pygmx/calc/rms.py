from ..dashboard import PlotBoxDF
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSF

# rms calculation from mda

@PlotBoxDF
def d(top, trr, selections):
    u = mda.Universe(top)

    ref = mda.Universe(top)
    
    rmsd = []
    
    for l in trr:
        u.trajectory =  mda.Universe(l).trajectory
    
        rmsd.extend([(i.data['step'], i.data['step']*i.dt/10000, *[align.alignto(u.select_atoms(s), ref.select_atoms(s))[1]\
                      for s in selections]) for i in u.trajectory])
        u.trajectory.close()
    
    cols = [' '.join([w.title() if w.islower() else w for w in s.split()]) for s in selections]
    df = pd.DataFrame(rmsd, columns=('Step', 'Time', *cols))
    return df.set_index(['Step', 'Time'])

#@PlotBoxDF(x_axis=['Name', 'ID'])
def f(u, selection, align_with='protein'):
    align.AlignTraj(u, u, select=align_with, in_memory=True).run()
    
    sele = u.select_atoms(selection)
    rmsfer = RMSF(sele, verbose=True).run()
    
    x = [(a.name,a.id) for a in rmsfer.atomgroup.atoms]
    df = pd.DataFrame(x, columns=('Name', 'ID'))
    df['RMSF'] = rmsfer.rmsf
    return df.set_index(['Name', 'ID'])
