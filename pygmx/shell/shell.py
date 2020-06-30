# This is part of MiMiCPy

"""

This module contains classes to wrap the functionalities of local/remote shell into a common framework

"""

from . import _local as local, _remote as remote
import os
from stat import S_ISDIR, S_ISREG
from shutil import copyfile
import re
from ..utils.errors import SlurmBatchError
from ..utils.strs import f
from .._global import _Global as _gbl
from abc import ABC, abstractmethod

class Base(ABC):
    """

    Base shell class, contains methods common to both local and remote shell
    Cannot be instantiated, depends on functions from local/remote

    """
    
    @abstractmethod
    def __init__(self, directory,  *loaders):
        """ Init cwd and loaders """
        
        if directory.strip() == '':
            directory = '.'
        
        ret = self.cd(directory, mkdir=True)
        
        if ret == -1: # dir is home directory
            pass
        else:
            if ret == 1: # dir doesn't exist and was created by cd()
                _gbl.logger.write('debug', f"{directory} not found, creating new directory..")
            # if ret == 0: --> dir already exists
            _gbl.logger.write('debug', f"Setting current working directory to {self.pwd()}..")
        
        self.loaders = [] # list of loaders
        self.loader_str = '' # stringified version of loaders list
        self.loader_out = '' # output of running loaders
        
        if loaders: # add loaders to self.loaders
            self.addLoaders(*loaders)
        
    def addLoaders(self, *loader):
         self.loaders.extend(loader) # add loaders
         if self.loader_str.strip() == '': # check if loaders already there
             self.loader_str += ' ; '.join(loader) # join all loaders with ';'
         else:
             self.loader_str = ' ; '.join(loader) # join all loaders with ';'
        
         # replace defaultHook in run(), so it doesn't check for errors 
         print_hook = lambda cmd,out: _gbl.logger.write('debug2', f'Running loader {cmd}..\n{out}')
        
         # get output of loaders, so we can remove it from all subsequent run() outputs
         self.loader_out = self.run(self.loader_str, hook = print_hook, fresh=True)
    
    #convenience functions
    def rename(self, a, b): self.hndl().rename(a, b)
    def rm(self, a): self.hndl().remove(a)
    def pwd(self): return self.hndl().getcwd()
    def isLocal(self):
        if self.name == 'localhost': return True
        else: return False
    
    @abstractmethod
    def hndl(self): pass
    
    def checkFile(self, file, throw=False):
        """ Check if file exists, for debugging """
        ret = self.fileExists(file)
        
        if ret == True:
            _gbl.logger.write('debug', f"{file} found!")
        elif throw:
            raise FileNotFoundError(f"{file} not found!")
    
    def read(self, file, asbytes=False):
        """ Read and return contents of full file """
        mode = 'rb' if asbytes else 'r'
        with self.open(file, mode) as f:
            out = f.read()
            if not asbytes:
                # sftp file objects, when read, returns only bytes
                # so need to decode it if self is remote
                try:
                    return out.decode()
                except (UnicodeDecodeError, AttributeError): # if error, means out is string
                    return out
            else:
                # mode should be rb
                # so o/p will be in bytes, just return
                return out
    
    def write(self, content, file, asbytes=False):
        """ Write file """
        mode = 'wb' if asbytes else 'w'
        with self.open(file, mode) as f: f.write(content)
        
    def mkdir(self, directory):
        """ Make directory """
        if not self.fileExists(directory):
            self.hndl().mkdir(directory)
            return 0
        else: return 1
        
    def cd(self, directory, mkdir=False):
        """
        Change directory
        If mkdir is True, cd() will make a new directory if not found
        Otherwise it will raise exception
        """
        if directory.strip() == '.':
            return -1
        
        if not self.fileExists(directory):
             if mkdir:       
                self.hndl().mkdir(directory)
                return 1
             else:
                 raise FileNotFoundError(f'Directory {directory} not found')
        else:
            self.hndl().chdir(directory+'/')
            return 0
    
    def sbatch(self, job, dirc=''):
        """
        Write jobscript to file
        And run sbatch <jobscript.sh>
        """
        jbs = f'{dirc}/{job.name}.sh'
        
        self.write(str(job), jbs)
        
        if job.noCommands():
            raise SlurmBatchError(jbs, "No commands found!")
        
        jid = 0
        def _sbatch_err(txt):
            if 'error' in txt.lower():
                raise SlurmBatchError(jbs, txt)
            else:
                match = re.search(r'Submitted batch job (\w*)', txt)
                if match:
                    global jid
                    jid = match.groups()[0]
                else: raise SlurmBatchError(jbs, txt)
        
        self.run(f'sbatch {job.name}.sh', hook=_sbatch_err, dirc=dirc, fresh=True)
        
        return jid
    
    def scancel(self, jobid):
        """Cancel slurm job"""
        return self.run(f'scancel {jobid}', fresh=True)
        
class Local(Base):
    """

    Class to handle localhost specific shell function

    """
    
    def __init__(self, directory, shell_path, *loaders):
        """Init cwd and shell executaable path"""
        self.name = 'localhost'
        if shell_path:
            self.shell_path = shell_path
        else: # by default uses $SHELL to get shell path
            self.shell_path = os.environ['SHELL'] # check env vars and get shell path, works only for UNIX
        super().__init__(directory, *loaders)
    
    # return handle for file io
    def hndl(self): return os
    
    def fileExists(self, file):
        """Check if file exists"""
        return os.path.isfile(file) or os.path.isdir(file)
    
    def open(self, file, mode):
        """Wrapper around open()"""
        if 'r' in mode: self.checkFile(file, throw=True)
        return open(file, mode)
        
    def run(self, cmd, stdin=None, hook=None, fresh=False, dirc=''):
        """Local shell specific run command"""
        replace = ''
        if not fresh and self.loader_str:
            cmd = self.loader_str + ' ; ' + cmd # add loader string
            replace = self.loader_out # remove loader_out cmd from output
        if dirc != '': cmd = f'cd {dirc} ; ' + cmd # move to dirc if specified
        return local.run(cmd, self.shell_path, replace, stdin=stdin, hook=hook) # call run from _local.py
    
    def runbg(self, cmd, hook=None, fresh=False, dirc='', query_rate=1):
        """Local shell specific run background command"""
        replace = ''
        if not fresh:
            cmd = self.loader_str + ' ; ' + cmd
            replace = self.loader_out
        if not dirc != '': cmd = f'cd {dirc} ; ' + cmd
        return local.runbg(cmd, self.shell_path, replace, hook=hook) # call runbg from _local.py

    def ls(self, dirc=None, file_eval=lambda a: True, dir_eval=lambda a: True):
        """Local shell specific ls command to return files/folders as list"""
        return local.ls(dirc, file_eval, dir_eval)
    
    # copy command
    def cp(self, f1, f2): copyfile(f1, f2)
   
    def join(self, *args): return os.path.join(*args)
    
    def dirname(self, file): return os.path.dirname(file)
    
    # provided to match remote.__del__()
    def __del__(self): pass

class Remote(remote.Shell, Base):
    
    def __init__(self, work_dir, ssh_config, *loaders):
        
        if ':' in work_dir:
            r = work_dir.split(':')
            server = r[0]
            dir_ = r[1]
            if dir_ == '': dir_ = '.'
        else:
            server = work_dir
            dir_ = '.'
            
        _gbl.logger.write('debug', f("Setting remote machine {r[0]} as host.."))
        
        remote.Shell.__init__(self, server, ssh_config)
        
        Base.__init__(self, dir_, *loaders)
    
    def hndl(self): return self.sftp
        
    def fileExists(self, file):
        try:
            self.sftp.stat(file)
            return True
        except FileNotFoundError:
           return False
    
    def ls(self, dirc=None, file_eval=lambda a: True, dir_eval=lambda a: True):
        
        files = []
        
        if dirc:
            dirs = self.sftp.listdir_attr(dirc)
        else:
            dirs = self.sftp.listdir_attr()
        
        for entry in dirs:
            if S_ISDIR(entry.st_mode) and dir_eval(entry.filename):
                files.append(entry.filename)
            elif S_ISREG(entry.st_mode) and file_eval(entry.filename):
                files.append(entry.filename)
            
        return files
    
    def open(self, name, s, prefetch=True):
        if 'r' in s: self.checkFile(name, throw=True)
        file = self.sftp.open(name, s)
        if 'w' not in s and prefetch: file.prefetch()
        return file
    
    def cp(self, f1, f2): self.run(f"cp {f1} {f2}")
    
    def join(self, *args):
        # both unix and windows based ssh uses forward slash
        # see discussion of issue 306 on the fabric github page
        return "/".join([l for l in args if l != ''])

    def dirname(self, file):
        # see comment in join()
        return "/".join(i for i in file.split('/')[:-1])