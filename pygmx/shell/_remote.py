# This is part of MiMiCPy

"""

This module contains Shell class to take care of low level
paramiko code for opening a remote shell and running commands

"""

from os.path import expanduser
import time
from ..utils.errors import MiMiCPyError, defaultHook

class Shell:
    def __init__(self, name, config=None):
        
        self.name = name
        
        if config == None:
            config = expanduser("~")+'/.ssh/config'
        
        self.config = config
        
        conf = self._lookup(name)
        
        if 'proxyjump' in conf:
            proxy_name, port = conf['proxyjump'].split(':')
            vmtransport = self._getssh(proxy_name).get_transport()
            dest_addr = (self._lookup(name)['hostname'], int(port))
            local_addr = (self._lookup(proxy_name)['hostname'], int(port))
            vmchannel = vmtransport.open_channel("direct-tcpip", dest_addr, local_addr)
            
            self.ssh = self._getssh(name, sock=vmchannel)
        
        else:
            self.ssh =  self._getssh(name)
            
        self.sftp = self.ssh.open_sftp()
        self.decoder = 'utf-8'
        self.is_open = True
        
        self.ssh_bg = None
        
    def _lookup(self, name):
        try:
            import paramiko
        except ImportError:
            raise MiMiCPyError("Paramiko python package not installed! Install it to remotely run commands")
        
        config = paramiko.SSHConfig()
        try:
            f = open(self.config)
        except FileNotFoundError:
            raise FileNotFoundError(f"No SSH config file found in {self.config}")
            
        config.parse(f)
        conf = config.lookup(name)
        
        return conf
    
    def _getssh(self, name, sock=None):
        try:
            import paramiko
        except ImportError:
            raise MiMiCPyError("Paramiko python package not installed! Install it to remote host")
            
        conf = self._lookup(name)
    
        if 'hostname' not in conf:
            raise MiMiCPyError(f"Hostname not found for {name} in SSH config file")
        if 'user' not in conf:
            raise MiMiCPyError(f"Username not found for {name} in SSH config file")
        
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        if sock is None:
            ssh.connect(hostname=conf['hostname'], username=conf['user'], look_for_keys=False)
        else:
            ssh.connect(hostname=conf['hostname'], username=conf['user'], look_for_keys=False, sock=sock)
        
        return ssh
    
    def __del__(self):
        if not hasattr(self, 'is_open'): return
        if self.is_open:
            self.is_open = False
            self.sftp.close()
            self.ssh.close()
            if self.ssh_bg: self.ssh_bg.__del__()
          
    def runbg(self, cmd, hook=None, dirc='', query_rate=3):
        if self.ssh_bg == None:
            self.ssh_bg = _ShellBGRun(self.ssh, self.decoder)
            
        self.ssh_bg.run(f'cd {self.pwd()}/{dirc}')
        self.ssh_bg.run(self.loader_str)
        if not hook:
            hook = defaultHook
            
        out = self.ssh_bg.run(cmd.split(';')[-1], hook, query_rate)
        
        return out
        # do not destory ssh_bg, as this stops the process
        # kill it in the deconstructor
    
    def run(self, cmd, stdin=None, hook=None, fresh=False, dirc=''):
        
        if not fresh and self.loader_str:
            cmd = self.loader_str + ' ; ' + cmd
        
        cmd = f'cd {self.pwd()}/{dirc} ; ' + cmd
        
        tran = self.ssh.get_transport()
        chan = tran.open_session()
        chan.get_pty()
        stdout = chan.makefile('rb')
        chan.exec_command(cmd)
        
        if stdin:
            sin = chan.makefile('wb')
            sin.channel.send(stdin+'\n')
            sin.channel.shutdown_write()
        
        out = stdout.read().decode(self.decoder)
        
        if not fresh:
            out = out.replace(self.loader_out, '')
    
        if not hook:
            hook = defaultHook
        
        hook(cmd.split(';')[-1], out)
        
        return out
        
    
    def __enter__(self): return self
        
    def __exit__(self, exc_type, exc_val, exc_tb): self.__del__()
    
class _ShellBGRun:
    
    def __init__(self, ssh, decoder):
        self.chan = ssh.invoke_shell()
        self.stdout = ''
        self.decoder = decoder
    
    def queryStdout(self, buff = 1024):
        
       while self.chan.recv_ready():
           self.stdout += self.chan.recv(buff).decode(self.decoder)
           
    def run(self, cmd, hook=None, query_rate=0.2):
        
        self.queryStdout()
        
        startout = self.stdout
        
        self.chan.send(cmd+'\n')
            
        prev = startout
        
        while True:
            time.sleep(query_rate)
                
            self.queryStdout()
                
            if prev != self.stdout:
                text = self.stdout.replace(prev, '')
                    
                prev = self.stdout
                    
                if hook:
                    if hook(cmd, text): break
                
            else: break
        
        lines = self.stdout.replace(startout, '').splitlines()[1:-1]
        
        return '\n'.join(lines)
    
    def __del__(self):
        self.chan.close()