from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
class Storage:
    def __init__(self):
        self._HTTP= HTTPRemoteProvider()
        self._FTP = FTPRemoteProvider()
    def http(self, *args, **kwargs):
        return self._HTTP.remote(*args, **kwargs)
    def ftp(self, *args, **kwargs):
        return self._FTP.remote(*args, **kwargs)

storage = Storage()