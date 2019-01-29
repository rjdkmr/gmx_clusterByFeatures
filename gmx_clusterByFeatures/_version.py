try:
    from .gmx_clusterByFeatures import gmx_version
    gmx_version = gmx_version()
except:
    pass

__version__ = "0.1.2"
