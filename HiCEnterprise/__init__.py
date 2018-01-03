from pkg_resources import get_distribution, DistributionNotFound

try:
    __version__ = get_distribution("HiCEnterprise").version
except DistributionNotFound:
    VERSION = "HiCEnterprise-(local)"
else:
    VERSION = "HiCEnterprise-" + __version__
