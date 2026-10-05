# -*- python -*-
from lsst.sconsUtils import scripts
from lsst.sconsUtils.state import env


# Python-only package
# Force shebang and policy to come first so the file first appears in the bin
# directory before it is used. This is required to run on macos.
PKG_ROOT = env.ProductDir("scene_modeler")

scripts.BasicSConstruct("scene_modeler", disableCc=True)