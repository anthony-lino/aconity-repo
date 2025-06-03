import pyslm
from libSLM import cli
solidpart = pyslm.Part("cubes")
solidpart.setGeometry("/CAD/cubes_Lydia_lab/flat_cubes.stl")
header = cli.Header()
header.filename = "flat_cubes.cli"
