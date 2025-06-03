# Example comment
import os
import pyslm
from libSLM import cli
print(os.listdir("/CAD"))
solidpart = pyslm.Part("cubes")
solidpart.setGeometry("/CAD/cubes_Lydia_lab/flat_cubes.stl")
header = cli.Header()
header.filename = "flat_cubes.cli"
