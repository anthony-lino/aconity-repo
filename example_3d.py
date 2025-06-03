"""
A simple example showing how to use PySLM for generating slices across a 3D model
"""
import pyslm
from pyslm import hatching as hatching
import numpy as np
import numpy as np

# Imports the part and sets the geometry to  an STL file (frameGuide.stl)
solidPart = pyslm.Part('inversePyramid')
solidPart.setGeometry('/CAD/examples/inversePyramid.stl')

solidPart.origin = [5.0, 10.0, 0.0]
solidPart.scaleFactor = 1.0
solidPart.rotation = [0, 0.0, 45]
solidPart.dropToPlatform()
print(solidPart.boundingBox)


# Create a BasicIslandHatcher object for performing any hatching operations
myHatcher = hatching.BasicIslandHatcher()
myHatcher.islandWidth = 3.0
myHatcher.stripeWidth = 5.0

# Set the base hatching parameters which are generated within Hatcher
myHatcher.hatchAngle = 10 # [°] The angle used for the islands
myHatcher.volumeOffsetHatch = 0.08 # [mm] Offset between internal and external boundary
myHatcher.spotCompensation = 0.06 # [mm] Additional offset to account for laser spot size
myHatcher.numInnerContours = 2
myHatcher.numOuterContours = 1
myHatcher.hatchSortMethod = hatching.AlternateSort()

# Set the layer thickness
layerThickness = 0.03 # [mm]

# Perform the slicing. Return coords paths should be set so they are formatted internally.
myHatcher.layerAngleIncrement = 66.7

#Perform the hatching operations
print('Hatching Started')

layers = []

for z in np.arange(0, solidPart.boundingBox[5], layerThickness):

    # Typically the hatch angle is globally rotated per layer by usually 66.7 degrees per layer
    #myHatcher.hatchAngle += 66.7
    # Slice the boundary
    geomSlice = solidPart.getVectorSlice(z,simplificationFactor=0.1)
    # Hatch the boundary using myHatcher
    layer = myHatcher.generateHatching(paths = geomSlice,hatchSpacing = .05)
    print(type(layer))
    
    # The layer height is set in integer increment of microns to ensure no rounding error during manufacturing

    layers.append(layer)

print('Completed Hatching')

# Plot the layer geometries using matplotlib
# Note: the use of python slices to get the arrays
pyslm.visualise.plotLayers(layers[0:-1:10])


