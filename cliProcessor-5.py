import numpy as np
from typing import List
import struct

import sys
import os
from copy import deepcopy

import pyslm
import pyslm.visualise
import pyslm.analysis
import pyslm.geometry
from pyslm import hatching as hatching
import numpy as np

"""
The structures necessary for creating a machine build file should be imported from the geometry submodule. 
Fallback python classes that are equivalent to those in libSLM are provided to ensure prototyping can take place.
"""
from pyslm import geometry as slm

def writeBuildStyle(fp, models, mid, bid):
    """
    If the cli-plus format is used, the associated laser parameters for each :class:`LayrerGeometry` are
    written inline into the .cli file for each geometry section
    """
    bstyle = pyslm.geometry.getBuildStyleById(models, mid, bid)

    fp.write('$$POWER/{:.1f}\n'.format(bstyle.laserPower))
    fp.write('$$SPEED/{:.1f}\n'.format(bstyle.laserSpeed))
    fp.write('$$FOCUS/{:.1f}\n'.format(bstyle.laserFocus))

    if bstyle.jumpSpeed > 0:
        fp.write('$$PARAM/jump_speed, double, {:.1f}\n'.format(bstyle.jumpSpeed))

    if bstyle.jumpDelay > 0:
        fp.write('$$PARAM/jump_delay, double, {:.1f}\n'.format(bstyle.jumpDelay))

    if hasattr(bstyle, 'markDelay'):
        fp.write('$$PARAM/mark_delay, double, {:.1f}\n'.format(bstyle.markDelay))

    if hasattr(bstyle, 'laserOnDelay'):
        fp.write('$$PARAM/laser_on_delay, double, {:.1f}\n'.format(bstyle.laserOnDelay))

    if hasattr(bstyle, 'laserOffDelay'):
        fp.write('$$PARAM/laser_off_delay, double, {:.1f}\n'.format(bstyle.laserOffDelay))

    if hasattr(bstyle, 'polygonDelay'):
        fp.write('$$PARAM/polygon_delay, double, {:.1f}\n'.format(bstyle.polygonDelay))

    if hasattr(bstyle, 'laserPowerRate'):
        fp.write('$$PARAM/laser_power_rate, double, {:.3f}\n'.format(bstyle.laserPowerRate))

    if hasattr(bstyle, 'spotDiameter'):
        fp.write('$$PARAM/spot_diameter, double, {:.3f}'.format(bstyle.spotDiameter))

    if hasattr(bstyle, 'disableDefocusCalibration'):
        fp.write('$$PARAM/disable_defocus_calculation, bool, {:d}\n'.format(bstyle.disableDefocusCalibration))

    if hasattr(bstyle, 'defocusSetPoint'):
        fp.write('$$PARAM/defocus_set_point, double, {:.3f}\n'.format(bstyle.defocusSetPoint))

    if hasattr(bstyle, 'markSpeed'):
        fp.write('$$PARAM/mark_speed, double, {:.3f}\n'.format(bstyle.markSpeed))

def checkParamSize(params, geom):
    if isinstance(geom, slm.HatchGeometry):
        numVecs = geom.coords.shape[0] / 2
    elif isinstance(geom, slm.ContourGeometry):
        numVecs = geom.coords.shape[0] - 1
    elif isinstance(geom, slm.PointsGeometry):
        numVecs = geom.coords.shape[0]
    else:
        raise Exception('Invalid LayerGeometry structures has been provided')

    if len(params) != numVecs:
        raise Exception('Parameters must match the size of the geometry')

    return True

def writeCLI(filename: str,
                      header: slm.Header,
                      models: List[slm.Model],
                      layers: List[slm.Layer],
                      isCliPlus: bool = False):
    """
    Writes an ASCII Cli Plus File
    :param filename:  Filename to write
    :param header:  Header information for the .cli file
    :param models: A list of  :class:`Model` describing Laser Parameters used if isCliPlus is `True`
    :param layers:  A list of :class:`Layer` objects containing the geometry to be written to the .cli file
    :param isCliPlus: If `True` uses the cli+ functionality for Aconitiy L-PBF systems
    """
    try:
        os.remove(filename)
    except:
        pass

    fp = open(filename, 'w')

    """ 
    Write heder information
    """
    fp.write('$$HEADERSTART\n')
    fp.write('$$ASCII\n')
    fp.write('$$UNITS/{:.5f}\n'.format(header.scaleFactor))
    fp.write('$$LAYERS/{:d}\n'.format(len(layers)))
    fp.write('$$HEADEREND\n')

    """
    Write the geometry section
    """
    fp.write('$$GEOMETRYSTART\n')

    for layer in layers:
        fp.write('$$LAYER/{:.4f}\n'.format((layer.z)*header.scaleFactor))

        for geom in layer.geometry:

            if isCliPlus:
                # Check if there are cusotm parameters embedded within the geometry
                writeBuildStyle(fp, models, geom.mid, geom.bid)

            #validParams = ['wait']
            #attr in dir(geom)

            if isinstance(geom, slm.HatchGeometry):
                # Write hatch geometry
                numHatches = int(geom.coords.shape[0] / 2)
                coords = geom.coords.ravel() / header.scaleFactor
                #coordsStr = np.array2string(coords.ravel(), precision=3 ,prefix='', suffix='', threshold=sys.maxsize,
                #                            formatter={'float_kind':lambda x: "%.3f" % x}, separator=',')[1:-1]
                coordsStr = ','.join('{:.3f}'.format(x) for x in coords.ravel())
                fp.write('\n$$HATCHES/')
                fp.write('{:d},{:d},{:s}\n'.format(geom.mid, numHatches, coordsStr))

            elif isinstance(geom, slm.ContourGeometry):
                numVecs = int(geom.coords.shape[0]) - 1
                coords = geom.coords.ravel() / header.scaleFactor
                #coordsStr = np.array2string(coords, precision=3 ,prefix='', suffix='',threshold=sys.maxsize,
                #                            formatter={'float_kind':lambda x: "%.3f" % x}, separator=',')[1:-1]
                coordsStr = ','.join('{:.3f}'.format(x) for x in coords.ravel())
                fp.write('$$POLYLINE/')
                fp.write('{:d},0,{:d},{:s}\n'.format(geom.mid, numVecs, coordsStr))

            elif isinstance(geom, slm.PointsGeometry):
                # Note that Point geometries are reprsented by zero length hatch vectors
                numPoints  = int(geom.coords.shape[0])

                coords = geom.coords / header.scaleFactor
                coordsCpy = np.tile(coords, (1,2))
                #coordsStr = np.array2string(coordsCpy.ravel(), precision=3 ,prefix='', suffix='',threshold=sys.maxsize,
                #                            formatter={'float_kind':lambda x: "%.3f" % x}, separator=',')[1:-1]
                coordsStr = ','.join('{:.3f}'.format(x) for x in coordsCpy.ravel())
                fp.write('$$HATCHES/')
                fp.write('{:d},{:d},{:s}\n'.format(geom.mid, numPoints, coordsStr))

    fp.write('$$GEOMETRYEND\n')

def writeBinaryCLI(filename: str,
                      header: slm.Header,
                      models: List[slm.Model],
                      layers: List[slm.Layer]):
    """
    Writes a binary .cli file

    :param filename: Filename
    :param header:  Header information for the .cli file
    :param models:  A list of  :class:`Model` describing Laser Parameters used if isCliPlus is `True`
    :param layers:  A list of :class:`Layer` objects containing the geometry to be written to the .cli file
    """
    def writeStr(fp, value):
        fp.write(value.encode('ascii'))  # write as an ascii encoded string

    def writeFloat(fp, value):
        if isinstance(value, float):
            fp.write(struct.pack('f', value))  # write a float
        elif isinstance(value, np.ndarray):
            from numpy import array
            a = array(value, 'float32')
            a.tofile(fp)

    def writeLong(fp, value):
        if isinstance(value, int):
            fp.write(struct.pack('i', value))  # write as a long integer

    def writeUInt(fp, value):
        if isinstance(value, int):
            fp.write(struct.pack('H', value))  # write an unsigned short integer

    try:
        os.remove(filename)
    except:
        pass

    fp = open(filename, 'wb')

    """ 
    Write heder information
    """
    writeStr(fp, '$$HEADERSTART\n')
    writeStr(fp, '$$BINARY\n')
    writeStr(fp, '$$VERSION/200\n')
    writeStr(fp, '$$UNITS/{:.5f}\n'.format(header.scaleFactor))
    writeStr(fp,'$$LAYERS/{:d}\n'.format(len(layers)))
    writeStr(fp,'$$HEADEREND')

    """
    Write the geometry section
    """
    for layer in layers:
        # Write layer
        writeUInt(fp, 127)
        z = float(layer.z)
        print('Writing z', z, z/1000.0)
        writeFloat(fp,z/1000.0)

        for geom in layer.geometry:

            if isinstance(geom, slm.HatchGeometry):
                # Write hatch geometry
                numHatches = int(geom.coords.shape[0] / 2)
                coords = geom.coords.ravel() / header.scaleFactor
                writeUInt(fp, 132)
                writeLong(fp, geom.mid)
                writeLong(fp, numHatches)
                writeFloat(fp, coords)

            elif isinstance(geom, slm.ContourGeometry):
                numVecs = int(geom.coords.shape[0])
                coords = geom.coords.ravel() / header.scaleFactor
                writeUInt(fp, 130)
                writeLong(fp, 0)  # note should be geom.mid
                writeLong(fp, 0)
                writeLong(fp, numVecs)
                writeFloat(fp, coords)

            elif isinstance(geom, slm.PointsGeometry):
                # Note that Point geometries are represented by zero length hatch vectors

                coords = geom.coords / header.scaleFactor
                coordsCpy = np.tile(coords, (1,2))
                numPoints = int(coordsCpy.shape[0] / 2)

                writeUInt(fp, 132)
                writeLong(fp, geom.mid)
                writeLong(fp, numPoints)
                writeFloat(fp, coords)

def processLayer(fp: object, units: object = 1.0) -> slm.Layer:
    """
    Internal function that parses the Layer and extracts both ContourGeometry
    and HatchGeometry  that is then stored  within a  :class:`Layer`

    :param fp:
    :param units: Scale factor to apply to the geometry within the processed layer
    :return: A :class:`Layer`
    """
    curPos = fp.tell()
    line = fp.readline()

    layer = slm.Layer()

    while line:

        if "$$HATCHES" in line.upper():
            geom = slm.HatchGeometry()
            tag, value= line.split('/')
            vec = value.split(',')
            geom.bid = int(vec[0])
            geom.mid = 1
            numHatches = int(vec[1])

            if numHatches * 4 != len(vec) -2:
                raise Exception('Malformed hatch line')

            coords = np.array([float(val) for val in vec[2:]]).reshape(-1, 2)
            coords *= units
            geom.coords = coords
            layer.geometry.append(geom)

        elif "$$POLYLINE" in line.upper():
            geom = slm.ContourGeometry()
            tag, value= line.split('/')
            vec = value.split(',')
            geom.bid = int(vec[0])
            geom.mid = 1
            numContours = int(vec[2])

            if numContours * 2 != len(vec) -3:
                raise Exception('Malformed contour (Polyline) line')

            coords = np.array([float(val) for val in vec[3:]]).reshape(-1,2)
            coords *= units
            geom.coords = coords
            layer.geometry.append(geom)
        else:
            break

        curPos = fp.tell()
        line = fp.readline()

    fp.seek(curPos)

    return layer

def readCLI(filepath: str)  -> List[pyslm.geometry.Layer]:
    """
    Reads an ASCII .cli file and returns a list of layers with the contour and hatch geometry.

    :param filepath: The path to the .cli file.
    :return: A list of `slm.Layer` objects containing the parsed geometry.
    """
    units = 1.0
    numLayers = 0
    layers = []

    with open(filepath) as fp:

        """ Process the Header """
        if "$$HEADERSTART" not in fp.readline():
            raise Exception('invalid .cli provided')

        if "$$ASCII" not in fp.readline():
            raise Exception('.cli is not ASCII format')
        else:
            print('.cli ASCII Header Preamble located')

        cnt = 2
        line = fp.readline()
        layerId = 1

        """ Process the header first"""
        while line:

            if "$$UNITS" in line.upper():
                tag, value = line.split('/')
                units = float(value)
            elif "$$UNITS" in line.upper():
                tag, value = line.split('/')
                units = float(value)
                print("Scale factor: ({:.3f})".format(units))
            elif "$$LAYERS" in line.upper():
                tag, value = line.split('/')
                numLayers = int(value)
                print("Num layers: ({:d})".format(numLayers))
            elif "$$HEADEREND" in line.upper():
                print("Processed header ({:s})".format(filepath))
                break

            line = fp.readline()
            cnt += 1

        print('.cli units: {:.3f}'.format(units))
        """ Process the geoemtry """
        while line:
            if "$$GEOMETRYSTART" in line.upper():
                print('Reading geometry sections')
            elif "$$GEOMETRYEND" in line.upper():
                print('Processed .cli file')
                break
            elif "$$LAYER" in line.upper():
                tag, value = line.split('/')
                layerZ = float(value)

                layer = processLayer(fp, units)
                layer.layerId = layerId
                layer.z = int(layerZ * 1000) # Multiply the z-unit value for the layer position
                layers.append(layer)
                layerId += 1

            line = fp.readline()

    return layers

def readBinaryCli(filepath: str) -> List[pyslm.geometry.Layer]:
    """
    Reads a binary .cli file and returns a list of layers with the contour and hatch geometry

    :param filepath:
    :return:
    """
    def readBinaryLine(fp):
        data = bytearray()
        while True:
            byte = fp.read(1)
            if not byte or byte == b'\n':
                break
            data.extend(byte)

        return data.decode('ascii')

    def readAsciiStr(fp, length):
        return fp.read(length).decode('ascii')

    def readFloat(fp):
        return struct.unpack('f', fp.read(4))[0]

    def readLong(fp):
        return struct.unpack('i', fp.read(4))[0]

    def readUInt(fp):
        return struct.unpack('H', fp.read(2))[0]

    layers = []

    with open(filepath, 'rb') as fp:

        line = readBinaryLine(fp)

        if line != '$$HEADERSTART':
            raise Exception('Invalid .cli file')

        line = readBinaryLine(fp)

        if line != '$$BINARY':
            raise Exception('.cli is not binary format')

        line = readBinaryLine(fp)

        units = 1.0
        numLayers = 0

        while line != '$$HEADEREND':

            if line.startswith('$$VERSION'):
                if line != '$$VERSION/200':
                    raise Exception('Unsupported .cli version')

            if line.startswith('$$UNITS'):
                units = float(line.split('/')[1])
            elif line.startswith('$$LAYERS'):
                numLayers = int(line.split('/')[1])

            """
            Store the current binary file position
            Read-ahead and determine whether  $$HEADEREND is reached
            Otherwise iterate and read the next line
            """

            curPos = fp.tell()
            tempLine = fp.read(11).decode('ascii')
            if tempLine == '$$HEADEREND':
                break
            else:
                fp.seek(curPos)
                line = readBinaryLine(fp).strip()

        """ Iterate and read across the total number of layers available """
        for _ in range(numLayers):
            code = readUInt(fp)

            if code != 127:
                raise Exception('Invalid Layer Start Code {}'.format(code))

            layerZ = readFloat(fp) * 1000

            layer = slm.Layer()
            layer.z = int(layerZ)

            while True:
                curPos = fp.tell()

                try:
                    geomType = readUInt(fp)
                except:
                    break

                if geomType == 132:  # Hatch geometry
                    mid = readLong(fp)
                    numHatches = readLong(fp)
                    coords = np.fromfile(fp, dtype=np.float32, count=numHatches * 4).reshape(-1, 2) * units
                    geom = slm.HatchGeometry()
                    geom.mid = mid
                    geom.coords = coords
                    layer.geometry.append(geom)
                elif geomType == 130:  # Contour geometry
                    mid = readLong(fp)
                    bid = readLong(fp)
                    numVecs = readLong(fp)
                    coords = np.fromfile(fp, dtype=np.float32, count=numVecs * 2).reshape(-1, 2) * units
                    geom = slm.ContourGeometry()
                    geom.mid = mid
                    geom.bid = bid
                    geom.coords = coords
                    layer.geometry.append(geom)
                elif geomType == 127 or geomType == 128:
                    fp.seek(curPos)
                    # next layer so go back to the start of the next layer
                    break;
                else:
                    raise Exception('Invalid Code {} was parsed'.format(geomType))
                    break

            layers.append(layer)

    return layers

""" Main execution code """

if __name__ == "__main__":
    # Imports the part and sets the geometry to  an STL file (frameGuide.stl)
    solidPart = pyslm.Part('cube')
    solidPart.setGeometry('/CAD/examples/cube.stl')

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
    
    # Robust layer data conversion with debugging
    layers = []

    def convert_to_2d_coords(data, description=""):
        """Helper function to convert various data formats to 2D coordinate arrays"""
        try:
            if isinstance(data, np.ndarray):
                if data.ndim == 2:
                    return data.astype(np.float32)
                elif data.ndim == 1:
                    if len(data) % 2 == 0:
                        return data.reshape(-1, 2).astype(np.float32)
                    else:
                        print(f"Warning: {description} length {len(data)} is not even, cannot reshape to pairs")
                        return None
                else:
                    print(f"Warning: {description} has {data.ndim} dimensions, expected 1 or 2")
                    return None
            elif isinstance(data, list):
                arr = np.array(data, dtype=np.float32)
                return convert_to_2d_coords(arr, description)
            else:
                print(f"Warning: {description} has unsupported type {type(data)}")
                return None
        except Exception as e:
            print(f"Error converting {description}: {e}")
            return None

    for i, z in enumerate(np.arange(0, solidPart.boundingBox[5], layerThickness)):
        myHatcher.hatchAngle += 66.7
        # Slice the boundary
        geomSlice = solidPart.getVectorSlice(z, simplificationFactor=0.1)
        # Hatch the boundary using myHatcher
        layerData = myHatcher.generateHatching(paths=geomSlice, hatchSpacing=.05)
        
        # Debug: Print the structure of layerData (only for first layer)
        if i == 0:
            print(f"Layer {i}, Z={z}")
            print(f"layerData type: {type(layerData)}")
            if hasattr(layerData, '__len__'):
                print(f"layerData length: {len(layerData)}")
            
            # If it's a list or array, check the first few elements
            if isinstance(layerData, (list, np.ndarray)) and len(layerData) > 0:
                print(f"First element type: {type(layerData[0])}")
                if hasattr(layerData[0], 'shape'):
                    print(f"First element shape: {layerData[0].shape}")
                elif isinstance(layerData[0], (list, tuple)):
                    print(f"First element length: {len(layerData[0])}")
        
        # Create an slm.Layer object
        layer = slm.Layer(id=i, z=int(z * 1000))  # Convert z to micrometers
        
        # Process the layer data
        try:
            # Case 1: layerData is a single 2D array (all hatches in one array)
            if isinstance(layerData, np.ndarray) and layerData.ndim == 2:
                coords_2d = convert_to_2d_coords(layerData, f"layerData for layer {i}")
                if coords_2d is not None:
                    hatchGeom = slm.HatchGeometry(mid=1, bid=1)
                    hatchGeom.coords = coords_2d
                    layer.appendGeometry(hatchGeom)
            
            # Case 2: layerData is a single 1D array (flattened coordinates)
            elif isinstance(layerData, np.ndarray) and layerData.ndim == 1:
                coords_2d = convert_to_2d_coords(layerData, f"layerData for layer {i}")
                if coords_2d is not None:
                    hatchGeom = slm.HatchGeometry(mid=1, bid=1)
                    hatchGeom.coords = coords_2d
                    layer.appendGeometry(hatchGeom)
            
            # Case 3: layerData is a list of paths/hatches
            elif isinstance(layerData, (list, tuple)):
                for j, hatchPath in enumerate(layerData):
                    coords_2d = convert_to_2d_coords(hatchPath, f"hatchPath {j} in layer {i}")
                    if coords_2d is not None:
                        hatchGeom = slm.HatchGeometry(mid=1, bid=1)
                        hatchGeom.coords = coords_2d
                        layer.appendGeometry(hatchGeom)
            
            # Case 4: layerData has .hatches or similar attribute (pyslm specific)
            elif hasattr(layerData, 'hatches'):
                for j, hatch in enumerate(layerData.hatches):
                    coords_2d = convert_to_2d_coords(hatch, f"hatch {j} in layer {i}")
                    if coords_2d is not None:
                        hatchGeom = slm.HatchGeometry(mid=1, bid=1)
                        hatchGeom.coords = coords_2d
                        layer.appendGeometry(hatchGeom)
            
            else:
                print(f"Unknown layerData structure for layer {i}: {type(layerData)}")
                if hasattr(layerData, '__dict__'):
                    print(f"Available attributes: {list(layerData.__dict__.keys())}")
        
        except Exception as e:
            print(f"Error processing layer {i}: {e}")
            import traceback
            traceback.print_exc()
        
        layers.append(layer)
                    
    print('Completed Hatching')
    
    # Create a model - FIXED: Properly set the topLayerId
    from pyslm.geometry import Model
    model = Model(mid=1, topSliceNum=len(layers))  # Set topSliceNum to number of layers
    
    # CRITICAL FIX: Set the topLayerId to the last layer's ID
    model.mid = 1
    model.topLayerId = len(layers) - 1  # Set to the ID of the last layer (0-indexed)
    
    # Also set the name to avoid the empty name issue
    model.name = "InversePyramid"
    
    from pyslm.geometry import BuildStyle

    bstyle = pyslm.geometry.BuildStyle()

    # Attach the build styles to the model
    bstyle.setStyle(
        bid=1,
        focus=0.5,
        power=200,
        pointExposureTime=20,
        pointExposureDistance=50,
        speed=800,
        laserId=1,
        laserMode=pyslm.geometry.LaserMode.Default
    )
    model.buildStyles.append(bstyle) 
    
    print(f"Model validation info: mid={model.mid}, topLayerId={model.topLayerId}, name='{model.name}'")
    print(f"Number of layers: {len(layers)}")
    
    """ Validate the input model """
    pyslm.geometry.ModelValidator.validateBuild([model], layers)

    """
    Generate the header for use in the .cli file
    """
    header = pyslm.geometry.Header()
    header.filename = "Cli File "
    header.version = (1, 0)
    header.zUnit = 1
    header.scaleFactor = .005

    """
    Note:  Both a model and header is required. The scale factor can be 1.0 
    as floating points are defined explicitly as textual literals. 
    
    """

    # Write ASCII Version
    writeCLI('filename.cli', header, [model], layers)
    print("Successfully generated CLI file: filename.cli")

    # Write Binary Version (uncomment if needed)
    # writeBinaryCLI('filename_binary.cli', header, [model], layers)