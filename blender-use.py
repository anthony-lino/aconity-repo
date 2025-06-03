import bpy
import sys
from mathutils import Vector

# Get STL path from CLI args
stl_path = sys.argv[-1]

# Reset Blender scene
bpy.ops.wm.read_factory_settings(use_empty=True)

# Import STL
bpy.ops.import_mesh.stl(filepath=stl_path)

obj = bpy.context.selected_objects[0]

# Convert bound_box points to Vector, then transform by matrix_world
bbox = [obj.matrix_world @ Vector(v) for v in obj.bound_box]

min_x = min(v.x for v in bbox)
max_x = max(v.x for v in bbox)
min_y = min(v.y for v in bbox)
max_y = max(v.y for v in bbox)
min_z = min(v.z for v in bbox)
max_z = max(v.z for v in bbox)

print(f"X: {min_x:.4f} to {max_x:.4f}")
print(f"Y: {min_y:.4f} to {max_y:.4f}")
print(f"Z: {min_z:.4f} to {max_z:.4f}")
