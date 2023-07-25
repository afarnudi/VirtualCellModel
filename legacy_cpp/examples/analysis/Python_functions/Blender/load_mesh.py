import bpy
import bmesh
import mathutils

import numpy as np


def build_mesh_from_xyz(xyz, name="AllSpheres", radius=1.0, colour=(0.8,0.7,0.2,1.0) ):
    MatAu = bpy.data.materials.new('Mat.Au')
    MatAu.diffuse_color = colour
    MatC  = bpy.data.materials.new('Mat.C')
    MatC.diffuse_color = (0.1,0.1,0.1,1.0)
    MatH  = bpy.data.materials.new('Mat.H')
    MatH.diffuse_color = (0.8,0.7,0.2,1.0)

    # create an empty mesh object and add it to the scene
    sphereMesh = bpy.data.meshes.new(name)
    sphereObj  = bpy.data.objects.new(name, sphereMesh)
    bpy.context.collection.objects.link(sphereObj)
    bpy.context.view_layer.objects.active = sphereObj

    # create slots for each material then add each material to the object
    #while len(sphereObj.material_slots) < 3:
    while len(sphereObj.material_slots) < 1:
        bpy.ops.object.material_slot_add()
    sphereObj.material_slots[0].material = MatAu
    #sphereObj.material_slots[1].material = MatC
    #sphereObj.material_slots[2].material = MatH
    
    # test data to be swapped for data in file
    # type, location, scale
    #data = [('C', (1.0, 1.0, 1.0), 0.2),
    #        ('H', (2.0, 2.0, 2.0), 0.4),
    #        ('Au',(3.0, 3.0, 3.0), 0.8),
    #        ('C', (4.0, 4.0, 4.0), 1.2),
    #        ('H', (5.0, 5.0, 5.0), 0.2),
    #        ('Au',(6.0, 6.0, 6.0), 0.8),
    #        ('C', (7.0, 7.0, 7.0), 0.4),
    #        ('H', (8.0, 8.0, 8.0), 0.8),
    #        ('Au',(9.0, 9.0, 9.0), 1.2),
    #]
    
    #bm = bmesh.new()
    
    u_seg = 3
    v_seg = 2
    num_verts = u_seg + v_seg
    

    count = 0
    for i in xyz[0]:
        bm = bmesh.new()
        for o in bpy.data.objects:
            if o.name.startswith(name):
                obj=o
        bpy.context.view_layer.objects.active = obj
        bpy.ops.object.mode_set(mode='EDIT')
        bm.from_mesh(sphereMesh)
        bpy.ops.object.mode_set(mode='OBJECT')
        
        
        locMatrix = mathutils.Matrix.Translation( (i[0], i[1], i[2]) )
        scaleMatrix = mathutils.Matrix.Scale(radius, 4)
        mesh = bmesh.ops.create_uvsphere(bm, u_segments=u_seg, v_segments=v_seg,
                        diameter=1.0, matrix=locMatrix @ scaleMatrix)
        #print(mesh)
        if count % 10 == 0:
            print('Particles {} out of {} [{:3.0f}% ]'.format(count,xyz[0].shape[0],100*count/xyz[0].shape[0]),end='\r')
        count+=1
        
        bm.to_mesh(sphereMesh)
        bm.free()
        
        
    bm.free()
    return num_verts

