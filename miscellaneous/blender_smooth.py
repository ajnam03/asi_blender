import bpy
import bmesh
import os, sys, tempfile
import subprocess as sp
import csv
import numpy as np
import imp

'''This functions is a smoothing routine to smooth your object meshes using Gamer: 
https://gamer.readthedocs.io/en/latest/tutorials/tutorials/blendertutorial.html'''

#not_analyzed - list of not analyzed objects (e.g. due to errors in mesh)
#filter -filters object in object list to include
#ex_filter -"exclusion filter" filters objects in object list NOT to include
def smooth(filter=[''], ex_filter=['foo']):
    print('smoothing')
    collection_name = 'Collection' #define name of collection
    obj_name_list = [obj.name for obj in bpy.data.collections[collection_name].all_objects if any(f1 in obj.name for f1 in filter) and any(f2 not in obj.name for f2 in ex_filter)]
    print('tot_num_objects = ', len(obj_name_list))
    non_smoothed = []
    i = 0
    for o in obj_name_list:
        print('obj name = ', o)
        print('i = ', i, 'out of tot objs = ', len(obj_name_list))
        i += 1
        #blender/selection settings
        bpy.ops.object.mode_set(mode ='OBJECT') #object mode
        obj = bpy.data.objects[o]
        bpy.context.view_layer.objects.active = obj #make active obj
        #selection setting
        bpy.ops.object.select_all(action='DESELECT') #deselect all objs
        obj.select_set(True)#select only object 'obj'
        
        print('smoothing meshes')

        #gamer version available
        gamer_version = bpy.context.scene.gamer.gamer_version[1]
        if gamer_version == '2':
            # We have found GAMer version 2
            gamer_smiprops = bpy.context.scene.gamer.surfmesh_improvement_properties
            gamer_smiprops.dense_rate = 2.5
            gamer_smiprops.dense_iter = 1
            gamer_smiprops.smooth_iter = 10
            gamer_smiprops.preserve_ridges = True   
        else:
            # We have found GAMer version 1
            gamer_mip = bpy.context.scene.gamer.mesh_improve_panel
            gamer_mip.dense_rate = 2.5
            gamer_mip.dense_iter = 1
            gamer_mip.max_min_angle = 20.0
            gamer_mip.smooth_iter = 10
            gamer_mip.preserve_ridges = True

        bpy.ops.object.mode_set(mode='EDIT')

        try:
            bpy.ops.mesh.beautify_fill(angle_limit=1.57)
            bpy.ops.mesh.subdivide(number_cuts=2)
            bpy.ops.object.mode_set(mode='OBJECT')
            
            #smoothing 
            bpy.ops.npt_gamer.smooth('INVOKE_DEFAULT')
            bpy.ops.npt_gamer.coarse_dense('INVOKE_DEFAULT')
            bpy.ops.npt_gamer.smooth('INVOKE_DEFAULT')
            bpy.ops.npt_gamer.coarse_dense('INVOKE_DEFAULT')
            bpy.ops.npt_gamer.smooth('INVOKE_DEFAULT')
            bpy.ops.npt_gamer.coarse_dense('INVOKE_DEFAULT')
            bpy.ops.npt_gamer.smooth('INVOKE_DEFAULT')
            bpy.ops.npt_gamer.coarse_dense('INVOKE_DEFAULT')
            bpy.ops.npt_gamer.smooth('INVOKE_DEFAULT')
            bpy.ops.npt_gamer.normal_smooth('INVOKE_DEFAULT')
            bpy.ops.npt_gamer.normal_smooth('INVOKE_DEFAULT')
            bpy.ops.npt_gamer.normal_smooth('INVOKE_DEFAULT')
            bpy.ops.npt_gamer.normal_smooth('INVOKE_DEFAULT')
        
        except RuntimeError:
            print('run time error')
            non_smoothed.append(o)
            obj.name = obj.name + 'nonSmooth'
            obj.data.name = obj.data.name + 'nonSmooth'
            continue

    not_analyzed = non_smoothed #any of the synapses that won't be analyzed
    not_analyzed = [i for i in not_analyzed if 'c' in i.rstrip('abc')] #drop any c-objects (it's ok if not cobject wasn't smoothed)

    return not_analyzed

