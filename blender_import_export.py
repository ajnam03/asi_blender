import bpy
import bmesh
import os, sys, tempfile
import subprocess as sp
import csv
import numpy as np
import imp


def export_objects(temp_dir, filter=[''], ex_filter =['foo'], not_analyzed = []):
    '''EXPORT'''
    #export objects
    collection_name = 'Collection' #define name of collection
    obj_name_list = [obj.name for obj in bpy.data.collections[collection_name].all_objects if any(o in obj.name for o in filter) and all(ef not in obj.name for ef in ex_filter)]
    # temp_dir = tempfile.TemporaryDirectory(dir='/tmp')
    # temp_dir = '/Volumes/LaCie/Harris_lab/brain_init23/tmp/'
    print('temp dir = ', temp_dir)
    for o in obj_name_list:
        if o not in not_analyzed: 
            print('obj name = ', o)
            #blender/selection settings

            obj = bpy.data.objects[o]
            bpy.context.view_layer.objects.active = obj #make active obj
            bpy.ops.object.mode_set(mode ='OBJECT') #object mode
            #selection setting
            bpy.ops.object.select_all(action='DESELECT') #deselect all objs
            obj.select_set(True)#select only object 'obj'

            #save obj to temp directory
            print('exporting obj to temp directory') 
            # filepath = temp_dir.name + '/f"{o}.ply'
            filepath = temp_dir + '/' + o + '.ply'
            bpy.ops.export_mesh.ply(filepath=str(filepath), use_selection=True)
        else:
            print('Synapse in not-analyzed list = ', o)

def import_objects(temp_dir, filter=['']):
    '''IMPORT'''
    #import cgal output meshes back into blender
    import_filepaths = os.listdir(temp_dir)
    # cgal_filter = ['_smooth_seg', '_remesh']
    import_filepaths = [i for i in import_filepaths if any(out in i for out in filter)]

    #make 'Collection' collection active
    col = bpy.context.view_layer.layer_collection.children['Collection']
    bpy.context.view_layer.active_layer_collection = col
    for mesh in import_filepaths:
        mesh_path = temp_dir + '/' + mesh        
        bpy.ops.import_mesh.ply(filepath=str(mesh_path))
        

def obj_move_collection(current_collection, new_collection_name, filter=[''], ex_filter=['foo']): #move objects to new collection
    #get sp, ax objects names (keep c objects)
    collection_name = current_collection
    obj_list = [obj for obj in bpy.data.collections[collection_name].all_objects if any(f in obj.name for f in filter) and all(ef not in obj.name for ef in ex_filter)]

    # collection_list =[c for c in bpy.data.collections]
    try: 
        col = bpy.data.collections[new_collection_name]
    except KeyError:
        col = bpy.data.collections.new(new_collection_name)
        bpy.context.scene.collection.children.link(col)
    for o in obj_list:
        # bpy.data.collections['Collection'].objects.unlink(o)
        bpy.data.collections[current_collection].objects.unlink(o)
        bpy.data.collections[new_collection_name].objects.link(o)



    

