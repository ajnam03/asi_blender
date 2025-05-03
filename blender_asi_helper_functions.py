import bpy
import bmesh
import math
from itertools import compress
from mathutils import Vector
import re

def extract_face_density(obj):
    bpy.ops.object.mode_set(mode='EDIT')
    #create bmesh
    bm = bmesh.new()
    bm.from_mesh(obj.data)
    #find num faces
    num_faces = len([f for f in bm.faces])
    #find volume
    vol = bm.calc_volume()
    density = num_faces /vol

    bm.free()
    return density

def min_face_density(col_name, filter, exclude_filter, min_density):
    obj_name_list = [obj.name for obj in bpy.data.collections[col_name].all_objects if any(j in obj.name for j in filter)]
    obj_name_list = [o for o in obj_name_list if any(j not in o for j in exclude_filter)]
    for o in obj_name_list:
        print('min face density object = ', o)
        obj = bpy.data.objects[o]
        density = extract_face_density(obj) #extract face density

        print('original density = ', density)
        while density < min_density: #subdivide if less than face density
            print('face density new = ', density)
            #blender settings
            bpy.ops.object.mode_set(mode='OBJECT')
            bpy.ops.object.select_all(action='DESELECT')
            obj.select_set(True)
            bpy.context.view_layer.objects.active = obj #set o to active object
            bpy.ops.object.mode_set(mode='EDIT')
            bpy.ops.mesh.subdivide(number_cuts=1)
            density = extract_face_density(obj) #update face density
            #print('updated face density = ', density)
        #print('final face density = ', density)

def check_obj_name(name): 
    if (not bpy.data.objects.get(name)) & name[-1].isalpha():
        name = name.rstrip('abcd')
    return name
        
#bounding box method
def extract_obj_center_bbox(obj_name):
    obj = bpy.data.objects[obj_name]
    #sum vectors from each corner of bbox and divide by # of corners (i.e. multiply 0.125)
    local_bbox_center = 0.125 * sum((Vector(b) for b in obj.bound_box), Vector())
    #convert to world coords
    global_bbox_center = obj.matrix_world @ local_bbox_center
    return local_bbox_center, global_bbox_center

def create_bmesh(obj):
    bpy.context.view_layer.objects.active = obj #make active obj
    bpy.ops.object.mode_set(mode ='EDIT') #change to edit mode
    bm = bmesh.new()
    bm.from_mesh(obj.data)  
    bm_faces = [f for f in bm.faces]
    bm_edges = [e for e in bm.edges]
    return bm, bm_faces, bm_edges

def check_obj_intersect(faulty_seg, valid_seg): #check two objects intersects (i.e. boolean difference = 0)
    valid = bpy.data.objects[valid_seg]
    faulty = bpy.data.objects[faulty_seg]
    bpy.context.view_layer.objects.active=faulty
    bpy.ops.object.mode_set(mode='OBJECT')
    mod1 = faulty.modifiers.new('bool_diff', 'BOOLEAN')
    mod1.operation = 'DIFFERENCE'
    mod1.object = valid

    #bmesh with modifiers applied (in object mode)
    depsgraph = bpy.context.evaluated_depsgraph_get()
    obj_eval = faulty.evaluated_get(depsgraph)
    me = obj_eval.to_mesh()
    bm = bmesh.new()
    bm.from_mesh(me)
    obj_eval.to_mesh_clear()
    volume = bm.calc_volume()
    faulty.modifiers.remove(mod1) #delete the modifier on faulty

    #check volume of difference boolean
    intersect=False
    if volume == 0.0:
        #apply the modifier BUT to the valid seg and delete the faulty 
        bpy.context.view_layer.objects.active = valid
        bpy.ops.object.mode_set(mode='OBJECT')
        mod2 = valid.modifiers.new('bool_union', 'BOOLEAN')
        mod2.operation='UNION'
        mod2.object =faulty
        bpy.ops.object.modifier_apply(modifier=mod2.name) #apply modifier to merge faulty with valid
        print('merged w/ valid seg = ', valid_seg)    
        intersect = True
    return intersect


#fixes if there are sharp edge selections that would create overlapping boundaryLoops
def fix_edgeSelection(obj): 
    
    #select any faces that are neighbored by two non selected faces-> creates problems w/ boundaryLoop
    bm, bm_faces, bm_edges = create_bmesh(obj)
#    bm.select_flush(True)
    selected_faces = [] #find all selected faces
    for f in bm_faces:
        if f.select == True:
            selected_faces.append(f.index)
    
    unselected_neighbor = []
    for f_index in selected_faces: #list of non-selected faces adjacent to selected
        linked_faces = set()
        for e in bm_faces[f_index].edges:
            linked_faces_e = [linked.index for linked in e.link_faces if linked.select == False]
            linked_faces.update(set(linked_faces_e))
        unselected_neighbor = unselected_neighbor + list(linked_faces)

    add_selection = list(set([i for i in unselected_neighbor if unselected_neighbor.count(i) == 2]))
    new_selection = selected_faces + add_selection
    for f in bm_faces:
        if f.index in new_selection:
            f.select = True
        else:
            f.select = False

    bpy.ops.object.mode_set(mode='OBJECT')
    bm.to_mesh(obj.data)
    bm.free()
    bpy.ops.object.mode_set(mode='EDIT')

    #first select more/less
    bpy.context.view_layer.objects.active = obj
    bpy.ops.object.mode_set(mode ='EDIT') 
    bpy.ops.mesh.select_more() #select more
    bpy.ops.mesh.select_less() #select less

def getLinkedEdges(bm_edges, edge_index, selected_edge_verts, selected_edges):
    
    #get verts
    edge_verts = [v.index for v in bm_edges[edge_index].verts]
    vert_filter = [set(edge_verts).intersection(set(i)) for i in selected_edge_verts] 
    
    linked_edges = list(compress(selected_edges, vert_filter))
    return linked_edges
    
def getBoundaryLoop(obj):
    bpy.ops.mesh.region_to_loop()
    bpy.ops.object.mode_set(mode='EDIT')
    bm = bmesh.from_edit_mesh(obj.data)
    bm_edges = [e for e in bm.edges]
    
    #get selected edges and verts
    selected_edges = []
    for e in bm_edges:
        if e.select == True:
            selected_edges.append(e.index)
    selected_edge_verts = [[v.index for v in bm_edges[e].verts] for e in selected_edges]
    
    linkedEdges = [getLinkedEdges(bm_edges, e, selected_edge_verts, selected_edges) for e in selected_edges]
    
    #construct list of all linked selected edges
    pooled = [set(subList) for subList in linkedEdges]
    merging = True
    while merging:
        merging=False
        for i,group in enumerate(pooled):
            merged = next((g for g in pooled[i+1:] if g.intersection(group)),None)
            if not merged: continue
            group.update(merged)
            pooled.remove(merged)
            merging = True
    bm.free()
    return pooled

def collapse_nested_bL_helper(obj, bl): #returns inner faces associated with that bL
    bpy.ops.mesh.select_all(action='DESELECT')
    bm, bm_faces, bm_edges = create_bmesh(obj)
    for f in bm_faces:
        f.select = False
    #select bL
    bpy.ops.mesh.select_mode(type='EDGE')
    for e in bm_edges:
        if e.index in bl:
            e.select = True
        else:
            e.select = False
    
    #visualize changes
    bpy.ops.object.mode_set(mode='OBJECT')
    bm.to_mesh(obj.data)
    bpy.ops.object.mode_set(mode='EDIT')
    bm.free()
    
    #get inner faces
    bpy.ops.mesh.loop_to_region()
    
    bl_innerFaces = []
    bm = bmesh.from_edit_mesh(obj.data)
    
    for f in bm.faces:
        if f.select == True:
            bl_innerFaces.append(f.index)
    
    bm.free()
    return bl_innerFaces

#collapse multiple boundary nested or outlier boundaryloops into the single main asi boundary loop
def collapse_nested_bL(obj, boundaryLoop_all):
    bl_innerFaces_all = []
    for b in boundaryLoop_all:
        bl_innerFaces = collapse_nested_bL_helper(obj, b)
        bl_innerFaces_all.append(bl_innerFaces)
        
    #sort bL_innerFaces_all list in reverse order based on list length (longest face list -> shortest)
    bl_innerFaces_all = sorted(bl_innerFaces_all, key=lambda l: len(l), reverse=True)
    
    #merge any overlapping innerFace sets        
    overlap_all = []
    overlap_all.append(set(bl_innerFaces_all[0])) #append first 
    for i in range(1, len(bl_innerFaces_all)):
        intersection = False
        for j in range(len(overlap_all)):
            if set(bl_innerFaces_all[i]).intersection(overlap_all[j]): #check if there is an intersection w/ existing face list
                overlap_all[j].update(set(bl_innerFaces_all[i]))
                intersection = True
        if not intersection:
            overlap_all.append(set()) #otherwise, append a new empty set
            overlap_all[-1].update(bl_innerFaces_all[i]) #add to new set

    #keep only the largest innerFace list
    #find max innerFace list
    bl_innerFaces_new_lens = [len(face_list) for face_list in overlap_all]
    max_bl_index = bl_innerFaces_new_lens.index(max(bl_innerFaces_new_lens)) #find index of longest face list
    bl_innerFaces_new = list(overlap_all[max_bl_index])

    return bl_innerFaces_new

def calc_loop_perimeter(edges, loop_idx_list, conversion):
    if loop_idx_list:
        loop_edge_lens = [e.calc_length() for e in edges if e.index in loop_idx_list] #length of each edge in loop idx list
        loopLength_convert = sum(loop_edge_lens) * conversion

        return loopLength_convert

#find midpoint of an edge
def get_edge_midpoint(edge, matrix_world):
    edge_v = [v for v in edge.verts]
    edge_v_midpoint_x = (edge_v[0].co[0] + edge_v[1].co[0]) * 0.5
    edge_v_midpoint_y = (edge_v[0].co[1] + edge_v[1].co[1]) * 0.5
    edge_v_midpoint_z = (edge_v[0].co[2] + edge_v[1].co[2]) * 0.5
    edge_v_midpoint_loc = Vector([edge_v_midpoint_x, edge_v_midpoint_y, edge_v_midpoint_z]) #calculate edge midpoint
    edge_v_midpoint_world = matrix_world @ edge_v_midpoint_loc

    return [edge_v_midpoint_loc, edge_v_midpoint_world]

#yes/no within a particular distance threshold
def extract_glia_appo(edge_midpoint, ax_matrix_world, g, glia_appo, conversion, contour:bool=False, mesh:bool=False):
    #get midpoint of the edge
    # edge_v_midpoint_loc, edge_v_midpoint_world = get_edge_midpoint(edge, ax_matrix_world)
    # print('calc distance')
    if contour: #if using glia contours
        astro_appo = False
        for astro in g:
            glia_dist = [math.dist(edge_midpoint[1], astro.matrix_world @ v.co) for v in astro.data.vertices]
            if min(glia_dist) < glia_appo:
                astro_appo = True
                return astro_appo
        return astro_appo

    if mesh: #if using glia mesh
        astro_appo = False
        # g_mesh = [i for i in g if 'contour' not in i.name]
        # g_contour = [j for j in g if 'contours' in j.name][0]
        g_mesh = g
        for astro in g_mesh:
            #find closest point on glia mesh --> uncomment if using glia mesh instead of contours
            bpy.ops.object.mode_set(mode ='OBJECT')
            g_hit, g_loc, g_norm, g_idx = astro.closest_point_on_mesh(edge_midpoint[0])
            # convert closest glia point into global coordinates
            glia_loc_world = astro.matrix_world @ g_loc
            if math.dist(glia_loc_world, edge_midpoint[1]) < glia_appo: #math.dist doesn't need to be converted
                astro_appo = True
                return astro_appo
        return astro_appo

#extract min appo dist
def extract_min_glia_appo(edge_midpoint, g, contour:bool=False, mesh:bool=False):
    if contour: #if using glia contours
        astro_appo_dist = []
        for astro in g:
            glia_dist = [math.dist(edge_midpoint[1], astro.matrix_world @ v.co) for v in astro.data.vertices]
            astro_appo_dist = astro_appo_dist + glia_dist 
        return min(astro_appo_dist)

    if mesh: #if using glia mesh
        astro_appo_dist = []
        g_mesh = g
        for astro in g_mesh:
            #find closest point on glia mesh --> uncomment if using glia mesh instead of contours
            bpy.ops.object.mode_set(mode ='OBJECT')
            g_hit, g_loc, g_norm, g_idx = astro.closest_point_on_mesh(edge_midpoint[0]) #edge_midpoint[0] is the local coordinates
            # convert closest glia point into global coordinates
            glia_loc_world = astro.matrix_world @ g_loc
            astro_appo_dist.append(math.dist(glia_loc_world, edge_midpoint[1]))
        return min(astro_appo_dist)

def fix_obj_name(col_name, filter_str, remove_str):
    #obj list
    obj_name_list = [obj.name for obj in bpy.data.collections[col_name].all_objects if filter_str in obj.name]
    for o in obj_name_list:
        print('o = ', o)
        obj = bpy.data.objects[o]
        try:
            new_name =obj.name.replace(remove_str + o.split(remove_str)[1], '')
        except IndexError:
            new_name = obj.name
        print('new name = ', new_name)
        obj.name = new_name
        obj.data.name = new_name

#general function to delete faulty meshes give obj name
def del_faulty_meshes_general(obj_name):
    print('obj_name = ', obj_name)
    obj = bpy.data.objects[obj_name]
    #blender settings
    bpy.ops.object.select_all(action='DESELECT')
    obj.select_set(True)
    #delete faulty mesh
    bpy.ops.object.delete()   

#function to delete faulty meshes with certain filterable name
def del_faulty_meshes_filter(col_name, delete_filter):
    # bpy.ops.object.mode_set(mode='OBJECT')
    #obj list
    obj_name_list = [obj.name for obj in bpy.data.collections[col_name].all_objects if all(j in obj.name for j in delete_filter)]
    
    for o in obj_name_list:
        print('o = ', o)
        try:
            del_faulty_meshes_general(o)
        except RuntimeError:
            bpy.ops.object.editmode_toggle()
        
def fix_obj_name(col_name, filter_str, delete_str):
    #obj list
    obj_name_list = [obj.name for obj in bpy.data.collections[col_name].all_objects if all(j in obj.name for j in filter_str)]
    obj_name_list = [o for o in obj_name_list if 'seg' not in o] #don't include seg objects
    #obj_name_list = [obj.name for obj in bpy.data.collections[col_name].all_objects]
    for o in obj_name_list:
        print('o = ', o)
        obj = bpy.data.objects[o]
        try:
            new_name =obj.name.replace(delete_str + o.split(delete_str)[1], '')
        except IndexError:
            new_name = obj.name
        print('new name = ', new_name)
        obj.name = new_name
        obj.data.name = new_name
        matnum = re.search("\d\d\d$", obj.name)
        if matnum:
            obj.name = (obj.name.rstrip(obj.name[-4:]))
            obj.data.name = (obj.name.rstrip(obj.name[-4:]))

#edges are edge references not edge idxs
def find_connected_segments(edges):
    visited = set()
    segments = []
    
    for edge in edges:
        if edge not in visited:
            # Start a new segment
            segment = []
            stack = [edge]
            
            while stack:
                current_edge = stack.pop()
                if current_edge not in visited:
                    visited.add(current_edge)
                    segment.append(current_edge)
                    
                    # Find connected edges
                    for vert in current_edge.verts:
                        for linked_edge in vert.link_edges:
                            if linked_edge in edges and linked_edge not in visited:
                                stack.append(linked_edge)
            
            segments.append(segment)
    return segments

#extract geometric centroid
def extract_local_geoCentroid(bm):
    
    #convert to edit mode
    bpy.ops.object.mode_set(mode='EDIT')

    #vert coordinates
    vertices = [v.co for v in bm.verts]

    # Calculate the geometric centroid in local coordinates
    if vertices:
        centroid_local = sum(vertices, Vector()) / len(vertices)
        return centroid_local
