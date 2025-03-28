import bpy
import bmesh
import math
from itertools import compress
from mathutils import Vector
import numpy as np
import csv
import sys
import blender_asi_helper_functions as ahf
import importlib

"""MAKE SURE YOU HAVE LOOP TOOLS"""

class asi:

    def __init__(self, c_name, glia_appo, asi_appo = 0.045):
        self.c_name = c_name
        self.sp_name = ahf.check_obj_name(c_name.split('c')[0] + 'sp' + c_name.split('c')[1].split('_')[0].rstrip('abcd')+'_remesh')
        self.ax_name = ahf.check_obj_name(c_name.split('c', 1)[0] + 'ax' + c_name.split('c', 1)[1])
        self.g_name = [g.name for g in bpy.data.collections['Collection'].objects if 'astro' in g.name]
        self.scale_name = 'scale'
        self.ref_name = '1um_rad_ref' #1um rad reference circle imported from PyReconstruct
        self.asi_appo = asi_appo
        self.glia_appo = glia_appo
    
    def check_objects_exist(self):
        c = bpy.data.objects[self.c_name]
        sp = bpy.data.objects[self.sp_name]
        ax = bpy.data.objects[self.ax_name]

    def calc_unit_conversion(self): 
        ref = bpy.data.objects[self.ref_name]
        self.ref_rad = (ref.dimensions.x / 2) 
        ref_bm, ref_bm_faces, ref_bm_edges = ahf.create_bmesh(ref)
        self.ref_area = (sum(f.calc_area() for f in ref_bm_faces) / 2) 

        self.conversion_length = 1/self.ref_rad

    #duplicate axon so don't affect original axon mesh
    def duplicateAxon(self):
        ax = bpy.data.objects[self.ax_name]
        ax_copy = ax.copy()
        ax_copy.data = ax.data.copy()
        ax_copy.name = self.ax_name + '_copy'
        ax_copy.data.name = self.ax_name + '_copy'
        bpy.data.collections["Collection"].objects.link(ax_copy)
        self.ax_name = self.ax_name + '_copy' #update self.ax_name so you only edit the copy

    def extract_syn_center(self):
        bpy.ops.object.mode_set(mode ='OBJECT')
        #get tri syn components

        # c = bpy.data.objects[self.c_name] #c object
        sp = bpy.data.objects[self.sp_name] #spine object
        ax = bpy.data.objects[self.ax_name] #axon object

        #synapse center -> c object center
        self.c_center_local, self.c_center_global = ahf.extract_obj_center_bbox(self.c_name) 
  
        #closest axon point to syn center
        bpy.context.view_layer.objects.active = ax
        ax_hit, ax_loc, ax_norm, ax_idx = ax.closest_point_on_mesh(self.c_center_local)
        self.ax_syn_idx = ax_idx #id of axon face closest to c_object
        ax_loc_global = ax.matrix_world @ ax_loc

        #closest spine point to syn center
        sp_hit, sp_loc, sp_norm, sp_idx = sp.closest_point_on_mesh(self.c_center_local)
        self.sp_syn_idx = sp_idx #id of spine face closest to c_object
        sp_loc_global = sp.matrix_world @ sp_loc

        self.syn_dist = math.dist(sp_loc_global, ax_loc_global) #synapse distance
        self.syn_dist_convert = self.syn_dist * self.conversion_length

    #segment spine and axon w/ cgal segments    
    def segmentSA(self):
        #check if segments exist
        segment = True
        try: 
            sp_seg0 = bpy.data.objects[self.sp_name + '_seg0']
        except KeyError:
            #returns a KeyError if seg0 doesn't exists (i.e. cgal program was not able to generate any segments)
            #only need to check _seg0 because ig segmentation was successful seg0 and seg1 were generated in all cases
            segment = False 

        if segment:
            #extract spine head segment -> should be largest sgment + closest to c-object otherwise do not segment
            self.sp_head_seg = ahf.extract_spineHead_seg(self.sp_name, self.c_name, self.c_center_global) #empty list if both conditions weren't met
            if self.sp_head_seg: 
                self.sp_name = ahf.segmentSpineHead(self.c_center_local, self.c_center_global, self.sp_name, self.sp_head_seg[0]) #segment + update self.sp_name to be used for asi analysis
 
    #extract asi face indices using axon face normals 
    def extract_asi_idx(self):
        asi_idx = []
        ax = bpy.data.objects[self.ax_name]
        sp = bpy.data.objects[self.sp_name]
        sp_bm, sp_bm_faces, sp_bm_edges = ahf.create_bmesh(sp)

        bpy.ops.mesh.select_mode(type="VERT")  ### NEED VERT MODE FOR ASI SELECTION
        # bpy.ops.object.mode_set(mode ='OBJECT')
        sp_bm.faces.ensure_lookup_table()
        #define asi by if spine face normal hits axon by ray casting w/ max dist = avg syn_dist * factor 
        for f1 in sp_bm_faces: 
            # ax_hit, ax_loc, ax_norm, ax_idx = ax.ray_cast(f.calc_center_median(), f.normal, distance=(self.syn_dist * asi.factor))
            ax_hit, ax_loc, ax_norm, ax_idx = ax.ray_cast(f1.calc_center_median(), f1.normal, distance=self.asi_appo) #0.045 default
            if ax_hit: 
                asi_idx.append(ax_idx)
    
        #remove dupliaces
        self.asi_idx = list(set(asi_idx)).copy()

    def select_asi_faces(self):
        ax = bpy.data.objects[self.ax_name]
        ax_bm, ax_bm_faces, ax_bm_edges = ahf.create_bmesh(ax)
        ax_bm.faces.ensure_lookup_table()
        for f in ax_bm.faces:
            if f.index in self.asi_idx:
                f.select = True
            else:
                f.select = False

        #update mesh and visualize changes
        bpy.ops.object.mode_set(mode ='OBJECT')
        ax_bm.to_mesh(ax.data)
        ax_bm.free()
            
    #clean asi by removing any outlier clusters and filling any selection holes
    def cleanUpASI(self): 
        #ax obj
        ax = bpy.data.objects[self.ax_name]

        #select more/less to improve selection holes
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_more()
        bpy.ops.mesh.select_more()
        bpy.ops.mesh.select_less()
        bpy.ops.mesh.select_less()

        #fixes if there are sharp edge selections that would create overlapping boundaryLoops
        print('fix edge selection')
        ahf.fix_edgeSelection(ax)

        #list of all boundaryloops from current pre-cleaned asi selection
        print('get all boundaryLoops')
        boundaryLoop_all = ahf.getBoundaryLoop(ax)

        #collapse multiple boundary nested or outlier boundaryloops into the single main asi boundary loop
        #record inner faces of this single boundaryloop 
        self.asi_idx = ahf.collapse_nested_bL(ax, boundaryLoop_all)

    def calc_asi_area(self):
        bpy.ops.object.mode_set(mode='EDIT')
        ax = bpy.data.objects[self.ax_name]
        ax_bm, ax_bm_faces, ax_bm_edges = ahf.create_bmesh(ax)
        ax_bm.faces.ensure_lookup_table()

        self.sp_area = sum([ax_bm_faces[i].calc_area() for i in range(len(ax_bm_faces))]) #to check units against PyReconstruct
        self.asi_area = 0
        for i in range(len(ax_bm_faces)):
            if ax_bm_faces[i].select == True:
                self.asi_area += ax_bm_faces[i].calc_area()
        

    #convert asi face selection to boundary loop -> calc perimeter
    def asi_boundaryLoop(self):
        bpy.ops.object.mode_set(mode='EDIT')
        ax = bpy.data.objects[self.ax_name]

        # ax_bm_asi = bmesh.from_edit_mesh(ax.data)
        # ax_bm_edges = [e for e in ax_bm_asi.edges]
        bpy.ops.mesh.select_mode(type="EDGE")

        #convert currently selected asi faces to boundary loop -> flatten/relax loop
        bpy.ops.mesh.region_to_loop() #boundary loop
        for i in range(10):
            bpy.ops.mesh.looptools_relax()  #relax fxn from looptools (smooths loop) x10

        # ax_bm, ax_bm_faces, ax_bm_edges = ahf.create_bmesh(ax)
        ax_bm = bmesh.from_edit_mesh(ax.data)
        ax_bm.edges.ensure_lookup_table()
        
        #add selected asi edges as attribute
        asi_bL_idx = [e.index for e in ax_bm.edges if e.select == True]
        self.asi_bL_idx =  asi_bL_idx.copy()#find selected -> add as asi boundary loop attribute

        bmesh.update_edit_mesh(ax.data)
        ax_bm.free()

    #select edges and show asi boundary loop w/ glia apposition
    def asi_glia_boundaryLoop(self, contour:bool =False, mesh:bool=False):

        #define bmeshes
        ax = bpy.data.objects[self.ax_name]
        g = [bpy.data.objects[g] for g in self.g_name]
        
        # ax_bm = bmesh.from_edit_mesh(ax.data)
        ax_bm, ax_bm_faces, ax_bm_edges = ahf.create_bmesh(ax)
        ax_bm.edges.ensure_lookup_table()

        #find midpoint of an edge
        asi_midpoints = [ahf.get_edge_midpoint(ax_bm.edges[edge], ax.matrix_world) for edge in self.asi_bL_idx]
        
        #extract asi_bL_idx edges that are apposed by glia 
        #@glia appo threshold = glia_appo[i]
        self.asiGlia_bL = list()
        for ig, glia_appo in enumerate(self.glia_appo):
            self.asiGlia_bL.append([self.asi_bL_idx[i] for i in range(len(self.asi_bL_idx)) if ahf.extract_glia_appo(asi_midpoints[i], 
                                                                                                                     ax.matrix_world, g, glia_appo, 
                                                                                                                     self.conversion_length, contour=contour, mesh=mesh)])
        
        bpy.ops.object.mode_set(mode ='EDIT')
        bmesh.update_edit_mesh(ax.data)
        ax_bm.free()

    #the blender e.calc_length() function always calculates in meters no matter what blender units are set to (at least for 3.6.5)
    def calc_peri_frac(self, asi:bool = False, asi_glia:bool=False):
        #ax bmesh
        ax = bpy.data.objects[self.ax_name]
        ax_bm, ax_bm_faces, ax_bm_edges = ahf.create_bmesh(ax)
        if asi: 
            self.asi_peri = ahf.calc_loop_perimeter(ax_bm_edges, self.asi_bL_idx, self.conversion_length) #perimeter of boundary loop
        if asi_glia:
            #glia_appo[0]
            self.asiGlia_peri = list()
            self.asiGlia_frac = list()
            for i, j in enumerate(self.glia_appo):
                if len(self.asiGlia_bL[i]) > 0:
                    self.asiGlia_peri.append(ahf.calc_loop_perimeter(ax_bm_edges, self.asiGlia_bL[i], self.conversion_length)) #perimeter of glia appo boundary loop
                    self.asiGlia_frac.append(self.asiGlia_peri[i]/self.asi_peri)
                else:
                    self.asiGlia_peri.append(0) 
                    self.asiGlia_frac.append(0)
       
    
    def select_bL(self, asi:bool=False, asi_glia:bool=False, appo_dist=0.01):
        #ax bmesh
        ax = bpy.data.objects[self.ax_name]
        ax_bm = bmesh.from_edit_mesh(ax.data)
        ax_bm.edges.ensure_lookup_table()

        #select asi boundary loop only
        if asi:
            for e in ax_bm.edges:
                if e.index in self.asi_idx:
                    e.select = True
                else:
                    e.select = False
        if asi_glia:
            idx = self.glia_appo.index(appo_dist)
            for e in ax_bm.edges:
                if e.index in self.asiGlia_bL[idx]:
                    e.select = True
                else:
                    e.select = False

        # update mesh and visualize changes
        bmesh.update_edit_mesh(ax.data)
        ax_bm.free()

    def extract_asi_idx_asiFix(self):
        asi_idx = []
        ax = bpy.data.objects[self.ax_name]
        ax_bm_asiFix = bmesh.new()
        ax_bm_asiFix.from_mesh(ax.data)

        ax_bm_asiFix_faces = [f for f in ax_bm_asiFix.faces]
        for f in ax_bm_asiFix_faces:
            if f.select == True:
                asi_idx.append(f.index)
        ax_bm_asiFix.free()
        self.asi_idx = list(set(asi_idx)).copy()

    def extract_asiPeri_idx(self):
        asi_bL_idx = []
        ax = bpy.data.objects[self.ax_name]

        #define bmeshs
        ax_bm = bmesh.new()
        ax_bm.from_mesh(ax.data)

        #extract edges
        ax_bm_edges = [e for e in ax_bm.edges]

        #add to asi_bL_idx list if edge is selected
        for e in ax_bm_edges:
            if e.select == True:
                asi_bL_idx.append(e.index)

        ax_bm.free()
        self.asi_bL_idx = list(set(asi_bL_idx)).copy()
    
    "PROBLEM WITH THIS FUNCTION!"
    #this might not be meaningful because some synapses had more glia traced so min dist could look farther out
    def extract_min_astro_appo(self, contour:bool=False, mesh:bool=False):

        #define bmeshes
        ax = bpy.data.objects[self.ax_name]
        g = [bpy.data.objects[g] for g in self.g_name]
        
        # ax_bm = bmesh.from_edit_mesh(ax.data)
        ax_bm, ax_bm_faces, ax_bm_edges = ahf.create_bmesh(ax)
        ax_bm.edges.ensure_lookup_table()

        #find midpoint of an edge
        asi_midpoints = [ahf.get_edge_midpoint(ax_bm.edges[edge], ax.matrix_world) for edge in self.asi_bL_idx]
        
        #min distance to astro for every asi peri vert
        asi_astro_appo_dist = [ahf.extract_min_glia_appo(i, g, contour=contour, mesh=mesh) for i in asi_midpoints]
        
        try:
            #min 
            min_astro_appo_dist = min(asi_astro_appo_dist)* self.conversion_length
        except ValueError:
            print('no nearest astro...')
            min_astro_appo_dist = np.nan

        self.min_astro_appo = min_astro_appo_dist
        ax_bm.free()

    def get_psd_obj(self):
        psd = bpy.data.objects[self.c_name.split('_remesh')[0]]

        psd_sp_name = self.c_name.split('_remesh')[0] + '_sp' 
        psd_ax_name = self.c_name.split('_remesh')[0] + '_ax'
        
        bpy.ops.object.select_all(action='DESELECT')
        psd.select_set(True)
        bpy.ops.object.duplicate()

        psd_copy = bpy.data.objects[self.c_name.split('_remesh')[0] + '.001']
        psd_copy.name = psd_ax_name
        psd_copy.data.name = psd_ax_name
        psd.name = psd_sp_name
        psd.data.name = psd_sp_name
        
        #get correct psd objects with updated names 
        psd_sp = bpy.data.objects[psd_sp_name]
        psd_ax = bpy.data.objects[psd_ax_name]

        bpy.context.view_layer.objects.active = psd_sp
        bpy.ops.object.modifier_add(type="SHRINKWRAP")
        shrinkwrap_mod = bpy.context.object.modifiers["Shrinkwrap"]
        shrinkwrap_mod.target = bpy.data.objects[self.sp_name]
        shrinkwrap_mod.wrap_method = 'NEAREST_SURFACEPOINT'
        shrinkwrap_mod.wrap_mode = 'ON_SURFACE'
        bpy.ops.object.modifier_move_to_index(modifier="Shrinkwrap", index=0)  # Move to the top
        bpy.ops.object.modifier_apply(modifier="Shrinkwrap")

        # Apply Shrinkwrap to the ax object
        bpy.context.view_layer.objects.active = psd_ax
        bpy.ops.object.modifier_add(type="SHRINKWRAP")
        shrinkwrap_mod = bpy.context.object.modifiers["Shrinkwrap"]
        shrinkwrap_mod.target = bpy.data.objects[self.ax_name]
        shrinkwrap_mod.wrap_method = 'NEAREST_SURFACEPOINT'
        shrinkwrap_mod.wrap_mode = 'ON_SURFACE'
        bpy.ops.object.modifier_move_to_index(modifier="Shrinkwrap", index=0)  # Move to the top
        bpy.ops.object.modifier_apply(modifier="Shrinkwrap")

   
    def extract_psd_area_centroid(self): #should be based on psd_sp_name
        try:
            print(self.psd_sp_name)
        except AttributeError:
            self.psd_sp_name = self.c_name.split('_remesh')[0] + '_sp'
            self.psd_ax_name  = self.c_name.split('_remesh')[0] + '_ax'

        # Switch to Object Mode if necessary
        if bpy.context.object.mode != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')
        
        psd_sp = bpy.data.objects[self.psd_sp_name]
        psd_ax = bpy.data.objects[self.psd_ax_name]
        
        #get bmesh
        psd_sp_bm, psd_sp_faces, psd_sp_edges = ahf.create_bmesh(psd_sp)
        psd_ax_bm, psd_ax_faces, psd_ax_edges = ahf.create_bmesh(psd_sp)
        psd_sp_bm.edges.ensure_lookup_table()
        psd_ax_bm.edges.ensure_lookup_table()
    
        # Get the surface area of psd shrinkwrapped to sp
        psd_area = sum(f.calc_area() for f in psd_sp_bm.faces)
        self.psd_area = psd_area * (self.conversion_length **2)
        
        #get psd centroid (in global coordinates)
        psd_sp_local_centroid = ahf.extract_local_geoCentroid(psd_sp_bm)
        psd_ax_local_centroid = ahf.extract_local_geoCentroid(psd_ax_bm)
        
        #convert to global coords
        self.psd_sp_centroid = psd_sp.matrix_world @ psd_sp_local_centroid
        self.psd_ax_centroid = psd_ax.matrix_world @ psd_ax_local_centroid

        #split into xyz coords
        self.psd_sp_centroid_x = self.psd_sp_centroid.x
        self.psd_sp_centroid_y = self.psd_sp_centroid.y
        self.psd_sp_centroid_z = self.psd_sp_centroid.z
        self.psd_ax_centroid_x = self.psd_ax_centroid.x
        self.psd_ax_centroid_y = self.psd_ax_centroid.y
        self.psd_ax_centroid_z = self.psd_ax_centroid.z

        psd_sp_bm.free()
        psd_ax_bm.free()

            #extract area of asi faces + centroid 
    def extact_asi_area_centroid(self):
        ax = bpy.data.objects[self.ax_name]

        #select asi peri
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.select_all(action='DESELECT') #deselect everything

        #temp mesh
        bpy.ops.object.mode_set(mode='EDIT')
        mesh = bmesh.from_edit_mesh(ax.data)
        mesh_edges = [e for e in mesh.edges]
        for idx in self.asi_bL_idx:
            mesh_edges[idx].select = True
        
        # Update the mesh
        bmesh.update_edit_mesh(ax.data)

        #switch back to Object Mode
        bpy.ops.object.mode_set(mode='OBJECT')
        bpy.context.view_layer.objects.active = ax
        bpy.ops.object.select_all(action='DESELECT')
        ax.select_set(True)

        #select inner faces
        bpy.ops.object.mode_set(mode='EDIT')
        bpy.ops.mesh.loop_to_region()

        asi_faces = []
        selected_verts = []
        mesh2 = bmesh.from_edit_mesh(ax.data)
        for f in mesh2.faces:
            if f.select ==True:
                asi_faces.append(f)
                for vert in f.verts: #also record asi face verts
                    selected_verts.append(vert.co)

        #calc asi area
        print(len(asi_faces))
        asi_area = sum(face.calc_area() for face in asi_faces)
        self.asi_area = asi_area * (self.conversion_length **2)

        print('length of asi bL = ', len(self.asi_bL_idx))
        print('asi area = ', self.asi_area)
        #calc centroid of asi
        if selected_verts:
            centroid = sum(selected_verts, Vector())/len(selected_verts)
        self.asi_centroid = ax.matrix_world @ centroid #convert to global coords
        self.asi_centroid_x = self.asi_centroid.x
        self.asi_centroid_y = self.asi_centroid.y
        self.asi_centroid_z = self.asi_centroid.z

        #select inner faces back to boundary loop
        bpy.ops.mesh.region_to_loop()

    #calc offset of psd center from asi center
    def extract_asi2psd_offset(self):
        psd_offset = math.dist(self.psd_ax_centroid, self.asi_centroid)
        self.psd_offset = psd_offset * self.conversion_length
        
    'THIS FUNCTION ASSUMES ONLY ONE GLIA NAME'
    #have to extract asiGlia_bL first for each appo Dist -> assume only one g_name!!! 
    def extract_asi2psd2glia_dist(self, contour:bool=False, mesh:bool=False):
        #define bmeshes
        ax = bpy.data.objects[self.ax_name] #SHOULD this be self.sp_name!!! --> SINCE THE PSD IS AGAINST THE POST SYN SIDE? OR Ax psd because it's asi
        psd = bpy.data.objects[self.psd_ax_name]
        ax_bm, ax_bm_faces, ax_bm_edges = ahf.create_bmesh(ax)
        ax_bm.edges.ensure_lookup_table()
        g = [bpy.data.objects[g] for g in self.g_name]

        #find midpoint of an edge
        asi_midpoints = [ahf.get_edge_midpoint(ax_bm.edges[edge], ax.matrix_world) for edge in self.asi_bL_idx]
        
        #for all asi vertices - get distance to closest point on psd (regardless of glia appo status)
        asi2psd_dist = [math.dist(i[0], psd.closest_point_on_mesh(i[0])[1]) for i in asi_midpoints]
        self.min_asi2psd_dist = min(asi2psd_dist) * self.conversion_length
        self.max_asi2psd_dist = max(asi2psd_dist) * self.conversion_length
        self.avg_asi2psd_dist = np.mean(asi2psd_dist) * self.conversion_length

        #for all asiGlia vertices - get distance to closest point on psd (for only asi peri vertices with glia within 120 nm)
        #last list of asiGlia_bL corresponds to vertices that had astro within 120 nm
        asiGlia_midpoints = [ahf.get_edge_midpoint(ax_bm.edges[edge], ax.matrix_world) for edge in self.asiGlia_bL[-1]]
        
        if(len(asiGlia_midpoints) >0):
            print('length of asiGlia modpoints = ', len(asiGlia_midpoints))
            #list of distances between asiGlia boundaryloop vertex and nearest point on glia
            asiGlia2glia_dist = [math.dist(i[0], g[0].closest_point_on_mesh(i[0])[1]) for i in asiGlia_midpoints]
            self.min_asiGlia2glia_dist = min(asiGlia2glia_dist) * self.conversion_length #min dist to glia 
            self.avg_asiGlia2glia_dist = np.mean(asiGlia2glia_dist) * self.conversion_length #avg dist to glia

            #list of distances between asiGlia boundaryloop vertex and nearest point on psd
            asiGlia2psd_dist = [math.dist(i[0], psd.closest_point_on_mesh(i[0])[1]) for i in asiGlia_midpoints] 
            self.min_asiGlia2psd_dist = min(asiGlia2psd_dist) * self.conversion_length #min dist to psd
            self.max_asiGlia2psd_dist = max(asiGlia2psd_dist) * self.conversion_length #max dist to psd
            self.avg_asiGlia2psd_dist = np.mean(asiGlia2psd_dist) * self.conversion_length #avg dist to psd

            #get for each asiGlia vertex that had glia within 120 nm - get ratio of distance to glia vs. distance to psd 
            asiGlia2glia_psd_distRatio = (np.array(asiGlia2glia_dist)/np.array(asiGlia2psd_dist)).tolist()
            self.min_glia2psd_distRatio = min(asiGlia2glia_psd_distRatio)
            self.max_glia2psd_distRatio = max(asiGlia2glia_psd_distRatio)
            self.avg_glia2psd_distRatio = np.mean(asiGlia2glia_psd_distRatio)

            #get avg dist from each asiGlia vertex to psd centroid
            asiGlia2psdCentroid_dist = [math.dist(i[1], self.psd_ax_centroid) for i in asiGlia_midpoints]#get dist for each vertex
            self.avg_asiGlia2psdCentroid_dist = np.mean(asiGlia2psdCentroid_dist) * self.conversion_length  #get vert avg

            #get avg dist from each asiGlia vertex to asi centroid (use global coordinates)
            asiGlia2asiCentroid_dist = [math.dist(i[1], self.asi_centroid) for i in asiGlia_midpoints] #get dist for each vertex
            self.avg_asiGlia2asiCentroid_dist = np.mean(asiGlia2asiCentroid_dist) * self.conversion_length  #take overall mean for all asiGlia verts

            #get avg dist ratio from each asiGlia vertex to psd centroid/ASI centroid
            asiGlia2asi_psdCentroid_distRatio = (np.array(asiGlia2psdCentroid_dist)/np.array(asiGlia2asiCentroid_dist)).tolist() #dist ratio per vertex
            self.avg_asiGlia2asi_psdCentroid_distRatio = np.mean(asiGlia2asi_psdCentroid_distRatio) #mean dist ratio for all asiGlia verts

        else:
            print('no asiGlia coverage for this synapse')
            self.min_asiGlia2glia_dist = np.nan
            self.avg_asiGlia2glia_dist = np.nan
            self.min_asiGlia2psd_dist = np.nan
            self.max_asiGlia2psd_dist = np.nan
            self.avg_asiGlia2psd_dist = np.nan
            self.min_glia2psd_distRatio = np.nan
            self.max_glia2psd_distRatio = np.nan
            self.avg_glia2psd_distRatio = np.nan
            self.avg_asiGlia2psdCentroid_dist = np.nan
            self.avg_asiGlia2asiCentroid_dist = np.nan
            self.avg_asiGlia2asi_psdCentroid_distRatio = np.nan

    #'THIS FUNCTION ASSUMES ONLY ONE GLIA NAME'
    #have to extract asiGlia_bL first for each appo Dist -> assume only one g_name!!! 
    def extract_glia_vector(self, contour:bool=False, mesh:bool=False):
        #define bmeshes
        ax = bpy.data.objects[self.ax_name]
        psd = bpy.data.objects[self.psd_ax_name]
        ax_bm, ax_bm_faces, ax_bm_edges = ahf.create_bmesh(ax)
        ax_bm.edges.ensure_lookup_table()
        g = [bpy.data.objects[g] for g in self.g_name]

        #for all asiGlia vertices - get distance to closest point on psd (for only asi peri vertices with glia within 120 nm)
        #last list of asiGlia_bL corresponds to vertices that had astro within 120 nm
        asiGlia_midpoints = [ahf.get_edge_midpoint(ax_bm.edges[edge], ax.matrix_world) for edge in self.asiGlia_bL[-1]]
        
        if(len(asiGlia_midpoints) >0):
            print('length of asiGlia modpoints = ', len(asiGlia_midpoints))
            #toward vector (if glia "attacts" psd toward) (target:ASI centroid - origin:glia vector)
            toward_vector = Vector()
            for i in asiGlia_midpoints: 
                #calc vector 
                vec = self.asi_centroid - i[1]
                #convert to unit vector
                unit_vec = vec/(math.dist(i[1], self.asi_centroid)) #divide by distance
                #add to 'toward vector
                toward_vector = toward_vector + unit_vec

            #away vector (if glia "repels" psd away) (target: glia vector - origin: ASI centroid)
            away_vector = Vector()
            for i in asiGlia_midpoints: 
                #calc vector 
                vec = i[1] - self.asi_centroid
                #convert to unit vector
                unit_vec = vec/(math.dist(i[1], self.asi_centroid)) #divide by distance
                #add to 'toward vector
                away_vector = away_vector + unit_vec
            
            #record vectors
            self.toward_vec_x = toward_vector.x
            self.toward_vec_y = toward_vector.y
            self.toward_vec_z = toward_vector.z
            self.away_vec_x = away_vector.x
            self.away_vec_y = away_vector.y
            self.away_vec_z = away_vector.z

        else:
            print('no asiGlia coverage for this synapse')
            self.toward_vec_x = np.nan
            self.toward_vec_y = np.nan
            self.toward_vec_z = np.nan
            self.away_vec_x = np.nan
            self.away_vec_y = np.nan
            self.away_vec_z = np.nan
    
        ax_bm.free()
    #number of asiGlia segments!
    def extract_asiGlia_segments(self):
        #define bmeshes
        ax = bpy.data.objects[self.ax_name]
        ax_bm, ax_bm_faces, ax_bm_edges = ahf.create_bmesh(ax)
        ax_bm.edges.ensure_lookup_table()

        #number of asiGlia continous segments using different appos
        asiGlia_seg_count = []
        asiGlia_seg_avgLength = []
        for i, j in enumerate(self.glia_appo):
            asiGlia_bL_idxs = self.asiGlia_bL[i] #list of asiGlia boundary loop edge idxs for appo dist i
            if len(asiGlia_bL_idxs)>0:
                asiGLia_bL_edges = [ax_bm_edges[edge_idx] for edge_idx in asiGlia_bL_idxs] #list of actual edge references
                asiGlia_bL_segments = ahf.find_connected_segments(asiGLia_bL_edges) #list of connected edge segments

                seg_count = len(asiGlia_bL_segments) #number of segments
                #average length of connected segments
                seg_lengths = [sum(e.calc_length() for e in segment) for segment in asiGlia_bL_segments]
                avg_seg_length = (sum(seg_lengths) / len(seg_lengths))*self.conversion_length
            else:
                #set default vals to 0 if no glia apposed edges
                seg_count = 0
                avg_seg_length = 0
            
            #append to list
            asiGlia_seg_count.append(seg_count)
            asiGlia_seg_avgLength.append(avg_seg_length)
        
        #record values
        self.asiGlia_seg_count = asiGlia_seg_count 
        self.asiGlia_seg_avgLength = asiGlia_seg_avgLength

        ax_bm.free()

    #extract vol and surface area of spines import from pyrecon
    def extract_pyrecon_spine_stats(self):
        pyrecon_sp_name = self.sp_name.split('_remesh')[0]
        pyrecon_sp = bpy.data.objects[pyrecon_sp_name]

        #convert to bm
        pyrecon_sp_bm, pyrecon_sp_bm_faces, pyrecon_sp_bm_edges = ahf.create_bmesh(pyrecon_sp)

        #ensure bmesh data is up-to-date
        pyrecon_sp_bm.normal_update()

        #calculate surface area
        sa = sum(f.calc_area() for f in pyrecon_sp_bm.faces)
        volume = pyrecon_sp_bm.calc_volume()

        self.pyrecon_sp_area = sa * (self.conversion_length **2)
        self.pyrecon_sp_vol = volume * (self.conversion_length **3)

        pyrecon_sp_bm.free()



