import bpy
import sys
import csv
import numpy as np

import blender_asi_functions as af
import blender_asi_helper_functions as ahf

if __name__== '__main__':
    #filepaths
    blend_args =4
    blend_path = sys.argv[blend_args+1]
    series_name =blend_path.split('/')[-1].split('.')[0]
    outpath = blend_path.split(blend_path.split('/')[-1])[0] + series_name + '_blender_output.csv'

    #open blend file
    bpy.ops.wm.open_mainfile(filepath=blend_path)

    headless = True
    
    collection_name = 'Collection' #define name of collection
    syn_list = [c.name for c in bpy.data.collections[collection_name].all_objects 
                if ('c' in c.name.split('_')[0].rstrip('abc'))] #'_' condition in case there are duplicate axons with 'copy' in name

    #define glia apposition distance (dag) that will be tested: 0.01 nm -> 0.12 nm
    glia_appo = list(np.arange(0.01, 0.13, 0.01))
    glia_appo = [round(i,2) for i in glia_appo]
    
    #define dataframe headers 
    asi_headers = ['syn_name', 
                   'ref_rad', 'ref_area',  #reference circle dimensions to calculate conversion factor
                   'syn_center', 'syn_dist',  #synapse
                   'sp_area', 'sp_vol', #spine area and vol
                   'asi_faces', 'asi_edges', 'asi_peri', 'asi_area', #asi size
                   'asi_centroid_x', 'asi_centroid_y', 'asi_centroid_z', #asi centroid xyz
                   'psd_sp_centroid_x', 'psd_sp_centroid_y', 'psd_sp_centroid_z', #psd centroid xyz (shrinkwrapped to spine)
                   'psd_ax_centroid_x', 'psd_ax_centroid_y', 'psd_ax_centroid_z', #psd centroid xyz (shrinkwrapped to axon)
                   'psd_offset,' #psd offset
                   'min_asi2psd_dist', 'max_asi2psd_dist', 'avg_asi2psd_dist', #dist from ASI to psd 
                   'min_asiGlia2psd_dist', 'max_asiGlia2psd_dist', 'avg_asiGlia2psd_dist', #dist frpm astro-apposed ASI to psd 
                   'min_asiGlia2glia_dist', 'avg_asiGlia2glia_dist'] #dist from astro-apposed ASI to nearest astro
    
    for i in glia_appo:
        asi_headers += ['asiGlia_edges_' + str(i)] #asiGlia boundary edges id headers
    for i in glia_appo:
        asi_headers += ['asiGlia_peri_' + str(i)] #asiGlia perimeter headers
    for i in glia_appo:
        asi_headers += ['asiGlia_frac_' + str(i)] #asiGlia fractions headers 

    asi_results = [] #analysis results
    omitted_syn = [] ## ADD LINE HERE TO RECORD OMITTED SYNAPSES

    '''Make sure you have loop tools addon activated'''

    for syn in syn_list:
        results = []
        print('syn = ', syn)        
        asi = af.asi(syn, glia_appo)
        try:
            asi.check_objects_exist() #check all three tri-syn components exist
        except KeyError:
            print('Missing tri-syn component')
            omitted_syn.append(syn) #append omitted syn
            continue
        asi.calc_unit_conversion()
        asi.duplicateAxon()
        print('extract synapse center')
        asi.extract_syn_center()
        print('get spine surface area and vol')
        asi.extract_spine_stats()
        print('extract asi indices')
        asi.extract_asi_idx()
        print('asi initial len = ', len(asi.asi_idx))
        if len(asi.asi_idx) > 0:
            print('select asi faces')
            asi.select_asi_faces()
            print('clean up ASI')
            asi.cleanUpASI()
            print('asi post cleanup len = ', len(asi.asi_idx))
            print('select update ASI faces')
            asi.select_asi_faces()
            print('convert asi to boundary loop')
            asi.asi_boundaryLoop()
            print('calculate asi perimeter')
            asi.calc_peri_frac(asi = True)
            print('asi_peri = ', asi.asi_peri)
            print('extract asiGlia_boundaryloop')
            asi.asi_glia_boundaryLoop(mesh=True)
            print('calculate asiGlia peri')
            asi.calc_peri_frac(asi_glia = True)
            for i, j in enumerate(glia_appo): 
                print('asiGlia peri ' + str(j) + ' :', asi.asiGlia_peri[i])
                print('asiGlia frac ' + str(j) + ' :', asi.asiGlia_frac[i])
            #in non-headless form blender calc_length() function doesn't seem to adjust units correctly
            if not headless: 
                print('Non-headless Blender form: update of perimeter length units needed')
                asi.asi_peri = asi.asi_peri() * 10**(-3)
                print('asi_peri updated = ', asi.asi_peri)
                for i, j in enumerate(glia_appo):
                    asi.asiGlia_peri[i] = asi.asiGlia_peri[i] * 10**(-3)
                    print('asiGlia peri updated' + str(j) + ' :', asi.asiGlia_peri[i])

            if len(asi.asi_bL_idx) > 0:
                #get min asiGlia to glia dist
                print('asiGlia to glia distance for appo dist = 120 nm')
                asi.extract_asiGlia2Glia_dist(mesh=False)
                # get asi area + centroid
                print('get asi area + geometric centroid')
                asi.extact_asi_area_centroid()
                #shrinkwrap psd obj
                print('shrinkwrap psd to spine and ax')
                asi.shrinkwrap_psd()
                # get psd area + centroid
                print('get psd area + centroid')
                asi.extract_psd_area_centroid()
                #get asi2psd centroid offset
                print('get asi to psd offset')
                asi.extract_psd_centroid_offset()
                #get asi to psd and asiGlia to psd distances
                print('get asi to psd and asiGlia to psd distances')
                asi.extract_asiGlia2psd_dist(mesh=True)
            else:
                print('problem calculating asi for synapse: ', syn)
                omitted_syn.append(syn) #append omitted syn
                continue
                
        else:
            print('problem calculating asi for synapse: ', syn)
            omitted_syn.append(syn) #append omitted syn
            continue

        #append results
        results.append(syn) #append synapse name
        results.append(asi.ref_rad) #append reference circle radius (as a reference to compare against PyReconstruct)
        results.append(asi.ref_area) #append reference circle area (as a reference to compare against PyReconstruct)
        results.append(asi.c_center_global) #append c obj center coordinates
        results.append(asi.syn_dist_convert) #append synapsse distance (i.e. syn cleft size)
        results.append(asi.sp_area) #append spine surface area 
        results.append(asi.sp_vol) #append spine volume
        results.append(asi.asi_idx) #append asi face ids
        results.append(asi.asi_bL_idx) #append asi boundary loop edge ids
        results.append(asi.asi_peri) #append asi perimeter
        results.append(asi.asi_area) #append asi area
        results.append(asi.asi_centroid_x) #geometric centroid (xyz) of asi 
        results.append(asi.asi_centroid_y)
        results.append(asi.asi_centroid_z)
        results.append(asi.psd_sp_centroid_x) #geometric centroid (xyz) of psd shrinkwrapped to spine
        results.append(asi.psd_sp_centroid_y)
        results.append(asi.psd_sp_centroid_z)
        results.append(asi.psd_ax_centroid_x) #geometric centroid (xyz) of psd shrinkwrapped to axon
        results.append(asi.psd_ax_centroid_y)
        results.append(asi.psd_ax_centroid_z)
        results.append(asi.psd_offset) #offset of psd centroid from asi centroid
        results.append(asi.min_asi2psd_dist) #min/max/avg distances between all asi peri verts + closest psd vert 
        results.append(asi.max_asi2psd_dist)
        results.append(asi.avg_asi2psd_dist)
        results.append(asi.min_asiGlia2psd_dist) #min/max/avg distances between just asiGlia peri verts + closest psd vert 
        results.append(asi.max_asiGlia2psd_dist)
        results.append(asi.avg_asiGlia2psd_dist)
        results.append(asi.min_asiGlia2glia_dist) #min/avg distances between asiGlia peri verts + closest glia vert
        results.append(asi.avg_asiGlia2glia_dist)

        for i, j in enumerate(glia_appo):
            results.append(asi.asiGlia_bL[i]) #append asiGlia boundary loop edges idxs based on each appo dist
        for i, j in enumerate(glia_appo):
            results.append(asi.asiGlia_peri[i]) #append asiGlia perimeter based on each appo dist
        for i, j in enumerate(glia_appo):
            results.append(asi.asiGlia_frac[i]) #append asiGlia frac based on each appo dist
        #append to larger results list
        asi_results.append(results)

    #append 
    for o in omitted_syn:
        o_results = [o]
        o_results += [np.nan for h in range(len(asi_headers)-1)]
        asi_results.append(o_results)

    #save blender output as csv file
    with open(outpath, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(asi_headers)
        writer.writerows(asi_results)
    #save
    bpy.ops.wm.save_mainfile(filepath=blend_path)
