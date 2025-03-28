import os, sys, tempfile
import subprocess as sp
import csv
import numpy as np


'''access cgal executable''' 
def cgal_remesh_smooth(bin_dir, temp_dir_path, filter=[''], ex_filter = ['foo'], remesh:bool=False, smooth_seg:bool=False):
    filepaths = os.listdir(temp_dir_path)
    filepaths = [i for i in filepaths if any(f in i for f in filter) and any(ef not in i for ef in ex_filter)]

    #remesh
    if remesh:
        print('run cgal script to remesh')
        for path in filepaths:
            cmd1 = [bin_dir + "cgal_remesh", temp_dir_path + '/' + path] #run cgal executable to smooth -> segment 
            proc = sp.call(cmd1)
            print('testing')
            if proc != 0:
                print('remesh error')
            else:
                print('remesh done')
    
    #smooth and segment w/ cgal
    if smooth_seg:
        print('run cgal script to segment mesh')
        for path in filepaths:
            cmd2 = [bin_dir + "cgal_seg",  temp_dir_path + '/' + path] #run cgal executable to smooth -> segment 
            proc = sp.call(cmd2)
            if proc != 0:
                print('segment error')
            else:
                print('segment done')


