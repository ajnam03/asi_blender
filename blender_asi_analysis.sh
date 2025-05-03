#!/usr/bin/env bash

blender_path="/Applications/Blender.app/Contents/MacOS/Blender"  # replace with the full path to Blender binary
script_path_asi="blender_asi_method.py"  # replace with the correct path to your script

# Declare an array of .blend file paths
declare -a blend_paths=(
    # "/path/to/file1.blend"
    # "/path/to/file2.blend"
)

echo 'Starting Blender analysis...'

for index in "${!blend_paths[@]}"; do 
    "$blender_path" --background --python "$script_path_asi" -- "${blend_paths[$index]}"
done