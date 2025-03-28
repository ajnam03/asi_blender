
blender_path="Blender.app/Contents/MacOS" #replace with path to blender
scripts_path_asi="asi_blender/blender_asi_method.py" #replace with path to blender_asi_method.py script

declare -a blend_paths=(
# PATH TO BLENDER FILES
)

echo 'blender analysis'
for index in "${!blend_paths[@]}"; do 
./Blender --background --python $script_path_asi -- ${blend_paths[$index]}
done
