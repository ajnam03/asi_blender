### Scripts for Automated 3D Analysis of Astroglial Apposition at the Axon-Spine Interface  

#### 🔧 Notes
1. **Blender version**: This script currently works with **Blender 3.6.5**.  
2. **LoopTools Add-on**: Ensure the LoopTools add-on is activated in Blender.  
   ➤ [LoopTools Documentation](https://docs.blender.org/manual/en/3.5/addons/mesh/looptools.html)  
3. **Naming conventions**:
   - `d##c##`: **c-object**
   - `d##sp##`: **spine object**
   - `d##ax##`: **axon object**
   - `astroAll`: **astroglia object**
   - `1um_rad_ref`: **1 µm reference circle**

#### ▶️ How to Run
1. Edit the path to your `.blend` file in `blender_asi_analysis.sh`.
2. Set the path to the Blender executable.
3. Make the script executable:
   ```bash
   chmod u+x blender_asi_analysis.sh