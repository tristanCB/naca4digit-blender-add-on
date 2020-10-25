# Python 3
# An addon for Blender 2.90 which generates points for the outline
# of a 4 digit series NACA airfoil
# Copyright 2020 Tristan C.B.
# License: MIT

bl_info = {
    "name": "Add 4-digit NACA profile outline",
    "blender": (2, 90, 0),
    "category": "Object",
}

import bpy
import bmesh
from bpy_extras.object_utils import AddObjectHelper
import numpy as np
import math

def fourDigitSeries(digits:int, numPanels: int, resType = "LINEAR"):
    """
    Simple functions which returnsNACA profiles as needed by aifoil.py

    returns  X coords, Y coords -> Both numpy arrays of shape (numPanels+1)
    The format is specifically for the vortex panel method.
    return X coords and Y coords describing a NACA foil. Equations are taken from Moran (2003) 

    http://web.stanford.edu/~cantwell/AA200_Course_Material/The%20NACA%20airfoil%20series.pdf
        > " The first digit specifies the maximum camber (m) in percentage of the chord (airfoil length),
        the second indicates the position of the maximum camber (p) in tenths of chord, and the last two
        numbers provide the maximum thickness (t) of the airfoil in percentage of chord. For example, the
        NACA 2415 airfoil has a maximum thickness of 15% with a camber of 2% located 40% back from
        the airfoil leading edge (or 0.4c)"

    references: 
        Moran, Jack. An Introduction to Theoretical and Computational Aerodynamics. Dover Publications, 2003.
        M. and Chowc-Y. Foundations of Aerodynamics: Bases of Aerodynamic Design — Fourth Edition. John Wieley & Sons, 1968.
    # """
    
    # Should handle error division by zeros
    # Add camber line
    # Add resolution type
    
    # Panels should be a multiple a 2.
    # assert numPanels % 2 == 0

    digits = str(digits).zfill(4)

    ## Definitions ##
    # maxCamber = location of the max camber
    # cpMaxCamber = chordwise position of the mmaximum camber
    # thicknessRatio =  airfoil having a thickness XX % (base chord) of the chord.
    # c = chord length = 1

    # Extracts digits describing specific physical params.
    maxCamber, cpMaxCamber, thicknessRatio = [int(i) for i in [digits[0], digits[1:2], digits[2:4]]]
    # Convert to format needed for computaion
    epsilon = maxCamber/100
    p = cpMaxCamber/10
    t = thicknessRatio/100
    c = 1

    length = int(numPanels/2+1)
    
    if resType == "COSINE":
    # Cosine spacing for discretization scheme applied to half circle b/c we are mirroring it, we want y_upper and y_lowwer for same x
        x = np.asarray(([0.5*(np.cos(i)+1) for i in np.linspace(0,np.pi,num=length)]))
    elif resType == "QUADRATIC":
    # Spacing increasing resolution near begining of profile.
        x = np.linspace(1,0,int(numPanels/2+1))**2 # power of 2 will increase resoltion near the leading edge
    elif resType == "LINEAR":
        x = np.linspace(1,0,int(numPanels/2+1))
    # elif resType == "LOGSPACE":
    #     x = np.logspace(0,10,num=int(numPanels/2+1))[::-1]

    meanCamberLine = np.zeros((x.shape))
    yt, dycOVERdx, x_upper, x_lower, y_upper, y_lower, theta  = [np.zeros((x.shape)) for i in range(7)]

    # Main loop to derive coords from equation.
    for i, ij in enumerate(x):
        if i >= 0 and ij <= p:
            meanCamberLine[i] =  epsilon/p**2 * (2*p*ij-ij**2)
            dycOVERdx[i] = 2 * epsilon/p**2 * (p - ij/c)
        else:
            meanCamberLine[i] = epsilon/(1-p)**2 * ((1-2*p)+2*p*ij-ij**2)
            dycOVERdx[i] = 2 * epsilon / (1-p)**2 * (p - ij/c)

        # yt[i] = t/0.2*(0.2969*np.sqrt(ij) - 0.1260*ij - 0.3516*ij**2 + 0.2843*ij**3 - 0.1015*ij**4)
        
        # Trailing edge adjustment to close it.
        yt[i] = t/0.2*(0.2969*np.sqrt(ij) - 0.1260*ij - 0.3516*ij**2 + 0.2843*ij**3 - 0.1036*ij**4)
        theta[i] = np.arctan(dycOVERdx[i])
         
        x_upper[i], y_upper[i] = (ij - yt[i]*np.sin(theta[i])), (meanCamberLine[i] + yt[i]*np.cos(theta[i]))
        x_lower[i], y_lower[i] = (ij + yt[i]*np.sin(theta[i])), (meanCamberLine[i] - yt[i]*np.cos(theta[i]))
    
    # Stitch the arrays together to get convention clockwise from leading edge
    # Takes the x_lower array in inverse order and removes coinciding point (0,0)
    XB = np.hstack((x_lower,x_upper[::-1][1:]))
    YB = np.hstack((y_lower,y_upper[::-1][1:]))
    
    # To run the airfoil simulation we need the same format at prescribed bellow.
    # The bellow array of for a NACA 2412, with slight modifications to make the trailing edge close
    YB[0], YB[-1] = 0, 0
    # YB[int(YB.shape[0]/2)] = 0

    # print(
    # f'''
    # ## NACA foil 4 digit series: {digits} ##
    # Physical params:
    # epsilon => {epsilon} 
    # p       => {p}
    # t       => {t}
    # ########################################
    # --- XB{XB.shape} ---
    #         {XB}
    # --- YB{YB.shape} ---
    #         {YB}
    # '''
    # )

    return XB, YB

def generate_NACA_mesh(width, camber, cpMaxCamber, thicknessRatio, resolution, resolutionType):
    """
    This function takes inputs and returns vertex and face arrays.
    no actual mesh data creation is done here.
    """
    if thicknessRatio <= 9:
        thicknessRatio = "0"+str(thicknessRatio)
        
    digitIdentifier = int(f"{camber}{cpMaxCamber}{thicknessRatio}")
    print(f"Plotting NACA{digitIdentifier}")
    
    
    XB, YB = fourDigitSeries(digitIdentifier,resolution, resType = resolutionType)
    # print(XB)
    
    verts = [(XB[i], YB[i], 0) for i, _ in enumerate(XB)]
    # faces = [([i for i, _ in enumerate(XB)])]
    # verts = [A
    #     (+1.0, +1.0, -1.0),
    #     (+1.0, -1.0, -1.0),
    #     (-1.0, -1.0, -1.0),
    #     (-1.0, +1.0, -1.0),
    #     (+1.0, +1.0, +1.0),
    #     (+1.0, -1.0, +1.0),
    #     (-1.0, -1.0, +1.0),
    #     (-1.0, +1.0, +1.0),
    # ]

    # faces = [
    #     (0, 1, 2, 3),
    #     (4, 7, 6, 5),
    #     (0, 4, 5, 1),
    #     (1, 5, 6, 2),
    #     (2, 6, 7, 3),
    #     (4, 0, 3, 7),
    # ]

    faces = [
    ]

    # apply size
    for i, v in enumerate(verts):
        verts[i] = v[0] * width, v[1], v[2]

    return verts, faces

from bpy.props import (
    BoolProperty,
    BoolVectorProperty,
    EnumProperty,
    FloatProperty,
    IntProperty,
    FloatVectorProperty,
)

class AddNACA(bpy.types.Operator):
    """Add a simple NACA 4 digit mesh"""
    bl_idname = "mesh.primitive_naca_add"
    bl_label = "Add NACA foil"
    bl_options = {'REGISTER', 'UNDO'}

    width: FloatProperty(
        name="Width",
        description="Box Width",
        min=0.01, max=100.0,
        default=1.0,
    )
    
    # First digit
    camber: IntProperty(
        name="Camber",
        description="Maximum camber as percentage of the chord",
        min=0, max=9,
        default=2,
    )
    # Second digit
    cpMaxCamber: IntProperty(
        name="cpMaxCamber",
        description="Distance of maximum camber from the airfoil leading edge in tenths of the chord",
        min=0, max=9,
        default=4,
    )
    # Last two digits
    thicknessRatio: IntProperty(
        name="thicknessRatio",
        description="Last two digits describing maximum thickness of the airfoil as percent of the chord",
        min=0, max=99,
        default=12,
    )
    
    
    resolution: IntProperty(
        name="Resolution",
        description="Amount of point plotted",
        min=8, max=9999,
        default=64,
    )    

    layers: BoolVectorProperty(
        name="Layers",
        description="Object Layers",
        size=20,
        options={'HIDDEN', 'SKIP_SAVE'},
    )

    # generic transform props
    align_items = (
        ('WORLD', "World", "Align the new object to the world"),
        ('VIEW', "View", "Align the new object to the view"),
        ('CURSOR', "3D Cursor", "Use the 3D cursor orientation for the new object")
    )
    align: EnumProperty(
        name="Align",
        items=align_items,
        default='WORLD',
        update=AddObjectHelper.align_update_callback,
    )
    location: FloatVectorProperty(
        name="Location",
        subtype='TRANSLATION',
    )
    rotation: FloatVectorProperty(
        name="Rotation",
        subtype='EULER',
    )
    resolutionType: EnumProperty(
        name="Resolution type",
        items=(
            ('COSINE', "Cosine", "Cosine spacing"),
            ('LINEAR', "Linear", "Numpy linespace"),
            ('QUADRATIC', "Quadratic", "Numpy linespace**2"),
            #('LOGSPACE', "Logspace", "Numpy logspace"),
                ),
        default='COSINE',
    )

    def execute(self, context):

        verts_loc, faces = generate_NACA_mesh(
            self.width,
            self.camber,
            self.cpMaxCamber,
            self.thicknessRatio,
            self.resolution,
            self.resolutionType
        )

        mesh = bpy.data.meshes.new("airfoil")

        bm = bmesh.new()

        for v_co in verts_loc:
            bm.verts.new(v_co)

        bm.verts.ensure_lookup_table()
        for f_idx in faces:
            bm.faces.new([bm.verts[i] for i in f_idx])

        bm.to_mesh(mesh)
        mesh.update()

        # add the mesh as an object into the scene with this utility module
        from bpy_extras import object_utils
        object_utils.object_data_add(context, mesh, operator=self)

        return {'FINISHED'}


def menu_func(self, context):
    self.layout.operator(AddNACA.bl_idname, icon='MESH_CUBE')

def register():
    bpy.utils.register_class(AddNACA)
    bpy.types.VIEW3D_MT_mesh_add.append(menu_func)

def unregister():
    bpy.utils.unregister_class(AddNACA)
    bpy.types.VIEW3D_MT_mesh_add.remove(menu_func)

if __name__ == "__main__":
    register()

    # test call
    bpy.ops.mesh.primitive_naca_add()
