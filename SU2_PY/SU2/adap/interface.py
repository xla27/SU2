#!/usr/bin/env python

## \file interface.py
#  \brief Wrapper functions for interfacing with the Inria AMG library
#  \author Victorien Menier, Brian Mungu\'ia
#  \version 8.0.1 "Harrier"
#
# SU2 Project Website: https://su2code.github.io
#
# The SU2 Project is maintained by the SU2 Foundation
# (http://su2foundation.org)
#
# Copyright 2012-2022, SU2 Contributors (cf. AUTHORS.md)
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

import sys, os, shutil
import subprocess
import csv
import numpy as np

# ------------------------------------------------------------
#  Setup
# ------------------------------------------------------------
MMG_RUN = os.environ["MMG_RUN"]
sys.path.append(MMG_RUN)
command_mmg2D = MMG_RUN+'/mmg2d_O3'
command_mmg3D = MMG_RUN+'/mmg3d_O3'

# ------------------------------------------------------------
#  SU2-MMG Interface Functions
# ------------------------------------------------------------

class MeshSolConverter():
    """
    Class to convert .su2 to .mesh files and viceversa.
    Class to convert .csv or .dat SS2 solution files to .sol
    Works only with 2D-triangular and 3D-tetrahedral unstructured meshes.
    """

    def __init__(self, verbose=False):
        self.verbose = verbose
        return
    
    def SetDim(self, dim):
        """
        Setting the number of dimensions.
        """
        self.dim = dim

    def GetDim(self):
        """
        Returning the number of dimensions.
        """
        return self.dim
    
    def SetMeshDict(self, mesh_dict):
        """
        Setting the mesh dictionary once read either from .su2 or .mesh
        """
        self.mesh_dict = mesh_dict

    def GetMeshDict(self):
        """
        Returning the mesh dictionary.
        """
        return self.mesh_dict
    
    def SetMetricDict(self, metric_dict):
        """
        Setting the dictionary for adaptation metric values.
        """
        self.metric_dict = metric_dict

    def GetMetricDict(self):
        """
        Returning the dictionary for adaptation metric values.
        """
        return self.metric_dict
        
    def SetSU2MeditMarkersMap(self, su2_markers_list):
        """
        Constructing a unique between SU2 markers and Medit colors.
        """
        self.markers_map = []
        for medit_tag, su2_tag in enumerate(su2_markers_list):
            # medit_tag starts from "1" since the tag "0" is left for the volume domain
            self.markers_map.append([str(medit_tag+1), su2_tag])

    def GetSU2MeditMarkersMap(self):
        """
        Returning the map bewteen SU2 markers and Medit colors.
        """
        if self.markers_map:
            return self.markers_map
        else:
            raise ValueError('Su2-Medit markers map not set!')  

    def GetSU2Marker(self, medit_tag):
        """
        Returning the SU2 markers correspondent to a Medit color.
        """
        for match in self.markers_map:
            if match[0] == medit_tag:
                return match[1]

    def GetMeditMarker(self, su2_tag):
        """
        Returning the Medit color correspondent to a SU2 marker.
        """
        for match in self.markers_map:
            if match[1] == su2_tag:
                return match[0]

    def ReadMeshSU2(self, su2_filename):
        """ 
        Reads a .su2 mesh file and returns node coordinates, elements, and boundary markers in a dictionary data structure. 
        """
        with open(su2_filename, "r") as f:
            lines = f.readlines()

        vertices = []
        elements = []
        boundaries = {}

        su2_markers_list = []
        
        i = 0
        while i < len(lines):
            line = lines[i].strip()

            if line.startswith("%"):
                i += 1
                continue

            if line.startswith("NDIME="):
                dim = int(line.split("=")[1].strip())
                self.SetDim(dim)

            elif line.startswith("NPOIN="):
                n_vertices = int(line.split("=")[1].strip())
                if dim == 2:
                    for j in range(n_vertices):
                        x, y = lines[i + 1 + j].split()[:dim]
                        vertices.append([float(x), float(y)])
                elif dim == 3:
                    for j in range(n_vertices):
                        x, y, z = lines[i + 1 + j].split()[:dim]
                        vertices.append([float(x), float(y), float(z)])
                i += n_vertices  # Move index past nodes

            elif line.startswith("NELEM="):
                n_elements = int(line.split("=")[1].strip())
                for j in range(n_elements):
                    elem_type, *id_vert = lines[i + 1 + j].split()[:(dim+2)]
                    elements.append([int(vert) for vert in id_vert])
                i += n_elements  # Move index past elements

            elif line.startswith("NMARK="):
                n_markers = int(line.split("=")[1].strip())
                for _ in range(n_markers):
                    i += 1
                    marker_tag = lines[i].split("=")[1].strip()
                    su2_markers_list.append(marker_tag)
                    i += 1
                    n_faces = int(lines[i].split("=")[1].strip())
                    boundaries[marker_tag] = []
                    for j in range(n_faces):
                        face_type, *id_vert = lines[i + 1 + j].split()
                        boundaries[marker_tag].append([int(vert) for vert in id_vert])
                    i += n_faces  # Move index past boundary elements

            i += 1

        mesh_dict = {'Dim': dim, 'Vertices': vertices}
        if int(elem_type) == 5:
            mesh_dict['Triangles'] = elements
        if int(elem_type) == 10:
            mesh_dict['Tetrahedra'] = elements

        # reordering markers in alphabetical orders
        #boundaries = {key: value for key, value in sorted(boundaries.items())}
        if int(face_type) == 3:
            mesh_dict['Edges'] =  boundaries
        if int(face_type) == 5:
            mesh_dict['Triangles'] =  boundaries

        self.SetMeshDict(mesh_dict)
        self.SetSU2MeditMarkersMap(su2_markers_list)

        return mesh_dict

    def ReadSolSU2(self, su2_filename):
        """
        Reads a .csv/.dat SU2 solution file to obtain the metric. 
        """
        metric_dict = {}

        if '.dat' in su2_filename:
            fieldnames, data = read_SU2_restart_binary(su2_filename)
        elif '.csv' in su2_filename:
            fieldnames, data = read_SU2_restart_ascii(su2_filename) 

        if "z" in fieldnames:
            dim = 3
        else:
            dim = 2

        metric_dict['Dim'] = dim

        if dim == 2:
            metric_dict ['Metric_xx'] = data[:,-3]
            metric_dict ['Metric_xy'] = data[:,-2]
            metric_dict ['Metric_yy'] = data[:,-1]
        elif dim == 3:
            metric_dict ['Metric_xx'] = data[:,-6]
            metric_dict ['Metric_xy'] = data[:,-5]
            metric_dict ['Metric_yy'] = data[:,-4]
            metric_dict ['Metric_xz'] = data[:,-3]
            metric_dict ['Metric_yz'] = data[:,-2]
            metric_dict ['Metric_zz'] = data[:,-1]      

        metric_dict['NumberVertices'] = data.shape[0]  
        
        self.SetMetricDict(metric_dict)

        return metric_dict
    
    def ReadMeshMedit(self, medit_filename):
        """ 
        Reads a .mesh file and returns node coordinates, elements, and boundary markers  in a dictionary data structure. 
        """
        with open(medit_filename, "r") as f:
            lines = f.readlines()

        dim = None
        vertices = []
        elements = []
        boundaries = {}

        i = 0
        while i < len(lines):
            line = lines[i].strip()

            if line.startswith("Dimension"):
                dim = int(line.split()[1])
                self.SetDim(dim)

            elif line.startswith("Vertices"):
                n_vert = int(lines[i+1])
                if dim == 2:
                    for j in range(n_vert):
                        x, y, id_dom = lines[i + 2 + j].split()
                        vertices.append([float(x), float(y)])
                elif dim == 3:
                    for j in range(n_vert):
                        x, y, z, id_dom = lines[i + 2 + j].split()
                        vertices.append([float(x), float(y), float(z)])
                i += n_vert  # Move index past nodes

            elif (line.startswith("Triangles") and dim == 2) or (line.startswith("Tetrahedra") and dim == 3):
                n_elements = int(lines[i+1])
                for j in range(n_elements):
                    elem_data = list(map(int, lines[i + 2 + j].split()))
                    elem_data = [elem-1 for elem in elem_data]
                    elements.append(elem_data[:-1])  # Last column is a region marker
                i += n_elements  # Move index past elements

            elif (line.startswith("Edges") and dim == 2) or (line.startswith("Triangles") and dim == 3):
                # These define boundary markers
                n_faces = int(lines[i+1])
                for j in range(n_faces):
                    face_data = list(map(int, lines[i + 2 + j].split()))
                    marker = str(face_data[-1])  # Last column is the boundary marker
                    face_data = [face - 1 for face in face_data[:-1]]
                    if marker not in boundaries.keys():
                        boundaries[marker] = []
                    boundaries[marker].append(face_data)  # Store only connectivity
                i += n_faces  # Move index past boundary elements

            i += 1

        mesh_dict = {'Dim': dim, 'Vertices': vertices}
        if dim == 2:
            mesh_dict['Triangles'] = elements
        if dim == 3:
            mesh_dict['Tetrahedra'] = elements

        # reordering markers in alphabetical orders
        if dim == 2:
            mesh_dict['Edges'] =  boundaries
        if dim == 3:
            mesh_dict['Triangles'] =  boundaries

        self.SetMeshDict(mesh_dict)

        return mesh_dict

    def WriteMeshSU2(self, su2_filename):
        """ 
        Writes a .su2 mesh file from given mesh data. 
        """
        mesh = self.GetMeshDict()
        dim = mesh['Dim']
        vertices = mesh["Vertices"]
        if dim == 2:
            elements = mesh['Triangles']
            elem_type = 5
            boundaries = mesh['Edges']
            face_type = 3
        if dim == 3:
            elements = mesh['Tetrahedra']
            elem_type = 10
            boundaries = mesh['Triangles']
            face_type = 5


        with open(su2_filename, "w") as f:
            f.write("NDIME= {}\n".format(dim))

            f.write("NELEM= {}\n".format(len(elements)))
            for elem in elements:
                f.write("{} {}\n".format(elem_type, " ".join(map(str, elem))))

            f.write("NPOIN= {}\n".format(len(vertices)))
            for i, node in enumerate(vertices):
                f.write("{}\n".format(" ".join(map(str, node))))

            f.write("NMARK= {}\n".format(len(boundaries)))
            for medit_tag in boundaries.keys():
                f.write("MARKER_TAG= {}\n".format(self.GetSU2Marker(medit_tag)))
                f.write("MARKER_ELEMS= {}\n".format(len(boundaries[medit_tag])))
                for face in boundaries[medit_tag]:
                    f.write("{} {}\n".format(face_type, " ".join(map(str, face))))

        return
    
    def WriteMeshMedit(self, medit_filename):
        """ 
        Writes a .mesh mesh file from given mesh data. 
        """
        mesh = self.GetMeshDict()
        dim = mesh['Dim']
        vertices = mesh['Vertices']
        if dim == 2:
            elements = mesh['Triangles']
            boundaries = mesh['Edges']
        if dim == 3:
            elements = mesh['Tetrahedra']
            boundaries = mesh['Triangles']

        with open(medit_filename, "w") as f:
            f.write("MeshVersionFormatted 2\n")
            f.write("Dimension {}\n".format(mesh['Dim']))
            
            # Write nodes
            f.write("\nVertices \n{}\n".format(len(mesh['Vertices'])))
            for vert in vertices:
                f.write(" ".join(map(str, vert[:dim])) + " 0\n")  # 0 is the default region ID
            
            # Write elements (assume triangles for 2D, tetrahedra for 3D)
            if dim == 2:
                f.write("\nTriangles \n{}\n".format(len(elements)))
            elif dim == 3:
                f.write("\nTetrahedra \n{}\n".format(len(elements)))
            
            for elem in elements:
                elem = [el + 1 for el in elem]
                f.write(" ".join(map(str, elem)) + " 0\n")  # Last value is a region ID
            
            # Write boundary elements correctly
            if boundaries:
                if dim == 2:
                    f.write("\nEdges \n{}\n".format(sum(len(faces) for faces in boundaries.values())))
                elif dim == 3:
                    f.write("\nTriangles \n{}\n".format(sum(len(faces) for faces in boundaries.values())))
                
                for su2_tag in boundaries.keys():
                    for face in boundaries[su2_tag]:
                        face = [fa + 1 for fa in face]
                        f.write(" ".join(map(str, face)) + " {}\n".format(self.GetMeditMarker(su2_tag)))

            f.write("\nEnd\n")
        
        return

    def WriteSolMedit(self, medit_filename):
        """ 
        Writes a .sol Medit file from given metric data. 
        """
        metric = self.GetMetricDict()

        dim = metric['Dim']
        numvert = metric['NumberVertices']

        with open(medit_filename, "w") as f:
            f.write("MeshVersionFormatted 2\n")
            f.write("Dimension {}\n".format(dim))
            f.write("SolAtVertices\n")
            f.write("{}\n".format(numvert))
            f.write("1 3\n")

            
            # Write metric per node
            if dim == 2:
                for vert in range(numvert):
                    met = [metric['Metric_xx'][vert],
                           metric['Metric_xy'][vert],
                           metric['Metric_yy'][vert],]
                    f.write("{}\n".format(" ".join(map(str, met))))

            if dim == 3:
                for vert in range(numvert):
                    met = [metric['Metric_xx'][vert],
                           metric['Metric_xy'][vert],
                           metric['Metric_yy'][vert],
                           metric['Metric_xz'][vert],
                           metric['Metric_yz'][vert],
                           metric['Metric_zz'][vert],]
                    f.write("{}\n".format(" ".join(map(str, met))))

            f.write("\nEnd\n")
        
        return
    
    def SU2ToMeditMesh(self, su2_filename, medit_filename):
        """
        Full mesh file conversion (reading-writing) from SU2 to Medit
        """
        self.ReadMeshSU2(su2_filename)
        self.WriteMeshMedit(medit_filename)
        if self.verbose:
            print(f"Converted {su2_filename} to {medit_filename}")

    def SU2ToMeditSol(self, su2_filename, medit_filename):
        """
        Full sol file conversion (reading-writing) from SU2 to Medit
        """
        self.ReadSolSU2(su2_filename)
        self.WriteSolMedit(medit_filename)
        if self.verbose:
            print(f"Converted {su2_filename} to {medit_filename}")

    def MeditToSU2Mesh(self, medit_filename, su2_filename):
        """
        Full mesh file conversion (reading-writing) from Medit to SU2
        """
        self.ReadMeshMedit(medit_filename)
        self.WriteMeshSU2(su2_filename)
        if self.verbose:
            print(f"Converted {medit_filename} to {su2_filename}")

    def WriteParamFile(self, config_mmg, mesh_filename):
        """
        Writing the .mmg2d/.mmg3d parameter file if required. 
        """
        param_required = isinstance(config_mmg['hausd'], dict)
        if param_required:
            mesh = self.GetMeshDict()
            dim = mesh['Dim']
            if dim == 2:
                boundaries = mesh['Edges']
                elem_type = 'Edges'
                mmg_ext = '.mmg2d'
            if dim == 3:
                boundaries = mesh['Triangles']
                elem_type = 'Triangles'
                mmg_ext = '.mmg3d'

            param_filename = mesh_filename + mmg_ext
            with open(param_filename, 'w') as f:
                f.write('Parameters\n')
                f.write(str(len(boundaries.keys()))+'\n')
                f.write('\n')

                if len(boundaries.keys()) != len(config_mmg['hausd'].keys()):
                    print('WARNING: Different number of markers between SU2 (%i) mesh and MMG parameters (%i). ' \
                                  'For unspecified markers, HAUSD = 0.01 is assumed.' %
                                   (len(boundaries.keys()), len(config_mmg['hausd'].keys())))
                
                for su2_tag in config_mmg['hausd'].keys():
                    f.write('%s %s %1.2e %1.2e %1.2e\n' % 
                            (self.GetMeditMarker(su2_tag), 
                            elem_type, 
                            config_mmg['hmin'], 
                            config_mmg['hmax'], 
                            config_mmg['hausd'][su2_tag]))       
        
        else:
            pass

        return      

def call_mmg(meshin, meshout, solfile, config_mmg):
    """Adapt mesh using pyamg module"""

    remesh_options = config_mmg
    remesh_options['meshfile'] = meshin
    remesh_options['meshoutfile'] = meshout
    remesh_options['solfile'] = solfile

    dim = config_mmg['dim']

    if dim == 2:
        the_Command = build_command(command_mmg2D, remesh_options)
    elif dim == 3:
        the_Command = build_command(command_mmg3D, remesh_options)
    else:
        raise KeyError('Wrong number of dimensions!')
    
    try:
        run_command(the_Command)
    except:
        raise RuntimeError("mmg failed.")

    return 

# ------------------------------------------------------------
#  Helper functions
# ------------------------------------------------------------

def build_command(command_mmg, options):
    """builds the mmg command with options"""

    the_Command = command_mmg 
    the_Command += ' -in '    + options['meshfile'] 
    the_Command += ' -sol '   + options['solfile'] 
    the_Command += ' -out '   + options['meshoutfile']
    the_Command += ' -hmin '  + str(options['hmin'])
    the_Command += ' -hmax '  + str(options['hmax'])
    if 'hgrad' in options.keys():
        the_Command += ' -hgrad ' + str(options['hgrad'])
    if 'hausd' in options.keys() and not isinstance(options['hausd'], dict):
        the_Command += ' -hausd ' + str(options['hausd'])
    the_Command += ' > ' + options['mmg_log'] + ' 2> ' + options['mmg_err']
    
    return the_Command


def run_command(Command):
    """runs os command with subprocess
    checks for errors from command
    """

    sys.stdout.flush()

    proc = subprocess.Popen(
        Command, shell=True, stdout=sys.stdout, stderr=subprocess.PIPE
    )
    return_code = proc.wait()
    message = proc.stderr.read().decode()

    return 


CGNS_STRING_SIZE = 33  # Fixed string size per CGNS standard

def read_SU2_restart_binary(filename):
    """
    Read SU2 binary restart file and return fields and data array.

    Returns:
        fields (List[str]): Field names including "Point_ID".
        data (np.ndarray): Data array of shape (nPoints, nFields-1).

    Note that the Point_ID column is implicit in the ordering
    """
    #filename += ".dat"
    fields = ["Point_ID"]  # Initialize with Point_ID as SU2 convention

    with open(filename, 'rb') as f:
        # Read 5 integers (magic number + metadata)
        header = np.fromfile(f, dtype=np.int32, count=5)
        if header.size != 5:
            raise RuntimeError("Error reading header from restart file.")
        
        magic_number, nFields, nPoints, _, _ = header

        # Check the magic number
        if magic_number != 535532:
            raise RuntimeError(f"{filename} is not a binary SU2 restart file.")

        # Read field names (each is CGNS_STRING_SIZE characters)
        for _ in range(nFields):
            name_bytes = f.read(CGNS_STRING_SIZE)
            name_str = name_bytes.decode('utf-8').strip('\x00').strip()
            fields.append(name_str)

        # Read restart data as a flat array of doubles
        data = np.fromfile(f, dtype=np.float64, count=nFields * nPoints)

        if data.size != nFields * nPoints:
            raise RuntimeError("Error reading restart data.")

        # Reshape to 2D: each row is a point, each column is a field
        data = data.reshape((nPoints, nFields))

    return fields, data


def read_SU2_restart_ascii(filename):
    """
    Read SU2 ASCII restart file and return fields and data array.

    Returns:
        fields (List[str]): Field names.
        data (np.ndarray): Data array of shape (nPoints, nFields).
    """

    # reading the first line to get the fields name
    try:
        with open(filename, "r") as f:
            line = f.readline()
    except:
        raise("The solution file must be in ASCII format!")

    fields = line.lstrip('"').rstrip('"\n')
    fields = fields.split('","')

    data = np.empty((0, len(fields)))
    with open(filename, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            data = np.vstack((data, np.array(line.split(','), dtype=np.float64)))

    return fields, data
