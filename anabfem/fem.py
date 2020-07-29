# -----------------------------------------------------------------------------
#     Anabfem
#     Finite element analysis of annular ablation experiments
#
#     Copyright (C) 2020 Alejandro Torres-Sanchez
#
#     This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#     the Free Software Foundation, either version 3 of the License, or
#     (at your option) any later version.
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>.
# -------------------------------------------------------------------------------

import pkg_resources
import os
import numpy as np
import vtk
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import factorized
from vtk.util import numpy_support
from matplotlib import pyplot as plt

class FEM2DActiveElastic:

    """
    This class computes a finite element approximation for a 2D active elastic material attached to a substrate
    with elastic bonds and with a free boundary.

    The class requires vtk for input and output of the mesh and displacement information.

    It depends on parameters [k, zeta_x, zeta_y, bar{K}] / K (parameters are re-normalised by shear modulus)

    This class also provides a method (compute_point_displacement) for the interpolation of the finite element solution
    in a cloud of points, which can also return gradient of point displacements with respect to material parameters.
    It can also return stretching and shear fields at those points (as well as their gradients wrt material parameters)
    """

    def __init__(self, filemesh='circle_0.vtk', internal_mesh=True, lintrans=np.array([[1.0,0.0],[0.0,1.0]])):

        self.Dbf = np.array([[-1.0, -1.0],
                             [1.0, 0.0],
                             [0.0, 1.0]])  # Derivatives of the basis functions in local coordinates of the triangle

        self.parameters = np.ones(4)  # 4 parameters (assuming K fixed), k, \zeta_x, \zeta_y, \bat{K}
        self.filemesh = filemesh

        if internal_mesh:
            self.filemesh = pkg_resources.resource_filename('anabfem', 'meshes/' + filemesh)
            if not os.path.isfile(self.filemesh):
                print("The file " + filemesh + " does not exists in the set of meshes stored by anabfem")
                print("Here is a list of available meshes:")
                filenames = np.array(pkg_resources.resource_listdir('anabfem','meshes/'))
                for f in filenames:
                    print(f)
                print()
                return

        # Load mesh
        self.x_nodes = np.array([])
        self.connec = np.array([])
        self.vtkmesh = None
        self.load_vtk()

        #Transform nodes
        self.x_nodes[:,0:2] = self.x_nodes[:,0:2].dot(lintrans.T)
        self.vtkmesh.GetPoints().SetData(numpy_support.numpy_to_vtk(self.x_nodes))

        # Compute terms in the linear system
        self.mat_K = None
        self.mat_bK = None
        self.mat_k = None
        self.rhs_zx = np.array([])
        self.rhs_zy = np.array([])
        self.linear_system()

        # Set displacements to zero
        self.disp = np.zeros([self.x_nodes.shape[0],2])

        # Initialise empty stretch and shear
        self.stretch = np.array([])
        self.shear = np.array([])

    def update_mesh_displacements(self, compute_shear_and_stretch=False, compute_parameter_gradients=False):

        """
        Update mesh displacements given model parameters
        """

        mat1 = self.mat_K
        mat2 = self.mat_bK
        mat3 = self.mat_k
        rhs1 = self.rhs_zx
        rhs2 = self.rhs_zy

        # [1] SOLVE THE PROBLEM USING FEM

        # Compute total matrix
        mat = mat1 + self.parameters[3] * mat2 + self.parameters[0] * mat3

        # Factorise matrix
        solve = factorized(mat)

        # Solve
        self.disp = solve(self.parameters[1] * rhs1 + self.parameters[2] * rhs2)

        # In case return_paramgrad, compute gradients of the solution wrt parameters
        if compute_parameter_gradients:
            self.dbK_disp = solve(-self.mat_bK * self.disp)
            self.dk_disp = solve(-self.mat_k * self.disp)
            self.dsx_disp = solve(self.rhs_zx)
            self.dsy_disp = solve(self.rhs_zy)

        # Reshape
        self.disp = self.disp.reshape([len(self.x_nodes), 2])

        if compute_parameter_gradients:
            self.dbK_disp = self.dbK_disp.reshape([len(self.x_nodes), 2])
            self.dk_disp = self.dk_disp.reshape([len(self.x_nodes), 2])
            self.dsx_disp = self.dsx_disp.reshape([len(self.x_nodes), 2])
            self.dsy_disp = self.dsy_disp.reshape([len(self.x_nodes), 2])

        if compute_shear_and_stretch:

            wpoints = np.zeros(self.disp.shape[0])
            self.stretch = np.zeros(self.disp.shape[0])
            self.shear = np.zeros([self.disp.shape[0], 2, 2])

            for e in self.connec:
                x_nbors = self.x_nodes[e, 0:2]
                T = x_nbors.T.dot(self.Dbf)
                det = T[0, 0] * T[1, 1] - T[1, 0] * T[0, 1]
                iT = np.array([[T[1, 1], -T[0, 1]], [-T[1, 0], T[0, 0]]]) / det
                Dbf_glo = self.Dbf.dot(iT)

                Du = np.einsum('Ij,Ii->ij', Dbf_glo, self.disp[e])

                self.stretch[e] += np.trace(Du)
                self.shear[e] += 0.5 * (Du + Du.T - np.trace(Du) * np.identity(2))

                wpoints[e] += 1

            self.stretch /= wpoints
            self.shear[:, 0, 0] /= wpoints
            self.shear[:, 0, 1] /= wpoints
            self.shear[:, 1, 0] /= wpoints
            self.shear[:, 1, 1] /= wpoints

    def compute_point_displacements(self, points, compute_shear_and_stretch=False, compute_parameter_gradients=False):

        """
        Compute point displacements (requires update_mesh_displacements to be called beforehand)
        """

        Dbf = np.array([[-1.0, -1.0], [1.0, 0.0],
                        [0.0, 1.0]])  # Derivatives of the basis functions in local coordinates of the triangle

        x_nodes = self.x_nodes
        connec = self.connec
        vtkmesh = self.vtkmesh

        nPts = points.shape[0]  # Number of points
        displ = np.zeros([nPts, 2])  # Displacements
        dbK_displ = np.zeros([nPts, 2])  # Derivative of displacements wrt \bar{K}
        dk_displ = np.zeros([nPts, 2])   # Derivative of displacements wrt k
        dsx_displ = np.zeros([nPts, 2])  # Derivative of displacements wrt \zeta_x
        dsy_displ = np.zeros([nPts, 2])  # Derivative of displacements wrt \zeta_y

        # If one needs to interpolate shear and stretch as well, then initialise equivalent variables
        if compute_shear_and_stretch:

            stretch = np.zeros(nPts)
            shear = np.zeros([nPts, 2, 2])

            if compute_parameter_gradients:
                dbK_stretch = np.zeros(nPts)
                dk_stretch = np.zeros(nPts)
                dsx_stretch = np.zeros(nPts)
                dsy_stretch = np.zeros(nPts)

                dbK_shear = np.zeros([nPts, 2, 2])
                dk_shear = np.zeros([nPts, 2, 2])
                dsx_shear = np.zeros([nPts, 2, 2])
                dsy_shear = np.zeros([nPts, 2, 2])

        # Use vtk to build a cell locator for interpolation of the data
        my_cell_locator = vtk.vtkCellLocator()
        my_cell_locator.SetDataSet(vtkmesh)  # reverse.GetOutput() --> vtkPolyData
        my_cell_locator.BuildLocator()

        # These are some auxiliary variables used in the locator
        cellId = vtk.reference(0)
        c = [0.0, 0.0, 0.0]
        subId = vtk.reference(0)
        d = vtk.reference(0.0)
        pcoords = [0.0, 0.0, 0.0]
        weights = np.array([0.0, 0.0, 0.0])

        # Loop and fill data
        n = 0
        for p in points:

            # Use vtk to find the element and the weights of a given point
            my_cell_locator.FindClosestPoint([p[0], p[1], 0.0], c, cellId, subId, d)
            vtkmesh.GetCell(cellId).EvaluatePosition([p[0], p[1], 0.0], c, subId, pcoords, d, weights)

            # With this we can directly compute the displacement
            displ[n] = weights.dot(self.disp[connec[cellId]])

            # If we need to compute parameter gradients, then use a similar approach to calculate them
            if compute_parameter_gradients:
                dbK_displ[n] = weights.dot(self.dbK_disp[connec[cellId]])
                dk_displ[n] = weights.dot(self.dk_disp[connec[cellId]])
                dsx_displ[n] = weights.dot(self.dsx_disp[connec[cellId]])
                dsy_displ[n] = weights.dot(self.dsy_disp[connec[cellId]])

            # For shear and stretchs, we calculate them in the element where the point lies
            if compute_shear_and_stretch:

                # Positions (x,y) of the nodes in the element
                x_nbors = x_nodes[connec[cellId], 0:2]

                # Transformation matrix T (converts barycentric coordinates (s1,s2) to (x,y))
                T = x_nbors.T.dot(Dbf)
                det = T[0, 0] * T[1, 1] - T[1, 0] * T[0, 1]
                iT = np.array([[T[1, 1], -T[0, 1]], [-T[1, 0], T[0, 0]]]) / det

                # Derivatives of the basis functions in (x,y)
                Dbf_glo = Dbf.dot(iT)

                # Gradient of the displacement field
                Du = np.einsum('Ij,Ii->ij', Dbf_glo, self.disp[connec[cellId]])

                if compute_parameter_gradients:
                    dbK_Du = np.einsum('Ij,Ii->ij', Dbf_glo, self.dbK_disp[connec[cellId]])
                    dk_Du = np.einsum('Ij,Ii->ij', Dbf_glo, self.dk_disp[connec[cellId]])
                    dsx_Du = np.einsum('Ij,Ii->ij', Dbf_glo, self.dsx_disp[connec[cellId]])
                    dsy_Du = np.einsum('Ij,Ii->ij', Dbf_glo, self.dsy_disp[connec[cellId]])

                stretch[n] = np.trace(Du)

                if compute_parameter_gradients:
                    dbK_stretch[n] = np.trace(dbK_Du)
                    dk_stretch[n] = np.trace(dk_Du)
                    dsx_stretch[n] = np.trace(dsx_Du)
                    dsy_stretch[n] = np.trace(dsy_Du)

                shear[n] = 0.5 * (Du + Du.T - stretch[n] * np.identity(2))

                if compute_parameter_gradients:
                    dbK_shear[n] = 0.5 * (dbK_Du + dbK_Du.T - dbK_stretch[n] * np.identity(2))
                    dk_shear[n] = 0.5 * (dk_Du + dk_Du.T - dk_stretch[n] * np.identity(2))
                    dsx_shear[n] = 0.5 * (dsx_Du + dsx_Du.T - dsx_stretch[n] * np.identity(2))
                    dsy_shear[n] = 0.5 * (dsy_Du + dsy_Du.T - dsy_stretch[n] * np.identity(2))

            n += 1

        if compute_shear_and_stretch:
            if compute_parameter_gradients:
                return displ, stretch, shear, dbK_displ, dk_displ, dsx_displ, dsy_displ, dbK_stretch, dk_stretch, \
                       dsx_stretch, dsy_stretch, dbK_shear, dk_shear, dsx_shear, dsy_shear
            else:
                return displ, stretch, shear
        else:
            return displ

    def linear_system(self):
        """
        Fill the matrices and rhs of the problem
        """

        x_nodes = self.x_nodes
        connec = self.connec

        # We use a fix 3-point Gauss integration rule
        wsamples = np.array([1. / 3., 1. / 3., 1. / 3.])
        xsamples = np.array([[2. / 3.0, 1. / 6.], [1. / 6., 2. / 3.], [1. / 6., 1. / 6.]])
        bf = np.array([[1. - x[0] - x[1], x[0], x[1]] for x in xsamples])
        Dbf = np.array([[-1.0, -1.0], [1.0, 0.0],
                        [0.0, 1.0]])  # Derivatives of the basis functions in local coordinates of the triangle

        nPts = len(x_nodes)  # Nymber of points in the discretisation

        # Auxiliary variables to construct the sparse matrices using coo_matrix
        rows = []  # Row labels for every element added to the matrix
        cols = []  # Column labels for every element added to the matrix
        vals1 = []  # Values for matrix A^K
        vals2 = []  # Values for matrix A^{\bar{K}}
        vals3 = []  # Values for matrix A^k

        # Right-hand sides
        rhs1 = np.zeros(nPts * 2)
        rhs2 = np.zeros(nPts * 2)

        # Loop in the elements of the mesh
        for element in connec:
            # [1] Initial computations

            # [1.1] Positions (x,y) of the nodes in the element
            x_nbors = x_nodes[element, 0:2]  # (3*2 matrix)

            # [1.2] Transformation matrix T (converts barycentric coordinates (s1,s2) to (x,y))
            T = x_nbors.T.dot(Dbf)
            det = T[0, 0] * T[1, 1] - T[1, 0] * T[0, 1]  # Jacobian of T (size of the triangle)
            iT = np.array([[T[1, 1], -T[0, 1]], [-T[1, 0], T[0, 0]]]) / det  # Inverse used to convert gradients

            Dbf_glo = Dbf.dot(iT)  # Derivatives of the basis functions in (x,y)

            # [2] Compute terms to matrices and rhs in this element

            # [2.1] A^K
            A1 = np.einsum('ij,kl->ikjl', Dbf_glo.dot(Dbf_glo.T),
                           np.identity(2))  # \nabla_i N^a \nabla_i N^b \delta_{jk}
            A1 += np.einsum('ij,kl->ilkj', Dbf_glo, Dbf_glo)  # \nabla_k N^a \nabla_j N^b
            A1 -= np.einsum('ij,kl->ijkl', Dbf_glo, Dbf_glo)  # -\nabla_j N^a \nabla_k N^b

            # [2.2] A^{\bar{K}}
            A2 = np.einsum('ij,kl->ijkl', Dbf_glo, Dbf_glo)  # \nabla_j N^a \nabla_k N^b

            # [2.3] A^k
            A3 = np.einsum('G,Gi,Gj,kl->ikjl', wsamples, bf, bf, np.identity(2))  # N^a N^b \delta_{jk}

            # [2.4] B^x
            B1 = -Dbf_glo.dot(np.array([[1.0, 0.0], [0.0, 0.0]]))  # \nabla_j N^a \delta_{j1} \delta_{k1}

            # [2.5] B^y
            B2 = -Dbf_glo.dot(np.array([[0.0, 0.0], [0.0, 1.0]]))  # \nabla_j N^a \delta_{j2} \delta_{k2}

            # [3] Assemble the element matrices and rhs

            # [3.1] Lists with rows and columns for each element in A1,A2,A3
            rr = np.array([2 * i + s for i in element for s in np.arange(2)])
            rows.append(np.outer(rr, np.ones(rr.size)))
            cols.append(np.outer(np.ones(rr.size), rr))

            # [3.2] Append values to matrices
            vals1.append(det * A1.flatten())
            vals2.append(det * A2.flatten())
            vals3.append(det * A3.flatten())

            # [3.3] Fill right hand side
            rhs1[rr] += det * B1.flatten()
            rhs2[rr] += det * B2.flatten()

        # After the loop, construct the matrices using coo_matrix
        rows = np.concatenate(rows).flatten()
        cols = np.concatenate(cols).flatten()
        vals1 = np.concatenate(vals1)
        vals2 = np.concatenate(vals2)
        vals3 = np.concatenate(vals3)

        mat1 = coo_matrix((vals1, (rows, cols)), shape=(nPts * 2, nPts * 2))
        mat2 = coo_matrix((vals2, (rows, cols)), shape=(nPts * 2, nPts * 2))
        mat3 = coo_matrix((vals3, (rows, cols)), shape=(nPts * 2, nPts * 2))

        self.mat_K = mat1.tocsc()
        self.mat_bK = mat2.tocsc()
        self.mat_k = mat3.tocsc()
        self.rhs_zx = rhs1
        self.rhs_zy = rhs2

    def interpolate_data(self, points, disp, stretch=None, shear=None, on_deformed=True, kernel_radius=3.0,
                         kernel_sharpness=1.0, remove_null=True, return_deformed=True):

        """
        Interpolate data on the mesh. This can be used for comparing experimental and simulation data together.
        on_deform selects whether to interpolate the data on the original or on the deformed configuration
        """

        if stretch is None:
            stretch = np.array([])

        if shear is None:
            shear = np.array([])

        # We need nPts*3 data for vtk
        points_ = points
        if points.shape[1] == 2:
            points_ = np.zeros([points.shape[0], 3])

        # Deform data if needed
        if on_deformed:
            points_[:,0:2] = points + disp

        # Generate vtk points structure with all data
        vtkpoints = vtk.vtkPoints()
        vtkpoints.SetData(numpy_support.numpy_to_vtk(points_))

        vtkdisp = vtk.vtkDoubleArray()
        vtkdisp.SetName("disp")
        vtkdisp.SetNumberOfComponents(disp.shape[1])
        vtkdisp.SetNumberOfTuples(disp.shape[0])
        vtkdisp.SetVoidArray(disp.flatten(), disp.size, 1)

        vtkstretch = vtk.vtkDoubleArray()
        vtkstretch.SetName("stretch")
        vtkstretch.SetNumberOfComponents(1)
        vtkstretch.SetNumberOfTuples(stretch.size)
        vtkstretch.SetVoidArray(stretch, stretch.size, 1)

        shear = shear.reshape([-1, 4])
        vtkshear = vtk.vtkDoubleArray()
        vtkshear.SetName("shear")
        vtkshear.SetNumberOfComponents(shear.shape[1])
        vtkshear.SetNumberOfTuples(shear.shape[0])
        vtkshear.SetVoidArray(shear, shear.size, 1)

        vtkpointset = vtk.vtkPolyData()
        vtkpointset.SetPoints(vtkpoints)
        vtkpointset.GetPointData().AddArray(vtkdisp)
        vtkpointset.GetPointData().AddArray(vtkstretch)
        vtkpointset.GetPointData().AddArray(vtkshear)

        # Build the locator for the interpolation
        locator = vtk.vtkStaticPointLocator()
        locator.SetDataSet(vtkpointset)
        locator.BuildLocator()

        # Build the Gaussian kernel
        kernel = vtk.vtkGaussianKernel()
        kernel.SetKernelFootprint(0)
        kernel.SetRadius(kernel_radius)
        kernel.SetSharpness(kernel_sharpness)

        # If the interpolation is performed on the deformed configuration, warp the data
        if on_deformed:
            disp_ = np.zeros(self.x_nodes.shape)
            disp_[:,0:2] = self.disp

            warpData = vtk.vtkDoubleArray()
            warpData.SetName("warp")
            warpData.SetNumberOfComponents(3)
            warpData.SetNumberOfTuples(self.x_nodes.shape[0])
            warpData.SetVoidArray(disp_, self.x_nodes.shape[0], 1)

            self.vtkmesh.GetPointData().AddArray(warpData)
            self.vtkmesh.GetPointData().SetActiveVectors(warpData.GetName())

            warpVector = vtk.vtkWarpVector()
            warpVector.SetInputData(self.vtkmesh)
            warpVector.Update()

            self.vtkmesh = warpVector.GetOutput()

        coarseInterpolator = vtk.vtkPointInterpolator()
        coarseInterpolator.SetSourceData(vtkpointset)
        coarseInterpolator.SetInputData(self.vtkmesh)
        coarseInterpolator.SetKernel(kernel)
        coarseInterpolator.SetLocator(locator)
        coarseInterpolator.PassPointArraysOff()
        coarseInterpolator.SetNullPointsStrategyToMaskPoints() # Get points with an invalid interpolation
        coarseInterpolator.Update()

        vtkmesh = coarseInterpolator.GetOutput()

        if return_deformed:
            self.x_nodes[:,0:2] += self.disp

        self.disp = numpy_support.vtk_to_numpy(vtkmesh.GetPointData().GetArray("disp"))
        self.stretch = numpy_support.vtk_to_numpy(vtkmesh.GetPointData().GetArray("stretch"))
        self.shear = numpy_support.vtk_to_numpy(vtkmesh.GetPointData().GetArray("shear"))

        if remove_null:
            not_null = (numpy_support.vtk_to_numpy(vtkmesh.GetPointData().GetArray(
                coarseInterpolator.GetValidPointsMaskArrayName()))).astype(bool)

            idlist_old = np.arange(self.x_nodes.shape[0])[not_null]

            self.x_nodes = self.x_nodes[not_null]
            self.disp = self.disp[not_null]
            self.stretch = self.stretch[not_null]
            self.shear = self.shear[not_null]

            idlist_new = np.arange(self.x_nodes.shape[0])

            c1 = np.in1d(self.connec[:, 0], idlist_old)
            c2 = np.in1d(self.connec[:, 1], idlist_old)
            c3 = np.in1d(self.connec[:, 2], idlist_old)
            self.connec = self.connec[np.logical_and(np.logical_and(c1,c2),c3)]

            for i in idlist_new:
                self.connec[self.connec == idlist_old[i]] = i

            points = vtk.vtkPoints()
            points.SetData(numpy_support.numpy_to_vtk(self.x_nodes))
            vtkmesh.SetPoints(points)

            cells = vtk.vtkCellArray()
            cells.SetCells(self.connec.shape[0], numpy_support.numpy_to_vtkIdTypeArray(self.connec))
            vtkmesh.SetPolys(cells)

        self.shear = self.shear.reshape([-1, 2, 2])

    def load_vtk(self):

        """
        Loads a vtk file. Returns x_nodes (nPts*3 array), connec (nElem*3 array) and a vtkmesh object
        """

        vtkreader = vtk.vtkPolyDataReader()

        vtkreader.SetFileName(self.filemesh)
        vtkreader.Update()
        vtkmesh = vtkreader.GetOutput()

        idList = vtk.vtkIdList()

        nPts = vtkmesh.GetNumberOfPoints()
        points = vtkmesh.GetPoints()
        x_nodes = np.zeros([nPts, 3])
        for i in range(nPts):
            x_nodes[i] = points.GetPoint(i)

        nElem = vtkmesh.GetNumberOfCells()
        connec = np.zeros([nElem, 3], dtype=int)

        for i in range(nElem):
            vtkmesh.GetCellPoints(i, idList)
            connec[i, 0] = idList.GetId(0)
            connec[i, 1] = idList.GetId(1)
            connec[i, 2] = idList.GetId(2)

        self.x_nodes = x_nodes
        self.connec = connec
        self.vtkmesh = vtkmesh

    def save_vtk(self, filename):

        """
        Save results in a vtk file
        """

        vtkwriter = vtk.vtkPolyDataWriter()

        vtkdisp = vtk.vtkDoubleArray()
        vtkdisp.SetNumberOfComponents(2)
        for d in self.disp:
            vtkdisp.InsertNextTuple(d)
        vtkdisp.SetName("Displacement")
        self.vtkmesh.GetPointData().AddArray(vtkdisp)

        if (self.stretch.size > 0):
            vtkdila = vtk.vtkDoubleArray()
            vtkdila.SetNumberOfComponents(1)
            for d in self.stretch:
                vtkdila.InsertNextValue(d)
            vtkdila.SetName("stretch")
            self.vtkmesh.GetPointData().AddArray(vtkdila)

        if (self.shear.size > 0):
            vtkshear = vtk.vtkDoubleArray()
            vtkshear.SetNumberOfComponents(4)
            for d in self.shear:
                vtkshear.InsertNextTuple(d.flatten())
            vtkshear.SetName("Shear")
            self.vtkmesh.GetPointData().AddArray(vtkshear)

        vtkwriter.SetFileName(filename)
        vtkwriter.SetInputData(self.vtkmesh)
        vtkwriter.Update()

    def plot(self, deformed=True, show_mesh=False, mesh_params=None, show_stretch=False, stretch_params=None,
             show_shear=False, shear_params=None):
        '''
        Add plot to matplotlib axis
        '''

        ax = plt.gca()

        if mesh_params is None:
            mesh_params = {"c": "gray", "linewidth": 1, "zorder": 1}

        if stretch_params is None:
            stretch_params = {"levels": 100, "zorder": 0}

        if shear_params is None:
            shear_params = {"headaxislength": 0, "headlength": 0, "pivot": "middle", "width": 0.01, "zorder": 2,
                            "scale": 10.0}


        x = np.copy(self.x_nodes)

        if deformed:
            x[:,0:2] += self.disp

        if show_mesh:
            im = ax.triplot(x[:,0], x[:,1], self.connec, **mesh_params)

        if show_stretch:
            im = ax.tricontourf(x[:, 0], x[:, 1], self.connec, self.stretch, **stretch_params)
            ax.tricontour(x[:, 0], x[:, 1], self.connec, self.stretch, **stretch_params, linewidths=0.3)

        if show_shear:
            S = 2.0 * np.sqrt(self.shear[:, 0, 0]**2 + self.shear[:, 0, 1]**2)
            nx = np.sqrt(self.shear[:, 0, 0] / S + 0.5)
            ny = self.shear[:, 0, 1]/(S * nx)

            nx *= S
            ny *= S

            ax.quiver(x[:,0], x[:,1], nx, ny, **shear_params)

        return im