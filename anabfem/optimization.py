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

import numpy as np
from scipy.optimize import minimize


class ParameterOptimizer:

    """
    This class optimises the parameters of a series of annular ablation experiments
    using the FEM2DActiveElastic class.
    """

    def __init__(self, fem, csvfiles, csv_delim=";", skiprows=1, kdisp=1.0, kstretch=0.0, kshear=0.0, stretch_shear=True,
                 shared_param=np.array([False, False, False, True], dtype=bool)):
        """
        Initialise the optimiser. This requires passing a FEM2DActiveElastic object and the files where experimental
        values are stored.

        Load data from csv with assuming:
             - column 0: original x position
             - column 1: final x position
             - column 2: original y position
             - column 3: final y position

        if stretch_shear:
             - column 4: original area
             - column 5: final area
             - column 6  : original Qxx elongation
             - column 7: final Qxx elongation
             - column 8: original Qxy elongation
             - column 9: final Qxy elongation
        """

        self.fem = fem

        self.params_per_file = 4
        self.shared_param = shared_param
        self.n_shared = shared_param.sum()
        self.n_notshared = self.params_per_file - self.n_shared

        # Map to build parameters for each file taking into account shared and not shared parameters among files
        self.map2params = np.zeros([len(csvfiles),4], dtype=int)
        for n in np.arange(len(csvfiles)):
            p_s = 0
            p_n = 0
            for p in np.arange(self.params_per_file):
                if self.shared_param[p]:
                    self.map2params[n,p] = -p_n-1
                    p_n += 1
                else:
                    self.map2params[n,p] = n * self.n_notshared + p_s
                    p_s += 1

        # 3 parameters per file (k, \zeta_x, \zeta_y) and a global parameter \bar{K}
        self.parameters = np.ones([self.n_notshared * len(csvfiles) + self.n_shared])

        self.stretch_shear = stretch_shear

        self.kdisp = kdisp
        self.kstretch = kstretch
        self.kshear = kshear

        if self.kstretch != 0 and (not self.stretch_shear):
            print("ParameterOptimizer: Warning! Boolean flag stretch_shear is false but kstretch is not 0 "
                  "(kstretch = {0:.2e}".format(self.kstretch))

        if self.kshear != 0 and (not self.stretch_shear):
            print("ParameterOptimizer: Warning! Boolean flag stretch_shear is false but kshear is not 0 "
                  "(kshear = {0:.2e}".format(kshear))

        self.x0 = []
        self.u = []
        self.stretch = []
        self.qxx = []
        self.qxy = []

        self.stretch_shear = stretch_shear

        for expfile in csvfiles:
            expdata = np.loadtxt(expfile, delimiter=csv_delim, skiprows=skiprows)

            self.x0.append(expdata[:, [0, 2]])
            self.u.append(expdata[:, [1, 3]] - self.x0[-1])

            if stretch_shear:
                self.stretch.append((expdata[:, 5] - expdata[:, 4]) / expdata[:, 4])
                self.qxx.append(expdata[:, 7] - expdata[:, 6])
                self.qxy.append(expdata[:, 9] - expdata[:, 8])


    def target_func(self, parameters, return_residuals=False):

        """
        Function to minimise (least square error in displacements, stretch and shear)
        """

        func = 0.0
        gradient = np.zeros(len(parameters))

        # If the function is used for bootstrapping, then it returns residuals and predicted values
        if return_residuals:
            pred_u = []
            pred_str = []
            pred_qxx = []
            pred_qxy = []

            res_u = []
            res_dil = []
            res_qxx = []
            res_qxy = []

        for n in np.arange(len(self.x0)):

            # Generate the list of parameters in right order for FEM evaluation
            parameters_ = parameters[self.map2params[n]]

            # Update the list in the FEM object, update vertex displacements and compute expected displacements
            # (and shear and dilatation if needed)
            self.fem.parameters = parameters_
            self.fem.update_mesh_displacements(True, True)
            disp, stretch, shear, dbK_displ, dk_displ, dsx_displ, dsy_displ, dbK_dilatation, dk_dilatation, \
                dsx_dilatation, dsy_dilatation, dbK_shear, dk_shear, dsx_shear, dsy_shear = \
                self.fem.compute_point_displacements(self.x0[n], True, True)

            # Compute differences between experiment and model
            res_u_ = self.u[n] - disp
            res_str_ = self.stretch[n] - stretch
            res_qxx_ = self.qxx[n] - shear[:, 0, 0]
            res_qxy_ = self.qxy[n] - shear[:, 0, 1]

            if return_residuals:
                pred_u.append(disp)
                pred_str.append(stretch)
                pred_qxx.append(shear[:, 0, 0])
                pred_qxy.append(shear[:, 0, 1])

                res_u.append(res_u_)
                res_dil.append(res_str_)
                res_qxx.append(res_qxx_)
                res_qxy.append(res_qxy_)

            res_u_ = res_u_.flatten()

            # Compute target function
            func += 0.5 * (self.kdisp * res_u_.dot(res_u_) + self.kstretch * res_str_.dot(res_str_)
                           + self.kshear * (res_qxx_.dot(res_qxx_) + res_qxy_.dot(res_qxy_)))

            # Compute gradients
            gradient[self.map2params[n]] -= self.kdisp * np.array(
                [res_u_.dot(dk_displ.flatten()), res_u_.dot(dsx_displ.flatten()),
                 res_u_.dot(dsy_displ.flatten()), res_u_.dot(dbK_displ.flatten())])

            gradient[self.map2params[n]] -= self.kstretch * np.array(
                [res_str_.dot(dk_dilatation.flatten()), res_str_.dot(dsx_dilatation.flatten()),
                 res_str_.dot(dsy_dilatation.flatten()), res_str_.dot(dbK_dilatation.flatten())])

            gradient[self.map2params[n]] -= self.kshear * np.array(
                [res_qxx_.dot(dk_shear[:, 0, 0].flatten()), res_qxx_.dot(dsx_shear[:, 0, 0].flatten()),
                 res_qxx_.dot(dsy_shear[:, 0, 0].flatten()), res_qxx_.dot(dbK_shear[:, 0, 0].flatten())])

            gradient[self.map2params[n]] -= self.kshear * np.array(
                [res_qxy_.dot(dk_shear[:, 0, 1].flatten()), res_qxy_.dot(dsx_shear[:, 0, 1].flatten()),
                 res_qxy_.dot(dsy_shear[:, 0, 1].flatten()), res_qxy_.dot(dbK_shear[:, 0, 1].flatten())])

        if return_residuals:
            return func, gradient, pred_u, pred_str, pred_qxx, pred_qxy, res_u, res_dil, res_qxx, res_qxy
        else:
            return func, gradient

    def optimize(self, method='L-BFGS-B', options=None, print_results=True):

        """
        Parameter optimisation

        This method provides an envelop for the minimize function in scipy.
        """

        if options is None:
            options = {'ftol': 1.E-14, 'gtol': 1E-14, 'maxiter': 5000, 'iprint': 50}

        res = minimize(self.target_func, self.parameters, jac=True, method=method, options=options)

        if print_results:
            print(res)

        self.parameters = res.x

    def bootstrap(self, nbootstrap, opt_method='L-BFGS-B', opt_options=None, opt_print_results=True):

        """
        Bootstrapping analysis

        This method performs a bootstrapping analysis to check the stability of the optimal parameters.
        For that, this functions uses the numerical solution obtained by FEM that optimizes the fit to the experimental
        data as the true solution and assumes that the measured residuals represent well the "population" of residuals.
        Then one can produce a large number of different samples by assigning to every point the FEM solution plus a
        residual from the list of residuals (with re-emplacement).
        """

        # Save the loaded data
        u = np.copy(self.u)
        stretch = np.copy(self.stretch)
        qxx = np.copy(self.qxx)
        qxy = np.copy(self.qxy)
        parameters = np.copy(self.parameters)

        # Optimise and obtain the residuals
        self.optimize()

        energy, gradient, pred_u, pred_stretch, pred_qxx, pred_qxy, \
            res_u, res_stretch, res_qxx, res_qxy = self.target_func( self.parameters, return_residuals=True)

        # Do the bootstrapping analysis
        param_list = []
        for _ in np.arange(nbootstrap):

            # Generate new residuals by randomly selecting them from the initial list of residuals
            for m in np.arange(len(res_u)):
                res_u_ = res_u[m][np.random.randint(0, res_stretch[m].size, size=res_stretch[m].size), :]
                res_str_ = res_stretch[m][np.random.randint(0, res_stretch[m].size, size=res_stretch[m].size)]
                res_qxx_ = res_qxx[m][np.random.randint(0, res_stretch[m].size, size=res_stretch[m].size)]
                res_qxy_ = res_qxy[m][np.random.randint(0, res_stretch[m].size, size=res_stretch[m].size)]

                # Generate a new sample using the previously calculated residuals
                self.u[m] = pred_u[m] + res_u_
                self.stretch[m] = pred_stretch[m] + res_str_
                self.qxx[m] = pred_qxx[m] + res_qxx_
                self.qxy[m] = pred_qxy[m] + res_qxy_

            self.optimize(method=opt_method, options=opt_options, print_results=opt_print_results)
            param_list.append(self.parameters)
            self.parameters = parameters

        # Restore experimental values from files
        self.u = u
        self.stretch = stretch
        self.qxx = qxx
        self.qxy = qxy

        return np.array(param_list)
