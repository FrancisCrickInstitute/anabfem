# Anabfem

![ScreenShot](doc/example.png?raw=true)

## Summary
This small python project provides tools for the **analysis of annular ablation experiments using a finite element approach**.

It includes the packages:

1. **fem**: containing the object **FEM2DActiveElastic**. This object uses a triangular mesh (either provided by the user or loaded from the set of meshes provided with the package) to calculate the tissue recoil after ablation assuming the tissue behaves as a 2D  active, linear elastic solid adhered to the substrate with elastic bonds. Model parameters (bulk modulus, active tensions and elastic constant for adhesion relative to the shear modulus) are provided by the user to compute the solution. The solution can be plotted with matplotlib with the plot method (see examples below). It can also be exported in a vtk format, which can be loaded in Paraview for visualisation and further manipulation. 
2. **optimization**: containing the object **ParameterOptimizer**. This object uses FEM2DActiveElastic and user-provided data to fit the model parameters via the . The data must be provided in a csv file and may contain data regarding tissue displacements (such as cell centre  displacements, or vertex displacements if the tissue has been triangulated), as well as tissue stretching and orientation. The

## Installation

```bash
pip install .
```

## Usage

1. Finite element solution given a set of parameters:
```python

# Compute FEM solution
from anabfem.fem import FEM2DActiveElastic
fem = FEM2DActiveElastic("circle_0.vtk", lintrans=np.array([[50.0, 0.0],[0.0, 50.0]]))
fem.parameters = [0.1, 0.6, -0.6, 0.1] # [k, \zeta_x, \zeta_y, \bar{K}]/K
fem.update_mesh_displacements(True, True) # Compute displacement (and stretch and shear)
fem.save_vtk("solution.vtk") # Save in vtk format (to open in Paraview, for instance)

# Plot it with matplotlib
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(3,3))
ax = fig.add_subplot(111)
fem.plot(deformed=True, show_mesh=False, show_stretch=True, show_shear=True)
plt.show()

```
2. Optimise parameters given an experimental data set
```python
from anabfem.fem import FEM2DActiveElastic
from anabfem.optimization import ParameterOptimizer

# Data fileds
expfiles = ["data.csv"]

# Initialise the fem object
fem = FEM2DActiveElastic("circle_0.vtk", lintrans=np.array([[25.0, 0.0],[0.0, 25.0]]))

# Initialise the optimiser
opt = ParameterOptimizer(fem, expfiles, kdisp=1.0, kshear=1.0, kstretch=1.0)

# Optimise
opt.optimize()

# Print the parameters
print(opt.parameters)

# Update the fem object with the best parameters
fem.parameters = opt.parameters
fem.update_mesh_displacements(True, True)

#Plot it in matplotlib
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(3,3))
ax = fig.add_subplot(111)
fem.plot(deformed=True, show_mesh=False, show_stretch=True, show_shear=True)
plt.show()
```
## Author and ackowledgement
This code has been produced by Alejandro Torres-Sanchez at the Francis Crick Institute
