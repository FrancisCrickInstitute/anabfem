# Anabfem
### Alejandro Torres-Sanchez (The Francis Crick Institute)

![ScreenShot](doc/example.png?raw=true)

## Summary
This small python project provides tools for the **analysis of annular ablation experiments using a finite element approach**.

It contains the packages:

1. **fem**: containing the object **FEM2DActiveElastic**. This object uses a triangular mesh (either provided by the user or loaded from the set of
meshes provided with the package) to calculate the tissue recoil after ablation assuming the tissue behaves as a 2D 
active, linear elastic solid. Model parameters (bulk modulus, active tensions and elastic constant for adhesion relative
 to the shear modulus) are provided by the user to compute the solution. The solution can be exported in a vtk format, 
 which can be loaded in Paraview for visualisation.
1. **optimization**: containing the object **ParameterOptimizer**. This object uses FEM2DActiveElastic and user-provided data to fit the model parameters. The
data must be provided in a csv file and may contain data regarding tissue displacements (such as cell centre 
displacements, or vertex displacements if the tissue has been triangulated), as well as tissue stretching and 
orientation. 

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
fem.update_mesh_displacements(True, True) #Compute displacement (and stretch and shear)
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
expfiles = ["data_16hAPF.csv", "data_21hAPF.csv", "data_26hAPF.csv","data_31hAPF.csv"]

fem = FEM2DActiveElastic("circle_0.vtk", lintrans=np.array([[45.54/2.0, 0.0],[0.0, 45.54/2.0]]))

opt = ParameterOptimizer(fem, expfiles, kdisp=1.0/(45.54*45.54), kshear=1.0, kstretch=1.0)
opt.parameters = np.array([ 3.6651338743344142e-02,  3.4416483373974804e-01,  3.9598881253401425e-01,
  3.1549244529775110e-02,  3.2515831957694152e-01,  2.7603108245392211e-01,
  9.1302426300089532e-03,  1.2222134424818480e-01,  6.8860193147434345e-01,
  3.7447641907539177e-03, -7.9378001681324278e-02,  6.4175591987936942e-01,
  3.2286999799468258e-01])

opt.optimize()
print(opt.parameters)

import matplotlib.pyplot as plt
fig = plt.figure(figsize=(3,3))
ax = fig.add_subplot(111)
fem.plot(deformed=True, show_mesh=False, show_stretch=True, show_shear=True)
plt.show()
```


