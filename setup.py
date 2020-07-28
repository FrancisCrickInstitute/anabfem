from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()

# Get the long description from the README file
long_description = (here / 'README.md').read_text(encoding='utf-8')

setup(
    name='anabfem',
    version='0.0.1',
    description='Finite element analysis for annular ablation experiments',
    long_description_content_type='text/markdown',
    url='https://github.com/pypa/sampleproject',  # FIXME: change this!
    author='Alejandro Torres-Sanchez',
    author_email='torres.sanchez.a@gmail.com',
    python_requires='>=3.5, <4',
    install_requires=['vtk', 'numpy', 'scipy'],
    packages=['anabfem'],
    package_dir={'anabfem': 'anabfem'},
    package_data={'anabfem': ['meshes/*.vtk']}
)
