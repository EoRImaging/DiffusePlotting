from setuptools import setup
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

extensions = [
    Extension("diffuseplotting.lic_internal", ["diffuseplotting/lic_internal.pyx"],
              include_dirs=[numpy.get_include()])
]

setup_args = {
    'name': 'diffuseplotting',
    'package_dir': {'diffuseplotting': 'diffuseplotting'},
    'packages': ['diffuseplotting'],
    'ext_modules': cythonize(extensions)
}

if __name__ == '__main__':
    apply(setup, (), setup_args)
