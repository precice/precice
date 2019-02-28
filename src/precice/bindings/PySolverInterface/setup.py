from setuptools import setup

setup(
    name='PySolverInterface',
    version='0.1',
    description='Python language bindings for preCICE coupling library',
    url='https://github.com/precice/precice',
    author='the preCICE developers',
    author_email='info@precice.org',
    license='LGPL-3.0',
    python_requires='>=3',
    install_requires=[
        'precice'
    ]
)
