from setuptools import setup

setup(
    name='precice',
    version='1.6.1',
    description='Python language bindings for preCICE coupling library',
    url='https://github.com/precice/precice',
    author='the preCICE developers',
    author_email='info@precice.org',
    license='LGPL-3.0',
    python_requires='>=3',
    packages=['precice'],
    install_requires=[
        'precice_future'
    ]
)
