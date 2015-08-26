from distutils.core import setup
import os
import sys

if __name__ == '__main__':

    try:
        from distutils.command.build_scripts import \
                build_scripts_2to3 as build_scripts
    except ImportError:
        from distutils.command.build_scripts import build_scripts

    scripts = ['pdb4amber']

    setup(name='pdb4amber',
          version='15.0',
          description='PDB file preparation script for Amber',
          author='Romain Wolf, Pawel Janowski, and Jason Swails',
          license='GPL v2 or later',
          cmdclass={'build_scripts' : build_scripts},
          scripts=scripts,
    )
