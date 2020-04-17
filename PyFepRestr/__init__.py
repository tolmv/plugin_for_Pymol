# coding=utf-8
"""
PyFepRestr
PyFepRestr is a plug-in designed to set restraints on ligand atoms in the active center of a protein,
as suggested by Boresch et al. and used in free energy calculations with the Gromacs molecular dynamics program.
The plug-in allows you to determine from the coordinates of three ligand atoms and three protein atoms,
the interatomic distance, the two angles between bonds and the three dihedral angles used in the Boresch formula
as well as write their values​to the Gromacs topology file.
(c) Alexander Lashkov, Ivan Tolmachev, Sergey Rubinsky
License:MIT License
"""
from __future__ import absolute_import

from .output import Output
from .tkinput import Restraints
from .wizard import RestraintWizard


def __init_plugin__(app):
    app.menuBar.addmenuitem('Plugin', 'command',
                            label='PyFepRestr',
                            command=lambda: main(app.root))


def main(parent):
    bondForceParams = {'T': None, 'K_r_aA': None,
                       'K_th_a': None, 'K_th_A': None,
                       'K_phi_ba': None, 'K_phi_aA': None, 'K_phi_AB': None,
                       'r_aA': None, 'th_a': None, 'th_A': None,
                       'phi_ba': None, 'phi_aA': None, 'phi_AB': None,
                       'index_a': None, 'index_b': None, 'index_c': None,
                       'index_A': None, 'index_B': None, 'index_C': None
                       }
    atoms_def = {
        'index_c': None,
        'index_b': None,
        'index_a': None,
        'index_A': None,
        'index_B': None,
        'index_C': None
    }
    Restraints(parent, bondForceParams, atoms_def)
