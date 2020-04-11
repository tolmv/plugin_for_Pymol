# coding=utf-8
from __future__ import absolute_import
from __future__ import print_function

from pymol import cmd, CmdException
from pymol.wizard import Wizard

from .output import Output

try:
    from Tkinter import Toplevel
except ImportError:
    from tkinter import Toplevel


def getAtomString(sel):
    atoms = cmd.get_model(sel)
    s = ""
    for at in atoms.atom:
        s += (at.name + '_' + at.resn + str(at.resi) + "/" + at.chain)
    return s


class RestraintWizard(Wizard):
    def __init__(self, parent, bondForceParams, atoms_def):
        Wizard.__init__(self)
        self.parent = parent
        self.atoms_def = atoms_def
        self.params_str = ['r_aA', 'th_a', "th_A1'", 'phi_ba', 'phi_aA', 'phi_AB']
        self.indexes_list = ['atom_c', 'atom_b', 'atom_a', "atom_A1", "atom_B1", "atom_C1"]
        self.bondForceParams = bondForceParams
        # some attributes to do with picking
        self.pick_count = 0
        self.object_prefix = "pw"
        self.error = None
        self.iswait = True
        self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
        cmd.set("mouse_selection_mode", 0)  # set selection mode to atomic
        cmd.deselect()

    def reset(self, reset=True):
        if reset:
            for par, ind in zip(self.params_str, self.indexes_list):
                cmd.delete(par)
                cmd.delete(ind)
        cmd.delete(self.object_prefix + "*")
        cmd.delete("sele*")
        cmd.delete("_indicate*")
        cmd.unpick()
        self.pick_count = 0
        cmd.refresh_wizard()

    def cleanup(self):
        cmd.set("mouse_selection_mode", self.selection_mode)  # restore selection mode
        self.reset(reset=False)

    def get_prompt(self):
        self.prompt = None
        if self.pick_count == 0:
            self.prompt = ['Please click on the first (c) atom...']
        elif self.pick_count == 1:
            self.prompt = ['Please click on the second (b) atom...']
        elif self.pick_count == 2:
            self.prompt = ['Please click on the third (a) atom...']
        elif self.pick_count == 3:
            self.prompt = ['Please click on the fourth (A) atom...']
        elif self.pick_count == 4:
            self.prompt = ['Please click on the fifth (B) atom...']
        elif self.pick_count == 5:
            self.prompt = ['Please click on the sixth (C) atom...']
        return self.prompt

    def do_select(self, name):
        # "edit" only this atom, and not others with the object prefix
        try:
            cmd.edit("%s and not %s*" % (name, self.object_prefix))
            self.do_pick(0)
        except CmdException as pmce:
            print(pmce)

    def pickNextAtom(self, atom_name):
        # transfer the click selection to a named selection
        cmd.select(atom_name, "(pk1)")
        # delete the click selection
        cmd.unpick()
        # using the magic of indicate, highlight stuff
        indicate_selection = "_indicate" + self.object_prefix
        cmd.select(indicate_selection, atom_name)
        cmd.enable(indicate_selection)
        view = cmd.get_view()
        cmd.create(self.indexes_list[self.pick_count], self.object_prefix + str(self.pick_count))
        cmd.set_view(view)
        cmd.set("sphere_scale", '0.3', self.indexes_list[self.pick_count])
        cmd.show_as('spheres', self.indexes_list[self.pick_count])
        if self.pick_count < 3:
            cmd.color("red", self.indexes_list[self.pick_count])
        else:
            cmd.color("green", self.indexes_list[self.pick_count])
        self.pick_count += 1
        self.error = None
        # necessary to force update of the prompt
        cmd.refresh_wizard()

    def do_pick(self, picked_bond):

        # this shouldn't actually happen if going through the "do_select"
        if picked_bond:
            self.error = "Error: please select an atom, not a bond."
            print(self.error)
            return
        atom_name = self.object_prefix + str(self.pick_count)
        if self.pick_count < 5:
            self.pickNextAtom(atom_name)
        else:
            self.pickNextAtom(atom_name)
            self.doFinish()

    def doFinish(self):
        r_aA = cmd.get_distance(atom1="pw2", atom2="pw3", state=0)
        cmd.distance(self.params_str[0], "pw2", "pw3")
        th_a = cmd.get_angle(atom1="pw1", atom2="pw2", atom3="pw3", state=0)
        cmd.angle(self.params_str[1], "pw1", "pw2", "pw3")
        th_A = cmd.get_angle(atom1="pw2", atom2="pw3", atom3="pw4", state=0)
        cmd.angle(self.params_str[2], "pw2", "pw3", "pw4")
        phi_ba = cmd.get_dihedral(atom1="pw0", atom2="pw1", atom3="pw2", atom4="pw3", state=0)
        cmd.dihedral(self.params_str[3], "pw0", "pw1", "pw2", "pw3")
        phi_aA = cmd.get_dihedral(atom1="pw1", atom2="pw2", atom3="pw3", atom4="pw4", state=0)
        cmd.dihedral(self.params_str[4], "pw1", "pw2", "pw3", "pw4")
        phi_AB = cmd.get_dihedral(atom1="pw2", atom2="pw3", atom3="pw4", atom4="pw5", state=0)
        cmd.dihedral(self.params_str[5], "pw2", "pw3", "pw4", "pw5")
        index_c = cmd.id_atom("pw0")
        index_c_name = getAtomString('pw0')
        index_b = cmd.id_atom("pw1")
        index_b_name = getAtomString('pw1')
        index_a = cmd.id_atom("pw2")
        index_a_name = getAtomString('pw2')
        index_A = cmd.id_atom("pw3")
        index_A_name = getAtomString('pw3')
        index_B = cmd.id_atom("pw4")
        index_B_name = getAtomString('pw4')
        index_C = cmd.id_atom("pw5")
        index_C_name = getAtomString('pw5')
        self.setBondForceParam(r_aA, th_a, th_A, phi_ba, phi_aA, phi_AB,
                               index_c, index_b, index_a, index_A, index_B, index_C)
        self.setAtomsDef(index_c_name, index_b_name, index_a_name, index_A_name, index_B_name, index_C_name)
        top = Toplevel(self.parent)
        Output(top, self.bondForceParams, self.atoms_def)
        cmd.set_wizard()

    def get_panel(self):
        return [
            [1, 'PyFepRestr', ''],
            [2, 'Reset', 'cmd.get_wizard().reset()'],
        ]

    def setAtomsDef(self, *args):
        self.atoms_def['index_c'] = args[0]
        self.atoms_def['index_b'] = args[1]
        self.atoms_def['index_a'] = args[2]
        self.atoms_def['index_A'] = args[3]
        self.atoms_def['index_B'] = args[4]
        self.atoms_def['index_C'] = args[5]

    def setBondForceParam(self, *args):
        self.bondForceParams['r_aA'] = args[0] / 10.0
        self.bondForceParams['th_a'] = args[1]
        self.bondForceParams['th_A'] = args[2]
        self.bondForceParams['phi_ba'] = args[3]
        self.bondForceParams['phi_aA'] = args[4]
        self.bondForceParams['phi_AB'] = args[5]
        self.bondForceParams['index_a'] = args[8]
        self.bondForceParams['index_b'] = args[7]
        self.bondForceParams['index_c'] = args[6]
        self.bondForceParams['index_A'] = args[9]
        self.bondForceParams['index_B'] = args[10]
        self.bondForceParams['index_C'] = args[11]
