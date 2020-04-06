# coding=utf-8
from __future__ import absolute_import
from __future__ import print_function

from pymol import cmd, CmdException
from pymol.wizard import Wizard


class RestraintWizard(object, Wizard):
    def __init__(self, bondForceParams):
        Wizard.__init__(self)
        self.bondForceParams = bondForceParams
        # some attributes to do with picking
        self.pick_count = 0
        self.object_prefix = "pw"
        self.error = None
        self.iswait = True
        self.selection_mode = cmd.get_setting_legacy("mouse_selection_mode")
        cmd.set("mouse_selection_mode", 0)  # set selection mode to atomic
        cmd.deselect()

    def reset(self):
        cmd.delete(self.object_prefix + "*")
        cmd.delete("sele*")
        cmd.delete("_indicate*")
        cmd.unpick()
        cmd.refresh_wizard()

    def cleanup(self):
        cmd.set("mouse_selection_mode", self.selection_mode)  # restore selection mode
        self.reset()

    def get_prompt(self):
        self.prompt = None
        if self.pick_count == 0:
            self.prompt = ['Please click on the first atom...']
        elif self.pick_count == 1:
            self.prompt = ['Please click on the second atom...']
        elif self.pick_count == 2:
            self.prompt = ['Please click on the third atom...']
        elif self.pick_count == 3:
            self.prompt = ['Please click on the fourth atom...']
        elif self.pick_count == 4:
            self.prompt = ['Please click on the fifth atom...']
        elif self.pick_count == 5:
            self.prompt = ['Please click on the sixth atom...']
        return self.prompt

    def do_select(self, name):
        # "edit" only this atom, and not others with the object prefix
        try:
            cmd.edit("%s and not %s*" % (name, self.object_prefix))
            self.do_pick(0)
        except CmdException, pmce:
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

            # point1 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "0"))
            # point2 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "1"))
            # point3 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "2"))
            r_aA = cmd.get_distance(atom1="pw2", atom2="pw3", state=0)
            th_a = cmd.get_angle(atom1="pw1", atom2="pw2", atom3="pw3", state=0)
            th_A = cmd.get_angle(atom1="pw2", atom2="pw3", atom3="pw4", state=0)
            phi_ba = cmd.get_dihedral(atom1="pw0", atom2="pw1", atom3="pw2", atom4="pw3", state=0)
            phi_aA = cmd.get_dihedral(atom1="pw1", atom2="pw2", atom3="pw3", atom4="pw4", state=0)
            phi_AB = cmd.get_dihedral(atom1="pw2", atom2="pw3", atom3="pw4", atom4="pw5", state=0)
            index_c = cmd.id_atom("pw0")
            index_b = cmd.id_atom("pw1")
            index_a = cmd.id_atom("pw2")
            index_A = cmd.id_atom("pw3")
            index_B = cmd.id_atom("pw4")
            index_C = cmd.id_atom("pw5")
            self.pick_count = 0
            self.reset()
            self.SetBondForceParam(r_aA, th_a, th_A, phi_ba, phi_aA, phi_AB,
                                   index_c, index_b, index_a, index_A, index_B, index_C)
            cmd.set_wizard()
            self.iswait = False

    def get_panel(self):
        return [
            [1, 'PyFepRestr', ''],
            [2, 'Reset', 'cmd.get_wizard().reset()'],
            # [2, 'Done', 'cmd.set_wizard()'],
        ]

    def SetBondForceParam(self, *args):
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
