# coding=utf-8
from pymol import cmd, CmdException
from pymol.wizard import Wizard


class RestraintWizard(object, Wizard):
    def __init__(self):
        Wizard.__init__(self)

        # some attributes to do with picking
        self.pick_count = 0
        self.object_prefix = "pw"

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

            point1 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "0"))
            point2 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "1"))
            point3 = cmd.get_atom_coords("(%s%s)" % (self.object_prefix, "2"))
            self.pick_count = 0
            self.reset()

    def get_panel(self):
        return [
            [1, 'Plane Wizard', ''],
            [2, 'Reset', 'cmd.get_wizard().reset()'],
            [2, 'Done', 'cmd.set_wizard()'],
        ]
