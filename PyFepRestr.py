# coding=utf-8
"""
PyFepRestr
(c) Alexander Lashkov, Ivan Tolmachev, Sergey Rubinsky
License:MIT License
"""
from __future__ import absolute_import

import math
import tempfile
import webbrowser
from pymol import cmd, CmdException
from pymol.wizard import Wizard
from sys import platform

try:
    from Tkinter import BooleanVar, Button, END, Entry, Label, Radiobutton, Text, Toplevel, W
    from tkFileDialog import askopenfilename
    from tkFont import Font
    from tkMessageBox import showerror, showinfo
except ImportError:
    from tkinter import BooleanVar, Button, END, Entry, Label, Radiobutton, Text, Toplevel, W
    from tkinter.filedialog import askopenfilename
    from tkinter.font import Font
    from tkinter.messagebox import showerror, showinfo

if platform == "darwin":
    button_font = label_font = radiobutton_font = Font(family='Arial', size=15)
else:
    radiobutton_font = Font(font=Radiobutton()["font"])
    label_font = Font(font=Label()["font"])
    button_font = Font(font=Button()["font"])

help_1 = """<html>
<title>Help</title>
<h1 align="left" style="color: Black">QuickStart</h1>
<h2 align="left" style="color: Black">Unit Button (kJ or kCal)</h2>

<body>
You can choose the desired unit of measurement.</body>

<h2 align="left" style="color: Black">Exit</h2>
<body>
Exit the program.
</body>

<h2 align="left" style="color: Black"> Text Fields</h2>

<body>
 K are the force constants of the harmonic restraints for the one distance (raA),<br>
 two angle (θa and θA), and three dihedral (φba, φaA, and φAB) restraints. Temp and K values must be > 0.<br>
 The value is in the range of 5-50 kCal/mol/rad<sup>2</sup> (or kCal/mol/Å<sup>2</sup> are acceptable.<br>
 Recommended value is 10 kCal/mol/rad<sup>2</sup> (or kCal/mol/Å<sup>2</sup>.<br>

 <h2 align="left" style="color: Black"> Next Button</h2>

 On the next step you choice six atoms (3 for ligand (c-b-a) and 3 for protein (A-B-C)).<br>
 Sequence of selection is c-b-a-A-B-C or C-B-A-a-b-c.
 </body>
</html>

"""

help_2 = """<html>
<title>Help</title>
<h1 align="left" style="color: Black">QuickStart</h1>

<h2 align="left" style="color: Black"> Gromacs topology</h2>

<h3 align="left" style="color: Black"> Preview</h4>

<body>
Preview an additional section in the topology file.
</body>

<h3 align="left" style="color: Black"> Write</h4>

<body>
Select the topology file and writes all the indices of the atoms, <br>
force constants, distances, angles and dihedral angles into it,<br>
where A is &lambda;<sub>restr</sub> =  0 and B is &lambda;<sub>restr</sub> =  1. <br>
The bonded-lambdas vector was interpolated <br>
between the force constant (and equilibrium posi- tions) in state A and B.
</body>

<h2 align="left" style="color: Black">Unit Button (kJ or kCal)</h2>

<body>
You can choose the desired unit of measurement,<br>
depending on what you need to output.</body>
<h2 align="left" style="color: Black">Exit</h2>
<body>
Exit the program.
</body>
<h2 align="left" style="color: Black">Free Energy</h2>
<body>
Here you can see value of &Delta;G.<br>
If you need some more info, see documentation.
</body>

</html>
"""


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


def getHelp(help):
    with tempfile.NamedTemporaryFile('w', delete=False, suffix='.html') as f:
        url = "file://" + f.name
        f.write(help)
    webbrowser.open(url)


def kCal_to_kJ(E):
    return E * 4.1868


def getAtomString(sel):
    atoms = cmd.get_model(sel)
    s = ""
    for at in atoms.atom:
        s += (at.name + '_' + at.resn + str(at.resi) + ("/" + at.chain if at.chain else ""))
    return s


def kJ_to_kCal(E):
    return E / 4.1868


def calc_dG(T, r_aA, th_a, th_A, K_r_aA, K_th_a, K_th_A, K_phi_ba, K_phi_aA, K_phi_AB):
    """BORESCH FORMULA - Calculate dG restraints off"""
    R = 8.314472 * 0.001  # Gas constant in kJ/mol/K
    V = 1.66  # standard volume in nm^3
    th_a = math.radians(th_a)  # convert angle from degrees to radians --> math.sin() wants radians
    th_A = math.radians(th_A)  # convert angle from degrees to radians --> math.sin() wants radians
    dG = - R * T * math.log(
        (8.0 * math.pi ** 2.0 * V) / (r_aA ** 2.0 * math.sin(th_a) * math.sin(th_A))
        *
        (
                ((K_r_aA * K_th_a * K_th_A * K_phi_ba * K_phi_aA * K_phi_AB) ** 0.5) / ((2.0 * math.pi * R * T) ** 3.0)
        )
    )
    return dG


class Restraints(object):
    def __init__(self, parent, bondForceParams, atoms_def):
        self.parent = parent
        self.main = Toplevel(self.parent)
        self.bondForceParams = bondForceParams
        self.atoms_def = atoms_def
        self.main.title('PyFepRestr')
        self.validated_values = []
        self.r_var = BooleanVar()
        self.r_var.set(1)

        rj1 = Radiobutton(self.main, text='kJ', variable=self.r_var, value=0, command=self.refresh,
                          font=radiobutton_font)
        rcal1 = Radiobutton(self.main, text="kCal", variable=self.r_var, value=1, command=self.refresh,
                            font=radiobutton_font)
        rj1.grid(row=0, column=0, padx=5, pady=5)
        rcal1.grid(row=0, column=1, padx=5, pady=5)

        labels = ['Temp',
                  'K raA',
                  u'K \u03b8a', u'K \u03b8A',
                  u'K \u03c6ba', u'K \u03c6aA', u'K \u03c6AB']

        label_all = []
        self.entry_all = []
        self.entry_all_get = []
        self.dimen_all = []
        for lab in labels:
            label_answer = Label(self.main, text=lab, anchor=W, font=label_font)
            label_all.append(label_answer)
            entry = Entry(self.main)
            self.entry_all.append(entry)
            self.entry_all_get.append(entry.get)
            dimen = Label(self.main, anchor=W, font=label_font)
            self.dimen_all.append(dimen)

        for i, (label, entry, dimen) in enumerate(zip(label_all, self.entry_all, self.dimen_all)):
            label.grid(row=i + 1, column=1, padx=5, pady=5, sticky=W)
            entry.grid(row=i + 1, column=2, padx=5, pady=5, sticky=W)
            dimen.grid(row=i + 1, column=3, padx=10, pady=5, sticky=W)

        self.dimen_all[0]['text'] = 'Kelvin'
        self.refresh()

        self.button_res = Button(self.main, text="Next -> ", command=self.validate, font=button_font)
        self.button_res.grid(row=11, column=2, padx=5, pady=5)

        self.destroyProgr = Button(self.main, text='Exit', bg='red', command=self.main.destroy, font=button_font)
        self.destroyProgr.grid(row=0, column=3, padx=5, pady=5)

        self.helpProgr = Button(self.main, text=' ? ', bg='#ffb3fe', command=lambda: getHelp(help_1), font=button_font)
        self.helpProgr.grid(row=12, column=0, padx=5, pady=5)

    def refresh(self):
        if self.r_var.get():
            self.dimen_all[1].configure(text=u'kCal/mol/\u212b\u00b2')
            for dimen in self.dimen_all[2:]:
                dimen.configure(text=u'kCal/mol/rad\u00b2')
        else:
            self.dimen_all[1].configure(text=u'kJ/mol/\u212b\u00b2')
            for dimen in self.dimen_all[2:]:
                dimen.configure(text=u'kJ/mol/rad\u00b2')
        for dimen in self.dimen_all:
            dimen.update()

    def validate(self):
        for i, k in enumerate(self.entry_all_get):
            try:
                f = float(k())
                if f <= 0:
                    raise ValueError
                self.entry_all[i]['bg'] = 'white'
            except ValueError:
                self.entry_all[i]['bg'] = "red"
                return
            self.validated_values.append(f)
        self.finish()

    def finish(self):
        if self.r_var.get():
            self.validated_values = list((self.validated_values[0],)) + \
                                    list(map(kCal_to_kJ, self.validated_values[1:]))

        self.bondForceParams['T'] = self.validated_values[0]  # Temperature (K)
        self.bondForceParams['K_r_aA'] = self.validated_values[1] * 100  # force constant for distance (kJ/mol/nm^2)
        self.bondForceParams['K_th_a'] = self.validated_values[2]  # force constant for angle (kJ/mol/rad^2)
        self.bondForceParams['K_th_A'] = self.validated_values[3]  # force constant for angle (kJ/mol/rad^2)
        self.bondForceParams['K_phi_ba'] = self.validated_values[4]  # force constant for dihedral (kJ/mol/rad^2)
        self.bondForceParams['K_phi_aA'] = self.validated_values[5]  # force constant for dihedral (kJ/mol/rad^2)
        self.bondForceParams['K_phi_AB'] = self.validated_values[6]  # force constant for dihedral (kJ/mol/rad^2)
        showinfo('Info', 'Now choose the atoms you need')
        wiz = RestraintWizard(self.parent, self.bondForceParams, self.atoms_def)
        cmd.set_wizard(wiz)
        cmd.refresh_wizard()
        self.main.withdraw()
        self.main.destroy()


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
            showerror("Error", str(pmce))

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
            showerror("Error", self.error)
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


class Output(object):
    def __init__(self, main, bondForceParams, atoms_def):
        self.bondForceParams = bondForceParams
        self.dG_off_kJ = calc_dG(
            bondForceParams['T'],
            bondForceParams['r_aA'],
            bondForceParams['th_a'],
            bondForceParams['th_A'],
            bondForceParams['K_r_aA'],
            bondForceParams['K_th_a'],
            bondForceParams['K_th_A'],
            bondForceParams['K_phi_ba'],
            bondForceParams['K_phi_aA'],
            bondForceParams['K_phi_AB']
        )
        self.dG_on_kJ = -self.dG_off_kJ
        self.dG_off_kCal = kJ_to_kCal(self.dG_off_kJ)
        self.dG_on_kCal = kJ_to_kCal(self.dG_on_kJ)
        self.atoms_def = atoms_def
        self.main = main
        self.main.title('PyFepRestr')
        self.r_var = BooleanVar()
        self.r_var.set(0)

        rj1 = Radiobutton(self.main, text='kJ', variable=self.r_var, value=0, command=self.refresh,
                          font=radiobutton_font)
        rcal1 = Radiobutton(self.main, text="kCal", variable=self.r_var, value=1, command=self.refresh,
                            font=radiobutton_font)
        rj1.grid(row=0, column=0, padx=5, pady=5)
        rcal1.grid(row=0, column=1, padx=5, pady=5)

        name0 = Label(self.main, text=u'\u0394G_off = ', font=label_font)
        name1 = Label(self.main, text=u'\u0394G_on = ', font=label_font)
        name0.grid(row=1, column=0, padx=5, pady=5)
        name1.grid(row=2, column=0, padx=5, pady=5)

        self.answer0 = Label(self.main, font=label_font)
        self.answer1 = Label(self.main, font=label_font)
        self.answer0.grid(row=1, column=1, padx=5, pady=5)
        self.answer1.grid(row=2, column=1, padx=5, pady=5)

        self.dimen0 = Label(self.main, font=label_font)
        self.dimen1 = Label(self.main, font=label_font)
        self.dimen0.grid(row=1, column=2, padx=5, pady=5)
        self.dimen1.grid(row=2, column=2, padx=5, pady=5)
        self.refresh()

        destroyProgr = Button(self.main, text='Exit', bg='red', command=self.main.destroy, font=button_font)
        destroyProgr.grid(row=0, column=3, padx=5, pady=5)

        helpProgr = Button(self.main, text=' ? ', bg='#ffb3fe', command=lambda: getHelp(help_2), font=button_font)
        helpProgr.grid(row=4, column=0, padx=5, pady=5)

        name3 = Label(self.main, text='Gromacs topology file:', font=label_font)
        name3.grid(row=3, column=0, padx=5, pady=5)

        previewButton = Button(self.main, text='Preview', bg='gray', command=self.ViewGromacsTopol, font=button_font)
        previewButton.grid(row=3, column=2, padx=5, pady=5)

        saveFileButton = Button(self.main, text='Save in...', bg='gray', command=self.writeTopolFile, font=button_font)
        saveFileButton.grid(row=3, column=3, padx=5, pady=5)

    def refresh(self):
        if self.r_var.get():
            self.dimen0.configure(text='kCal/mol')
            self.dimen1.configure(text='kCal/mol')
            self.answer0.configure(text='{:>.3f}'.format(self.dG_off_kCal))
            self.answer1.configure(text='{:>.3f}'.format(self.dG_on_kCal))
        else:
            self.dimen0.configure(text='kJ/mol')
            self.dimen1.configure(text='kJ/mol')
            self.answer0.configure(text='{:>.3f}'.format(self.dG_off_kJ))
            self.answer1.configure(text='{:>.3f}'.format(self.dG_on_kJ))
        self.dimen0.update()
        self.dimen0.update()
        self.answer0.update()
        self.answer1.update()

    def ViewGromacsTopol(self):
        top = Toplevel(self.main)
        restraints = self.createGronacsRestr()
        tx = Text(top, width=130, height=20)
        tx.grid(row=0, column=0, pady=5, padx=5)
        tx.insert(END, restraints)
        tx.configure(state='disabled')
        closeTopolPrevB = Button(top, text='Exit', bg='red', command=top.destroy, font=button_font)
        closeTopolPrevB.grid(row=1, column=0, pady=5)

    def writeTopolFile(self):
        topolFile = askopenfilename(initialdir="/", title="Select file",
                                    filetypes=(("Topology files", "*.top"), ("all files", "*.*")))
        if topolFile is None:
            showerror("Error", "Topology file is not selected!")
            return
        restraints = self.createGronacsRestr()
        try:
            with open(topolFile, 'at') as f:
                f.write("\n\n" + restraints)
        except IOError:
            showerror("Error", "Topology file {:s} is not accessible for writing!".format(topolFile))

    def createGronacsRestr(self):
        restraints = ("[ intermolecular_interactions ]\n"
                      "[ bonds ]\n"
                      "; ai     aj    type   bA      kA     bB      kB\n"
                      " {12:5d}  {13:5d}  6    {0:.3f}   0.0  {0:.3f} {6:.1f} ; {20:s} - {21:s}\n"
                      " \n"
                      "[ angles ]\n"
                      "; ai     aj    ak     type    thA      fcA        thB      fcB\n"
                      " {14:5d}  {12:5d}  {13:5d}   1     {1:>6.2f}    0.0     {1:>6.2f}   {7:.2f} ; {19:s} - {20:s} - {21:s}\n"
                      " {12:5d}  {13:5d}  {15:5d}   1     {2:>6.2f}    0.0     {2:>6.2f}   {8:.2f} ; {20:s} - {21:s} - {22:s}\n"
                      "\n"
                      "[ dihedrals ]\n"
                      "; ai     aj    ak    al    type     thA      fcA       thB      fcB\n"
                      " {16:5d}  {14:5d}  {12:5d}  {13:5d}   2  {3:>7.2f}   0.0   {3:>7.2f}     {9:.2f} ; {18:s} - {19:s} - {20:s} - {21:s}\n"
                      " {14:5d}  {12:5d}  {13:5d}  {15:5d}   2  {4:>7.2f}   0.0   {4:>7.2f}   {10:>7.2f} ; {19:s} - {20:s} - {21:s} - {22:s}\n"
                      " {12:5d}  {13:5d}  {15:5d}  {17:5d}   2  {5:>7.2f}   0.0   {5:>7.2f}   {11:>7.2f} ; {20:s} - {21:s} - {22:s} - {23:s}\n"
                      "; T = {24:.1f} K\n"
                      "; dG_off = {25:.3f} kJ/mol ({27:.3f} kCal/mol)\n"
                      "; dG_on = {26:.3f} kJ/mol ({28:.3f} kCal/mol)\n"
                      ).format(
            self.bondForceParams['r_aA'],  # 0
            self.bondForceParams['th_a'],  # 1
            self.bondForceParams['th_A'],  # 2
            self.bondForceParams['phi_ba'],  # 3
            self.bondForceParams['phi_aA'],  # 4
            self.bondForceParams['phi_AB'],  # 5
            self.bondForceParams['K_r_aA'],  # 6
            self.bondForceParams['K_th_a'],  # 7
            self.bondForceParams['K_th_A'],  # 8
            self.bondForceParams['K_phi_ba'],  # 9
            self.bondForceParams['K_phi_aA'],  # 10
            self.bondForceParams['K_phi_AB'],  # 11
            self.bondForceParams['index_a'],  # 12
            self.bondForceParams['index_A'],  # 13
            self.bondForceParams['index_b'],  # 14
            self.bondForceParams['index_B'],  # 15
            self.bondForceParams['index_c'],  # 16
            self.bondForceParams['index_C'],  # 17
            self.atoms_def['index_c'],  # 18
            self.atoms_def['index_b'],  # 19
            self.atoms_def['index_a'],  # 20
            self.atoms_def['index_A'],  # 21
            self.atoms_def['index_B'],  # 22
            self.atoms_def['index_C'],  # 23
            self.bondForceParams['T'],  # 24
            self.dG_off_kJ,  # 25
            self.dG_on_kJ,  # 26
            self.dG_off_kCal,  # 27
            self.dG_on_kCal  # 28
        )
        return restraints
