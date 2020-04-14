# coding=utf-8
from __future__ import absolute_import

import tempfile
import webbrowser
from pymol import cmd

try:
    from Tkinter import BooleanVar, Radiobutton, Entry, Label, Button, Toplevel, W
    from tkMessageBox import showinfo
except ImportError:
    from tkinter import BooleanVar, Radiobutton, Entry, Label, Button, Toplevel, W
    from tkinter.messagebox import showinfo

from .wizard import RestraintWizard

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
 </body>
</html>

"""


def kCal_to_kJ(E):
    return E * 4.1868


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
        rj1 = Radiobutton(self.main, text='kJ', variable=self.r_var, value=0, command=self.refresh)
        rcal1 = Radiobutton(self.main, text="kCal", variable=self.r_var, value=1, command=self.refresh)
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
            label_answer = Label(self.main, text=lab, font=15, anchor=W)
            label_all.append(label_answer)
            entry = Entry(self.main)
            self.entry_all.append(entry)
            self.entry_all_get.append(entry.get)
            dimen = Label(self.main, font=15, anchor=W)
            self.dimen_all.append(dimen)

        for i, (label, entry, dimen) in enumerate(zip(label_all, self.entry_all, self.dimen_all)):
            label.grid(row=i + 1, column=1, padx=5, pady=5, sticky=W)
            entry.grid(row=i + 1, column=2, padx=5, pady=5, sticky=W)
            dimen.grid(row=i + 1, column=3, padx=10, pady=5, sticky=W)

        self.dimen_all[0]['text'] = 'Kelvin'
        self.refresh()

        self.button_res = Button(self.main, text="Next -> ", command=self.validate)
        self.button_res.grid(row=11, column=2, padx=5, pady=5)

        self.destroyProgr = Button(self.main, text='Exit', bg='red', command=self.main.destroy)
        self.destroyProgr.grid(row=0, column=3, padx=5, pady=5)

        self.helpProgr = Button(self.main, text=' ? ', bg='#ffb3fe', command=self.getHelp)
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

    @staticmethod
    def getHelp():
        with tempfile.NamedTemporaryFile('w', delete=False, suffix='.html') as f:
            url = "file://" + f.name
            f.write(help_1)
        webbrowser.open(url)
