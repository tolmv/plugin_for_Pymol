# coding=utf-8
from __future__ import absolute_import
from __future__ import print_function

import tempfile
import webbrowser
from pymol import cmd

try:
    from Tkinter import BooleanVar, Radiobutton, Entry, Label, Button
    from tkMessageBox import showinfo
except ImportError:
    from tkinter import BooleanVar, Radiobutton, Entry, Label, Button
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
    def __init__(self, parent, main, bondForceParams):
        self.help = []
        self.parent = parent
        self.main = main
        self.main.title('PyFepRestr')
        self.main.protocol('WM_DELETE_WINDOW', self.exit)
        self.labels = ['Temp',
                       'K raA',
                       u'K \u03b8a', u'K \u03b8A',
                       u'K \u03c6ba', u'K \u03c6aA', u'K_\u03c6AB']
        self.r_var = BooleanVar()
        self.r_var.set(1)
        self.rj1 = Radiobutton(main, text='kJ', variable=self.r_var, value=0, command=self.refresh)
        self.rcal1 = Radiobutton(main, text="kCal", variable=self.r_var, value=1, command=self.refresh)
        self.rj1.grid(row=0, column=0, padx=5, pady=5)
        self.rcal1.grid(row=0, column=1, padx=5, pady=5)

        self.entry0 = Entry(main)
        self.entry1 = Entry(main)
        self.entry2 = Entry(main)
        self.entry3 = Entry(main)
        self.entry4 = Entry(main)
        self.entry5 = Entry(main)
        self.entry6 = Entry(main)

        self.label_answer0 = Label(main, font=15)
        self.label_answer1 = Label(main, font=15)
        self.label_answer2 = Label(main, font=15)
        self.label_answer3 = Label(main, font=15)
        self.label_answer4 = Label(main, font=15)
        self.label_answer5 = Label(main, font=15)
        self.label_answer6 = Label(main, font=15)

        self.dimen0 = Label(main, font=15)
        self.dimen1 = Label(main, font=15)
        self.dimen2 = Label(main, font=15)
        self.dimen3 = Label(main, font=15)
        self.dimen4 = Label(main, font=15)
        self.dimen5 = Label(main, font=15)
        self.dimen6 = Label(main, font=15)

        self.dimen_all = [self.dimen0, self.dimen1, self.dimen2, self.dimen3, self.dimen4,
                          self.dimen5, self.dimen6]

        self.entry_all = [self.entry0, self.entry1, self.entry2, self.entry3, self.entry4,
                          self.entry5, self.entry6]

        self.label_all = [self.label_answer0, self.label_answer1,
                          self.label_answer2, self.label_answer3, self.label_answer4,
                          self.label_answer5, self.label_answer6]

        for i, (label, entry, dimen) in enumerate(zip(self.label_all, self.entry_all, self.dimen_all)):
            label.grid(row=i + 1, column=1, padx=5, pady=5)
            entry.grid(row=i + 1, column=2, padx=5, pady=5)
            dimen.grid(row=i + 1, column=3, padx=5, pady=5)

        self.dimen_all[0]['text'] = 'Kelvin'
        self.refresh()

        for i in range(len(self.label_all)):
            self.label_all[i]['text'] = self.labels[i]

        self.entry_all_get = [self.entry0.get, self.entry1.get, self.entry2.get,
                              self.entry3.get, self.entry4.get,
                              self.entry5.get, self.entry6.get]
        self.bondForceParams = bondForceParams
        self.button_res = Button(main, text="Next -> ", command=self.validate)
        self.button_res.grid(row=11, column=2, padx=5, pady=5)

        self.destroyProgr = Button(main, text='Exit', bg='red', command=self.exit)
        self.destroyProgr.grid(row=0, column=3, padx=5, pady=5)

        self.helpProgr = Button(main, text=' ? ', bg='#ffb3fe', command=self.getHelp)
        self.helpProgr.grid(row=12, column=0, padx=5, pady=5)

    def exit(self):
        self.main.destroy()

    def refresh(self):

        if self.r_var.get():
            self.dimen_all[1].configure(text=u'kCal/mol/\u212b\u00b2')
            self.dimen_all[2].configure(text=u'kCal/mol/rad\u00b2')
            self.dimen_all[3].configure(text=u'kCal/mol/rad\u00b2')
            self.dimen_all[4].configure(text=u'kCal/mol/rad\u00b2')
            self.dimen_all[5].configure(text=u'kCal/mol/rad\u00b2')
            self.dimen_all[6].configure(text=u'kCal/mol/rad\u00b2')
        else:
            self.dimen_all[1].configure(text=u'kJ/mol/\u212b\u00b2')
            self.dimen_all[2].configure(text=u'kJ/mol/rad\u00b2')
            self.dimen_all[3].configure(text=u'kJ/mol/rad\u00b2')
            self.dimen_all[4].configure(text=u'kJ/mol/rad\u00b2')
            self.dimen_all[5].configure(text=u'kJ/mol/rad\u00b2')
            self.dimen_all[6].configure(text=u'kJ/mol/rad\u00b2')

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
            self.help.append(f)

        if self.r_var.get():
            self.help = list((self.help[0],)) + list(map(kCal_to_kJ, self.help[1:]))

        self.bondForceParams['T'] = self.help[0]  # Temperature (K)
        self.bondForceParams['K_r_aA'] = self.help[1] * 100  # force constant for distance (kJ/mol/nm^2)
        self.bondForceParams['K_th_a'] = self.help[2]  # force constant for angle (kJ/mol/rad^2)
        self.bondForceParams['K_th_A'] = self.help[3]  # force constant for angle (kJ/mol/rad^2)
        self.bondForceParams['K_phi_ba'] = self.help[4]  # force constant for dihedral (kJ/mol/rad^2)
        self.bondForceParams['K_phi_aA'] = self.help[5]  # force constant for dihedral (kJ/mol/rad^2)
        self.bondForceParams['K_phi_AB'] = self.help[6]  # force constant for dihedral (kJ/mol/rad^2)
        showinfo('Info', 'Now choose the atoms you need')
        atoms_def = {
            'index_c': None,
            'index_b': None,
            'index_a': None,
            'index_A': None,
            'index_B': None,
            'index_C': None
        }
        wiz = RestraintWizard(self.parent, self.bondForceParams, atoms_def)
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
