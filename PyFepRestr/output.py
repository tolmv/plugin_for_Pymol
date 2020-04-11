# coding=utf-8
from __future__ import absolute_import
from __future__ import print_function

import math
import tempfile
import webbrowser

try:
    from Tkinter import BooleanVar, Radiobutton, Label, Button, Frame, Text, Scrollbar, W, END, WORD, Y, LEFT, RIGHT, \
        Toplevel
    from tkFileDialog import askopenfilename
    from tkMessageBox import showerror
except ImportError:
    from tkinter import BooleanVar, Radiobutton, Label, Button, Frame, Text, Scrollbar, W, END, WORD, Y, LEFT, RIGHT, \
        Toplevel
    from tkinter.filedialog import askopenfilename
    from tkinter.messagebox import askyesno, showerror

R = 8.314472 * 0.001  # Gas constant in kJ/mol/K
V = 1.66  # standard volume in nm^3

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


def kJ_to_kCal(E):
    return E / 4.1868


def calc_dG(T, r0, thA, thB, K_r, K_thA, K_thB, K_phiA, K_phiB, K_phiC):
    """BORESCH FORMULA - Calculate dG restraints off"""
    thA = math.radians(thA)  # convert angle from degrees to radians --> math.sin() wants radians
    thB = math.radians(thB)  # convert angle from degrees to radians --> math.sin() wants radians

    arg = (
            (8.0 * math.pi ** 2.0 * V) / (r0 ** 2.0 * math.sin(thA) * math.sin(thB))
            *
            (
                    ((K_r * K_thA * K_thB * K_phiA * K_phiB * K_phiC) ** 0.5) / ((2.0 * math.pi * R * T) ** 3.0)
            )
    )

    dG = - R * T * math.log(arg)
    return dG


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
        self.topolFile = None
        self.main = main
        self.main.title('PyFepRestr')
        self.r_var = BooleanVar()
        self.r_var.set(0)
        rj1 = Radiobutton(main, text='kJ', variable=self.r_var, value=0, command=self.refresh)
        rcal1 = Radiobutton(main, text="kCal", variable=self.r_var, value=1, command=self.refresh)
        rj1.grid(row=0, column=0, padx=5, pady=5)
        rcal1.grid(row=0, column=1, padx=5, pady=5)

        name0 = Label(main, text=u'\u0394G_off = ', font=15)
        name1 = Label(main, text=u'\u0394G_on = ', font=15)
        name0.grid(row=1, column=0, padx=5, pady=5)
        name1.grid(row=2, column=0, padx=5, pady=5)

        self.answer0 = Label(main, font=15)
        self.answer1 = Label(main, font=15)
        self.answer0['text'] = '{:>.3f}'.format(self.dG_off_kJ)
        self.answer1['text'] = '{:>.3f}'.format(self.dG_on_kJ)
        self.answer0.grid(row=1, column=1, padx=5, pady=5)
        self.answer1.grid(row=2, column=1, padx=5, pady=5)

        self.dimen0 = Label(main, font=15)
        self.dimen1 = Label(main, font=15)
        self.refresh()
        self.dimen0.grid(row=1, column=2, padx=5, pady=5)
        self.dimen1.grid(row=2, column=2, padx=5, pady=5)

        self.tx = None

        destroyProgr = Button(main, text='Exit', bg='red', command=main.destroy)
        destroyProgr.grid(row=0, column=3, padx=5, pady=5)

        helpProgr = Button(main, text=' ? ', bg='#ffb3fe', command=self.getHelp)
        helpProgr.grid(row=4, column=0, padx=5, pady=5)

        name0 = Label(main, text='Gromacs topology file:', font=15)
        name0.grid(row=3, column=0, padx=5, pady=5)

        previewButton = Button(main, text='Preview', bg='gray', command=self.ViewGromacsTopol)
        previewButton.grid(row=3, column=2, padx=5, pady=5)

        saveFileButton = Button(main, text='Save', bg='gray', command=self.writeTopolFile)
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

    @staticmethod
    def getHelp():
        with tempfile.NamedTemporaryFile('w', delete=False, suffix='.html') as f:
            url = "file://" + f.name
            f.write(help_2)
        webbrowser.open(url)

    def ViewGromacsTopol(self):
        top = Toplevel(self.main)
        restraints = self.createGronacsRestr()
        self.tx = Text(top, width=130, height=20, wrap=WORD)
        self.tx.grid(row=0, column=0, pady=5, padx=5)
        self.tx.insert(END, restraints)
        self.tx.configure(state='disabled')
        closeTopolPrevB = Button(top, text='Exit', bg='red', command=top.destroy)
        closeTopolPrevB.grid(row=1, column=0, pady=5)

    def writeTopolFile(self):
        topolFile = askopenfilename(initialdir="/", title="Select file",
                                    filetypes=(("Topology files", "*.top"), ("all files", "*.*")))
        if topolFile is None:
            showerror("Error", "Topology file is not selected!")
            return
        restraints = self.createGronacsRestr()
        with open(topolFile, 'at') as f:
            f.write("\n\n" + restraints)

    def createGronacsRestr(self):
        restraints = ("[ intermolecular_interactions ]\n"
                      "[ bonds ]\n"
                      "; ai     aj    type   bA      kA     bB      kB\n"
                      " {12:d}    {13:d}  6      {0:.3f}     0.0    {0:.3f}   {6:.1f} ; {20:s} - {21:s}\n"
                      " \n"
                      "[ angles ]\n"
                      "; ai     aj    ak     type    thA      fcA        thB      fcB\n"
                      " {14:d}   {12:d}   {13:d}   1       {1:.2f}     0.0        {1:.2f}    {7:.2f} ; {19:s} - {20:s} - {21:s}\n"
                      " {12:d}   {13:d}   {15:d}   1       {2:.2f}     0.0        {2:.2f}    {8:.2f} ; {20:s} - {21:s} - {22:s}\n"
                      "\n"
                      "[ dihedrals ]\n"
                      "; ai     aj    ak    al    type     thA      fcA       thB      fcB\n"
                      " {16:d}   {14:d}   {12:d}   {13:d}   2     {3:.2f}    0.0     {3:.2f}     {9:.2f} ; {18:s} - {19:s} - {20:s} - {21:s}\n"
                      " {14:d}   {12:d}   {13:d}   {15:d}   2     {4:.2f}    0.0     {4:.2f}    {10:.2f} ; {19:s} - {20:s} - {21:s} - {22:s}\n"
                      " {12:d}   {13:d}   {15:d}   {17:d}   2     {5:.2f}    0.0     {5:.2f}    {11:.2f} ; {20:s} - {21:s} - {22:s} - {23:s}\n"
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
