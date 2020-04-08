# coding=utf-8
from __future__ import absolute_import
from __future__ import print_function

import math
import tempfile
import webbrowser
from Tkinter import BooleanVar, Radiobutton, Label, Button
from tkFileDialog import askopenfilename
from tkMessageBox import showerror

R = 8.314472 * 0.001  # Gas constant in kJ/mol/K
V = 1.66  # standard volume in nm^3

help_2 = """<html>
<title>Help</title>
<h1 align="left" style="color: Black">QuickStart</h1>

<h2 align="left" style="color: Black"> Gromacs topology</h2>

<h3 align="left" style="color: Black"> Select</h4>

<body>
Select the topology file that you need.
</body>

<h3 align="left" style="color: Black"> Write</h4>

<body>
This button writes all the indices of the atoms, <br>
force constants, distances, angles and dihedral angles in your topology file,<br>
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
                    ((K_r * K_thA * K_thB * K_phiA * K_phiB * K_phiC) ** 0.5) / ((2.0 * math.pi * R * T) ** (3.0))
            )
    )

    dG = - R * T * math.log(arg)
    return dG


class Output(object):
    def __init__(self, main, bondForceParams):
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
        self.topolFile = None
        self.main = main
        self.main.title('PyFepRestr')
        self.r_var = BooleanVar()
        self.r_var.set(0)
        self.rj1 = Radiobutton(main, text='kJ', variable=self.r_var, value=0, command=self.refresh)
        self.rcal1 = Radiobutton(main, text="kCal", variable=self.r_var, value=1, command=self.refresh)
        self.rj1.grid(row=0, column=0, padx=5, pady=5)
        self.rcal1.grid(row=0, column=1, padx=5, pady=5)

        self.name0 = Label(main, text=u'\u0394G_off = ', font=15)
        self.name1 = Label(main, text=u'\u0394G_on = ', font=15)
        self.name0.grid(row=1, column=0, padx=5, pady=5)
        self.name1.grid(row=2, column=0, padx=5, pady=5)

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

        self.destroyProgr = Button(main, text='Exit', bg='red', command=main.destroy)
        self.destroyProgr.grid(row=0, column=3, padx=5, pady=5)

        self.helpProgr = Button(main, text=' ? ', bg='#ffb3fe', command=self.getHelp)
        self.helpProgr.grid(row=4, column=0, padx=5, pady=5)

        self.name0 = Label(main, text='Gromacs topology file:', font=15)
        self.name0.grid(row=3, column=0, padx=5, pady=5)

        self.openFileButton = Button(main, text='Select', bg='gray', command=self.selectFile)
        self.openFileButton.grid(row=3, column=2, padx=5, pady=5)

        self.openFileButton = Button(main, text='Write', bg='green', command=self.writeTopolFile)
        self.openFileButton.grid(row=3, column=3, padx=5, pady=5)

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

    def selectFile(self):
        self.topolFile = askopenfilename(initialdir="/", title="Select file",
                                         filetypes=(("Topology files", "*.top"), ("all files", "*.*")))

    def writeTopolFile(self):
        if self.topolFile is None:
            showerror("Error", "Topology file is not selected")
            return
        restraints = ("\n\n[ intermolecular_interactions ]\n"
                      "[ bonds ]\n"
                      "; ai     aj    type   bA      kA     bB      kB\n"
                      " {12:d}    {13:d}  6      {0:.3f}     0.0    {0:.3f}   {6:.1f}\n"
                      " \n"
                      "[ angles ]\n"
                      "; ai     aj    ak     type    thA      fcA        thB      fcB\n"
                      " {14:d}   {12:d}   {13:d}   1       {1:.2f}     0.0        {1:.2f}    {7:.2f}\n"
                      " {12:d}   {13:d}   {15:d}   1       {2:.2f}     0.0        {2:.2f}    {8:.2f}\n"
                      "\n"
                      "[ dihedrals ]\n"
                      "; ai     aj    ak    al    type     thA      fcA       thB      fcB\n"
                      " {16:d}   {14:d}   {12:d}   {13:d}   2       149.2    0.0       149.2    {9:.2f}\n"
                      " {14:d}   {12:d}   {13:d}   {15:d}   2        76.2    0.0        76.2    {10:.2f}\n"
                      " {12:d}   {13:d}   {15:d}   {17:d}   2        71.1    0.0        71.1    {11:.2f}\n"
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
            self.bondForceParams['index_C'])  # 17
        print(restraints)
        with open(self.topolFile, 'at') as f:
            f.write(restraints)
