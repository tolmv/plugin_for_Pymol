# coding=utf-8
from __future__ import absolute_import
from __future__ import print_function

import math
import tempfile
import tkFileDialog
import tkMessageBox
import webbrowser
from Tkinter import BooleanVar, Radiobutton, Label, Button, Tk

R = 8.314472 * 0.001  # Gas constant in kJ/mol/K
V = 1.66  # standard volume in nm^3

help_2 = """<html>
<title>Help</title>
<body>Help</body>
</html>
"""

bondForceParams = {'T': 300.0,
                   'K_r': 4184.0, 'K_thA': 41.84, 'K_thB': 41.84,
                   'K_phiA': 41.84, 'K_phiB': 41.84, 'K_phiC': 41.84,
                   'r0': 0.50, 'thA': 69.7, 'thB': 48.1,
                   'phiA': 132.2, 'phiB': 123.2, 'phiC': -12.3,
                   'index_a': 1, 'index_b': 2, 'index_c': 3,
                   'index_A': 4, 'index_B': 5, 'index_C': 6}


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
    def __init__(self, main):
        self.dG_off_kJ = calc_dG(
            bondForceParams['T'],
            bondForceParams['r0'],
            bondForceParams['thA'],
            bondForceParams['thB'],
            bondForceParams['K_r'],
            bondForceParams['K_thA'],
            bondForceParams['K_thB'],
            bondForceParams['K_phiA'],
            bondForceParams['K_phiB'],
            bondForceParams['K_phiC']
        )
        self.dG_on_kJ = -self.dG_off_kJ
        self.dG_off_kCal = kJ_to_kCal(self.dG_off_kJ)
        self.dG_on_kCal = kJ_to_kCal(self.dG_on_kJ)
        self.topolFile = None
        self.main = main
        self.r_var = BooleanVar()
        self.r_var.set(0)
        self.rj1 = Radiobutton(text='kJ', variable=self.r_var, value=0, command=self.refresh)
        self.rcal1 = Radiobutton(text="kCal", variable=self.r_var, value=1, command=self.refresh)
        self.rj1.grid(row=0, column=0)
        self.rcal1.grid(row=0, column=1)

        self.name0 = Label(main, text='dG_off = ', font=15)
        self.name1 = Label(main, text='dG_on = ', font=15)
        self.name0.grid(row=1, column=0)
        self.name1.grid(row=2, column=0)

        self.answer0 = Label(main, font=15)
        self.answer1 = Label(main, font=15)
        self.answer0['text'] = '{:>.3f}'.format(self.dG_off_kJ)
        self.answer1['text'] = '{:>.3f}'.format(self.dG_on_kJ)
        self.answer0.grid(row=1, column=1)
        self.answer1.grid(row=2, column=1)

        self.dimen0 = Label(main, font=15)
        self.dimen1 = Label(main, font=15)
        self.dimen0['text'] = 'kJ'
        self.dimen1['text'] = 'kJ'
        self.dimen0.grid(row=1, column=2)
        self.dimen1.grid(row=2, column=2)

        self.destroyProgr = Button(main, text='Exit', bg='red', command=main.destroy)
        self.destroyProgr.grid(row=0, column=3)

        self.helpProgr = Button(main, text=' ? ', bg='#ffb3fe', command=self.getHelp)
        self.helpProgr.grid(row=4, column=0)

        self.name0 = Label(main, text='Gromacs topology file:', font=15)
        self.name0.grid(row=3, column=0)

        self.openFileButton = Button(main, text='Select', bg='gray', command=self.selectFile)
        self.openFileButton.grid(row=3, column=2)

        self.openFileButton = Button(main, text='Write', bg='green', command=self.writeTopolFile)
        self.openFileButton.grid(row=3, column=3)

    def refresh(self):

        if self.r_var.get():
            self.dimen0.configure(text='kCal')
            self.dimen1.configure(text='kCal')
            self.answer0.configure(text='{:>.3f}'.format(self.dG_off_kCal))
            self.answer1.configure(text='{:>.3f}'.format(self.dG_on_kCal))
        else:
            self.dimen0.configure(text='kJ')
            self.dimen1.configure(text='kJ')
            self.answer0.configure(text='{:>.3f}'.format(self.dG_off_kJ))
            self.answer1.configure(text='{:>.3f}'.format(self.dG_on_kJ))

        self.dimen0.update()
        self.dimen0.update()
        self.answer0.update()
        self.answer1.update()

    def getHelp(self):
        with tempfile.NamedTemporaryFile('w', delete=False, suffix='.html') as f:
            url = "file://" + f.name
            f.write(help_2)
        webbrowser.open(url)

    def selectFile(self):
        self.topolFile = tkFileDialog.askopenfilename(initialdir="/", title="Select file",
                                                      filetypes=(("Topology files", "*.top"), ("all files", "*.*")))

    def writeTopolFile(self):
        if self.topolFile is None:
            tkMessageBox.showerror("Error", "Topology file is not selected")
            return
        restraints = ("[ intermolecular_interactions ]\n"
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
            bondForceParams['r0'],  # 0
            bondForceParams['thA'],  # 1
            bondForceParams['thB'],  # 2
            bondForceParams['phiA'],  # 3
            bondForceParams['phiB'],  # 4
            bondForceParams['phiC'],  # 5
            bondForceParams['K_r'],  # 6
            bondForceParams['K_thA'],  # 7
            bondForceParams['K_thB'],  # 8
            bondForceParams['K_phiA'],  # 9
            bondForceParams['K_phiB'],  # 10
            bondForceParams['K_phiC'],  # 11
            bondForceParams['index_a'],  # 12
            bondForceParams['index_A'],  # 13
            bondForceParams['index_b'],  # 14
            bondForceParams['index_B'],  # 15
            bondForceParams['index_c'],  # 16
            bondForceParams['index_C'])  # 17
        print(restraints)
        with open(self.topolFile, 'at') as f:
            f.write(restraints)


def main():
    root = Tk()
    app = Output(root)
    root.mainloop()


if __name__ == '__main__':
    main()
