import math
import webbrowser
import tempfile
import tkFileDialog
import tkMessageBox
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
                   'r0':  0.50, 'thA': 69.7, 'thB': 48.1,
                    'phiA' : None, 'phiA' : None, 'phiA' : None,
                   'index_a': None, 'index_b': None, 'index_c': None,
                   'index_A': None, 'index_B': None, 'index_C': None}

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

class Output():

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
            self.answer0.configure(text = '{:>.3f}'.format(self.dG_off_kCal))
            self.answer1.configure(text = '{:>.3f}'.format(self.dG_on_kCal))
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
        restraints = """[ intermolecular_interactions]
[ bonds ]
; ai     aj    type   bA      kA     bB      kB
 6573    1629  6      0.5     0.0    0.5   4184.0
 
[ angles ]
; ai     aj    ak     type    thA      fcA        thB      fcB        
 6573   1629   1626   1       69.7     0.0        69.7    41.84
 6574   6573   1629   1       48.1     0.0        48.1    41.84 ; 

[ dihedrals ]
; ai     aj    ak    al    type     thA      fcA       thB      fcB        
 6573    1629  1626  1623   2       149.2    0.0       149.2    41.84
 6574    6573  1629  1626   2        76.2    0.0        76.2    41.84
 6576    6574  6573  1629   2        71.1    0.0        71.1    41.84
        """.format()



def main():
    root = Tk()
    app = Output(root)
    root.mainloop()


if __name__ == '__main__':
    main()

