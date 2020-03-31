# coding=utf-8
from __future__ import absolute_import
from __future__ import print_function

import tempfile
import webbrowser
from Tkinter import BooleanVar, Radiobutton, Entry, Label, Button, Tk, Toplevel

help_1 = """<html>
<title>Help</title>
<body>Help</body>
</html>
"""


def kCal_to_kJ(E):
    return E * 4.1868


bondForceParams = {'T': None,
                   'K_r': None, 'K_thA': None, 'K_thB': None,
                   'K_phiA': None, 'K_phiB': None, 'K_phiC': None,
                   'r0': None, 'thA': None, 'thB': None,
                   'phiA': None, 'phiB': None, 'phiC': None,
                   'index_a': None, 'index_b': None, 'index_c': None,
                   'index_A': None, 'index_B': None, 'index_C': None}


class Restraints(object):
    def __init__(self, main):
        self.help = []
        self.main = main
        self.labels = ['Temp', 'K_r', 'K_thA',
                       'K_thB', 'K_phiA', 'K_phiB', 'K_phiC']
        self.r_var = BooleanVar()
        self.r_var.set(0)
        self.rj1 = Radiobutton(text='kJ', variable=self.r_var, value=0, command=self.refresh)
        self.rcal1 = Radiobutton(text="kCal", variable=self.r_var, value=1, command=self.refresh)
        self.rj1.grid(row=0, column=0)
        self.rcal1.grid(row=0, column=1)

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

        for i in range(7):
            self.label_all[i].grid(row=i + 1, column=1)
            self.entry_all[i].grid(row=i + 1, column=2)

        for i in range(len(self.dimen_all)):
            self.dimen_all[i].grid(row=i + 1, column=3)

        self.dimen_all[0]['text'] = 'Kelvin'
        self.dimen_all[1]['text'] = 'kJ/mol/nm^2'
        self.dimen_all[2]['text'] = 'kJ/mol/rad^2'
        self.dimen_all[3]['text'] = 'kJ/mol/rad^2'
        self.dimen_all[4]['text'] = 'kJ/mol/rad^2'
        self.dimen_all[5]['text'] = 'kJ/mol/rad^2'
        self.dimen_all[6]['text'] = 'kJ/mol/rad^2'

        for i in range(len(self.label_all)):
            self.label_all[i]['text'] = self.labels[i]

        self.entry_all_get = [self.entry0.get, self.entry1.get, self.entry2.get,
                              self.entry3.get, self.entry4.get,
                              self.entry5.get, self.entry6.get]

        self.button_res = Button(main, text="Next -> ", command=self.validate)
        self.button_res.grid(row=11, column=2)

        self.destroyProgr = Button(main, text='Exit', bg='red', command=main.destroy)
        self.destroyProgr.grid(row=0, column=3)

        self.helpProgr = Button(main, text=' ? ', bg='#ffb3fe', command=self.getHelp)
        self.helpProgr.grid(row=12, column=0)

    def refresh(self):

        if self.r_var.get():
            self.dimen_all[1].configure(text='kCal/mol/nm^2')
            self.dimen_all[2].configure(text='kCal/mol/rad^2')
            self.dimen_all[3].configure(text='kCal/mol/rad^2')
            self.dimen_all[4].configure(text='kCal/mol/rad^2')
            self.dimen_all[5].configure(text='kCal/mol/rad^2')
            self.dimen_all[6].configure(text='kCal/mol/rad^2')
        else:
            self.dimen_all[1].configure(text='kJ/mol/nm^2')
            self.dimen_all[2].configure(text='kJ/mol/rad^2')
            self.dimen_all[3].configure(text='kJ/mol/rad^2')
            self.dimen_all[4].configure(text='kJ/mol/rad^2')
            self.dimen_all[5].configure(text='kJ/mol/rad^2')
            self.dimen_all[6].configure(text='kJ/mol/rad^2')

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
            self.help = list((self.help[1],)) + list(map(kCal_to_kJ, self.help[1:]))

        bondForceParams['T'] = self.help[0]  # Temperature (K)
        bondForceParams['K_r'] = self.help[1]  # force constant for distance (kJ/mol/nm^2)
        bondForceParams['K_thA'] = self.help[2]  # force constant for angle (kJ/mol/rad^2)
        bondForceParams['K_thB'] = self.help[3]  # force constant for angle (kJ/mol/rad^2)
        bondForceParams['K_phiA'] = self.help[4]  # force constant for dihedral (kJ/mol/rad^2)
        bondForceParams['K_phiB'] = self.help[5]  # force constant for dihedral (kJ/mol/rad^2)
        bondForceParams['K_phiC'] = self.help[6]  # force constant for dihedral (kJ/mol/rad^2)

        self.rt = App(self.main)

    def getHelp(self):
        with tempfile.NamedTemporaryFile('w', delete=False, suffix='.html') as f:
            url = "file://" + f.name
            f.write(help_1)
        webbrowser.open(url)


class App(object):
    def __init__(self, main):
        self.res_top = Toplevel(main)
        self.now_do = Label(self.res_top, font=15)
        self.now_do['text'] = 'Now choose the atoms you need'
        self.now_do.config(bd=20, bg='#aaffff')
        self.now_do.pack()


def main():
    root = Tk()
    app = Restraints(root)
    root.mainloop()
    print(bondForceParams)


if __name__ == '__main__':
    main()
