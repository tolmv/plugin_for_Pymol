import sys
import math
from Tkinter import *
import tkMessageBox as mb
import numpy as np


class Restraints():

    def __init__(self, main):
        self.help = []
        self.bol = True
        self.error = 0
        self.dG = 0
        self.main = main
        self.K = 8.314472*0.001     # Gas constant in kJ/mol/K
        self.V = 1.66               # standard volume in nm^3

        self.labels = ['K_r', 'K_thA',
                       'K_thB', 'K_phiA', 'K_phiB', 'K_phiC']

        self.r_var = BooleanVar()
        self.r_var.set(1)
        self.rcal1 = Radiobutton(text="kCal", variable=self.r_var, value=0)
        self.rj1 = Radiobutton(text='kJ', variable=self.r_var, value=1)

        self.rj1.grid(row=0, column=0)
        self.rcal1.grid(row=0, column=1)


        self.entry0 = Entry(main)
        self.entry1 = Entry(main)
        self.entry2 = Entry(main)
        self.entry3 = Entry(main)
        self.entry4 = Entry(main)
        self.entry5 = Entry(main)


        self.label_answer0 = Label(main, font=15)
        self.label_answer1 = Label(main, font=15)
        self.label_answer2 = Label(main, font=15)
        self.label_answer3 = Label(main, font=15)
        self.label_answer4 = Label(main, font=15)
        self.label_answer5 = Label(main, font=15)


        self.dimen0 = Label(main, font=15)
        self.dimen1 = Label(main, font=15)
        self.dimen2 = Label(main, font=15)
        self.dimen3 = Label(main, font=15)
        self.dimen4 = Label(main, font=15)
        self.dimen5 = Label(main, font=15)


        self.dimen_all = [self.dimen0, self.dimen1, self.dimen2, self.dimen3, self.dimen4,
                          self.dimen5]

        self.entry_all = [self.entry0, self.entry1, self.entry2, self.entry3, self.entry4,
                          self.entry5]

        self.label_all = [self.label_answer0, self.label_answer1,
                          self.label_answer2, self.label_answer3, self.label_answer4,
                          self.label_answer5]

        self.button_res = Button(main, text="Next -> ")

        for i in range(10):
            self.label_all[i].grid(row=i + 1, column=1)
            self.entry_all[i].grid(row=i + 1, column=2)

        for i in range(len(self.dimen_all)):
            self.dimen_all[i].grid(row=i + 1, column=3)

        self.dimen_all[0]['text'] = 'kJ/mol/rad^2'
        self.dimen_all[1]['text'] = 'kJ/mol/rad^2'
        self.dimen_all[2]['text'] = 'kJ/mol/rad^2'
        self.dimen_all[3]['text'] = 'kJ/mol/rad^2'
        self.dimen_all[4]['text'] = 'kJ/mol/rad^2'
        self.dimen_all[5]['text'] = 'kJ/mol/rad^2'

        if self.r_var == 0:
            self.dimen_all[0]['text'] = 'kCal/mol/nm^2'
            self.dimen_all[1]['text'] = 'kCal/mol/rad^2'
            self.dimen_all[2]['text'] = 'kCal/mol/rad^2'
            self.dimen_all[3]['text'] = 'kCal/mol/rad^2'
            self.dimen_all[4]['text'] = 'kCal/mol/rad^2'
            self.dimen_all[5]['text'] = 'kCal/mol/rad^2'



        for i in range(len(self.label_all)):
            self.label_all[i]['text'] = self.labels[i]


        self.entry_all_get = [self.entry0.get, self.entry1.get, self.entry2.get,
                              self.entry3.get, self.entry4.get,
                              self.entry5.get]

        self.button_res.grid(row=11, column=2)

        self.button_res.bind('<Button-1>', self.__calculate)

        self.destroyProgr = Button(main, text='Exit', bg='black', command=main.destroy)
        self.destroyProgr.grid(row=0, column=3)

        self.helpProgr = Button(main, text=' ? ', bg='#ffb3fe')
        self.helpProgr.grid(row=12, column=0)

    def __calculate(self, event):

        for i, k in enumerate(self.entry_all_get):
            try:
                f = float(k())
                if f <= 0:
                    raise ValueError
            except ValueError:
                self.entry_all[i]['bg'] = "red"
                return
            self.help.append(f)

        self.K_r = self.help[0]     # force constant for distance (kJ/mol/nm^2)
        self.K_thA = self.help[1]   # force constant for angle (kJ/mol/rad^2)
        self.K_thB = self.help[2]   # force constant for angle (kJ/mol/rad^2)
        self.K_phiA = self.help[3]  # force constant for dihedral (kJ/mol/rad^2)
        self.K_phiB = self.help[4]  # force constant for dihedral (kJ/mol/rad^2)
        self.K_phiC = self.help[5]  # force constant for dihedral (kJ/mol/rad^2)
        if self.r_var == 0:

            self.K_r = self.help[0] * 0.238846     # force constant for distance (kCal/mol/nm^2)
            self.K_thA = self.help[1] * 0.238846   # force constant for angle (kCal/mol/rad^2)
            self.K_thB = self.help[2] * 0.238846   # force constant for angle (kCal/mol/rad^2)
            self.K_phiA = self.help[3] * 0.238846  # force constant for dihedral (kCal/mol/rad^2)
            self.K_phiB = self.help[4] * 0.238846  # force constant for dihedral (kCal/mol/rad^2)
            self.K_phiC = self.help[5] * 0.238846  # force constant for dihedral (kCal/mol/rad^2)


        self.rt = App(self.main)


class App():
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


if __name__ == '__main__':
    main()
