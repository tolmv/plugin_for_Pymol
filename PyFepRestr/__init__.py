# coding=utf-8
from __future__ import absolute_import
from __future__ import print_function

from Tkinter import Toplevel
from pymol import cmd

from .output import Output
from .tkinput import Restraints
from .wizard import RestraintWizard


def __init_plugin__(app):
    app.menuBar.addmenuitem('Plugin', 'command',
                            label='PyFepRestr',
                            command=lambda: main(app.root))


def main(parrent):
    top = Toplevel(parrent)
    bondForceParams = {'T': None,
                       'K_r': None, 'K_thA': None, 'K_thB': None,
                       'K_phiA': None, 'K_phiB': None, 'K_phiC': None,
                       'r0': None, 'thA': None, 'thB': None,
                       'phiA': None, 'phiB': None, 'phiC': None,
                       'index_a': None, 'index_b': None, 'index_c': None,
                       'index_A': None, 'index_B': None, 'index_C': None}
    restr = Restraints(top, bondForceParams)
    top.transient(parrent)
    top.grab_set()
    top.focus_set()
    top.wait_window()
    if restr.isexit:
        del top
        return
    del top
    wiz = RestraintWizard(bondForceParams)
    cmd.set_wizard(wiz)
    top = Toplevel(parrent)
    Output(top, bondForceParams)
    top.transient(parrent)
    top.grab_set()
    top.focus_set()
    top.wait_window()
    del top


if __name__ == '__main__':
    from Tkinter import Tk

    root = Tk()
    main(root)
    root.mainloop()
