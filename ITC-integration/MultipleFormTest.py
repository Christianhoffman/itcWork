from Tkinter import *
from PIL import Image, ImageTk

master = Tk()

def f1():
    print("Change to form A")

def f2():
    print("Change to form B")

title = Label(text="Swithing Menus or Forms")
title.grid(row=0, column=1)

blank = Label(text="\n")
blank.grid(row=1, column=1)


b = Button(master, text="Submenu A", command=f1)
b.grid(row=2, column=1)

b2 = Button(master, text="Submenu B", command=f2)
b2.grid(row=2, column=2)
mainloop()