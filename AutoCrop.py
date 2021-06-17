from PIL import Image
import os
import tkinter as tk
from tkinter import ttk



def autocrop(folderPath ,searchstring,Newfolderpath,testmachine):
    directory = os.listdir(folderPath)
    for fname in directory:
        if os.path.isfile(folderPath+'/'+fname):
            f = open(folderPath + '/' + fname, 'r')
            if searchstring in fname:
                # load image
                image = Image.open((folderPath + '/' + fname))
                # create a cropped image
                if str(testmachine)=='RNS':
                    xx=7
                    xy=44
                    yx=479
                    yy=377
                if testmachine == 'test2':
                    xx =0
                    xy =0
                    yx =0
                    yy =0
                if testmachine == 'test3':
                    xx =0
                    xy =0
                    yx =0
                    yy =0
                cropped = image.crop((int(xx), int(xy), int(yx), int(yy)))
                cropped.show()
                # starts at 0:0 and then goes down to 960:540
                # topright to bottom right
                # show cropped image
                fnamenew = ' '
                for svalues in fname:
                    if svalues == '.':
                        fnamenew = fnamenew + '_cr' + svalues
                    else:
                        fnamenew = fnamenew + svalues
                cropped.save(Newfolderpath+'/'+fnamenew)
            else:
                f.close()


class button_work(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)

        self.entry1 = tk.Entry(self)
        self.entry1.grid(row=0, column=1)

        self.label1 = tk.Label(self, text="File Path").grid(row=0)

        self.entry2 = tk.Entry(self)
        self.entry2.grid(row=1,column=1)

        self.label2 = tk.Label(self, text="Search Linking String").grid(row=1)

        self.entry3 = tk.Entry(self)
        self.entry3.grid(row=2,column=1)

        self.label3 = tk.Label(self, text="new folder path").grid(row=2)

        self.label4 = tk.Label(self, text="select test type").grid(row=3)

        self.n = tk.StringVar()
        self.testmachine = ttk.Combobox(self, width=27, textvariable=self.n)

        # Adding combobox drop down list
        self.testmachine['values'] = ('RNS',' Test 2 ',' Test 3')

        self.testmachine.grid(column=1, row=3)
        self.testmachine.current()

        self.button = tk.Button(self, text="autocrop", command=self.getinputs)
        self.button.grid(column = 2 , row = 5)

    def getinputs(self):
        autocrop(self.entry1.get(),self.entry2.get(),self.entry3.get(),self.testmachine.get())

guirun = button_work()
guirun.mainloop()

