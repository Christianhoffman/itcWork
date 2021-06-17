#interactive leave word document to be sent to HR
import tkinter as tk
from tkinter import ttk



def instantiate(name,address,contactnum,datefrom,dateto,numdays,comment,typeofleave,emailname,emailpass):
    tf =open(name+' Leave request'+'.txt',"w+")
    tf.write(name)
    tf = open(name+' Leave request'+'.txt',"a")
    tf.write('\n')
    tf.write(address)






class GUI(tk.Tk):
    def __init__(self):
        tk.Tk.__init__(self)

        self.entname = tk.Entry()
        self.entname.grid(row=0, column=1)

        self.label1 = tk.Label(self, text="Name of employee").grid(row=0)

        self.entaddress = tk.Entry(self)
        self.entaddress.grid(row=1, column=1)

        self.label2 = tk.Label(self, text="Address during leave").grid(row=1)

        self.entcontactnum = tk.Entry()
        self.entcontactnum.grid(row=2, column=1)

        self.label3 = tk.Label(self, text="Contact Number").grid(row=2)

        self.entdatefrom = tk.Entry()
        self.entdatefrom.grid(row=1, column=3)

        self.label4 = tk.Label(self, text="Date From ").grid(row=0,column=3)

        self.entdateto= tk.Entry()
        self.entdateto.grid(row=1, column=4)

        self.label5 = tk.Label(self, text="Date to ").grid(row=0,column=4)

        self.entnumdays = tk.Entry()
        self.entnumdays.grid(row=1, column=5)

        self.label6 = tk.Label(self, text="no. working days ").grid(row=0,column=5)

        self.entcomment = tk.Entry()
        self.entcomment.grid(row=4, column=1)

        self.label7 = tk.Label(self, text="Comment about leave").grid(row=4,column=0)

        self.label13 = tk.Label(self, text="type of leave").grid(row=3,column=0)

        self.n = tk.StringVar()
        self.cmbtypeofleave = ttk.Combobox(self, width=27, textvariable=self.n)

        # Adding combobox drop down list
        self.cmbtypeofleave['values'] = ('Special leave', 'Annual leave', 'Sick Leave','unpaid leave')

        self.cmbtypeofleave.grid(column=1, row=3)
        self.cmbtypeofleave.current()

        self.button = tk.Button(self, text="Submit", command=self.accessor)
        self.button.grid(column = 0 , row = 7)

        self.entemailname = tk.Entry()
        self.entemailname.grid(row=5, column=1)

        self.label8 = tk.Label(self, text="email of sender").grid(row=5, column=0)

        self.entemailpass = tk.Entry()
        self.entemailpass.grid(row=6, column=1)

        self.label9 = tk.Label(self, text="password of sender").grid(row=6, column=0)

    def accessor(self):
        instantiate(self.entname.get(),self.entaddress.get(),self.entcontactnum.get(),self.entdatefrom.get(),self.entdateto.get(),self.entnumdays.get(),self.entcomment.get(),self.cmbtypeofleave.get(),self.entemailname.get(),self.entemailpass.get())




GUIrun = GUI()
GUIrun.mainloop()







