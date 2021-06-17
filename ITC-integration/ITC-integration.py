#this program sends emails containing python scripts
import email, smtplib, ssl

from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

#email server Hosting command
#python -m smtpd -c DebuggingServer -n localhost:1025

def emailsend(sender_email,password,name):
    subject = "An email with attachment from Python"
    body = "Leave Form submitted by:"
    #sender_email = "pythontesterssmtp@gmail.com"
    receiver_email = "pythontesterssmtp@gmail.com"
    #password = input("Type your password and press enter:")

    # Create a multipart message and set headers
    message = MIMEMultipart()
    message["From"] = sender_email
    message["To"] = receiver_email
    message["Subject"] = subject
    message["Bcc"] = receiver_email  # Recommended for mass emails

    # Add body to email
    message.attach(MIMEText(body, "plain"))

    filename =name+'_Leave_request'+'.txt'  # In same directory as script

    # Open PDF file in binary mode
    with open(filename, "rb") as attachment:
     # Add file as application/octet-stream
        # Email client can usually download this automatically as attachment
         part = MIMEBase("application", "octet-stream")
         part.set_payload(attachment.read())

    # Encode file in ASCII characters to send by email
    encoders.encode_base64(part)

    # Add header as key/value pair to attachment part
    part.add_header(
    "Content-Disposition",
    f"attachment; filename= {filename}",)

    # Add attachment to message and convert message to string
    message.attach(part)
    text = message.as_string()

    # Log in to server using secure context and send email
    context = ssl.create_default_context()
    with smtplib.SMTP_SSL("smtp.gmail.com", 465, context=context) as server:
         server.login(sender_email, password)
         server.sendmail(sender_email, receiver_email, text)



#------------------------------------------------------------------------------------------------------------------------
#interactive leave text document to be sent to HR

import tkinter as tk
from tkinter import ttk


def instantiate(name,address,contactnum,datefrom,dateto,numdays,comment,typeofleave):
    tf =open(name+'_Leave_request'+'.txt',"w+")
    tf.write('Name of employee  :' +name)
    tf = open(name+'_Leave_request'+'.txt',"a")
    tf.write('\n')
    tf.write('Address during leave  :' +address)
    tf.write('\n')
    tf.write("Contact Number :" +contactnum)
    tf.write('\n')
    tf.write("Type of leave :" +typeofleave)
    tf.write('\n')
    tf.write("From(including) :"+datefrom)
    tf.write('\n')
    tf.write("To(including) :"+ dateto)
    tf.write('\n')
    tf.write("No. of working days  :"+ numdays)
    tf.write('\n')
    tf.write("Remarks  :"+ comment)


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
        instantiate(self.entname.get(),self.entaddress.get(),self.entcontactnum.get(),self.entdatefrom.get(),self.entdateto.get(),self.entnumdays.get(),self.entcomment.get(),self.cmbtypeofleave.get(),)
        emailsend(self.entemailname.get(),self.entemailpass.get(),self.entname.get())



GUIrun = GUI()
GUIrun.mainloop()


