tfDatabasename = "D:\ITC work\Database Python\DataBase.txt"





tf2 = open("organizedDB.txt", "w+")
tf3 = open("organizedDBCodes.txt", "w+")
tf4 = open("organizedDBText.txt", "w+")
iloop=0
iloop2=0
arr2D=[[0]*2]*icount
arr2D2=[[0]*2]*icount

for iloop in range(icount):
    arr2D[iloop][0]=arrofstr1
    arr2D[iloop][1]=arrofstr2
    tf2.write(arrofstr1[iloop]+"@"+arrofstr2[iloop])
    tf3.write(arrofstr1[iloop])
    tf4.write(arrofstr2[iloop])

tf2.close()
tf3.close()
tf4.close()



