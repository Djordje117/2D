import os


def extract():
    Temp = []
    Coor = [] 
    Names = []

    if os.path.exists('POSCAR'):
        POS = open('POSCAR','r').readlines()
        for i in range(8,len(POS)):
            data = POS[i].split()
            Temp.append(float(data[0]))
            Temp.append(float(data[1]))
            Temp.append(float(data[2]))
            Names.append(data[3])
            Coor.append(Temp)
            Temp = []
    yield Coor
    yield Names
                    
def majority_minority_carriers(Coor = [], Names =[]):
    Min_car = Names[0]
    Max_car = Names[len(Names)-1]
    Min = 1
    Max = 0
    for i in range(1,len(Names)):
        if Names[i] == Min_car:
            Min+=1
        else:
            Max +=1
        if Max<Min:
            temp = Max
            Max = Min
            Min = temp
            Temp = Max_car
            Max_car = Min_car
            Min_car = Temp
    yield Min_car
    yield Max_car
    yield Max
    yield Min
    
def Reduce( Max_cor, Min_cor,Coor = [], Names = []):
    Red_coor = []
    Red_names = []
    difference = 1.2
    bond_z_lenght = 1.2
    check = True
    c = 0
    
    while check == True:
         if Names[0] == Min_cor:
            mark = True
            for i in range(len(Coor)):
                if Names[i]==Min_cor:
                    dif = float(Coor[i][2])-float(Coor[0][2])
                    if abs(dif) < difference:
                        c+=1
         elif Names[len(Names)-1] == Min_cor:
            mark = False
            for i in range(1,len(Coor)):
                if Names[i]==Min_cor:
                    dif = float(Coor[i][2])-float(Coor[len(Names)-1][2])
                    if abs(dif) < difference:
                        c+=1    
         if c!=2:
            print "Minority is ",c
            c=0
            difference -=0.01
         else:
            check = False
            
         if difference <= -0.01:
            raise SystemExit
            
    for i in range(1,len(Coor)):
        if Names[0] == Min_cor and abs(float(Coor[i][2])-float(Coor[0][2]))<difference and Names[i] == Names[0]:
            Red_coor.append(Coor[i])
            Red_coor.append(Coor[0])
            Red_names.append(Names[i])
            Red_names.append(Names[0])
            break
            
        elif Names[len(Names)-1] == Min_cor and abs(float(Coor[i][2])-float(Coor[len(Names)-1][2]))<difference and Names[i] == Names[len(Names)-1]:
            Red_coor.append(Coor[i])
            Red_coor.append(Coor[len(Names)-1])
            Red_names.append(Names[i])
            Red_names.append(Names[len(Names)-1])
            break

    check = True
    c=0
#Extracting majority carriers   
    while check == True:
        for i in range(len(Coor)):
            if Names[i]!=Min_cor:
                dif = float(Coor[i][2])-float(Red_coor[0][2])
                if abs(dif) < bond_z_lenght:

                        c+=1
        if c!=6:
            print "Majority is: ", c
            c=0
            bond_z_lenght -=0.01
        else:
            break

        if bond_z_lenght <= 0.01:
            raise SystemExit
    
    for i in range(len(Coor)):
        if Red_names[0] !=Names[i]:
            dif = float(Coor[i][2])-float(Red_coor[0][2])
            if abs(dif) < bond_z_lenght:
                Red_coor.append(Coor[i])
                Red_names.append(Names[i])
    
    yield Red_coor
    yield Red_names

def write_POSCAR(Min_car, Red_coor=[],Red_names = []):
    POS = open('POSCAR','r').readlines()
    Write = open('POSCAR','w')
    for i in range(0,6):
        Write.write(POS[i])
        
    elem = POS[5].split()
    if elem[0] ==Min_car:
        string = '2 6\n'
    else:
        string = '6 2\n'
    Write.write(string)
    Write.write(POS[7]) 
    for i in range(len(Red_coor)):
        string = str(Red_coor[i][0])+' '+str(Red_coor[i][1])+' '+str(Red_coor[i][2])+' '+Red_names[i]+'\n'
        Write.write(string)
    Write.close()
