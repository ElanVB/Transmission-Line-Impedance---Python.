# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 14:33:58 2015

@author: Elan van Biljon - 18384439
"""
from pylab import *

#a method to format matrices into strings in a nice way
def printComplexMatrix(m):
    stringMat = "["
    N, M = m.shape
    for i in range(N):
        for j in range(M):
            real = m[i, j].real
            imag = m[i, j].imag
            
            if (j == 0):
                stringMat += "["
            if (j != M - 1):
                if(imag >= 0):
                    stringMat += ("%.3g+%.3gj, " % (real, imag))
                else:
                     stringMat += ("%.3g %.3gj, " % (real, imag))
            else:
                if(imag >= 0):
                    stringMat += ("%.3g+%.3gj" % (real, imag))
                else:
                    stringMat += ("%.3g %.3gj" % (real, imag))
            
            if (j == M - 1):
                stringMat += "]"
        
        if(i != N - 1):
            stringMat += ",\n"
        else:
            stringMat += "]"
        
    return stringMat

#an oject with a x and y co-ordinate - usefull for working out the geometry of the problem
class point:
    __x = 0;
    __y = 0;
    
    #constructor
    def __init__(self, x, y):
        self.__x = x
        self.__y = y
        
    #returns the x cooord if the point
    def getX(self):
        return self.__x
        
    #returns the y coord of the point
    def getY(self):
        return self.__y

#calculates and returns the distance between two points
def D(p1, p2):
    x = p1.getX() - p2.getX()
    y = p1.getY() - p2.getY()
    return sqrt(x**2 + y**2)
    
#Calculates the distance between two conductors taking inductive earth images into account
def D_dash(p1, p2):
    x = p1.getX() - p2.getX()
    y1 = p1.getY() if (p1.getY() >= 0) else (abs(p1.getY()) - Dk_dash)
    y2 = p2.getY() if (p2.getY() >= 0) else (abs(p2.getY()) - Dk_dash)
    y = y1 - y2
    return sqrt(x**2 + y**2)
    
#some geometric properties of the phase conductor
phase_c_space = 15.8 #horizontal spacing
phase_c_hight = 37.3 #vertical height
bundle_r = 0.027 /2 #radius in metres
bundle_spacing = 0.32 #spacing between bundled strands
bundle_Rk = (72e-6 / 6)  * (128.0/117) * (740.0/711) #Resistance of line
bundle_GMR = 9.040168e-3 #effective radius for inductive calculations

#some geometric properties of the neutral conductor
neutral_c_space = 13.9  #horizontal spacing
neutral_c_hight = phase_c_hight + 11.25 #vertical height
neutral_r = 0.01 /2 #neutral radius in metres
neutral_Rk = 220e-6 * (128.0/117) * (214.0/213) #Resistance of line
neutral_GMR = 6.824412e-3 #effective radius for inductive calculations

p = 1000.0 #resistivity of the ground
f = 50.0 #frequency
w = 2 * math.pi * f #the frequency in radians per second
Deq = 0.0 #this variable will be used to store the Dkm' / Dkm terms
epsilon = 8.854e-12 #the constant epsilon
length = 290e3 #length of the transmission line in metres
Vr = Vs = 765e3 #kilo-Volts

#The following are paramaters if the transformation matrix and the matrix it's self
a = 1 * exp(1j * deg2rad(120)) #IN RADIANS
a2 = a**2
A = np.matrix([[1, 1, 1], [1, a2, a], [1, a, a2]])

Dk_dash = 658.5 * sqrt(p/f) #the efective inductive distance between earth conductors 
Rk_dash = 9.869e-7 * f #the effective resistance of the earth conductors

bundle_r_dash = bundle_r * exp(-0.25) #equivilent radius of the bundle strands
bundle_DsL = (6 * (bundle_spacing ** 5) * bundle_r_dash) ** (1.0/6) #the equivilent radius of the bundled conductor
bundle_Dsc = (6 * (bundle_spacing ** 5) * bundle_r) ** (1.0/6) #the capacitive radius of the bundled conductor

neutral_r_dash = bundle_r * exp(-0.25) #equivilent radius of the neutral strands#neutral_DsL = ((neutral_c_space*2)**2 * (neutral_r_dash**2))**0.5 #the equivilent radius of the neutral conductor
neutral_Dsc = bundle_r #the capacitive radius of the neutral conductor

#here I set the points' x and y equal to the conductors positions
#You will notice that there are 10 points - the extra 5 are the reflected earth conductors
points = range(0, 10)
points[0] = point(0, phase_c_hight)
points[1] = point(phase_c_space, phase_c_hight)
points[2] = point(phase_c_space * 2, phase_c_hight)
points[3] = point(phase_c_space - neutral_c_space, neutral_c_hight)
points[4] = point(phase_c_space + neutral_c_space, neutral_c_hight)
points[5] = point(0, -phase_c_hight)
points[6] = point(phase_c_space, -phase_c_hight)
points[7] = point(phase_c_space * 2, -phase_c_hight)
points[8] = point(phase_c_space - neutral_c_space, -neutral_c_hight)
points[9] = point(phase_c_space + neutral_c_space, -neutral_c_hight)

#here I create a few arrays or matracies
distances_for_capacitance = np.zeros((10,10)) #distance between conductors as seen by the capacitive calculations
distances_for_inductance = np.zeros((10,10)) #distance between conductors as seen by the inductive calculations
R_array = np.zeros((5,5)) #The R matrix [Sec 4.7 - equ 4.7.8]
L_array = np.zeros((5,5), dtype = complex) #the L matrix [Sec 4.7 - equ 4.7.8]
P_array = np.zeros((5,5)) #the P matrix [Sec 4.10]

#work out the distances between the conductors for the capacitive calculations
#if the counters k and m are equal the distance is the equivilent radius of the k conductor
for k in range(0,10):
    for m in range(0,10):
        distances_for_capacitance[k][m] = D(points[k], points[m]) if (k != m) else (bundle_Dsc if (k < 3) else neutral_Dsc)

#work out the distances between the conductors for the inductive calculations
#if the counters i and j are equal the distance is the equvilent radius of the i conductor
#else if either of the conductors are a reflected earth conductor the distance is the Dkk' value
#else the distance is the geometric distance between the conductors
for i in range(0, 10):
    for j in range(0, 10):
        distances_for_inductance[i][j] = (bundle_GMR if (i < 3) else neutral_r_dash) if(i == j) else D_dash(points[i], points[j])
        
#poulating the R matrix
#if the counters i and j are equal then the value is Rk' plus the resistivity of the conductor
#else the value is just Rkk'
for i in range(0, 5):
    for j in range(0, 5):
        R_array[i][j] = ((bundle_Rk if(i < 3) else neutral_Rk) + Rk_dash) if (i == j) else Rk_dash

#Populating the L matrix
#the if statements determine what Deq will be --> Deq = numerator / denominator
    #if the counters i and j are equal then Deq is the numerator
    #else the inductive distance between the conductors is the numerator
        #if the counters are lower than 3 (i or j is a phase conductor) then the bundle's DsL is the denominator
        #else the neutral's GMR is the denominator
for i in range(0, 5):
    for j in range(0, 5):
        Deq = distances_for_inductance[i][j + 5] / distances_for_inductance[i][j]
        L_array[i, j] = 2e-7 * log(Deq) * 1j * w #the formula to work out inductance
        
#now we get the impedence matrix by summing the R and the L matrix
Z_matrix = np.matrix(R_array + L_array, dtype = complex) 

#we then split the impedence matrix up into 4 parts
ZA = np.matrix(Z_matrix[:3, :3]) #(3 x 3) matrix
ZB = np.matrix(Z_matrix[:3, 3:]) #(3 x (N-3)) matrix
ZC = np.matrix(Z_matrix[3:, :3]) #((N-3) x 3) matrix
ZD = np.matrix(Z_matrix[3:, 3:]) #((N-3) x (N-3) matrix

#I now create the Zp matrix from [Sec 4.7 - equ 4.7.19]
Zp = np.matrix(ZA - (ZB.dot(inv(ZD))).dot(ZC), dtype = complex) 

#calculate the Zp_hat matrix paramaters 
Zaaeq_hat = (1.0/3)*(Zp[0, 0] + Zp[1, 1] + Zp[2, 2]) #[Sec 4.7 - equ 4.7.22]
Zabeq_hat = (1.0/3)*(Zp[0, 1] + Zp[0, 2] + Zp[1, 2]) #[Sec 4.7 - equ 4.7.23]

#create the Zp_hat matrix with Zabeq_hat as every element and the fill the diagonal wiht Zaaeq_hat
Zp_hat = np.matrix([Zabeq_hat for i in range(9)], dtype = complex).reshape(3, 3)
np.fill_diagonal(Zp_hat, Zaaeq_hat)

#small z is then calculated from 
z = Zaaeq_hat - Zabeq_hat #z = Zaaeq - Zabeq
print("z: %.3g+%.3gj" % (z.real, z.imag))

#the transformation matrix can be used as well
Zs = np.matrix((inv(A).dot(Zp)).dot(A))
temp_z = Zs[1, 1]
print("temp z: %.3g+%.3gj" % (temp_z.real, temp_z.imag))
#populate the P matrix
#Deq here is calculated as (the distance from conductor i to the earth conductor j) / (the actual geometric distance between conductors i and j) 
for i in range(0, 5):
    for j in range(0, 5):
        Deq = distances_for_capacitance[i][j + 5] / distances_for_capacitance[i][j]
        P_array[i, j] = (1.0 / (2 * math.pi * epsilon)) * log(Deq)  #the formula to get P

#we then split the P matrix up into 4 parts
PA = np.matrix(P_array[:3, :3]) #(3 x 3) matrix
PB = np.matrix(P_array[:3, 3:]) #(3 x (N-3)) matrix
PC = np.matrix(P_array[3:, :3]) #((N-3) x 3) matrix
PD = np.matrix(P_array[3:, 3:]) #((N-3) x (N-3) matrix

#I now create the Cp matrix from [Sec 4.10 - equ ]
Cp = np.matrix(inv(PA - (PB.dot(inv(PD))).dot(PC)), dtype = complex)

#calculate the Cp_hat matrix paramaters 
Caaeq_hat = (1.0/3)*(Cp[0, 0] + Cp[1, 1] + Cp[2, 2]) #[Sec 4.11 - equ 4.11.14]
Cabeq_hat = (1.0/3)*(Cp[0, 1] + Cp[0, 2] + Cp[1, 2]) #[Sec 4.11 - equ 4.11.15]

#create the Cp_hat matrix with Cabeq_hat as every element and the fill the diagonal wiht Caaeq_hat
Cp_hat = np.matrix([Cabeq_hat for i in range(9)], dtype = complex).reshape(3, 3)
np.fill_diagonal(Cp_hat, Caaeq_hat)

#get Yp and Yp_hat 
Yp = np.matrix(w * 1j * Cp, dtype = complex)
Yp_hat = np.matrix(w * 1j * Cp_hat, dtype = complex)

#the transformation matrix can be used as well
Ys = np.matrix((inv(A).dot(Yp)).dot(A))

#small y is then calculated from
y = Yp_hat[0, 0] - Yp_hat[0, 1] #y = Yaaeq - Yabeq
print("y: %.3g+%.3gj" % (y.real, y.imag))
temp_y = Ys[1, 1]
print("temp y: %.3g+%.3gj" % (temp_y.real, temp_y.imag))

#Gamma can now be calculated from [Sec 5.2 - equ 5.2.12]
gamma = sqrt(z * y)

#the characteristic impedence can now be calculated from [Sec 5.2 - equ 5.2.16]
Zc = sqrt(z / y)

#The term that goes into the hyperbolic trig functions can now be obtained
gl = gamma * length

#the ABCD paramaters can be calculated [Sec 5.2 - equ 5.2.34-36]
A_Par = cosh(gl)
print("A: %.3g+%.3gj" % (A_Par.real, A_Par.imag))
B_Par = Zc *  sinh(gl)
print("B: %.3g+%.3gj" % (B_Par.real, B_Par.imag))
C_Par = (1.0 / Zc) *  sinh(gl)
print("C: %.3g+%.3gj" % (C_Par.real, C_Par.imag))

#the circuit paramaters Z and Y can now be calculated - not necassary
Z = B_Par
Y = (A_Par - 1) / B_Par

#Caclulation of theoretical max Power can now be done
Pr = (Vs * Vr) / abs(Z) - ( abs(A_Par) * Vr**2) * cos(angle(Z) - angle(A_Par)) /  abs(Z) #[Sec 5.5 - equ 5.5.3]

#to check the validity of my paramater calculations, the following must be as close to zero as possible
Param_accuracy = (A_Par**2 - B_Par*C_Par)
print("AD-BC: %.3g%.3gj" % (Param_accuracy.real, Param_accuracy.imag))
