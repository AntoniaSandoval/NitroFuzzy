# -*- coding: utf-8 -*-
"""
Created on Mon May 21 14:46:09 2018

@author: Usuario
modulo fuzzy2

Advantages of FLSs
- Mathematical concepts within fuzzy reasoning are very simple
- You can modify a FLS by just adding or deleting rules due to flexibility of fuzzy logic
- Fuzzy logic Systems can take imprecise, distorted, noisy input information
- FLSs are easy to construct and understand
- Fuzzy logic is a solution to complex problems in all fields of life, including medicine as it resembles human reasoning and decision making

Disadvantages of FLSs
- There is nos sytematic approach to fuzzy system designing
- They are understandable only when simple
- They are suitable for the problems which do not need high accuracy

Applications
- Automatic control of dam gates for hydroelectric power plants
- Simplified control of robots
- Camera-aiming for the telecast of sporting events
- Efficient and stable control of car engines
- Cruise-control for automobiles
- Substitution of an expert for the assesment of stock exchange activities
- Optimised planning of bus timetables
- Archiving system for documents
- Prediction system for early recognition of earthquakes
- Medicine technoology: cancer diagnosis
- Recognition of motives in picturesn with video cameras
- Automatic motor-control for vacuum cleaners with recgnition of a surface condition and a degree of soiling
- Back-light control for camcoders

Fuzzy controllers applications
- Consumer products
* washing machines
* microwave ovens
* rice cookers
* vacuum cleaners
* camcoders
* TVs and VCRs
* thermal rugs
* word translations

- Systems
* elevators
* train
* cranes
* automotive (engines, transmissions, brakes)
* traffic control

- Software
* medical diagnosis
* securities
* data compression

"""

import numpy as np

"""===========================================================================
   Clase 1: Construccion de funciones de membresia
   ========================================================================"""
   
class Fuzzificador:

    def __init__(self,D,dx,alpha,numeroEtiquetas):
        
        self.D = float(D)
        self.dx = float(dx)
        self.alpha = float(alpha)
        self.numeroEtiquetas = int(numeroEtiquetas)-1
      
#=============================================================================
# Metodo 1. Calculo de coeficientes a y b (y = ax + b) Matriz con Netiq-1 cols 
#=============================================================================
        
    def calcularCoeficientes(self):
         
        coeficiente = np.empty((4,self.numeroEtiquetas)) 
        
        for etiqueta in range(self.numeroEtiquetas):
            
            coeficiente[0][etiqueta] = (self.D*(1+etiqueta)-self.dx+self.alpha)/(self.D-self.dx)
            coeficiente[1][etiqueta] = -1/(self.D-self.dx)
            coeficiente[2][etiqueta] = -(self.dx+self.D*etiqueta+self.alpha)/(self.D-self.dx)
            coeficiente[3][etiqueta] = 1/(self.D-self.dx)
            
        return coeficiente
    
    
#=============================================================================
# Metodo 2. Calculo de dominios de rectas y (Matriz con Netiq cols)
#============================================================================= 
        
    def calcularDominios(self):
        
        dominio = np.empty((4,self.numeroEtiquetas))
        
        for etiqueta in range(self.numeroEtiquetas):
            
            dominio[0][etiqueta] = self.D*etiqueta + self.alpha
            dominio[1][etiqueta] = self.D*(1+etiqueta) + self.alpha - self.dx
            dominio[2][etiqueta] = self.D*etiqueta + self.dx + self.alpha
            dominio[3][etiqueta] = self.D*(1+etiqueta) + self.alpha
            
        return dominio   

#=============================================================================
# Metodo 3. Calculo de grados de pertenencia sobre funciones de membresia
#=============================================================================

    def gpSfm(self,G,valorGP):

        etiquetas = np.zeros(self.numeroEtiquetas+1)
        
        for i in range(2):
            for etiqueta in range(self.numeroEtiquetas):
                if G[i][etiqueta] == 1:
                    etiquetas[i+etiqueta] = valorGP[i]
                    
        return etiquetas
        
     
#=============================================================================
# Metodo 4. Calculo de grados de pertenencia 
#=============================================================================     
    
    def calcularGP(self,dato):
        
        coeficiente = self.calcularCoeficientes()
        dominio = self.calcularDominios()
                
        G = np.zeros((4,self.numeroEtiquetas))        
        valorGP = np.zeros((2,1))
                
        for etiqueta in range(self.numeroEtiquetas):
            for i in range(3):
                if (dato >= dominio[i][etiqueta] and dato < dominio[i+1][etiqueta]):
                    
                    if (dominio[i+1][etiqueta] > dominio[i][etiqueta]):

                        G[i][etiqueta] = 1
                        G[i+1][etiqueta] = 1
                        
                        if i<2:
                            valorGP[0] = coeficiente[i][etiqueta] + coeficiente[i+1][etiqueta]*dato
                        elif i>=2:
                            valorGP[1] = coeficiente[i][etiqueta] + coeficiente[i+1][etiqueta]*dato
        
        G = np.vstack((G[0][:],G[2][:])) 
        etiquetas = self.gpSfm(G,valorGP)

        return etiquetas
         
  
#=============================================================================
# Metodo 5. Calculo de puntos de coordenadas de los puntos de cruce entre 
#           etiquetas
#=============================================================================
        
    def calcularCoordCruces(self):
        
        coeficiente = self.calcularCoeficientes()
        
        coord_x = np.zeros(self.numeroEtiquetas)
        coord_y = np.zeros(self.numeroEtiquetas)
        coord_xy = np.zeros((2,self.numeroEtiquetas))
        
        for etiqueta in range(self.numeroEtiquetas):
            coord_x[etiqueta] = (coeficiente[0][etiqueta]-coeficiente[2][etiqueta])/(coeficiente[3][etiqueta]-coeficiente[1][etiqueta])
            coord_y[etiqueta] = (coeficiente[0][etiqueta]+coeficiente[1][etiqueta]*coord_x[etiqueta])
            
        coord_xy = np.vstack((coord_x,coord_y))
        
        return coord_xy


"""===========================================================================
   Clase 2: Defuzzificador
   ========================================================================"""

class Defuzzificador:
    
    def __init__(self,numeroEtiquetas,GP,PC):
        
        self.numeroEtiquetas = numeroEtiquetas
        self.GP = GP
        self.PC = PC
 
#=============================================================================
# Metodo 1. Calculo de coeficientes a y b (y = ax + b) Matriz con Netiq cols 
#=============================================================================  
        
    def calcularCoeficientes2(self,coeficiente):
        
        matrizColZero = np.zeros((2,1))
        filasSuperior = np.hstack((matrizColZero,coeficiente[2:4][:]))
        filasInferior = np.hstack((coeficiente[0:2][:],matrizColZero))
        coeficiente2 = np.vstack((filasSuperior,filasInferior))
      
        return coeficiente2

    
#=============================================================================
# Metodo 2.1. Formula 1
#=============================================================================     
    
    def centroArea_1(self,i,sProd,sArea,pc,M):
        
        dx1 = self.PC[0][i]
        dx2 = (self.GP[i]-M[2][i])/M[3][i]
        dx3 = self.PC[0][i]-(self.GP[i]-M[2][i])/M[3][i]
        
        A1 = dx1*pc
        A2 = dx2*(self.GP[i]-pc)
        A3 = dx3*(self.GP[i]-pc)/2
        A = A1 + A2 + A3
        
        xg1 = 0.5*dx1
        xg2 = 0.5*dx2
        xg3 = 1/3*dx3
        
        xA1 = xg1*A1
        xA2 = xg2*A2
        xA3 = xg3*A3
        xA = xA1 + xA2 + xA3
        
        sProd = sProd + xA
        sArea = sArea + A
        
        return sProd, sArea
      
#=============================================================================
# Metodo 2.2. Formula 2
#=============================================================================     

    def centroArea_2(self,i,sProd,sArea):
        
        dx = self.PC[0][i]
        Ab = dx*self.GP[i]
        xg = 0.5*self.PC[0][i]
        xAb = xg*Ab
        
        sProd = sProd + xAb
        sArea = sArea + Ab  
        
        return sProd, sArea
  
#=============================================================================
# Metodo 2.3. Formula 3
#=============================================================================     

    def centroArea_3(self,i,sProd,sArea,pc,M):
        
        dx1 = self.PC[0][i]-(self.GP[i]-M[0][i+1])/M[1][i+1]
        A1 = dx1*(pc-self.GP[i])/2
        xg1 = (self.GP[i]-M[0][i+1])/M[1][i+1]+2/3*dx1
        xA1 = xg1*A1
        
        sProd = sProd + xA1
        sArea = sArea + A1        

        return sProd, sArea
    
#=============================================================================
# Metodo 2.4. Formula 4
#=============================================================================  

    def centroArea_4(self,i,sProd,sArea,M):
        
        dx1 = (self.GP[i+1]-M[0][i+1])/M[1][i+1]-(self.GP[i]-M[0][i+1])/M[1][i+1]
        dx2 = (self.GP[i+1]-M[0][i+1])/M[1][i+1]-self.PC[0][i]
        
        A1 = dx1*(self.GP[i+1]-self.GP[i])/2
        A2 = dx2*(self.GP[i+1]-self.GP[i])
        A = A1 + A2
        
        xg1 = (self.GP[i]-M[0][i+1])/M[1][i+1] + 2/3*dx1
        xg2 = 0.5*(self.PC[0][i]+(self.GP[i+1]-M[0][i+1])/M[1][i+1])
        
        xA1 = xg1*A1
        xA2 = xg2*A2
        xA = xA1 + xA2
        
        sProd = sProd + xA
        sArea = sArea + A

        return sProd, sArea        

#=============================================================================
# Metodo 2.5. Formula 5
#============================================================================= 

    def centroArea_5(self,i,sProd,sArea,pc,M):
        
        dx1 = (self.PC[0][i-1]+self.PC[0][0])-self.PC[0][i-1]
        dx2 = (self.PC[0][i-1]+self.PC[0][0])-(self.GP[i]-M[0][i])/M[1][i]
        dx3 = (self.GP[i]-M[0][i])/M[1][i]-self.PC[0][i-1]
        
        A1 = dx1*pc
        A2 = dx2*(self.GP[i]-pc)
        A3 = dx3*(self.GP[i]-pc)/2
        A = A1 + A2 + A3
        
        xg1 = self.PC[0][i-1]+0.5*self.PC[0][0]
        xg2 = 0.5*(self.PC[0][i-1]+self.PC[0][0]+(self.GP[i]-M[0][i])/M[1][i])
        xg3 = self.PC[0][i-1]+2/3*dx3
        
        xA1 = xg1*A1
        xA2 = xg2*A2
        xA3 = xg3*A3
        xA = xA1 + xA2 + xA3
        
        sProd = sProd + xA
        sArea = sArea + A

        return sProd,sArea

#=============================================================================
# Metodo 2.6. Formula 6
#=============================================================================

    def centroArea_6(self,i,sProd,sArea):
        
        dx = self.PC[0][i-1]+self.PC[0][0]-self.PC[0][i-1]
        Ab = dx*self.GP[i]
        xg = self.PC[0][i-1]+0.5*self.PC[0][0]
        xAb = xg*Ab
        
        sProd = sProd + xAb
        sArea = sArea + Ab  
        
        return sProd, sArea

#=============================================================================
# Metodo 2.7. Formula 7
#=============================================================================

    def centroArea_7(self,i,sProd,sArea,pc,M):
        
        dx1 = (self.GP[i]-M[2][i-1])/M[3][i-1]-self.PC[0][i-1]
        A1 = dx1*(pc-self.GP[i])/2
        xg1 = self.PC[0][i-1]+1/3*dx1
        xA1 = xg1*A1
                            
        sProd = sProd + xA1
        sArea = sArea + A1  
        
        return sProd, sArea

#=============================================================================
# Metodo 2.8. Formula 8
#=============================================================================

    def centroArea_8(self,i,sProd,sArea,M):
        
        dx1 = (self.GP[i]-M[2][i-1])/M[3][i-1]-(self.GP[i-1]-M[2][i-1])/M[3][i-1]
        dx2 = (self.GP[i-1]-M[2][i-1])/M[3][i-1]-self.PC[0][i-1]
        
        A1 = dx1*(self.GP[i-1]-self.GP[i])/2
        A2 = dx2*(self.GP[i-1]-self.GP[i])
        A = A1 + A2
        
        xg1 = (self.GP[i-1]-M[2][i-1])/M[3][i-1] + 1/3*dx1
        xg2 = 0.5*((self.GP[i-1]-M[2][i-1])/M[3][i-1]+self.PC[0][i-1])
        
        xA1 = xg1*A1
        xA2 = xg2*A2
        xA = xA1 + xA2
        
        sProd = sProd + xA
        sArea = sArea + A

        return sProd,sArea        

#=============================================================================
# Metodo 2.9. Formula 9
#=============================================================================

    def centroArea_9(self,i,sProd,sArea,pc,M):
        
        x1 = (self.GP[i]-M[0][i])/M[1][i]
        x2 = (self.GP[i]-M[2][i])/M[3][i]
        
        dx1 = self.PC[0][i]-self.PC[0][i-1]
        dx2 = x2-x1
        dx3 = x1-self.PC[0][i-1]
          
        A1 = pc*dx1
        A2 = dx2*(self.GP[i]-pc)
        A3 = dx3*(self.GP[i]-pc)
        A = A1 + A2 + A3
        
        xg = 0.5*(self.PC[0][i]+self.PC[0][i-1])
        xA = xg*A
        
        sProd = sProd + xA
        sArea = sArea + A

        return sProd,sArea 

#=============================================================================
# Metodo 2.10. Formula 10
#=============================================================================

    def centroArea_10(self,i,sProd,sArea):
        
        dx = self.PC[0][i]-self.PC[0][i-1]
        Ab = dx*self.GP[i]
        xg = 0.5*(self.PC[0][i-1]+self.PC[0][i])
        xAb = xg*Ab 
        
        sProd = sProd + xAb
        sArea = sArea + Ab        

        return sProd, sArea

#=============================================================================
# Metodo 2.11. Formula 11
#=============================================================================

    def centroArea_11(self,i,sProd,sArea,pc,M):
        
        dx1 = self.PC[0][i]-(self.GP[i]-M[0][i+1])/M[1][i+1]
        A1 = (pc-self.GP[i])*dx1/2
        xg1 = (self.GP[i]-M[0][i+1])/M[1][i+1]+2/3*dx1
        xA1 = xg1*A1
        
        sProd = sProd + xA1
        sArea = sArea + A1 
        
        return sProd,sArea

#=============================================================================
# Metodo 2.12. Formula 13
#=============================================================================

    def centroArea_13(self,i,sProd,sArea,pc,M):
        
        dx1 = (self.GP[i]-M[2][i-1])/M[3][i-1]-self.PC[0][i-1]
        A1 = dx1*(pc-self.GP[i])/2
        xg1 = self.PC[0][i-1] + 1/3*dx1
        xA1 = xg1*A1
        
        sProd = sProd + xA1
        sArea = sArea + A1        

        return sProd,sArea

#=============================================================================
# Metodo 2.13. Formula 14
#=============================================================================

    def centroArea_14(self,i,sProd,sArea,M):
        
        dx1 = (self.GP[i]-M[2][i-1])/M[3][i-1]-(self.GP[i-1]-M[2][i-1])/M[3][i-1]
        dx2 = (self.GP[i-1]-M[2][i-1])/M[3][i-1]-self.PC[0][i-1]
        
        A1 = dx1*(self.GP[i-1]-self.GP[i])/2
        A2 = dx2*(self.GP[i-1]-self.GP[i])
        A = A1 + A2
        
        xg1 = (self.GP[i-1]-M[2][i-1])/M[3][i-1]+1/3*dx1
        xg2 = 0.5*((self.GP[i-1]-M[2][i-1])/M[3][i-1]+self.PC[0][i-1])
        xA1 = xg1*A1
        xA2 = xg2*A2
        xA = xA1 + xA2
        
        sProd = sProd + xA
        sArea = sArea + A
        
        return sProd,sArea

#=============================================================================
# Metodo 2.14. Calcular salida difusa 
#============================================================================= 
        
    def salidaCentroArea(self,coef):
        
        pc = self.PC[1][0]        
        sumaProd = np.zeros(self.numeroEtiquetas)
        sumaArea = np.zeros(self.numeroEtiquetas)
        ajusteActuador = 0
        
        for etiqueta in range(self.numeroEtiquetas):
            
            if etiqueta == 0:

                if self.GP[etiqueta]>pc:
                
                    sumaProd[etiqueta], sumaArea[etiqueta] = self.centroArea_1(etiqueta,sumaProd[etiqueta],sumaArea[etiqueta],pc,self.calcularCoeficientes2(coef))

                else:
                                      
                    sumaProd[etiqueta], sumaArea[etiqueta] = self.centroArea_2(etiqueta,sumaProd[etiqueta],sumaArea[etiqueta])
                    
                    if self.GP[etiqueta]<=self.GP[etiqueta+1]:
                        
                        if self.GP[etiqueta+1]>pc:
                            
                            sumaProd[etiqueta], sumaArea[etiqueta] = self.centroArea_3(etiqueta,sumaProd[etiqueta],sumaArea[etiqueta],pc,self.calcularCoeficientes2(coef))
                         
                        else:
                            
                            sumaProd[etiqueta], sumaArea[etiqueta] = self.centroArea_4(etiqueta,sumaProd[etiqueta],sumaArea[etiqueta],self.calcularCoeficientes2(coef))
                                                
            elif etiqueta == self.numeroEtiquetas-1:

                if self.GP[etiqueta]>pc:
                    
                    sumaProd[etiqueta], sumaArea[etiqueta] = self.centroArea_5(etiqueta,sumaProd[etiqueta],sumaArea[etiqueta],pc,self.calcularCoeficientes2(coef))
                    
                else:
                    
                    sumaProd[etiqueta],sumaArea[etiqueta] = self.centroArea_6(etiqueta,sumaProd[etiqueta],sumaArea[etiqueta])
                                       
                    if self.GP[etiqueta]<=self.GP[etiqueta-1]:
                        
                        if self.GP[etiqueta-1]>pc:
                            
                            sumaProd[etiqueta], sumaArea[etiqueta] = self.centroArea_7(etiqueta,sumaProd[etiqueta],sumaArea[etiqueta],pc,self.calcularCoeficientes2(coef))

                        else:
                            
                            sumaProd[etiqueta],sumaArea[etiqueta] = self.centroArea_8(etiqueta,sumaProd[etiqueta],sumaArea[etiqueta],self.calcularCoeficientes2(coef))                   
                                     
            else:  
                                   
                if self.GP[etiqueta]>pc:
                    
                    sumaProd[etiqueta],sumaArea[etiqueta] = self.centroArea_9(etiqueta,sumaProd[etiqueta],sumaArea[etiqueta],pc,self.calcularCoeficientes2(coef))

                else:
                                        
                    sumaProd[etiqueta],sumaArea[etiqueta] = self.centroArea_10(etiqueta,sumaProd[etiqueta],sumaArea[etiqueta])
                                                           
                    if self.GP[etiqueta]<=self.GP[etiqueta+1]:
                        
                        if self.GP[etiqueta+1]>pc:
                            
                            sumaProd[etiqueta],sumaArea[etiqueta] = self.centroArea_11(etiqueta,sumaProd[etiqueta],sumaArea[etiqueta],pc,self.calcularCoeficientes2(coef))
                       
                        else:
                            
                            sumaProd[etiqueta],sumaArea[etiqueta] = self.centroArea_4(etiqueta,sumaProd[etiqueta],sumaArea[etiqueta],self.calcularCoeficientes2(coef))
                                                            
                    if self.GP[etiqueta]<=self.GP[etiqueta-1]:
                        
                        if self.GP[etiqueta-1]>pc:
                            
                            sumaProd[etiqueta],sumaArea[etiqueta] = self.centroArea_13(etiqueta,sumaProd[etiqueta],sumaArea[etiqueta],pc,self.calcularCoeficientes2(coef))

                        else:
                            
                            sumaProd[etiqueta],sumaArea[etiqueta] = self.centroArea_14(etiqueta,sumaProd[etiqueta],sumaArea[etiqueta],self.calcularCoeficientes2(coef))
                            
        salidaCentroArea = sum(sumaProd)/sum(sumaArea)-ajusteActuador
        
        return salidaCentroArea

"""===========================================================================
   Clase 3: Regla
   ========================================================================"""

class Regla:
    
    def __init__(self,rdir):
        self.rdir = rdir
    
#=============================================================================
# Metodo 3.1. extraccion de reglas difusas
#============================================================================= 

    def extraerReglas(self):
        
        f = open(self.rdir,"r")
        
        cont = 0
        rd = []
        
        for linea in f.readlines():
            
            rd.append(linea)
            rd[cont] = rd[cont].split(";")
            
            try:
                rd[cont][len(rd[cont])-1] = rd[cont][len(rd[cont])-1].replace("\r\n","")
            
            except:
                pass
            
            cont += 1
        
        f.close()
        
        return rd
    

#=============================================================================
# Metodo 3.2. Obtiene el numero de antecedentes
#============================================================================= 
    
    def obtenerNumAnt(self,rd,tokens):
      
        numAntecedentes = [None for i in range(len(rd))]
        
        for i in range(len(rd)):
            for j in range(len(rd[i])):
                
                if rd[i][j] == 'then':
                    numAntecedentes[i] = int(j - j/2)
        
        return numAntecedentes


#=============================================================================
# Metodo 3.3. Obtiene el numero de consecuentes
#============================================================================= 
   
    def obtenerNumCon(self,rd,tokens):
        
        numConsecuentes = [None for i in range(len(rd))]
    
        for i in range(len(rd)):
            for j in range(len(rd[i])):
                if rd[i][j] == 'then':
                    numConsecuentes[i] = int((len(rd[i])-j)/2)
            
        return numConsecuentes    

    
