# -*- coding: utf-8 -*-
#############################   Code description   ############################


## This code defines a class to read .csv files and allow its separation in
## columns, rows and a matrix, instead of using directly a data frame.
## Created by: Joel Alpízar Castillo.
## TU Delft
## Version: 2.1
##  - 1.1: the data2matrix() function was added.
##  - 2.0: the for loops in the data2cols(), data2rows() and data2matrix() are
##      now nested.
##  - 2.1: the data2array() function was added.



###############################################################################
###########################   Imported modules   ##############################

from itertools import product

############################   Classes definition  ############################

## The class read_data() allows to read the data from a .csv file and store it.

## For the class read_data, the description of the parameters is the next:
##  csv: (string) holds the name of the .csv file.
##  delim: a string with the separator of the values in the .csv file. If not
##      indicated, a semicolon would be assumed as default.
##  address: (string) holds the route of the .csv file.

## Considerations:
##  - The name of the file, when indicated in the variable "csv" must include 
##    the extension, e.g. "csv='file.csv'".
##  - If not indicated, the default address for the files is a folder named
##    "CSVs".
##  - The index of the data is: self.cols[column][row]


class read_data:
    def __init__(self, csv='', delim=';', address='CSVs\\'):
        from pandas import read_csv
        self.csv=csv
        self.delim=delim
        self.doc=''.join([address,csv])
        self.df=read_csv(self.doc, sep=delim, header=None)
        self.df.values        
        self.nrows=len(self.df)
        self.ncols=int(self.df.size/self.nrows)
    ## End of def __init__(self, csv='', delim=';', address='CSVs\\'):

## The function data2cols reads the data frame stored in self.df and create a
## list whose elements are lists of the columns. i.e. the first element is a
## list that stores the elements of the first column of the data. The resulting
## list is stored in self.col.    

    def data2cols(self):
        self.col=[]                                            
        for __i in range(self.ncols):
            self.col.append([0])
        ## End of for __i in range(__ncols):
        
        for __i, __j in product(range(self.nrows), range(self.ncols)):
            self.col[__j].append(self.df[__j][__i])
        ## End of for __i, __j in product(range(self.nrows), range(self.ncols)):
                
        for __i in range(self.ncols):
            self.col[__i].pop(0)               
        ## End of for __i in range(__ncols):
        
    ## End of def data2cols(self):

## The function data2rows reads the data frame stored in self.df and create a
## list whose elements are lists of the rows. i.e. the first element is a
## list that stores the elements of the first row of the data. The resulting
## list is stored in self.row.  
    
    
    def data2rows(self):
        self.row=[]        
        for __i in range(self.nrows):
            self.row.append([0])
        ## End of for __i in range(__nrows):
        
        for __i, __j in product(range(self.nrows), range(self.ncols)):
            self.row[__i].append(self.df[__j][__i])
        ## End of for __i, __j in product(range(self.nrows), range(self.ncols)):
                
        for __i in range(self.nrows):
            self.row[__i].pop(0)               
        ## End of for __i in range(__nrows):

    ## End def data2rows(self):
    
## The function data2matrix reads the data frame stored in self.df and create
## a matrix with the data stored in self.df The resulting matrix is stored
## in self.mat.
    
    def data2matrix(self):
        from numpy import matrix, zeros       
        self.mat=matrix(zeros((self.nrows,self.ncols)))        
        
        for __i, __j in product(range(self.nrows), range (self.ncols)):
            self.mat[__i,__j]=self.df[__j][__i]
        ## End of for i, j in product(range(self.nrows),range (self.ncols)):
        
    ## End of def data2matrix(self):

## The function data2array reads the data frame stored in self.df and create
## an array with the data stored in self.df The resulting array is stored
## in self.ar.
    
    def data2array(self):
        from numpy import array, zeros       
        self.ar=array(zeros((self.nrows,self.ncols)))        
        
        for __i, __j in product(range(self.nrows), range (self.ncols)):
            self.ar[__i,__j]=self.df[__j][__i]
        ## End of for i, j in product(range(self.nrows),range (self.ncols)):
        
    ## End of def data2matrix(self):


## End of class read_data:   

###############################################################################