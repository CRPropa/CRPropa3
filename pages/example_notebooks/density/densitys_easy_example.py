
# coding: utf-8

# # first example code for use of densitys
# 
# #### all densitys are seperated in atomic hydrogyn (HI), ioninised hydrogyn (HII) and molecular hydrogn (H2)
# 
# ### this example contains 
# * constant densitys
# * superposition of different typses 
# * superposition of same types
# 
# Further on there is another example for all implemented density models of the Milky-Way

# In[4]:

from crpropa import *


# ### first we start with the use of a constant density 
# the first double set the HI, the second the HII and the third the H2 density-number in SI units

# In[5]:

CD = ConstantDensity(1,2,0.5)


# to see the output option we check the density at a random position
# 
# the output options are:
# * getDensity: for the sum of all densitys
# * getHIDensity: for the HI part
# * getHIIDensity: for the HII part
# * getH2Density: for the H2 part
# * getNucleonDensity: for the sum of nuclei ($n_{HI} + n_{HII} + 2\cdot n_{H2}$)
# 

# In[6]:

position = Vector3d(2,1,3)

n_tot = CD.getDensity(position)
n_HI = CD.getHIDensity(position)
n_HII = CD.getHIIDensity(position)
n_H2 = CD.getH2Density(position)
n_nucl = CD.getNucleonDensity(position)


print('total density n = %f' %n_tot)
print('HI density n_HI = %f' %n_HI)
print('HII density n_HII = %f' %n_HII)
print('H2 density n_H2 = %f' %n_H2)
print('nucleon density n_nucl = %f' %n_nucl)


# #### the ConstantDensity can be adjusted to new values and the useage of several parts can be chosen
# therefore are methodes to change and activate (set-function) and methodes to see actual configuration (getisfor-functions)

# In[7]:

#see the actual configuration
print('HI: %s, HII: %s, H2: %s, \n \n' %(CD.getIsForHI(), CD.getIsForHII(), CD.getIsForH2()))

# change activity
CD.setHI(False)

# change activity and density number
CD.setHII(False, 1.5)

# change density number
CD.setH2(1.3)

# see the changes in the Description
print(CD.getDescription())


# #### the output of the getDensity and getNucleonDensity depends on the activity of the types. Only acitvated types are used for summing up

# In[8]:

n_tot = CD.getDensity(position)
n_HI = CD.getHIDensity(position)
n_HII = CD.getHIIDensity(position)
n_H2 = CD.getH2Density(position)
n_nucl = CD.getNucleonDensity(position)


print('total density n = %f' %n_tot)
print('HI density n_HI = %f' %n_HI)
print('HII density n_HII = %f' %n_HII)
print('H2 density n_H2 = %f' %n_H2)
print('nucleon density n_nucl = %f' %n_nucl)


# ## customize a density with different density models
# To customize a density use the DensityList.
# In a superposition of globel models of density distribution of the Milkyway keep care of normalisation. Therfore you can just add components by deactivating the other  

# In[9]:

CD1 = ConstantDensity(0,2,0)     # for use HII
CD2 = ConstantDensity(3,1,2.5)   # for use HI, H2

CUS = DensityList()

# first deactivate not wanted parts

CD1.setHI(False)
CD1.setH2(False)

CD2.setHII(False)  

# add density 

CUS.addDensity(CD1)
CUS.addDensity(CD2)

# get output

n_tot = CUS.getDensity(position)
n_HI = CUS.getHIDensity(position)
n_HII = CUS.getHIIDensity(position)
n_H2 = CUS.getH2Density(position)
n_nucl = CUS.getNucleonDensity(position)

print('total density n = %f' %n_tot)
print('HI density n_HI = %f' %n_HI)
print('HII density n_HII = %f' %n_HII)
print('H2 density n_H2 = %f' %n_H2)
print('nucleon density n_nucl = %f' %n_nucl)


# ### DensityList
# you can also superposition total models without deactivating several components

# In[10]:

#set wanted density
CD1 = ConstantDensity(1,3,4)
CD2 = ConstantDensity(1.5,2,3.3)

DL = DensityList()

# add density to list

DL.addDensity(CD1)
DL.addDensity(CD2)

# see output
n_tot = DL.getDensity(position)
n_HI = DL.getHIDensity(position)
n_HII = DL.getHIIDensity(position)
n_H2 = DL.getH2Density(position)
n_nucl = DL.getNucleonDensity(position)

print('total density n = %f' %n_tot)
print('HI density n_HI = %f' %n_HI)
print('HII density n_HII = %f' %n_HII)
print('H2 density n_H2 = %f' %n_H2)
print('nucleon density n_nucl = %f' %n_nucl)


# In[ ]:




# In[ ]:



