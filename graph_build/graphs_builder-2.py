#!/usr/bin/env python
# coding: utf-8

# In[33]:


import numpy as np
import matplotlib.pyplot as plt


# In[170]:


def split_coords(str):
    tmp = str.split()
    x = []
    y = []
    for i in range(len(tmp)):
        if i % 2 == 0:
            x.append(float(tmp[i]))
        else:
            y.append(float(tmp[i]))
    return x,y


# In[174]:


dih = """-0.33 0.33
-0.190687 0.42689
-0.0726713 0.51732
0.00810188 0.58507
0.0538657 0.626146
0.0768309 0.647623
0.0876762 0.657986
0.0926611 0.662798
0.0949254 0.664995
0.0959487 0.66599
0.09641 0.666439
0.0966177 0.666641
0.0967112 0.666732
0.0967533 0.666773
0.0967722 0.666792
0.0967808 0.6668
0.0967846 0.666804
"""

gold = """"""
fib = """"""

alpha = """"""
dih_x, dih_y = split_coords(dih)
gold_x,gold_y = split_coords(gold)
fib_x,fib_y = split_coords(fib)
a_x,a_y = split_coords(alpha)


# In[179]:


#-4x - 2y + x^2 + y^2 + 5


delta = 0.01
x = np.arange(-1/2, 1/2, delta)
y = np.arange(0,1, delta)
X, Y = np.meshgrid(x, y)

Z =  np.sin(1/2*X**2 - 1/4*Y**2 +3)*np.cos(2*x + 1 + np.exp(y))

fig, ax = plt.subplots()
CS = ax.contour(X, Y, Z,levels = 5,colors ='black')
ax.clabel(CS, inline=True, fontsize=10)
fig.set_figwidth(10)    #  ширина и
fig.set_figheight(10)
plt.plot(gold_x,gold_y,color = 'red',label = 'golden')
plt.plot(fib_x,fib_y,color = 'blue', label = 'fib')
plt.plot(dih_x, dih_y,color = 'green', label = 'dih')
plt.scatter(dih_x, dih_y, color = 'green')
plt.plot(a_x,a_y,color = 'black', label = 'alpha_const')
plt.legend()
plt.show()


# In[ ]:





# In[ ]:




