
# coding: utf-8

# In[66]:

#README
#Down the line we will be looking for peaks for gamma spectroscopy purposes. So to prepare for this, I began the 
#framework for a script that will find the peaks of graphs. 


# In[67]:

#import numpy and matplotlib
import numpy as np
get_ipython().magic('matplotlib inline')
import matplotlib.pyplot as plt


# In[68]:

#Lets start with a quadratic 
#One way to plot functions is to go completely discrete
#we use the np.linspace function to create a list from -2 to 2 with 50 entries
x = np.linspace(-2, 2, 50)
def test_func_1(x): 
    '''
    This just defines the quadratic
    '''
    return -x**2 + 1
y = [] 
i = 0
#we iterate through all of x and make a new list with the values of the function for each entry in x
while i < 50:
    y.append(test_func_1(x[i]))
    i+=1
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')


# In[69]:

def peak_finder(y):
    '''
    This is the main function that needs development. 
    Essentially it looks through three consequative values of y at the time and checks in the middle value is a max. 
    If it is a max, it saves it to a peak list.
    '''
    peaks = []
    i, j, k = 0, 1, 2
    counter = 0
    first, second, third = y[i], y[j], y[k]
    while k < len(y): 
        if second >= third: 
            if second >= first: 
                peaks.append([j, second])
        
        i += 1
        j += 1
        k += 1
        if k < len(y):
            first, second, third = y[i], y[j], y[k]
        else: 
            pass
        

    return peaks

peak = peak_finder(y)
print(peak)


# In[70]:

def plotting_max(x, peak): 
    '''
    This function just automates part of the process
    '''
    length = len(peak)
    i = 0
    horizontal = []
    vertical = []
    argument = []
    while i < length: 
        horizontal.append(x[peak[i][0]])
        vertical.append(peak[i][1])
        i+=1
    return horizontal, vertical
print(plotting_max(x, peak))
argument = plotting_max(x, peak)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(argument[0], argument[1], 'ro')
plt.axis([-2.5, 2.5, -3, 1.5])
#The Red circles are maxes. Why do we get two?


# In[71]:

#Lets use np.random.randint
#Look up documentation to see what it does
y = np.random.randint(10, size=50)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')


# In[72]:

peak = peak_finder(y)
argument = plotting_max(x, peak)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(argument[0], argument[1], 'ro')
#We have some double peaks once again. What should we do with those? Also some of these peaks are pretty small, 
#How should we deal with the small peaks?


# In[73]:

x = np.linspace(-2, 2, 1000)
y = np.random.randint(100, size=1000)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')


# In[74]:

peak = peak_finder(y)
argument = plotting_max(x, peak)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(argument[0], argument[1], 'ro')
#It can get chaotic.....


# In[75]:

x = np.linspace(-2, 2, 1000)
y = np.random.randint(10, size=1000)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')


# In[76]:

peak = peak_finder(y)
argument = plotting_max(x, peak)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(argument[0], argument[1], 'ro')


# In[77]:

x = np.linspace(-2, 2, 500)
y = np.random.randint(10, size=500)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')


# In[78]:

peak = peak_finder(y)
argument = plotting_max(x, peak)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(argument[0], argument[1], 'ro')


# In[79]:

x = np.linspace(-2, 2, 250)
y = np.random.randint(10, size=250)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')


# In[80]:

peak = peak_finder(y)
argument = plotting_max(x, peak)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(argument[0], argument[1], 'ro')


# In[81]:

x = np.linspace(-2, 2, 120)
y = np.random.randint(10, size=120)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')


# In[82]:

peak = peak_finder(y)
argument = plotting_max(x, peak)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(argument[0], argument[1], 'ro')


# In[83]:

#More accurate Spectra generator

def random_spectra_generator(spikes, length): 
    '''
    Generates a more accurate spectra for testing purposes. 
    There is length entries. 
    Spikes is how many peaks there should be. 
    '''
    higher = length//3
    i = 0 
    spectra = []
    while i < higher: 
        spectra.append(np.random.randint(low=45, high=55))
        i+=1
    while i < length: 
        spectra.append(np.random.randint(low=0, high=10))
        i+=1
    j = 0 
    while j < spikes: 
        spectra[np.random.randint(low=0, high=length-1)] = np.random.randint(low=40, high=100)
        j+=1
    return spectra


# In[84]:

x = np.linspace(0, 1800, 120)
y = random_spectra_generator(1, 120)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')


# In[85]:

peak = peak_finder(y)
argument = plotting_max(x, peak)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(argument[0], argument[1], 'ro')


# In[86]:

#So the prior graph is a demonstration of what a spectra might look like. 
#We have all this noise which leads to the current iteration of the peak finder to add extra peaks. 
#We have to deal with this noise somehow. 


# In[87]:

#I will do a couple more examples then there will be instructions how to run the code in its current state. 
x = np.linspace(0, 1800, 120)
y = random_spectra_generator(5, 120)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')


# In[88]:

peak = peak_finder(y)
argument = plotting_max(x, peak)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(argument[0], argument[1], 'ro')


# In[89]:

length = 1000
x = np.linspace(0, 1800, length)
y = random_spectra_generator(5, length)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')


# In[90]:

peak = peak_finder(y)
argument = plotting_max(x, peak)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(argument[0], argument[1], 'ro')


# In[91]:

length = 10000
x = np.linspace(0, 1800, length)
y = random_spectra_generator(10, length)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')


# In[92]:

peak = peak_finder(y)
argument = plotting_max(x, peak)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(argument[0], argument[1], 'ro')


# In[93]:

#Instructions how to run: 
#There are three main functions that I already designed for you
#You can modify them how you see fit
def peak_finder(y):
    '''
    This is the main function that needs development. 
    Essentially it looks through three consequative values of y at the time and checks in the middle value is a max. 
    If it is a max, it saves it to a peak list.
    '''
    peaks = []
    i, j, k = 0, 1, 2
    counter = 0
    first, second, third = y[i], y[j], y[k]
    while k < len(y): 
        if second >= third: 
            if second >= first: 
                peaks.append([j, second])
        
        i += 1
        j += 1
        k += 1
        if k < len(y):
            first, second, third = y[i], y[j], y[k]
        else: 
            pass
        

    return peaks

def plotting_max(x, peak): 
    '''
    This function just automates part of the process
    '''
    length = len(peak)
    i = 0
    horizontal = []
    vertical = []
    argument = []
    while i < length: 
        horizontal.append(x[peak[i][0]])
        vertical.append(peak[i][1])
        i+=1
    return horizontal, vertical

def random_spectra_generator(spikes, length): 
    '''
    Generates a more accurate spectra for testing purposes. 
    There is length entries. 
    Spikes is how many peaks there should be. 
    '''
    higher = length//3
    i = 0 
    spectra = []
    while i < higher: 
        spectra.append(np.random.randint(low=45, high=55))
        i+=1
    while i < length: 
        spectra.append(np.random.randint(low=0, high=10))
        i+=1
    j = 0 
    while j < spikes: 
        spectra[np.random.randint(low=0, high=length-1)] = np.random.randint(low=40, high=100)
        j+=1
    return spectra


# In[94]:

#After you have your functions defined, you need to first generate data
length = 1000
x = np.linspace(0, 1800, length)
y = random_spectra_generator(7, length)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
#This will generate a sample spectra with a length that you can define. You can also define how many spikes that will 
#occur. Each spikes intensity will be random, and the location is also random. I also plot it for good measure. 
#For this example I am just going to use length=1000 and spikes=7


# In[95]:

#Now is when we use the peak finder. After we get the peaks, we also need to add the circles to the graph. 
peak = peak_finder(y)
argument = plotting_max(x, peak)
plt.plot(x, y)
plt.xlabel('x')
plt.ylabel('y')
plt.plot(argument[0], argument[1], 'ro')


# In[96]:

#We should only get 7 red circles but the randomly generated noise causes more circles to appear. 
#Also there should be some nice way to demonstrate which peak corresponds to which values