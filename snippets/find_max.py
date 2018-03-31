from matplotlib.pyplot import plot, show
from numpy import array
import operator
from detect_peaks import detect_peaks
from numpy.random import random


freq = array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19])
power = array([10, 8, 3, 2, 1, 30, 31, 32, 35, 30, 29, 25, 50, 9, int(random()*25), int(random()*25), int(random()*25), 40, int(random()*25)])
print power

# index = detect_peaks(power, mpd=4, show=True)
# value = power[index]
index, value = max(enumerate(power), key=operator.itemgetter(1))
print (index, value)
print freq[index]
print "Average power"
print (power[index] + power[index+1] + power[index-1])/3.

# plot(freq, power)
# plot(index+1, value, '+')
# show()
