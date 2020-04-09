from Fourier_Series_Class import *

"""
This is the example file. It contains an example of the the Fourier Series calculator calculating the Fourier series of exp(x) over the period [-2, 2].
"""

# Defining the function (Note: the calculator imported sympy as sp. You can overwrite this with importing it yourself)
func1 = Fourier_Func(sp.exp(x), 2)

# coefficents is a method that does what the name describes, prints out the indicated coefficent in its general form
def coefficent_example():
  func1.coefficients('b', printing = True)

# series is a method that returns the Fourier series of the indicated function up to the indicated number of terms
def series_example():
  func1.series(2, printing = True)

# plot is the method that plots the calculated Fourier series
def plot_example():
  func1.plot(2, printing = True)