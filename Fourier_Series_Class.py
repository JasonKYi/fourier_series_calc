import sympy as sp
import matplotlib.pyplot as plt
from sympy.abc import o
from progress.bar import ChargingBar
from halo import Halo
from os import system, name


def clear():

    if name == 'nt':
        _ = system('cls')

    else:
        _ = system('clear')


clear()

x = sp.Symbol('x')


class Fourier_Func:
    """A class of functions with methods for computing their Fourier Series.

    Parameters:
        func (sympy.core.power.Pow):
            func is a sympy function on which the Fourier coefficient (and thus the Fourier series) is evaluated on.

        period (float):
            The parameter period indicates the range on which the function func is periodic.
            As an example, given period = pi, then the function is periodic on [-pi, pi]

    """

    def __init__(self, func, period=sp.pi):
        n = sp.Symbol('n')

        self.func = func
        self.period = period
        self.n = n

    def coefficients(self, tpe, position_n=None, printing=False):
        """Returns the corresponding Fourier coefficients depending on the parameters.

        Parameters:
            tpe (string):
                For tpe == 'a' or 'b', the real Fourier coefficient of a or b at the specified position is printed / returned respectively.
                For tpe == 'c', the imaginary Fourier coefficient is printed / returned.

            position_n (int):
                position_n is the parameter that specifies where the Fourier coefficient is evaluated at.
                Leave empty for the general formulae.

            printing (bool):
                printing == True will print out the specified coefficient.
                printing == False with return the coefficient instead.        
        """

        if tpe != 'a' and tpe != 'b':
            raise NameError('Please enter a valid tpe a or b as strings!')

        if tpe == 'a':
            with Halo(text='Calculating a ', placement='right'):
                self.constant_a = (1 / self.period) * sp.integrate(sp.cos(
                    (self.n * sp.pi * x) / self.period) * self.func, (x, - self.period, self.period))
            clear()

            if position_n == None:
                if printing == True:
                    print(self.constant_a)
                elif printing == False:
                    return self.constant_a
            elif type(position_n) == int:
                if printing == True:
                    print(self.constant_a.subs(self.n, position_n).evalf())
                elif printing == False:
                    return self.constant_a.subs(self.n, position_n).evalf()

        if tpe == 'b':
            with Halo(text='Calculating b ', placement='right'):
                self.constant_b = (1 / self.period) * sp.integrate(sp.sin(
                    (self.n * sp.pi * x) / self.period) * self.func, (x, - self.period, self.period))
            clear()

            if position_n == None:
                if printing == True:
                    print(self.constant_b)
                elif printing == False:
                    return self.constant_b
            elif type(position_n) == int:
                if printing == True:
                    print(self.constant_b.subs(self.n, position_n).evalf())
                elif printing == False:
                    return self.constant_b.subs(self.n, position_n).evalf()

    def series(self, terms, printing=True):
        """Returns the Fourier series up to the specified number of terms.

        Parameter:
            terms (int):
                The parameter terms specifies the number of terms the Fourier series will have.
                As an example, if terms == 3, then the series method will output a Fourier series with three terms.

            printing (bool):
                printing == True will print out the calculated Fourier series.
                printing == False will not print out othe calculated Fourier series.

        """

        # Calculating if the function is even
        even_odd = 0
        for i in range(1, 100):
            if abs(self.func.subs(x, self.period * i / 100) + self.func.subs(x, - self.period * i / 100)) == 0:
                even_odd += 1
            elif abs(self.func.subs(x, self.period * i / 100) - self.func.subs(x, - self.period * i / 100)) == 0:
                even_odd -= 1
            else:
                break

        # Calculating the coefficients
        constant_a, constant_b = 0 * x, 0 * x

        if even_odd != 99:
            with Halo(text='Calculating a ', placement='right'):
                constant_a = (1 / self.period) * sp.integrate(sp.cos((self.n * sp.pi *
                                                                      x) / self.period) * self.func, (x, - self.period, self.period))
            clear()

        if even_odd != -99:
            with Halo(text='Calculating b ', placement='right'):
                constant_b = (1 / self.period) * sp.integrate(sp.sin((self.n * sp.pi *
                                                                      x) / self.period) * self.func, (x, - self.period, self.period))
            clear()

        self.four_func = 0.5 * constant_a.subs(self.n, 0)

        with ChargingBar('Appending series', max=terms) as bar:
            for k in range(1, terms + 1):
                self.four_func = self.four_func + constant_a.subs(self.n, k) * sp.cos(
                    k * sp.pi * x / self.period) + constant_b.subs(self.n, k) * sp.sin(k * sp.pi * x / self.period)
                bar.next()
        clear()

        if printing == True:
            print(self.four_func)

        return self.four_func

    def plot(self, terms, printing=True):
        """Plot the Fourier series.

        Parameter:
            terms (int):
                The parameter terms specifies the number of terms the Fourier series will have.
                As an example, if terms == 3, then the plot method will plot a Fourier series with three terms.

            printing (bool):
                printing == True will print out the calculated Fourier series.
                printing == False will not print out othe calculated Fourier series.

        """

        sp.plot(self.series(terms, printing=printing), show=True)
