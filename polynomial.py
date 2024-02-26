def flatten(data):
    """Function to flatten nested tuples

    Parameters
    ----------
    data : tuple
        a tuple consisting of nested tuples
    
    Returns
    -------
    tuple :
        a flatten version of a nested tuple
    """
    if isinstance(data, tuple):
        for x in data:
            yield from flatten(x)
    else:
        yield data


def zip_longest(iter1, iter2, fillvalue=None):
    """Works similarly to zip for two parameters, but does not stop
    when one of the iterators is exhausted, but uses "fillvalue" instead.
    Alternatively use:
        from itertools import zip_longest

    Parameters
    ----------
    iter1 : tuple
        first parameter to zip

    iter2 : tuple
    
    fillvalue : int
        Default = None

    Returns
    -------
    tuple
        of two zipped parameters
    """
    for i in range(max(len(iter1), len(iter2))):
        if i >= len(iter1):
            yield (fillvalue, iter2[i])
        elif i >= len(iter2):
            yield (iter1[i], fillvalue)
        else:
            yield (iter1[i], iter2[i])
        i += 1


class Polynomial:
    """Class to define and calculate poylnomials of the type:
       p(x) = an*x^n + a(n-1)*x^(n-1) + ... + a1*x^1 + a0*x^0
            = \sum\lim_{k=0}^{n}\left(a_{k}x^{k}\right)
       and the derivative of:
       p'(x) = n*an*x^(n-1) + (n-1)*a(n-1)*x^(n-2) + ... + 2*a2*x^1+ a1*x^0
             = \sum\lim_{k=0}^{n}\left(a_{k}x^{k}\right)

    Reference
    ---------
    Bronstein et. al.: Taschenbuch der Mathematik. 7th edition, 2008,
                       Harri Deutsch GmbH.
    """


    def __init__(self, *coeff):
        """
        Parameters
        ----------
        coeff : array
            coefficients of the polynomial in the order from
            the largest to the smallest coefficient:
            [an, a(n-1), ..., a2, a1, a0]
        """
        self.coeff = list(coeff)


    def __repr__(self):
        """Method to return the canonical string representation 
        of a polynomial.
        
        Returns
        -------
        tuple :
            of the polynomial coefficients
        """
        return "Polynomial" + str(self.coeff)


    def degree(self):
        """Method to return the degree of the polynomial

        Returns
        -------
        int :
            degree of the polynomial
        """
        return len(self.coeff)


    def __str__(self):
        """Method for returning the polynomial as a string for better
        visualisation and understanding.

        Returns
        -------
        str :
            a consistent sequence of coefficients and their degree of x,
            zero coefficients are neglected.
        """
        
        def x_expr(degree):
            """function to add the "x^",
               into the polynomial string representation
            
            Parameters
            ----------
            degree : int
                the degree of the polynomial

            Returns
            -------
            str :
                an empty string, x or x^ in order of the degree
            """
            if degree == 0:
                res = ""
            elif degree == 1:
                res = "x"
            else:
                res = "x^" + str(degree)
                
            return res


        degree = len(self.coeff) - 1
        res = ""

        for i in range(0, degree+1):
            # nothing has to be done if coeff is 0:
            if abs(self.coeff[i]) == 1 and i < degree:
                # 1 in front of x shouldn't occur, e.g. x instead of 1x
                # but we need the plus or minus sign:
                res += f"{'+' if self.coeff[i]>0 else '-'}{x_expr(degree-i)}"  
            elif self.coeff[i] != 0:
                res += f"{self.coeff[i]:+g}{x_expr(degree-i)}" 

        return res.lstrip("+")  # removing leading '+'


    def __call__(self, x):
        """Calculates the Polynomial Function at the point x acc. to
        Bronstein et. al. (2008), see Horner's method.
        
        Parameters
        ----------
        x : float
            Point at which the polynomial function is to be calculated.
        
        Returns
        -------
        float :
            for the value of the polynomial function at point x.
        """
        res = 0
        for coefficient in self.coeff:
            res = res * x + coefficient
        #res = [res * x + c for c in self.coeff]
        return res
    

    def __add__(self, other):
        """Method to add two polynoms

        Parameters
        ----------
        other : tuple
            another polynom of the type:
            p2 = Polynomial(...)
        
        Returns
        -------
        list :
            fo the combined polynomials
        """
        p1 = self.coeff
        p2 = other.coeff
        res = [sum(f) for f in zip_longest(p1, p2, fillvalue=0)]
        return Polynomial(*res)

    
    def __sub__(self, other):
        """Method to subtract two polynoms

        Parameters
        ----------
        other : tuple
            another polynom of the type:
            p2 = Polynomial(...)
        
        Returns
        -------
        list :
            fo the subtracted polynomials
        """
        p1 = self.coeff
        p2 = other.coeff
        
        res = [f-g for f, g in zip_longest(p1, p2, fillvalue=0)]
        return Polynomial(*res)
    

    def derivative(self):
        """Method to calculate the derivative of the polynomial acc. to
        Bronstein et. al. (2008), see Horner's method.

        Returns
        -------

        """
        derived_coeffs = []
        exponent = len(self.coeff) - 1
        for i in range(exponent):
            derived_coeffs.append(self.coeff[i] * exponent)
            exponent -= 1
        return Polynomial(*derived_coeffs)