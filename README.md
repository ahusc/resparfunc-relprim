# resparfunc-relprim
Sage Code for the Restricted Partition Function and Fourier-Dedekind Sums

This is a sagemath implementation of two algorithms to construct the restricted partition function for a given list of pairwise relatively prime positive numbers (also known als Sylvester's denumerant) and of an algorithm for the calculation of Fourier-Dedekind sums. 

The code implements the algorithms outlined in the paper "Restricted Partition Function and Fourier-Dedekind Sums".

For the construction of the restricted partition function there are two different algorithms. The first one constructs the partition function iteratively, starting with the construction of the partition function for the list that contains only the first element of the initial list. Then it constructs the partition function for the list containing the first two elements using the previously constructed partition function and so on.
The second one is based on the ability to construct the periodic functions which are part of the partition function independently. Having obtained the periodic functions first, it then calculates the function values for the partition function for several consecutive integers using a brute force method and finally obtains the polynomial part of the partition function using Lagrange's interpolation formula.


Fourier-Dedekind sums are number-theoretic objects that are closely
related to the periodical functions which are part of the representation
of the restricted partition function for the special case of pairwise
relatively prime numbers. A part of the second algorithm for the construction of the restricted partition function has been extracted and adapted for the calculation of Fourier-Dedekind sums.

# Authors

The code has been written by Anne Huschitt based on ideas of Alfred Seiler. These ideas are described in the paper "Restricted Partition Function and Fourier-Dedekind sums" (to be published). 

# License

The code is released under the GNU General Public License, version 2, or any later version as published by the Free Software Foundation. 

# Installation

A prerequisite for using this code is a properly working installation of [sage] (http://www.sagemath.org). After downloading the code, you may import one of the modules for the restricted partition function by typing:

```python
sage: load("resparfunc-relprim-algo2.sage")
```
or
```python
sage: load("resparfunc-relprim-algo1.sage")
```
The two versions of the algorithm cannot be used in the same sage session.

For the function to calculate Fourier-Dedekind sums, type

```python
sage: load("fourier-dedekind-sum.sage")
```

# Usage

After importing the code, type
```python
sage: RestrictedPartitionFunctionRelprim?
```
or, respectively
```python
sage: fourier_dedekind_sum?
```
to get information about the usage of the package, including various examples. The documentation can also be found on the [project website] (http://ahusc.github.io/resparfunc-relprim).

For the restricted partition function preferably use the second algorithm as it is considerably more performant in most cases than the first one.


