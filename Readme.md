# MPAWL for MatLab =
This `MatLab` Library is an implementation of the periodic Wavelet Transform based on an integral regular matrix __M__ and its factorization into dilation matrices. Introducing the multivariate de la Vallée Poussin means, this Library provides many scaling functions for the levels of decomposition. This package  was transcribed from the Mathematica version, which is available at [https://github.com/kellertuer/MPAWL](https://github.com/kellertuer/MPAWL)

The underlying theory of patterns, its generating groups and the fast Fourier transform on these patterns is also implemented in this Library yielding a fast Wavelet Transform, when computing in Fourier coefficients. Further, several functions to work with the translates of a function with respect to the pattern

For the dyadic case, i.e. |det __J__<sub>k</sub>|=2 for all factor matrices, this Library also provides the construction of corresponding wavelets and an algorithm to decompose on several different factorizations at the same time.

Further, for the dyadic two-dimensional case, several visualization methods are given for the pattern, the wavelet and scaling functions and the obtained fractions of a function sampled on a pattern and decomposed with respect to the wavelets. 

## Dependencies

The box spline evaluation was written by Leif Kobbelt and is available at [http://www.netlib.org/numeralgo/na11](http://www.netlib.org/numeralgo/na11). For a copyright note to that, see `/helpers/READ.ME`

## License
    MPAWL is free software : you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    MPAWL is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
  
    You should have received a copy of the GNU General Public License
    along with the MPAWL. If not, see <http://www.gnu.org/licenses/>.
