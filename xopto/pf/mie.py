# -*- coding: utf-8 -*-
################################ Begin license #################################
# Copyright (C) Laboratory of Imaging technologies,
#               Faculty of Electrical Engineering,
#               University of Ljubljana.
#
# This file is part of PyXOpto.
#
# PyXOpto is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# PyXOpto is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with PyXOpto. If not, see <https://www.gnu.org/licenses/>.
################################# End license ##################################

import numpy as np
from scipy.special import jv, yv

from .pfbase import PfBase


def _Mie_eff(a, b, x):
    '''
    Computation of Mie Efficiencies for given
    complex refractive-index ratio m=m'+im"
    and size parameter x=k0*a, where k0= wave number in ambient
    medium, a=sphere radius, using complex Mie Coefficients
    an and bn for n=1 to nmax,
    s. Bohren and Huffman (1983) BEWI:TDD122, p. 103,119-122,477.
    Result: m', m", x, efficiencies for extinction (qext),
    scattering (qsca), absorption (qabs), backscattering (qb),
    asymmetry parameter (asy=<costeta>) and (qratio=qb/qsca).
    Uses the function "Mie_ab" for an and bn, for n=1 to nmax.
    C. M�tzler, May 2002, revised July 2002.
    '''
    if x == 0.0: # To avoid a singularity at x=0
        return 0.0, 0.0, 0.0, 0.0, 0.0, 1.5

    elif x > 0.0: # This is the normal situation
        nmax = int(np.round(2.0 + x + 4.0*x**(1.0/3.0)))
        n1 = nmax-1

        n = np.arange(1, nmax + 1)
        cn = 2.0*n + 1.0
        c1n = n*(n + 2.0)/(n + 1.0)
        c2n = cn/n/(n + 1.0)

        x2 = x*x

        # f=Mie_ab(m,x);

        anp, anpp = np.real(a), np.imag(a)

        bnp, bnpp = np.real(b), np.imag(b)

        g1 = np.zeros([4, nmax])    # displaced numbers used for
        g1[0, :n1] = anp[1:nmax]        # asymmetry parameter, p. 120
        g1[1, :n1] = anpp[1:nmax]
        g1[2, :n1] = bnp[1:nmax]
        g1[3, :n1] = bnpp[1:nmax]

        dn = cn*(anp+bnp)
        q = dn.sum()
        qext = 2.0*q/x2

        en = cn*(anp*anp+anpp*anpp+bnp*bnp+bnpp*bnpp)
        q = np.sum(en)
        qsca = 2.0*q/x2

        qabs = qext - qsca

        fn = (a - b)*cn
        gn = (-1)**n
        q = (fn*gn).sum()
        qb = np.dot(q, q)/x2

        asy1 = c1n*(anp*g1[0] + anpp*g1[1] + \
            bnp*g1[2] + bnpp*g1[3])
        asy2 = c2n*(anp*bnp + anpp*bnpp)
        asy = 4.0/x2*(asy1 + asy2).sum()/qsca

        qratio = qb/qsca

        return qext, qsca, qabs, qb, asy, qratio


def _Mie_S12(a, b, x, costheta):

    # Computation of Mie Scattering functions S1 and S2
    # for complex refractive index m=m'+im",
    # size parameter x=k0*a, and u=cos(scattering angle),
    # where k0=vacuum wave number, a=sphere radius;
    # s. p. 111-114, Bohren and Huffman (1983) BEWI:TDD122
    # C. M�tzler, May 2002

    Nmax = int(np.round(2.0 + x + 4.0*x**(1.0/3.0)))

    #ab = Mie_ab(m,x);

    p, t = _Mie_pt(costheta, Nmax)
    pin = p
    tin = t

    n = np.arange(1.0, Nmax + 1.0)
    n2 = (2.0*n + 1.0)/(n*(n + 1.0))
    n2 = np.tile(n2, [costheta.size, 1])

    pin = (n2*pin).transpose()
    tin = (n2*tin).transpose()

    S1 = np.dot(a, pin) + np.dot(b, tin)
    S2 = np.dot(a, tin) + np.dot(b, pin)

    return S1, S2

def _Mie_ab(m, x):
    '''
    Computes a matrix of Mie Coefficients, an, bn,
    of orders n=1 to nmax, for given complex refractive-index
    ratio m=m'+im" and size parameter x=k0*a where k0= wave number in ambient
    medium for spheres of radius a;
    Eq. (4.88) of Bohren and Huffman (1983), BEWI:TDD122
    using the recurrence relation (4.89) for Dn on p. 127 and
    starting conditions as described in Appendix A.
    C. M�tzler, July 2002
    '''

    z = m*x
    nmax = int(np.round(2 + x + 4*x**(1.0/3.0)))
    nmx = int(np.round(max(nmax, np.abs(z)) + 16.0))

    n = np.arange(nmax)
    nu = n + 1.5

    sx = np.sqrt(0.5*np.pi*x)
    px = sx*jv(nu, x)
    p1x = np.hstack([np.sin(x), px[:nmax-1]])

    chx = -sx*yv(nu, x)
    ch1x = np.hstack([np.cos(x), chx[:nmax-1]])
    gsx = px - 1.0j*chx
    gs1x = p1x - 1.0j*ch1x

    dnx = np.zeros([nmx], dtype=np.complex)
    for j in range(nmx, 1, -1):      # Computation of Dn(z) according to (4.89) of B+H (1983)
        dnx[j-2] = j/z - 1.0/(dnx[j-1] + j/z)


    dn = dnx[n]          # Dn(z), n=1 to nmax
    da = dn/m + (n + 1)/x
    db = m*dn + (n + 1)/x

    an = (da*px - p1x)/(da*gsx - gs1x)
    bn = (db*px - p1x)/(db*gsx - gs1x)

    return an, bn


def _Mie_pt(costheta, Nmax):
    # pi_n and tau_n, -1 <= u <= 1, n1 integer from 1 to Nmax
    # angular functions used in Mie Theory
    # Bohren and Huffman (1983), p. 94 - 95

    Nmax = int(Nmax)
    n = costheta.size
    p = np.zeros([n, Nmax])
    t = np.zeros([n, Nmax])

    p[:, 0] = np.ones([n])
    t[:, 0] = costheta

    p[:, 1] = 3.0*costheta
    t[:, 1] = 3.0*np.cos(2.0*np.arccos(costheta))

    for n1 in range(3, Nmax + 1):

        p1 = (2*n1 - 1.0)/(n1 - 1.0)*p[:, n1 - 2]*costheta
        p2 = n1/(n1 - 1.0)*p[:, n1 - 3]
        p[:, n1 - 1] = p1 - p2

        t1 = n1*costheta*p[:, n1 - 1]
        t2 = (n1 + 1.0)*p[:, n1 - 2]
        t[:, n1 - 1] = t1 - t2

    return p, t


class Mie(PfBase):
    def __init__(self, nsphere: float or complex, nmedium: float or complex,
                 diameter: float, wavelength: float):
        '''
        Mie scattering phase function of spherical particles.

        Parameters
        ----------
        nsphere:  float or complex
            Refractive index (complex) of the spherical paricle.
        nmedium: float or complex
            Refractive index (complex) of the surrounding medium.
        diameter: float
            Diameter of the spherical particle (m).
        wavelength: float
            Wavelength of light (m).

        Examples
        --------
        Mie scattering phase function for a 1 um spherical particle
        (refractive index 1.6) at 550 nm in water (refractive index 1.33).

        >>> import numpy as np
        >>> from matplotlib import pyplot as pp
        >>>
        >>> cos_theta = np.linspace(-1.0, 1.0, 1000)
        >>>
        >>> pp.figure()
        >>> for d in [0.1, 0.5, 1, 2, 5]:
        >>>   pf = Mie(1.6, 1.33, d*1e-6, 550e-9)
        >>>   pp.semilogy(cos_theta, pf(cos_theta), label='Monodisperse({} um)'.format(d))
        >>> pp.legend()
        '''
        self._nmedium = nmedium
        self._nsphere = nsphere
        self._wavelength = wavelength
        self._diameter = diameter

        super().__init__()

        self._m = nsphere/nmedium
        # light wavelength in medium
        wavelength_medium = wavelength/complex(nmedium).real
        # size parameter
        self._x = np.pi*diameter/wavelength_medium

        self._a, self._b = _Mie_ab(self._m, self._x)
        self._result = _Mie_eff(self._a, self._b, self._x)
        self._qext = self._result[0]
        self._qsca = self._result[1]
        self._qabs = self._result[2]
        self._g1 = self._result[4]
        s = np.pi*(0.5*diameter)**2
        self._sigma_scs = self._qsca*s
        self._sigma_ecs = self._qext*s

    def g(self, n: int, **kwargs) -> float:
        '''
        Overloads the :py:meth:`PfBase.g` method of the base class.
        Computes the n-th Legendre moment of the scattering phase function.

        Note
        ----
            If n is 1, a precalculated g1 is returned.
        '''
        if n == 1:
            return self._g1
        else:
            return PfBase.g(self, n, **kwargs)

    def fastg(self, n: int, *args, **kwargs) -> float:
        '''
        Overloads the :py:meth:`PfBase.fastg` method of the base class.

        Note
        ----
            If n is 1, a precalculated g1 is returned.
        '''
        if n == 1:
            return self._g1
        else:
            return PfBase.fastg(self, n, *args, **kwargs)

    def scs(self) -> float:
        '''
        Returns the scattering cross section.
        '''
        return self._sigma_scs

    def ecs(self) -> float:
        '''
        Returns the extinction cross section.
        '''
        return self._sigma_ecs

    def acs(self) -> float:
        '''
        Returns the absorption cross section.
        '''
        return self._sigma_ecs - self._sigma_scs

    def __call__(self, costheta: float or np.ndarray, scsmul: bool = False) \
                -> float or np.ndarray:
        '''
        Call method of the Mie scattering phase function.

        Parameters
        ----------
        costheta: float or np.ndarray
            Scattering angle cosines at which the scattering phase function is
            evaluated.
        scsmul: bool
            If nonzero, the scattering phase function is multiplied by the
            scattering cross section.

        Returns
        -------
        pf: float or np.ndarray
            Scattering phase function at the specified scattering angle cosines.
        '''
        S1, S2 = _Mie_S12(self._a, self._b, self._x, np.asarray(costheta))
        pf = 2*np.pi*(np.abs(S1)**2+np.abs(S2)**2)/ \
            (2*np.pi*self._x**2.0*self._qsca)

        if scsmul:
            return pf*self._sigma_scs
        else:
            return pf

    def __repr__(self):
        return 'Mie(nsphere={}, nmedium={}, diameter={}, '\
                    'wavelength={})'.format(self._nsphere, self._nmedium,
                                            self._diameter, self._wavelength)
