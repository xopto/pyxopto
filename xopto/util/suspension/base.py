from typing import Callable, Tuple
import io

from xopto import pf
from xopto.materials import ri
from xopto.materials import density
from xopto.mcbase.mcutil import cache

import xopto.mcbase.mcpf

import numpy as np

DEFAULT_TEMPERATURE = 293.15


class Suspension:
    def __init__(
            self,
            pd: Callable[[float], float] or 'Suspension',
            particle_ri: float or Callable[[float, float], float] = None,
            medium_ri: float or Callable[[float, float], float] = None,
            particle_density: float or Callable[[float], float] = None,
            medium_density: float or Callable[[float], float] = None,
            particle_mua: float or Callable[[float, float], float] = 0.0,
            medium_mua: float or Callable[[float, float], float] = 0.0,
            solid_content: float = 10.0, nd: int or None = 100):
        '''
        Suspension of spherical particles that follows the provided size
        distribution function.

        pd: Callable[[float], float] or Suspension
            Particle size distribution that takes diameter and returns
            probability.
            This object must be hashable and implement several methods.
            The object must implement raw_moment method
            that computes the requested statistical moment of the distribution.
            The object must implement range method that returns a
            tuple with full range of size/diameter covered by the distribution
            as (dmin, dmax). 
            The object must implement todict/fromdict methods. Fromdict
            exports the object to a dict from which it can be constructed by
            calling static method fromdict.
            Use one of the distributions implemented in
            :py:mod:`xopto.util.distribution` or derive a custom distribution
            by following these implementations.
            If this parameter is a Suspension instance, a new copy is
            created, which inherits all the properties of the original
            suspension including the cache.
        particle_ri: float or Callable[[float, float], float]
            Refractive index of the suspended particles.
            A callable that takes wavelength and temperature or a fixed
            floating-point value that is independent of wavelength and
            temperature.
            Default implementation of the polystyrene refractive index
            is used by default.
        medium_ri: float or Callable[[float, float], float]
            Refractive index of the medium surrounding the suspended particles.
            A callable that takes wavelength and temperature or a fixed
            floating-point value that is independent of wavelength and
            temperature.
            Default implementation of the water refractive index
            is used by default.
        particle_density: float or Callable[[float], float]
            Density of the suspended particles (kg/m3).
            A callable that takes temperature or a fixed
            floating-point value that is independent of temperature.
            Default implementation of the polystyrene density
            is used by default.
        medium_density: float or Callable[[float], float]
            Density of the medium surrounding the particles (kg/m3).
            A callable that takes temperature or a fixed
            floating-point value that is independent of temperature.
            Default implementation of the water density
            is used by default.
        particle_mua: float or Callable[[float, float], float]:
            Absorption coefficient (1/m) of the suspended particles.
            A callable that takes wavelength and temperature or a fixed
            floating-point value that is independent of wavelength and
            temperature.
            Default implementation is 0.0 1/m.
        medium_mua: float or Callable[[float, float], float]:
            Absorption coefficient (1/m) of the medium surrounding
            the suspended particles.
            A callable that takes wavelength and temperature or a fixed
            floating-point value that is independent of wavelength and
            temperature.
            Default implementation is 0.0 1/m.
        solid_content: float
            Solid content of the suspension given as % wt/v
            (1% wt/v equals 1 g/100 ml, i.e. 1 g of particles per 100 ml
            of suspension). Defaults to 10% wt/v.
        nd: int or None
            If an integer, the number of nodes to use when integrating
            the Mie scattering phase function over the range of the size
            distribution `pd`.
            If None, adaptive step integration is used (slow).

        Note
        ----
        The scattering properties of the suspension can be adjusted through
        the :py:meth:`~Suspension.set_mus`, :py:meth:`~Suspension.set_musr`,
        :py:meth:`~Suspension.set_number_density` or
        :py:meth:`~Suspension.set_solid_content` methods.
        This class uses :py:class:`~xopto.pf.miepd.MiePd` scattering phase
        function to make the required computations.

        If the absorption coefficients of the medium or particles are nonzero,
        the corresponding refractive index that is used to compute the
        scattering phase function and related properties becomes complex
        :math:`n + ik`. Imaginary part :math:`k` is related to the absorption
        coefficient :math:`\\mu_{a}` as
        :math:`\\mu_{a} = 4 \\pi k / \\lambda_0`, where :math:`\\lambda_0`
        is the wavelength of light in vacuum.
        '''
        if isinstance(pd, Suspension):
            obj = pd
            pd = obj._pd
            particle_ri = obj._particle_ri
            medium_ri = obj._medium_ri
            particle_density = obj._particle_density
            particle_mua = obj._particle_mua
            medium_mua = obj._medium_mua
            solid_content = obj.solid_content(293.15)
            nd = obj._nd
            self._pf_cache = obj._pf_cache
            self._mcpf_lut_cache = obj._mcpf_lut_cache
        else:
            self._pf_cache = cache.ObjCache()
            self._mcpf_lut_cache = cache.LutCache()

        if nd is not None:
            nd = int(nd)

        self._number_density = None
        self._nd = nd

        if particle_ri is None:
            particle_ri = ri.polystyrene.default

        if medium_ri is None:
            medium_ri = ri.water.default

        if particle_density is None:
            particle_density = density.polystyrene.default

        if medium_density is None:
            medium_density = density.water.default

        if isinstance(particle_ri, (float, int)):
            particle_ri_value = float(particle_ri)
            particle_ri = lambda w, t: particle_ri_value

        if isinstance(medium_ri, (float, int)):
            medium_ri_value = float(medium_ri)
            medium_ri = lambda w, t: medium_ri_value

        if isinstance(particle_mua, (float, int)):
            particle_mua_value = float(particle_mua)
            particle_mua = lambda w, t: particle_mua_value

        if isinstance(medium_mua, (float, int)):
            medium_mua_value = float(medium_mua)
            medium_mua = lambda w, t: medium_mua_value

        if isinstance(particle_density, (float, int)):
            particle_density_value = float(particle_density)
            particle_density = lambda t: particle_density_value

        if isinstance(medium_density, (float, int)):
            medium_density_value = float(medium_density)
            medium_density = lambda t: medium_density_value

        self._particle_ri = particle_ri
        self._medium_ri = medium_ri
        self._particle_mua = particle_mua
        self._medium_mua = medium_mua
        self._particle_density = particle_density
        self._medium_density = medium_density
        self._pd = pd

        # this will initialize self._number_density
        self.set_solid_content(solid_content, 293.15)

    def set_musr(self, musr: float, wavelength: float,
                 temperature: float = 293.15):
        '''
        Set the scattering properties of the suspension by specifying the
        reduced scattering coefficient (1/m) at the given wavelength (m).

        Parameters
        ----------
        musr: float
            Reduced scattering coefficient (1/m) at the specified wavelength of
            light.
        wavelength: float
            Wavelength of light (m) at which reduced scattering coefficient is
            specified.
        temperature: float
            Suspension temperature (K).
        '''
        pf_obj = self.pf(wavelength, temperature)
        g, scs = pf_obj.g(1), pf_obj.scs()
        self._number_density = musr/(1.0 - g)/scs

    def set_mus(self, mus: float, wavelength: float,
                temperature: float = 293.15):
        '''
        Set the scattering properties of the suspension by specifying the
        scattering coefficient (1/m) at the given wavelength (m).

        Parameters
        ----------
        mus: float
            Scattering coefficient (1/m) at the specified wavelength of light.
        wavelength: float
            Wavelength of light (m) at which reduced scattering coefficient is
            specified.
        temperature: float
            Suspension temperature (K).
        '''
        self._number_density = mus/self.pf(wavelength, temperature).scs()

    def set_number_density(self, nd: float):
        '''
        Set the scattering properties of the suspension by specifying the
        number density of the particles.

        Parameters
        ----------
        nd: float
            Number density of particles (per cubic metre).
        '''
        self._number_density = nd

    def set_solid_content(self, sc: float, temperature: float = 293.15):
        '''
        Set the scattering properties of the suspension by specifying the
        solid content of the particles in  the suspension.

        Parameters
        ----------
        sc: float
            Solid content of the suspension given as % wt/v
            (1% wt/v equals 1 g/100 ml, i.e. 1 g of particles per 100 ml
            of suspension).
        temperature: float
            Suspension temperature (K).

        Note
        -----
        Solid content of 1% wt/v equals 1 g/100 ml, which equals 10 g/l or
        10 kg/m3. 
        '''
        average_particle_volume = np.pi/6.0*self._pd.raw_moment(3)
        average_particle_weight = \
            average_particle_volume*self._particle_density(temperature)
        # 1 % g/ml ~ 0.01 g/ml ~ 10 g/l ~ 10 kg/m3
        self._number_density = sc*10.0/average_particle_weight

    def mua(self, wavelength: float, temperature: float = 293.15) -> float:
        '''
        Computes the absorption coefficient of the suspension (contributions
        from both the medium and particles) at the given wavelength and
        temperature.

        Parameters
        ----------
        wavelength: float
            Wavelength of light (m).
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        mua: float
            Absorption coefficient (1/m) of the suspension
            at the given wavelength and temperature.

        Note
        ----
        The absorption coefficient of the suspension depends on the absorption
        coefficients of the medium and absorption coefficient of the particles.
        The absorption coefficient of the suspension is computed as the
        sum of the effective absorption coefficient of the medium
        and the effective absorption coefficient of the particles.
        The effective absorption coefficient of the medium is computed as
        the absorption coefficient of the medium weighted by the
        volume fraction of the medium.
        The effective absorption coefficient of the particles is computed
        using the Mie solution for the absorption cross section multiplied
        by the number density of the particles in the suspension.
        '''
        # solid content 1% wt/v equals 1 g/ 100 ml or 10 g/l or 10 kg/m3
        particle_volume_fraction = self.particle_volume_fraction(temperature)

        # absorption coefficient of a pure/bulk medium
        medium_mua = self.medium_mua(wavelength, temperature)

        # absorption coefficient of bulk particles (for comparison)
        #particle_mua = self.particle_mua(wavelength, temperature)

        # computing the absorption cross section and from there the effective
        # absorption coefficient of the suspended particles
        pf_obj = self.pf(wavelength, temperature)
        mua_particle_effective = pf_obj.acs()*self.number_density()

        return  mua_particle_effective + \
            (1.0 - particle_volume_fraction)*medium_mua

    def medium_mua(self, wavelength: float,
                   temperature: float = 293.15) -> float:
        '''
        Computes the absorption coefficient of the suspension medium at
        the given wavelength and temperature.

        Parameters
        ----------
        wavelength: float
            Wavelength of light (m).
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        mua: float
            Absorption coefficient (1/m) of the suspension medium
            at the given wavelength and temperature.
        '''
        return self._medium_mua(wavelength, temperature)

    def particle_mua(self, wavelength: float,
                     temperature: float = 293.15) -> float:
        '''
        Computes the absorption coefficient of the suspended particles at
        the given wavelength and temperature.

        Parameters
        ----------
        wavelength: float
            Wavelength of light (m).
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        mua: float
            Absorption coefficient (1/m) of the suspended particles
            at the given wavelength and temperature.
        '''
        return self._particle_mua(wavelength, temperature)

    def mus(self, wavelength: float, temperature: float = 293.15) -> float:
        '''
        Computes the scattering coefficient of the suspension at
        the given wavelength and temperature.

        Parameters
        ----------
        wavelength: float
            Wavelength of light (m).
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        musr: float
            Scattering coefficient (1/m) of the suspension
            at the given wavelength and temperature.
        '''
        pf_obj = self.pf(wavelength, temperature)
        return self.number_density()*pf_obj.scs()

    def musr(self, wavelength: float, temperature: float = 293.15) -> float:
        '''
        Computes the reduced scattering coefficient of the suspension at
        the given wavelength and temperature.

        Parameters
        ----------
        wavelength: float
            Wavelength of light (m).
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        musr: float
            Reduced scattering coefficient (1/m) of the suspension
            at the given wavelength and temperature.
        '''
        g = self.g(wavelength, temperature)
        return self.mus(wavelength, temperature)*(1.0 - g)

    def particle_volume_fraction(self, temperature: float = 293.15) -> float:
        '''
        Computes and returns the volume fraction of particles in the suspension.
        (a value from 0.0, to 1.0).

        Parameters
        ----------
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        fraction: float
            Volume fraction of particles in the suspension, i.e. total
            volume of particles divided by the total volume of suspension.

        Note
        ----
        The volume fractions of particles and medium returned by
        the :py:meth:`~Suspension.particle_volume_fraction` and
        :py:meth:`Suspension.medium_volume_fraction` sum up to 1 if computed
        at the same temperature.
        '''
        return self.solid_content(temperature)*10.0/ \
            self._particle_density(temperature)

    def medium_volume_fraction(self, temperature: float = 293.15) -> float:
        '''
        Computes and returns the volume fraction of medium in the suspension.
        (a value from 0.0, to 1.0).

        Parameters
        ----------
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        fraction: float
            Volume fraction of medium in the suspension, i.e. total
            volume of medium divided by the total volume of suspension.

        Note
        ----
        The volume fractions of particles and medium returned by
        the :py:meth:`~Suspension.particle_volume_fraction` and
        :py:meth:`Suspension.medium_volume_fraction` sum up to 1 if computed
        at the same temperature.
        '''
        return 1.0 - self.particle_volume_fraction(temperature)


    def number_density(self) -> float:
        '''
        Returns
        -------
        nd: float
            Number density of the particles (1/m3) in the suspension.
        '''
        return self._number_density

    def solid_content(self, temperature: float = 293.15) -> float:
        '''
        Computes and returns the solid content of suspension at the given
        temperature (K)

        Parameters
        ----------
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        sc: float
            Solid content of the suspension given as % wt/v
            (1% wt/v equals 1 g/100 ml, i.e 1 g of particles per 100 ml
            of suspension).

        Note
        -----
        Solid content of 1% wt/v equals 1 g/100 ml, which equals 10 g/l or
        10 kg/m3.
        '''
        average_particle_volume = np.pi/6.0*self._pd.raw_moment(3)
        average_particle_weight = \
            average_particle_volume*self._particle_density(temperature)
        # kg/m3 ~ g/l ~ 0.001 g/ml ~ 0.1 % g/ml
        return self.number_density()*average_particle_weight*0.1

    def pf(self, wavelength: float, temperature: float = 293.15) -> pf.MiePd:
        '''
        Computes and returns an instance of the scattering phase function
        that describes the scattering of suspended particles at the given
        wavelength and temperature.

        Parameters
        -----------
        wavelength: float
            Wavelength of light.
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        pf_obj: pf.MiePd
            Scattering phase function instance at the given wavelength
            and temperature.
        '''
        nd = int(self._nd) if self._nd is not None else None
        return self._pf_cache.get(
            pf.MiePd,
            float(self._particle_ri(wavelength, temperature)) +
                1j*float(self._particle_mua(wavelength, temperature))*
                    wavelength/(4.0*np.pi),
            float(self._medium_ri(wavelength, temperature)) +
                1j*float(self._medium_mua(wavelength, temperature))*
                    wavelength/(4.0*np.pi),
            float(wavelength),
            self._pd, self._pd.range, nd
        )

    def mcpf(self, wavelength: float, temperature: float = 293.15) \
            -> xopto.mcbase.mcpf.PfBase:
        '''
        Returns a Monte Carlo simulator-compatible scattering phase function
        instance at the given wavelength and temperature.

        Parameters
        -----------
        wavelength: float
            Wavelength of light (m).
        temperature: float
            Suspension temperature (K).

        Returns
        -------
            Monte Carlo simulator-compatible scattering phase function
            instance at the given wavelength and temperature.
        '''
        nd = int(self._nd) if self._nd is not None else None
        return self._mcpf_lut_cache.get(
            pf.MiePd,
            (
                float(self._particle_ri(wavelength, temperature)) + 
                    1j*float(self._particle_mua(wavelength, temperature))*
                        wavelength/(4.0*np.pi),
                float(self._medium_ri(wavelength, temperature)) +
                    1j*float(self._medium_mua(wavelength, temperature))*
                        wavelength/(4.0*np.pi),
                float(wavelength),
                self._pd, self._pd.range, nd
            )
        )

    def scs(self, wavelength: float, temperature: float = 293.15) -> float:
        '''
        Computes and returns the scattering cross section of the suspended
        particles at the given wavelength and temperature.

        Parameters
        -----------
        wavelength: float
            Wavelength of light.
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        scs: float
            Scattering cross section (m2) of the suspended particles
            at the given wavelength and temperature.
        '''
        return self.pf(wavelength, temperature).scs()

    def g(self, wavelength: float, temperature: float = 293.15, n: int = 1):
        '''
        Computes and returns the specified Legendre moment of the scattering
        phase function that describes the scattering of the suspended particles
        at the given wavelength and temperature.

        Parameters
        -----------
        wavelength: float
            Wavelength of light.
        temperature: float
            Suspension temperature (K).
        n: int
            Legendre moment (defaults to 1 - anisotropy of the scattering
            phase function).

        Returns
        -------
        scs: float
            The specified Legendre moment of the scattering
            phase function that describes the scattering of
            the suspended particles at the given wavelength and temperature.
        '''
        return self.pf(wavelength, temperature).g(n)

    def particle_ri(self, wavelength: float,
                    temperature: float = 293.15) -> float:
        '''
        Computes and returns the refractive index of suspended particles at
        the given wavelength and temperature.

        Parameters
        -----------
        wavelength: float
            Wavelength of light.
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        ri: float
            Refractive index of the suspended particles at the given
            wavelength and temperature.
        '''
        return self._particle_ri(wavelength, temperature)

    def medium_ri(self, wavelength: float,
                  temperature: float = 293.15) -> float:
        '''
        Computes and returns the refractive index of medium that surrounds
        the suspended particles at the given wavelength and temperature.

        Parameters
        -----------
        wavelength: float
            Wavelength of light.
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        ri: float
            Refractive index of the medium that surrounds the suspended
            particles at the given wavelength and temperature.
        '''
        return self._medium_ri(wavelength, temperature)

    def ri(self, wavelength: float,
                  temperature: float = 293.15) -> float:
        '''
        Computes and returns the refractive index of suspension that equals
        the refractive index of the medium that surrounds
        the suspended particles at the given wavelength and temperature.

        Parameters
        -----------
        wavelength: float
            Wavelength of light.
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        ri: float
            Refractive index of the suspension/medium that surrounds the
            suspended particles at the given wavelength and temperature.
        '''
        return self.medium_ri(wavelength, temperature)

    def particle_density(self, temperature: float = 293.15) -> float:
        '''
        Computes and returns the density (kg/m3) of the suspended particles at
        the given temperature.

        Parameters
        -----------
        wavelength: float
            Wavelength of light.
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        ri: float
            Density (kg/m3) of the suspended particles at the given temperature.
        '''
        return self._particle_density(temperature)

    def medium_density(self, temperature: float = 293.15) -> float:
        '''
        Computes and returns the density (kg/m3) of the medium at the given
        temperature.

        Parameters
        -----------
        wavelength: float
            Wavelength of light.
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        ri: float
            Density (kg/m3) of the medium at the given temperature.
        '''
        return self._medium_density(temperature)

    def make_mus(self, mus: float, volume: float,
                 wavelength: float, temperature: float = 293.15) \
                    -> Tuple[float, 'Suspension']:
        '''
        Compute the volume that needs to be taken from this suspension to
        prepare a diluted suspension with target volume and target scattering
        coefficient at the given wavelength and temperature.

        Parameters
        ----------
        mus: float
            Target scattering coefficient (1/m) of the diluted suspension.
        volume: float
            Volume of the target suspension (m3).
        wavelength: float
            Wavelength of light (m).
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        required_volume: float
            The required volume of this suspension (m3).
        diluted_suspension: Suspension
            Diluted suspension.
        '''
        if mus > self.mus(wavelength, temperature):
            raise ValueError(
                'The scattering coefficient of the diluted '
                'suspension exceeds the value of this suspension!')

        diluted_suspension = Suspension(self)
        diluted_suspension.set_mus(mus, wavelength, temperature)
        # solid content units % wt/v or 1 g/100 ml or 10 g/l or 10 kg/m3
        solid_mass = diluted_suspension.solid_content(temperature)*10.0*volume
        required_volume = solid_mass/(self.solid_content(temperature)*10.0)

        return required_volume, diluted_suspension

    def make_mus_from(self, mus: float, volume: float,
                      wavelength: float, temperature: float = 293.15) \
                        -> Tuple[float, 'Suspension']:
        '''
        Compute the volume of a new suspension that is obtained by
        taking the specified volume of this suspension and diluting it to
        the specified target scattering coefficient at the given wavelength
        and temperature.

        Parameters
        ----------
        mus: float
            Target scattering coefficient (1/m) of the diluted suspension.
        volume: float
            Volume taken from this suspension (m3).
        wavelength: float
            Wavelength of light (m).
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        total_volume: float
            Total volume of the diluted suspension (m3).
        diluted_suspension: Suspension
            Diluted suspension.
        '''
        if mus > self.mus(wavelength, temperature):
            raise ValueError(
                'The scattering coefficient of the diluted '
                'suspension exceeds the value of this suspension!')

        diluted_suspension = Suspension(self)
        diluted_suspension.set_mus(mus, wavelength, temperature)
        # solid content units % wt/v or 1 g/100 ml or 10 g/l or 10 kg/m3
        solid_mass = self.solid_content(temperature)*10.0*volume
        total_volume = \
            solid_mass/(diluted_suspension.solid_content(temperature)*10.0)

        return total_volume, diluted_suspension

    def make_musr(self, musr: float, volume: float,
                  wavelength: float, temperature: float = 293.15) \
                    -> Tuple[float, 'Suspension']:
        '''
        Compute the volume that needs to be taken from this suspension to
        prepare a diluted suspension with target volume and target
        reduced scattering coefficient at the given wavelength and temperature.

        Parameters
        ----------
        musr: float
            Target reduced scattering coefficient (1/m) of the diluted
            suspension.
        volume: float
            Volume of the target suspension (m3).
        wavelength: float
            Wavelength of light (m).
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        required_volume: float
            The required volume of this suspension.
        diluted_suspension: Suspension
            Diluted suspension.
        '''
        if musr > self.musr(wavelength, temperature):
            raise ValueError(
                'The reduced scattering coefficient of the diluted '
                'suspension exceeds the value of this suspension!')

        diluted_suspension = Suspension(self)
        diluted_suspension.set_musr(musr, wavelength, temperature)
        # solid content units % wt/v or 1 g/100 ml or 10 g/l or 10 kg/m3
        solid_mass = diluted_suspension.solid_content(temperature)*10.0*volume
        required_volume = solid_mass/(self.solid_content(temperature)*10.0)

        return required_volume, diluted_suspension

    def make_musr_from(self, musr: float, volume: float,
                       wavelength: float, temperature: float = 293.15) \
                        -> Tuple[float, 'Suspension']:
        '''
        Compute the volume of a new suspension that is obtained by
        taking the specified volume of this suspension and diluting it to
        the specified target reduced scattering coefficient at the given
        wavelength and temperature.

        Parameters
        ----------
        musr: float
            Target reduced scattering coefficient (1/m) of the diluted
            suspension.
        volume: float
            Volume taken from this suspension (m3).
        wavelength: float
            Wavelength of light (m).
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        total_volume: float
            Total volume of the diluted suspension (m3).
        diluted_suspension: Suspension
            Diluted suspension.
        '''
        if musr > self.musr(wavelength, temperature):
            raise ValueError(
                'The reduced scattering coefficient of the diluted '
                'suspension exceeds the value of this suspension!')

        diluted_suspension = Suspension(self)
        diluted_suspension.set_musr(musr, wavelength, temperature)
        # solid content units % wt/v or 1 g/100 ml or 10 g/l or 10 kg/m3
        solid_mass = self.solid_content(temperature)*10.0*volume
        total_volume = \
            solid_mass/(diluted_suspension.solid_content(temperature)*10.0)

        return total_volume, diluted_suspension

    def make_solid_content(self, solid_content: float, volume: float,
                           temperature: float = 293.15) \
                            -> Tuple[float, 'Suspension']:
        '''
        Compute the volume that needs to be taken from this suspension to
        prepare a diluted suspension with target volume and target solid
        content at the given temperature.

        Parameters
        ----------
        solid_content: float
            Target solid content of the diluted suspension given as % wt/v
            (1% wt/v equals 1 g/100 ml, i.e. 1 g of particles per 100 ml
            of suspension). Defaults to 10% wt/v.
        volume: float
            Volume of the target suspension (m3).
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        required_volume: float
            The required volume of this suspension (m3).
        diluted_suspension: Suspension
            Diluted suspension.
        '''
        if solid_content > self.solid_content(temperature):
            raise ValueError(
                'The solid content of the diluted '
                'suspension exceeds the solid content of this suspension!')

        diluted_suspension = Suspension(self)
        diluted_suspension.set_solid_content(solid_content, temperature)
        # solid content units % wt/v or 1 g/100 ml or 10 g/l or 10 kg/m3
        solid_mass = solid_content*10.0*volume
        required_volume = solid_mass/(self.solid_content(temperature)*10.0)

        return required_volume, diluted_suspension

    def make_solid_content_from(self, solid_content: float, volume: float,
                                temperature: float = 293.15) \
                                    -> Tuple[float, 'Suspension']:
        '''
        Compute the volume of a new suspension that is obtained by
        taking the specified volume of this suspension and diluting it to
        the specified target solid content at the given temperature.

        Parameters
        ----------
        solid_content: float
            Target solid content of the diluted suspension given as % wt/v
            (1% wt/v equals 1 g/100 ml, i.e. 1 g of particles per 100 ml
            of suspension). Defaults to 10% wt/v.
        volume: float
            Volume taken from this suspension (m3).
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        total_volume: float
            Total volume of the diluted suspension (m3).
        diluted_suspension: Suspension
            Diluted suspension
        '''
        if solid_content > self.solid_content(temperature):
            raise ValueError(
                'The solid content of the diluted '
                'suspension exceeds the solid content of this suspension!')

        diluted_suspension = Suspension(self)
        diluted_suspension.set_solid_content(solid_content, temperature)
        # solid content units % wt/v or 1 g/100 ml or 10 g/l or 10 kg/m3
        solid_mass = solid_content*10.0*volume
        total_volume = solid_mass/ \
            (diluted_suspension.solid_content(temperature)*10.0)

        return total_volume, diluted_suspension

    def dilute_volume(self, take: float, dilute: float,
                      temperature: float = 293.15) -> 'Suspension':
        '''
        Dilute the given volume of this suspension (m3) to the target
        volume (m3).

        Parameters
        ----------
        take: float
            Volume (m3) of this suspension that will be diluted to the
            specified volume.
        dilute: float
            Diluted volume (m3). Must be greater than the taken volume.
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        suspension: Suspension
            Diluted suspension.
        '''
        if dilute < take:
            raise ValueError('Diluted suspension volume must not be smaller '
                             'than the take volume!')
        # solid content units % wt/v or 1 g/100 ml or 10 g/l or 10 kg/m3
        solid_mass = self.solid_content(temperature)*10.0*take
        diluted_solid_content = solid_mass/(dilute*10.0)

        diluted_suspension = Suspension(self)
        diluted_suspension.set_solid_content(diluted_solid_content, temperature)

        return diluted_suspension

    def dilute_mass(self, take: float, dilute: float,
                    temperature: float = 293.15) -> 'Suspension':
        '''
        Dilute the given mass of this suspension (kg) to the target
        mass (kg).

        Parameters
        ----------
        take: float
            Mass (kg) of this suspension that will be diluted to the
            specified mass.
        dilute: float
            Diluted mass (kg). Must be greater than the taken mass.
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        suspension: Suspension
            Diluted suspension.
        '''
        if dilute < take:
            raise ValueError('Diluted suspension mass must not be smaller '
                             'than the take mass!')
        # solid content units % wt/v or 1 g/100 ml or 10 g/l or 10 kg/m3
        v_take = self.volume(take, temperature)
        solid_mass = self.solid_content(temperature)*10.0*v_take

        v_dilute = solid_mass/self._particle_density(temperature) + \
                  (dilute - solid_mass)/self._medium_density(temperature)

        return self.dilute_volume(v_take, v_dilute)

    def sequential_dilution_recipe(
            self,
            dilutions: Tuple['Suspension', ...],
            volumes: Tuple[float, ...],
            temperature: float = 293.15) -> Tuple[
                Tuple[Tuple[float, float], ...],
                Tuple[Tuple[float, float], ...],
                Tuple[int, ...]
            ]:
        '''
        Create a sequential dilution recipe for the given tuple of diluted
        suspensions.

        Parameters
        ----------
        dilutions: Tuple[Suspension, ...]
            Diluted suspensions. The suspensions can be listed in an arbitrary
            order, but the solid content of any diluted suspensions must not
            exceed the solid content of this suspension.
            The length of the dilutions and volumes tuple must be the same.
        volumes: Tuple[float, ..]
            Required volumes (m3) of the diluted suspensions.
            The length of the dilutions and volumes tuple must be the same.
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        volume_recipe: Tuple[Tuple[float, float], ...]
            Dilution recipe as a tuple of volumes ((take, dilute), ...) (m3),
            where take is the volume (m3) to take from the previous suspension
            in the sequence (starting with this suspension) and dilution
            the volume (m3) up to which the taken suspension sample should be
            diluted.
        mass_recipe: Tuple[Tuple[float, float], ...]
            Dilution recipe as a tuple of masses ((take, dilute), ...) (kg),
            where take is the mass (kg) to take from the previous suspension
            in the sequence (starting with this suspension) and dilution
            the mass (kg) up to which the taken suspension sample should be
            diluted.
        argsort: Tuple[ind, ...]
            Sorting indices into the dilutions argument tuple that produces an
            ordered sequence of diluted suspensions matching
            the volume_recipe and mass_recipe arguments.

        Note
        ----
        Dilution starts with this suspension.
        '''
        sortedDilutions = list(zip(dilutions, volumes))
        sortedDilutions.sort(key=lambda x: x[0].solid_content(temperature))
        volumeRecipe = []
        massRecipe = []
        diluted, v_diluted = sortedDilutions[0]
        suspensionSequence = sortedDilutions[1:]
        suspensionSequence.append((self, 0.0))
        for fromSuspension, volume in suspensionSequence:
            v_take, _ = fromSuspension.make_solid_content(
                diluted.solid_content(), v_diluted, temperature=temperature)
            m_take = fromSuspension.mass(v_take, temperature=temperature)
            m_diluted = diluted.mass(v_diluted, temperature=temperature)
            volumeRecipe.append((v_take, v_diluted))
            massRecipe.append((m_take, m_diluted))
            diluted, v_diluted = fromSuspension, volume + v_take

        return tuple(volumeRecipe[::-1]), tuple(massRecipe[::-1]), \
               tuple(dilutions.index(s[0]) for s in sortedDilutions[::-1])

    def apply_sequential_volume_dilution_recipe(
            self, recipe: Tuple[Tuple[float, float], ...],
            temperature: float = 293.15) \
                -> Tuple[Tuple['Suspension', ...], Tuple[float, ...]]:
        '''
        Apply a sequential volume dilution recipe.

        Parameters
        -------
        recipe: Tuple[Tuple[float, float], ...]
            Dilution recipe as a tuple of volumes ((take, dilute), ...) (m3),
            where take is the volume (m3) to take from the previous suspension
            in the sequence (starting with this suspension) and dilution
            the volume (m3) up to which the taken suspension sample should be
            diluted.
        temperature: float
            Suspension temperature (K).

        Returns
        ----------
        dilutions: Tuple[Suspension, ...]
            Diluted suspensions produced by the recipe.
        volumes: Tuple[float, ..]
            Final available volumes (m3) of each diluted suspensions.

        Note
        ----
        Dilution starts with this suspension.
        '''
        dilutions = []
        volumes = []
        suspension = self
        for v_take, v_dilute in recipe:
            suspension = suspension.dilute_volume(
                    v_take, v_dilute, temperature=temperature)
            dilutions.append(suspension)
            volumes.append(v_dilute)

        for index, (r1, r2) in enumerate(zip(recipe[:-1], recipe[1:])):
            _, v_dilute1 = r1
            v_take2, _ = r2
            volumes[index] = v_dilute1 - v_take2

        return tuple(dilutions), tuple(volumes)

    def apply_sequential_mass_dilution_recipe(
            self, recipe: Tuple[Tuple[float, float], ...],
            temperature: float = 293.15) \
                -> Tuple[Tuple['Suspension', ...], Tuple[float, ...]]:
        '''
        Apply a sequential volume dilution recipe.

        Parameters
        -------
        recipe: Tuple[Tuple[float, float], ...]
            Dilution recipe as a tuple of masses ((take, dilute), ...) (kg),
            where take is the mass (kg) to take from the previous suspension
            in the sequence (starting with this suspension) and dilution
            the mass (kg) up to which the taken suspension sample should be
            diluted.
        temperature: float
            Suspension temperature (K).

        Returns
        ----------
        dilutions: Tuple[Suspension, ...]
            Diluted suspensions produced by the recipe.
        masses: Tuple[float, ..]
            Final available masses (m3) of each diluted suspensions.

        Note
        ----
        Dilution starts with this suspension.
        '''
        dilutions = []
        masses = []
        suspension = self
        for m_take, m_dilute in recipe:
            suspension = suspension.dilute_mass(
                    m_take, m_dilute, temperature=temperature)
            dilutions.append(suspension)
            masses.append(m_dilute)

        for index, (r1, r2) in enumerate(zip(recipe[:-1], recipe[1:])):
            _, m_dilute1 = r1
            m_take2, _ = r2
            masses[index] = m_dilute1 - m_take2

        return tuple(dilutions), tuple(masses)

    def mass(self, volume: float, temperature: float = 293.15) -> float:
        '''
        Compute mass (kg) of the given suspension volume (m3).

        Parameters
        ----------
        volume: float
            Suspension volume (m3).
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        mass: float
            Mass (kg) of the given suspension volume,
        '''
        k_medium = self.medium_volume_fraction(temperature)

        return volume*k_medium*self._medium_density(temperature) + \
               volume*(1.0 - k_medium)*self._particle_density(temperature)

    def volume(self, mass: float, temperature: float = 293.15) -> float:
        '''
        Compute volume (m3) of the given suspension mass (kg).

        Parameters
        ----------
        mass: float
            Suspension volume (m3).
        temperature: float
            Suspension temperature (K).

        Returns
        -------
        volume: float
            Volume (m3) of the given suspension mass,
        '''
        k_medium = self.medium_volume_fraction(temperature)

        return mass/(k_medium*self._medium_density(temperature) + 
                     (1.0 - k_medium)*self._particle_density(temperature))

    def save_cache(self, filename: str or io.IOBase):
        '''
        Save scattering phase function and corresponding MC lookup table
        cache to a file.

        Parameters
        ----------
        filename: str or io.IOBase
            Target filename or a file like object. If filename is a string,
            ".pkl" extension is added automatically. If filename is a file
            like object, it must support binary write operation.
        '''
        if isinstance(filename, str):
            if not filename.endswith(('.pkl', '.PKL')):
                filename = filename + '.pkl'

        if isinstance(filename, str):
            with open(filename, 'wb') as fid:
                self._pf_cache.save(fid)
                self._mcpf_lut_cache.save(fid)
        else:
            self._pf_cache.save(fid)
            self._mcpf_lut_cache.save(fid)

    def load_cache(self, filename: str or io.IOBase):
        '''
        Load scattering phase function and corresponding MC lookup table
        cache from a file.

        Parameters
        ----------
        filename: str or io.IOBase
            Source filename or a file like object. If filename is a string,
            ".pkl" extension is added automatically. If filename is a file
            like object, it must support binary read operation.
        '''
        if isinstance(filename, str):
            if not filename.endswith(('.pkl', '.PKL')):
                filename = filename + '.pkl'

        if isinstance(filename, str):
            with open(filename, 'rb') as fid:
                self._pf_cache = cache.ObjCache.load(fid)
                self._mcpf_lut_cache = cache.LutCache.load(fid)
        else:
            self._pf_cache = cache.ObjCache.load(fid)
            self._mcpf_lut_cache = cache.LutCache.load(fid)

    def _get_cache(self) -> Tuple[cache.ObjCache, cache.LutCache]:
        return self._pf_cache, self._mcpf_lut_cache
    cache = property(_get_cache, None, None,
                    'Returns the scattering phase function and '
                    'corresponding MC lookup table cache objects as tuple '
                    '(pf_cache, mcpf_cache).')
