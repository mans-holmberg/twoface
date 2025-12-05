from . import _twoface
import numpy as np

deg_to_rad = np.pi / 180.

def quadraticlimbdarkening(r, u, answer):
    mask = (r <= 1.0)
    mu = np.sqrt(1.0 - np.power(r[mask], 2))
    answer[mask] = 1.0 - u[0] * (1. - mu) - u[1] * (1. - mu)**2


def nonlinearlimbdarkening(r, u, answer):
    mask = (r <= 1.0)
    mu = np.sqrt(1.0 - np.power(r[mask], 2))
    answer[mask] = 1.0 - u[0] * (1. - np.sqrt(mu)) - u[1] * (1. - mu) - u[2] * (1. - np.power(mu, 1.5)) - u[3] * (1. - np.power(mu, 2.0))


class TransitParams:
    def __init__(self):
        self.t0 = np.nan  # time of inferior conjunction
        self.per = np.nan  # orbital period
        self.rp = np.nan  # planet radius (in units of stellar radii)
        self.a = np.nan  # semi-major axis (in units of stellar radii)
        self.inc = np.nan  # orbital inclination (in degrees)
        self.ecc = 0.  # eccentricity
        self.w = 0.  # longitude of periastron (in degrees)
        self.u = []  # limb darkening coefficients
        self.limb_dark = np.nan

        self.num_spots = 0
        self.spot_mu = []
        self.spot_phi = []
        self.spot_radius = []
        self.spot_f = []


class TransitModel:
    def __init__(self, params, t, num_rings=1500, limb_dark_function=None, use_spotrod=True):
        self.params = params
        self.t = t
        self.num_ints = t.shape[0]
        self.num_rings = num_rings
        self.num_spots = params.num_spots

        self.r = np.linspace(1.0 / (2 * num_rings), 1.0 - 1.0 / (2 * num_rings), num_rings)
        self.planetangle = np.zeros((self.num_ints, self.num_rings))
        self.model = np.zeros((self.num_ints))

        self.f = np.zeros((self.num_rings))
        self.spotcenterdistance = np.zeros((self.num_spots + self.num_spots * self.num_rings + self.num_rings))
        self.bounds = np.zeros((8 + 2 * self.num_spots))

        self.planetz =  np.zeros((self.num_ints))
        self.planetx =  np.zeros((self.num_ints))
        self.planety =  np.zeros((self.num_ints))
        self.psi =      np.zeros((self.num_ints))

        self.limb_dark_function = limb_dark_function
        self.use_spotrod = use_spotrod


    def light_curve(self, params=None):

        if params is None:
            params = self.params

        t0 = params.t0  # time of inferior conjunction
        per = params.per  # orbital period
        rp = params.rp  # planet radius (in units of stellar radii)
        a = params.a  # semi-major axis (in units of stellar radii)
        inc = params.inc  # orbital inclination (in degrees)
        ecc = params.ecc  # eccentricity
        w = params.w  # longitude of periastron (in degrees)
        u = params.u  # limb darkening coefficients
        limb_dark = params.limb_dark
        rp2 = params.rp2
        phi = params.phi

        if self.num_spots != params.num_spots:
            self.num_spots = params.num_spots
            self.spotcenterdistance = np.zeros((self.num_spots + self.num_spots * self.num_rings + self.num_rings))
            self.bounds = np.zeros((8 + 2 * self.num_spots))

        if self.num_spots == 0:
            spot_mu = np.array([0.0])
            spot_phi = np.array([0.0])
            spot_radius = np.array([0.0])
            spot_f = np.array([1.0])
        else:
            spot_mu = np.array(params.spot_mu)
            spot_phi = np.array(params.spot_phi)
            spot_radius = np.array(params.spot_radius)
            spot_f = np.array(params.spot_f)

        spotx = np.sqrt(1.0 - spot_mu**2) * np.cos(spot_phi)
        spoty = np.sqrt(1.0 - spot_mu**2) * np.sin(spot_phi)
        spotradius = spot_radius
        spotcontrast = spot_f

        if limb_dark == 'quadratic':
            quadraticlimbdarkening(self.r, u, self.f)
        elif limb_dark == 'nonlinear':
            nonlinearlimbdarkening(self.r, u, self.f)
        elif limb_dark == 'custom':
            self.limb_dark_function(self.r, u, self.f)
        else:
            raise Exception('Only supports quadratic, nonlinear, or custom limb-darkening laws.')
        
        _twoface._sky(self.num_ints, self.t, t0, per, a, inc * deg_to_rad, ecc, w * deg_to_rad, self.planetz, self.planetx, self.planety, self.psi)
        self.psi += phi * deg_to_rad - np.pi/2

        if (rp == rp2) and self.use_spotrod:
            _twoface._circleangle(self.r, rp, self.planetz, self.planetangle)
            _twoface._integratetransit(self.num_ints, self.num_rings, self.num_spots, self.planetx, self.planety, self.planetz, rp, self.r, self.f,
                                             spotx, spoty, spotradius, spotcontrast, self.planetangle, self.model)
        else:
            _twoface._integratetransit_asymmetric(self.num_ints, self.num_rings, self.num_spots, self.r, self.f, self.planetx, self.planety, self.planetz, self.psi, rp2, rp, 
                                         spotx, spoty, spotradius, spotcontrast, self.spotcenterdistance, self.bounds, self.model)
        return np.copy(self.model)
