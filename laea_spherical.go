package crs

import (
	"math"
)

type LambertAzimuthalEqualAreaSpherical struct {
	Lonf, Latf, Eastf, Northf float64
}

func (cs LambertAzimuthalEqualAreaSpherical) String() string {
	return build("lambert_azimuthal_equal_area_spherical").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs LambertAzimuthalEqualAreaSpherical) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	r := s.A
	phi0 := radian(cs.Latf)
	lam0 := radian(cs.Lonf)
	phi := radian(lat)
	lam := radian(lon)
	sinPhi0, cosPhi0 := math.Sincos(phi0)
	sinPhi, cosPhi := math.Sincos(phi)
	dLam := lam - lam0
	denom := 1 + sinPhi0*sinPhi + cosPhi0*cosPhi*math.Cos(dLam)

	if denom < 1e-14 {
		return cs.Eastf, cs.Northf, h
	}

	k := math.Sqrt(2 / denom)
	east := cs.Eastf + r*k*cosPhi*math.Sin(dLam)
	north := cs.Northf + r*k*(cosPhi0*sinPhi-sinPhi0*cosPhi*math.Cos(dLam))

	return east, north, h
}

func (cs LambertAzimuthalEqualAreaSpherical) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	r := s.A
	phi0 := radian(cs.Latf)
	lam0 := radian(cs.Lonf)
	x := (east - cs.Eastf) / r
	y := (north - cs.Northf) / r
	rho := math.Hypot(x, y)

	if rho < 1e-14 {
		return cs.Lonf, cs.Latf, h
	}

	c := 2 * math.Asin(clamp(rho/2, -1, 1))
	sinC, cosC := math.Sincos(c)
	sinPhi0, cosPhi0 := math.Sincos(phi0)
	phi := math.Asin(clamp(cosC*sinPhi0+y*sinC*cosPhi0/rho, -1, 1))
	lam := lam0 + math.Atan2(x*sinC, rho*cosPhi0*cosC-y*sinPhi0*sinC)

	return degree(lam), degree(phi), h
}
