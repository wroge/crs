package crs

import (
	"math"
)

type AzimuthalEquidistant struct {
	Lonf, Latf, Eastf, Northf float64
}

func (cs AzimuthalEquidistant) String() string {
	return build("azimuthal_equidistant").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs AzimuthalEquidistant) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	a := s.A
	phi0 := radian(cs.Latf)
	lam0 := radian(cs.Lonf)
	phi := radian(lat)
	lam := radian(lon)
	sinPhi0, cosPhi0 := math.Sincos(phi0)
	sinPhi, cosPhi := math.Sincos(phi)
	dLam := lam - lam0
	cosC := sinPhi0*sinPhi + cosPhi0*cosPhi*math.Cos(dLam)
	c := math.Acos(clamp(cosC, -1, 1))

	if math.Abs(c) < 1e-14 {
		return cs.Eastf, cs.Northf, h
	}

	k := c / math.Sin(c)
	east := cs.Eastf + a*k*cosPhi*math.Sin(dLam)
	north := cs.Northf + a*k*(cosPhi0*sinPhi-sinPhi0*cosPhi*math.Cos(dLam))

	return east, north, h
}

func (cs AzimuthalEquidistant) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	a := s.A
	phi0 := radian(cs.Latf)
	lam0 := radian(cs.Lonf)
	x := (east - cs.Eastf) / a
	y := (north - cs.Northf) / a
	rho := math.Hypot(x, y)

	if rho < 1e-14 {
		return cs.Lonf, cs.Latf, h
	}

	c := rho
	sinC, cosC := math.Sincos(c)
	sinPhi0, cosPhi0 := math.Sincos(phi0)
	phi := math.Asin(clamp(cosC*sinPhi0+y*sinC*cosPhi0/rho, -1, 1))
	lam := lam0 + math.Atan2(x*sinC, rho*cosPhi0*cosC-y*sinPhi0*sinC)

	return degree(lam), degree(phi), h
}
