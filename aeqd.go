package crs

import "math"

// AzimuthalEquidistant is EPSG method 1125.
// On an ellipsoid the map is the geodesic distance and azimuth from the
// natural origin (Vincenty). A sphere keeps the spherical formula.
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
	if s.F() == 0 {
		return cs.sphereFrom(s, lon, lat, h)
	}

	dist, azi, ok := geodesicInverse(s, cs.Lonf, cs.Latf, lon, lat)
	if !ok || dist < 1e-8 {
		return cs.Eastf, cs.Northf, h
	}

	return cs.Eastf + dist*math.Sin(azi), cs.Northf + dist*math.Cos(azi), h
}

func (cs AzimuthalEquidistant) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	if s.F() == 0 {
		return cs.sphereTo(s, east, north, h)
	}

	de := east - cs.Eastf
	dn := north - cs.Northf
	dist := math.Hypot(de, dn)
	if dist < 1e-8 {
		return cs.Lonf, cs.Latf, h
	}

	lon, lat, ok := geodesicDirect(s, cs.Lonf, cs.Latf, math.Atan2(de, dn), dist)
	if !ok {
		return cs.Lonf, cs.Latf, h
	}
	return lon, lat, h
}

func (cs AzimuthalEquidistant) sphereFrom(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	a := s.A
	phi0 := radian(cs.Latf)
	phi := radian(lat)
	sinPhi0, cosPhi0 := math.Sincos(phi0)
	sinPhi, cosPhi := math.Sincos(phi)
	dLam := deltaLon(lon, cs.Lonf)
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

func (cs AzimuthalEquidistant) sphereTo(s Spheroid, east, north, h float64) (float64, float64, float64) {
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
