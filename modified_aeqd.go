package crs

import (
	"math"
)

type ModifiedAzimuthalEquidistant struct {
	Lonf, Latf, Eastf, Northf float64
}

func (cs ModifiedAzimuthalEquidistant) String() string {
	return build("modified_azimuthal_equidistant").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs ModifiedAzimuthalEquidistant) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	a, e, e2 := s.A, s.E(), s.E2()
	phi0 := radian(cs.Latf)
	lam0 := radian(cs.Lonf)
	phi := radian(lat)
	dLam := radian(lon) - lam0
	sinPhi0, cosPhi0 := math.Sincos(phi0)
	sinPhi, cosPhi := math.Sincos(phi)
	nu0 := a / math.Sqrt(1-e2*sinPhi0*sinPhi0)
	nu := a / math.Sqrt(1-e2*sinPhi*sinPhi)
	psi := math.Atan((1-e2)*math.Tan(phi) + e2*nu0*sinPhi0/(nu*cosPhi))
	alpha := math.Atan2(math.Sin(dLam), cosPhi0*math.Tan(psi)-sinPhi0*math.Cos(dLam))
	sinAlpha, cosAlpha := math.Sincos(alpha)
	g := e * sinPhi0 / math.Sqrt(1-e2)
	hCoeff := e * cosPhi0 * cosAlpha / math.Sqrt(1-e2)

	var sAng float64
	if math.Abs(sinAlpha) < 1e-14 {
		sAng = math.Asin(clamp(cosPhi0*math.Sin(psi)-sinPhi0*math.Cos(psi), -1, 1))
		if cosAlpha < 0 {
			sAng = -sAng
		}
	} else {
		sAng = math.Asin(clamp(math.Sin(dLam)*math.Cos(psi)/sinAlpha, -1, 1))
	}

	s2 := sAng * sAng
	s3 := s2 * sAng
	s4 := s2 * s2
	s5 := s4 * sAng
	h2 := hCoeff * hCoeff
	c := nu0 * sAng * (1 - s2*h2*(1-h2)/6 + s3/8*g*hCoeff*(1-2*h2) +
		s4/120*(h2*(4-7*h2)-3*g*g*(1-7*h2)) - s5/48*g*hCoeff)

	return cs.Eastf + c*sinAlpha, cs.Northf + c*cosAlpha, h
}

func (cs ModifiedAzimuthalEquidistant) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	a, e2 := s.A, s.E2()
	phi0 := radian(cs.Latf)
	sinPhi0, cosPhi0 := math.Sincos(phi0)
	nu0 := a / math.Sqrt(1-e2*sinPhi0*sinPhi0)
	x := east - cs.Eastf
	y := north - cs.Northf
	c := math.Hypot(x, y)

	if c < 1e-14 {
		return cs.Lonf, cs.Latf, h
	}

	alpha := math.Atan2(x, y)
	sinAlpha, cosAlpha := math.Sincos(alpha)
	aa := -e2 * cosPhi0 * cosPhi0 * cosAlpha * cosAlpha / (1 - e2)
	b := 3 * e2 * (1 - aa) * sinPhi0 * cosPhi0 * cosAlpha / (1 - e2)
	d := c / nu0
	j := d - aa*(1+aa)*d*d*d/6 - b*(1+3*aa)*d*d*d*d/24
	k := 1 - aa*j*j/2 - b*j*j*j/6
	psi := math.Asin(clamp(sinPhi0*math.Cos(j)+cosPhi0*math.Sin(j)*cosAlpha, -1, 1))
	phi := math.Atan((1 - e2*k*sinPhi0/math.Sin(psi)) * math.Tan(psi) / (1 - e2))
	lam := radian(cs.Lonf) + math.Asin(clamp(sinAlpha*math.Sin(j)/math.Cos(psi), -1, 1))

	return degree(lam), degree(phi), h
}

type GuamProjection struct {
	Lonf, Latf, Eastf, Northf float64
}

func (cs GuamProjection) String() string {
	return build("guam_projection").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs GuamProjection) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	phi := radian(lat)
	dLam := radian(lon) - radian(cs.Lonf)
	sinPhi, cosPhi := math.Sincos(phi)

	t := math.Sqrt(1 - s.E2()*sinPhi*sinPhi)
	x := s.A * dLam * cosPhi / t
	m0 := meridianDistance(s, radian(cs.Latf))
	m := meridianDistance(s, phi)
	east := cs.Eastf + x
	north := cs.Northf + m - m0 + 0.5*x*x*math.Tan(phi)*t/s.A

	return east, north, h
}

func (cs GuamProjection) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	m0 := meridianDistance(s, radian(cs.Latf))
	x := east - cs.Eastf
	x2 := 0.5 * x * x
	phi := radian(cs.Latf)

	for range 3 {
		sinPhi := math.Sin(phi)
		t := math.Sqrt(1 - s.E2()*sinPhi*sinPhi)
		phi = footpointLatitude(s, m0+(north-cs.Northf)-x2*math.Tan(phi)*t/s.A)
	}

	sinPhi, cosPhi := math.Sincos(phi)
	t := math.Sqrt(1 - s.E2()*sinPhi*sinPhi)
	lam := radian(cs.Lonf) + x*t/(s.A*cosPhi)

	return degree(lam), degree(phi), h
}
