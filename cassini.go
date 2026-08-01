package crs

import (
	"math"
)

type CassiniSoldner struct {
	Lonf, Latf, Eastf, Northf float64
}

func (cs CassiniSoldner) String() string {
	return build("cassini_soldner").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func meridianDistance(s Spheroid, phi float64) float64 {
	e2 := s.E2()
	e4 := e2 * e2
	e6 := e4 * e2

	return s.A * ((1-e2/4-3*e4/64-5*e6/256)*phi -
		(3*e2/8+3*e4/32+45*e6/1024)*math.Sin(2*phi) +
		(15*e4/256+45*e6/1024)*math.Sin(4*phi) -
		(35*e6/3072)*math.Sin(6*phi))
}

func footpointLatitude(s Spheroid, m float64) float64 {
	e2 := s.E2()
	e4 := e2 * e2
	e6 := e4 * e2
	e1 := (1 - math.Sqrt(1-e2)) / (1 + math.Sqrt(1-e2))
	e12 := e1 * e1
	e13 := e12 * e1
	e14 := e12 * e12

	mu := m / (s.A * (1 - e2/4 - 3*e4/64 - 5*e6/256))

	return mu +
		(3*e1/2-27*e13/32)*math.Sin(2*mu) +
		(21*e12/16-55*e14/32)*math.Sin(4*mu) +
		(151*e13/96)*math.Sin(6*mu) +
		(1097*e14/512)*math.Sin(8*mu)
}

func (cs CassiniSoldner) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	phi := radian(lat)
	lam := radian(lon)
	phi0 := radian(cs.Latf)
	lam0 := radian(cs.Lonf)

	sinPhi, cosPhi := math.Sincos(phi)
	e2 := s.E2()
	nu := s.A / math.Sqrt(1-e2*sinPhi*sinPhi)
	t := math.Tan(phi)
	t2 := t * t
	c := e2 * cosPhi * cosPhi / (1 - e2)
	a := (lam - lam0) * cosPhi
	a2 := a * a
	a3 := a2 * a
	a4 := a2 * a2
	a5 := a4 * a

	east := cs.Eastf + nu*(a-t2*a3/6-(8-t2+8*c)*t2*a5/120)
	north := cs.Northf + (meridianDistance(s, phi) - meridianDistance(s, phi0)) +
		nu*t*(a2/2+(5-t2+6*c)*a4/24)

	return east, north, h
}

func (cs CassiniSoldner) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	phi0 := radian(cs.Latf)
	lam0 := radian(cs.Lonf)
	e2 := s.E2()

	m1 := meridianDistance(s, phi0) + (north - cs.Northf)
	phi1 := footpointLatitude(s, m1)

	if math.Abs(math.Abs(phi1)-math.Pi/2) < 1e-12 {
		return degree(lam0), degree(phi1), h
	}

	sinPhi1, cosPhi1 := math.Sincos(phi1)
	nu1 := s.A / math.Sqrt(1-e2*sinPhi1*sinPhi1)
	rho1 := s.A * (1 - e2) / math.Pow(1-e2*sinPhi1*sinPhi1, 1.5)
	t1 := math.Tan(phi1)
	t12 := t1 * t1
	d := (east - cs.Eastf) / nu1
	d2 := d * d
	d3 := d2 * d
	d4 := d2 * d2
	d5 := d4 * d

	phi := phi1 - (nu1*t1/rho1)*(d2/2-(1+3*t12)*d4/24)
	lam := lam0 + (d-t12*d3/3+(1+3*t12)*t12*d5/15)/cosPhi1

	return degree(lam), degree(phi), h
}
