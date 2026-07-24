package crs

import (
	"math"
)

type BonneSouthOrientated struct {
	Lonf, Latf, Eastf, Northf float64
}

func (cs BonneSouthOrientated) String() string {
	return build("bonne_south_orientated").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs BonneSouthOrientated) am1(s Spheroid) float64 {
	phi1 := radian(cs.Latf)
	sin1 := math.Sin(phi1)
	cos1 := math.Cos(phi1)
	return cos1 / (math.Sqrt(1-s.E2()*sin1*sin1) * sin1)
}

func (cs BonneSouthOrientated) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	phi1 := radian(cs.Latf)
	phi := radian(lat)
	lam := radian(lon) - radian(cs.Lonf)
	am1 := cs.am1(s)
	m1 := meridianDistance(s, phi1) / s.A
	m := meridianDistance(s, phi) / s.A
	sinPhi, cosPhi := math.Sincos(phi)
	rh := am1 + m1 - m

	var x, y float64

	if math.Abs(rh) > 1e-14 {
		e := cosPhi * lam / (rh * math.Sqrt(1-s.E2()*sinPhi*sinPhi))
		x = rh * math.Sin(e)
		y = am1 - rh*math.Cos(e)
	}

	return cs.Eastf - s.A*x, cs.Northf - s.A*y, h
}

func (cs BonneSouthOrientated) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	am1 := cs.am1(s)
	x := (cs.Eastf - east) / s.A
	y := am1 - (cs.Northf-north)/s.A
	phi1 := radian(cs.Latf)
	rh := math.Copysign(math.Hypot(x, y), phi1)
	m1 := meridianDistance(s, phi1) / s.A
	phi := footpointLatitude(s, s.A*(am1+m1-rh))
	sinPhi, cosPhi := math.Sincos(phi)

	var lam float64

	if math.Abs(math.Abs(phi)-math.Pi/2) <= 1e-14 {
		lam = 0
	} else {
		lm := rh * math.Sqrt(1-s.E2()*sinPhi*sinPhi) / cosPhi
		if phi1 > 0 {
			lam = lm * math.Atan2(x, y)
		} else {
			lam = lm * math.Atan2(-x, -y)
		}
	}

	return degree(lam + radian(cs.Lonf)), degree(phi), h
}
