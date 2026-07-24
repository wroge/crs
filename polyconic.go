package crs

import (
	"math"
)

type AmericanPolyconic struct {
	Lonf, Latf, Eastf, Northf float64
}

func (cs AmericanPolyconic) String() string {
	return build("american_polyconic").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs AmericanPolyconic) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	phi := radian(lat)
	phi0 := radian(cs.Latf)
	dLam := radian(lon - cs.Lonf)
	m0 := meridianDistance(s, phi0)

	if math.Abs(phi) <= 1e-10 {
		return cs.Eastf + s.A*dLam, cs.Northf - m0, h
	}

	sinPhi, cosPhi := math.Sincos(phi)
	ms := cosPhi / (math.Sqrt(1-s.E2()*sinPhi*sinPhi) * sinPhi)
	e := dLam * sinPhi
	east := cs.Eastf + s.A*ms*math.Sin(e)
	north := cs.Northf + (meridianDistance(s, phi) - m0) + s.A*ms*(1-math.Cos(e))

	return east, north, h
}

func (cs AmericanPolyconic) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	a := s.A
	e2 := s.E2()
	oneEs := 1 - e2
	m0 := meridianDistance(s, radian(cs.Latf))
	x := (east - cs.Eastf) / a
	y := (north - cs.Northf) / a
	y += m0 / a

	if math.Abs(y) <= 1e-10 {
		return degree(radian(cs.Lonf) + x), 0, h
	}

	r := y*y + x*x
	phi := y
	
	for range 20 {
		sinPhi, cosPhi := math.Sincos(phi)

		if math.Abs(cosPhi) < 1e-12 {
			return cs.Lonf, degree(math.Copysign(math.Pi/2, phi)), h
		}

		s2ph := sinPhi * cosPhi
		mlp := math.Sqrt(1 - e2*sinPhi*sinPhi)
		c := sinPhi * mlp / cosPhi
		ml := meridianDistance(s, phi) / a
		mlb := ml*ml + r
		mlp = oneEs / (mlp * mlp * mlp)
		dPhi := (ml + ml + c*mlb - 2*y*(c*ml+1)) /
			(e2*s2ph*(mlb-2*y*ml)/c + 2*(y-ml)*(c*mlp-1/s2ph) - mlp - mlp)
		phi += dPhi

		if math.Abs(dPhi) <= 1e-12 {
			break
		}
	}

	sinPhi := math.Sin(phi)
	if math.Abs(sinPhi) < 1e-12 {
		return degree(radian(cs.Lonf) + x), 0, h
	}

	lam := math.Asin(clamp(x*math.Tan(phi)*math.Sqrt(1-e2*sinPhi*sinPhi), -1, 1)) / sinPhi
	
	return degree(lam + radian(cs.Lonf)), degree(phi), h
}
