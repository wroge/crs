package crs

import (
	"math"
)

// LambertConformalConic1SPVariantB is EPSG method 1102.
// Same cone as LCC 1SP at Lat0 (natural origin); FE/FN apply at Latf (false origin).
type LambertConformalConic1SPVariantB struct {
	Lonf, Lat0, Latf, Scale, Eastf, Northf float64
}

func (cs LambertConformalConic1SPVariantB) String() string {
	return build("lambert_conformal_conic_1sp_variant_b").addAll(
		"lonf", cs.Lonf,
		"lat0", cs.Lat0,
		"latf", cs.Latf,
		"scale", cs.Scale,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs LambertConformalConic1SPVariantB) consts(s Spheroid) (n, f, rf float64) {
	k0 := cs.Scale

	if k0 == 0 {
		k0 = 1
	}

	phi0 := radian(cs.Lat0)
	phif := radian(cs.Latf)
	t0 := math.Tan(math.Pi/4-phi0/2) / math.Pow((1-s.E()*math.Sin(phi0))/(1+s.E()*math.Sin(phi0)), s.E()/2)
	tf := math.Tan(math.Pi/4-phif/2) / math.Pow((1-s.E()*math.Sin(phif))/(1+s.E()*math.Sin(phif)), s.E()/2)
	m0 := math.Cos(phi0) / math.Sqrt(1-s.E2()*sin2(phi0))
	n = math.Sin(phi0)
	f = k0 * m0 / (n * math.Pow(t0, n))
	rf = s.A * f * math.Pow(tf, n)

	return n, f, rf
}

func (cs LambertConformalConic1SPVariantB) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	n, f, rf := cs.consts(s)
	phi := radian(lat)
	t := math.Tan(math.Pi/4-phi/2) / math.Pow((1-s.E()*math.Sin(phi))/(1+s.E()*math.Sin(phi)), s.E()/2)
	r := s.A * f * math.Pow(t, n)
	theta := n * (radian(lon) - radian(cs.Lonf))
	east := cs.Eastf + r*math.Sin(theta)
	north := cs.Northf + rf - r*math.Cos(theta)

	return east, north, h
}

func (cs LambertConformalConic1SPVariantB) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	n, f, rf := cs.consts(s)
	ri := math.Hypot(east-cs.Eastf, rf-(north-cs.Northf))
	ti := math.Pow(ri/(s.A*f), 1/n)
	theta := math.Atan2(east-cs.Eastf, rf-(north-cs.Northf))
	phi := math.Pi/2 - 2*math.Atan(ti)

	for range 6 {
		next := math.Pi/2 - 2*math.Atan(ti*math.Pow((1-s.E()*math.Sin(phi))/(1+s.E()*math.Sin(phi)), s.E()/2))
		if math.Abs(next-phi) < 1e-12 {
			phi = next
			break
		}
		phi = next
	}
	
	return degree(theta/n + radian(cs.Lonf)), degree(phi), h
}
