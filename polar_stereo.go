package crs

import (
	"math"
)

type PolarStereographicA struct {
	Lonf, Latf, Scale, Eastf, Northf float64
}

func (cs PolarStereographicA) String() string {
	return build("polar_stereographic_a").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"scale", cs.Scale,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs PolarStereographicA) south() bool {
	return cs.Latf < 0
}

func (cs PolarStereographicA) k0() float64 {
	if cs.Scale == 0 {
		return 1
	}
	return cs.Scale
}

func (cs PolarStereographicA) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	return polarStereoFromGeographic(s, lon, lat, h, cs.Lonf, cs.Eastf, cs.Northf, cs.k0(), cs.south())
}

func (cs PolarStereographicA) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	return polarStereoToGeographic(s, east, north, h, cs.Lonf, cs.Eastf, cs.Northf, cs.k0(), cs.south())
}

type PolarStereographicB struct {
	Lonf, Latf, Eastf, Northf float64
}

func (cs PolarStereographicB) String() string {
	return build("polar_stereographic_b").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs PolarStereographicB) south() bool {
	return cs.Latf < 0
}

func (cs PolarStereographicB) k0(s Spheroid) float64 {
	e := polarStereoE(s)
	mF, tF := polarStereoMF(s, radian(cs.Latf), cs.south())
	return mF * polarStereoC(e) / (2 * tF)
}

func (cs PolarStereographicB) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	return polarStereoFromGeographic(s, lon, lat, h, cs.Lonf, cs.Eastf, cs.Northf, cs.k0(s), cs.south())
}

func (cs PolarStereographicB) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	return polarStereoToGeographic(s, east, north, h, cs.Lonf, cs.Eastf, cs.Northf, cs.k0(s), cs.south())
}

// PolarStereographicC is EPSG method 9830.
// Latf is the standard parallel (latitude of false origin); FE/FN apply there, not at the pole.
type PolarStereographicC struct {
	Lonf, Latf, Eastf, Northf float64
}

func (cs PolarStereographicC) String() string {
	return build("polar_stereographic_c").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs PolarStereographicC) south() bool {
	return cs.Latf < 0
}

func (cs PolarStereographicC) k0(s Spheroid) float64 {
	e := polarStereoE(s)
	mF, tF := polarStereoMF(s, radian(cs.Latf), cs.south())
	return mF * polarStereoC(e) / (2 * tF)
}

func (cs PolarStereographicC) rhoF(s Spheroid) float64 {
	mF, _ := polarStereoMF(s, radian(cs.Latf), cs.south())
	return s.A * mF
}

// northfAtPole shifts NF so shared A/B helpers (origin at pole) match variant C.
func (cs PolarStereographicC) northfAtPole(s Spheroid) float64 {
	rhoF := cs.rhoF(s)
	if cs.south() {
		return cs.Northf - rhoF
	}
	return cs.Northf + rhoF
}

func (cs PolarStereographicC) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	return polarStereoFromGeographic(s, lon, lat, h, cs.Lonf, cs.Eastf, cs.northfAtPole(s), cs.k0(s), cs.south())
}

func (cs PolarStereographicC) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	return polarStereoToGeographic(s, east, north, h, cs.Lonf, cs.Eastf, cs.northfAtPole(s), cs.k0(s), cs.south())
}

func polarStereoE(s Spheroid) float64 {
	return math.Sqrt(s.E2())
}

func polarStereoC(e float64) float64 {
	return math.Exp(0.5 * ((1+e)*math.Log(1+e) + (1-e)*math.Log(1-e)))
}

func polarStereoT(phi, e float64, south bool) float64 {
	sinPhi := math.Sin(phi)
	es := e * sinPhi
	f := math.Pow((1+es)/(1-es), e/2)
	if south {
		return math.Tan(math.Pi/4+phi/2) / f
	}
	return math.Tan(math.Pi/4-phi/2) * f
}

func polarStereoMF(s Spheroid, phiF float64, south bool) (mF, tF float64) {
	e := polarStereoE(s)
	sinF := math.Sin(phiF)
	mF = math.Cos(phiF) / math.Sqrt(1-s.E2()*sinF*sinF)
	tF = polarStereoT(phiF, e, south)
	return mF, tF
}

func polarStereoFromGeographic(s Spheroid, lon, lat, h, lonf, eastf, northf, k0 float64, south bool) (float64, float64, float64) {
	e := polarStereoE(s)
	c := polarStereoC(e)
	phi := radian(lat)
	dLam := radian(lon - lonf)
	t := polarStereoT(phi, e, south)
	rho := 2 * s.A * k0 * t / c
	sinD := math.Sin(dLam)
	cosD := math.Cos(dLam)
	east := eastf + rho*sinD
	north := northf + rho*cosD
	if !south {
		north = northf - rho*cosD
	}
	return east, north, h
}

func polarStereoToGeographic(s Spheroid, east, north, h, lonf, eastf, northf, k0 float64, south bool) (float64, float64, float64) {
	e2 := s.E2()
	e4 := e2 * e2
	e6 := e4 * e2
	e8 := e4 * e4
	c := polarStereoC(polarStereoE(s))

	de := east - eastf
	dn := north - northf
	rho := math.Hypot(de, dn)
	t := rho * c / (2 * s.A * k0)

	var chi float64
	if south {
		chi = 2*math.Atan(t) - math.Pi/2
	} else {
		chi = math.Pi/2 - 2*math.Atan(t)
	}

	phi := chi +
		(e2/2+5*e4/24+e6/12+13*e8/360)*math.Sin(2*chi) +
		(7*e4/48+29*e6/240+811*e8/11520)*math.Sin(4*chi) +
		(7*e6/120+81*e8/1120)*math.Sin(6*chi) +
		(4279*e8/161280)*math.Sin(8*chi)

	var lon float64
	if de == 0 {
		lon = lonf
	} else if south {
		lon = lonf + degree(math.Atan2(de, dn))
	} else {
		lon = lonf + degree(math.Atan2(de, -dn))
	}
	return lon, degree(phi), h
}
