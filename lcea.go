package crs

import (
	"math"
)

type LambertCylindricalEqualArea struct {
	Lonf, Sp1, Eastf, Northf float64
}

func (cs LambertCylindricalEqualArea) String() string {
	return build("lambert_cylindrical_equal_area").addAll(
		"lonf", cs.Lonf,
		"sp1", cs.Sp1,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs LambertCylindricalEqualArea) k0(s Spheroid) float64 {
	phi1 := radian(cs.Sp1)
	sin1 := math.Sin(phi1)

	return math.Cos(phi1) / math.Sqrt(1-s.E2()*sin1*sin1)
}

func (cs LambertCylindricalEqualArea) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	k0 := cs.k0(s)
	east := cs.Eastf + s.A*k0*(radian(lon)-radian(cs.Lonf))
	north := cs.Northf + s.A*0.5*authalicQ(math.Sin(radian(lat)), s.E(), s.E2())/k0

	return east, north, h
}

func (cs LambertCylindricalEqualArea) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	k0 := cs.k0(s)
	qp := authalicQ(1, s.E(), s.E2())
	lam := radian(cs.Lonf) + (east-cs.Eastf)/(s.A*k0)
	sinBeta := clamp(2*(north-cs.Northf)*k0/(s.A*qp), -1, 1)
	beta := math.Asin(sinBeta)

	return degree(lam), degree(authalicToGeodetic(beta, s)), h
}

type LambertCylindricalEqualAreaSpherical struct {
	Lonf, Sp1, Eastf, Northf float64
}

func (cs LambertCylindricalEqualAreaSpherical) String() string {
	return build("lambert_cylindrical_equal_area_spherical").addAll(
		"lonf", cs.Lonf,
		"sp1", cs.Sp1,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs LambertCylindricalEqualAreaSpherical) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	r := s.A
	k0 := math.Cos(radian(cs.Sp1))
	east := cs.Eastf + r*k0*(radian(lon)-radian(cs.Lonf))
	north := cs.Northf + r*math.Sin(radian(lat))/k0

	return east, north, h
}

func (cs LambertCylindricalEqualAreaSpherical) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	r := s.A
	k0 := math.Cos(radian(cs.Sp1))
	lam := radian(cs.Lonf) + (east-cs.Eastf)/(r*k0)
	phi := math.Asin(clamp((north-cs.Northf)*k0/r, -1, 1))

	return degree(lam), degree(phi), h
}
