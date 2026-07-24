package crs

import (
	"math"
)

// ColombiaUrban is EPSG method 1052 (PROJ +proj=col_urban).
type ColombiaUrban struct {
	Lonf, Latf, Eastf, Northf, H0 float64
}

func (cs ColombiaUrban) String() string {
	return build("colombia_urban").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
		"h0", cs.H0,
	).String()
}

func (cs ColombiaUrban) consts(s Spheroid) (A, B, C, D, rho0, h0a float64) {
	h0a = cs.H0 / s.A
	phi0 := radian(cs.Latf)
	sin0 := math.Sin(phi0)
	nu0 := 1 / math.Sqrt(1-s.E2()*sin0*sin0)
	A = 1 + h0a/nu0
	rho0 = (1 - s.E2()) / math.Pow(1-s.E2()*sin0*sin0, 1.5)
	B = math.Tan(phi0) / (2 * rho0 * nu0)
	C = 1 + h0a
	D = rho0 * (1 + h0a/(1-s.E2()))

	return A, B, C, D, rho0, h0a
}

func (cs ColombiaUrban) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	A, B, _, _, rho0, h0a := cs.consts(s)
	phi := radian(lat)
	dLam := radian(lon - cs.Lonf)
	sinPhi := math.Sin(phi)
	nu := 1 / math.Sqrt(1-s.E2()*sinPhi*sinPhi)
	lamNuCos := dLam * nu * math.Cos(phi)
	east := cs.Eastf + s.A*A*lamNuCos

	sinPhiM := math.Sin(0.5 * (phi + radian(cs.Latf)))
	rhoM := (1 - s.E2()) / math.Pow(1-s.E2()*sinPhiM*sinPhiM, 1.5)
	G := 1 + h0a/rhoM
	north := cs.Northf + s.A*G*rho0*((phi-radian(cs.Latf))+B*lamNuCos*lamNuCos)

	return east, north, h
}

func (cs ColombiaUrban) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	_, B, C, D, _, _ := cs.consts(s)
	de := (east - cs.Eastf) / s.A
	dn := (north - cs.Northf) / s.A
	phi := radian(cs.Latf) + dn/D - B*(de/C)*(de/C)
	sinPhi := math.Sin(phi)
	nu := 1 / math.Sqrt(1-s.E2()*sinPhi*sinPhi)
	lon := radian(cs.Lonf) + de/(C*nu*math.Cos(phi))

	return degree(lon), degree(phi), h
}
