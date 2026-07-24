package crs

import (
	"math"
)

// EqualEarth is EPSG method 1078 (PROJ +proj=eqearth).
// Formulas: IOGP Guidance Note 7-2 §3.4.4.
type EqualEarth struct {
	Lonf, Eastf, Northf float64
}

func (cs EqualEarth) String() string {
	return build("equal_earth").addAll(
		"lonf", cs.Lonf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

var equalEarthA = [...]float64{1.340264, -0.081106, 0.000893, 0.003796}

func (cs EqualEarth) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	e := s.E()
	e2 := s.E2()
	qp := authalicQ(1, e, e2)
	rqda := math.Sqrt(0.5 * qp)

	sbeta := math.Sin(radian(lat))
	if e >= 1e-12 {
		sbeta = authalicQ(sbeta, e, e2) / qp
		sbeta = clamp(sbeta, -1, 1)
	}

	m := math.Sqrt(3) / 2
	psi := math.Asin(m * sbeta)
	psi2 := psi * psi
	psi6 := psi2 * psi2 * psi2
	a1, a2, a3, a4 := equalEarthA[0], equalEarthA[1], equalEarthA[2], equalEarthA[3]
	lam := radian(lon - cs.Lonf)
	x := lam * math.Cos(psi) / (m * (a1 + 3*a2*psi2 + psi6*(7*a3+9*a4*psi2)))
	y := psi * (a1 + a2*psi2 + psi6*(a3+a4*psi2))

	return cs.Eastf + s.A*rqda*x, cs.Northf + s.A*rqda*y, h
}

func (cs EqualEarth) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	e := s.E()
	e2 := s.E2()
	qp := authalicQ(1, e, e2)
	rqda := math.Sqrt(0.5 * qp)
	a1, a2, a3, a4 := equalEarthA[0], equalEarthA[1], equalEarthA[2], equalEarthA[3]
	m := math.Sqrt(3) / 2

	x := (east - cs.Eastf) / (s.A * rqda)
	y := (north - cs.Northf) / (s.A * rqda)
	yc := clamp(y, -1.3173627591574, 1.3173627591574)

	for range 12 {
		y2 := yc * yc
		y6 := y2 * y2 * y2
		f := yc*(a1+a2*y2+y6*(a3+a4*y2)) - y
		fder := a1 + 3*a2*y2 + y6*(7*a3+9*a4*y2)
		dy := f / fder
		yc -= dy
		if math.Abs(dy) < 1e-11 {
			break
		}
	}

	y2 := yc * yc
	y6 := y2 * y2 * y2
	lam := m * x * (a1 + 3*a2*y2 + y6*(7*a3+9*a4*y2)) / math.Cos(yc)

	sbeta := math.Sin(yc) / m
	beta := math.Asin(clamp(sbeta, -1, 1))
	phi := beta

	if e >= 1e-12 {
		phi = authalicToGeodetic(beta, s)
	}

	return degree(lam + radian(cs.Lonf)), degree(phi), h
}
