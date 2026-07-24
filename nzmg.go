package crs

import (
	"math"
)

type NewZealandMapGrid struct {
	Lonf, Latf, Eastf, Northf float64
}

func (cs NewZealandMapGrid) String() string {
	return build("new_zealand_map_grid").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

const (
	nzmgSec5ToRad = 0.4848136811095359935899141023
	nzmgRadToSec5 = 2.062648062470963551564733573
)

var nzmgBf = [...]complex128{
	complex(0.7557853228, 0),
	complex(0.249204646, 0.003371507),
	complex(-0.001541739, 0.041058560),
	complex(-0.10162907, 0.01727609),
	complex(-0.26623489, -0.36249218),
	complex(-0.6870983, -1.1651967),
}

var nzmgTpsi = [...]float64{
	0.6399175073, -0.1358797613, 0.063294409, -0.02526853, 0.0117879,
	-0.0055161, 0.0026906, -0.001333, 0.00067, -0.00034,
}

var nzmgTphi = [...]float64{
	1.5627014243, 0.5185406398, -0.03333098, -0.1052906, -0.0368594,
	0.007317, 0.01220, 0.00394, -0.0013,
}

func nzmgZpoly1(z complex128) complex128 {
	a := nzmgBf[5]

	for i := 4; i >= 0; i-- {
		a = nzmgBf[i] + z*a
	}

	return z * a
}

func nzmgZpolyd1(z complex128) (f, fp complex128) {
	a := nzmgBf[5]
	b := a
	for i := 4; i >= 0; i-- {
		b = a + z*b
		a = nzmgBf[i] + z*a
	}
	b = a + z*b

	return z * a, b
}

func (cs NewZealandMapGrid) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	dPhi := (radian(lat) - radian(cs.Latf)) * nzmgRadToSec5
	p := nzmgTpsi[9]

	for i := 8; i >= 0; i-- {
		p = nzmgTpsi[i] + dPhi*p
	}

	p *= dPhi
	z := nzmgZpoly1(complex(p, radian(lon)-radian(cs.Lonf)))

	return cs.Eastf + s.A*imag(z), cs.Northf + s.A*real(z), h
}

func (cs NewZealandMapGrid) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	z := complex((north-cs.Northf)/s.A, (east-cs.Eastf)/s.A)
	p := z

	for range 20 {
		f, fp := nzmgZpolyd1(p)
		f -= z
		den := real(fp)*real(fp) + imag(fp)*imag(fp)
		dp := complex(-(real(f)*real(fp)+imag(f)*imag(fp))/den, -(imag(f)*real(fp)-real(f)*imag(fp))/den)
		p += dp
		if math.Abs(real(dp))+math.Abs(imag(dp)) <= 1e-10 {
			break
		}
	}

	lam := imag(p)
	phiNom := nzmgTphi[8]

	for i := 7; i >= 0; i-- {
		phiNom = nzmgTphi[i] + real(p)*phiNom
	}

	phi := radian(cs.Latf) + real(p)*phiNom*nzmgSec5ToRad
	
	return degree(lam + radian(cs.Lonf)), degree(phi), h
}
