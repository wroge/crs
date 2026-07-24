package crs

import (
	"math"
)

type LocalOrthographic struct {
	Lonf, Latf, Azimuth, Scale, Eastf, Northf float64
}

func (cs LocalOrthographic) String() string {
	return build("local_orthographic").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"azimuth", cs.Azimuth,
		"scale", cs.Scale,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

func (cs LocalOrthographic) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	k0 := cs.Scale

	if k0 == 0 {
		k0 = 1
	}

	e2 := s.E2()
	phi0 := radian(cs.Latf)
	phi := radian(lat)
	lam := radian(lon) - radian(cs.Lonf)
	sinPhi0, cosPhi0 := math.Sincos(phi0)
	sinPhi, cosPhi := math.Sincos(phi)
	sinLam, cosLam := math.Sincos(lam)
	nu0 := 1 / math.Sqrt(1-e2*sinPhi0*sinPhi0)
	nu := 1 / math.Sqrt(1-e2*sinPhi*sinPhi)
	xp := nu * cosPhi * sinLam
	yp := nu*(sinPhi*cosPhi0-cosPhi*sinPhi0*cosLam) + e2*(nu0*sinPhi0-nu*sinPhi)*cosPhi0
	sinA, cosA := math.Sincos(radian(cs.Azimuth))
	east := cs.Eastf + s.A*k0*(cosA*xp-sinA*yp)
	north := cs.Northf + s.A*k0*(sinA*xp+cosA*yp)

	return east, north, h
}

func (cs LocalOrthographic) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	k0 := cs.Scale

	if k0 == 0 {
		k0 = 1
	}

	e2 := s.E2()
	sinA, cosA := math.Sincos(radian(cs.Azimuth))
	xf := (east - cs.Eastf) / (s.A * k0)
	yf := (north - cs.Northf) / (s.A * k0)
	x := cosA*xf + sinA*yf
	y := -sinA*xf + cosA*yf

	phi0 := radian(cs.Latf)
	sinPhi0, cosPhi0 := math.Sincos(phi0)
	nu0 := 1 / math.Sqrt(1-e2*sinPhi0*sinPhi0)
	yShift := e2 * nu0 * sinPhi0 * cosPhi0
	yScale := 1 / math.Sqrt(1-e2*cosPhi0*cosPhi0)
	yRec := (y - yShift) / yScale

	rh := math.Hypot(x, yRec)
	sinc := rh

	if sinc > 1 {
		sinc = 1
	}

	cosc := math.Sqrt(math.Max(0, 1-sinc*sinc))

	var phi, lam float64
	if rh < 1e-14 {
		phi, lam = phi0, 0
	} else {
		phi = math.Asin(clamp(cosc*sinPhi0+yRec*sinc*cosPhi0/rh, -1, 1))
		lam = math.Atan2(x*sinc, rh*cosPhi0*cosc-yRec*sinPhi0*sinc)
	}

	for range 20 {
		sinPhi, cosPhi := math.Sincos(phi)
		sinLam, cosLam := math.Sincos(lam)
		oneMinus := 1 - e2*sinPhi*sinPhi
		nu := 1 / math.Sqrt(oneMinus)
		xp := nu * cosPhi * sinLam
		yp := nu*(sinPhi*cosPhi0-cosPhi*sinPhi0*cosLam) + e2*(nu0*sinPhi0-nu*sinPhi)*cosPhi0
		rho := (1 - e2) * nu / oneMinus
		j11 := -rho * sinPhi * sinLam
		j12 := nu * cosPhi * cosLam
		j21 := rho * (cosPhi*cosPhi0 + sinPhi*sinPhi0*cosLam)
		j22 := nu * sinPhi0 * cosPhi * sinLam
		d := j11*j22 - j12*j21
		dx := x - xp
		dy := y - yp
		dPhi := (j22*dx - j12*dy) / d
		dLam := (-j21*dx + j11*dy) / d
		phi += dPhi
		lam += dLam
		if math.Abs(dPhi) < 1e-12 && math.Abs(dLam) < 1e-12 {
			break
		}
	}

	return degree(lam + radian(cs.Lonf)), degree(phi), h
}
