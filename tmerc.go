package crs

import (
	"math"
)

type TransverseMercator struct {
	Lonf      float64
	Latf      float64
	Scale     float64
	Eastf     float64
	Northf    float64
	ZoneWidth float64
}

func (cs TransverseMercator) String() string {
	return build("transverse_mercator").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"scale", cs.Scale,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
		"zone_width", cs.ZoneWidth,
	).String()
}

func (cs TransverseMercator) zone(lon float64) int {
	w := cs.ZoneWidth
	d := lon - cs.Lonf

	for d < 0 {
		d += 360
	}
	for d >= 360 {
		d -= 360
	}

	z := 1 + int(math.Floor(d/w))
	nZones := int(math.Round(360 / w))

	if z > nZones {
		z = nZones
	}

	if z < 1 {
		z = 1
	}

	return z
}

func (cs TransverseMercator) centralMeridian(zone int) float64 {
	return cs.Lonf + (float64(zone)-0.5)*cs.ZoneWidth
}

// fixedLonf returns a single-zone TM with the given central meridian.
func (cs TransverseMercator) fixedLonf(lonf float64) TransverseMercator {
	cs.Lonf = lonf
	cs.ZoneWidth = 0

	return cs
}

func (cs TransverseMercator) FromGeographic(s Spheroid, lon, lat, h float64) (float64, float64, float64) {
	if cs.ZoneWidth > 0 {
		return cs.fixedLonf(cs.centralMeridian(cs.zone(lon))).FromGeographic(s, lon, lat, h)
	}

	phiO := radian(cs.Latf)
	lambdaO := radian(cs.Lonf)
	n := s.F() / (2 - s.F())
	n2 := n * n
	n3 := n * n * n
	n4 := n * n * n * n
	B := (s.A / (1 + n)) * (1 + n2/4 + n4/64)
	h1 := n/2.0 - (2/3.0)*n2 + (5/16.0)*n3 + (41/180.0)*n4
	h2 := (13/48.0)*n2 - (3/5.0)*n3 + (557/1440.0)*n4
	h3 := (61/240.0)*n3 - (103/140.0)*n4
	h4 := (49561 / 161280.0) * n4

	var MO float64

	switch phiO {
	case 0:
		MO = 0
	case math.Pi / 2:
		MO = B * (math.Pi / 2)
	case -math.Pi / 2:
		MO = B * (-math.Pi / 2)
	default:
		Q0 := math.Asinh(math.Tan(phiO)) - (s.E() * math.Atanh(s.E()*math.Sin(phiO)))
		xi00 := math.Atan(math.Sinh(Q0))
		xi01 := h1 * math.Sin(2*xi00)
		xi02 := h2 * math.Sin(4*xi00)
		xi03 := h3 * math.Sin(6*xi00)
		xi04 := h4 * math.Sin(8*xi00)
		xi0 := xi00 + xi01 + xi02 + xi03 + xi04
		MO = B * xi0
	}

	phi := radian(lat)
	lambda := radian(lon)

	Q := math.Asinh(math.Tan(phi)) - s.E()*math.Atanh(s.E()*math.Sin(phi))
	beta := math.Atan(math.Sinh(Q))
	eta0 := math.Atanh(math.Cos(beta) * math.Sin(lambda-lambdaO))
	xi0 := math.Asin(math.Sin(beta) * math.Cosh(eta0))

	xi1 := h1 * math.Sin(2*xi0) * math.Cosh(2*eta0)
	xi2 := h2 * math.Sin(4*xi0) * math.Cosh(4*eta0)
	xi3 := h3 * math.Sin(6*xi0) * math.Cosh(6*eta0)
	xi4 := h4 * math.Sin(8*xi0) * math.Cosh(8*eta0)

	xi := xi0 + xi1 + xi2 + xi3 + xi4

	eta1 := h1 * math.Cos(2*xi0) * math.Sinh(2*eta0)
	eta2 := h2 * math.Cos(4*xi0) * math.Sinh(4*eta0)
	eta3 := h3 * math.Cos(6*xi0) * math.Sinh(6*eta0)
	eta4 := h4 * math.Cos(8*xi0) * math.Sinh(8*eta0)

	eta := eta0 + eta1 + eta2 + eta3 + eta4

	east := cs.Eastf + cs.Scale*B*eta
	north := cs.Northf + cs.Scale*(B*xi-MO)

	return east, north, h
}

func (cs TransverseMercator) ToGeographic(s Spheroid, east, north, h float64) (float64, float64, float64) {
	if cs.ZoneWidth > 0 {
		nZones := int(math.Round(360 / cs.ZoneWidth))
		bestDist := math.Inf(1)
		bestLon, bestLat, bestH := 0.0, 0.0, h
		found := false

		for z := 1; z <= nZones; z++ {
			lon, lat, hOut := cs.fixedLonf(cs.centralMeridian(z)).ToGeographic(s, east, north, h)
			minLon := cs.Lonf + float64(z-1)*cs.ZoneWidth
			maxLon := minLon + cs.ZoneWidth
			lonN := lon
			for lonN < cs.Lonf {
				lonN += 360
			}
			for lonN >= cs.Lonf+360 {
				lonN -= 360
			}
			if lonN < minLon || lonN >= maxLon {
				continue
			}
			dist := math.Abs(lon - cs.centralMeridian(z))
			if dist < bestDist {
				bestDist = dist
				bestLon, bestLat, bestH = lon, lat, hOut
				found = true
			}
		}

		if found {
			return bestLon, bestLat, bestH
		}
		
		return cs.fixedLonf(cs.centralMeridian(1)).ToGeographic(s, east, north, h)
	}

	phiO := radian(cs.Latf)
	lambdaO := radian(cs.Lonf)
	n := s.F() / (2 - s.F())
	n2 := n * n
	n3 := n * n * n
	n4 := n * n * n * n
	B := (s.A / (1 + n)) * (1 + n2/4 + n4/64)
	h1 := n/2.0 - (2/3.0)*n2 + (5/16.0)*n3 + (41/180.0)*n4
	h2 := (13/48.0)*n2 - (3/5.0)*n3 + (557/1440.0)*n4
	h3 := (61/240.0)*n3 - (103/140.0)*n4
	h4 := (49561 / 161280.0) * n4
	e := math.Sqrt(s.E2())

	var MO float64

	switch phiO {
	case 0:
		MO = 0
	case math.Pi / 2:
		MO = B * (math.Pi / 2)
	case -math.Pi / 2:
		MO = B * (-math.Pi / 2)
	default:
		Q0 := math.Asinh(math.Tan(phiO)) - (e * math.Atanh(e*math.Sin(phiO)))
		xi00 := math.Atan(math.Sinh(Q0))
		xi01 := h1 * math.Sin(2*xi00)
		xi02 := h2 * math.Sin(4*xi00)
		xi03 := h3 * math.Sin(6*xi00)
		xi04 := h4 * math.Sin(8*xi00)
		xi0 := xi00 + xi01 + xi02 + xi03 + xi04
		MO = B * xi0
	}

	h1i := n/2.0 - (2/3.0)*n2 + (37/96.0)*n3 + (1/360.0)*n4
	h2i := (1/48.0)*n2 - (1/15.0)*n3 + (437/1440.0)*n4
	h3i := (17/480.0)*n3 - (37/840.0)*n4
	h4i := (4397 / 161280.0) * n4

	etai := (east - cs.Eastf) / (B * cs.Scale)
	xii := ((north - cs.Northf) + cs.Scale*MO) / (B * cs.Scale)

	xi1i := h1i * math.Sin(2*xii) * math.Cosh(2*etai)
	xi2i := h2i * math.Sin(4*xii) * math.Cosh(4*etai)
	xi3i := h3i * math.Sin(6*xii) * math.Cosh(6*etai)
	xi4i := h4i * math.Sin(8*xii) * math.Cosh(8*etai)

	xi0i := xii - (xi1i + xi2i + xi3i + xi4i)

	eta1i := h1i * math.Cos(2*xii) * math.Sinh(2*etai)
	eta2i := h2i * math.Cos(4*xii) * math.Sinh(4*etai)
	eta3i := h3i * math.Cos(6*xii) * math.Sinh(6*etai)
	eta4i := h4i * math.Cos(8*xii) * math.Sinh(8*etai)

	eta0i := etai - (eta1i + eta2i + eta3i + eta4i)

	betai := math.Asin(math.Sin(xi0i) / math.Cosh(eta0i))
	Qi := math.Asinh(math.Tan(betai))
	Qii := Qi + (s.E() * math.Atanh(s.E()*math.Tanh(Qi)))

	for range 15 {
		newQ := Qi + (s.E() * math.Atanh(s.E()*math.Tanh(Qii)))
		if math.Abs(newQ-Qii) < 1e-14 {
			break
		}

		Qii = newQ
	}

	phi := math.Atan(math.Sinh(Qii))
	lambda := lambdaO + math.Asin(math.Tanh(eta0i)/math.Cos(betai))

	return degree(lambda), degree(phi), h
}
