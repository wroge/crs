package crs

import (
	"fmt"
	"math"
)

type Spheroid struct {
	Name string
	A    float64
	Fi   float64
}

func (s Spheroid) String() string {
	return fmt.Sprintf("a=%s fi=%s", formatFloat(s.A), formatFloat(s.Fi))
}

func (s Spheroid) F() float64 {
	if s.Fi == 0 {
		return 0
	}

	return 1 / s.Fi
}

func (s Spheroid) A2() float64 {
	return s.A * s.A
}

func (s Spheroid) F2() float64 {
	f := s.F()

	return f * f
}

func (s Spheroid) B() float64 {
	return s.A * (1 - s.F())
}

func (s Spheroid) E2() float64 {
	if s.Fi == 0 {
		return 0
	}
	return 2/s.Fi - s.F2()
}

func (s Spheroid) E() float64 {

	return math.Sqrt(s.E2())

}

func (s Spheroid) E4() float64 {
	e2 := s.E2()

	return e2 * e2
}

func (s Spheroid) E6() float64 {
	e2 := s.E2()

	return e2 * e2 * e2
}

func (s Spheroid) Ei() float64 {
	e2 := s.E2()
	t := math.Sqrt(1 - e2)

	return (1 - t) / (1 + t)
}

func (s Spheroid) Ei2() float64 {
	ei := s.Ei()

	return ei * ei
}

func (s Spheroid) Ei3() float64 {
	ei := s.Ei()

	return ei * ei * ei
}

func (s Spheroid) Ei4() float64 {
	ei := s.Ei()

	return ei * ei * ei * ei
}

func (s Spheroid) GeographicToGeocentric(lon, lat, h float64) (x, y, z float64) {
	n := s.A / math.Sqrt(1-s.E2()*intPow(math.Sin(radian(lat)), 2))

	x = (n + h) * math.Cos(radian(lon)) * math.Cos(radian(lat))
	y = (n + h) * math.Cos(radian(lat)) * math.Sin(radian(lon))
	z = (n*intPow(s.A*(1-s.F()), 2)/(s.A2()) + h) * math.Sin(radian(lat))

	return x, y, z
}

func (s Spheroid) GeocentricToGeographic(x, y, z float64) (lon, lat, h float64) {
	sd := math.Sqrt(x*x + y*y)
	T := math.Atan(z * s.A / (sd * s.B()))
	B := math.Atan((z + s.E2()*(s.A2())/s.B()*
		intPow(math.Sin(T), 3)) / (sd - s.E2()*s.A*intPow(math.Cos(T), 3)))
	n := s.A / math.Sqrt(1-s.E2()*intPow(math.Sin(B), 2))
	h = sd/math.Cos(B) - n
	lon = degree(math.Atan2(y, x))
	lat = degree(B)

	return lon, lat, h
}
