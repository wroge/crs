package crs

import (
	"math"
	"testing"
)

func TestAzimuthalEquidistantEllipsoid(t *testing.T) {
	// cs2cs EPSG:4326 → EPSG:27701 / EPSG:27705, WGS84, metres.
	cases := []struct {
		code           int
		lon, lat, e, n float64
	}{
		{27701, -5.67747881, 1.63138717, 2591027.07343396, 5310164.44164079},
		{27701, 23.9998013, 14.95484472, 5890925.66637192, 6705768.98464155},
		{27701, 0.99913049, 20.31608451, 3468215.72952056, 7382267.91512091},
		{27701, -4.14953775, -18.13248039, 2808531.26850418, 3039436.32614133},
		{27705, -162, 12.081, 1225374.59214, 3265885.58203},
		{27705, 162, 85.899, 7778192.97096, 9182141.03158},
	}

	for _, c := range cases {
		cs, err := Load(c.code)
		if err != nil {
			t.Fatal(err)
		}
		east, north, _ := cs.Conversion.FromGeographic(cs.Datum.Spheroid, c.lon, c.lat, 0)
		if d := math.Hypot(east-c.e, north-c.n); d > 0.01 {
			t.Errorf("EPSG:%d (%.4f, %.4f): %.4f m off, got %.4f %.4f", c.code, c.lon, c.lat, d, east, north)
		}
		lon, lat, _ := cs.Conversion.ToGeographic(cs.Datum.Spheroid, east, north, 0)
		if math.Abs(math.Remainder(lon-c.lon, 360)) > 2e-5 || math.Abs(lat-c.lat) > 2e-5 {
			t.Errorf("EPSG:%d inverse (%.8f, %.8f)", c.code, lon, lat)
		}
	}
}

func TestAzimuthalEquidistantSouthPole(t *testing.T) {
	cs, err := Load(27702)
	if err != nil {
		t.Fatal(err)
	}
	s := cs.Datum.Spheroid
	east, north, _ := cs.Conversion.FromGeographic(s, -40, -70, 0)
	lon, lat, _ := cs.Conversion.ToGeographic(s, east, north, 0)
	if math.Abs(math.Remainder(lon-(-40), 360)) > 1e-8 || math.Abs(lat-(-70)) > 1e-8 {
		t.Fatalf("pole round-trip (%.10f, %.10f)", lon, lat)
	}
	if math.Hypot(east-3714266.977, north-3402016.506) < 1 {
		t.Fatal("south-pole origin should move a point at -70")
	}
}
