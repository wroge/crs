package crs

import (
	"math"
	"testing"
)

func TestVincentyFlinders(t *testing.T) {
	wgs84 := Spheroid{A: 6378137, Fi: 298.257223563}
	lon1 := 144 + 25.0/60 + 29.52440/3600
	lat1 := -(37 + 57.0/60 + 3.72030/3600)
	lon2 := 143 + 55.0/60 + 35.38390/3600
	lat2 := -(37 + 39.0/60 + 10.15610/3600)
	wantAzi := radian(306 + 52.0/60 + 5.37/3600)

	dist, azi, ok := geodesicInverse(wgs84, lon1, lat1, lon2, lat2)
	if !ok {
		t.Fatal("inverse did not converge")
	}
	if math.Abs(dist-54972.271) > 1e-3 {
		t.Fatalf("distance %.6f", dist)
	}
	if math.Abs(wrapPi(azi-wantAzi)) > 1e-7 {
		t.Fatalf("azimuth %.10f rad", azi)
	}

	gotLon, gotLat, ok := geodesicDirect(wgs84, lon1, lat1, azi, dist)
	if !ok {
		t.Fatal("direct did not converge")
	}
	if math.Abs(gotLon-lon2) > 1e-8 || math.Abs(gotLat-lat2) > 1e-8 {
		t.Fatalf("direct (%.10f, %.10f)", gotLon, gotLat)
	}
}

func TestGeodesicSphereRoundTrip(t *testing.T) {
	sphere := Spheroid{A: 6378137}
	lon, lat, ok := geodesicDirect(sphere, 10, -20, radian(35), 250000)
	if !ok {
		t.Fatal("direct")
	}
	dist, azi, ok := geodesicInverse(sphere, 10, -20, lon, lat)
	if !ok {
		t.Fatal("inverse")
	}
	if math.Abs(dist-250000) > 1e-6 || math.Abs(wrapPi(azi-radian(35))) > 1e-10 {
		t.Fatalf("dist %.6f azi %.10f", dist, azi)
	}
}
