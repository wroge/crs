package crs

import (
	"math"
	"testing"
)

// TestGN72Examples checks forward results published in IOGP Guidance Note 7-2
// (September 2019, the edition in docs/IOGP_373-07-02_2019.pdf). The July 2026
// edition is not in the tree. Methods that only cite the GIGS test dataset
// are covered by TestProjSamples instead.
//
// Coordinates are the geographic CRS of the same datum, so the check is the
// projection and not a datum-grid hop. Longitude is relative to that CRS.
func TestGN72Examples(t *testing.T) {
	cases := []struct {
		name        string
		geo, proj   int
		lon, lat    float64
		east, north float64
		tol         float64
	}{
		// §3.1.1.1 Lambert Conic Conformal (2SP), NAD27 / Texas South Central. US survey feet.
		{"lcc_2sp", 4267, 32040, -96, 28.5, 2963503.91, 254759.80, 0.02},
		// §3.1.1.2 Lambert Conic Conformal (1SP), JAD69 / Jamaica National Grid.
		{"lcc_1sp", 4242, 24200, -76.94368333333334, 17.932166666666667, 255966.58, 142493.51, 0.02},
		// §3.1.1.4 Lambert Conic Conformal (2SP Belgium), Belge 1972 / Belge Lambert 72.
		{"lcc_belgium", 4313, 31300, 5.807370277777778, 50.6795725, 251763.20, 153034.13, 0.05},
		// §3.1.1.6 Lambert Conic Near-Conformal, Deir ez Zor / Levant Zone.
		{"lcca", 4227, 22700, 34.13646972222222, 37.5215625, 15707.96, 623165.96, 0.05},
		// §3.1.2.2 Krovak (North Orientated). The note's positive westing and southing
		// change sign for east-north axes (X = −westing, Y = −southing). Longitude is east of Ferro.
		{"krovak", 4818, 5221, 34.51643861111111, 50.20901166666667, -568991.00, -1050538.64, 0.05},
		// §3.1.2.4 Krovak Modified (North Orientated), same sign convention.
		{"krovak_modified", 5229, 5225, 34.51643861111111, 50.20901166666667, -5568990.91, -6050538.71, 0.05},
		// §3.2.1 Mercator (1SP), Makassar / NEIEZ.
		{"mercator_a", 4257, 3002, 120, -3, 5009726.58, 569150.82, 0.02},
		// §3.2.1 Mercator variant B, Pulkovo 1942 / Caspian Sea Mercator.
		{"mercator_b", 4284, 3388, 53, 53, 165704.29, 5171848.07, 0.05},
		// §3.2.1.2 Web Mercator, WGS 84 / Pseudo-Mercator.
		{"web_mercator", 4326, 3857, -100.33333333333333, 24.381786944444446, -11169055.58, 2800000.00, 0.02},
		// §3.2.2 Cassini-Soldner, Trinidad 1903 / Trinidad Grid. Clarke's links.
		{"cassini", 4302, 30200, -62, 10, 66644.94, 82536.22, 0.05},
		// §3.2.3.1 Transverse Mercator, OSGB 1936 / British National Grid.
		{"tmerc", 4277, 27700, 0.5, 50.5, 577274.99, 69740.50, 0.02},
		// §3.2.3.3 Transverse Mercator (South Orientated), Hartebeesthoek94 / Lo29.
		// The note quotes westing and southing. This CRS is easting/northing, so the
		// signs are opposite and the false origin is zero.
		{"tmerc_south", 4148, 2053, 28.282633055555556, -25.732028333333334, -71984.49, -2847342.74, 0.05},
		// §3.2.5 Equidistant Cylindrical, WGS 84 / World Equidistant Cylindrical.
		{"eqc", 4326, 4087, 10, 55, 1113194.91, 6097230.31, 0.02},
		// §3.3.1.3 Polar Stereographic (variant B), WGS 84 / Australian Antarctic.
		{"polar_b", 4326, 3032, 120, -75, 7255380.79, 7053389.56, 0.05},
		// §3.3.1.4 Polar Stereographic (variant C), Petrels 1972 / Terre Adelie.
		{"polar_c", 4636, 2985, 140.0714, -66.60522777777778, 303169.52, 244055.72, 0.05},
		// §3.3.2 Lambert Azimuthal Equal Area, ETRS89 / ETRS-LAEA.
		{"laea", 4258, 3035, 5, 50, 3962799.45, 2999718.85, 0.05},
		// §3.4.4 Equal Earth, WGS 84 / Equal Earth Americas.
		{"equal_earth", 4326, 8858, -117.19676361111111, 34.05754694444444, -2390749.042, 4242849.758, 0.05},
	}

	for _, c := range cases {
		t.Run(c.name, func(t *testing.T) {
			geo, err := Load(c.geo)
			if err != nil {
				t.Fatalf("Load(%d): %v", c.geo, err)
			}
			proj, err := Load(c.proj)
			if err != nil {
				t.Fatalf("Load(%d): %v", c.proj, err)
			}
			if geo.Datum.Name != proj.Datum.Name {
				t.Fatalf("datum %s and %s differ; this example is projection-only", geo.Datum.Name, proj.Datum.Name)
			}
			// FromGeographic skips the area-of-use gate. The Deir ez Zor example
			// sits just outside the published bbox and is still a valid formula check.
			east, north, h := proj.Conversion.FromGeographic(proj.Datum.Spheroid, c.lon, c.lat, 0)
			east, north, _ = proj.Unit.fromInternal(east, north, h, false)
			if math.Hypot(east-c.east, north-c.north) > c.tol {
				t.Errorf("got east %.4f north %.4f, want %.4f %.4f (tol %.3f)", east, north, c.east, c.north, c.tol)
			}
		})
	}
}
