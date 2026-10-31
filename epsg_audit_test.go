package crs

import (
	"fmt"
	"math"
	"path"
	"sort"
	"strconv"
	"strings"
	"testing"
)

// TestProjectionRoundTrip projects several in-bounds coordinates through every
// EPSG code of every conversion and checks that the inverse comes back.
// One failing code fails the projection's subtest.
func TestProjectionRoundTrip(t *testing.T) {
	byConv := map[string][]int{}

	entries, err := epsgDir.ReadDir("epsg")
	if err != nil {
		t.Fatal(err)
	}
	for _, e := range entries {
		name := e.Name()
		if !strings.HasSuffix(name, ".txt") {
			continue
		}
		code, err := strconv.Atoi(strings.TrimSuffix(name, ".txt"))
		if err != nil {
			t.Errorf("bad filename %s", name)
			continue
		}
		raw, err := epsgDir.ReadFile(path.Join("epsg", name))
		if err != nil {
			t.Errorf("EPSG:%d: %v", code, err)
			continue
		}
		conv := fieldOf(string(raw), "conversion")
		if conv == "" {
			t.Errorf("EPSG:%d: missing conversion", code)
			continue
		}
		byConv[conv] = append(byConv[conv], code)
	}

	for _, name := range requiredCoordinateSystems {
		if len(byConv[name]) == 0 {
			t.Errorf("no EPSG codes exercise %s", name)
		}
	}

	names := make([]string, 0, len(byConv))
	for name := range byConv {
		names = append(names, name)
	}
	sort.Strings(names)

	for _, name := range names {
		codes := append([]int(nil), byConv[name]...)
		sort.Ints(codes)
		t.Run(name, func(t *testing.T) {
			var fails []string
			for _, code := range codes {
				msg := checkEPSG(code)
				if msg != "" {
					fails = append(fails, fmt.Sprintf("EPSG:%d: %s", code, msg))
				}
			}
			t.Logf("%d codes", len(codes))
			if len(fails) == 0 {
				return
			}
			shown := fails
			if len(shown) > 8 {
				shown = shown[:8]
			}
			t.Errorf("%d/%d codes failed\n%s", len(fails), len(codes), strings.Join(shown, "\n"))
		})
	}
}

// conversion=utm is a loader alias of transverse Mercator and has no EPSG file.
func TestUTMAliasRoundTrip(t *testing.T) {
	defs := []string{
		"conversion=utm zone=32 datum=wgs84 bbox=6,0,12,84",
		"conversion=utm zone=1 southern datum=wgs84 bbox=-180,-80,-174,0",
		"conversion=utm zone=60 datum=wgs84 bbox=174,0,180,84",
	}
	for _, def := range defs {
		c, err := Load(def)
		if err != nil {
			t.Errorf("Load %q: %v", def, err)
			continue
		}
		if msg := roundTripCRS(c); msg != "" {
			t.Errorf("%s: %s", def, msg)
		}
	}
}

func checkEPSG(code int) (msg string) {
	defer func() {
		if r := recover(); r != nil {
			msg = fmt.Sprintf("panic: %v", r)
		}
	}()

	c, err := Load(code)
	if err != nil {
		return "load: " + err.Error()
	}
	return roundTripCRS(c)
}

func fieldOf(txt, key string) string {
	for _, tok := range strings.Fields(txt) {
		k, v, ok := strings.Cut(tok, "=")
		if ok && k == key {
			return v
		}
	}
	return ""
}

func roundTripCRS(c CoordinateReferenceSystem) string {
	pts := samplePoints(c.BoundingBox)
	if len(pts) == 0 {
		return "no sample point"
	}

	// Area-of-use boxes are Greenwich longitudes. Projection parameters and
	// geographic coordinates are relative to the datum prime meridian.
	pm := c.Datum.primeMeridianLongitude()
	for i := range pts {
		pts[i][0] -= pm
	}

	s := c.Datum.Spheroid
	if s.A == 0 || math.IsNaN(s.A) || math.IsInf(s.A, 0) {
		return fmt.Sprintf("bad spheroid a=%v fi=%v", s.A, s.Fi)
	}

	return roundTripPoints(c, s, pts)
}

func roundTripPoints(c CoordinateReferenceSystem, s Spheroid, pts [][2]float64) string {
	const (
		tolDeg = 2e-5 // about 2 m
		height = 42.0
	)

	var worst string
	worstErr := 0.0
	for _, p := range pts {
		for _, hIn := range []float64{0, height} {
			east, north, h := c.Conversion.FromGeographic(s, p[0], p[1], hIn)
			if !finite3(east, north, h) {
				return fmt.Sprintf("forward non-finite at lon=%.6f lat=%.6f h=%.1f -> %v %v %v", p[0], p[1], hIn, east, north, h)
			}
			lon, lat, h2 := c.Conversion.ToGeographic(s, east, north, h)
			if !finite3(lon, lat, h2) {
				return fmt.Sprintf("inverse non-finite at lon=%.6f lat=%.6f h=%.1f", p[0], p[1], hIn)
			}
			dLon := math.Abs(lon - p[0])
			if dLon > 180 {
				dLon = 360 - dLon
			}
			dLat := math.Abs(lat - p[1])
			errDeg := math.Hypot(dLon, dLat)
			if errDeg > tolDeg && errDeg > worstErr {
				worstErr = errDeg
				worst = fmt.Sprintf("round-trip %.6f deg (%.1f m) at lon=%.6f lat=%.6f h=%.1f got lon=%.6f lat=%.6f",
					errDeg, errDeg*111320, p[0], p[1], hIn, lon, lat)
			}
			if math.Abs(h2-hIn) > 1e-4 {
				return fmt.Sprintf("height drift %g at lon=%.6f lat=%.6f", h2-hIn, p[0], p[1])
			}
		}
	}
	return worst
}

func finite3(a, b, c float64) bool {
	return !math.IsNaN(a) && !math.IsNaN(b) && !math.IsNaN(c) &&
		!math.IsInf(a, 0) && !math.IsInf(b, 0) && !math.IsInf(c, 0)
}

// samplePoints returns coordinates strictly inside b, including boxes that
// cross the antimeridian. Longitude is Greenwich, matching EPSG areas of use.
func samplePoints(b BoundingBox) [][2]float64 {
	if b.MaxLat < b.MinLat {
		return nil
	}

	fractions := []float64{0.2, 0.5, 0.8}
	lats := mix(b.MinLat, b.MaxLat, fractions)

	var lons []float64
	if b.MinLon <= b.MaxLon {
		lons = mix(b.MinLon, b.MaxLon, fractions)
	} else {
		span := (180 - b.MinLon) + (b.MaxLon + 180)
		for _, f := range fractions {
			lon := b.MinLon + f*span
			if lon > 180 {
				lon -= 360
			}
			lons = append(lons, lon)
		}
	}

	var out [][2]float64
	seen := map[[2]float64]bool{}
	for _, lon := range lons {
		for _, lat := range lats {
			if lat > 89.5 {
				lat = 89.5
			}
			if lat < -89.5 {
				lat = -89.5
			}
			p := [2]float64{lon, lat}
			if seen[p] {
				continue
			}
			seen[p] = true
			out = append(out, p)
		}
	}
	return out
}

func mix(min, max float64, fractions []float64) []float64 {
	span := max - min
	out := make([]float64, len(fractions))
	for i, f := range fractions {
		out[i] = min + f*span
	}
	return out
}

func TestUTMZonedGridPrefix(t *testing.T) {
	cases := []struct {
		grid, zone int
		lon, lat   float64
		zoneNum    float64
	}{
		{32600, 32632, 9, 50, 32},
		{32600, 32601, -177, 10, 1},
		{32600, 32660, 177, 10, 60},
		{32700, 32732, 9, -40, 32},
		{32700, 32701, -179, -20, 1},
	}

	for _, tc := range cases {
		grid, err := Load(tc.grid)
		if err != nil {
			t.Fatal(err)
		}
		zone, err := Load(tc.zone)
		if err != nil {
			t.Fatal(err)
		}

		ge, gn, _ := grid.Conversion.FromGeographic(grid.Datum.Spheroid, tc.lon, tc.lat, 0)
		ze, zn, _ := zone.Conversion.FromGeographic(zone.Datum.Spheroid, tc.lon, tc.lat, 0)
		if math.Abs(ge-(tc.zoneNum*1e6+ze)) > 0.01 || math.Abs(gn-zn) > 0.01 {
			t.Errorf("EPSG:%d lon=%.1f lat=%.1f: grid (%.3f, %.3f) zone (%.3f, %.3f)",
				tc.grid, tc.lon, tc.lat, ge, gn, ze, zn)
		}

		lon, lat, _ := grid.Conversion.ToGeographic(grid.Datum.Spheroid, ge, gn, 0)
		if math.Abs(lon-tc.lon) > 1e-6 || math.Abs(lat-tc.lat) > 1e-6 {
			t.Errorf("EPSG:%d inverse: got (%.8f, %.8f) want (%.1f, %.1f)", tc.grid, lon, lat, tc.lon, tc.lat)
		}
	}
}
