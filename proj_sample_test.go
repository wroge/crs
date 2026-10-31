package crs

import (
	"math"
	"strings"
	"testing"
)

// projSamples and projUncompared are filled by proj_sample_data_test.go,
// which tools/projcheck regenerates. go test does not call PROJ.
var (
	projSamples    []projSample
	projUncompared map[string]string
)

type projSample struct {
	name         string
	conv         string
	from, to     int
	inA, inB     float64
	wantA, wantB float64
	wantC        float64
	tol          float64
	deg          bool
	edge         bool
}

// TestProjSamples checks several EPSG codes of each conversion against cs2cs
// values taken inside the area of use. Edge samples belong to series methods
// and carry the measured accuracy budget instead of the 1 cm interior limit.
func TestProjSamples(t *testing.T) {
	if len(projSamples) == 0 {
		t.Fatal("no PROJ samples; regenerate with go run ./tools/projcheck")
	}

	seen := map[string]int{}
	for _, s := range projSamples {
		seen[s.conv]++
		t.Run(s.name, func(t *testing.T) {
			from, err := Load(s.from)
			if err != nil {
				t.Fatalf("Load(%d): %v", s.from, err)
			}
			to, err := Load(s.to)
			if err != nil {
				t.Fatalf("Load(%d): %v", s.to, err)
			}
			from = from.WithUnit(Degree)
			if s.deg {
				to = to.WithUnit(Degree)
			} else {
				to = to.WithUnit(Metre)
			}
			fn, err := from.TransformTo(to)
			if err != nil {
				t.Fatalf("TransformTo: %v", err)
			}
			a, b, c, err := fn(s.inA, s.inB, 0)
			if err != nil {
				t.Fatalf("transform: %v", err)
			}
			da, db := a-s.wantA, b-s.wantB
			if s.deg {
				da = math.Remainder(da, 360)
				db = math.Remainder(db, 360)
			}
			if math.Hypot(da, db) > s.tol {
				t.Errorf("got (%.6f, %.6f), want (%.6f, %.6f), tol %g", a, b, s.wantA, s.wantB, s.tol)
			}
			if math.Abs(c-s.wantC) > 0.01 {
				t.Errorf("height got %.6f, want %.6f", c, s.wantC)
			}
		})
	}

	for _, conv := range requiredCoordinateSystems {
		if seen[conv] > 0 {
			continue
		}
		reason := projUncompared[conv]
		if reason == "" {
			t.Errorf("no PROJ sample for %s", conv)
			continue
		}
		if strings.HasPrefix(reason, "PROJ") {
			t.Logf("no PROJ sample for %s: %s", conv, reason)
			continue
		}
		t.Errorf("no PROJ sample for %s: %s", conv, reason)
	}
}

// TestSeriesEdgeBudget records the largest cs2cs gap of the series methods
// at 5% and 95% of the area of use. The per-point check lives in TestProjSamples.
func TestSeriesEdgeBudget(t *testing.T) {
	need := []string{
		"cassini_soldner",
		"hyperbolic_cassini_soldner",
		"hotine_oblique_mercator_a",
		"laborde_oblique_mercator",
		"new_zealand_map_grid",
		"american_polyconic",
	}
	maxTol := map[string]float64{}
	n := map[string]int{}
	for _, s := range projSamples {
		if !s.edge {
			continue
		}
		n[s.conv]++
		if s.tol > maxTol[s.conv] {
			maxTol[s.conv] = s.tol
		}
	}
	for _, conv := range need {
		if n[conv] == 0 {
			if reason := projUncompared[conv]; reason != "" {
				t.Errorf("%s: %s", conv, reason)
				continue
			}
			t.Errorf("no edge sample for %s", conv)
			continue
		}
		t.Logf("%s edge budget %.4g m (%d points)", conv, maxTol[conv], n[conv])
	}
}
