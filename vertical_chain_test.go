package crs

import (
	"math"
	"strings"
	"testing"
)

// Catania 1965 height is not a geoid. EPSG publishes ETRS89 → Genoa 1942 height
// through the ITAL2005 grid, then Genoa → Catania as a vertical offset of
// +0.141 m (EPSG:9726, accuracy 0.01 m). The offset is stored negated on the
// Catania datum because ToTarget runs back toward Genoa.

func TestVerticalChainGraph(t *testing.T) {
	t.Parallel()

	offsets := []struct {
		datum  string
		target string
		dh     float64
		acc    float64
	}{
		{"catania1965", "genoa1942", -0.141, 0.01},
		{"etrs89_catania1965", "etrs89_genoa1942", -0.141, 0.01},
		{"kiunga", "egm2008geoid", 3, 0.11},
		{"kumul34", "egm96geoid", 0.87, 1},
		{"portmoresby2008", "egm2008geoid", 0.93, 0.11},
		{"portmoresby1996", "egm96geoid", 1.58, 1},
	}
	for _, c := range offsets {
		t.Run(c.datum, func(t *testing.T) {
			d, err := loadDatum(c.datum)
			if err != nil {
				t.Fatal(err)
			}
			if len(d.Transformations) != 1 {
				t.Fatalf("%s has %d operations, want the offset only", c.datum, len(d.Transformations))
			}
			tr := d.Transformations[0]
			off, ok := tr.Operation.(VerticalOffset)
			if !ok {
				t.Fatalf("%s operation is %T, want a vertical offset", c.datum, tr.Operation)
			}
			if tr.Target == nil || tr.Target.Name != c.target {
				got := ""
				if tr.Target != nil {
					got = tr.Target.Name
				}
				t.Fatalf("target %q, want %q", got, c.target)
			}
			if math.Abs(off.Dh-c.dh) > 1e-12 {
				t.Fatalf("dh %v, want %v", off.Dh, c.dh)
			}
			if math.Abs(tr.Accuracy-c.acc) > 1e-12 {
				t.Fatalf("accuracy %v, want %v", tr.Accuracy, c.acc)
			}
		})
	}

	for _, name := range []string{"genoa1942", "etrs89_genoa1942"} {
		t.Run(name, func(t *testing.T) {
			d, err := loadDatum(name)
			if err != nil {
				t.Fatal(err)
			}
			var grids int
			for _, tr := range d.Transformations {
				if _, ok := tr.Operation.(VerticalOffset); ok {
					t.Fatalf("%s carries a vertical offset; that step belongs to Catania", name)
				}
				if g, ok := tr.Operation.(VerticalGrid); ok {
					grids++
					if !strings.Contains(string(g), "geo_igm_mar06.grd") {
						t.Fatalf("grid %s", g)
					}
				}
			}
			if grids != 1 {
				t.Fatalf("%s has %d ITAL2005 grids, want 1", name, grids)
			}
		})
	}
}

func TestVerticalOffsetHeight(t *testing.T) {
	t.Parallel()

	// cs2cs EPSG:9723 EPSG:9724 of (lat 37.5206217, lon 13.63490925, h 100)
	// is height 100.141. The same translation is EPSG:9726.
	sicily := [][2]float64{
		{13.63490925, 37.5206217},
		{13.03081594, 37.24498089},
		{13.89407518, 37.55446393},
		{14.3459583, 37.43302753},
		{12.5, 36.7},
		{15.5, 38.2},
	}
	heights := []float64{0, -5, 100, 250, 1234.5}

	chains := []struct {
		name     string
		from, to int
		dh       float64
		lonLat   [][2]float64
	}{
		{"genoa_to_catania_compound", 9723, 9724, 0.141, sicily},
		{"genoa_to_catania_vertical", 5214, 9721, 0.141, sicily},
		{"egm2008_to_kiunga", 3855, 7652, -3, [][2]float64{{142, -7}, {141, -9}, {143.5, -5.2}}},
		{"egm96_to_kumul34", 5773, 7651, -0.87, [][2]float64{{143, -7}, {142.5, -6}, {144.5, -5.7}}},
		{"egm2008_to_portmoresby2008", 3855, 7841, -0.93, [][2]float64{{147, -9}, {145, -10}, {149, -7}}},
		{"egm96_to_portmoresby1996", 5773, 7832, -1.58, [][2]float64{{147, -9}, {145, -10}, {149, -7}}},
	}

	for _, c := range chains {
		t.Run(c.name, func(t *testing.T) {
			fwd := makeTransform(t, c.from, c.to, 0)
			rev := makeTransform(t, c.to, c.from, 0)
			for _, ll := range c.lonLat {
				for _, h := range heights {
					lon, lat := ll[0], ll[1]
					gotLon, gotLat, gotH, err := fwd(lon, lat, h)
					if err != nil {
						t.Fatalf("forward (%.4f, %.4f, %v): %v", lon, lat, h, err)
					}
					if math.Abs(gotLon-lon) > 1e-12 || math.Abs(gotLat-lat) > 1e-12 {
						t.Fatalf("horizontal moved to (%.8f, %.8f)", gotLon, gotLat)
					}
					if math.Abs(gotH-(h+c.dh)) > 1e-9 {
						t.Fatalf("height %v %+v m: got %v", h, c.dh, gotH)
					}
					backLon, backLat, backH, err := rev(gotLon, gotLat, gotH)
					if err != nil {
						t.Fatal(err)
					}
					if math.Abs(backLon-lon) > 1e-12 || math.Abs(backLat-lat) > 1e-12 || math.Abs(backH-h) > 1e-9 {
						t.Fatalf("round trip got (%.8f, %.8f, %v)", backLon, backLat, backH)
					}
				}
			}
		})
	}
}

func TestEllipsoidToCataniaNeedsGeoid(t *testing.T) {
	t.Parallel()

	// PROJ without geo_igm_mar06.grd leaves the height unchanged and marks the
	// operation ballpark. The library has no readable grid, so the transform
	// fails instead of applying the 0.141 m offset on its own.
	const lon, lat, h = 13.63490925, 37.5206217, 100.0
	blocked := []struct{ from, to int }{
		{4937, 9724},
		{9724, 4937},
		{4937, 9723},
		{9723, 4937},
		{4979, 9721},
		{9721, 4979},
		{4979, 5214},
		{5214, 4979},
	}
	for _, c := range blocked {
		_, err := datumHops(c.from, c.to, lon, lat, 0)
		if err == nil {
			t.Errorf("EPSG:%d→%d found a datum path; the ITAL2005 geoid is not a readable grid", c.from, c.to)
		}
		f := makeTransform(t, c.from, c.to, 0)
		_, _, gotH, err := f(lon, lat, h)
		if err == nil {
			t.Errorf("EPSG:%d→%d height %v; want an error, not a ballpark or the bare offset", c.from, c.to, gotH)
		}
	}
}

func TestCataniaOffsetStaysInsideItsArea(t *testing.T) {
	t.Parallel()

	// Genoa itself is inside the Genoa height area and outside the Catania area.
	f := makeTransform(t, 9723, 9724, 0)
	_, _, _, err := f(8.9, 44.4, 100)
	if err == nil {
		t.Fatal("Genoa city was accepted by the Catania height area")
	}
	if !strings.Contains(err.Error(), "out of bounds") {
		t.Fatalf("got %v", err)
	}
}
