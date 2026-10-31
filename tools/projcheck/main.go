// Command projcheck writes proj_sample_test.go.
//
// For every conversion it asks PROJ (cs2cs) for a same-datum geographic →
// projected sample inside the area of use: one ordinary code, one southern
// code, and one antimeridian code, when such a code exists. Series methods
// are also sampled at 5% and 95% of the area of use; the measured gap is the
// regression budget. EPSG:32600 and EPSG:32700 are compared with the
// per-zone CRS plus the zone prefix, because this PROJ build does not
// implement the zoned grid method.
//
//	go run ./tools/projcheck
package main

import (
	"bytes"
	"context"
	"encoding/json"
	"fmt"
	"go/format"
	"math"
	"os"
	"os/exec"
	"path/filepath"
	"sort"
	"strconv"
	"strings"
	"time"

	"github.com/wroge/crs"
)

const (
	interiorTolM   = 0.01
	interiorTolDeg = 2e-5
	acceptM        = 0.01
	acceptDeg      = 2e-5
)

var seriesMethods = map[string]bool{
	"cassini_soldner":            true,
	"hyperbolic_cassini_soldner": true,
	"hotine_oblique_mercator_a":  true,
	"laborde_oblique_mercator":   true,
	"new_zealand_map_grid":       true,
	"american_polyconic":         true,
}

var required = []string{
	"albers_equal_area",
	"american_polyconic",
	"azimuthal_equidistant",
	"bonne_south_orientated",
	"cassini_soldner",
	"colombia_urban",
	"equal_earth",
	"equidistant_cylindrical",
	"geocentric",
	"geographic",
	"guam_projection",
	"hotine_oblique_mercator_a",
	"hyperbolic_cassini_soldner",
	"krovak",
	"krovak_modified",
	"laborde_oblique_mercator",
	"lambert_azimuthal_equal_area",
	"lambert_azimuthal_equal_area_spherical",
	"lambert_conformal_conic",
	"lambert_conformal_conic_1sp_variant_b",
	"lambert_conformal_conic_2sp",
	"lambert_conformal_conic_2sp_belgium",
	"lambert_conformal_conic_2sp_michigan",
	"lambert_conic_near_conformal",
	"lambert_cylindrical_equal_area",
	"lambert_cylindrical_equal_area_spherical",
	"local_orthographic",
	"mercator_a",
	"mercator_b",
	"modified_azimuthal_equidistant",
	"new_zealand_map_grid",
	"oblique_stereographic",
	"polar_stereographic_a",
	"polar_stereographic_b",
	"polar_stereographic_c",
	"swiss_oblique_mercator",
	"transverse_mercator",
	"tunisia_mining_grid",
	"web_mercator",
}

type rec struct {
	code  int
	conv  string
	datum string
	unit  string
	bbox  crs.BoundingBox
	zoned bool
}

type sample struct {
	name                string
	conv                string
	from, to            int
	inA, inB            float64
	wantA, wantB, wantC float64
	tol                 float64
	deg                 bool
	edge                bool
	measured            float64
}

type axes struct {
	name    [3]string
	dir     [3]string
	factor  [3]float64
	angular [3]bool
	n       int
}

func main() {
	root, err := os.Getwd()
	if err != nil {
		fatal(err)
	}
	byConv, geoByDatum, err := loadDefs(filepath.Join(root, "epsg"))
	if err != nil {
		fatal(err)
	}

	ver := projVersion()
	var samples []sample
	uncompared := map[string]string{}

	for _, conv := range required {
		list := byConv[conv]
		if len(list) == 0 {
			uncompared[conv] = "no EPSG file"
			continue
		}
		var got []sample
		var last string
		if reason := unsupported(list); reason != "" {
			uncompared[conv] = reason
			fmt.Fprintf(os.Stderr, "%s: %s\n", conv, reason)
			continue
		}
		for _, kind := range []string{"normal", "south", "wrap"} {
			cands := ranked(list, kind)
			s, reason := firstGood(cands, geoByDatum, conv, kind)
			if reason != "" {
				last = reason
			}
			got = append(got, s...)
		}
		if seriesMethods[conv] {
			cands := ranked(list, "largest")
			s, reason := firstEdge(cands, geoByDatum, conv)
			if reason != "" && len(s) == 0 {
				last = reason
			}
			got = append(got, s...)
		}
		if len(got) == 0 {
			s, reason := firstEdge(ranked(list, "largest"), geoByDatum, conv)
			if reason != "" {
				last = reason
			}
			got = append(got, s...)
		}
		if len(got) == 0 {
			if last == "" {
				last = "no same-datum geographic CRS that cs2cs could project"
			}
			uncompared[conv] = last
			fmt.Fprintf(os.Stderr, "%s: %s\n", conv, last)
			continue
		}
		samples = append(samples, got...)
		fmt.Fprintf(os.Stderr, "%s: %d samples\n", conv, len(got))
	}

	zoned, zerr := zonedSamples()
	if zerr != "" {
		fmt.Fprintf(os.Stderr, "zoned grid: %s\n", zerr)
	}
	samples = append(samples, zoned...)

	if err := writeTest(filepath.Join(root, "proj_sample_data_test.go"), ver, samples, uncompared); err != nil {
		fatal(err)
	}
	fmt.Fprintf(os.Stderr, "wrote %d samples, %d uncompared\n", len(samples), len(uncompared))
}

func fatal(err error) {
	fmt.Fprintln(os.Stderr, err)
	os.Exit(1)
}

func loadDefs(dir string) (map[string][]rec, map[string][]rec, error) {
	entries, err := os.ReadDir(dir)
	if err != nil {
		return nil, nil, err
	}
	byConv := map[string][]rec{}
	geo := map[string][]rec{}
	for _, e := range entries {
		if !strings.HasSuffix(e.Name(), ".txt") {
			continue
		}
		code, err := strconv.Atoi(strings.TrimSuffix(e.Name(), ".txt"))
		if err != nil {
			continue
		}
		b, err := os.ReadFile(filepath.Join(dir, e.Name()))
		if err != nil {
			return nil, nil, err
		}
		txt := string(b)
		r := rec{
			code:  code,
			conv:  field(txt, "conversion"),
			datum: field(txt, "datum"),
			unit:  field(txt, "unit"),
			zoned: strings.Contains(txt, "zone_width="),
		}
		if r.conv == "" || r.datum == "" {
			continue
		}
		r.bbox, err = parseBBox(field(txt, "bbox"))
		if err != nil {
			continue
		}
		byConv[r.conv] = append(byConv[r.conv], r)
		if r.conv == "geographic" {
			geo[r.datum] = append(geo[r.datum], r)
		}
	}
	return byConv, geo, nil
}

func field(txt, key string) string {
	line, _, _ := strings.Cut(txt, "\n")
	for _, tok := range strings.Fields(line) {
		k, v, ok := strings.Cut(tok, "=")
		if ok && k == key {
			return v
		}
	}
	return ""
}

func parseBBox(s string) (crs.BoundingBox, error) {
	var b crs.BoundingBox
	_, err := fmt.Sscanf(s, "%f,%f,%f,%f", &b.MinLon, &b.MinLat, &b.MaxLon, &b.MaxLat)
	return b, err
}

func unsupported(list []rec) string {
	if len(list) == 0 {
		return ""
	}
	cmd := exec.Command("projinfo", "-o", "PROJ", "-q", fmt.Sprintf("EPSG:%d", list[0].code))
	out, _ := cmd.CombinedOutput()
	text := string(out)
	const marker = "Unsupported conversion method:"
	i := strings.Index(text, marker)
	if i < 0 {
		return ""
	}
	line, _, _ := strings.Cut(strings.TrimSpace(text[i:]), "\n")
	return "PROJ 9.8 " + line
}

func ranked(list []rec, kind string) []rec {
	var out []rec
	for _, r := range list {
		if r.zoned {
			continue
		}
		south := midLat(r.bbox) < 0
		wrap := r.bbox.MinLon > r.bbox.MaxLon
		switch kind {
		case "south":
			if !south || wrap {
				continue
			}
		case "wrap":
			if !wrap {
				continue
			}
		case "normal":
			if south || wrap {
				continue
			}
		case "largest":
		default:
			continue
		}
		out = append(out, r)
	}
	sort.Slice(out, func(i, j int) bool {
		si, sj := score(out[i]), score(out[j])
		if kind == "largest" {
			ai, aj := area(out[i].bbox), area(out[j].bbox)
			if ai != aj {
				return ai > aj
			}
		}
		if si != sj {
			return si > sj
		}
		return out[i].code < out[j].code
	})
	if len(out) > 12 {
		out = out[:12]
	}
	return out
}

func score(r rec) int {
	s := 0
	if r.unit == "m" || r.unit == "degree" {
		s += 4
	}
	if r.datum == "wgs84" {
		s += 2
	}
	if area(r.bbox) < 80*80 {
		s++
	}
	return s
}

func area(b crs.BoundingBox) float64 {
	lat := b.MaxLat - b.MinLat
	lon := b.MaxLon - b.MinLon
	if b.MinLon > b.MaxLon {
		lon = (180 - b.MinLon) + (b.MaxLon + 180)
	}
	if lat < 0 {
		lat = 0
	}
	if lon < 0 {
		lon = 0
	}
	return lat * lon
}

func midLat(b crs.BoundingBox) float64 {
	return (b.MinLat + b.MaxLat) / 2
}

func pickGeo(datum string, geos map[string][]rec) (rec, bool) {
	list := geos[datum]
	if len(list) == 0 {
		return rec{}, false
	}
	if datum == "wgs84" {
		for _, g := range list {
			if g.code == 4326 {
				return g, true
			}
		}
	}
	sort.Slice(list, func(i, j int) bool {
		di := list[i].unit == "degree"
		dj := list[j].unit == "degree"
		if di != dj {
			return di
		}
		return list[i].code < list[j].code
	})
	return list[0], true
}

func firstGood(cands []rec, geos map[string][]rec, conv, kind string) ([]sample, string) {
	var last string
	for _, r := range cands {
		s, reason := projectPair(r, geos, conv, kind, [][2]float64{{0.50, 0.50}, {0.35, 0.65}}, false)
		if reason != "" {
			last = reason
			continue
		}
		if len(s) > 0 {
			return s, ""
		}
	}
	return nil, last
}

func firstEdge(cands []rec, geos map[string][]rec, conv string) ([]sample, string) {
	var last string
	for _, r := range cands {
		s, reason := projectPair(r, geos, conv, "edge", [][2]float64{{0.05, 0.05}, {0.95, 0.95}}, true)
		if reason != "" {
			last = reason
			continue
		}
		if len(s) > 0 {
			return s, ""
		}
	}
	return nil, last
}

func projectPair(r rec, geos map[string][]rec, conv, kind string, fracs [][2]float64, edge bool) ([]sample, string) {
	geo, ok := pickGeo(r.datum, geos)
	if !ok {
		return nil, "no geographic CRS for datum " + r.datum
	}
	if geo.code == r.code {
		return nil, "geographic CRS is the target"
	}
	src, err := crs.Load(geo.code)
	if err != nil {
		return nil, err.Error()
	}
	dst, err := crs.Load(r.code)
	if err != nil {
		return nil, err.Error()
	}
	pm := prime(src.Datum)
	srcIn := src.WithUnit(crs.Degree)
	deg := dst.Unit.IsAngular()
	dstOut := dst.WithUnit(crs.Metre)
	if deg {
		dstOut = dst.WithUnit(crs.Degree)
	}
	fn, err := srcIn.TransformTo(dstOut)
	if err != nil {
		return nil, err.Error()
	}
	srcAxes, err := projAxes(geo.code)
	if err != nil {
		return nil, err.Error()
	}
	dstAxes, err := projAxes(r.code)
	if err != nil {
		return nil, err.Error()
	}
	latFirst := srcAxes.n > 0 && (srcAxes.dir[0] == "north" || srcAxes.dir[0] == "south")

	var out []sample
	for _, fr := range fracs {
		lon, lat, ok := pointAt(r.bbox, fr[0], fr[1])
		if !ok {
			continue
		}
		local := wrap180(lon - pm)
		la, lb, lc, err := fn(local, lat, 0)
		if err != nil {
			return nil, err.Error()
		}
		if !finite(la, lb, lc) {
			continue
		}
		inLon := toAngular(local, src.Unit)
		inLat := toAngular(lat, src.Unit)
		ca, cb, cc, err := runCS2CS(geo.code, r.code, latFirst, inLon, inLat)
		if err != nil {
			return nil, err.Error()
		}
		pa, pb, pc, err := align(dstAxes, ca, cb, cc, deg)
		if err != nil {
			return nil, err.Error()
		}
		da, db := la-pa, lb-pb
		if deg {
			da, db = wrap180(da), wrap180(db)
		}
		errM := math.Hypot(da, db)
		tol := interiorTolM
		limit := acceptM
		if deg {
			tol = interiorTolDeg
			limit = acceptDeg
			errM = math.Max(math.Abs(da), math.Abs(db))
		}
		if !edge && math.Abs(lc-pc) > acceptM {
			fmt.Fprintf(os.Stderr, "  skip %s EPSG:%d %s height err=%.4g\n", conv, r.code, kind, lc-pc)
			continue
		}
		if edge {
			if !finite(errM) {
				continue
			}
			tol = budget(errM)
		} else if errM > limit {
			fmt.Fprintf(os.Stderr, "  skip %s EPSG:%d %s err=%.4g (limit %.4g)\n", conv, r.code, kind, errM, limit)
			continue
		}
		tag := kind
		if edge {
			tag = fmt.Sprintf("edge%.0f", fr[0]*100)
		}
		out = append(out, sample{
			name:     fmt.Sprintf("%s/%d/%s", conv, r.code, tag),
			conv:     conv,
			from:     geo.code,
			to:       r.code,
			inA:      local,
			inB:      lat,
			wantA:    pa,
			wantB:    pb,
			wantC:    pc,
			tol:      tol,
			deg:      deg,
			edge:     edge,
			measured: errM,
		})
		_ = lc
	}
	if len(out) == 0 {
		return nil, fmt.Sprintf("EPSG:%d interior outside %.4g", r.code, acceptM)
	}
	return out, ""
}

func zonedSamples() ([]sample, string) {
	pts := []struct {
		grid int
		lon  float64
		lat  float64
		tag  string
	}{
		{32600, 9, 50, "zone32"},
		{32600, -177, 10, "zone1"},
		{32600, 177, 10, "zone60"},
		{32700, 9, -40, "zone32"},
		{32700, -179, -20, "zone1"},
	}
	var out []sample
	for _, p := range pts {
		src, err := crs.Load(4326)
		if err != nil {
			return nil, err.Error()
		}
		dst, err := crs.Load(p.grid)
		if err != nil {
			return nil, err.Error()
		}
		fn, err := src.WithUnit(crs.Degree).TransformTo(dst.WithUnit(crs.Metre))
		if err != nil {
			return nil, err.Error()
		}
		east, north, h, err := fn(p.lon, p.lat, 0)
		if err != nil {
			return nil, err.Error()
		}
		zone := int(math.Floor(east / 1e6))
		if zone < 1 || zone > 60 {
			return nil, fmt.Sprintf("EPSG:%d lon %.1f gave zone %d", p.grid, p.lon, zone)
		}
		zoneCode := p.grid + zone
		ca, cb, cc, err := runCS2CS(4326, zoneCode, true, p.lon, p.lat)
		if err != nil {
			return nil, err.Error()
		}
		axes, err := projAxes(zoneCode)
		if err != nil {
			return nil, err.Error()
		}
		pa, pb, pc, err := align(axes, ca, cb, cc, false)
		if err != nil {
			return nil, err.Error()
		}
		wantE := float64(zone)*1e6 + pa
		errM := math.Hypot(east-wantE, north-pb)
		if errM > acceptM || math.Abs(h-pc) > acceptM {
			return nil, fmt.Sprintf("EPSG:%d %s prefix err %.4g m (easting lib %.3f proj %.3f)", p.grid, p.tag, errM, east, wantE)
		}
		out = append(out, sample{
			name:     fmt.Sprintf("transverse_mercator/%d/%s", p.grid, p.tag),
			conv:     "transverse_mercator",
			from:     4326,
			to:       p.grid,
			inA:      p.lon,
			inB:      p.lat,
			wantA:    wantE,
			wantB:    pb,
			wantC:    pc,
			tol:      interiorTolM,
			measured: errM,
		})
	}
	return out, ""
}

func pointAt(b crs.BoundingBox, flon, flat float64) (float64, float64, bool) {
	if b.MaxLat < b.MinLat {
		return 0, 0, false
	}
	lat := b.MinLat + flat*(b.MaxLat-b.MinLat)
	if lat > 89.5 {
		lat = 89.5
	}
	if lat < -89.5 {
		lat = -89.5
	}
	var lon float64
	if b.MinLon <= b.MaxLon {
		if b.MaxLon == b.MinLon && b.MaxLat == b.MinLat {
			return 0, 0, false
		}
		lon = b.MinLon + flon*(b.MaxLon-b.MinLon)
	} else {
		span := (180 - b.MinLon) + (b.MaxLon + 180)
		lon = b.MinLon + flon*span
		if lon > 180 {
			lon -= 360
		}
	}
	return lon, lat, true
}

func prime(d crs.Datum) float64 {
	for _, t := range d.Transformations {
		if lr, ok := t.Operation.(crs.LongitudeRotation); ok {
			return lr.Lon
		}
	}
	return 0
}

func toAngular(deg float64, u crs.Unit) float64 {
	if u.Name == "grad" {
		return deg * 200 / 180
	}
	return deg
}

func wrap180(d float64) float64 {
	return math.Remainder(d, 360)
}

func finite(v ...float64) bool {
	for _, x := range v {
		if math.IsNaN(x) || math.IsInf(x, 0) {
			return false
		}
	}
	return true
}

func budget(err float64) float64 {
	if err < 0.005 {
		return 0.01
	}
	padded := err * 1.25
	if padded < 1 {
		return math.Ceil(padded*1000) / 1000
	}
	return math.Ceil(padded*100) / 100
}

var axisCache = map[int]axes{}

func projAxes(code int) (axes, error) {
	if a, ok := axisCache[code]; ok {
		return a, nil
	}
	ctx, cancel := context.WithTimeout(context.Background(), 20*time.Second)
	defer cancel()
	cmd := exec.CommandContext(ctx, "projinfo", "-o", "PROJJSON", "-q", fmt.Sprintf("EPSG:%d", code))
	out, err := cmd.Output()
	if err != nil {
		if ee, ok := err.(*exec.ExitError); ok {
			return axes{}, fmt.Errorf("projinfo EPSG:%d: %s", code, bytes.TrimSpace(ee.Stderr))
		}
		return axes{}, err
	}
	var doc struct {
		CoordinateSystem struct {
			Axis []struct {
				Name      string          `json:"name"`
				Direction string          `json:"direction"`
				Unit      json.RawMessage `json:"unit"`
			} `json:"axis"`
		} `json:"coordinate_system"`
	}
	if err := json.Unmarshal(out, &doc); err != nil {
		return axes{}, fmt.Errorf("projinfo EPSG:%d: %w", code, err)
	}
	var a axes
	for i, ax := range doc.CoordinateSystem.Axis {
		if i >= 3 {
			break
		}
		f, ang, err := unitFactor(ax.Unit)
		if err != nil {
			return axes{}, fmt.Errorf("EPSG:%d axis %s: %w", code, ax.Direction, err)
		}
		a.name[i] = ax.Name
		a.dir[i] = ax.Direction
		a.factor[i] = f
		a.angular[i] = ang
		a.n++
	}
	if a.n < 2 {
		return axes{}, fmt.Errorf("EPSG:%d: %d axes", code, a.n)
	}
	axisCache[code] = a
	return a, nil
}

func unitFactor(raw json.RawMessage) (float64, bool, error) {
	if len(raw) == 0 {
		return 0, false, fmt.Errorf("missing unit")
	}
	var name string
	if raw[0] == '"' {
		if err := json.Unmarshal(raw, &name); err != nil {
			return 0, false, err
		}
	} else {
		var o struct {
			Name   string  `json:"name"`
			Factor float64 `json:"conversion_factor"`
		}
		if err := json.Unmarshal(raw, &o); err != nil {
			return 0, false, err
		}
		name = o.Name
		if o.Factor != 0 {
			return o.Factor, name == "degree" || name == "grad", nil
		}
	}
	switch name {
	case "degree":
		return math.Pi / 180, true, nil
	case "grad":
		return math.Pi / 200, true, nil
	case "metre", "meter":
		return 1, false, nil
	default:
		return 0, false, fmt.Errorf("unsupported unit %q", name)
	}
}

func runCS2CS(from, to int, latFirst bool, lon, lat float64) (float64, float64, float64, error) {
	line := fmt.Sprintf("%.10f %.10f\n", lon, lat)
	if latFirst {
		line = fmt.Sprintf("%.10f %.10f\n", lat, lon)
	}
	ctx, cancel := context.WithTimeout(context.Background(), 30*time.Second)
	defer cancel()
	cmd := exec.CommandContext(ctx, "cs2cs", "-d", "8", fmt.Sprintf("EPSG:%d", from), fmt.Sprintf("EPSG:%d", to))
	cmd.Stdin = strings.NewReader(line)
	out, err := cmd.Output()
	if err != nil {
		msg := ""
		if ee, ok := err.(*exec.ExitError); ok {
			msg = string(bytes.TrimSpace(ee.Stderr))
		}
		return 0, 0, 0, fmt.Errorf("cs2cs EPSG:%d EPSG:%d: %v %s", from, to, err, msg)
	}
	return parseXYZ(string(out))
}

func parseXYZ(s string) (float64, float64, float64, error) {
	f := strings.Fields(s)
	if len(f) < 2 {
		return 0, 0, 0, fmt.Errorf("cs2cs: %q", strings.TrimSpace(s))
	}
	vals := make([]float64, 0, 3)
	for _, tok := range f[:min(3, len(f))] {
		if strings.Contains(tok, "*") || tok == "inf" || tok == "-inf" {
			return 0, 0, 0, fmt.Errorf("cs2cs: %s", strings.TrimSpace(s))
		}
		v, err := strconv.ParseFloat(tok, 64)
		if err != nil {
			return 0, 0, 0, fmt.Errorf("cs2cs: %q", strings.TrimSpace(s))
		}
		vals = append(vals, v)
	}
	h := 0.0
	if len(vals) == 3 {
		h = vals[2]
	}
	return vals[0], vals[1], h, nil
}

func align(a axes, c0, c1, c2 float64, wantDeg bool) (float64, float64, float64, error) {
	v0, err := canon(c0, a.factor[0], a.angular[0], wantDeg)
	if err != nil {
		return 0, 0, 0, err
	}
	v1, err := canon(c1, a.factor[1], a.angular[1], wantDeg)
	if err != nil {
		return 0, 0, 0, err
	}
	h := 0.0
	if a.n > 2 && !a.angular[2] {
		h = c2 * a.factor[2]
	} else if a.n > 2 && wantDeg {
		h = c2
	} else {
		h = c2
	}
	// Library order is east, north or lon, lat. Axis direction is not the
	// order: polar stereographic names its easting "north" because of the
	// meridian of the axis. Names decide the order.
	if librarySwaps(a) {
		return v1, v0, h, nil
	}
	return v0, v1, h, nil
}

func librarySwaps(a axes) bool {
	if a.n < 2 {
		return false
	}
	n0 := strings.ToLower(a.name[0])
	n1 := strings.ToLower(a.name[1])
	if strings.Contains(n0, "latitude") || strings.Contains(n0, "northing") {
		if strings.Contains(n1, "longitude") || strings.Contains(n1, "easting") {
			return true
		}
	}
	if a.angular[0] && (a.dir[0] == "north" || a.dir[0] == "south") {
		return true
	}
	return false
}

func canon(v, factor float64, angular, wantDeg bool) (float64, error) {
	if factor == 0 {
		return 0, fmt.Errorf("zero unit factor")
	}
	if angular {
		if !wantDeg {
			return 0, fmt.Errorf("angular PROJ axis for a linear CRS")
		}
		// factor is radians per CRS unit. Value is in that unit.
		return v * factor / (math.Pi / 180), nil
	}
	if wantDeg {
		return v, nil
	}
	return v * factor, nil
}

func projVersion() string {
	cmd := exec.Command("cs2cs")
	out, _ := cmd.CombinedOutput()
	line, _, _ := strings.Cut(string(out), "\n")
	line = strings.TrimSpace(line)
	if line == "" {
		return "PROJ"
	}
	return line
}

func writeTest(path, ver string, samples []sample, uncompared map[string]string) error {
	var b strings.Builder
	fmt.Fprintf(&b, "// Generated by tools/projcheck (%s). DO NOT EDIT.\n", ver)
	fmt.Fprintf(&b, "// Same-datum geographic to projected samples against cs2cs.\n")
	fmt.Fprintf(&b, "// Interior points use 1 cm or 2e-5 degrees. Edge points of series\n")
	fmt.Fprintf(&b, "// methods use the measured gap, padded, as the regression budget.\n")
	fmt.Fprintf(&b, "// EPSG:32600 and EPSG:32700 are the zone CRS plus the zone prefix.\n")
	fmt.Fprintf(&b, "// go test does not require PROJ.\n")
	edgeMax := map[string]float64{}
	for _, s := range samples {
		if s.edge && s.measured > edgeMax[s.conv] {
			edgeMax[s.conv] = s.measured
		}
	}
	if len(edgeMax) > 0 {
		b.WriteString("//\n// Series-method gap at 5% and 95% of the area of use (metres):\n")
		ek := make([]string, 0, len(edgeMax))
		for k := range edgeMax {
			ek = append(ek, k)
		}
		sort.Strings(ek)
		for _, k := range ek {
			fmt.Fprintf(&b, "//   %s  max %.4g\n", k, edgeMax[k])
		}
	}
	fmt.Fprintf(&b, "package crs\n\n")
	b.WriteString("func init() {\n")
	b.WriteString("\tprojSamples = []projSample{\n")
	for _, s := range samples {
		fmt.Fprintf(&b, "\t\t{name: %q, conv: %q, from: %d, to: %d, inA: %s, inB: %s, wantA: %s, wantB: %s, wantC: %s, tol: %s, deg: %t, edge: %t},\n",
			s.name, s.conv, s.from, s.to,
			num(s.inA), num(s.inB), num(s.wantA), num(s.wantB), num(s.wantC),
			num(s.tol), s.deg, s.edge)
	}
	b.WriteString("\t}\n")
	b.WriteString("\tprojUncompared = map[string]string{\n")
	keys := make([]string, 0, len(uncompared))
	for k := range uncompared {
		keys = append(keys, k)
	}
	sort.Strings(keys)
	for _, k := range keys {
		fmt.Fprintf(&b, "\t\t%q: %q,\n", k, uncompared[k])
	}
	b.WriteString("\t}\n")
	b.WriteString("}\n")

	src, err := format.Source([]byte(b.String()))
	if err != nil {
		return fmt.Errorf("format: %w\n%s", err, b.String())
	}
	return os.WriteFile(path, src, 0o644)
}

func num(v float64) string {
	return strconv.FormatFloat(v, 'g', 12, 64)
}
