package crs

import (
	"fmt"
	"slices"
	"strconv"
	"strings"
)

var (
	datumStartKeys = []string{"accuracy"}
	crsEmbedKeys   = []string{"datum", "spheroid", "accuracy"}
)

func flattenDSL(txt string) string {
	var b strings.Builder

	for line := range strings.SplitSeq(txt, "\n") {
		line = stripComment(line)
		if line == "" {
			continue
		}
		if b.Len() > 0 {
			b.WriteByte(' ')
		}
		b.WriteString(line)
	}

	return b.String()
}

func (pp parts) indexOfKey(keys ...string) int {
	for i, p := range pp {
		if slices.Contains(keys, p.key) {
			return i
		}
	}

	return -1
}

func (pp parts) hasKey(key string) bool {
	return pp.indexOfKey(key) >= 0
}

func (pp parts) segmentByStartKeys(keys ...string) []parts {
	if len(pp) == 0 {
		return nil
	}

	starts := []int{0}
	for i := 1; i < len(pp); i++ {
		if slices.Contains(keys, pp[i].key) {
			starts = append(starts, i)
		}
	}

	out := make([]parts, 0, len(starts))
	for i, start := range starts {
		end := len(pp)
		if i+1 < len(starts) {
			end = starts[i+1]
		}
		out = append(out, pp[start:end])
	}

	return out
}

func parseSpheroid(txt string) (Spheroid, error) {
	return asParts(txt).asSpheroid()
}

func (pp parts) asSpheroid() (Spheroid, error) {
	var s Spheroid

	if name, ok := pp.getPart("spheroid").asString(); ok {
		loaded, err := loadSpheroid(name)
		if err != nil {
			return s, err
		}

		return loaded, nil
	}

	a, hasA := pp.float("a")
	fi, hasFi := pp.float("fi")
	if hasA {
		s.A = a
	}
	if hasFi {
		s.Fi = fi // fi=0 is a sphere (authalic / planetary)
	}

	if !hasA || s.A == 0 || !hasFi {
		return s, fmt.Errorf("spheroid requires spheroid=<name> or a= and fi=")
	}

	return s, nil
}

func parseDatum(txt string) (Datum, error) {
	return asParts(txt).asDatum()
}

func (pp parts) asDatum() (Datum, error) {
	var d Datum

	if len(pp) == 0 {
		return d, fmt.Errorf("empty datum")
	}

	headerEnd := len(pp)
	if i := pp.indexOfKey(datumStartKeys...); i >= 0 {
		headerEnd = i
	}
	header := pp[:headerEnd]

	if name, ok := header.getPart("datum").asString(); ok {
		d.Name = strings.ToLower(name)
	}

	if name, ok := header.getPart("spheroid").asString(); ok {
		s, err := loadSpheroid(name)
		if err != nil {
			return d, err
		}
		d.Spheroid = s
	}

	if a, ok := header.float("a"); ok {
		d.Spheroid.A = a
	}

	if fi, ok := header.float("fi"); ok {
		d.Spheroid.Fi = fi
	}

	if d.Spheroid.A == 0 {
		if d.Name == "" {
			return d, fmt.Errorf("datum missing spheroid (need spheroid=<name> or a= and fi=)")
		}

		loaded, err := loadDatum(d.Name)
		if err != nil {
			return d, err
		}

		if headerEnd == len(pp) {
			return loaded, nil
		}

		d.Spheroid = loaded.Spheroid
	}

	for _, seg := range pp[headerEnd:].segmentByStartKeys(datumStartKeys...) {
		if !seg.hasKey("accuracy") {
			continue
		}

		t, err := seg.asTransformation()
		if err != nil {
			return d, err
		}

		if t.Operation != nil {
			d.Transformations = append(d.Transformations, t)
		}
	}

	sortTransformations(d.Transformations)

	return d, nil
}

func sortTransformations(ts []Transformation) {
	slices.SortFunc(ts, func(a, b Transformation) int {
		if a.Accuracy < b.Accuracy {
			return -1
		}

		if a.Accuracy > b.Accuracy {
			return 1
		}

		if a.BoundingBox.Area() < b.BoundingBox.Area() {
			return -1
		}

		if a.BoundingBox.Area() > b.BoundingBox.Area() {
			return 1
		}

		return 0
	})
}

func parseCoordinateReferenceSystem(txt string) (CoordinateReferenceSystem, error) {
	return asParts(txt).asCRS()
}

func (pp parts) asCRS() (CoordinateReferenceSystem, error) {
	var c CoordinateReferenceSystem
	if len(pp) == 0 {
		return c, fmt.Errorf("empty crs")
	}

	if !pp.hasKey("conversion") {
		return c, fmt.Errorf("missing conversion")
	}

	expanded := pp.hasKey("accuracy")

	if !expanded {
		return pp.asCRSCompact()
	}

	return pp.asCRSExpanded()
}

func (pp parts) asCRSCompact() (CoordinateReferenceSystem, error) {
	var c CoordinateReferenceSystem

	fields := make(map[string]string, len(pp))
	for _, p := range pp {
		if p.key == "" {
			continue
		}
		fields[p.key] = p.value
	}

	csName, ok := fields["conversion"]
	if !ok {
		return c, fmt.Errorf("missing conversion")
	}

	if bboxStr, ok := fields["bbox"]; ok {
		bbox, ok := parseBoundingBox(bboxStr)
		if !ok {
			return c, fmt.Errorf("invalid bbox: %s", bboxStr)
		}
		c.BoundingBox = bbox
	} else {
		c.BoundingBox = World
	}

	cs, err := buildConversion(csName, fields)
	if err != nil {
		return c, err
	}

	c.Conversion = cs

	switch {
	case fields["datum"] != "":
		c.Datum, err = loadDatum(fields["datum"])
		if err != nil {
			return c, fmt.Errorf("datum %s: %w", fields["datum"], err)
		}
	case fields["spheroid"] != "":
		s, err := loadSpheroid(fields["spheroid"])
		if err != nil {
			return c, fmt.Errorf("spheroid %s: %w", fields["spheroid"], err)
		}
		c.Datum = Datum{Spheroid: s}
	default:
		a, okA := parseFloatField(fields, "a")
		fi, okFi := parseFloatField(fields, "fi")

		if !okA || !okFi {
			return c, fmt.Errorf("missing datum or spheroid (need datum=, spheroid=, or a= and fi=)")
		}
		c.Datum = Datum{Spheroid: Spheroid{A: a, Fi: fi}}
	}

	c.Datum, err = c.Datum.Intersects(c.BoundingBox)
	if err != nil {
		return c, fmt.Errorf("datum: %w", err)
	}

	return c, nil
}

func parseFloatField(fields map[string]string, key string) (float64, bool) {
	v, ok := fields[key]
	if !ok || v == "" {
		return 0, false
	}

	f, err := strconv.ParseFloat(v, 64)
	if err != nil {
		return 0, false
	}

	return f, true
}

func (pp parts) asCRSExpanded() (CoordinateReferenceSystem, error) {
	var c CoordinateReferenceSystem

	split := pp.indexOfKey(crsEmbedKeys...)
	if split < 0 {
		return c, fmt.Errorf("expanded crs missing datum/spheroid/accuracy")
	}

	header := pp[:split]
	rest := pp[split:]

	fields := make(map[string]string, len(header))
	for _, p := range header {
		if p.key == "" {
			continue
		}
		fields[p.key] = p.value
	}

	csName, ok := fields["conversion"]
	if !ok {
		return c, fmt.Errorf("missing conversion")
	}

	if bboxStr, ok := fields["bbox"]; ok {
		bbox, ok := parseBoundingBox(bboxStr)
		if !ok {
			return c, fmt.Errorf("invalid bbox: %s", bboxStr)
		}
		c.BoundingBox = bbox
	} else {
		c.BoundingBox = World
	}

	cs, err := buildConversion(csName, fields)
	if err != nil {
		return c, err
	}

	c.Conversion = cs

	c.Datum, err = rest.asDatum()
	if err != nil {
		return c, err
	}

	c.Datum, err = c.Datum.Intersects(c.BoundingBox)
	if err != nil {
		return c, fmt.Errorf("datum: %w", err)
	}

	return c, nil
}
