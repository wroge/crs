// Package crs provides coordinate reference systems, map projections, and
// datum transformations.
package crs

import (
	"fmt"
)

type Func func(float64, float64, float64) (float64, float64, float64, error)

func (f Func) Round(dec int) Func {
	return f.RoundAxis(dec, dec, dec)
}

func (f Func) RoundAxis(decA, decB, decC int) Func {
	return func(a, b, c float64) (float64, float64, float64, error) {
		a, b, c, err := f(a, b, c)

		return round(a, decA), round(b, decB), round(c, decC), err
	}
}

type Conversion interface {
	ToGeographic(s Spheroid, a, b, c float64) (lon, lat, h float64)
	FromGeographic(s Spheroid, lon, lat, h float64) (a, b, c float64)
}

type CoordinateReferenceSystem struct {
	Conversion  Conversion
	Datum       Datum
	BoundingBox BoundingBox
}

func (crs CoordinateReferenceSystem) String() string {
	return build("").addAll(
		"conversion", crs.Conversion,
		"bbox", crs.BoundingBox,
		"datum", crs.Datum.Name,
	).String()
}

func (crs CoordinateReferenceSystem) TransformTo(to CoordinateReferenceSystem) (Func, error) {
	var err error

	crs, err = crs.Intersects(to.BoundingBox)
	if err != nil {
		return nil, err
	}

	to, err = to.Intersects(crs.BoundingBox)
	if err != nil {
		return nil, err
	}

	return func(a, b, c float64) (float64, float64, float64, error) {
		lon, lat, h := crs.Conversion.ToGeographic(crs.Datum.Spheroid, a, b, c)

		// EPSG may publish different areas of use for geographic vs projected CRS
		// that share a datum (e.g. Guam 1963 geographic vs Yap Islands). Same-datum
		// hops are allowed when the point lies in either CRS extent.
		inFrom := crs.BoundingBox.Contains(lon, lat)
		inTo := to.BoundingBox.Contains(lon, lat)
		if crs.Datum.Name == to.Datum.Name {
			if !inFrom && !inTo {
				return 0, 0, 0, fmt.Errorf("out of bounds: [%f,%f]", lon, lat)
			}
		} else if !inFrom || !inTo {
			return 0, 0, 0, fmt.Errorf("out of bounds: [%f,%f]", lon, lat)
		}

		if crs.Datum.Name != to.Datum.Name {
			var err error

			lon, lat, h, err = crs.Datum.TransformTo(to.Datum, lon, lat, h)
			if err != nil {
				return 0, 0, 0, err
			}
		}

		a1, b1, c1 := to.Conversion.FromGeographic(to.Datum.Spheroid, lon, lat, h)

		return a1, b1, c1, nil
	}, nil
}

func (crs CoordinateReferenceSystem) Intersects(bbox ...BoundingBox) (CoordinateReferenceSystem, error) {
	datum, err := crs.Datum.Intersects(bbox...)
	if err != nil {
		return crs, err
	}

	crs.Datum = datum

	return crs, nil
}

// AtEpoch evaluates time-dependent datum operations at the given coordinate epoch
// (decimal year), mirroring Intersects as a pre-Transform datum rewrite.
// Time-specific Helmerts are only kept when epoch matches their transformation epoch.
func (crs CoordinateReferenceSystem) AtEpoch(epoch float64) CoordinateReferenceSystem {
	crs.Datum = crs.Datum.AtEpoch(epoch)

	return crs
}

type CRS interface {
	CoordinateReferenceSystem | int | string
}

func Load[C CRS](crs C, intersects ...BoundingBox) (CoordinateReferenceSystem, error) {
	switch c := any(crs).(type) {
	case int:
		out, err := loadEPSG(c)
		if err != nil {
			return CoordinateReferenceSystem{}, err
		}

		return out.Intersects(intersects...)
	case string:
		out, err := parseCoordinateReferenceSystem(c)
		if err != nil {
			return CoordinateReferenceSystem{}, err
		}

		return out.Intersects(intersects...)
	case CoordinateReferenceSystem:
		return c.Intersects(intersects...)
	}

	return CoordinateReferenceSystem{}, fmt.Errorf("invalid crs")
}

func Transform[F, T CRS](from F, to T, intersects ...BoundingBox) (Func, error) {
	fromCRS, err := Load(from, intersects...)
	if err != nil {
		return nil, err
	}

	toCRS, err := Load(to, intersects...)
	if err != nil {
		return nil, err
	}

	return fromCRS.TransformTo(toCRS)
}

// TransformAt is Transform with both datums evaluated at epoch (decimal year).
func TransformAt[F, T CRS](from F, to T, epoch float64, intersects ...BoundingBox) (Func, error) {
	fromCRS, err := Load(from, intersects...)
	if err != nil {
		return nil, err
	}

	toCRS, err := Load(to, intersects...)
	if err != nil {
		return nil, err
	}

	return fromCRS.AtEpoch(epoch).TransformTo(toCRS.AtEpoch(epoch))
}
