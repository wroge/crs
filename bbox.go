package crs

import (
	"fmt"
	"math"
)

var World = BoundingBox{
	MinLon: -180,
	MinLat: -90,
	MaxLon: 180,
	MaxLat: 90,
}

type BoundingBox struct {
	MinLon float64
	MinLat float64
	MaxLon float64
	MaxLat float64
}

func (b BoundingBox) String() string {
	return fmt.Sprintf("%s,%s,%s,%s", formatFloat(b.MinLon), formatFloat(b.MinLat), formatFloat(b.MaxLon), formatFloat(b.MaxLat))
}

func (b BoundingBox) Contains(lon, lat float64) bool {
	if b.MaxLat < b.MinLat {
		return false
	}

	if lat < b.MinLat || lat > b.MaxLat {
		return false
	}

	if b.MinLon <= b.MaxLon {
		return lon >= b.MinLon && lon <= b.MaxLon
	}

	return lon >= b.MinLon || lon <= b.MaxLon
}

func (b BoundingBox) Area() float64 {
	lonSpan := b.MaxLon - b.MinLon
	if lonSpan < 0 {
		lonSpan = (180 - b.MinLon) + (b.MaxLon + 180)
	}

	latSpan := b.MaxLat - b.MinLat
	if lonSpan <= 0 || latSpan <= 0 {
		return math.MaxFloat64
	}

	return lonSpan * latSpan
}

func lonIntervals(minLon, maxLon float64) [][2]float64 {
	if minLon <= maxLon {
		return [][2]float64{{minLon, maxLon}}
	}

	return [][2]float64{{minLon, 180}, {-180, maxLon}}
}

func intervalsOverlap(a, b [2]float64) bool {
	return a[0] <= b[1] && b[0] <= a[1]
}

func (b BoundingBox) Intersects(bbox BoundingBox) bool {
	if b.MaxLat < b.MinLat || bbox.MaxLat < bbox.MinLat {
		return false
	}

	if b.MaxLat < bbox.MinLat || bbox.MaxLat < b.MinLat {
		return false
	}

	for _, bi := range lonIntervals(b.MinLon, b.MaxLon) {
		for _, bj := range lonIntervals(bbox.MinLon, bbox.MaxLon) {
			if intervalsOverlap(bi, bj) {
				return true
			}
		}
	}

	return false
}
