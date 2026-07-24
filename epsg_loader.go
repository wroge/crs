package crs

import (
	"embed"
	"fmt"
	"io"
	"strconv"
	"strings"
	"sync"
)

//go:embed epsg/*.txt
var epsgDir embed.FS

var epsgStore sync.Map

func RegisterEPSG(code int, c CoordinateReferenceSystem) {
	epsgStore.Store(code, c)
}

func loadEPSG(code int) (c CoordinateReferenceSystem, err error) {
	v, ok := epsgStore.Load(code)
	if ok {
		return v.(CoordinateReferenceSystem), nil
	}

	defer func() {
		if err == nil {
			epsgStore.Store(code, c)
		}
	}()

	file, err := epsgDir.Open(fmt.Sprintf("epsg/%d.txt", code))
	if err != nil {
		return c, fmt.Errorf("epsg not found: %d", code)
	}

	data, err := io.ReadAll(file)
	if err != nil {
		return c, err
	}

	return Parse(string(data))
}

func parseBoundingBox(s string) (BoundingBox, bool) {
	parts := strings.Split(s, ",")
	if len(parts) != 4 {
		return BoundingBox{}, false
	}

	vals := make([]float64, 4)
	for i, p := range parts {
		v, err := strconv.ParseFloat(strings.TrimSpace(p), 64)
		if err != nil {
			return BoundingBox{}, false
		}
		vals[i] = v
	}

	return BoundingBox{
		MinLon: vals[0],
		MinLat: vals[1],
		MaxLon: vals[2],
		MaxLat: vals[3],
	}, true
}

type epsgParams struct {
	str map[string]string
	f   map[string]float64
}

func newEPSGParams(fields map[string]string) epsgParams {
	p := epsgParams{
		str: fields,
		f:   make(map[string]float64),
	}

	for k, v := range fields {
		if fv, err := strconv.ParseFloat(v, 64); err == nil {
			p.f[k] = fv
		}
	}

	return p
}

func (p epsgParams) float(keys ...string) float64 {
	for _, k := range keys {
		if v, ok := p.f[k]; ok {
			return v
		}
	}

	return 0
}

func buildConversion(name string, fields map[string]string) (Conversion, error) {
	p := newEPSGParams(fields)

	switch name {
	case "geographic":
		return Geographic{}, nil
	case "geocentric":
		return Geocentric{}, nil
	case "transverse_mercator":
		scale := p.float("scale")
		if scale == 0 {
			scale = 1
		}

		return TransverseMercator{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Scale: scale,
			Eastf: p.float("eastf"), Northf: p.float("northf"),
			ZoneWidth: p.float("zone_width"),
		}, nil
	case "utm":
		zone := int(p.float("zone"))
		if zone < 1 || zone > 60 {
			return nil, fmt.Errorf("utm zone must be 1..60, got %d", zone)
		}

		northf := 0.0
		if _, ok := p.str["southern"]; ok {
			northf = 1e7
		}

		return TransverseMercator{
			Lonf: float64(zone)*6 - 183, Scale: 0.9996,
			Eastf: 500000, Northf: northf,
		}, nil
	case "web_mercator":
		return WebMercator{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Scale: p.float("scale"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "lambert_conformal_conic":
		scale := p.float("scale")
		if scale == 0 {
			scale = 1
		}

		return LambertConformalConic{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Scale: scale,
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "lambert_conformal_conic_1sp_variant_b":
		scale := p.float("scale")
		if scale == 0 {
			scale = 1
		}

		return LambertConformalConic1SPVariantB{
			Lonf: p.float("lonf"), Lat0: p.float("lat0"), Latf: p.float("latf"), Scale: scale,
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "lambert_conformal_conic_2sp":
		return LambertConformalConic2SP{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Sp1: p.float("sp1"), Sp2: p.float("sp2"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "lambert_conformal_conic_2sp_michigan":
		return LambertConformalConic2SPMichigan{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Sp1: p.float("sp1"), Sp2: p.float("sp2"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "lambert_conformal_conic_2sp_belgium":
		return LambertConformalConic2SPBelgium{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Sp1: p.float("sp1"), Sp2: p.float("sp2"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "lambert_conic_near_conformal":
		scale := p.float("scale")
		if scale == 0 {
			scale = 1
		}

		return LambertConicNearConformal{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Scale: scale,
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "lambert_azimuthal_equal_area":
		return LambertAzimuthalEqualArea{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "lambert_azimuthal_equal_area_spherical":
		return LambertAzimuthalEqualAreaSpherical{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "lambert_cylindrical_equal_area":
		return LambertCylindricalEqualArea{
			Lonf: p.float("lonf"), Sp1: p.float("sp1", "latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "lambert_cylindrical_equal_area_spherical":
		return LambertCylindricalEqualAreaSpherical{
			Lonf: p.float("lonf"), Sp1: p.float("sp1", "latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "albers_equal_area":
		return AlbersEqualArea{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Sp1: p.float("sp1"), Sp2: p.float("sp2"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "krovak":
		return Krovak{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Alpha: p.float("alpha"), Scale: p.float("scale"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "krovak_modified", "krovak_modified_north_orientated": // alias: always east,north
		return KrovakModified{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Alpha: p.float("alpha"), Scale: p.float("scale"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "hotine_oblique_mercator_a":
		return HotineObliqueMercatorA{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Alpha: p.float("azimuth", "alpha"), Gamma: p.float("gamma"),
			Scale: p.float("scale"), Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "swiss_oblique_mercator":
		return SwissObliqueMercator{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Alpha: p.float("azimuth", "alpha"), Gamma: p.float("gamma"),
			Scale: p.float("scale"), Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "laborde_oblique_mercator":
		return LabordeObliqueMercator{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Alpha: p.float("azimuth", "alpha"),
			Scale: p.float("scale"), Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "cassini_soldner":
		return CassiniSoldner{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "hyperbolic_cassini_soldner":
		return HyperbolicCassiniSoldner{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "polar_stereographic_a":
		return PolarStereographicA{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Scale: p.float("scale"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "polar_stereographic_b":
		return PolarStereographicB{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "polar_stereographic_c":
		return PolarStereographicC{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "oblique_stereographic":
		return ObliqueStereographic{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Scale: p.float("scale"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "mercator_a":
		return MercatorA{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Scale: p.float("scale"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "mercator_b":
		return MercatorB{
			Lonf: p.float("lonf"), Sp1: p.float("sp1"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "azimuthal_equidistant":
		return AzimuthalEquidistant{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "modified_azimuthal_equidistant":
		return ModifiedAzimuthalEquidistant{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "american_polyconic":
		return AmericanPolyconic{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "bonne_south_orientated":
		return BonneSouthOrientated{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "equal_earth":
		return EqualEarth{
			Lonf:  p.float("lonf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "equidistant_cylindrical":
		return EquidistantCylindrical{
			Lonf: p.float("lonf"), Latf: p.float("latf"), Sp1: p.float("sp1", "latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "colombia_urban":
		return ColombiaUrban{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"), H0: p.float("h0"),
		}, nil
	case "guam_projection":
		return GuamProjection{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "local_orthographic":
		scale := p.float("scale")
		if scale == 0 {
			scale = 1
		}

		return LocalOrthographic{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Azimuth: p.float("azimuth", "alpha"), Scale: scale,
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "tunisia_mining_grid":
		return TunisiaMiningGrid{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	case "new_zealand_map_grid":
		return NewZealandMapGrid{
			Lonf: p.float("lonf"), Latf: p.float("latf"),
			Eastf: p.float("eastf"), Northf: p.float("northf"),
		}, nil
	}

	return nil, fmt.Errorf("unknown conversion: %s", name)
}
