# crs

[![Coverage](https://wroge.github.io/crs/badges/coverage.svg)](https://wroge.github.io/crs/)
[![go.dev reference](https://img.shields.io/badge/go.dev-reference-007d9c?logo=go&logoColor=white)](https://pkg.go.dev/github.com/wroge/crs@v0.0.3)

```sh
go get github.com/wroge/crs@v0.0.3
```

Go library for coordinate reference systems, map projections, and datum
transformations. Supports the majority of EPSG codes (thousands of embedded
definitions).

```go
transform, err := crs.Transform(4326, 32632)
if err != nil {
	log.Fatal(err)
}

east, north, h, err := transform(10.0, 50.0, 0)
```

Conventions: geographic `lon, lat, h` (degrees and longitude relative to the CRS
prime meridian), projected `east, north, h`, geocentric `x, y, z`.

## Geodetic grids

Some datum (nad83, osgb36, ...) transforms need grid files. Provide them via CDN and/or local sources:

```go
crs.SetGridCDN("https://cdn.proj.org/")             // download on demand
crs.SetGridCDN("https://cdn.proj.org/", "./.cache") // optional disk cache

crs.RegisterGridDir("/path/to/grids") // search a directory
crs.RegisterGridFS("", embedFS)       // or an fs.FS (e.g. embed.FS)
```

Lookup order: registered directories, then filesystems, then CDN.

## Time-dependent transforms

This also enables the use of time_specific_* transformations.

```go
transform, err := crs.TransformAt(9989, 9069, 2010.0) // coordinate epoch (decimal year)

// etrs89: accuracy=0.105 bbox=3.34,54.36,8.88,55.92 operation=time_specific_position_vector tx=0.054 ty=0.051 tz=-0.085 rx=0.0021 ry=0.012600000000000002 rz=-0.0204 ds=0.0025 epoch=2014.81 inverse
transform, err := crs.TransformAt(4326, 4258, 2014.81) // regional time_specific transformation
```

## Development

This library succeeds [github.com/wroge/wgs84](https://github.com/wroge/wgs84).

Much of the conversion and transformation math comes from *EPSG Guidance Note 7, part 2* (Coordinate Conversions and Transformations including Formulas, revised September 2019) implemented by hand with help from AI tools. Where results diverged from [PROJ](https://proj.org/), those differences were often found (and sometimes fixed) with the same tooling.

I am documenting that explicitly: AI was used heavily while building this. That does not replace checking the code. Implementations are compared against each other and there is an extensive test suite aimed at catching disagreements with PROJ.

This package is not a full PROJ drop-in, and PROJ still covers more ground. For many Go applications, a self-contained implementation like this is enough.

## License

Licensed under the [MIT License](LICENSE).
