package crs

type TunisiaMiningGrid struct {
	Lonf, Latf, Eastf, Northf float64
}

func (cs TunisiaMiningGrid) String() string {
	return build("tunisia_mining_grid").addAll(
		"lonf", cs.Lonf,
		"latf", cs.Latf,
		"eastf", cs.Eastf,
		"northf", cs.Northf,
	).String()
}

const tunisiaLonPerKm = 0.012185 // grads per kilometre easting

func tunisiaLatPerKm(aboveOrigin bool) float64 {
	if aboveOrigin {
		return 0.010015
	}

	return 0.01002
}

func (cs TunisiaMiningGrid) FromGeographic(_ Spheroid, lon, lat, h float64) (float64, float64, float64) {
	latG := lat / 0.9
	lonG := lon / 0.9
	lat0G := cs.Latf / 0.9
	lon0G := cs.Lonf / 0.9
	a := tunisiaLatPerKm(latG > lat0G)
	east := cs.Eastf + (lonG-lon0G)/tunisiaLonPerKm*1000
	north := cs.Northf + (latG-lat0G)/a*1000

	return east, north, h
}

func (cs TunisiaMiningGrid) ToGeographic(_ Spheroid, east, north, h float64) (float64, float64, float64) {
	lat0G := cs.Latf / 0.9
	lon0G := cs.Lonf / 0.9
	eKm := (east - cs.Eastf) / 1000
	nKm := (north - cs.Northf) / 1000
	a := tunisiaLatPerKm(nKm > 0)
	lonG := lon0G + eKm*tunisiaLonPerKm
	latG := lat0G + nKm*a

	return lonG * 0.9, latG * 0.9, h
}
