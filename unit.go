package crs

import (
	"fmt"
	"math"
	"strings"
)

// Unit is an axis I/O unit. Projection parameters stay in metres/degrees internally;
// Unit only scales caller coordinates on Transform entry/exit.
//
// For length units, ToSI is metres per unit. For angular units, ToSI is radians per unit.
type Unit struct {
	Name string
	ToSI float64
}

func (u Unit) String() string {
	if u.Name == "" {
		return ""
	}

	return u.Name
}

func (u Unit) IsAngular() bool {
	return u.Name == "degree" || u.Name == "grad"
}

func (u Unit) IsLength() bool {
	return u.Name != "" && !u.IsAngular()
}

var (
	Degree = Unit{Name: "degree", ToSI: math.Pi / 180}
	Grad   = Unit{Name: "grad", ToSI: math.Pi / 200}

	Metre            = Unit{Name: "m", ToSI: 1}
	Kilometre        = Unit{Name: "km", ToSI: 1000}
	Foot             = Unit{Name: "ft", ToSI: 0.3048}
	USSurveyFoot     = Unit{Name: "us-ft", ToSI: 1200.0 / 3937.0}
	GermanLegalMetre = Unit{Name: "german_legal_metre", ToSI: 1.0000135965}
	ClarkeFoot       = Unit{Name: "clarke_foot", ToSI: 0.3047972654}
	GoldCoastFoot    = Unit{Name: "gold_coast_foot", ToSI: 0.3047997101815088}
	SearsFoot        = Unit{Name: "sears_foot", ToSI: 0.3047994715386762}
	BritishFoot1936  = Unit{Name: "british_foot_1936", ToSI: 0.3048007491}
	Link             = Unit{Name: "link", ToSI: 0.201168}
	ClarkeLink       = Unit{Name: "clarke_link", ToSI: 0.201166195164}
	IndianYard       = Unit{Name: "indian_yard", ToSI: 0.9143985307444408}
	SearsYard        = Unit{Name: "sears_yard", ToSI: 0.9143984146160287}
	ClarkeYard       = Unit{Name: "clarke_yard", ToSI: 0.9143917962}
	SearsChain       = Unit{Name: "sears_chain", ToSI: 20.116765121552632}
	SearsChainTrunc  = Unit{Name: "sears_chain_truncated", ToSI: 20.116756}
)

var unitsByName = map[string]Unit{
	Degree.Name:           Degree,
	Grad.Name:             Grad,
	Metre.Name:            Metre,
	"metre":               Metre,
	"meter":               Metre,
	Kilometre.Name:        Kilometre,
	Foot.Name:             Foot,
	"foot":                Foot,
	USSurveyFoot.Name:     USSurveyFoot,
	GermanLegalMetre.Name: GermanLegalMetre,
	ClarkeFoot.Name:       ClarkeFoot,
	GoldCoastFoot.Name:    GoldCoastFoot,
	SearsFoot.Name:        SearsFoot,
	BritishFoot1936.Name:  BritishFoot1936,
	Link.Name:             Link,
	ClarkeLink.Name:       ClarkeLink,
	IndianYard.Name:       IndianYard,
	SearsYard.Name:        SearsYard,
	ClarkeYard.Name:       ClarkeYard,
	SearsChain.Name:       SearsChain,
	SearsChainTrunc.Name:  SearsChainTrunc,
}

func parseUnit(name string) (Unit, error) {
	name = strings.ToLower(strings.TrimSpace(name))
	if name == "" {
		return Unit{}, fmt.Errorf("empty unit")
	}

	if u, ok := unitsByName[name]; ok {
		return u, nil
	}

	return Unit{}, fmt.Errorf("unknown unit: %s", name)
}

func defaultUnit(cs Conversion) Unit {
	switch cs.(type) {
	case Geographic:
		return Degree
	case Geocentric:
		return Metre
	default:
		return Metre
	}
}

// toInternal scales caller axis values into library units (degrees / metres).
func (u Unit) toInternal(a, b, c float64, scaleHeight bool) (float64, float64, float64) {
	if u.Name == "" {
		return a, b, c
	}

	if u.IsAngular() {
		factor := u.ToSI / Degree.ToSI // unit → degrees

		a, b = a*factor, b*factor

		return a, b, c
	}

	a, b = a*u.ToSI, b*u.ToSI

	if scaleHeight {
		c *= u.ToSI
	}

	return a, b, c
}

// fromInternal scales library units (degrees / metres) into caller axis units.
func (u Unit) fromInternal(a, b, c float64, scaleHeight bool) (float64, float64, float64) {
	if u.Name == "" {
		return a, b, c
	}

	if u.IsAngular() {
		factor := Degree.ToSI / u.ToSI // degrees → unit

		a, b = a*factor, b*factor

		return a, b, c
	}

	inv := 1 / u.ToSI

	a, b = a*inv, b*inv

	if scaleHeight {
		c *= inv
	}

	return a, b, c
}
