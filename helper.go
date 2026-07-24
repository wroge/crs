package crs

import (
	"errors"
	"fmt"
	"math"
	"reflect"
	"strconv"
	"strings"
)

func formatFloat(v float64) string {
	return strconv.FormatFloat(v, 'f', -1, 64)
}

func degree(r float64) float64 {
	return r * 180 / math.Pi
}

func radian(g float64) float64 {
	return g * math.Pi / 180
}

func intPow(val float64, times int) float64 {
	result := 1.0

	for range times {
		result *= val
	}

	return result
}

func round(val float64, dec int) float64 {
	factor := math.Pow(10, float64(dec))

	r := math.Round(val*factor) / factor

	if r == -0 {
		return 0
	}

	return r
}

func sin2(r float64) float64 {
	return intPow(math.Sin(r), 2)
}

func sign(x float64) float64 {
	if x < 0 {
		return -1
	}
	return 1
}

func clamp(x, min, max float64) float64 {
	if x < min {
		return min
	}
	if x > max {
		return max
	}
	return x
}

func authalicQ(sinPhi, e, e2 float64) float64 {
	if e < 1e-12 {
		return 2 * sinPhi
	}
	es := e * sinPhi
	return (1 - e2) * (sinPhi/(1-e2*sinPhi*sinPhi) - (1/(2*e))*math.Log((1-es)/(1+es)))
}

func authalicToGeodetic(beta float64, s Spheroid) float64 {
	e2, e4, e6 := s.E2(), s.E4(), s.E6()
	return beta +
		(e2/3+31*e4/180+517*e6/5040)*math.Sin(2*beta) +
		(23*e4/360+251*e6/3780)*math.Sin(4*beta) +
		(761*e6/45360)*math.Sin(6*beta)
}

func build(txt string) builder {
	b := &strings.Builder{}

	if len(txt) > 0 {
		b.WriteString(txt)
	}

	return builder{
		s: b,
	}
}

type builder struct {
	s   *strings.Builder
	err error
}

func (b builder) addKey(key string) builder {
	b.s.WriteString(" " + key)

	return b
}

func (b builder) addErr(err error) builder {
	b.err = errors.Join(b.err, err)

	return b
}

func (b builder) addAll(kv ...any) builder {
	for i := 0; i+1 < len(kv); i += 2 {
		k := kv[i]
		kk, ok := k.(string)
		if !ok {
			return b.addErr(fmt.Errorf("invalid key: %v", k))
		}

		b = b.add(kk, kv[i+1])
	}

	return b
}

func (b builder) add(key string, val any) builder {
	switch v := val.(type) {
	case nil:
		return b.addKey(key)
	case string:
		if v == "" {
			return b
		}
	case fmt.Stringer:
		val = v.String()
		if val == "" {
			return b
		}
	case float64:
		if v == 0 {
			return b
		}

		val = formatFloat(v)
	}

	if key == "" {
		if val == nil {
			return b
		}

		r := reflect.ValueOf(val)
		if r.Kind() == reflect.Slice && r.Type().Elem().Implements(reflect.TypeFor[fmt.Stringer]()) {
			for l := range r.Len() {
				fmt.Fprintf(b.s, " %v", r.Index(l))
			}

			return b
		}

		b.s.WriteString(" " + r.String())

		return b
	}

	fmt.Fprintf(b.s, " %s=%v", key, val)

	return b
}

func (b builder) String() string {
	if b.err != nil {
		return b.err.Error()
	}

	return strings.TrimSpace(b.s.String())
}
