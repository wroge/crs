package crs

import "math"

// timeSpecificEpochTol is the max |coordEpoch − TransformationEpoch| (years)
// for a time-specific Helmert to participate in path selection.
const timeSpecificEpochTol = 1e-6

type TimeSpecificPositionVector struct {
	Tx, Ty, Tz, Rx, Ry, Rz, Ds float64
	TransformationEpoch        float64
}

// timeSpecificEpoch returns the transformation epoch when op is (or wraps)
// TimeSpecificPositionVector.
func timeSpecificEpoch(op Operation) (epoch float64, ok bool) {
	switch o := op.(type) {
	case TimeSpecificPositionVector:
		return o.TransformationEpoch, true
	case Inverse:
		return timeSpecificEpoch(o.Operation)
	default:
		return 0, false
	}
}

// timeSpecificAllowed reports whether a time-specific op may be used for the
// given coordinate epoch. Non-time-specific ops always return true.
// Without a coordinate epoch (Transform), time-specific edges are excluded.
func timeSpecificAllowed(op Operation, hasCoordEpoch bool, coordEpoch float64) bool {
	te, isTS := timeSpecificEpoch(op)
	if !isTS {
		return true
	}
	if !hasCoordEpoch {
		return false
	}
	return math.Abs(coordEpoch-te) <= timeSpecificEpochTol
}

func (t TimeSpecificPositionVector) String() string {
	return build("time_specific_position_vector").addAll(
		"tx", t.Tx,
		"ty", t.Ty,
		"tz", t.Tz,
		"rx", t.Rx,
		"ry", t.Ry,
		"rz", t.Rz,
		"ds", t.Ds,
		"epoch", t.TransformationEpoch,
	).String()
}

func (t TimeSpecificPositionVector) asPositionVector() PositionVector {
	return PositionVector{
		Tx: t.Tx, Ty: t.Ty, Tz: t.Tz,
		Rx: t.Rx, Ry: t.Ry, Rz: t.Rz,
		Ds: t.Ds,
	}
}

func (t TimeSpecificPositionVector) ToTarget(source, target Spheroid, lon, lat, h float64) (float64, float64, float64, error) {
	return t.asPositionVector().ToTarget(source, target, lon, lat, h)
}

func (t TimeSpecificPositionVector) FromTarget(source, target Spheroid, lon0, lat0, h0 float64) (float64, float64, float64, error) {
	return t.asPositionVector().FromTarget(source, target, lon0, lat0, h0)
}

type TimeDependentPositionVector struct {
	Tx, Ty, Tz, Rx, Ry, Rz, Ds                             float64
	TxRate, TyRate, TzRate, RxRate, RyRate, RzRate, DsRate float64
	ReferenceEpoch                                         float64
}

func (t TimeDependentPositionVector) String() string {
	return build("time_dependent_position_vector").addAll(
		"tx", t.Tx,
		"ty", t.Ty,
		"tz", t.Tz,
		"rx", t.Rx,
		"ry", t.Ry,
		"rz", t.Rz,
		"ds", t.Ds,
		"dtx", t.TxRate,
		"dty", t.TyRate,
		"dtz", t.TzRate,
		"drx", t.RxRate,
		"dry", t.RyRate,
		"drz", t.RzRate,
		"dds", t.DsRate,
		"epoch", t.ReferenceEpoch,
	).String()
}

func (t TimeDependentPositionVector) at(epoch float64) PositionVector {
	dt := epoch - t.ReferenceEpoch

	return PositionVector{
		Tx: t.Tx + t.TxRate*dt,
		Ty: t.Ty + t.TyRate*dt,
		Tz: t.Tz + t.TzRate*dt,
		Rx: t.Rx + t.RxRate*dt,
		Ry: t.Ry + t.RyRate*dt,
		Rz: t.Rz + t.RzRate*dt,
		Ds: t.Ds + t.DsRate*dt,
	}
}

// ToTarget evaluates Helmert parameters at ReferenceEpoch (dt=0). Use
// Datum.AtEpoch / TransformAt so rates apply at a coordinate epoch.
func (t TimeDependentPositionVector) ToTarget(source, target Spheroid, lon, lat, h float64) (float64, float64, float64, error) {
	return t.at(t.ReferenceEpoch).ToTarget(source, target, lon, lat, h)
}

// FromTarget evaluates Helmert parameters at ReferenceEpoch (dt=0). Use
// Datum.AtEpoch / TransformAt so rates apply at a coordinate epoch.
func (t TimeDependentPositionVector) FromTarget(source, target Spheroid, lon0, lat0, h0 float64) (float64, float64, float64, error) {
	return t.at(t.ReferenceEpoch).FromTarget(source, target, lon0, lat0, h0)
}
