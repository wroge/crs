package crs

import (
	"flag"
	"fmt"
	"math"
	"os"
	"strings"
	"testing"
)

var cache = flag.Bool("cache", false, "cache grids")

func TestMain(m *testing.M) {
	flag.Parse()

	if *cache {
		SetGridCDN("https://cdn.proj.org/", "./.cache")
	} else {
		SetGridCDN("https://cdn.proj.org/")
	}

	os.Exit(m.Run())
}

const transformTol = 0.01
const transformTolDeg = 2e-5

type transformPt struct{ A, B, C float64 }

type transformCase struct {
	Name     string
	FromEPSG int
	ToEPSG   int
	CoversCS string
	CoversOp string
	Epoch    float64 // 0 → Transform; otherwise TransformAt(..., Epoch)
	In       [4]transformPt
	Want     [4]transformPt
}

var transformCases = []transformCase{
	{
		Name:     "conv_geographic",
		FromEPSG: 4326,
		ToEPSG:   4258,
		CoversCS: "geographic",
		CoversOp: "",
		In: [4]transformPt{
			{A: 15.48163044, B: 44.32638214, C: 0.0},
			{A: 3.65110185, B: 50.44719402, C: 0.0},
			{A: 18.63227444, B: 64.45183357, C: 0.0},
			{A: 23.68750185, B: 46.23884503, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 15.48163044, B: 44.32638214, C: 0.0},
			{A: 3.65110185, B: 50.44719402, C: 0.0},
			{A: 18.63227444, B: 64.45183357, C: 0.0},
			{A: 23.68750185, B: 46.23884503, C: 0.0},
		},
	},
	{
		Name:     "conv_geocentric",
		FromEPSG: 4326,
		ToEPSG:   4978,
		CoversCS: "geocentric",
		CoversOp: "",
		In: [4]transformPt{
			{A: -16.86488695, B: -50.7819003, C: 0.0},
			{A: -60.77419744, B: 0.57837112, C: 0.0},
			{A: -102.26823055, B: -32.52553373, C: 0.0},
			{A: 32.37503856, B: 4.85367991, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 3867066.28567992, B: -1172316.40364588, C: -4918239.7200331},
			{A: 3113985.14541967, B: -5565933.0480599, C: 63951.90291349},
			{A: -1143815.46459941, B: -5260022.30647536, C: -3409710.96305441},
			{A: 5367539.52181719, B: 3403064.87021672, C: 536063.3112849},
		},
	},
	{
		Name:     "conv_transverse_mercator",
		FromEPSG: 4326,
		ToEPSG:   32632,
		CoversCS: "transverse_mercator",
		CoversOp: "",
		In: [4]transformPt{
			{A: 7.99358624, B: 46.49899047, C: 0.0},
			{A: 10.11394964, B: 17.12753749, C: 0.0},
			{A: 10.10094931, B: 51.98622551, C: 0.0},
			{A: 8.42490186, B: 24.63616679, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 422774.7445097, B: 5149983.18014845, C: 0.0},
			{A: 618498.12624472, B: 1894003.22961724, C: 0.0},
			{A: 575602.8462278, B: 5760078.48038194, C: 0.0},
			{A: 441796.64488126, B: 2724783.76313037, C: 0.0},
		},
	},
	{
		Name:     "conv_web_mercator",
		FromEPSG: 4326,
		ToEPSG:   3857,
		CoversCS: "web_mercator",
		CoversOp: "",
		In: [4]transformPt{
			{A: 12.74327843, B: 49.01956727, C: 0.0},
			{A: 7.55647506, B: 47.58029826, C: 0.0},
			{A: 12.0849662, B: 50.62235619, C: 0.0},
			{A: 11.84276964, B: 51.37839072, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 1418575.26622526, B: 6278182.202867, C: 0.0},
			{A: 841182.95590254, B: 6037313.29273976, C: 0.0},
			{A: 1345292.28342472, B: 6554762.81796277, C: 0.0},
			{A: 1318331.08586746, B: 6688501.57924964, C: 0.0},
		},
	},
	{
		Name:     "conv_lambert_conformal_conic",
		FromEPSG: 4807,
		ToEPSG:   27561,
		CoversCS: "lambert_conformal_conic",
		CoversOp: "",
		In: [4]transformPt{
			{A: 1.9647528, B: 50.49160838, C: 0.0},
			{A: 0.7252802, B: 49.73367314, C: 0.0},
			{A: 4.26912066, B: 49.85333555, C: 0.0},
			{A: 4.52301624, B: 49.77923386, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 739417.87526692, B: 312105.82514745, C: 0.0},
			{A: 652282.00028793, B: 226238.10594434, C: 0.0},
			{A: 906827.12674424, B: 247988.96797489, C: 0.0},
			{A: 925548.01703143, B: 240827.25956769, C: 0.0},
		},
	},
	{
		Name:     "conv_lambert_conformal_conic_1sp_variant_b",
		FromEPSG: 4326,
		ToEPSG:   9549,
		CoversCS: "lambert_conformal_conic_1sp_variant_b",
		CoversOp: "",
		In: [4]transformPt{
			{A: 6.66146022, B: 45.10204452, C: 0.0},
			{A: 5.73766686, B: 45.25110543, C: 0.0},
			{A: 5.45063685, B: 45.21646802, C: 0.0},
			{A: 5.49174077, B: 45.24411985, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 137783.23857711, B: 40976.772433, C: 0.0},
			{A: 65289.32734678, B: 58090.60925413, C: 0.0},
			{A: 42692.70280965, B: 54577.51904063, C: 0.0},
			{A: 45971.00551625, B: 57597.54579158, C: 0.0},
		},
	},
	{
		Name:     "conv_lambert_conformal_conic_2sp",
		FromEPSG: 4326,
		ToEPSG:   2154,
		CoversCS: "lambert_conformal_conic_2sp",
		CoversOp: "",
		In: [4]transformPt{
			{A: 1.90775189, B: 45.51074179, C: 0.0},
			{A: -1.31652234, B: 44.54058091, C: 0.0},
			{A: -2.56982133, B: 49.08234455, C: 0.0},
			{A: 2.05774172, B: 47.03663226, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 614718.86071981, B: 6490730.9830579, C: 0.0},
			{A: 357241.97718452, B: 6391760.70695503, C: 0.0},
			{A: 293429.50692894, B: 6901289.77477103, C: 0.0},
			{A: 628451.34549122, B: 6660026.50210425, C: 0.0},
		},
	},
	{
		Name:     "conv_lambert_conformal_conic_2sp_belgium",
		FromEPSG: 4313,
		ToEPSG:   31300,
		CoversCS: "lambert_conformal_conic_2sp_belgium",
		CoversOp: "",
		In: [4]transformPt{
			{A: 3.68046444, B: 50.78132692, C: 0.0},
			{A: 3.66236184, B: 50.35962326, C: 0.0},
			{A: 5.59548464, B: 50.67383971, C: 0.0},
			{A: 4.5832624, B: 50.72764479, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 101547.37411825, B: 163590.78584881, C: 0.0},
			{A: 99825.25528493, B: 116695.86330683, C: 0.0},
			{A: 236800.30170348, B: 152127.48061708, C: 0.0},
			{A: 165235.26922413, B: 157417.09973715, C: 0.0},
		},
	},
	{
		Name:     "conv_lambert_conformal_conic_2sp_michigan",
		FromEPSG: 4326,
		ToEPSG:   6201,
		CoversCS: "lambert_conformal_conic_2sp_michigan",
		CoversOp: "",
		In: [4]transformPt{
			{A: -83.67964358, B: 45.21107189, C: 0.0},
			{A: -85.44371584, B: 44.26483151, C: 0.0},
			{A: -85.19538794, B: 44.56456639, C: 0.0},
			{A: -85.49563531, B: 45.42338116, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 660949.48708435, B: 210716.10609631, C: 0.0},
			{A: 520939.95443241, B: 105966.53490951, C: 0.0},
			{A: 541121.12498766, B: 139030.0210911, C: 0.0},
			{A: 518638.48112929, B: 234755.2120703, C: 0.0},
		},
	},
	{
		Name:     "conv_lambert_conic_near_conformal",
		FromEPSG: 4227,
		ToEPSG:   22700,
		CoversCS: "lambert_conic_near_conformal",
		CoversOp: "",
		In: [4]transformPt{
			{A: 38.9032823, B: 34.16642019, C: 0.0},
			{A: 38.3730528, B: 34.31213742, C: 0.0},
			{A: 38.99491422, B: 34.42593333, C: 0.0},
			{A: 37.4357124, B: 34.04392951, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 443169.96688565, B: 247481.59290632, C: 0.0},
			{A: 394135.45325473, B: 263013.67101362, C: 0.0},
			{A: 451145.49773539, B: 276387.62738511, C: 0.0},
			{A: 307912.19328097, B: 232799.38564307, C: 0.0},
		},
	},
	{
		Name:     "conv_lambert_azimuthal_equal_area",
		FromEPSG: 4326,
		ToEPSG:   3035,
		CoversCS: "lambert_azimuthal_equal_area",
		CoversOp: "",
		In: [4]transformPt{
			{A: 7.585767, B: 46.10519175, C: 0.0},
			{A: 8.70593568, B: 69.01765399, C: 0.0},
			{A: -0.22852323, B: 44.53865435, C: 0.0},
			{A: 28.62919936, B: 55.00868962, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 4134151.43163056, B: 2557719.0695705, C: 0.0},
			{A: 4268676.27646318, B: 5099267.82196468, C: 0.0},
			{A: 3509347.39858626, B: 2436493.51468078, C: 0.0},
			{A: 5497973.50677255, B: 3698108.93335958, C: 0.0},
		},
	},
	{
		Name:     "conv_lambert_azimuthal_equal_area_spherical",
		FromEPSG: 4326,
		ToEPSG:   3408,
		CoversCS: "lambert_azimuthal_equal_area_spherical",
		CoversOp: "",
		In: [4]transformPt{
			{A: -94.90908705, B: 44.56539651, C: 0.0},
			{A: -94.68421044, B: 51.5293525, C: 0.0},
			{A: -86.49504763, B: 49.0659196, C: 0.0},
			{A: -95.23766753, B: 48.57943144, C: 0.0},
		},
		Want: [4]transformPt{
			{A: -4902889.04228191, B: 421109.18207208, C: 0.0},
			{A: -4183965.35670889, B: 342823.78453394, C: 0.0},
			{A: -4447301.23965291, B: -272394.43984522, C: 0.0},
			{A: -4487461.60882722, B: 411365.76876755, C: 0.0},
		},
	},
	{
		Name:     "conv_lambert_cylindrical_equal_area",
		FromEPSG: 4326,
		ToEPSG:   6933,
		CoversCS: "lambert_cylindrical_equal_area",
		CoversOp: "",
		In: [4]transformPt{
			{A: 11.90691313, B: 0.34937214, C: 0.0},
			{A: 11.30588106, B: 4.32935643, C: 0.0},
			{A: -11.72445547, B: 2.64866183, C: 0.0},
			{A: 4.36104886, B: 0.44364396, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 1148853.75677538, B: 44570.28420093, C: 0.0},
			{A: 1090862.40870157, B: 551798.84768238, C: 0.0},
			{A: -1131249.09659296, B: 337781.38344636, C: 0.0},
			{A: 420781.38216941, B: 56596.5691843, C: 0.0},
		},
	},
	{
		Name:     "conv_lambert_cylindrical_equal_area_spherical",
		FromEPSG: 4326,
		ToEPSG:   3410,
		CoversCS: "lambert_cylindrical_equal_area_spherical",
		CoversOp: "",
		In: [4]transformPt{
			{A: -5.59619544, B: 1.69154158, C: 0.0},
			{A: -9.32274783, B: -0.78281699, C: 0.0},
			{A: -1.11063105, B: 5.44579113, C: 0.0},
			{A: 9.02047057, B: -2.83933139, C: 0.0},
		},
		Want: [4]transformPt{
			{A: -538919.650686, B: 217164.81912794, C: 0.0},
			{A: -897790.66131902, B: -100511.69921992, C: 0.0},
			{A: -106954.96659794, B: 698195.05482214, C: 0.0},
			{A: 868681.03504275, B: -364425.02982483, C: 0.0},
		},
	},
	{
		Name:     "conv_albers_equal_area",
		FromEPSG: 4269,
		ToEPSG:   5070,
		CoversCS: "albers_equal_area",
		CoversOp: "",
		In: [4]transformPt{
			{A: -95.82964547, B: 32.08056247, C: 0.0},
			{A: -81.5202604, B: 42.44610921, C: 0.0},
			{A: -102.84960928, B: 38.97674133, C: 0.0},
			{A: -92.0656825, B: 31.69383792, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 16004.45843208, B: 1000772.68561012, C: 0.0},
			{A: 1179018.22056504, B: 2250796.36588922, C: 0.0},
			{A: -587388.33971401, B: 1793480.6823246, C: 0.0},
			{A: 371298.62224231, B: 965365.78576565, C: 0.0},
		},
	},
	{
		Name:     "conv_krovak",
		FromEPSG: 4156,
		ToEPSG:   5514,
		CoversCS: "krovak",
		CoversOp: "",
		In: [4]transformPt{
			{A: 18.97409285, B: 49.4736793, C: 0.0},
			{A: 19.07533154, B: 49.45564664, C: 0.0},
			{A: 14.18759265, B: 49.0436638, C: 0.0},
			{A: 14.3063529, B: 50.25233904, C: 0.0},
		},
		Want: [4]transformPt{
			{A: -424218.03464014, B: -1145772.79202197, C: 0.0},
			{A: -417053.41936915, B: -1148328.05736525, C: 0.0},
			{A: -775634.41618565, B: -1155722.51542267, C: 0.0},
			{A: -748515.76605619, B: -1023801.18510333, C: 0.0},
		},
	},
	{
		Name:     "conv_krovak_modified",
		FromEPSG: 4326,
		ToEPSG:   5515,
		CoversCS: "krovak_modified",
		CoversOp: "",
		In: [4]transformPt{
			{A: 17.01336827, B: 50.31351831, C: 0.0},
			{A: 14.69312238, B: 49.16219265, C: 0.0},
			{A: 17.01047499, B: 50.48506077, C: 0.0},
			{A: 13.79192432, B: 49.79915381, C: 0.0},
		},
		Want: [4]transformPt{
			{A: -5556076.62449541, B: -6040110.75594175, C: 0.0},
			{A: -5737189.49154148, B: -6147621.58248851, C: 0.0},
			{A: -5554332.76223784, B: -6021105.54506083, C: 0.0},
			{A: -5792035.41710906, B: -6068414.11965381, C: 0.0},
		},
	},
	{
		Name:     "conv_hotine_oblique_mercator_a",
		FromEPSG: 4326,
		ToEPSG:   3375,
		CoversCS: "hotine_oblique_mercator_a",
		CoversOp: "",
		In: [4]transformPt{
			{A: 101.56028556, B: 3.31488721, C: 0.0},
			{A: 102.1039057, B: 5.39156372, C: 0.0},
			{A: 100.52138698, B: 3.7153589, C: 0.0},
			{A: 104.23729899, B: 5.07071126, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 396048.98430906, B: 366873.21547655, C: 0.0},
			{A: 456912.76753137, B: 596345.21518618, C: 0.0},
			{A: 280751.80746508, B: 411539.12919044, C: 0.0},
			{A: 693448.66268388, B: 560746.1338808, C: 0.0},
		},
	},
	{
		Name:     "conv_swiss_oblique_mercator",
		FromEPSG: 4326,
		ToEPSG:   2056,
		CoversCS: "swiss_oblique_mercator",
		CoversOp: "",
		In: [4]transformPt{
			{A: 8.05601323, B: 46.83109101, C: 0.0},
			{A: 7.19034145, B: 46.4796368, C: 0.0},
			{A: 7.78297358, B: 46.91597046, C: 0.0},
			{A: 7.48821322, B: 46.47426086, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 2647103.10993301, B: 1186845.44389839, C: 0.0},
			{A: 2580932.03171158, B: 1147620.99414099, C: 0.0},
			{A: 2626230.32030269, B: 1196153.83423511, C: 0.0},
			{A: 2603807.26797883, B: 1146994.11162607, C: 0.0},
		},
	},
	{
		Name:     "conv_laborde_oblique_mercator",
		FromEPSG: 4326,
		ToEPSG:   8441, // tananarive (Greenwich); 29701 is tananariveparis
		CoversCS: "laborde_oblique_mercator",
		CoversOp: "",
		In: [4]transformPt{
			{A: 44.97035738, B: -17.6834006, C: 0.0},
			{A: 45.66975422, B: -15.42028489, C: 0.0},
			{A: 48.46246555, B: -22.30542686, C: 0.0},
			{A: 45.70988452, B: -17.37093333, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 244451.21284211, B: 934134.73648334, C: 0.0},
			{A: 317729.07277262, B: 1184977.74002525, C: 0.0},
			{A: 608656.79967132, B: 421810.84404126, C: 0.0},
			{A: 322760.88157918, B: 969164.07375348, C: 0.0},
		},
	},
	{
		Name:     "conv_cassini_soldner",
		FromEPSG: 4326,
		ToEPSG:   24500,
		CoversCS: "cassini_soldner",
		CoversOp: "",
		In: [4]transformPt{
			{A: 103.7477002, B: 1.22499162, C: 0.0},
			{A: 103.9554281, B: 1.31449279, C: 0.0},
			{A: 103.82212926, B: 1.35806236, C: 0.0},
			{A: 103.91855914, B: 1.23684362, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 18473.57705848, B: 23077.56717961, C: 0.0},
			{A: 41592.0848522, B: 32974.14064568, C: 0.0},
			{A: 26757.27437176, B: 37791.62588036, C: 0.0},
			{A: 37489.14470825, B: 24387.95918432, C: 0.0},
		},
	},
	{
		Name:     "conv_hyperbolic_cassini_soldner",
		FromEPSG: 4326,
		ToEPSG:   3139,
		CoversCS: "hyperbolic_cassini_soldner",
		CoversOp: "",
		In: [4]transformPt{
			{A: 179.20652679, B: -16.62482157, C: 0.0},
			{A: 179.36331774, B: -16.61187112, C: 0.0},
			{A: 179.50995641, B: -16.53758876, C: 0.0},
			{A: 179.6323993, B: -16.74456957, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 238589.584338, B: 292896.50001061, C: 0.0},
			{A: 255318.95600776, B: 294333.70305188, C: 0.0},
			{A: 270973.1678762, B: 302546.10995721, C: 0.0},
			{A: 284009.1162004, B: 279624.95067609, C: 0.0},
		},
	},
	{
		Name:     "conv_polar_stereographic_a",
		FromEPSG: 4326,
		ToEPSG:   5041,
		CoversCS: "polar_stereographic_a",
		CoversOp: "",
		In: [4]transformPt{
			{A: 17.66291077, B: 74.03581563, C: 0.0},
			{A: 28.68014087, B: 73.491938, C: 0.0},
			{A: 12.5650138, B: 74.69168129, C: 0.0},
			{A: 18.12515936, B: 73.67127087, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 2541142.18976398, B: 300580.11361035, C: 0.0},
			{A: 2885480.39502905, B: 381305.27382312, C: 0.0},
			{A: 2371867.15356652, B: 331577.15868671, C: 0.0},
			{A: 2567669.12327082, B: 265793.25922307, C: 0.0},
		},
	},
	{
		Name:     "conv_polar_stereographic_b",
		FromEPSG: 4326,
		ToEPSG:   3031,
		CoversCS: "polar_stereographic_b",
		CoversOp: "",
		In: [4]transformPt{
			{A: -39.00696788, B: -72.4604064, C: 0.0},
			{A: -32.04729318, B: -72.83190537, C: 0.0},
			{A: -28.18828875, B: -77.69647002, C: 0.0},
			{A: -12.02583114, B: -72.98383449, C: 0.0},
		},
		Want: [4]transformPt{
			{A: -1208565.53356152, B: 1492083.08512075, C: 0.0},
			{A: -996966.56776788, B: 1592553.41950768, C: 0.0},
			{A: -633808.76207812, B: 1182628.71922717, C: 0.0},
			{A: -387954.5048643, B: 1821144.84065904, C: 0.0},
		},
	},
	{
		Name:     "conv_polar_stereographic_c",
		FromEPSG: 4636,
		ToEPSG:   2985,
		CoversCS: "polar_stereographic_c",
		CoversOp: "",
		In: [4]transformPt{
			{A: 140.58139775, B: -66.37325394, C: 0.0},
			{A: 140.50921744, B: -66.510064, C: 0.0},
			{A: 140.29138468, B: -66.50152549, C: 0.0},
			{A: 140.24062418, B: -66.52944563, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 326071.41218788, B: 269842.53722165, C: 0.0},
			{A: 322698.88070274, B: 254586.39156227, C: 0.0},
			{A: 312993.73179225, B: 255608.1240706, C: 0.0},
			{A: 310717.08185416, B: 252499.58889688, C: 0.0},
		},
	},
	{
		Name:     "conv_oblique_stereographic",
		FromEPSG: 4289,
		ToEPSG:   28992,
		CoversCS: "oblique_stereographic",
		CoversOp: "",
		In: [4]transformPt{
			{A: 4.91808316, B: 53.08399665, C: 0.0},
			{A: 4.64366978, B: 52.72780497, C: 0.0},
			{A: 5.10148018, B: 52.08872325, C: 0.0},
			{A: 6.31305015, B: 53.10189816, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 123541.39531833, B: 566332.08004981, C: 0.0},
			{A: 104748.5781735, B: 526856.15419413, C: 0.0},
			{A: 135390.2968795, B: 455536.4759189, C: 0.0},
			{A: 216972.88195587, B: 568619.34919786, C: 0.0},
		},
	},
	{
		Name:     "conv_mercator_a",
		FromEPSG: 4326,
		ToEPSG:   3395,
		CoversCS: "mercator_a",
		CoversOp: "",
		In: [4]transformPt{
			{A: 10.33460994, B: 6.31044965, C: 0.0},
			{A: 7.92878095, B: 3.78024695, C: 0.0},
			{A: 12.81225619, B: 5.47508174, C: 0.0},
			{A: 10.25317121, B: 6.48785336, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 1150443.51611618, B: 699207.31771383, C: 0.0},
			{A: 882627.85814846, B: 418305.72391995, C: 0.0},
			{A: 1426253.83496317, B: 606339.015668, C: 0.0},
			{A: 1141377.79791816, B: 718948.23848198, C: 0.0},
		},
	},
	{
		Name:     "conv_mercator_b",
		FromEPSG: 4326,
		ToEPSG:   3388,
		CoversCS: "mercator_b",
		CoversOp: "",
		In: [4]transformPt{
			{A: 48.58540816, B: 42.64587308, C: 0.0},
			{A: 50.4519374, B: 44.19589922, C: 0.0},
			{A: 49.00532826, B: 44.81961583, C: 0.0},
			{A: 48.68150682, B: 40.34658167, C: 0.0},
		},
		Want: [4]transformPt{
			{A: -199954.71684862, B: 3892011.27289067, C: 0.0},
			{A: -45310.73057851, B: 4068204.3230576, C: 0.0},
			{A: -165161.7254104, B: 4140420.05650878, C: 0.0},
			{A: -191994.87873024, B: 3638588.44852621, C: 0.0},
		},
	},
	{
		Name:     "conv_azimuthal_equidistant",
		FromEPSG: 4326,
		ToEPSG:   27701,
		CoversCS: "azimuthal_equidistant",
		CoversOp: "",
		In: [4]transformPt{
			{A: 30.13819767, B: 6.14799563, C: 0.0},
			{A: 6.17775696, B: -21.40394554, C: 0.0},
			{A: 49.79845167, B: -15.13627173, C: 0.0},
			{A: 30.1038414, B: 3.37799425, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 6577730.5905754, B: 5738511.56313955, C: 0.0},
			{A: 3959330.83521181, B: 2655513.1382068, C: 0.0},
			{A: 8749501.83464612, B: 3366643.72600085, C: 0.0},
			{A: 6578798.62537088, B: 5428978.81984501, C: 0.0},
		},
	},
	{
		Name:     "conv_modified_azimuthal_equidistant",
		FromEPSG: 4675,
		ToEPSG:   3295,
		CoversCS: "modified_azimuthal_equidistant",
		CoversOp: "",
		In: [4]transformPt{
			{A: 138.12676094, B: 9.54353457, C: 0.0},
			{A: 138.13732384, B: 9.57934004, C: 0.0},
			{A: 138.10483444, B: 9.55705156, C: 0.0},
			{A: 138.10834597, B: 9.52437016, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 35390.62163991, B: 59649.26829084, C: 0.0},
			{A: 36550.68488897, B: 63609.16489965, C: 0.0},
			{A: 32983.58253582, B: 61144.59117028, C: 0.0},
			{A: 33368.46720897, B: 57530.02101084, C: 0.0},
		},
	},
	{
		Name:     "conv_american_polyconic",
		FromEPSG: 4326,
		ToEPSG:   5472,
		CoversCS: "american_polyconic",
		CoversOp: "",
		In: [4]transformPt{
			{A: -79.51236732, B: 8.11139557, C: 0.0},
			{A: -80.76021804, B: 8.7973303, C: 0.0},
			{A: -81.61537367, B: 8.35167742, C: 0.0},
			{A: -78.36542491, B: 9.16807441, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 1078350.02212542, B: 984376.79466789, C: 0.0},
			{A: 940772.55283694, B: 1059943.4535778, C: 0.0},
			{A: 846609.62550406, B: 1010702.16081649, C: 0.0},
			{A: 1203950.20429134, B: 1101998.01258345, C: 0.0},
		},
	},
	{
		Name:     "conv_bonne_south_orientated",
		FromEPSG: 4666,
		ToEPSG:   5017,
		CoversCS: "bonne_south_orientated",
		CoversOp: "",
		In: [4]transformPt{
			{A: -8.73786682, B: 38.65832038, C: 0.0},
			{A: -8.34976476, B: 40.90936882, C: 0.0},
			{A: -7.10489264, B: 40.74059878, C: 0.0},
			{A: -8.13881623, B: 38.4851166, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 52736.99599153, B: 111755.12461984, C: 0.0},
			{A: 18352.43643939, B: -137997.1733018, C: 0.0},
			{A: -86733.52072543, B: -119731.05169191, C: 0.0},
			{A: 602.84022718, B: 131157.90568704, C: 0.0},
		},
	},
	{
		Name:     "conv_equal_earth",
		FromEPSG: 4326,
		ToEPSG:   8857,
		CoversCS: "equal_earth",
		CoversOp: "",
		In: [4]transformPt{
			{A: 8.00987891, B: 2.4424791, C: 0.0},
			{A: 2.68026638, B: 5.84679676, C: 0.0},
			{A: 3.69543163, B: -5.90612271, C: 0.0},
			{A: 7.61049924, B: -2.40745497, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 767014.51681117, B: 313780.48854927, C: 0.0},
			{A: 256136.24310527, B: 750587.95930756, C: 0.0},
			{A: 353131.43297284, B: -758190.5313808, C: 0.0},
			{A: 728779.42133867, B: -309282.34146451, C: 0.0},
		},
	},
	{
		Name:     "conv_equidistant_cylindrical",
		FromEPSG: 4326,
		ToEPSG:   4087,
		CoversCS: "equidistant_cylindrical",
		CoversOp: "",
		In: [4]transformPt{
			{A: 3.92132916, B: 5.26716005, C: 0.0},
			{A: -8.77701325, B: -4.61485595, C: 0.0},
			{A: -9.43113653, B: 0.63868369, C: 0.0},
			{A: -5.4636429, B: 1.25795792, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 436520.36523339, B: 582428.85546717, C: 0.0},
			{A: -977052.64618334, B: -510295.42184326, C: 0.0},
			{A: -1049869.31667514, B: 70622.0159413, C: 0.0},
			{A: -608209.94599903, B: 139098.01091278, C: 0.0},
		},
	},
	{
		Name:     "conv_colombia_urban",
		FromEPSG: 4326,
		ToEPSG:   6247,
		CoversCS: "colombia_urban",
		CoversOp: "",
		In: [4]transformPt{
			{A: -74.08713548, B: 4.57886335, C: 0.0},
			{A: -74.1016426, B: 4.59335614, C: 0.0},
			{A: -74.12699546, B: 4.74728076, C: 0.0},
			{A: -74.06477795, B: 4.55215163, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 98935.16587869, B: 98079.10764833, C: 0.0},
			{A: 97324.616844, B: 99682.26887471, C: 0.0},
			{A: 94509.75441882, B: 116710.23163772, C: 0.0},
			{A: 101417.43154001, B: 95124.34167636, C: 0.0},
		},
	},
	{
		Name:     "conv_guam_projection",
		FromEPSG: 4675,
		ToEPSG:   3993,
		CoversCS: "guam_projection",
		CoversOp: "",
		In: [4]transformPt{
			{A: 144.76353637, B: 13.40980244, C: 0.0},
			{A: 144.70053185, B: 13.49880146, C: 0.0},
			{A: 144.79556701, B: 13.40715195, C: 0.0},
			{A: 144.81118464, B: 13.45930248, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 51601.36677254, B: 43067.67132008, C: 0.0},
			{A: 44779.56609807, B: 52913.92060815, C: 0.0},
			{A: 55070.51350428, B: 42774.88603382, C: 0.0},
			{A: 56760.54510658, B: 48544.56246734, C: 0.0},
		},
	},
	{
		Name:     "conv_local_orthographic",
		FromEPSG: 4326,
		ToEPSG:   10622,
		CoversCS: "local_orthographic",
		CoversOp: "",
		In: [4]transformPt{
			{A: -122.38803715, B: 37.60234811, C: 0.0},
			{A: -122.40283976, B: 37.63379183, C: 0.0},
			{A: -122.368035, B: 37.62164125, C: 0.0},
			{A: -122.37094701, B: 37.62297034, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 1838.85835951, B: -2370.62941328, C: 0.0},
			{A: -944.45911115, B: 107.36484183, C: 0.0},
			{A: 2402.26690605, B: 347.22983093, C: 0.0},
			{A: 2106.06352548, B: 357.78253858, C: 0.0},
		},
	},
	{
		Name:     "conv_tunisia_mining_grid",
		FromEPSG: 4816,
		ToEPSG:   22300,
		CoversCS: "tunisia_mining_grid",
		CoversOp: "",
		In: [4]transformPt{
			{A: 9.07771254, B: 34.19116828, C: 0.0},
			{A: 9.26991002, B: 35.34847223, C: 0.0},
			{A: 9.85534677, B: 35.29105387, C: 0.0},
			{A: 9.97870956, B: 34.31511481, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 454808.967554, B: 499169.94274402, C: 0.0},
			{A: 472334.83973156, B: 627566.67588771, C: 0.0},
			{A: 525718.9409443, B: 621196.41349103, C: 0.0},
			{A: 536967.99925517, B: 512921.15243908, C: 0.0},
		},
	},
	{
		Name:     "conv_new_zealand_map_grid",
		FromEPSG: 4326,
		ToEPSG:   27200,
		CoversCS: "new_zealand_map_grid",
		CoversOp: "",
		In: [4]transformPt{
			{A: 170.65754076, B: -43.86802423, C: 0.0},
			{A: 174.56053508, B: -37.66573874, C: 0.0},
			{A: 171.81131208, B: -39.75718893, C: 0.0},
			{A: 169.95889436, B: -37.3026045, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 2321745.06165214, B: 5701743.34937189, C: 0.0},
			{A: 2647906.4927081, B: 6391850.74313471, C: 0.0},
			{A: 2408109.56053257, B: 6160318.18041594, C: 0.0},
			{A: 2240274.34412839, B: 6429588.385578, C: 0.0},
		},
	},
	{
		Name:     "op_geocentric_translations",
		FromEPSG: 4326,
		ToEPSG:   2043,
		CoversCS: "",
		CoversOp: "geocentric_translations",
		In: [4]transformPt{
			{A: -6.73402748, B: 9.35791735, C: 0.0},
			{A: -6.81833149, B: 8.99108071, C: 0.0},
			{A: -8.04918456, B: 8.43050451, C: 0.0},
			{A: -7.56779756, B: 9.18225748, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 748843.40468147, B: 1034759.19201818, C: 0.0},
			{A: 739826.20694532, B: 994115.4753244, C: 0.0},
			{A: 604638.74894814, B: 931553.87111874, C: 0.0},
			{A: 657323.52055924, B: 1014851.81157883, C: 0.0},
		},
	},
	{
		Name:     "op_position_vector",
		FromEPSG: 4326,
		ToEPSG:   2137,
		CoversCS: "",
		CoversOp: "position_vector",
		In: [4]transformPt{
			{A: 0.22309898, B: 4.74792302, C: 0.0},
			{A: 0.25318808, B: 3.07798876, C: 0.0},
			{A: 0.17058152, B: 2.63423537, C: 0.0},
			{A: 0.47023741, B: 4.73262673, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 635616.99274659, B: 524594.58741211, C: 0.0},
			{A: 639220.51042334, B: 339970.13564543, C: 0.0},
			{A: 630089.46166298, B: 290900.09934345, C: 0.0},
			{A: 663034.77996488, B: 522956.52529761, C: 0.0},
		},
	},
	{
		Name:     "op_coordinate_frame",
		FromEPSG: 4326,
		ToEPSG:   20249,
		CoversCS: "",
		CoversOp: "coordinate_frame",
		In: [4]transformPt{
			{A: 110.8206053, B: -23.59249196, C: 0.0},
			{A: 111.50138786, B: -29.92868515, C: 0.0},
			{A: 112.46027882, B: -30.89009191, C: 0.0},
			{A: 110.25172761, B: -31.31712192, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 481557.84813298, B: 7390730.64384564, C: 0.0},
			{A: 548252.76988067, B: 6688869.11548245, C: 0.0},
			{A: 639433.26825717, B: 6581523.95959043, C: 0.0},
			{A: 428663.03409921, B: 6534869.70897356, C: 0.0},
		},
	},
	{
		Name:     "op_coordinate_frame_full_matrix",
		FromEPSG: 4326,
		ToEPSG:   10744,
		CoversCS: "",
		CoversOp: "coordinate_frame_full_matrix",
		In: [4]transformPt{
			{A: -62.98251728, B: 17.532164, C: 0.0},
			{A: -62.91737731, B: 17.47247075, C: 0.0},
			{A: -62.95056886, B: 17.4847672, C: 0.0},
			{A: -62.91592273, B: 17.498694, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 2373.52905379, B: 7974.29279227, C: 0.0},
			{A: 9285.32165561, B: 1361.18185091, C: 0.0},
			{A: 5761.04359478, B: 2725.09108038, C: 0.0},
			{A: 9441.99562137, B: 4263.38845706, C: 0.0},
		},
	},
	{
		Name:     "op_molodensky_badekas",
		FromEPSG: 4326,
		ToEPSG:   22991,
		CoversCS: "",
		CoversOp: "molodensky_badekas",
		In: [4]transformPt{
			{A: 36.01599202, B: 24.49783543, C: 0.0},
			{A: 36.08984945, B: 24.85405108, C: 0.0},
			{A: 36.07120633, B: 25.34363749, C: 0.0},
			{A: 34.04691404, B: 26.29633222, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 402810.25217394, B: 490680.76973445, C: 0.0},
			{A: 409983.7978105, B: 530199.6386995, C: 0.0},
			{A: 407669.46035977, B: 584424.234617, C: 0.0},
			{A: 204648.30969504, B: 689888.81534154, C: 0.0},
		},
	},
	{
		Name:     "op_molodensky_badekas_pv",
		FromEPSG: 4326,
		ToEPSG:   5456,
		CoversCS: "",
		CoversOp: "molodensky_badekas_pv",
		In: [4]transformPt{
			{A: -83.77828299, B: 10.1860688, C: 0.0},
			{A: -84.03078493, B: 10.38658298, C: 0.0},
			{A: -84.48695663, B: 10.45266028, C: 0.0},
			{A: -84.75625275, B: 10.58670828, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 560605.85980797, B: 240998.05970091, C: 0.0},
			{A: 532917.25849163, B: 263138.70288212, C: 0.0},
			{A: 482968.64051781, B: 270435.82925984, C: 0.0},
			{A: 453506.10705847, B: 285289.62615735, C: 0.0},
		},
	},
	{
		Name:     "op_geographic_offset",
		FromEPSG: 4326,
		ToEPSG:   2000,
		CoversCS: "",
		CoversOp: "geographic_offset",
		In: [4]transformPt{
			{A: -63.15969557, B: 18.27617592, C: 0.0},
			{A: -63.06307864, B: 18.24896476, C: 0.0},
			{A: -63.02644899, B: 18.24252296, C: 0.0},
			{A: -63.09444014, B: 18.16323654, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 277301.11792751, B: 2021306.44450358, C: 0.0},
			{A: 287496.14645938, B: 2018233.71860039, C: 0.0},
			{A: 291364.48889808, B: 2017498.89177616, C: 0.0},
			{A: 284124.01093583, B: 2008768.71872123, C: 0.0},
		},
	},
	{
		Name:     "op_longitude_rotation",
		FromEPSG: 4326,
		ToEPSG:   4807,
		CoversCS: "",
		CoversOp: "longitude_rotation",
		In: [4]transformPt{
			{A: 3.80886786, B: 45.22351981, C: 0.0},
			{A: 0.76106611, B: 48.27759412, C: 0.0},
			{A: 4.29186209, B: 45.04730074, C: 0.0},
			{A: 0.72077656, B: 45.68470138, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 1.47225883, B: 45.2235405, C: 0.0},
			{A: -1.57538486, B: 48.27766984, C: 0.0},
			{A: 1.95523277, B: 45.04731595, C: 0.0},
			{A: -1.61570963, B: 45.68474772, C: 0.0},
		},
	},
	{
		Name:     "op_horizontal_grid",
		FromEPSG: 4326,
		ToEPSG:   27200,
		CoversCS: "",
		CoversOp: "horizontal_grid",
		In: [4]transformPt{
			{A: 171.78205725, B: -42.337089, C: 0.0},
			{A: 169.75832905, B: -41.346497, C: 0.0},
			{A: 175.73931516, B: -39.30745015, C: 0.0},
			{A: 175.46303759, B: -39.7980426, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 2409633.98754907, B: 5873724.32692024, C: 0.0},
			{A: 2238649.31202354, B: 5979314.8603128, C: 0.0},
			{A: 2746237.22613603, B: 6207181.28144994, C: 0.0},
			{A: 2720903.97654087, B: 6153443.90373303, C: 0.0},
		},
	},
	{
		Name:     "op_vertical_grid",
		FromEPSG: 4326,
		ToEPSG:   9518,
		CoversCS: "",
		CoversOp: "vertical_grid",
		In: [4]transformPt{
			{A: -42.99482709, B: 5.17721902, C: 100.0},
			{A: -107.91231703, B: -23.01331858, C: 100.0},
			{A: -15.1441596, B: 8.63835637, C: 100.0},
			{A: 33.41641472, B: -3.78127545, C: 100.0},
		},
		Want: [4]transformPt{
			{A: -42.99482709, B: 5.17721902, C: 77.83541578},
			{A: -107.91231703, B: -23.01331858, C: 100.0},
			{A: -15.1441596, B: 8.63835637, C: 85.04866638},
			{A: 33.41641472, B: -3.78127545, C: 47.60622887},
		},
	},
	{
		Name:     "op_vertical_offset",
		FromEPSG: 4326,
		ToEPSG:   9724,
		CoversCS: "",
		CoversOp: "vertical_offset",
		In: [4]transformPt{
			{A: 13.9187412, B: 37.16766868, C: 100.0},
			{A: 13.98110423, B: 37.89364695, C: 100.0},
			{A: 14.63000977, B: 37.12119411, C: 100.0},
			{A: 13.20043903, B: 37.48631732, C: 100.0},
		},
		Want: [4]transformPt{
			{A: 13.9187412, B: 37.16766868, C: 99.859},
			{A: 13.98110423, B: 37.89364695, C: 99.859},
			{A: 14.63000977, B: 37.12119411, C: 99.859},
			{A: 13.20043903, B: 37.48631732, C: 99.859},
		},
	},
	{
		Name:     "op_time_dependent_position_vector",
		FromEPSG: 4326,
		ToEPSG:   7912,
		CoversCS: "",
		CoversOp: "time_dependent_position_vector",
		In: [4]transformPt{
			{A: 6.31905805, B: 56.00278238, C: 0.0},
			{A: 6.76421631, B: 56.62670721, C: 0.0},
			{A: 6.41470961, B: 55.836961, C: 0.0},
			{A: 5.27791184, B: 55.53663808, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 6.31905805, B: 56.00278238, C: 0.0},
			{A: 6.76421631, B: 56.62670721, C: 0.0},
			{A: 6.41470961, B: 55.836961, C: 0.0},
			{A: 5.27791184, B: 55.53663808, C: 0.0},
		},
	},
	{
		Name:     "op_time_dependent_position_vector_atepoch",
		FromEPSG: 9989,
		ToEPSG:   9069,
		CoversCS: "",
		CoversOp: "time_dependent_position_vector",
		Epoch:    2010.0, // ITRF2020→ETRF2014 reference epoch is 2015.0
		In: [4]transformPt{
			{A: 10.0, B: 55.0, C: 0.0},
			{A: 8.0, B: 56.0, C: 0.0},
			{A: 6.0, B: 55.5, C: 0.0},
			{A: 12.0, B: 57.0, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 9.99999405, B: 54.99999706, C: -0.00421603},
			{A: 7.99999415, B: 55.99999702, C: -0.00418564},
			{A: 5.99999432, B: 55.49999699, C: -0.00421132},
			{A: 11.99999378, B: 56.99999709, C: -0.00412674},
		},
	},
	{
		Name:     "op_time_specific_position_vector",
		FromEPSG: 4326,
		ToEPSG:   5554,
		CoversCS: "",
		CoversOp: "time_specific_position_vector",
		In: [4]transformPt{
			{A: 5.91203583, B: 51.38717035, C: 100.0},
			{A: 5.96647639, B: 51.17957937, C: 100.0},
			{A: 5.92729969, B: 51.46695894, C: 100.0},
			{A: 5.90749917, B: 51.50131879, C: 100.0},
		},
		Want: [4]transformPt{
			{A: 702611.86777109, B: 5696905.87514295, C: 100.0},
			{A: 707332.13516969, B: 5673978.50858792, C: 100.0},
			{A: 703319.03872321, B: 5705819.06106673, C: 100.0},
			{A: 701792.22489856, B: 5709584.40737409, C: 100.0},
		},
	},
	{
		Name:     "op_velocity_grid",
		FromEPSG: 4326,
		ToEPSG:   9062,
		CoversCS: "",
		CoversOp: "velocity_grid",
		In: [4]transformPt{
			{A: 5.98650521, B: 55.86597666, C: 0.0},
			{A: 6.37453923, B: 55.50831723, C: 0.0},
			{A: 6.60231474, B: 56.65506928, C: 0.0},
			{A: 5.05580951, B: 56.13771929, C: 0.0},
		},
		Want: [4]transformPt{
			{A: 5.98649888, B: 55.86597238, C: -0.02258583},
			{A: 6.37453287, B: 55.50831295, C: -0.02181285},
			{A: 6.60231474, B: 56.65506928, C: 0.00007308},
			{A: 5.05580951, B: 56.13771929, C: 0.00007220},
		},
	},
}

func isGeographicEPSG(code int) bool {
	b, err := epsgDir.ReadFile(fmt.Sprintf("epsg/%d.txt", code))
	if err != nil {
		return false
	}
	line, _, _ := strings.Cut(string(b), "\n")
	return strings.Contains(line, "conversion=geographic")
}

func checkHeight(op string) bool {
	return op == "vertical_grid" || op == "vertical_offset" || op == "velocity_grid"
}

func heightTol(op string) float64 {
	if op == "vertical_grid" {
		return 50.0
	}
	return transformTol
}

func transformTolForEPSG(code int) float64 {
	if isGeographicEPSG(code) {
		return transformTolDeg
	}

	return transformTol
}

func makeTransform(t *testing.T, from, to int, epoch float64) Func {
	t.Helper()
	var (
		f   Func
		err error
	)
	if epoch != 0 {
		f, err = TransformAt(from, to, epoch)
	} else {
		f, err = Transform(from, to)
	}
	if err != nil {
		t.Fatalf("Transform(%d,%d): %v", from, to, err)
	}
	return f
}

func checkTransformPt(t *testing.T, label string, i int, gotA, gotB, gotC float64, want transformPt, tol, htol float64, checkC bool) {
	t.Helper()
	bad := math.Abs(gotA-want.A) > tol || math.Abs(gotB-want.B) > tol
	if checkC && math.Abs(gotC-want.C) > htol {
		bad = true
	}
	if bad {
		t.Errorf("%s point %d: got (%.8f, %.8f, %.8f) want (%.8f, %.8f, %.8f)",
			label, i, gotA, gotB, gotC, want.A, want.B, want.C)
	}
}

func TestEPSGTransform(t *testing.T) {
	for _, tc := range transformCases {
		t.Run(tc.Name, func(t *testing.T) {
			fwd := makeTransform(t, tc.FromEPSG, tc.ToEPSG, tc.Epoch)
			rev := makeTransform(t, tc.ToEPSG, tc.FromEPSG, tc.Epoch)

			fwdTol := transformTolForEPSG(tc.ToEPSG)
			revTol := transformTolForEPSG(tc.FromEPSG)
			htol := heightTol(tc.CoversOp)
			checkC := checkHeight(tc.CoversOp) || tc.Epoch != 0

			for i, in := range tc.In {
				want := tc.Want[i]

				gotA, gotB, gotC, err := fwd(in.A, in.B, in.C)
				if err != nil {
					t.Errorf("forward point %d: %v", i, err)
					continue
				}
				checkTransformPt(t, "forward", i, gotA, gotB, gotC, want, fwdTol, htol, checkC)

				// Round-trip on the computed forward result (Hin und zurück).
				backA, backB, backC, err := rev(gotA, gotB, gotC)
				if err != nil {
					t.Errorf("reverse point %d: %v", i, err)
					continue
				}
				checkTransformPt(t, "round-trip", i, backA, backB, backC, in, revTol, htol, checkC)
			}
		})
	}
}
