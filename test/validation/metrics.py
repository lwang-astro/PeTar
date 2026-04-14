import math
from typing import Dict, List


def check_max_threshold(value: float, threshold: float) -> Dict[str, object]:
    passed = abs(value) <= threshold
    return {
        "passed": passed,
        "value": value,
        "threshold": threshold,
        "message": f"|value|={abs(value):.6e} <= {threshold:.6e}" if passed else f"|value|={abs(value):.6e} > {threshold:.6e}",
    }


def check_convergence_ratio(
    coarse: float,
    fine: float,
    expected_ratio: float,
    ratio_tolerance: float,
) -> Dict[str, object]:
    if fine == 0.0:
        return {
            "passed": False,
            "value": math.inf,
            "expected": expected_ratio,
            "tolerance": ratio_tolerance,
            "message": "fine error is zero, cannot evaluate convergence ratio",
        }

    ratio = abs(coarse / fine)
    lower = expected_ratio * (1.0 - ratio_tolerance)
    upper = expected_ratio * (1.0 + ratio_tolerance)
    passed = lower <= ratio <= upper
    return {
        "passed": passed,
        "value": ratio,
        "expected": expected_ratio,
        "tolerance": ratio_tolerance,
        "message": (
            f"ratio={ratio:.4f} within [{lower:.4f}, {upper:.4f}]"
            if passed
            else f"ratio={ratio:.4f} outside [{lower:.4f}, {upper:.4f}]"
        ),
    }


def check_loglog_slope(
    x_values: List[float],
    y_values: List[float],
    expected_slope: float,
    slope_tolerance: float,
) -> Dict[str, object]:
    if len(x_values) != len(y_values) or len(x_values) < 2:
        return {
            "passed": False,
            "value": math.nan,
            "expected": expected_slope,
            "tolerance": slope_tolerance,
            "message": "need at least two valid points for slope fit",
        }

    if any(x <= 0.0 for x in x_values) or any(y <= 0.0 for y in y_values):
        return {
            "passed": False,
            "value": math.nan,
            "expected": expected_slope,
            "tolerance": slope_tolerance,
            "message": "log-log slope requires positive x/y values",
        }

    lx = [math.log10(v) for v in x_values]
    ly = [math.log10(v) for v in y_values]
    mx = sum(lx) / len(lx)
    my = sum(ly) / len(ly)
    denom = sum((x - mx) ** 2 for x in lx)
    if denom == 0.0:
        return {
            "passed": False,
            "value": math.nan,
            "expected": expected_slope,
            "tolerance": slope_tolerance,
            "message": "zero variance in x values",
        }

    slope = sum((x - mx) * (y - my) for x, y in zip(lx, ly)) / denom
    lower = expected_slope - slope_tolerance
    upper = expected_slope + slope_tolerance
    passed = lower <= slope <= upper
    return {
        "passed": passed,
        "value": slope,
        "expected": expected_slope,
        "tolerance": slope_tolerance,
        "message": (
            f"slope={slope:.4f} within [{lower:.4f}, {upper:.4f}]"
            if passed
            else f"slope={slope:.4f} outside [{lower:.4f}, {upper:.4f}]"
        ),
    }


def summarize_results(results: List[Dict[str, object]]) -> Dict[str, object]:
    failed = [item for item in results if not item.get("passed", False)]
    return {
        "passed": len(failed) == 0,
        "n_total": len(results),
        "n_failed": len(failed),
    }