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


def summarize_results(results: List[Dict[str, object]]) -> Dict[str, object]:
    failed = [item for item in results if not item.get("passed", False)]
    return {
        "passed": len(failed) == 0,
        "n_total": len(results),
        "n_failed": len(failed),
    }