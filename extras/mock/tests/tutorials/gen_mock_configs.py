#!/usr/bin/env python3
"""Generate per-tutorial precice-mock-config.xml with bounded random read-data
so solvers that divide by mock data don't hit FPE/NaN with zero buffers."""
import re
import sys
from pathlib import Path

TUT = Path(__file__).resolve().parents[4] / "examples" / "tutorials"

# Tutorials that fail with numerical errors (FPE 136 / NaN 134) -> the cfg dir
TARGETS = [
    "breaking-dam-2d",
    "elastic-tube-3d",
    "flow-around-controlled-moving-cylinder",
    "flow-over-heated-plate-nearest-projection",
    "flow-over-heated-plate-partitioned-flow",
    "multiple-perpendicular-flaps",
    "partitioned-backwards-facing-step",
    "quickstart",
    "turek-hron-fsi3",
    "perpendicular-flap",
    "elastic-tube-1d",
]


# Hand-tuned overrides where the generic heuristics produce unphysical boundary
# data: (tutorial, mesh, data) -> (lower, upper).
OVERRIDES = {
    # buoyantPimpleFoam runs at an absolute pressure of ~103500 Pa; feeding the
    # coupled outlet a pressure near 0 lets the equation of state divide by zero
    # (SIGFPE). Keep the mocked pressure tight around the operating point.
    ("flow-over-heated-plate-partitioned-flow", "Fluid1-Fluid-Mesh", "Pressure"): (
        103499.999,
        103500.001,
    ),
    # Fully developed outflow: near-zero gradients keep the outlet BC well posed.
    (
        "flow-over-heated-plate-partitioned-flow",
        "Fluid1-Fluid-Mesh",
        "VelocityGradient",
    ): (-1e-6, 1e-6),
    (
        "flow-over-heated-plate-partitioned-flow",
        "Fluid1-Fluid-Mesh",
        "FlowTemperatureGradient",
    ): (-1e-6, 1e-6),
    # Spatially random heat flux slowly destabilizes the buoyant solver
    # (blowup around t=0.7 with +-10 W/m^2); keep the plate near-adiabatic.
    ("flow-over-heated-plate-partitioned-flow", "Fluid1-Solid-Mesh", "Heat-Flux"): (
        -1e-3,
        1e-3,
    ),
    # DisplacementDelta accumulates over windows on a 5 mm radius tube; white
    # +-1e-5 noise per window roughens the wall until cells invert (GAMG FPE).
    ("elastic-tube-3d", "Fluid-Mesh-Nodes", "DisplacementDelta"): (-1e-8, 1e-8),
    # Incompressible pimpleFoam: spatially white Dirichlet data on the coupled
    # outlet makes the pressure equation unsolvable (fluxes don't balance).
    # Near-zero values turn the outlet into a plain p=0 / zero-gradient outflow.
    ("partitioned-backwards-facing-step", "Fluid1-Mesh", "Pressure"): (-1e-6, 1e-6),
    ("partitioned-backwards-facing-step", "Fluid1-Mesh", "VelocityBack"): (-1e-6, 1e-6),
    ("partitioned-backwards-facing-step", "Fluid1-Mesh", "VelocityGradient"): (
        -1e-6,
        1e-6,
    ),
}


def bounds_for(name):
    """Return (lower, upper) for a data name, by keyword. Order matters."""
    n = name.lower()
    if "crosssectionlength" in n:
        return (0.9, 1.1)  # must stay positive ~1 (area term)
    if "gradient" in n:
        return (-1.0, 1.0)  # any gradient: small
    if "displacement" in n:
        return (-1e-5, 1e-5)  # mesh motion: tiny to avoid tangling
    if "temperature" in n:
        return (290.0, 310.0)  # physical, nonzero
    if "heat-flux" in n or "heatflux" in n:
        return (-10.0, 10.0)
    if "stress" in n:
        return (-1e-2, 1e-2)
    if "force" in n:
        return (-1e-2, 1e-2)
    if "pressure" in n:
        return (-1.0, 1.0)
    if "velocity" in n:
        return (-0.1, 0.1)
    return (-1e-3, 1e-3)  # conservative default


def main():
    # optional CLI filter: regenerate only the named tutorials
    targets = sys.argv[1:] if len(sys.argv) > 1 else TARGETS
    for t in targets:
        cfg = TUT / t / "precice-config.xml"
        if not cfg.exists():
            print(f"!! no config for {t}")
            continue
        text = cfg.read_text()
        pairs = sorted(
            set(re.findall(r'read-data\s+name="([^"]+)"\s+mesh="([^"]+)"', text))
        )
        if not pairs:
            print(f"!! no read-data in {t}")
            continue
        lines = [
            '<?xml version="1.0" encoding="UTF-8" ?>',
            "<mock-config>",
            "  <!-- Auto-generated bounded random read-data so the solver does not",
            "       divide by zero / produce NaN when run against the mock alone. -->",
            '  <logging-mode mode="mock" />',
            '  <max-iterations-override value="2" />',
        ]
        seed = 1
        for name, mesh in pairs:
            lo, hi = OVERRIDES.get((t, mesh, name)) or bounds_for(name)
            lines.append(f'  <mocked-data mesh="{mesh}" data="{name}" mode="random">')
            lines.append(f'    <bounds lower="{lo}" upper="{hi}" />')
            lines.append(f'    <seed value="{seed}" />')
            lines.append("  </mocked-data>")
            seed += 1
        lines.append("</mock-config>")
        out = TUT / t / "precice-mock-config.xml"
        out.write_text("\n".join(lines) + "\n")
        print(f"wrote {out}  ({len(pairs)} data items)")


if __name__ == "__main__":
    main()
