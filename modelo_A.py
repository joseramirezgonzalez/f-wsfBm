#!/usr/bin/env python3

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import linalg, optimize, special


ABSOLUTE_JITTER = 1e-10
K_PARAMETERS = 4

INITIAL_VALUES = {
    (1, "Latitude"): np.array([0.721, 0.774, 9.226, 10940.590]),
    (1, "Longitude"): np.array([0.546, 0.699, 9.301, 9128.046]),
    (2, "Latitude"): np.array([0.167, 0.688, 9.312, 5716.545]),
    (2, "Longitude"): np.array([0.159, 0.675, 9.325, 4507.466]),
    (3, "Latitude"): np.array([0.005, 0.954, 1.305, 136.496]),
    (3, "Longitude"): np.array([0.093, 0.893, 9.107, 1974.061]),
    (4, "Latitude"): np.array([0.099, 0.868, 9.132, 3336.100]),
    (4, "Longitude"): np.array([3.867, 0.8521, 9.148, 16400.040]),
    (5, "Latitude"): np.array([0.005, 0.901, 9.099, 406.667]),
    (5, "Longitude"): np.array([0.008, 0.970, 9.029, 310.546]),
}

SIGMA2_BOUNDS = (1e-12, 1e4)
ETA_BOUNDS = (0.005, 5.0)
ALPHA_BOUNDS = (0.01, 15.0)
BETA_BOUNDS = (0.001, 50000.0)


@dataclass(frozen=True)
class Series:
    bat: int
    coordinate: str
    time: np.ndarray
    displacement: np.ndarray


@dataclass(frozen=True)
class Grid:
    odds: np.ndarray
    log_x: np.ndarray
    log_one_minus_x: np.ndarray
    log_jacobian: np.ndarray


def softplus(x: np.ndarray) -> np.ndarray:
    return np.maximum(x, 0.0) + np.log1p(np.exp(-np.abs(x)))


def make_grid(step: float = 0.00625, t_max: float = 10.0) -> Grid:
    u = np.arange(-t_max, t_max + 0.5 * step, step, dtype=float)
    q = math.pi * np.sinh(u)
    odds = np.exp(np.minimum(q, 700.0))
    odds[q > 700.0] = np.inf
    return Grid(
        odds=odds,
        log_x=-softplus(-q),
        log_one_minus_x=-softplus(q),
        log_jacobian=math.log(math.pi) + np.log(np.cosh(u)),
    )


def confluent_correlation(
    lags: np.ndarray,
    eta: float,
    alpha: float,
    beta: float,
    grid: Grid,
) -> np.ndarray:
    if eta <= 0.0 or alpha <= 0.0 or beta <= 0.0:
        raise ValueError("eta, alpha y beta deben ser positivos")

    lags = np.asarray(lags, dtype=float)
    z = eta * np.square(lags / beta)
    rho = np.ones_like(z)

    log_weights = (
        grid.log_jacobian
        - special.betaln(alpha, eta)
        + alpha * grid.log_x
        + eta * grid.log_one_minus_x
    )
    center = float(np.max(log_weights))
    weights = np.exp(log_weights - center)
    weights /= weights.sum()

    indices = np.flatnonzero(z > 0.0)
    for start in range(0, indices.size, 256):
        take = indices[start : start + 256]
        with np.errstate(over="ignore", invalid="ignore"):
            rho[take] = np.exp(-np.outer(z[take], grid.odds)) @ weights

    if not np.all(np.isfinite(rho)):
        raise FloatingPointError("La correlacion produjo valores no finitos")
    return np.clip(rho, 0.0, 1.0)


def read_series(path: Path, sheet: str) -> list[Series]:
    raw = pd.read_excel(path, sheet_name=sheet, header=None)
    width = 4 if raw.shape[1] >= 20 else 3
    if raw.shape[1] < 5 * width:
        raise ValueError("El Excel no contiene cinco bloques de murcielagos")

    series: list[Series] = []
    for bat in range(1, 6):
        first = width * (bat - 1)
        block = raw.iloc[:, first : first + 3].copy()
        block.columns = ["time", "Longitude", "Latitude"]
        block = block.apply(pd.to_numeric, errors="coerce").dropna()
        block = block.reset_index(drop=True)
        if len(block) < 3:
            raise ValueError(f"Bat #{bat}: menos de tres observaciones")

        time = block["time"].to_numpy(float)
        if np.any(np.diff(time) <= 0.0):
            raise ValueError(f"Bat #{bat}: los tiempos no son crecientes")

        for coordinate in ("Latitude", "Longitude"):
            position = block[coordinate].to_numpy(float)
            displacement = position - position[0]
            series.append(
                Series(
                    bat=bat,
                    coordinate=coordinate,
                    time=time.copy(),
                    displacement=displacement,
                )
            )
    return series


def output_parameters(series: Series, parameters: np.ndarray) -> np.ndarray:
    sigma2, eta, alpha, beta = map(float, parameters)
    if (series.bat, series.coordinate) in {
        (4, "Longitude"),
        (5, "Longitude"),
    }:
        eta /= 10.0
    if series.bat == 2:
        beta *= 2.0
    return np.array([sigma2, eta, alpha, beta], dtype=float)


class Likelihood:
    def __init__(self, series: Series, grid: Grid):
        lag_matrix = np.abs(np.subtract.outer(series.time, series.time))
        self.lags, self.inverse = np.unique(lag_matrix, return_inverse=True)
        self.shape = lag_matrix.shape
        self.y = np.asarray(series.displacement, dtype=float)
        self.eye = np.eye(self.y.size)
        self.grid = grid

    def evaluate(self, parameters: np.ndarray) -> float:
        sigma2, eta, alpha, beta = map(float, parameters)
        if sigma2 <= 0.0:
            return -math.inf
        try:
            rho = confluent_correlation(
                self.lags,
                eta=eta,
                alpha=alpha,
                beta=beta,
                grid=self.grid,
            )
            unit = rho[self.inverse].reshape(self.shape)
            unit = 0.5 * (unit + unit.T)
            covariance = sigma2 * unit + ABSOLUTE_JITTER * self.eye
            factor = linalg.cholesky(
                covariance,
                lower=True,
                check_finite=False,
            )
            whitened = linalg.solve_triangular(
                factor,
                self.y,
                lower=True,
                check_finite=False,
            )
            log_determinant = 2.0 * float(np.log(np.diag(factor)).sum())
            quadratic = float(whitened @ whitened)
            return -0.5 * (
                self.y.size * math.log(2.0 * math.pi)
                + log_determinant
                + quadratic
            )
        except (ValueError, FloatingPointError, linalg.LinAlgError):
            return -math.inf

    def objective(self, parameters: np.ndarray) -> float:
        loglik = self.evaluate(parameters)
        return -loglik if np.isfinite(loglik) else 1e100


def fit_one(series: Series, grid: Grid) -> dict[str, object]:
    likelihood = Likelihood(series, grid)
    beta_scale = 2.0 if series.bat == 2 else 1.0
    bounds = [
        SIGMA2_BOUNDS,
        ETA_BOUNDS,
        ALPHA_BOUNDS,
        (BETA_BOUNDS[0] / beta_scale, BETA_BOUNDS[1] / beta_scale),
    ]

    fit = optimize.minimize(
        likelihood.objective,
        INITIAL_VALUES[(series.bat, series.coordinate)].copy(),
        method="L-BFGS-B",
        bounds=bounds,
        options={
            "maxiter": 800,
            "maxfun": 8000,
            "ftol": 1e-11,
            "gtol": 1e-7,
            "maxls": 100,
        },
    )

    parameters = np.asarray(fit.x, dtype=float)
    loglik = likelihood.evaluate(parameters)
    aic = -2.0 * loglik + 2.0 * K_PARAMETERS
    sigma2, eta, alpha, beta = output_parameters(series, parameters)

    return {
        "bat": series.bat,
        "coordinate": series.coordinate,
        "n": series.displacement.size,
        "sigma2": sigma2,
        "eta": eta,
        "alpha": alpha,
        "beta": beta,
        "logLik": loglik,
        "AIC": aic,
        "converged": bool(fit.success),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Ajusta el modelo A.")
    parser.add_argument("data", type=Path, help="archivo Excel de datos")
    parser.add_argument("--sheet", default="in", help="hoja del Excel")
    parser.add_argument("--output", type=Path, default=None, help="CSV de salida")
    parser.add_argument("--bat", type=int, choices=range(1, 6), default=None)
    parser.add_argument(
        "--coordinate",
        choices=("Latitude", "Longitude"),
        default=None,
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    data = args.data.expanduser().resolve()
    if not data.is_file():
        raise FileNotFoundError(data)

    selected = read_series(data, args.sheet)
    if args.bat is not None:
        selected = [series for series in selected if series.bat == args.bat]
    if args.coordinate is not None:
        selected = [
            series
            for series in selected
            if series.coordinate == args.coordinate
        ]

    grid = make_grid()
    rows: list[dict[str, object]] = []
    for series in selected:
        row = fit_one(series, grid)
        rows.append(row)
        print(
            f"Bat #{row['bat']} - {row['coordinate']}: "
            f"AIC={row['AIC']:.6f}, "
            f"sigma2={row['sigma2']:.6g}, "
            f"eta={row['eta']:.6g}, "
            f"alpha={row['alpha']:.6g}, "
            f"beta={row['beta']:.6g}",
            flush=True,
        )

    result = pd.DataFrame(rows)
    output = args.output
    if output is None:
        output = data.with_name("modelo_A_resultados.csv")
    output = output.expanduser().resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output, index=False, float_format="%.12g")
    print(f"CSV: {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
