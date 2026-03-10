import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd
import pathlib, os, sys, requests
from datetime import datetime, timedelta

path = pathlib.Path().resolve()
os.chdir(path)

# ──────────────────────────────────────────────────────────────────────────────
# The following functions represent S'(t), I'(t) and R'(t) of the SIR model. 
# These are numerically solved using the 4th-order Runge-Kutta model.
# ──────────────────────────────────────────────────────────────────────────────

def dS(s: float, i: float, beta: float) -> float:
    """S'(t) = -β*S*I"""
    return -beta * s * i

def dI(s: float, i: float, beta: float, gamma: float) -> float:
    """I'(t) = β*S*I - γ*I"""
    return beta * s * i - gamma * i

def dR(i: float, gamma: float) -> float:
    """R'(t) = γ*I"""
    return gamma * i

# ──────────────────────────────────────────────────────────────────────────────
# COUPLED RK4 SOLVER
# ──────────────────────────────────────────────────────────────────────────────

def solve_sir_rk4(
    S0: float, I0: float, R0: float,
    beta: float, gamma: float,
    days: int = 900, h: float = 0.05
) -> pd.DataFrame:
    """
    Solve the SIR system with RK4 on the coupled vector [S, I, R].

    Parameters
    ----------
    S0, I0, R0 : initial conditions
    beta       : transmission rate (per person per day)
    gamma      : recovery rate (per day)
    days       : forecast horizon in days
    h          : RK4 step size (days)

    Returns
    -------
    DataFrame with integer-day snapshots indexed by step number.
    """
    n_steps = int(days / h)
    S = np.zeros(n_steps)
    I = np.zeros(n_steps)
    R = np.zeros(n_steps)

    S[0], I[0], R[0] = S0, I0, R0

    for n in range(n_steps - 1):
        sn, i_n, rn = S[n], I[n], R[n]

        # --- Stage 1 ---
        ks1 = dS(sn, i_n, beta)
        ki1 = dI(sn, i_n, beta, gamma)
        kr1 = dR(i_n, gamma)

        # --- Stage 2  (use intermediate estimates for ALL three variables) ---
        s2 = sn + 0.5 * h * ks1
        i2 = i_n + 0.5 * h * ki1
        ks2 = dS(s2, i2, beta)
        ki2 = dI(s2, i2, beta, gamma)
        kr2 = dR(i2, gamma)

        # --- Stage 3 ---
        s3 = sn + 0.5 * h * ks2
        i3 = i_n + 0.5 * h * ki2
        ks3 = dS(s3, i3, beta)
        ki3 = dI(s3, i3, beta, gamma)
        kr3 = dR(i3, gamma)

        # --- Stage 4 ---
        s4 = sn + h * ks3
        i4 = i_n + h * ki3
        ks4 = dS(s4, i4, beta)
        ki4 = dI(s4, i4, beta, gamma)
        kr4 = dR(i4, gamma)

        S[n + 1] = sn + (h / 6) * (ks1 + 2*ks2 + 2*ks3 + ks4)
        I[n + 1] = i_n + (h / 6) * (ki1 + 2*ki2 + 2*ki3 + ki4)
        R[n + 1] = rn + (h / 6) * (kr1 + 2*kr2 + 2*kr3 + kr4)

    steps = np.arange(n_steps)
    day_vals = steps * h

    df = pd.DataFrame({"day": day_vals, "Susceptible": S, "Infected": I, "Recovered": R})
    # Keep only integer-day rows
    df = df[np.isclose(df["day"] % 1, 0)].copy()
    df["day"] = df["day"].astype(int)

    # Conservation check: S + I + R ≈ N
    N = S0 + I0 + R0
    drift = abs((df["Susceptible"] + df["Infected"] + df["Recovered"]).max() - N)
    if drift > 1.0:
        print(f"  ⚠  Population conservation drift detected: {drift:.2f} people")
    else:
        print(f"  ✓  Population conserved (max drift = {drift:.4f})")

    return df

# ──────────────────────────────────────────────────────────────────────────────
# DATA LOADERS
# ──────────────────────────────────────────────────────────────────────────────

def load_local_data(csv_path: str) -> dict:
    """Load the original Querétaro CSV and extract SIR parameters."""
    df = pd.read_csv(csv_path).fillna(0)
    df["Date"] = pd.to_datetime(df["Date"], format="%d/%m/%Y")
    df["Beta"]  = df["Beta"].astype(str).str.replace(",", ".").astype(float)
    df["Gamma"] = df["Gamma"].astype(str).str.replace(",", ".").astype(float)
    df = df.sort_values("Date").reset_index(drop=True)

    return {
        "beta":  df["Beta"].mean(),
        "gamma": df["Gamma"].mean(),
        "S0":    df["Susceptible"].iloc[0],
        "I0":    df["Infected"].iloc[0],
        "R0":    df["Recovered_SIR"].iloc[0],
        "start_date": df["Date"].iloc[0],
        "source": "Local CSV – Querétaro",
    }

def load_mexico_dge(state_name: str = "QUERETARO") -> dict:
    """
    Load Mexico's official DGE open COVID data.
    Dataset: https://datos.gob.mx/busca/dataset/informacion-referente-a-casos-covid-19-en-mexico

    The daily summary CSV (~10 MB) is hosted on the DGE S3 bucket.
    state_name: Mexican state name in upper-case (e.g. 'JALISCO', 'CDMX').

    Returns a dict compatible with solve_sir_rk4().
    """
    # Official DGE summary (updated daily, publicly accessible)
    DGE_URL = (
        "https://datosabiertos.salud.gob.mx/gobmx/salud/datos_abiertos/"
        "historicos/2024/COVID19MEXICO2024.csv.zip"
    )
    print(f"  Downloading DGE dataset … (this may take a moment)")
    df_raw = pd.read_csv(DGE_URL, compression="zip",
                         usecols=["ENTIDAD_RES", "FECHA_SINTOMAS",
                                  "RESULTADO_LAB", "FECHA_DEF"],
                         encoding="latin-1",
                         low_memory=False)

    # State mapping (ENTIDAD_RES is a 2-digit INEGI code)
    STATE_CODES = {
        "AGUASCALIENTES": "01", "BAJA CALIFORNIA": "02", "BAJA CALIFORNIA SUR": "03",
        "CAMPECHE": "04", "COAHUILA": "05", "COLIMA": "06", "CHIAPAS": "07",
        "CHIHUAHUA": "08", "CDMX": "09", "DURANGO": "10", "GUANAJUATO": "11",
        "GUERRERO": "12", "HIDALGO": "13", "JALISCO": "14", "MEXICO": "15",
        "MICHOACAN": "16", "MORELOS": "17", "NAYARIT": "18", "NUEVO LEON": "19",
        "OAXACA": "20", "PUEBLA": "21", "QUERETARO": "22", "QUINTANA ROO": "23",
        "SAN LUIS POTOSI": "24", "SINALOA": "25", "SONORA": "26", "TABASCO": "27",
        "TAMAULIPAS": "28", "TLAXCALA": "29", "VERACRUZ": "30", "YUCATAN": "31",
        "ZACATECAS": "32",
    }
    code = int(STATE_CODES.get(state_name.upper(), "22"))
    df_state = df_raw[df_raw["ENTIDAD_RES"] == code].copy()

    df_state["FECHA_SINTOMAS"] = pd.to_datetime(
        df_state["FECHA_SINTOMAS"], errors="coerce"
    )
    # Confirmed positives (RESULTADO_LAB == 1)
    positives = (
        df_state[df_state["RESULTADO_LAB"] == 1]
        .groupby("FECHA_SINTOMAS")
        .size()
        .sort_index()
        .cumsum()
    )
    # Deaths
    deaths = (
        df_state[df_state["FECHA_DEF"].notna() & (df_state["FECHA_DEF"] != "9999-99-99")]
        .groupby("FECHA_SINTOMAS")
        .size()
        .sort_index()
        .cumsum()
    )

    # State populations (2020 census, INEGI)
    STATE_POP = {
        "22": 2_368_467, "09": 9_209_944, "14": 8_348_151, "19": 5_784_442,
        "15": 16_992_418, "21": 6_583_278, "30": 8_112_505, "11": 6_166_934,
    }
    population = int(STATE_POP.get(str(code).zfill(2), 2_000_000))

    latest_cases  = int(positives.iloc[-1]) if len(positives) else 1
    latest_deaths = int(deaths.iloc[-1])    if len(deaths)    else 0
    # Assume 80% of confirmed cases resolved (proxy for recoveries)
    latest_rec    = int(latest_cases * 0.80)

    S0 = population - latest_cases
    I0 = max(latest_cases - latest_rec - latest_deaths, 1)
    R0 = latest_rec + latest_deaths

    # Rough β/γ estimates from the final 30-day window
    daily_new = positives.diff().clip(lower=0).fillna(0).iloc[-30:]
    active_ts  = (positives - deaths * 0.8 - positives * 0.8).clip(lower=0).iloc[-30:]
    S_ts       = population - positives.iloc[-30:]
    beta_est   = float(
        (daily_new / ((S_ts / population) * active_ts).replace(0, np.nan))
        .dropna().median()
    )
    gamma_est  = 1 / 14  # typical COVID infectious period ≈ 14 days
    beta_est   = np.clip(beta_est, 0.0, 2.0)

    return {
        "beta":  beta_est,
        "gamma": gamma_est,
        "S0": S0, "I0": I0, "R0": R0,
        "start_date": positives.index[-1] if len(positives) else datetime.today(),
        "source": f"DGE datos abiertos – {state_name.title()}",
    }

# ──────────────────────────────────────────────────────────────────────────────
# PLOTTING
# ──────────────────────────────────────────────────────────────────────────────

def plot_sir(df_sir: pd.DataFrame, params: dict, location: str = ""):
    """Plot the SIR curves with peak annotation."""
    peak = df_sir.loc[df_sir["Infected"].idxmax()]
    N    = params["S0"] + params["I0"] + params["R0"]
    final= df_sir.loc[df_sir["Susceptible"].idxmin()]
    percentage = (1 - final["Susceptible"] / N) * 100

    fig, ax = plt.subplots(figsize=(12, 6))

    ax.plot(df_sir.index, df_sir["Susceptible"] / 1e6, color="green",
            label="Susceptible", linewidth=2)
    ax.plot(df_sir.index, df_sir["Infected"]    / 1e6, color="red",
            label="Infected",    linewidth=2)
    ax.plot(df_sir.index, df_sir["Recovered"]   / 1e6, color="blue",
            label="Recovered",   linewidth=2)

    ax.axvline(peak.name, color="orange", linestyle="--", linewidth=1.5,
               label=f"Peak: {peak.name.strftime('%d %b %Y')}")

    # Annotations
    ax.annotate(
        f"Peak infections\n{round(peak['Infected'] / 1e3):,}k people\n{peak.name.strftime('%d %b %Y')}",
        xy=(peak.name, peak["Infected"] / 1e6),
        xytext=(20, 20), textcoords="offset points",
        arrowprops=dict(arrowstyle="->", color="black"),
        fontsize=9, bbox=dict(boxstyle="round,pad=0.3", fc="lightyellow")
    )

    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")

    title = f"SIR Model – {location or params['source']}"
    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("Date")
    ax.set_ylabel("Millions of People")
    ax.legend()
    ax.grid(alpha=0.3)

    print(f"\n{'─'*60}")
    print(f"  Source  : {params['source']}")
    print(f"  β       : {params['beta']:.6f}  |  γ : {params['gamma']:.6f}")
    print(f"  R₀ = β/γ: {params['beta']/params['gamma']:.2f}")
    print(f"  Peak    : {peak.name.strftime('%d %b %Y')} "
          f"({round(peak['Infected']):,} infected)")
    print(f"  Total infected by end: {percentage:.2f}% of population")
    print(f"{'─'*60}\n")

    plt.tight_layout()
    plt.show()
    return fig


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def run_sir(params: dict, days: int = 900, h: float = 0.05, location: str = ""):
    """End-to-end: solve SIR and plot, given a params dict."""
    print(f"\n→ Running SIR for: {params['source']}")

    df_sir = solve_sir_rk4(
        S0=params["S0"], I0=params["I0"], R0=params["R0"],
        beta=params["beta"], gamma=params["gamma"],
        days=days, h=h
    )

    # Attach real dates to the integer-day index
    start = pd.Timestamp(params["start_date"])
    df_sir.index = [start + timedelta(days=int(d)) for d in df_sir["day"]]
    df_sir.index.name = "Date"
    df_sir = df_sir.drop("day", axis=1)

    fig = plot_sir(df_sir, params, location=location)
    return df_sir, fig


# ──────────────────────────────────────────────────────────────────────────────
# USAGE
# ──────────────────────────────────────────────────────────────────────────────
# ── A) Original local data (Querétaro CSV) ──────────────────────────────
LOCAL_CSV = r"data\Covid-19 Queretaro.csv"
if pathlib.Path(LOCAL_CSV).exists():
    params_local = load_local_data(LOCAL_CSV)
    run_sir(params_local, location="Querétaro (local CSV)")

# ── D) Mexico official DGE data — any state ──────────────────────────────
# params_jal = load_mexico_dge("JALISCO")
# run_sir(params_jal, location="Jalisco (DGE)")
#
# params_cdmx = load_mexico_dge("CDMX")
# run_sir(params_cdmx, location="CDMX (DGE)")