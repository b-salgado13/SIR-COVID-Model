# Metadata for the Dataset of COVID-19 Querétaro

---
*March 9th, 2026 (originally obtained on November 25th, 2020)*

## General Information

- **Name of the document:** `Covid-19 Queretaro.csv`

- **Description:** The `.csv` file contains the general information and calculations to approximate the value of $\beta$ and $\gamma$ parameters for the SIR model based on the data of susceptible, infected and recovered people due to COVID-19 pandemic in Querétaro, México.

- **Source:** The dataset was curated from the [public archive](https://gobqro.tumblr.com/archive/2020/10) of the Government of the State of Querétaro, covering 6 months and 22 days (April 4, 2020 - October 27, 2020).

- **Initial Population:** According to INEGI, the 2020 projected population for the state of Querétaro is 2279950.

---

## File structure
- **Number of rows:** 204 excluding the header.
- **Number of fields:** 15 columns.
- **Field delimiter:** Comma `,`.
- **Index:** The `Date` column.

## Field descriptions

|Field Name|Type|Description|
|:---|:---|:---|
|Date|DATE|Date in format `DD/MM/YYYY`.|
|Susceptible|INT|The number of susceptible people. <br/>Defined as `2279950 - Cases` for each day, since 2279950 is considered the initial population.|
|Cases|INT|The accumulated cases.|
|Hospital|INT|Number of infected individuals currently hospitalized.|
|Deaths|INT|Cumulative number of fatalities due to COVID-19.|
|Home|INT|Number of infected individuals in home isolation.|
|ICU|INT|Number of infected people who are in the Intensive Care Unit|
|Discharged|INT|Daily recovery count.<br/>For a given day `d` it is defined as `Recovered(d)-Recovered(d-1)`|
|Recovered|INT|Number of people who have recovered from COVID-19|
|Increased_Cases|INT|New cases since the previous day. For a given day `d` it is defined as `Cases(d)-Cases(d-1)`|
|Infected|INT|Number of people currently infected. Defined as `Hospital + Home + ICU`|
|Beta|FLOAT|Daily transmission rate estimated via finite differences. <br/>For a given day `d` it is defined as $\beta\approx\frac{S(d)-S(d+1)}{S(d)\cdot I(d)}$|
|Gamma|FLOAT|Daily recovery/removal rate estimated via finite differences. <br/>For a given day `d` it is defined as $\gamma\approx\frac{R(d+1)-R(d)}{I(d)}$|
|Recovered_SIR|INT|The recovered persons for SIR model. Defined as `Recovered + Deaths`| 
|Death_rate|FLOAT|Percentage of people who died compared to those who recovered `100*Deaths/Recovered_SIR`|

## Data Limitations
It is important to point out some clear limitations for the usage of this dataset:

### 1. The "Closed Population" Assumption

The SIR model assumes $N$ is constant ($N = S + I + R$).

The population of **2,279,950** does not account for births, natural deaths, or migration into and out of Querétaro during those 6 months. In a real-world scenario, mobility between states can significantly alter the susceptible pool.

### 2. Sensitivity of Finite Difference Estimations

The calculation of $\beta$ and $\gamma$ is based on daily snapshots:


$$\beta \approx \frac{S_t - S_{t+1}}{S_t \cdot I_t}$$

This method is highly sensitive to reporting noise. In 2020, government data often suffered from "weekend effects" (lower reporting on Saturdays/Sundays and spikes on Mondays). This can lead to non-physical fluctuations or even negative values in $\beta$ and $\gamma$ estimates if the cumulative counts weren't perfectly monotonic.

### 3. Under-reporting and Testing Bias

The `Cases` and `Infected` numbers only reflect individuals who sought medical attention and tested positive. Early in the pandemic (2020), testing capacity in Mexico was limited. The "true" number of infected individuals was likely higher, which means your $S$ (Susceptible) count is likely an overestimate.

### 4. Constant vs. Time-Varying Parameters

The standard SIR model often assumes $\beta$ and $\gamma$ are constants. However, this dataset calculates them daily. While this shows the evolution, a simple RK4 prediction using a single average $\beta$ may fail to account for behavioral changes (e.g., lockdowns, mask mandates) that occurred during that period in Querétaro.

### 5. Recovery Definition ($R$)

The `Recovered_SIR` includes deaths. While mathematically correct for the "Removed" category in SIR, it assumes that "Recovered" individuals have permanent immunity. In reality, waning immunity and reinfections (though less common in 2020 than with later variants) are not captured by a basic SIR model.