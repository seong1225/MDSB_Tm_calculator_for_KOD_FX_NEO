import streamlit as st
import numpy as np
import pandas as pd

# ----------------------------
# Thermodynamic parameters
# ----------------------------
nn_params = {
    'AA': (-9100, -24.0), 'TT': (-9100, -24.0),
    'AT': (-8600, -23.9),
    'TA': (-6000, -16.9),
    'CA': (-5800, -12.9), 'TG': (-5800, -12.9),
    'GT': (-6500, -17.3), 'AC': (-6500, -17.3),
    'CT': (-7800, -20.8), 'AG': (-7800, -20.8),
    'GA': (-5600, -13.5), 'TC': (-5600, -13.5),
    'CG': (-11900, -27.8),
    'GC': (-11100, -26.7),
    'GG': (-11000, -26.6), 'CC': (-11000, -26.6),
}

# ----------------------------
# Constants
# ----------------------------
R = 1.987  # cal/K·mol
default_na = 0.05         # 50 mM
default_conc = 0.5e-6     # 0.5 µM
initiation_penalty = -10.8  # cal/K·mol

# ----------------------------
# Streamlit UI
# ----------------------------
st.title("🧬 MDSB Primer Tm Calculator (Nearest neighbor method)")
st.markdown("여러 개의 시퀀스를 줄 단위로 입력하세요. (각 줄 = 하나의 primer)")
st.markdown("Breslauer et al. (1986) ΔH/ΔS 파라미터 기반 계산 (Na⁺ 50 mM, Oligo 0.5 µM for KOD FX NEO polymerase)")
st.markdown("⚠️주의⚠️")
st.markdown("사용하고자 하는 polymerase 제품마다 계산에 반영하는 parameter, 고려하는 요소가 다를 수 있습니다. 본 calculator를 활용하기 전에, 제품의 manual을 꼭 확인하시기 바랍니다 :)")

# Input area
seq_input = st.text_area("Enter one DNA sequence per line", height=200,
                         value="TGCGGCTAGCTAGCATAACCCCTTG\nGCGCCAGGGAGTGTCCAACTTATC")

na_conc = st.number_input("Na⁺ Concentration (M)", value=default_na, step=0.01, format="%.2f")
oligo_conc = st.number_input("Oligonucleotide Concentration (M)", value=default_conc,
                             min_value=1e-9, step=1e-7, format="%.1e")

# ----------------------------
# Tm calculation function
# ----------------------------
def calculate_tm(seq: str, na: float, conc: float) -> float:
    seq = seq.upper()
    delta_h = 0.0
    delta_s = 0.0

    for i in range(len(seq) - 1):
        pair = seq[i:i+2]
        if pair in nn_params:
            dh, ds = nn_params[pair]
        else:
            rc_pair = pair.translate(str.maketrans('ATGC', 'TACG'))[::-1]
            dh, ds = nn_params.get(rc_pair, (0, 0))
        delta_h += dh
        delta_s += ds

    delta_s += initiation_penalty
    na_correction = 16.6 * np.log10(na)
    tm_k = delta_h / (delta_s + R * np.log(conc / 4))
    tm_c = tm_k - 273.15 + na_correction
    return round(tm_c, 2)

# ----------------------------
# Output
# ----------------------------
if st.button("Calculate Tm for all sequences"):
    if oligo_conc <= 0:
        st.error("❗ Oligonucleotide concentration must be greater than 0.")
    else:
        results = []
        lines = seq_input.strip().splitlines()
        for idx, line in enumerate(lines, 1):
            seq = ''.join([c for c in line.strip().upper() if c in "ATGC"])
            if len(seq) < 2:
                results.append((f"Seq {idx}", seq, "Invalid (too short)"))
            else:
                try:
                    tm = calculate_tm(seq, na_conc, oligo_conc)
                    results.append((f"Seq {idx}", seq, f"{tm} °C"))
                except Exception as e:
                    results.append((f"Seq {idx}", seq, f"Error: {str(e)}"))

        df = pd.DataFrame(results, columns=["ID", "Sequence", "Calculated Tm"])
        st.markdown(f"**입력 조건:** Na⁺ = `{na_conc}` M, Oligo = `{oligo_conc:.2e}` M")
        st.dataframe(df, use_container_width=True)


st.markdown("")
st.markdown("POSTECH MDSB lab., Developed by SH Nam / Reviewed by DY Baek")
