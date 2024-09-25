import streamlit as st
import numpy as np
import pandas as pd
from pyXSteam.XSteam import XSteam

# Configuración global
steamTable = XSteam(XSteam.UNIT_SYSTEM_MKS)  # Sistema métrico
atmosPress = 13.2 / 14.50377  # Presión atmosférica Cali en bar
max_iter = 100

# Funciones auxiliares
def calcLiquidEntalpy(temp, mass_frac, purity=100):
    # Calcula la entalpía del líquido
    Cp_sol = 1 - (0.6 - 0.0018 * temp + 0.0008 * (100 - purity)) * mass_frac
    Cp_sol = Cp_sol * 4.184  # Conversión a kJ/kg°C
    return Cp_sol * temp

def calcBPE(temp, mass_frac, purity=100):
    # Calcula el aumento del punto de ebullición (BPE)
    B = mass_frac * 100
    return 2 * B / (100 - B)

def calcU(effect, temp, mass_frac, purity=100):
    return 2.500 - effect * 1.500 / num_effects

# Interfaz de usuario
st.title("Simulación de Evaporador Multiefecto")
st.write("Aplicación para el cálculo de evaporadores en trenes de múltiple efecto en flujo en cocorriente para concentrar una solución agua-azúcar.")
st.image("Imagen1.png", caption="Evaporador multiple con alimentación hacia adelante")

# Sección de entrada de parámetros
st.sidebar.header("Parámetros del evaporador")

num_effects = st.sidebar.slider("Número de efectos", min_value=2, max_value=10, value=2)
press_last_effect_vac = st.sidebar.number_input("Presión del último efecto (vacío)", value=20)
P_units = st.sidebar.selectbox("Unidades de presión", ["inHg", "bar", "psi", "Pa"])

feed_flow = st.sidebar.number_input("Flujo de alimentación (kg/s)", value=100.0)
feed_mass_frac = st.sidebar.number_input("Fracción másica de sólidos en alimentación", value=0.1)
feed_temp = st.sidebar.number_input("Temperatura de alimentación (°C)", value=25.0)

prod_mass_frac = st.sidebar.number_input("Fracción másica de sólidos en producto", value=0.5)

steam_temp = st.sidebar.number_input("Temperatura del vapor de calderas (°C)", value=125.0)

# Conversión de unidades de presión
if P_units == "psi":
    press_last_effect_vac = press_last_effect_vac / 14.5038 * 1E5
elif P_units == "Pa":
    press_last_effect_vac = press_last_effect_vac / 1E5
elif P_units == "inHg":
    press_last_effect_vac = press_last_effect_vac / 29.52999

press_last_effect = atmosPress - press_last_effect_vac  # Presión del último efecto en bar

def calcular_evaporadores():

    # Cálculos iniciales
    feed_entalpy = calcLiquidEntalpy(feed_temp, feed_mass_frac)
    steam_entalpy = steamTable.hV_t(steam_temp)
    steam_lambda = steam_entalpy - steamTable.hL_t(steam_temp)

    L = []
    V = []

    for i in range(0,num_effects+1):
        # alimentación hacia adelante
        if i == 0:
            L.append({"name":"F", "mass_flow":feed_flow, "mass_frac": feed_mass_frac, "temp":feed_temp, "entalpy":feed_entalpy})
            V.append({"name":"S", "mass_flow":None, "temp":steam_temp, "entalpy":steam_entalpy, "entalpy_sat":steam_entalpy, "lambda":steam_lambda})
        else:
            L.append({"name":"L"+str(i), "mass_flow":None, "mass_frac": None, "temp":None, "BPE":None, "entalpy":None, "entalpy_sat":None, "U": None})
            V.append({"name":"V"+str(i), "mass_flow":None, "temp":None, "entalpy":None, "entalpy_sat":None, "lambda":None})

    L[num_effects]["mass_frac"] = prod_mass_frac
    L[num_effects]["mass_flow"] = L[0]["mass_flow"]*L[0]["mass_frac"]/L[num_effects]["mass_frac"]
    L[num_effects]["temp_sat_pure"] = steamTable.tsat_p(press_last_effect)
    dTtotal = V[0]["temp"] - L[num_effects]["temp_sat_pure"]
    dTicorrector = [1]*num_effects

    total_vapor_produced = L[0]['mass_flow'] - L[num_effects]['mass_flow']

    ############### INICIO DEL PROCESO DE CÁLCULO ###################

    # Estimación inicial de flujos de vapor
    for i in range(1,num_effects+1):
        V[i]["mass_flow"] = total_vapor_produced/num_effects

    finish = False
    current_iter = 1
    while not finish:

        dTneto = dTtotal
        ufactor = 0
        for i in range(1,num_effects+1):
            L[i]["mass_flow"] = L[i-1]["mass_flow"] - V[i]["mass_flow"]
            L[i]["mass_frac"] = L[i-1]["mass_flow"] * L[i-1]["mass_frac"] / L[i]["mass_flow"]

            if current_iter==1:
                L[i]["BPE"] = calcBPE(L[0]["temp"], L[i]['mass_frac'])
                L[i]["U"] = calcU(i, L[0]["temp"], L[i]['mass_frac'])
            else:
                L[i]["BPE"] = calcBPE(L[i]["temp"], L[i]['mass_frac'])
                L[i]["U"] = calcU(i, L[i]["temp"], L[i]['mass_frac'])

            dTneto -= L[i]["BPE"]
            ufactor += 1/L[i]["U"]

        dTi = [0]*(num_effects+1)
        for i in range(1,num_effects+1):
            dTi[i] = dTicorrector[i-1]*dTneto*(1/L[i]["U"])/ufactor
            L[i]["temp"] = V[i-1]["temp"] - dTi[i]
            L[i]["entalpy"] = calcLiquidEntalpy(L[i]["temp"], L[i]['mass_frac'])
            L[i]["entalpy_sat"] = steamTable.hL_t(L[i]["temp"])
            V[i]["temp"] = L[i]["temp"]
            V[i]["entalpy_sat"] = steamTable.hV_t(V[i]["temp"])
            V[i]["entalpy"] = V[i]["entalpy_sat"] + steamTable.CpV_t(V[i]["temp"])*L[i]["BPE"]
            V[i]["lambda"] = V[i]["entalpy_sat"] - L[i]["entalpy_sat"]

        matrix_model = np.zeros((num_effects*2, num_effects*2))
        ind_vector = np.zeros((num_effects*2, 1))

        for k in range(0,num_effects):
            if k==0:
                matrix_model[k,k] = 1
                matrix_model[k,k+num_effects] = 1
                ind_vector[k,0] = L[k]["mass_flow"]
                matrix_model[k+num_effects,k] = L[k+1]["entalpy"]
                matrix_model[k+num_effects,k+num_effects-1] = -V[k]["lambda"]
                matrix_model[k+num_effects,k+num_effects] = V[k+1]["entalpy"]
                ind_vector[k+num_effects,0] = L[k]["entalpy"]*L[k]["mass_flow"]
            elif k<num_effects-1:
                matrix_model[k,k] = 1
                matrix_model[k,k-1] = -1
                matrix_model[k,k+num_effects] = 1
                ind_vector[k,0] = 0
                matrix_model[k+num_effects,k] = L[k+1]["entalpy"]
                matrix_model[k+num_effects,k-1] = -L[k]["entalpy"]
                matrix_model[k+num_effects,k+num_effects-1] = -V[k]["lambda"]
                matrix_model[k+num_effects,k+num_effects] = V[k+1]["entalpy"]
            else:
                matrix_model[k,k-1] = -1
                matrix_model[k,k+num_effects] = 1
                ind_vector[k,0] = -L[k+1]["mass_flow"]
                matrix_model[k+num_effects,k-1] = -L[k]["entalpy"]
                matrix_model[k+num_effects,k+num_effects-1] = -V[k]["lambda"]
                matrix_model[k+num_effects,k+num_effects] = V[k+1]["entalpy"]
                ind_vector[k+num_effects,0] = L[k+1]["entalpy"]*L[k+1]["mass_flow"]

        results = np.linalg.solve(matrix_model, ind_vector)

        vap_error = []
        for k in range(0,num_effects):

            if k<num_effects-1:
                L[k+1]["mass_flow"] = results[k,0]
            else:
                V[0]["mass_flow"] = results[k,0]

            vap_error.append(100 * abs(V[k+1]["mass_flow"] - results[k+num_effects,0]) / V[k+1]["mass_flow"])
            V[k+1]["mass_flow"] = results[k+num_effects,0]

        if np.max(np.array(vap_error))<10.0:

            areas = []

            for k in range(1, num_effects+1):
                areas.append(V[k]["lambda"]*V[k]["mass_flow"]/(L[k]["U"]*dTi[k]))

            areas_error = []
            areas_mean = np.mean(areas)

            for area in areas:
                areas_error.append(100 * abs(area - areas_mean) / areas_mean)

            if np.max(np.array(areas_error))<10.0:
                finish = True
            else:
                dTicorrector = areas/areas_mean
        
        current_iter += 1
        if current_iter > max_iter:
            finish = True
            st.warning("El cálculo no ha convergido tras las iteraciones máximas permitidas.")

    # Resultados
    df1 = pd.DataFrame(V)
    df2 = pd.DataFrame(L)
    df = pd.concat([df2, df1], ignore_index=True)

    st.subheader("Resultados")
    st.dataframe(df)

    areas = [V[k]["lambda"] * V[k]["mass_flow"] / (L[k]["U"] * dTi[k]) for k in range(1, num_effects + 1)]
    areas_mean = np.mean(areas)

    st.write(f"Áreas de los efectos: {areas_mean:.1f} m²")
    steam_economy = total_vapor_produced / V[0]["mass_flow"]
    st.write(f"Economía del vapor: {steam_economy:.2f}")

if st.button('Realizar cálculos'):
    calcular_evaporadores()