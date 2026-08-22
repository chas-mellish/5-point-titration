import math

def log(x):
    return math.log(x) / math.log(10)

def ten_to(y):
    return math.exp(y * math.log(10))

def perH2CO3(ph, pk1, pk2):
    return 1 / (1 + ten_to(-pk1) / ten_to(-ph) + ten_to(-pk1) * ten_to(-pk2) / ten_to(-2 * ph))

def perHCO3(ph, pk1, pk2):
    return 1 / (ten_to(-ph) / ten_to(-pk1) + 1 + ten_to(-pk2) / ten_to(-ph))

def perCO3(ph, pk2):
    return 1 / (ten_to(-2 * ph) / (ten_to(-pk1) * ten_to(-pk2)) + ten_to(-ph) / ten_to(-pk2) + 1)

def dH2CO3alk(pHf, pHs, pk1, pk2):
    return perHCO3(pHf, pk1, pk2) - perHCO3(pHs, pk1, pk2) + 2 * (perCO3(pHf, pk2) - perCO3(pHs, pk2))

def dHAcalk(pHf, pHs, pKa):
    return perHCO3(pHf, pKa) - perHCO3(pHs, pKa)

def calculate_titration(vx1, vx2, vx3, vx4, ca, pH1, pH2, pH3, pH4, dil):
    # Placeholder for other parameters
    global pK1, pK2, pK11, pK22, pKa, TDS

    # Example values for pK constants (these need to be defined as per your requirements)
    pK1 = 6.35
    pK11 = pK1 + 0.1  # Example adjustment
    pK2 = 10.33
    pK22 = pK2 + 0.1  # Example adjustment
    pKa = 4.76  # Example pKa for acetic acid
    TDS = 3300  # Total Dissolved Solids

    A1 = (vx2 - vx1) * ca * (perH2CO3(pH1, pK11, pK22) - perH2CO3(pH2, pK11, pK22))
    B1 = dH2CO3alk(pH1, pH2, pK1, pK2)
    
    Ct1 = (A1 / B1) / dil * 50000
    
    print(f"H2CO3* alkalinity (Ct1): {Ct1:.2f} mg/l")

# Example usage
vx1 = 1.0  # Volume added for first titration step
vx2 = 2.0  # Volume added for second titration step
vx3 = 3.0  # Volume added for third titration step
vx4 = 4.0  # Volume added for fourth titration step
pH1 = 6.75  # Initial pH after adding Vx1
pH2 = 5.95  # pH after adding Vx2
pH3 = 5.18  # pH after adding Vx3
pH4 = 4.29  # pH after adding Vx4
ca = 0.0728  # Normality of titrant
dil = 1.0  # Dilution factor

# Run titration calculation
calculate_titration(vx1, vx2, vx3, vx4, ca, pH1, pH2, pH3, pH4, dil)