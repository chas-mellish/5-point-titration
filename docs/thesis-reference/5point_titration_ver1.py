import math

def log(x):
    return math.log(x) / math.log(10)

def ten_to(y):
    return math.exp(y * math.log(10))

def input_titration_data():
    print("Enter titration volumes (in mL):")
    values = []
    for i in range(1, 5):
        value = float(input(f"Enter volume for titration {i}: "))
        values.append(value)
    return values

def allocate_data(values):
    global pH0, pH1, pH2, pH3, pH4, vx1, vx2, vx3, vx4, ca, vsundil, vsdil, TDS, ionstr, Nt, Pt
    pH0 = values[0]
    pH1 = values[1]
    pH2 = values[2]
    pH3 = values[3]
    pH4 = values[4]
    vx1 = values[5]
    vx2 = values[6]
    vx3 = values[7]
    vx4 = values[8]
    ca = values[9]
    vsundil = values[10]
    vsdil = values[11]
    TDS = values[12]
    ionstr = values[13]
    Nt = values[14]
    Pt = values[15]

def perH2CO3(ph, pk11, pk22):
    return 1 / (1 + ten_to(-pk11) / ten_to(-ph) + ten_to(-pk11) * ten_to(-pk22) / ten_to(-2 * ph))

def perHCO3(ph, pk11, pk22):
    return 1 / (ten_to(-pk11) + 1 + ten_to(-pk22) / ten_to(-ph))

def perCO3(ph, pk22):
    return 1 / (ten_to(-2 * ph) / (ten_to(-pk22) + ten_to(-pk22) * ten_to(-pk22)) + 1)

def dH2CO3_alk(pHf, pHs, pk11, pk22):
    return perHCO3(pHf, pk11, pk22) - perHCO3(pHs, pk11, pk22) + 2 * (perCO3(pHf, pk22) - perCO3(pHs, pk22))

def dHAcalk(pHf, pHs, pkaa):
    return per(pHf, pkaa) - per(pHs, pkaa)

def calculate_titration(vx_values, ca, pH_values, dil):
    # Placeholder for pK constants (these should be defined as per your requirements)
    pK11 = 6.35
    pK22 = 10.33
    pKa = 4.76  # Example pKa for acetic acid

    # Unpack the titration volumes and pH values
    vx1, vx2, vx3, vx4 = vx_values
    pH1, pH2, pH3, pH4 = pH_values

    A1 = (vx2 - vx1) * ca * (perH2CO3(pH1, pK11, pK22) - perH2CO3(pH2, pK11, pK22))
    B1 = dH2CO3_alk(pH1, pH2, pK11, pK22)
    
    Ct1 = (A1 / B1) / dil * 50000
    
    print(f"H2CO3* alkalinity (Ct1): {Ct1:.2f} mg/l")

def main():
    # Get titration values from user input
    vx_values = input_titration_data()
    
    # Get pH values from user input
    pH_values = []
    for i in range(1, 5):
        pH = float(input(f"Enter pH after adding Vx{i}: "))
        pH_values.append(pH)

    # Get other parameters from user input
    values = []
    for i in range(1, 16):  # Assuming 15 parameters based on allocate_data
        value = float(input(f"Enter value for parameter {i}: "))
        values.append(value)

    # Allocate data based on input values
    allocate_data(values)

    # Run titration calculation
    dil = 1.0  # Dilution factor
    calculate_titration(vx_values, ca, pH_values, dil)

if __name__ == "__main__":
    main()

#c:\Users\charl\OneDrive\My%20Documents\Wastewater\5point_titration_ver2.py