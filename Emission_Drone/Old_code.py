# def motor_U_I(M, KV0, Um0, Im0, N, Rm):
#     U = eq.motorVoltage(M, KV0, Um0, Im0, N, Rm)
#     I = eq.motorCurrent( M, KV0, Um0, Im0, Rm)
#     print(U, I)
#     return U, I

# def check_propellers(propeller_matrix, total_mass, labda,dzeta,K0,eta,alpha0, **kwargs):
#     """
#     Calculates whether propeller options can provide enough thrust for T/W = 2
#     propeller_matrix should have a row for each option: ["prop_name", diameter (m), pitch (m), max rpm]
#     """
#     required_thrust = total_mass * 9.81 / 2
#     result_mat = []
#     for row in propeller_matrix:
#         name = row[0]
#         T = propellerThrust(labda,dzeta,K0,eta,row[2],alpha0,0.9*row[3],row[1])
#         if T > required_thrust:
#             result_mat.append([name, "yes"])
#         else:
#             result_mat.append([name, "no"])
#
#     return print(result_mat) # return matrix with [[prop name, yes/no],[prop name, yes/no],[prop name, yes/no]]


# check_propellers(usefullmatrix[1:28,1:5],3.028,labda,dzeta,K0,eta,alpha0)


# test_mat = [["APC 6×4.1SF", 0.1524, 0.10414, 20000],["T-Motor SW 13x5", 0.3302, 0.127, 9600],\
#             ["APC 11x12E", 0.2794, 0.3048, 13636.36364], ["Mezjlik 14x4.5", 0.3356,0.1143, 12900],["T-Motor SW 11x4.2",0.2794,0.10668,11000],\
#             ["T-Motor SW 15x5.6",0.381,0.14224,],["T-Motor CF 14x4.8",0.3556,0.12192,],["APC 10x4.6SF", 0.254, 0.11684, 6500],["APC 11x12WE", 0.254, 0.3048, 15000]]


# def flight_time(battery_w,  battery_cap, frame_w, no_propellers, prop_eff):
#     return ((battery_cap*battery_w)/(frame_w+battery_w))*prop_eff*((frame_w+battery_w)/no_propellers)

# test_mat = [["APC 6×4.1SF", 0.1524, 0.10414, 20000],["T-Motor SW 13x5", 0.3302, 0.127, 9600],\
#             ["APC 11x12E", 0.2794, 0.3048, 13636.36364]]
#
# testmat=np.array([["APC 6×4.1SF", 0.1778, 0.10414, 20000],["APC 11x4.6SF", 0.2794, 0.11684, 15000],["APC 11x12E", 0.2794, 0.3048, 13636.36364],["T-Motor SW 11x4.2",0.2794,0.10668,11000],["DJI Mavic 3", 0.239,0.135,13000]])
