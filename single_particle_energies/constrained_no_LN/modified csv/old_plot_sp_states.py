import matplotlib.pyplot as plt
import numpy as np
# z126_n184
# proton hole
#   0       1       2       3       4       5       6
#   1h_9/2  2f_7/2  2f_5/2  3p_3/2  3p_1/2  1i_13/2
# proton particle
#   2g_9/2  3d_5/2  1i_11/2 2g_7/2  4s_1/2  3d_3/2  1j_15/2

# neutron hole
#   0       1       2       3       4       5       6
# 2g_9/2    3d_5/2  1i_11/2 2g_7/2  4s_1/2  3d_3/2  1j_15/2
# neutron particle
# 2h_11/2   1j_13/2 3f_7/2  2h_9/2  4p_3/2  3f_5/2  1k_17/2

# *edf*[] and *edf*_sp[] has 4 elements, in the order of proton hole,particle ; neutron hole,particle

# UNEDF0
une0_gs = -2143.176664; une0 = []; une0_sp = []
une0.append(np.array([-2139.442616,-2140.749218,-2142.735031,-2142.762948,
-2144.240317,-2140.617194]))
une0.append(np.array([-2139.798686,-2136.677350,-2140.371333,-2137.263897,
-2135.149205,-2135.552172,-2139.986715]))
une0.append(np.array([-2132.130595,-2134.500307,-2131.849284,-2134.149479,
-2134.165359,-2135.464007,-2132.944024]))
une0.append(np.array([-2148.462664,-2148.436894,-2146.554833,-2146.623799,
-2145.633081,-2145.575093,-2148.043498]))
une0_sp.append(une0_gs-une0[0])
une0_sp.append(une0[1]-une0_gs)
une0_sp.append(une0_gs-une0[2])
une0_sp.append(une0[3]-une0_gs)

# UNEDF1
une1_gs = -2139.088146; une1 = []; une1_sp = []
une1.append(np.array([-2134.139618,-2136.692253,-2138.218004,-2139.309369,
-2139.872184,-2136.612877]))
une1.append(np.array([-2135.052265,-2131.866758,-2137.048595,-2133.104192,
-2130.361171,-2131.041439,-2135.237678]))
une1.append(np.array([-2127.843748,-2127.878587,-2127.189077,-2130.188906,
-2131.538420,-2131.672476,-2128.438509]))
une1.append(np.array([-2144.305687,-2143.886032,-2141.910187,-2141.801136,
-2140.746490,-2140.713737,-2143.857272]))
une1_sp.append(une1_gs-une1[0])
une1_sp.append(une1[1]-une1_gs)
une1_sp.append(une1_gs-une1[2])
une1_sp.append(une1[3]-une1_gs)

# UNEDF2
une2_gs = -2137.378883; une2 = []; une2_sp = []
une2.append(np.array([-2132.391579,-2134.992666,-2136.474817,-2137.705276
,-2138.20299,-2135.230836]))
une2.append(np.array([-2133.165895,-2129.833451,-2135.230836,-2131.225746
,-2128.170269,-2128.977881,-2133.705968]))
une2.append(np.array([-2126.153528,-2129.049387,-2125.534738,-2128.338592
,-2129.954793,-2130.024453,-2126.149876]))
une2.append(np.array([-2142.516121,-2142.211823,-2139.658127,-2139.875110
,-2138.617465,-2138.599913,-2142.607373]))
une2_sp.append(une2_gs-une2[0])
une2_sp.append(une2[1]-une2_gs)
une2_sp.append(une2_gs-une2[2])
une2_sp.append(une2[3]-une2_gs)

# SV-MIN
sv_min_gs = -2132.975829; sv_min = []; sv_min_sp = []
sv_min.append(np.array([-2128.433629,-2130.556336,-2132.282819,-2133.275787
,-2133.996848,-2130.616578]))
sv_min.append(np.array([-2129.027849,-2125.834923,-2130.490668,-2126.836260
,-2124.357874,-2124.867621,-2129.043923]))
sv_min.append(np.array([-2121.309098,-2121.342538,-2120.689783,-2123.813823
,-2125.214466,-2125.324267,-2122.145721]))
sv_min.append(np.array([-2138.020069,-2137.992197,-2135.581668,-2135.474666
,-2134.337252,-2134.303517,-2137.479031]))
sv_min_sp.append(sv_min_gs-sv_min[0])
sv_min_sp.append(sv_min[1]-sv_min_gs)
sv_min_sp.append(sv_min_gs-sv_min[2])
sv_min_sp.append(sv_min[3]-sv_min_gs)

# SLY4
sly4_gs =  -2129.59854; sly4 = []; sly4_sp = []
sly4.append(np.array([-2123.595088,-2126.113250,-2128.225801,-2129.321279
,-2129.315758,-2126.580546]))
sly4.append(np.array([-2125.656835,-2121.995711,-2127.259929,-2123.020937
,-2120.396112,-2120.789699,-2125.355418]))
sly4.append(np.array([-2117.320296,-2120.596606,-2116.204474,-2120.454623
,-2122.001715,-2122.124458,-2118.628565]))
sly4.append(np.array([-2134.159963,-2134.126655,-2131.351752,-2130.985805
,-2129.876787,-2129.836043,-2133.142787]))
sly4_sp.append(sly4_gs-sly4[0])
sly4_sp.append(sly4[1]-sly4_gs)
sly4_sp.append(sly4_gs-sly4[2])
sly4_sp.append(sly4[3]-sly4_gs)



# proton
#   0       1       2       3       4       5       6
#   1h_9/2  2f_7/2  2f_5/2  3p_3/2  3p_1/2  1i_13/2
# proton particle
#   2g_9/2  3d_5/2  1i_11/2 2g_7/2  4s_1/2  3d_3/2  1j_15/2

# delta_j = 3 states index are:
# proton::
# hole [0] - [3]; [1] - [4]; [5] - [1];  # 3 pairs

# particle [6] - [0]; [0] - [5]; [3] - [4]; [2] - [1]  # 4 pairs

# neutron::
# hole [6] - [0]; [0] - [5]; [3] - [4]; [2] - [1]  # 4 pairs

# particle  [6] - [0]; [0] - [5]; [1] - [2]; [3] - [4]; # 4 pairs

# proton hole
#   0       1       2       3       4       5       6
#   1h_9/2  2f_7/2  2f_5/2  3p_3/2  3p_1/2  1i_13/2
# proton particle
#   0       1       2       3       4       5       6
#   2g_9/2  3d_5/2  1i_11/2 2g_7/2  4s_1/2  3d_3/2  1j_15/2

#une0_sp = une1_sp
plt.title("z126_n184 neutron")
plt.suptitle("UNEDF0")
plt.ylabel(" s.p. levels (MeV)")
# delta_j pairs:
xx = np.array([0,0.5])
# proton
## hole
#plt.plot(xx,np.array([une0_sp[0][0],une0_sp[0][0]]),"r-")
#plt.plot(xx,np.array([une0_sp[0][3],une0_sp[0][3]]),"r-")
#plt.plot(xx+1,np.array([une0_sp[0][1],une0_sp[0][1]]),"b-")
#plt.plot(xx+1,np.array([une0_sp[0][4],une0_sp[0][4]]),"b-")
#plt.plot(xx+2,np.array([une0_sp[0][5],une0_sp[0][5]]),"g-")
#plt.plot(xx+2,np.array([une0_sp[0][1],une0_sp[0][1]]),"g-")
## particle
#plt.plot(xx+3,np.array([une0_sp[1][6],une0_sp[1][6]]),"r-")
#plt.plot(xx+3,np.array([une0_sp[1][0],une0_sp[1][0]]),"r-")
#plt.plot(xx+4,np.array([une0_sp[1][0],une0_sp[1][0]]),"b-")
#plt.plot(xx+4,np.array([une0_sp[1][5],une0_sp[1][5]]),"b-")
#plt.plot(xx+5,np.array([une0_sp[1][3],une0_sp[1][3]]),"g-")
#plt.plot(xx+5,np.array([une0_sp[1][4],une0_sp[1][4]]),"g-")
#plt.plot(xx+6,np.array([une0_sp[1][2],une0_sp[1][2]]),"y-")
#plt.plot(xx+6,np.array([une0_sp[1][1],une0_sp[1][1]]),"y-")




# hole
plt.plot(xx,np.array([une0_sp[2][6],une0_sp[2][6]]),"r-")
plt.plot(xx,np.array([une0_sp[2][0],une0_sp[2][0]]),"r-")
plt.plot(xx+1,np.array([une0_sp[2][0],une0_sp[2][0]]),"b-")
plt.plot(xx+1,np.array([une0_sp[2][5],une0_sp[2][5]]),"b-")
plt.plot(xx+2,np.array([une0_sp[2][3],une0_sp[2][3]]),"g-")
plt.plot(xx+2,np.array([une0_sp[2][4],une0_sp[2][4]]),"g-")
plt.plot(xx+3,np.array([une0_sp[2][2],une0_sp[2][2]]),"y-")
plt.plot(xx+3,np.array([une0_sp[2][1],une0_sp[2][1]]),"y-")

# particle
plt.plot(xx+4,np.array([une0_sp[3][6],une0_sp[3][6]]),"r-")
plt.plot(xx+4,np.array([une0_sp[3][0],une0_sp[3][0]]),"r-")
plt.plot(xx+5,np.array([une0_sp[3][0],une0_sp[3][0]]),"b-")
plt.plot(xx+5,np.array([une0_sp[3][5],une0_sp[3][5]]),"b-")
plt.plot(xx+6,np.array([une0_sp[3][1],une0_sp[3][1]]),"g-")
plt.plot(xx+6,np.array([une0_sp[3][2],une0_sp[3][2]]),"g-")
plt.plot(xx+7,np.array([une0_sp[3][3],une0_sp[3][3]]),"y-")
plt.plot(xx+7,np.array([une0_sp[3][4],une0_sp[3][4]]),"y-")


# neutron
# neutron hole
#   0       1       2       3       4       5       6
# 2g_9/2    3d_5/2  1i_11/2 2g_7/2  4s_1/2  3d_3/2  1j_15/2
# neutron particle
#   0       1       2       3       4       5       6
# 2h_11/2   1j_13/2 3f_7/2  2h_9/2  4p_3/2  3f_5/2  1k_17/2
# neutron::
# hole [6] - [0]; [0] - [5]; [3] - [4]; [2] - [1]  # 4 pairs

# particle  [6] - [0]; [0] - [5]; [1] - [2]; [3] - [4]; # 4 pairs

# proton ticks
plt.xticks([0.25,1.25,2.25,3.25,4.25,5.25,6.25],
["1h$_{9/2}$\n3p$_{3/2}$","2f$_{7/2}$\n3p$_{1/2}$","1i$_{13/2}$\n2f$_{7/2}$",
"1k$_{17/2}$\n2g$_{9/2}$","2g$_{9/2}$\n3d$_{3/2}$","2g$_{7/2}$\n4s$_{1/2}$","1i$_{11/2}$\n3d$_{5/2}$",])
# neutron ticks
plt.xticks([0.25,1.25,2.25,3.25,4.25,5.25,6.25,7.25],
["$_{}$\n$_{}$","$_{}$\n$_{}$","$_{}$\n$_{}$","$_{}$\n$_{}$","$_{}$\n$_{}$","$_{}$\n$_{}$","$_{}$\n$_{}$","$_{}$\n$_{}$"])
#"1k$_{17/2}$\n2g$_{9/2}$","2g$_{9/2}$\n3d$_{3/2}$","2g$_{7/2}$\n4s$_{1/2}$","1i$_{11/2}$\n3d$_{5/2}$"
#])

plt.show()
