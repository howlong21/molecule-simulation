from functools import reduce
import math
import numpy as np
import matplotlib.pyplot as plt

spath = "./no100dump"
file = open(spath)
arr_multi = np.array([])
clays = []
co2s = []
cations = []
arr_box = []
bound = np.array([])
arr_multi = np.array([])
arr_in = np.array([])
arr_dipole = np.array([])
arr_co2 = np.array([])
arr_mass = np.array([])
clay_ele = ["st", "ao", "mgo", "ob", "obos", "oh", "ohs", "ho"]
co2_ele = ["co", "oc"]


def read_box():
    global file
    box = [0, 0, 0]
    file.readline()
    file.readline()
    file.readline()
    file.readline()
    line_1 = file.readline()
    dd_1 = line_1.split()
    box[0] = float(dd_1[1]) - float(dd_1[0])
    line_1 = file.readline()
    dd_1 = line_1.split()
    box[1] = float(dd_1[1]) - float(dd_1[0])
    line_1 = file.readline()
    dd_1 = line_1.split()
    box[2] = float(dd_1[1]) - float(dd_1[0])
    file.readline()
    return [box[0], box[1], box[2]]


def read_xyz(water, cation, clay, co2):
    global file
    global clay_ele
    global co2_ele
    num_1 = 0
    num_co2 = 0
    water_o = [1, 2, 3]
    co2_o = [1, 2, 3]
    while 1:
        line = file.readline()
        if line:
            linesep = line.split()
            if linesep[0] != "ITEM:":
                if linesep[0] == "ow" or linesep[0] == "hw" or linesep[0] == "ehw" or linesep[0] == "eow":
                    water_o[num_1 % 3] = [float(linesep[1]), float(linesep[2]), float(linesep[3])]
                    num_1 += 1
                    if num_1 % 3 == 0:
                        water.append([water_o[0], water_o[1], water_o[2]])
                elif linesep[0] in clay_ele:
                    clay_atom = [linesep[0], float(linesep[1]), float(linesep[2]), float(linesep[3])]
                    clay.append(clay_atom)
                elif linesep[0] in co2_ele:
                    co2_o[num_co2 % 3] = [float(linesep[1]), float(linesep[2]), float(linesep[3])]
                    num_co2 += 1
                    if num_co2 % 3 == 0:
                        co2.append([co2_o[0], co2_o[1], co2_o[2]])
                elif linesep[0] == "na" or linesep[0] == "ca" or linesep[0] == "cl":
                    clay_cation = [linesep[0], float(linesep[1]), float(linesep[2]), float(linesep[3])]
                    cation.append(clay_cation)
            else:
                return line
        else:
            return line


def readcoor():  # read the coordinate from the dump file,then save water in the array arr_multi,
    # save cations in list cations,save clay_atoms in list clays,
    # the information of box size in array arr_box
    global arr_box
    global arr_multi
    global arr_co2
    global clays
    global cations
    global co2s
    global file
    ymax = 74
    multi = []
    boxs = []
    water = []
    clay = []
    cation = []
    co2 = []
    file.readline()
    while 1:
        box1 = read_box()
        line = read_xyz(water, cation, clay, co2)
        arr_wat = np.array(water, dtype=np.float32)
        arr_cla = np.array(np.array(clay)[:, 3], dtype=np.float32)
        water_min = np.min(arr_wat[:, :, 2])
        clay_min = np.min(arr_cla).item()
        if water_min < clay_min:
            for waat in range(len(water)):
                for ohh in range(3):
                    if water[waat][ohh][2] < clay_min + 3:
                        water[waat][ohh][2] += box1[2]
            for lac1 in range(len(co2)):
                for lac2 in range(len(co2[0])):
                    if co2[lac1][lac2][2] < clay_min + 3:
                        co2[lac1][lac2][2] += box1[2]
            for caion in range(len(cation)):
                if cation[caion][3] < clay_min + 3:
                    cation[caion][3] += box1[2]
            for clas in range(len(clay)):
                if clay[clas][3] < clay_min + 3:
                    clay[clas][3] += box1[2]
        arr_wat = np.array(water, dtype=np.float32)
        arr_cla = np.array(np.array(clay)[:, 3], dtype=np.float32)
        water_max = np.max(arr_wat[:, :, 2])
        clay_max = np.max(arr_cla).item()
        if water_max > clay_max - 2:
            for clas in range(len(clay)):
                if clay[clas][3] < clay_min + 3:
                    # print("before", clay[clas][2])
                    clay[clas][3] += box1[2]
                    # print("after:", clas, clay[clas][2])
            for clas in range(len(water)):
                for wid in range(3):
                    if water[clas][wid][2] < clay_min + 3:
                        water[clas][wid][2] += box1[2]
            for clas in range(len(co2)):
                for wid in range(3):
                    if co2[clas][wid][2] < clay_min + 3:
                        co2[clas][wid][2] += box1[2]
            for clas in range(len(cation)):
                if cation[clas][3] < clay_min + 3:
                    cation[clas][3] += box1[2]
        for ix in range(len(water)):
            for iy in range(3):
                if water[ix][iy][1] > ymax:
                    water[ix][iy][1] -= box1[1]
        for ix in range(len(co2)):
            for iy in range(3):
                if co2[ix][iy][1] > ymax:
                    co2[ix][iy][1] -= box1[1]
        for ix in range(len(cation)):
            if cation[ix][2] > ymax:
                cation[ix][2] -= box1[1]
        for ix in range(len(clay)):
            if clay[ix][2] > ymax:
                clay[ix][2] -= box1[1]
        boxs.append(box1)
        multi.append(water)
        clays.append(clay)
        cations.append(cation)
        co2s.append(co2)
        if line:
            water = []
            clay = []
            cation = []
            co2 = []
        else:
            arr_multi = np.array(multi)
            arr_box = np.array(boxs)
            arr_co2 = np.array(co2s)
            break
    file.close()
    print("all coordinates has been read...")


def writecor():
    global arr_multi
    global clays
    global cations
    global co2s
    file2 = open("./xx", "w")
    num = len(clays[0])+3*arr_multi.shape[1]+3*len(co2s[0])+len(cations[0])
    print(num, file=file2)
    print(num, file=file2)
    for i_a in range(arr_multi.shape[1]):
        print("%s %.5f %.5f %.5f" % ("ow", arr_multi[0, i_a, 0, 0], arr_multi[0, i_a, 0, 1],
                                     arr_multi[0, i_a, 0, 2]), sep=" ", file=file2)
        print("%s %.5f %.5f %.5f" % ("hw", arr_multi[0, i_a, 1, 0], arr_multi[0, i_a, 1, 1],
                                     arr_multi[0, i_a, 1, 2]), sep=" ", file=file2)
        print("%s %.5f %.5f %.5f" % ("hw", arr_multi[0, i_a, 2, 0], arr_multi[0, i_a, 2, 1],
                                     arr_multi[0, i_a, 2, 2]), sep=" ", file=file2)
    for ix in range(len(clays[0])):
        print("%s %.5f %.5f %.5f" % (clays[0][ix][0], clays[0][ix][1], clays[0][ix][2], clays[0][ix][3]),
              sep=" ", file=file2)
    for ix in range(len(cations[0])):
        print("%s %.5f %.5f %.5f" % (cations[0][ix][0], cations[0][ix][1], cations[0][ix][2], cations[0][ix][3]),
              sep=" ", file=file2)
    for ix in range(len(co2s[0])):
        print("%s %.5f %.5f %.5f" % ("oc", co2s[0][ix][0][0], co2s[0][ix][0][1], co2s[0][ix][0][2]),
              sep=" ", file=file2)
        print("%s %.5f %.5f %.5f" % ("co", co2s[0][ix][1][0], co2s[0][ix][1][1], co2s[0][ix][1][2]),
              sep=" ", file=file2)
        print("%s %.5f %.5f %.5f" % ("oc", co2s[0][ix][2][0], co2s[0][ix][2][1], co2s[0][ix][2][2]),
              sep=" ", file=file2)
    file2.close()


def calmass(thewater):  # calculate the center of mass on one water
    # return [(thewater[0][0] * 15.9994 + thewater[1][0] * 1.008 + thewater[2][0] * 1.008) / 18.0154,
    #         (thewater[0][1] * 15.9994 + thewater[1][1] * 1.008 + thewater[2][1] * 1.008) / 18.0154,
    #         (thewater[0][2] * 15.9994 + thewater[1][2] * 1.008 + thewater[2][2] * 1.008) / 18.0154]
    return [thewater[0][0], thewater[0][1], thewater[0][2]]


def allmass():  # calculate center of mass on all waters
    global arr_multi
    global arr_mass
    mass_temp = []
    for tt in range(arr_multi.shape[0]):
        mass1 = []
        for watind in range(arr_multi.shape[1]):
            mass1.append(calmass(arr_multi[tt][watind]))
        mass_temp.append(mass1)
    arr_mass = np.array(mass_temp)


def xshift(watt1, watt2):
    return (watt1[0] - watt2[0])*(watt1[0] - watt2[0])


def msdxy(whe_contin, max_time, re_left, re_right, arr_inclay):
    global arr_mass
    global arr_multi
    fin_dip_up = []
    fin_dip_down = []

    for time in range(0, min([int(arr_multi.shape[0] / 2) + 1, max_time])):
        fin_dip_1 = []
        fin_dip_2 = []
        for numb in range(arr_multi.shape[0] - time):
            diplist_1 = []
            diplist_2 = []
            if whe_contin:
                for wat_ind in arr_inclay:
                    whe_in = 1
                    for intime in range(numb, numb + time + 1):
                        if (arr_mass[intime][wat_ind][1] > re_right) or \
                                (arr_mass[intime][wat_ind][1] < re_left):
                            whe_in = 0
                            break
                    if whe_in:
                        diplist_1.append(xshift(watt1=arr_mass[numb][wat_ind], watt2=arr_mass[numb + time][wat_ind]))
                        diplist_2.append(1)
            else:
                for wat_ind in arr_inclay:
                    # print("every_bound", every_bound, "arr_mass", arr_mass[numb][wat_ind][1])
                    if re_left < arr_mass[numb][wat_ind][1] < re_right and \
                                        re_left < arr_mass[numb + time][wat_ind][1] < re_right:
                        diplist_1.append(xshift(watt1=arr_mass[numb][wat_ind], watt2=arr_mass[numb + time][wat_ind]))
                        diplist_2.append(1)
            if diplist_1:
                fin_dip_1.append(reduce(lambda x, y: x + y, diplist_1))
            if diplist_2:
                fin_dip_2.append(reduce(lambda x, y: x + y, diplist_2))
        if fin_dip_1:
            fin_dip_up.append(reduce(lambda x, y: x + y, fin_dip_1))
        else:
            fin_dip_up.append(0)
        if fin_dip_2:
            fin_dip_down.append(reduce(lambda x, y: x + y, fin_dip_2))
        else:
            fin_dip_down.append(0)
    return fin_dip_up, fin_dip_down


def msdxy_all(whe_con, max_time, re_left, re_right):  # calculate the corelation
    # of msd in xy direction of multi layers via invoking function msdxy
    global bound
    global arr_multi
    allco2 = list(range(arr_multi.shape[1]))
    output = msdxy(whe_con, max_time, re_left, re_right, allco2)
    return output[0], output[1]


def cal_msdfun(npara, z_left, z_right):
    global arr_multi
    mid_fun = []
    for timeind in range(arr_multi.shape[0]):
        recwater = []
        for waternum in range(arr_multi.shape[1]):
            if z_left < arr_multi[timeind][waternum][0][1] < z_right:
                recwater.append(math.sin(npara * math.pi * (arr_multi[timeind][waternum][0][1] -
                                                            z_left) / (z_right - z_left)))
            else:
                recwater.append(0)
        mid_fun.append(recwater)
    return np.array(mid_fun)


def midmsd(whe_contin, max_time, re_left, re_right, mid_fun, arr_inclay):
    global arr_multi
    fin_dip_up = []
    fin_dip_down = []
    for time in range(0, min([int(arr_multi.shape[0] / 2) + 1, max_time])):
        fin_dip_1 = []
        fin_dip_2 = []
        for numb in range(arr_multi.shape[0] - time):
            diplist_1 = []
            diplist_2 = []
            if whe_contin:
                for wat_ind in arr_inclay:
                    whe_in = 1
                    if re_left < arr_multi[numb][wat_ind][0][1] < re_right:
                        diplist_2.append(1)
                    else:
                        continue
                    for intime in range(numb + 1, numb + time + 1):
                        if (arr_multi[intime][wat_ind][0][1] > re_right) or\
                                (arr_multi[intime][wat_ind][0][1] < re_left):
                            whe_in = 0
                            break
                    if whe_in:
                        diplist_1.append(mid_fun[numb][wat_ind] *
                                         mid_fun[numb + time][wat_ind])
            else:
                for wat_ind in arr_inclay:
                    if re_left < arr_multi[numb][wat_ind][0][1] < re_right:
                        diplist_2.append(1)
                    else:
                        continue
                    if re_left < arr_multi[numb + time][wat_ind][0][1] < re_right:
                        diplist_1.append(mid_fun[numb][wat_ind] *
                                         mid_fun[numb + time][wat_ind])
            if diplist_1:
                fin_dip_1.append(reduce(lambda x, y: x + y, diplist_1))
            if diplist_2:
                fin_dip_2.append(reduce(lambda x, y: x + y, diplist_2))
        if fin_dip_1:
            fin_dip_up.append(reduce(lambda x, y: x + y, fin_dip_1))
        else:
            fin_dip_up.append(0)
        if fin_dip_2:
            fin_dip_down.append(reduce(lambda x, y: x + y, fin_dip_2))
        else:
            fin_dip_down.append(0)
    return fin_dip_up, fin_dip_down


def msdmid_all(npara, whe_con, max_time, re_left, re_right):  # calculate the corelation
    # of msd in z direction of multi layers via invoking function msdmid
    #  be aware of that the density of the vater must be a constant
    mid_fun = cal_msdfun(npara, re_left, re_right)
    allmulti = list(range(arr_multi.shape[1]))
    output = midmsd(whe_con, max_time, re_left, re_right, mid_fun, allmulti)
    return output[0], output[1]


def prepare():
    readcoor()  # read the coordinate from the dump file,then save water in the array arr_multi,
    # save cations in list cations,save clay_atoms in list clays,
    # the information of box size in array arr_box
    writecor()  # write the coordinate read to a file to confirm whether the lists are true or not
    # cal_dipole()  # calculate the dipole on each water
    allmass()


def msdout(whe_con, max_time, re_left, re_right):
    namelat = "_" + str(whe_con) + "_" + str(re_left) + "_" + str(re_right)
    fname = open("msd" + namelat, "w")
    fgh = msdxy_all(whe_con, max_time, re_left, re_right)
    for io in range(len(fgh[0])):
        print("%d %.5f %.5f " % (io, fgh[0][io], fgh[1][io]), sep=" ", file=fname)
    fname.close()


def midmsdout(npara, whe_con, max_time, re_left, re_right):
    namelat = "_" + str(whe_con) + "_" + str(re_left) + "_" + str(re_right)
    fname = open("midmsd" + namelat, "w")
    fgh = msdmid_all(npara, whe_con, max_time, re_left, re_right)
    for io in range(len(fgh[0])):
        print("%d %.5f %.5f " % (io, fgh[0][io], fgh[1][io]), sep=" ", file=fname)
    fname.close()

readcoor()
writecor()
# tmptmp = arr_co2
# arr_co2 = arr_multi
# arr_multi = tmptmp
allmass()
# midmsdout(1, 0, 500, 15, 50)
# prepare()
fgh = msdxy_all(0, 250, 65, 70)
# fgh = msdmid_all(1, 0, 250, 65, 70)
listx = []
listy = []

for i in range(len(fgh[0])):
    listx.append((i+1)*0.2)
    listy.append(fgh[0][i]/fgh[1][i])
    # listy.append(math.log2(fgh[0][i]/fgh[1][i]))
plt.plot(listx, listy)
plt.show()
