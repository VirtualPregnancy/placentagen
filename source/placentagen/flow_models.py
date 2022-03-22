import numpy as np
from scipy import special


def bisection_method_diam(a, b, fit_passive_params, fit_myo_params, fit_flow_params, fixed_flow_params, pressure,verbose):


    if (tension_balance(fit_passive_params, fit_myo_params, fit_flow_params, fixed_flow_params, a, pressure) * tension_balance(fit_passive_params,fit_myo_params, fit_flow_params, fixed_flow_params, b, pressure)) > 0:
        #if verbose:
        #    print ('Bisection interval will not work')
        #    print(a, tension_balance(fit_passive_params, fit_myo_params, fit_flow_params,fixed_flow_params, a, pressure), pressure)
        #    print(b, tension_balance(fit_passive_params, fit_myo_params, fit_flow_params, fixed_flow_params, b, pressure), pressure)
        diam = 0
    else:
        a1 = a
        b1 = b
        eps = 1.0e-12
        err = 1.0
        while err >= eps:
            x0 = (a1 + b1) / 2.
            if (tension_balance(fit_passive_params, fit_myo_params, fit_flow_params, fixed_flow_params, a1,
                                pressure) * tension_balance(fit_passive_params, fit_myo_params, fit_flow_params, fixed_flow_params,x0,
                                                            pressure)) < 0:
                b1 = x0
                err = abs(a1 - x0)
            else:
                a1 = x0
                err = abs(b1 - x0)
        diam = x0

    return diam


def calc_plug_resistance(mu,Dp,porosity,radius, length):
    #Calculates resistance of a completely plugged segment of tube (artery)
    K = (Dp ** 2. / 180.) * ((porosity ** 3.) / (1. - porosity) ** 2.) #permeability
    gamma = 1. / (1. + 2.5 * (1. - porosity))
    area = np.pi * radius ** 2.
    resistance = (mu * length) / ((K * area) * (1. - 2. * (np.sqrt(K/gamma))/radius* special.iv(1,np.sqrt(gamma/K) * radius) / (special.iv(0, np.sqrt( gamma/K) * radius))))
    return resistance

def calc_plug_shear(mu,Dp,porosity,radius, length,flow,resistance,r_val):
    # Calculates shear stress a completely plugged segment of tube (artery)
    K = (Dp ** 2. / 180.) * ((porosity ** 3.) / (1. - porosity) ** 2.)  # permeability
    gamma = 1. / (1. + 2.5 * (1. - porosity))
    shear = np.zeros(len(r_val))
    for j in range(0,len(r_val)):
        shear[j] = (np.sqrt(gamma)*flow*resistance/length)* (special.iv(1,np.sqrt(gamma/K) * r_val[j])/special.iv(0,np.sqrt(gamma/K) * radius))

    shear_at_wall = (np.sqrt(gamma)*flow*resistance/length)* (special.iv(1,np.sqrt(gamma/K) * radius)/special.iv(0,np.sqrt(gamma/K) * radius))
    return shear,shear_at_wall

def calc_tube_resistance(mu,radius,length):
    #Calculates resistance of a simple Poisseuille tube
    resistance = (8.* mu * length) / (np.pi * radius**4.);
    return resistance

def calc_tube_shear(mu, radius,flow):
    # Calculates shear in a simple Poisseuille tube
    shear = 4.* mu * flow/(np.pi * radius**3.)
    return shear


def calc_funnel_resistance(mu, radius_a, radius_b, length_a,length_b):
    #Calculates resistance of a funnel shaped segment of tube
    #Radius A is the smaller inlet radius, radius b is the larger outlet radius
    #length_a is distance down the main axis of the vessel that the funnel starts
    #length_b is distance down the main axis of the vessel that the funnel ends

    rate_increase_c = (radius_b -radius_a)/(length_b-length_a)

    resistance =(8.* mu) / (np.pi * radius_a**4.)*(radius_a/(3.*rate_increase_c)-radius_a**4./(3.*rate_increase_c*(radius_a+rate_increase_c*(length_b-length_a))**3.))

    return resistance

def calc_total_tension(fit_passive_params, fit_myo_params, fit_flow_params, fixed_flow_params,diameter, pressure):

    include_passive = True
    include_flow = True
    include_myo = True
    if np.array(fit_myo_params).all() == 0:
        include_myo = False
    if np.array(fit_flow_params).all() == 0:
        include_flow = False

    #Defining parameters of importance in passive model
    D0 = fit_passive_params[0]
    Cpass = fit_passive_params[1]
    Cpassdash = fit_passive_params[2]
    if include_myo:
        #Defining parameters relavent to myogenic model
        Cact = fit_myo_params[0]
        Cactdash = fit_myo_params[1]
        Cactdashdash = fit_myo_params[2]
        Cmyo = fit_myo_params[3]
        Cdashdashtone = fit_myo_params[4]
    if include_flow:
        #Defining parameters related to blood flow through vessel
        Cshear = fit_flow_params[0]
        Cshear2 = -fit_flow_params[1]
        tau1 = fit_flow_params[2]
        tau2 = fit_flow_params[3]
        mu = fixed_flow_params[0]
        length = fixed_flow_params[1]
        dp_blood = fixed_flow_params[2]
        resistance = calc_tube_resistance(mu,diameter/2.,length)
        flow = dp_blood/resistance
        tau = calc_tube_shear(mu,diameter/2.,flow)

    if not include_myo and not include_flow:
        Tmaxact=0.
        total_tension = Cpass * np.exp(Cpassdash * (diameter / D0 - 1.))
    else:
        if not include_flow:
            Stone = Cmyo * pressure * diameter / 2. + Cdashdashtone
        else:
            if tau < tau1:
                Stone = Cmyo * pressure * diameter / 2. + Cdashdashtone
            elif tau < tau2:
                Stone = Cmyo * pressure * diameter / 2. + Cdashdashtone + Cshear * (tau - tau1)
            else:
                Stone = Cmyo * pressure * diameter / 2. + Cdashdashtone + Cshear2 * (
                            tau - (Cshear / Cshear2 * (tau1 - tau2) + tau2))


        A = 1. / (1 + np.exp(-Stone))

        Tmaxact = Cact * np.exp(-((diameter / D0 - Cactdash) / Cactdashdash) ** 2.)
        total_tension = Cpass * np.exp(Cpassdash * (diameter / D0 - 1.))  + A * Tmaxact

    return Tmaxact,total_tension

def diameter_from_pressure(fit_passive_params,fit_myo_params,fit_flow_params,fixed_flow_params, pressure,verbose):
        #Dp_blood is driving pressure (mmHg)
        #pressure is transmural pressure (kPa)
        # calculates a passive diameter under zero flow conditions
        include_passive = True
        include_flow = True
        include_myo = True
        if np.array(fit_myo_params).all() == 0:
            include_myo = False
        if np.array(fit_flow_params).all() == 0:
            include_flow = False

        if not include_myo and not include_flow:
            diameter= bisection_method_diam(5.,500.,fit_passive_params,fit_myo_params,fit_flow_params,fixed_flow_params,pressure,verbose)
        else:
            diameter_pass = bisection_method_diam(5.,500.,fit_passive_params,[0.,0.],[0.,0.],[0.,0.],pressure,verbose)
            D0 = fit_passive_params[0]
            dp_blood = fixed_flow_params[2]
            reference_diameter = np.max([diameter_pass,D0])
            if(pressure>12.) or (dp_blood>0):
                #Looks for a diameter between 5 um and 10% larger than the passive diameter at that transmural pressure as a possible root
                lowest_sign = find_possible_roots(10.,reference_diameter*1.10, fit_passive_params, fit_myo_params,fit_flow_params,fixed_flow_params, pressure,verbose)
            else:
                lowest_sign=np.array([10.,reference_diameter*1.10])

            diameter = bisection_method_diam(lowest_sign[0], lowest_sign[1], fit_passive_params,fit_myo_params,fit_flow_params,fixed_flow_params, \
                                                                                                     pressure,verbose)
            #if verbose:
            if diameter<10:
                print('zero',pressure,dp_blood,lowest_sign,diameter_pass,reference_diameter)
                lowest_sign = find_possible_roots(0.,reference_diameter*1.10, fit_passive_params, fit_myo_params,fit_flow_params,fixed_flow_params, pressure,verbose)
                diameter = bisection_method_diam(lowest_sign[0], lowest_sign[1], fit_passive_params,fit_myo_params,fit_flow_params,fixed_flow_params, \
                                                                                                     pressure,verbose)
                print(diameter,lowest_sign)

        return diameter

def find_possible_roots(low_diam,high_diam, fit_passive_params, fit_myo_params, fit_flow_params,fixed_flow_params, pressure,verbose):

    discretise = 500
    diameter_range = np.linspace(high_diam,low_diam,discretise) #finding root closest to the passive diameter
    ten_resid = np.zeros(len(diameter_range)) #tension residual
    lowest_signchange = np.zeros(2)
    i = 0
    ten_resid[i] = tension_balance(fit_passive_params,fit_myo_params, fit_flow_params,fixed_flow_params,diameter_range[i], pressure)
    samesign = True
    #for every diameter in the range calculate the tension residual
    while samesign:
        i=i+1
        ten_resid[i] = tension_balance(fit_passive_params,fit_myo_params, fit_flow_params,fixed_flow_params,diameter_range[i],pressure)
        if np.sign(ten_resid[i]) != np.sign(ten_resid[i-1]):
            samesign = False
        if(i==(discretise -1)) and samesign:
            samesign = False


    if i==(discretise-1) and samesign: #should be searching for a higher diameter
        diameter_range = np.linspace(high_diam,high_diam + (high_diam-low_diam), discretise)
        ten_resid = np.zeros(len(diameter_range))
        i = 0
        ten_resid[i] = tension_balance(fit_passive_params,fit_myo_params, fit_flow_params,fixed_flow_params,diameter_range[i], pressure)
        samesign = True
        while samesign:
            i=i+1
            ten_resid[i] = tension_balance(fit_passive_params,fit_myo_params, fit_flow_params,fixed_flow_params,diameter_range[i],pressure)
            if np.sign(ten_resid[i]) != np.sign(ten_resid[i-1]):
                samesign = False
        if verbose:
            print('had to go again')
    lowest_signchange = np.zeros(2)
    lowest_signchange[0] = diameter_range[i]
    lowest_signchange[1] = diameter_range[i-1]

    return lowest_signchange

def tension_balance(fit_passive_params,fit_myo_params, fit_flow_params,fixed_flow_params, diameter, pressure):
    [Tmaxact,total_tension] = calc_total_tension(fit_passive_params,fit_myo_params, fit_flow_params,fixed_flow_params, diameter, pressure)
    f = total_tension - pressure * diameter / 2.
    return f


def human_total_resistance(mu,Dp,porosity,vessels,terminals,boundary_conds):
    # Calculates total resistance of the uterine arteries, outputs this resistance and a venous equivalent resistance (half of arterial resistance)
    resistance = np.zeros(np.size(vessels))
    flow = np.zeros(np.size(vessels) + 1)
    total_resistance = 0.0
    static_inlet_pressure = boundary_conds['inlet_p']  # Units of Pascals
    static_inlet_flow = boundary_conds['inlet_q'] * 1000. / 60.  # ml/min converted to mm^3/s
    pressure_out = np.zeros(np.size(vessels) + 1)

    anast_index = 0  # initialise
    radial_index = 0
    spiral_index = 0
    for i in range(0,np.size(vessels)):
        if (vessels['vessel_type'][i] == 'Anastomose'):
            anast_index = i
        elif (vessels['vessel_type'][i] == 'Radial'):
            radial_index = i
        elif (vessels['vessel_type'][i] == 'Spiral_tube'):
            spiral_index = i

    print("==============================")
    print("Resistance of each vessel type")
    print("==============================")

    for i in range(0, np.size(vessels)):
        # Poiseille resistance of each vessels
        # Units of resistance are Pa.s/mm^3
        if vessels['vessel_type'][i]=='Spiral_plug':
            resistance[i] =calc_plug_resistance(mu,Dp,porosity,vessels['radius'][i],vessels['length'][i])/vessels['number'][i]
        elif vessels['vessel_type'][i]=='Spiral_funnel':
            resistance[i] = calc_funnel_resistance(mu,vessels['radius'][spiral_index], vessels['radius'][i], 0.,vessels['length'][i]) / \
                            vessels['number'][i]
        else:
            resistance[i] = calc_tube_resistance(mu,vessels['radius'][i],vessels['length'][i])/vessels['number'][i]

        print(vessels['vessel_type'][i], resistance[i])

    print("=====================================")
    print("Flow and pressure in each vessel type")
    print("=====================================")
    if anast_index != 0:
        for i in range(0, np.size(vessels)):
            if i != anast_index:
                if (vessels['generation'][i] < vessels['generation'][anast_index]):
                    if i == 0:
                        flow[i] = static_inlet_flow / vessels['number'][
                            i]  # Assuming all arteries of a specific level are the same we just divide the flow
                        pressure_out[i] = static_inlet_pressure - resistance[i] * flow[i]
                    else:
                        flow[i] = flow[i - 1] * vessels['number'][i - 1] / vessels['number'][i]
                        pressure_out[i] = pressure_out[i - 1] - resistance[i] * static_inlet_flow

                    total_resistance = total_resistance + resistance[i]
                    print(str(vessels['vessel_type'][i]) + ' flow: ' + str(flow[i]) + ' mm^3/s ' + str(
                        flow[i] * 60. / 1000.) + ' ml/min ' + str(flow[i] * 60.) + ' ul/min ')
                    print(str(vessels['vessel_type'][i]) + ' pressure out: ' + str(pressure_out[i]) + ' Pa ' + str(
                        pressure_out[i] / 133.) + ' mmHg ')
    else:
        for i in range(0, np.size(vessels)):
            if i == 0:
                flow[i] = static_inlet_flow / vessels['number'][
                    i]  # Assuming all arteries of a specific level are the same we just divide the flow
                pressure_out[i] = static_inlet_pressure - resistance[i] * flow[i]
            else:
                flow[i] = flow[i - 1] * vessels['number'][i - 1] / vessels['number'][i]
                pressure_out[i] = pressure_out[i - 1] - resistance[i] * static_inlet_flow

            total_resistance = total_resistance + resistance[i]
            print(str(vessels['vessel_type'][i]) + ' flow: ' + str(flow[i]) + ' mm^3/s ' + str(
                flow[i] * 60. / 1000.) + ' ml/min ' + str(flow[i] * 60.) + ' ul/min ')
            print(str(vessels['vessel_type'][i]) + ' pressure out: ' + str(pressure_out[i]) + ' Pa ' + str(
                pressure_out[i] / 133.) + ' mmHg ')

    beyond_anast_resistance = 0.0
    if (vessels['length'][anast_index] == 0.0) or (anast_index == 0):
        venous_resistance = 0.0
        # Only IVS contribution to resistance (as in Mo et al. A transmission line modelling approach to the interpretation of uterine doppler waveforms. Ultrasound in Medicine and Biology, 1988. 14(5): p. 365-376.)
        parallel_resistance = terminals[0]  # Pa.s/mm^3
        flow[np.size(vessels)] = flow[np.size(vessels) - 2]
        pressure_out[np.size(vessels)] = pressure_out[np.size(vessels) - 2] - parallel_resistance * flow[
            np.size(vessels)] * vessels['number'][np.size(vessels) - 2]

        print('Terminal  flow: ' + str(flow[np.size(vessels)]) + ' mm^3/s ' + str(
            flow[np.size(vessels)] * 60. / 1000.) + ' ml/min ' + str(flow[np.size(vessels)] * 60.) + ' ul/min ')
        print('Terminal pressure out: ' + str(pressure_out[np.size(vessels)]) + ' Pa ' + str(
            pressure_out[np.size(vessels)] / 133.) + ' mmHg ')
    else:
        venous_resistance = total_resistance / 2.0  # Assuming veins are half as resistive as arteries
        for j in range(0, np.size(vessels)):
            if (vessels['generation'][j] > vessels['generation'][anast_index]):
                beyond_anast_resistance = beyond_anast_resistance + resistance[j]

        flow[anast_index] = (1. - (resistance[anast_index] + venous_resistance) / (
                    resistance[anast_index] + venous_resistance + terminals[2] * terminals[
                0] + beyond_anast_resistance)) * flow[anast_index - 1]


        flow[np.size(vessels)] = flow[anast_index - 1] - flow[anast_index]
        parallel_resistance = 1.0 / (1.0 / (resistance[anast_index] + venous_resistance) + terminals[2] * 1.0 / (
                    terminals[0] + beyond_anast_resistance * terminals[2]))

        kount_extra = 0
        for j in range(0, np.size(vessels)):
            if (vessels['generation'][j] > vessels['generation'][anast_index]):
                kount_extra = kount_extra + 1
                flow[j] = flow[np.size(vessels)] * terminals[2] / vessels['number'][j]
                pressure_out[j] = pressure_out[j - 2] - resistance[j] * flow[j]
                print(str(vessels['vessel_type'][j]) + ' flow: ' + str(flow[j]) + ' mm^3/s ' + str(
                    flow[i] * 60. / 1000.) + ' ml/min ' + str(flow[i] * 60.) + ' ul/min ')
                print(str(vessels['vessel_type'][j]) + ' pressure out: ' + str(pressure_out[j]) + ' Pa ' + str(
                    pressure_out[i] / 133.) + ' mmHg ')

        if kount_extra == 0:
            pressure_out[np.size(vessels)] = pressure_out[anast_index - 1] - terminals[2] * terminals[0] * flow[
                np.size(vessels)] * vessels['number'][anast_index - 1]
        else:
            pressure_out[np.size(vessels)] = pressure_out[anast_index + kount_extra] - terminals[2] * terminals[0] * \
                                             flow[np.size(vessels)] * vessels['number'][anast_index - 1]

        pressure_out[anast_index] = pressure_out[anast_index - 1] - (resistance[anast_index] + venous_resistance) * \
                                    flow[anast_index] * vessels['number'][anast_index]
        print('Anastomosis flow: ' + str(flow[anast_index]) + ' mm^3/s ' + str(
            flow[anast_index] * 60. / 1000.) + ' ml/min ' + str(flow[anast_index] * 60.) + ' ul/min ')
        print('Anastomosis pressure out: ' + str(pressure_out[anast_index]) + ' Pa ' + str(
            pressure_out[anast_index] / 133.) + ' mmHg ')
        print('Terminal  flow: ' + str(flow[np.size(vessels)]) + ' mm^3/s ' + str(
            flow[np.size(vessels)] * 60. / 1000.) + ' ml/min ' + str(flow[np.size(vessels)] * 60.) + ' ul/min ')
        print('Terminal pressure out: ' + str(pressure_out[np.size(vessels)]) + ' Pa ' + str(
            pressure_out[np.size(vessels)] / 133.) + ' mmHg ')
    total_resistance = total_resistance + parallel_resistance

    print("===============")
    print("Shear Stresses")
    print("===============")
    shear = np.zeros(np.size(vessels))
    for i in range(0, np.size(vessels)):
        if vessels['vessel_type'][i]=='Spiral_plug':
            if vessels['length'][i]>0.:
                [dummy,shear[i]]= calc_plug_shear(mu, Dp, porosity, vessels['radius'][i], vessels['length'][i], flow[i], resistance[i], np.linspace(0,vessels['radius'][i],10))
            else:
                shear[i]=0
        elif vessels['vessel_type'][i] == 'Spiral_funnel':
            shear[i]=(calc_tube_shear(mu,vessels['radius'][i],flow[i]) + calc_tube_shear(mu,vessels['radius'][spiral_index],flow[i]) )/2.#mean shear should be mid-way between inlet and outlet shear, but not accuarate
        else:
            shear[i] = calc_tube_shear(mu,vessels['radius'][i],flow[i])
        print(str(vessels['vessel_type'][i]) + ' shear: ' + str(shear[i]) + ' Pa ' + str(shear[i] / 10) + ' dyne/cm3')

    print("================")
    print("Total Resistance")
    print("================")
    print(str(total_resistance) + "Pa.mm^3/s")
    print("=============================")
    print("Total estimated pressure drop")
    print("=============================")
    print(str(total_resistance * static_inlet_flow) + " Pa, " + str(
        total_resistance * static_inlet_flow / 133.) + " mmHg")

    return [total_resistance, venous_resistance, shear,resistance,flow,pressure_out]