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
    
def calc_channel_resistance(mu,Dp, porosity, length,channel_rad,outer_rad):
    K = (Dp ** 2. / 180.) * ((porosity ** 3.) / (1. - porosity) ** 2.)  # permeability
    gamma = 1. / (1. + 2.5 * (1. - porosity))
    delta = gamma / K
    # special.iv(0,np.sqrt(delta) * radius) I0
    # special.iv(1,np.sqrt(delta) * radius) I1
    # special.kv(0,np.sqrt(delta) * radius) K0
    # special.kv(1,np.sqrt(delta) * radius) K1
    C2 = (-2. * K * np.sqrt(delta) * special.kv(1, np.sqrt(delta) * channel_rad) - channel_rad * special.kv(0, np.sqrt(
        delta) * outer_rad)) * special.iv(0, np.sqrt(delta) * channel_rad) / (2. * np.sqrt(delta) * (
                special.kv(0, np.sqrt(delta) * outer_rad) * special.iv(1, np.sqrt(delta) * channel_rad) + special.kv(1,
                                                                                                                     np.sqrt(
                                                                                                                         delta) * channel_rad) * special.iv(
            0, np.sqrt(delta) * outer_rad))) \
         - (2. * K * np.sqrt(delta) * special.iv(1, np.sqrt(delta) * channel_rad) - channel_rad * special.iv(0, np.sqrt(
        delta) * outer_rad)) * special.kv(0, np.sqrt(delta) * channel_rad) / (2. * np.sqrt(delta) * (
                special.kv(0, np.sqrt(delta) * outer_rad) * special.iv(1, np.sqrt(delta) * channel_rad) + special.kv(1,
                                                                                                                     np.sqrt(
                                                                                                                         delta) * channel_rad) * special.iv(
            0, np.sqrt(delta) * outer_rad))) + K + channel_rad ** 2 / 4.
    C2 = C2
    C3 = (-2. * K * np.sqrt(delta) * special.kv(1, np.sqrt(delta) * channel_rad) - channel_rad * special.kv(0, np.sqrt(
        delta) * outer_rad)) / (2. * np.sqrt(delta) * (
                special.kv(0, np.sqrt(delta) * outer_rad) * special.iv(1, np.sqrt(delta) * channel_rad) + special.kv(1,
                                                                                                                     np.sqrt(
                                                                                                                         delta) * channel_rad) * special.iv(
            0, np.sqrt(delta) * outer_rad)))
    C4 = -1. * (2. * K * np.sqrt(delta) * special.iv(1, np.sqrt(delta) * channel_rad) - channel_rad * special.iv(0,
                                                                                                                 np.sqrt(
                                                                                                                     delta) * outer_rad)) / (
                     2. * np.sqrt(delta) * (special.kv(0, np.sqrt(delta) * outer_rad) * special.iv(1, np.sqrt(
                 delta) * channel_rad) + special.kv(1, np.sqrt(delta) * channel_rad) * special.iv(0, np.sqrt(
                 delta) * outer_rad)))

    Q1 = ((2 * np.pi) / (mu * length)) * (C2 * channel_rad ** 2 / 2 - channel_rad ** 4 / 16)
    Q2 = (2 * np.pi / (np.sqrt(delta) * mu * length)) * (C3 * (
                outer_rad * special.iv(1, np.sqrt(delta) * channel_rad) - channel_rad * special.iv(1, np.sqrt(
            delta) * channel_rad)) - C4 * (outer_rad * special.kv(1, np.sqrt(
        delta) * outer_rad) - channel_rad * special.kv(1, np.sqrt(delta) * channel_rad)) + K * (
                                                                     outer_rad ** 2 - channel_rad ** 2) / 2)

    QT = Q1 + Q2
    resistance_channel = 1 / QT

    return resistance_channel


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
    #pressure is transmural pressure (kPa)
    #lengths are in um
    #everything else in SI units
    include_passive = True
    include_flow = True
    include_myo = True
    if np.array(fit_myo_params).all() == 0:
        include_myo = False
    if np.array(fit_flow_params).all() == 0:
        include_flow = False

    #Defining parameters of importance in passive model
    D0 = fit_passive_params[0] #um
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
        if(fixed_flow_params[4]==1):
           #Defining parameters related to blood flow through vessel, presribed flow
           Cshear = fit_flow_params[0]
           Cshear2 = -fit_flow_params[1]
           tau1 = fit_flow_params[2]
           tau2 = fit_flow_params[3]
           mu = fixed_flow_params[0] #viscosity Pa.s
           length = fixed_flow_params[1] #length of artery um
           flow = fixed_flow_params[2] #blood flow in artery !m3/s
           tau = calc_tube_shear(mu,diameter/2000000.,flow) #In Pa  
        else:
           #Defining parameters related to blood flow through vessel, prescribed pressure drop
           Cshear = fit_flow_params[0]
           Cshear2 = -fit_flow_params[1]
           tau1 = fit_flow_params[2]
           tau2 = fit_flow_params[3]
           mu = fixed_flow_params[0] #viscosity
           length = fixed_flow_params[1] #length of artery
           dp_blood = fixed_flow_params[2] #blood pressure drop in artery #in Pa
           system_resistance = fixed_flow_params[3] #any system resistance in myography, would be zero in a flow network model
           resistance = calc_tube_resistance(mu,diameter/2000000.,length/1000000.) + system_resistance# conversions take um to m #Pa.s/m3
           flow = dp_blood/resistance ## m3/s
           tau = calc_tube_shear(mu,diameter/2000000.,flow) #In Pa
        

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

        if not include_myo and not include_flow: #passive model
            diameter= bisection_method_diam(5.,fit_passive_params[0]*1.10,fit_passive_params,fit_myo_params,fit_flow_params,fixed_flow_params,pressure,verbose)
            if diameter <10.:
                 lowest_sign = find_possible_roots(0.,fit_passive_params[0]*2.0, fit_passive_params, fit_myo_params,fit_flow_params,fixed_flow_params, pressure,verbose)
                 diameter= bisection_method_diam(lowest_sign[0],lowest_sign[1],fit_passive_params,fit_myo_params,fit_flow_params,fixed_flow_params,pressure,verbose)
        else:
            diameter_pass = bisection_method_diam(5.,fit_passive_params[0]*1.10,fit_passive_params,[0.,0.],[0.,0.],[0.,0.],pressure,verbose)
            D0 = fit_passive_params[0]
            dp_blood = fixed_flow_params[2]
            reference_diameter = np.max([diameter_pass,D0])
            if(pressure>12.) or (dp_blood>0):
                #Looks for a diameter between 5 um and 10% larger than the passive diameter at that transmural pressure as a possible root
                lowest_sign = find_possible_roots(10.,reference_diameter*2.0, fit_passive_params, fit_myo_params,fit_flow_params,fixed_flow_params, pressure,verbose)
            else:
                lowest_sign=np.array([10.,reference_diameter*2.0])

            diameter = bisection_method_diam(lowest_sign[0], lowest_sign[1], fit_passive_params,fit_myo_params,fit_flow_params,fixed_flow_params, \
                                                                                                     pressure,verbose)
            #if verbose:
            if diameter<10:
                print('zero',pressure,dp_blood,lowest_sign,diameter_pass,reference_diameter)
                lowest_sign = find_possible_roots(0.,reference_diameter*2.0, fit_passive_params, fit_myo_params,fit_flow_params,fixed_flow_params, pressure,verbose)
                diameter = bisection_method_diam(lowest_sign[0], lowest_sign[1], fit_passive_params,fit_myo_params,fit_flow_params,fixed_flow_params, \
                                                                                                     pressure,verbose)
        return diameter

def find_possible_roots(low_diam,high_diam, fit_passive_params, fit_myo_params, fit_flow_params,fixed_flow_params, pressure,verbose):

    discretise = 1000
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
    if i==(discretise-1): #should be searching for a higher diameter
        #Has not changed sign in this range 
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
            if(i==(discretise -1)) and samesign:
                samesign = False

    lowest_signchange = np.zeros(2)
    lowest_signchange[0] = diameter_range[i]
    lowest_signchange[1] = diameter_range[i-1]

    return lowest_signchange

def tension_balance(fit_passive_params,fit_myo_params, fit_flow_params,fixed_flow_params, diameter, pressure):
    [Tmaxact,total_tension] = calc_total_tension(fit_passive_params,fit_myo_params, fit_flow_params,fixed_flow_params, diameter, pressure)
    f = total_tension - pressure * diameter / 2.
    return f


def human_total_resistance(mu,Dp,porosity,vessels,terminals,boundary_conds,channel_rad):
    # Calculates total resistance of the uterine arteries, outputs this resistance and a venous equivalent resistance (half of arterial resistance)
    resistance = np.zeros(np.size(vessels))
    flow = np.zeros(np.size(vessels) + 1)
    total_resistance = 0.0
    total_resistance_plus_myo = 0.
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
        elif (vessels['vessel_type'][i] == 'Uterine'):
            uterine_index = i
        elif (vessels['vessel_type'][i] == 'Arcuate'):
            arcuate_index = i

    print("==============================")
    print("Resistance of each vessel type")
    print("==============================")

    above_anast_resistance = 0.0
    beyond_anast_resistance = 0.0
    for i in range(0, np.size(vessels)):
        # Poiseille resistance of each vessels
        # Units of resistance are Pa.s/mm^3
        if vessels['vessel_type'][i]=='Spiral_plug':
            resistance[i] = calc_plug_resistance(mu,Dp,porosity,vessels['radius'][i],vessels['length'][i])/vessels['number'][i]
        elif vessels['vessel_type'][i]=='Spiral_funnel':
            resistance[i] = calc_funnel_resistance(mu,vessels['radius'][spiral_index], vessels['radius'][i], 0.,vessels['length'][i]) / \
                            vessels['number'][i]
        elif vessels['vessel_type'][i]=='Spiral_channel':
            if porosity > 0.15:
                resistance[i] = calc_channel_resistance(mu,Dp, porosity,vessels['length'][i],channel_rad,vessels['radius'][i])/vessels['number'][i]
            else: #calculate the equivalent channel resistance
                resistance[i] = calc_tube_resistance(mu,channel_rad,vessels['length'][i])/vessels['number'][i]
        else:
            resistance[i] = calc_tube_resistance(mu,vessels['radius'][i],vessels['length'][i])/vessels['number'][i]

        print(vessels['vessel_type'][i], resistance[i])

        if anast_index != 0:
            for i in range(0, np.size(vessels)):
                if i != anast_index:
                    if (vessels['generation'][i] < vessels['generation'][anast_index]):
                        above_anast_resistance = above_anast_resistance + resistance[i]
                    else:
                        beyond_anast_resistance = beyond_anast_resistance + resistance[i]
                else:
                    anast_resistance = resistance[i]

    print('Resistance, anastomosis',anast_resistance, 'above anast ',above_anast_resistance, 'beyond anast', beyond_anast_resistance)

    print("=====================================")
    print("------Calculating resistances--------")
    print("=====================================")
    uterine_resistance = resistance[uterine_index]
    feed_placenta_resistance = 0.0
    myometrial_resistance = terminals[3]

    if anast_index != 0:
        for i in range(0, np.size(vessels)):
            if i != anast_index:
                if (vessels['generation'][i] < vessels['generation'][anast_index]): #'branch 1, carries inlet flow'
                    if i != uterine_index:#exlclude uterine artery
                        feed_placenta_resistance = feed_placenta_resistance + resistance[i]
    else: #There is no anastomosis
        for i in range(0, np.size(vessels)):
            if i != uterine_index:#exlclude uterine artery
                    feed_placenta_resistance = feed_placenta_resistance + resistance[i]

    beyond_anast_resistance = 0.0
    if (vessels['length'][anast_index] == 0.0) or (anast_index == 0): #No anastomosis
        venous_resistance = 0.0
        # Only IVS contribution to resistance (as in Mo et al. A transmission line modelling approach to the interpretation of uterine doppler waveforms. Ultrasound in Medicine and Biology, 1988. 14(5): p. 365-376.)
        parallel_resistance = terminals[0]  # Pa.s/mm^3
    else:
        venous_resistance = feed_placenta_resistance / 2.0  # Assuming veins are half as resistive as arteries that lie upstream of the anastomoisi
        for j in range(0, np.size(vessels)):
            if (vessels['generation'][j] > vessels['generation'][anast_index]):
                beyond_anast_resistance = beyond_anast_resistance + resistance[j]
        venous_beyond_anast = beyond_anast_resistance / 2.
        parallel_resistance = 1.0 / (1.0 / (resistance[anast_index]) + terminals[2] * 1.0 / (
                    terminals[0] + (beyond_anast_resistance + venous_beyond_anast) * terminals[2]))


    feed_placenta_resistance = feed_placenta_resistance + parallel_resistance + venous_resistance
    total_resistance = uterine_resistance + 1./(1./feed_placenta_resistance + 1./myometrial_resistance) + uterine_resistance/2.
    print(uterine_resistance,feed_placenta_resistance,myometrial_resistance,total_resistance)
    print(beyond_anast_resistance,terminals[0]/terminals[2],venous_beyond_anast, resistance[anast_index],parallel_resistance)

    print("=====================================")
    print("------Flow divisions (no units)-----")
    print("=====================================")

    flow_prop_to_placenta = myometrial_resistance/(feed_placenta_resistance+myometrial_resistance)
    print("Percentage flow to placenta = " + str(flow_prop_to_placenta*100.))
    flow_prop_to_myometrium = 1. - flow_prop_to_placenta
    print("Percentage flow to myometrium = " + str(flow_prop_to_myometrium*100.))

    temp_inlet_flow = 1.

    if anast_index != 0:
        for i in range(0, np.size(vessels)):
            if i != anast_index:
                if (vessels['generation'][i] < vessels['generation'][anast_index]): #'branch 1, carries inlet flow'
                    if i == 0:
                        flow[i] = temp_inlet_flow / vessels['number'][i]  # Assuming all arteries of a specific level are the same we just divide the flow
                    else:
                        if i == arcuate_index:
                            flow[i] = flow_prop_to_placenta* flow[i - 1] * vessels['number'][i - 1] / vessels['number'][i]
                        else:
                            flow[i] = flow[i - 1] * vessels['number'][i - 1] / vessels['number'][i] #If vesels branch you divide the flow

                    print(str(vessels['vessel_type'][i]) + ' flow division: ' + str(flow[i]))
    else: #There is no anastomosis
        for i in range(0, np.size(vessels)):
            if i == 0:
                flow[i] = temp_inlet_flow / vessels['number'][i]  # Assuming all arteries of a specific level are the same we just divide the flow
            else:
                flow[i] = flow[i - 1] * vessels['number'][i - 1] / vessels['number'][i]

            print(str(vessels['vessel_type'][i]) + ' flow division: ' + str(flow[i]))

    if (vessels['length'][anast_index] == 0.0) or (anast_index == 0): #No anastomosis
        flow[np.size(vessels)] = flow[np.size(vessels) - 2] #need to check, all flow should go through the IVS
        print('Terminal  flow: ' + str(flow[np.size(vessels)]))
    else:
        flow[anast_index] = (1. - (resistance[anast_index]) / (resistance[anast_index]  + terminals[2] * terminals[
                0] + beyond_anast_resistance+venous_beyond_anast)) * flow[anast_index - 1]
        flow[np.size(vessels)] = flow[anast_index - 1] - flow[anast_index]
        for j in range(0, np.size(vessels)):
            if (vessels['generation'][j] > vessels['generation'][anast_index]):
                flow[j] = flow[np.size(vessels)] * terminals[2] * vessels['number'][anast_index - 1]/ vessels['number'][j]
                print(str(vessels['vessel_type'][j]) + ' flow division: ' + str(flow[j]))# + ' mm^3/s ' + str(

        print('Anastomosis flow division: ' + str(flow[anast_index]))# + ' mm^3/s ' + str(
        print('Terminal  flow division: ' + str(flow[np.size(vessels)]))# + ' mm^3/s ' + str(


    if boundary_conds['bc_type'] == 'flow':
        print(boundary_conds['bc_type'])
        static_inlet_pressure = boundary_conds['inlet_p']  # Units of Pascals
        static_inlet_flow = boundary_conds['inlet_q'] * 1000. / 60.  # ml/min converted to mm^3/s
        static_outlet_pressure = static_inlet_pressure - static_inlet_flow *total_resistance
    elif boundary_conds['bc_type'] == 'pressure':
        print(boundary_conds['bc_type'])
        static_inlet_pressure = boundary_conds['inlet_p']  # Units of Pascals
        static_outlet_pressure = boundary_conds['outlet_p']  # Units of Pascals
        static_inlet_flow = (static_inlet_pressure-static_outlet_pressure)/total_resistance

    print(flow)
    flow = flow*static_inlet_flow

    print(flow)
    print("==================")
    print("Flows (real units)")
    print("==================")

    for j in range(0, np.size(vessels)):
        print(str(vessels['vessel_type'][j]) + ' flow: ' + str(flow[j]) + ' mm^3/s ' + str(
                    flow[j] * 60. / 1000.) + ' ml/min ' + str(flow[j] * 60.) + ' ul/min ')
        print(str(vessels['vessel_type'][j]) + ' max velocity ' + str(2.*flow[j]/(np.pi*vessels['radius'][j]**2.)/1000.) + 'm/s')
    myometrial_flow = static_inlet_flow -flow[arcuate_index]*vessels['number'][arcuate_index]

    print('Flow feeding myometrium ' + str(myometrial_flow) + ' mm^3/s ' + str(
                    myometrial_flow * 60. / 1000.) + ' ml/min ' + str(myometrial_flow * 60.) + ' ul/min ')
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

        print(str(vessels['vessel_type'][i]) + ' shear: ' + str(shear[i]) + ' Pa ')

    print("================")
    print("Total Resistance")
    print("================")
    print(str(total_resistance) + "Pa.s/mm3")
    print("=============================")
    print("Total estimated pressure drop")
    print("=============================")
    print(str(total_resistance * static_inlet_flow) + " Pa, " + str(
        total_resistance * static_inlet_flow / 133.) + " mmHg")

    print("============")
    print("Total flow")
    print("============")
    print( str(static_inlet_flow) + " mm3/s, " + str(static_inlet_flow *60./1000) + " ml/min")



    return [total_resistance, venous_resistance, shear,resistance,flow,pressure_out]
    
    
def rat_total_resistance(mu,NumberPlacentae,vessels,terminals,boundary_conds,printme):
    #Calculates total resistance of the uterine arteries, outputs this resistance and a venous equivalent resistance (half of arterial resistance)
    resistance=np.zeros(np.size(vessels))
    flow = np.zeros(np.size(vessels)+1)
    total_resistance=0.0
    static_inlet_pressure = boundary_conds['inlet_p']  # Units of Pascals
    static_inlet_flow = boundary_conds['inlet_q'] * 1000. / 60.  # ml/min converted to mm^3/s
    pressure_out = np.zeros(np.size(vessels)+1)

    
    for i in range(0,np.size(vessels)):
        # Poiseille resistance of each vessels
        # Units of resistance are Pa.s/mm^3
        resistance[i]=81.0*mu*vessels['length'][i]/(8.0*np.pi* vessels['radius'][i]**4.0) /vessels['number'][i]
        #print(vessels['vessel_type'][i], resistance[i])
        if(vessels['vessel_type'][i]=='InUterine'):
            inlet_index = i
        elif(vessels['vessel_type'][i]=='Uterine'):
            uterine_index = i
        elif(vessels['vessel_type'][i]=='Arcuate'):
            arcuate_index = i
        elif (vessels['vessel_type'][i] == 'Radial'):
            radial_index = i

    ut_unit_resistance = 0.
    #We have a uterine segment in parallel to an arcuate from which branches everything else
    rad_sp_can_resistance = 0.
    if(mu*vessels['length'][arcuate_index]==0.):
        #No arcuate
        uterine_segment_resistance = resistance[uterine_index]/2.
        terminal_resistance = terminals[0]#Pa.s/mm^3
    else:
        arcuate_segment_resistance = resistance[arcuate_index] / 2.  # vessels['number'][radial_index]
        terminal_resistance = terminals[0]  # Pa.s/mm^3


    for i in range(arcuate_index+1, np.size(vessels)):
        rad_sp_can_resistance = rad_sp_can_resistance + resistance[i]

    if(mu*vessels['length'][arcuate_index]==0.):
        #No arcuate
        arc_and_ut_unit_resistance = uterine_segment_resistance + 1. / (
                    1. / uterine_segment_resistance + 1. / (rad_sp_can_resistance + terminal_resistance))
    else:
        total_arc_resistance =  arcuate_segment_resistance+ 1./(1./arcuate_segment_resistance + 1./(rad_sp_can_resistance + terminal_resistance))
        arc_and_ut_unit_resistance = 1./(1./total_arc_resistance + 1./resistance[uterine_index]) #add arcuate unit in parallel with uterine artery segment

    if	printme:
    	print("=====================================")
    	print("Flow and Pressure in each vessel type")
    	print("=====================================")
    for i in range(0, np.size(vessels)):
        if vessels['vessel_type'][i]=='InUterine':
            flow[i] = static_inlet_flow
            pressure_out[i] = static_inlet_pressure - flow[i]*resistance[i]
        elif vessels['vessel_type'][i] == 'Uterine':
            if (mu*vessels['length'][arcuate_index]==0.):
                flow[i] = (rad_sp_can_resistance + terminal_resistance)/(rad_sp_can_resistance + terminal_resistance+uterine_segment_resistance)*flow[inlet_index]
            else:
                flow[i] = (total_arc_resistance)/(total_arc_resistance+resistance[uterine_index])*flow[inlet_index]
            pressure_out[i] = pressure_out[i-1] - flow[i] * resistance[i]
        elif vessels['vessel_type'][i] == 'Arcuate':
            if (mu*vessels['length'][arcuate_index]==0.):
                flow[i] = (uterine_segment_resistance)/(rad_sp_can_resistance + terminal_resistance+uterine_segment_resistance)*flow[inlet_index]
                pressure_out[i] = pressure_out[inlet_index] - flow[inlet_index] * uterine_segment_resistance
            else:
                flow[i] = (resistance[uterine_index])/(total_arc_resistance+resistance[uterine_index])*flow[inlet_index]
                pressure_out[i] = pressure_out[inlet_index]-total_arc_resistance*flow[i]
        else:
            flow[i] = flow[arcuate_index]/vessels['number'][i] #flow per vessel
            if(i == radial_index):
                if (mu * vessels['length'][arcuate_index] == 0.):
                    pressure_in = pressure_out[arcuate_index]
                else:
                    pressure_in = pressure_out[inlet_index] - arcuate_segment_resistance*flow[arcuate_index]
                pressure_out[i] = pressure_in - flow[arcuate_index]*resistance[radial_index] #Assumption: Total arcuate flow goes through all the vessels below it
            else:
                pressure_out[i] = pressure_out[i-1] - resistance[i] * flow[arcuate_index]
        if mu * vessels['length'][arcuate_index] == 0 and vessels['vessel_type'][i] == 'Arcuate':
            print('No arcuate')
        else:
            if printme:
            	print(str(vessels['vessel_type'][i]) + ' flow: ' + str(flow[i]) + ' mm^3/s ' + str(flow[i]*60./1000.) + ' ml/min ' + str(flow[i]*60.) + ' ul/min ')
            	print(str(vessels['vessel_type'][i]) + ' pressure out: ' + str(pressure_out[i]) + ' Pa ' + str(
                pressure_out[i]/133.) + ' mmHg ')
    flow[np.size(vessels)]=flow[arcuate_index]
    pressure_out[np.size(vessels)] = pressure_out[np.size(vessels)-1]-terminal_resistance * flow[np.size(vessels)]
    if printme:
    	print('Terminals ' + ' flow: ' + str(flow[np.size(vessels)]) + ' mm^3/s ' + str(
        	flow[i] * 60. / 1000.) + ' ml/min ' + str(flow[np.size(vessels)] * 60.) + ' ul/min ')
    	print('Terminals '+ ' pressure out: ' + str(pressure_out[np.size(vessels)]) + ' Pa ' + str(
        pressure_out[np.size(vessels)] / 133.) + ' mmHg ')
    venous_resistance = 0.

    total_resistance = resistance[inlet_index] + arc_and_ut_unit_resistance * NumberPlacentae
    if printme:
    	print("===============")
    	print("Shear Stresses")
    	print("===============")
    shear = np.zeros(np.size(vessels))
    for i in range(0, np.size(vessels)):
        shear[i] = 4. * mu * flow[i] / (np.pi * vessels['radius'][i] ** 3.)
        if printme:
        	print(str(vessels['vessel_type'][i]) + ' shear: ' + str(shear[i]) + ' Pa ' + str(shear[i]/10) + ' dyne/cm3')

    if printme:
    	print("=====================================")
    	print("Total pressure drop along uterine horn")
    	print("======================================")
    dPress = static_inlet_flow *total_resistance #Pa
    if printme:
    	print(str(dPress) + ' Pa,' + str(dPress/133.) + ' mmHg' )

    total_resistance = resistance[inlet_index] + arc_and_ut_unit_resistance * NumberPlacentae

    if printme:
    	print("==========================================================")
    	print("Total pressure drop along single section of uterine artery")
    	print("==========================================================")
    dPress = static_inlet_flow *arc_and_ut_unit_resistance #Pa
    if printme:
        print(str(dPress) + ' Pa,' + str(dPress/133.) + ' mmHg' )
    
    if printme:
        print("================")
        print("Total Resistance")
        print("================")
        print(str(total_resistance) + "Pa.mm^3/s")

    return [total_resistance,venous_resistance,shear,resistance,flow,pressure_out] 


