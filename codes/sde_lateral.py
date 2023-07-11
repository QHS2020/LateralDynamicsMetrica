# -*- coding: utf-8 -*-

"""


"""

from RequiredModules import *



def Transformed_erf_integrator(y, lw = 3.5, rho = 0, alpha = 3.0, y_base = 0):
    """
    
    """
    
    tmp1 = -((y_base+rho*lw-y+lw/2.0)/(lw/alpha))**2
    tmp2 = -((y_base+rho*lw-y-lw/2.0)/(lw/alpha))**2
    
    return np.exp(tmp1) - np.exp(tmp2)
    
    pass


class SingleTrackModel:
    """
    
    """
    @classmethod
    def GetBeta(self, δ, lr, lf):
    
        #**************************Constrain the domain
        β = atan(lr*tan(δ) / (lr+lf))

        #**************************Constrain the domain
    #     if β > radians(60):
    #         β = radians(60)
    #     elif β < -radians(60):
    #         β = -radians(60)
        return β   
    
    
    @classmethod
    def GetSpeedPosition(self, v_last, β, ψ, vψ, a_long, dt, sx_last, sy_last):
        '''
            Calculate vehicle's velocity, acceleration, and position
        '''
        #*****************Speed
        v = a_long*dt + v_last
        vx = v*cos(β+ψ)
        vy = v*sin(β+ψ)

        #*****************Acceleration
        ax = a_long
        if sqrt((a_long)**2 + (v*vψ)**2) <=amax:
            ay = v*vψ
        else:
            ay = 0

        #*****************Position
        sx = vx*dt + 0.5*ax*dt**2 + sx_last
        sy = vy*dt + 0.5*ay*dt**2 + sy_last

        return sx, sy, v, vx, vy
    

    @classmethod
    def GetSteeringVelocity(self, δ_last, δ_lower, δ_upper, vδ_d, vδ_lower, vδ_upper, dt):
        '''
            Calculate the steering angle velocity of the vehicle: "vδ, δ"

            ----------------------------------------
            @input: 
                δ: steering angle
                δ_lower: the lower side of steering angle
                δ_upper: the upper side of steering angle

                vδ_d: desired steering angle velocity
                vδ_lower: the lower side of steering angle velocity
                vδ_upper: the upper side of steering angle velocity
                dt: simulation time step
        '''
        #*****************Calculate the steering angle velocity
        vδ = 0
        C1 = ((δ_last<=δ_lower) & (vδ_d<=0)) | ((δ_last>=δ_upper) & (vδ_d>=0))
        if C1:
            vδ = 0
        elif (not C1) & (vδ_d<=vδ_lower):
            vδ = vδ_lower
        elif (not C1) & (vδ_d>=vδ_upper):
            vδ = vδ_upper
        else:
            vδ = vδ_d

        #*****************Calculate the steering angle acceleration
        δ = δ_last + vδ*dt


        return vδ, δ

    
    @classmethod
    def GetLongitudinalAcceleration(self, a_long_d, vs, amax, v, v_lower, v_upper):
        '''
            Calculate the longitudinal acceleration of the vehicle: "a_long"

            ----------------------------------------
            @input: 
                a_long_d: desired acceleration
                vs: switching velocity
                amax: maximum acceleration

                v: speed
                v_lower: the lower side of speed
                v_upper: the upper side of speed
        '''
        #*****************Find the interval of the acceleration
        a_lower = 0
        a_upper = lambda v: amax*vs/v if v>vs else amax


        #*****************Calculate the acceleration
        a_long = 0

        C2 = ((v<=v_lower) & (a_long_d<=0)) | ((v>=v_upper) & (a_long_d>=0))
        if C2:
            a_long = 0
        elif (not C2) & (a_long_d<=a_lower):
            a_long = a_lower
        elif (not C2) & (a_long_d>=a_upper(v)):
            a_long = a_upper(v)
        else:
            a_long = a_long_d

        return a_long

    
    @classmethod
    def GetDirection(self, u, m, Iz, lf, lr, Csf, Csr, g, δ, a_long, hcg, β, v, ψ_last, vψ_last, dt):
        '''
            Calculate the parameter related to vehicle's direction: "ψ, vψ, aψ"
        '''
        #*****************Calculate the acceleration
        aψ = (u*m/(Iz*(lr+lf))) * ( lf*Csf*(g*lr-a_long*hcg)*δ + (lr*Csr*(g*lf+a_long*hcg) - lf*Csf*(g*lr-a_long*hcg))*β - 
                                   (lf**2*Csf*(g*lr-a_long*hcg) + lr**2*Csr*(g*lf+a_long*hcg)) * (vψ_last/v))

        #*****************Velocity
        vψ = vψ_last + aψ*dt

        #*****************Heading angle
        ψ = ψ_last + vψ_last*dt + 0.5*aψ*dt**2

        return ψ, vψ, aψ


    @classmethod
    def GetSlipAngleVelocity(self, u, v, lf, lr, Csf, Csr, g, δ, a_long, hcg, β, vψ, dt):
        '''
            Calculate the slip angle velocity: "vβ"
        '''
        #*****************Calculate the Velocity
        vβ = (u/(v*(lr+lf))) * ( Csf*(g*lr-a_long*hcg)*δ - (Csr*(g*lf+a_long*hcg) + Csf*(g*lr-a_long*hcg))*β +
                                  (Csr*(g*lf+a_long*hcg)*lr - Csf*(g*lr-a_long*hcg)*lf) * (vψ/v)) - vψ

        return vβ
    
    
    
    @classmethod
    def simulate(self, S, args, ts = np.linspace(0, 1 ,100)):
        """
        ----------------------------
        
        @input: S
            
            a list:
            
            S = [sx0, sy0, δ0, v0, ψ0, β0, vψ0, aψ0]
        
        
        ==============================================================================
        
        @input: args
            
            a dict:
            
            args['m'] is the mass of the vehicle. 
            
            args['lr'] and args['lf'] are the rear and front axis legnth. 
            ...
            
        
        ---------------------
        
        
        """
        #~~~~~~~~~~~~~~~~~~~~Convert 'dict' into varibles
        for i in range(len(args)):
            globals()[list(args.keys())[i]] = list(args.values())[i]

        
        #~~~~~~~~~~~~~~~~~~~~Simulation time
        ts = np.around(ts,3)
        dt = np.diff(ts)[0]
        
        
        #~~~~~~~~~~~~~~~~~~~~Compute the initial status
        sx0, sy0, δ0, v0, ψ0, β0, vψ0, aψ0 = S    ## the input initial values

        a_long_0 = SingleTrackModel.GetLongitudinalAcceleration(a_long_d, vs, amax, v0, v_lower, v_upper)
        vδ0 = SingleTrackModel.GetSteeringVelocity(δ0, δ_lower, δ_upper, vδ_d, vδ_lower, vδ_upper, dt)[0]
        vx0 = v0*cos(β0+ψ0) 
        vy0 = v0*sin(β0+ψ0)
        vβ0 = SingleTrackModel.GetSlipAngleVelocity(u, v0, lf, lr, Csf, Csr, g, δ0, a_long_0, hcg, β0, vψ0, dt)
        
        Status = {}
        Status[0] = [sx0, sy0, δ0, vδ0, a_long_0, v0, vx0, vy0, ψ0, vψ0, aψ0, β0, vβ0]


        #~~~~~~~~~~~~~~~~~~~~Simulation
        for i in range(1, len(ts)):

            #--------------------------------------------------------Get the last vehicle's status
            [sx_last, sy_last, δ_last, vδ_last, a_long_last, v_last, vx_last, vy_last, ψ_last, vψ_last, aψ_last, β_last, vβ_last] = Status[ts[i-1]]

            #--------------------------------------------------------Get new acceleration: a_long
            a_long = SingleTrackModel.GetLongitudinalAcceleration(a_long_d, vs, amax, v_last, v_lower, v_upper)

            #--------------------------------------------------------Get steering angle velocity: β
            vδ, δ = SingleTrackModel.GetSteeringVelocity(δ_last, δ_lower, δ_upper, vδ_d, vδ_lower, vδ_upper, dt)

            #--------------------------------------------------------Get slip angle: β
            β = SingleTrackModel.GetBeta(δ, lf, lr)

            #--------------------------------------------------------Get vehicle's direction: ψ, vψ, aψ
            ψ, vψ, aψ = SingleTrackModel.GetDirection(u, m, Iz, lf, lr, Csf, Csr, g, δ, a_long, hcg, β, v_last, ψ_last, vψ_last, dt)

            #--------------------------------------------------------Get slip angle velocity: vβ
            vβ = SingleTrackModel.GetSlipAngleVelocity(u, v_last, lf, lr, Csf, Csr, g, δ, a_long, hcg, β, vψ, dt)

            #--------------------------------------------------------Get position, velocity, and acceleration
            sx, sy, v, vx, vy = SingleTrackModel.GetSpeedPosition(v_last, β, ψ, vψ, a_long, dt, sx_last, sy_last)

            #--------------------------------------------------------Update vehicle's status
            Status[ts[i]] = [sx, sy, δ, vδ, a_long_0, v, vx, vy, ψ, vψ, aψ, β, vβ]


        #~~~~~~~~~~~~~~~~~~~~Convert results into dataframe
        results = pd.DataFrame.from_dict(Status).T
        results.columns = ['sx', 'sy', 'δ', 'vδ', 'a_long', 'v', 'vx', 'vy', 'ψ', 'vψ', 'aψ', 'β', 'vβ']     
        
        return results
    
    @classmethod
    def F(self, S, args, ts = np.linspace(0, 1 ,100)):
        """
        ----------------------------
        
        @input: S
            
            a list:
            
            S = [sx0, sy0, δ0, v0, ψ0, β0, vψ0, aψ0]
        
        
        ==============================================================================
        
        @input: args
            
            a dict:
            
            args['m'] is the mass of the vehicle. 
            
            args['lr'] and args['lf'] are the rear and front axis legnth. 
            ...
            
        
        ---------------------
        
        
        """
        #~~~~~~~~~~~~~~~~~~~~Convert 'dict' into varibles
        for i in range(len(args)):
            globals()[list(args.keys())[i]] = list(args.values())[i]

        
        #~~~~~~~~~~~~~~~~~~~~Simulation time
        ts = np.around(ts,3)
        dt = np.diff(ts)[0]
        
        
        #~~~~~~~~~~~~~~~~~~~~Compute the initial status
        sx0, sy0, δ0, v0, ψ0, β0, vψ0, aψ0 = S    ## the input initial values

        a_long_0 = SingleTrackModel.GetLongitudinalAcceleration(a_long_d, vs, amax, v0, v_lower, v_upper)
        vδ0 = SingleTrackModel.GetSteeringVelocity(δ0, δ_lower, δ_upper, vδ_d, vδ_lower, vδ_upper, dt)[0]
        vx0 = v0*cos(β0+ψ0) 
        vy0 = v0*sin(β0+ψ0)
        vβ0 = SingleTrackModel.GetSlipAngleVelocity(u, v0, lf, lr, Csf, Csr, g, δ0, a_long_0, hcg, β0, vψ0, dt)
        
        Status = {}
        Status[0] = [sx0, sy0, δ0, vδ0, a_long_0, v0, vx0, vy0, ψ0, vψ0, aψ0, β0, vβ0]


        #~~~~~~~~~~~~~~~~~~~~Simulation
        for i in range(1, len(ts)):

            #--------------------------------------------------------Get the last vehicle's status
            [sx_last, sy_last, δ_last, vδ_last, a_long_last, v_last, vx_last, vy_last, ψ_last, vψ_last, aψ_last, β_last, vβ_last] = Status[ts[i-1]]

            #--------------------------------------------------------Get new acceleration: a_long
            a_long = SingleTrackModel.GetLongitudinalAcceleration(a_long_d, vs, amax, v_last, v_lower, v_upper)

            #--------------------------------------------------------Get steering angle velocity: β
            vδ, δ = SingleTrackModel.GetSteeringVelocity(δ_last, δ_lower, δ_upper, vδ_d, vδ_lower, vδ_upper, dt)

            #--------------------------------------------------------Get slip angle: β
            β = SingleTrackModel.GetBeta(δ, lf, lr)

            #--------------------------------------------------------Get vehicle's direction: ψ, vψ, aψ
            ψ, vψ, aψ = SingleTrackModel.GetDirection(u, m, Iz, lf, lr, Csf, Csr, g, δ, a_long, hcg, β, v_last, ψ_last, vψ_last, dt)

            #--------------------------------------------------------Get slip angle velocity: vβ
            vβ = SingleTrackModel.GetSlipAngleVelocity(u, v_last, lf, lr, Csf, Csr, g, δ, a_long, hcg, β, vψ, dt)

            #--------------------------------------------------------Get position, velocity, and acceleration
            sx, sy, v, vx, vy = SingleTrackModel.GetSpeedPosition(v_last, β, ψ, vψ, a_long, dt, sx_last, sy_last)

            #--------------------------------------------------------Update vehicle's status
            Status[ts[i]] = [sx, sy, δ, vδ, a_long_0, v, vx, vy, ψ, vψ, aψ, β, vβ]


        #~~~~~~~~~~~~~~~~~~~~Convert results into dataframe
        results = pd.DataFrame.from_dict(Status).T
        results.columns = ['sx', 'sy', 'δ', 'vδ', 'a_long', 'v', 'vx', 'vy', 'ψ', 'vψ', 'aψ', 'β', 'vβ']   
        Data = results[['sx', 'sy', 'δ', 'v', 'ψ', 'β', 'vψ', 'aψ']].values[-1,:]
        
        return Data
    

class PlotVehicle:
    
    @classmethod  
    def Arrow(self, x, y, theta, L, c, alpha):
        '''
            Plot the vehicle's direction
            -----------------------------
            
            @input:
                
                x: x-position
                y: y-position
                theta: 
                c: color
                alpha: transparency coefficient
                
        '''
        PI = np.pi
        angle = np.deg2rad(30)
        d = 0.5 * L
        w = 2

        x_start = x
        y_start = y
        x_end = x + L * np.cos(theta)
        y_end = y + L * np.sin(theta)

        theta_hat_L = theta + PI - angle
        theta_hat_R = theta + PI + angle

        x_hat_start = x_end
        x_hat_end_L = x_hat_start + d * np.cos(theta_hat_L)
        x_hat_end_R = x_hat_start + d * np.cos(theta_hat_R)

        y_hat_start = y_end
        y_hat_end_L = y_hat_start + d * np.sin(theta_hat_L)
        y_hat_end_R = y_hat_start + d * np.sin(theta_hat_R)

        plt.plot([x_start, x_end], [y_start, y_end], color=c, linewidth=w, alpha=alpha)
        plt.plot([x_hat_start, x_hat_end_L],
                 [y_hat_start, y_hat_end_L], color=c, linewidth=w, alpha=alpha)
        plt.plot([x_hat_start, x_hat_end_R],
                 [y_hat_start, y_hat_end_R], color=c, linewidth=w, alpha=alpha)
    
    
    @classmethod
    def Car(self, x, y, yaw, w, L, alpha):
        '''
            Plot the counter of vehicle
            ----------------------------
            
            @input:
                
                x: x-position
                y: y-position
                yaw: the yaw angle
                w: vehicle's width
                L: vehicle's length
                alpha: transparency coefficient
            
            
        '''
        PI = np.pi
        
        theta_B = PI + yaw

        xB = x + L / 2 * np.cos(theta_B)
        yB = y + L / 2 * np.sin(theta_B)

        theta_BL = theta_B + PI / 2
        theta_BR = theta_B - PI / 2

        x_BL = xB + w / 2 * np.cos(theta_BL)        # Bottom-Left vertex
        y_BL = yB + w / 2 * np.sin(theta_BL)
        x_BR = xB + w / 2 * np.cos(theta_BR)        # Bottom-Right vertex
        y_BR = yB + w / 2 * np.sin(theta_BR)

        x_FL = x_BL + L * np.cos(yaw)               # Front-Left vertex
        y_FL = y_BL + L * np.sin(yaw)
        x_FR = x_BR + L * np.cos(yaw)               # Front-Right vertex
        y_FR = y_BR + L * np.sin(yaw)

        plt.plot([x_BL, x_BR, x_FR, x_FL, x_BL],
                 [y_BL, y_BR, y_FR, y_FL, y_BL],
                 linewidth=1.5, color='black', alpha=alpha)
        
        #---------------------- plot the arrow
        PlotVehicle.Arrow(x, y, yaw, L/2, 'black', alpha)






def Transformed_erf_integrator_assymetry(y, lw = 3.5, rho = 0, alpha_l = 3.0, alpha_r = 4.0,y_base = 0):
    """
    
    Difference:
    
        - Transformed_erf_integrator(), it is symmetric. 
        - Transformed_erf_integrator_assymetry(), )The non-symmetry relationship.
    
    ----------------------------------------------
    The assymetric relationship is given by alpha_l and alpha_r. alpha_l correspond to the negative-y case and alpha_r correspons to the positive-y case. 
    
    
    """
    
    tmp1 = -((y_base+rho*lw-y+lw/2.0)/(lw/alpha_l))**2
    tmp2 = -((y_base+rho*lw-y-lw/2.0)/(lw/alpha_r))**2
    
    return np.exp(tmp1) - np.exp(tmp2)
    
    pass


def mc_single_point_sum(log_probs):
    #print(np.format_float_scientific(log_probs[0]).split('+')[1])
    #given many points from monte carlao. calculate the log of the mean of the prob. 
    #
    digitsNs = []
    for i in log_probs:
        if i<=-1:
            #print(i, np.format_float_scientific(i).split('+'), '--')
            #print(int(np.format_float_scientific(i).split('+')[1]))
            digitsNs.append(int(np.format_float_scientific(i).split('+')[1]))
        else:
            #when i<=-1, then np.format_float_scientific(-.9) return '-9.e-01'
            digitsNs.append(0)
    #digitsNs = [int(np.format_float_scientific(i).split('+')[1]) for i in log_probs]
    digitsMIN = min(digitsNs)
    
    #the condition >=-10 is to omit the smaller number. 
    sample_logs_added = [i for i in 10**(digitsMIN)+log_probs if i>=-10]
    
    print(sample_logs_added)
    
    return -10**digitsMIN+np.log(np.mean(np.exp(sample_logs_added)))






def NOTES():
    """
    GET the max log(p) using Monte carla:
    
    
    lw = 3.5
    sizemonte = 1e7                                                                                                                                                                                                                                                                                 
        
    
    
    
    
    """
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    pass



class dataanalysis():
    """
    
    """
    
    
    @classmethod
    def LaneKeepingMove2Zero_NGSIM(self, observationsRAW, lanemarks_coors = np.array([0, 11.98, 24.05, 35.84, 48.15, 60.11])/3.2808, moment_col = 'Frame_ID', lateral_col = 'Local_X', lw_normalized = 3.5, vid_col = 'Vehicle_ID', mirror_trajectory = True, filter_UP = 6, filter_LW = -2, lat_to_meter_coefficient = 1/3.2808):
        """
        feet convert to meter: ft/3.2808
        
        @input: observationsRAW
            observationsRAW[vid] = data, a dataframe. 
            
        @input: mirror_trajectory
            make sure that the vehicle turn to the positive lateral coordinate. 
        
        """
        #----------------------------------
        #keys are vehicle ids. 
        observations = {}
        
        #----------------------------------
        #lanesCoors, boundary of each lane. laneid is labeled from 0~N
        lanecoors = {}
        for i in range(1, len(lanemarks_coors)):
            lanecoors[i-1] = (lanemarks_coors[i-1], lanemarks_coors[i])
        
        #------------------------------------------
        def findlaneid(loc, lanecoors):
            #find the lane id the vehicle is in. 
            for k in lanecoors.keys():
                if loc>=min(lanecoors[k]) and loc<max(lanecoors[k]):
                    return k
            return False
        #------------------------------------------
        
        #
        for vid in observationsRAW.keys():
            vdata = observationsRAW[vid]
            
            #-------get the lateral data in a Series, 
            lateraldata  = (vdata[lateral_col].values)*lat_to_meter_coefficient
            tdata = vdata[moment_col].values
            tdata = tdata-min(tdata)
            #
            observation = pd.Series(lateraldata, index = tdata).sort_index()
            
            #-----shift to the laterl = 0.
            laneid = findlaneid(observation.iloc[0], lanecoors)
            if laneid==False:continue
            #
            shift = (min(lanecoors[laneid]) + max(lanecoors[laneid]))/2.0
            lw =  max(lanecoors[laneid])-  min(lanecoors[laneid])
            #print(laneid,shift)
            observation0 = observation-shift
            observation0 = observation0*lw_normalized/lw
            
            #
            if mirror_trajectory:
                if observation0.iloc[-1]<-11.98/3.2808/2:
                    observation0 = -1*observation0
            #save and assign. 
            if not isinstance(filter_UP, bool):
                if max(observation0)>filter_UP:continue
            if not isinstance(filter_LW, bool):
                if min(observation0)<filter_LW:continue
            
            
            observations[vid] = observation0
        
        return observations

    @classmethod
    def Convertobservation(self, observation0, lanes_bounds = {1:[8.09, 11.59], 2:[ 11.59, 14.5], 3:[20.45, 23.97], 4:[23.97, 27]}):
        """
        Convert the observation that start from lane with middle line==0. and lane changing to positive lat coordinates. 
        
        
        [8.09, 11.59, 14.5, 20.45, 23.97, 27]
        
        lane_bounds = ( [8.09, 11.59], [ 11.59, 14.5], [20.45, 23.97], [23.97, 27])
        
        
        @input: observation0
         
         a pd.series. The index are treated as the time. 
        
        @input: 
        
        """
        #determine which lane the loc is in. 
        def determin_lane(loc, lanes = [1,2,3,4], lanes_bounds = ( [8.09, 11.59], [ 11.59, 14.5], [20.45, 23.97], [23.97, 27])):
            for i,(s,e) in zip(lanes, lanes_bounds):
                #print(loc)
                if loc>=min(s,e) and loc<max(s,e):
                    return i
            return False
        
        #
        observation = observation0.sort_index()
        
        #
        #print(observation.iloc[0],'---')
        lanes = list(lanes_bounds.keys())
        lanes_boundstmp = [lanes_bounds[l] for l in lanes]
        startlane = determin_lane(observation.iloc[0], lanes_bounds=lanes_boundstmp,lanes = lanes)
        endlane = determin_lane(observation.iloc[-1], lanes_bounds=lanes_boundstmp,lanes = lanes)
        
        if startlane==endlane or startlane==False or endlane==False:
            return False,False
        #
        #print(startlane, endlane, observation.iloc[0], observation.iloc[-1])
        lanelw, laneup =  lanes_bounds[startlane][0],lanes_bounds[startlane][1]
        observation = observation - (lanelw + laneup)/2.0
        
        if observation.iloc[-1]<0:
            observation= -1*observation
            
        return True,observation

    @classmethod
    def Convertobservations(self, raw_observations, lanes_bounds = {1:[8.09, 11.59], 2:[ 11.59, 14.5], 3:[20.45, 23.97], 4:[23.97, 27]}, filter_up = False, filter_lw = False, ):
        """
        ---------------------------------------------------------------
        For ngsim data. the callback is:
            
            #
            - Convertobservations(self, raw_observations, lanes_bounds = {1:[8.09, 11.59], 2:[ 11.59, 14.5], 3:[20.45, 23.97], 4:[23.97, 27]}, filter_up = False, filter_lw = False)
        
        
        ---------------------------------------------------------------
        @input: raw_observations
            
            a dict. 
            
            rach value is a pd.Series. 
        
        @input: filter_up and filter_lw
            
            the extra condition to filter the observation. 
            
            They should be either bool or float. 
            
        """
        
        observations = {}
        for vid in raw_observations.keys():
            
            successornot,o= self.Convertobservation( observation0 = raw_observations[vid], lanes_bounds = lanes_bounds, )
            if not successornot:continue
            if not isinstance(filter_up, bool):
                if max(o)>=filter_up:continue
            if not isinstance(filter_lw, bool):
                if min(o)<=filter_lw:continue
            observations[vid] = o
            pass
        
        
        return observations
    
    @classmethod
    def USELESS_aggregate_trs2dp_of_lateral(self, trs_dict_vids, lateral_col = 'Local_X', moment_col = 'Frame_ID', min_horizon_sec = 100):
        """
        aggregate the data to dataframe. 
        
        """
        #filter the data so that only long enough data is used
        trs_dict_vids_filtered = {}
        for vid in trs_dict_vids:
            if max(trs_dict_vids[vid][moment_col])<=min_horizon_sec:continue
            trs_dict_vids_filtered[vid] = trs_dict_vids[vid]
        
        #find the common moments
        allmoments = set()
        for vid in trs_dict_vids_filtered:
            allmoments.update(trs_dict_vids_filtered[vid][moment_col])
        
        commonmoments = copy.deepcopy(allmoments)
        for vid in trs_dict_vids_filtered:
            commonmoments.intersection_update(trs_dict_vids_filtered[vid][moment_col])
        
        #
        datas = []
        for vid in trs_dict_vids_filtered:
            partialdata0 = trs_dict_vids_filtered[vid][lateral_col]
            partialdata0.index = trs_dict_vids_filtered[vid][moment_col].values
            
            builtins.tmp = partialdata0,commonmoments,allmoments
            datas.append(partialdata0[commonmoments])
        
        
        return datas
        return pd.DataFrame(datas)
        
    
    @classmethod
    def normalize_lateral_to_zero_given_dicts_of_df(self, trs_dict_vids, moment_col = 'Frame_ID', lateral_col = 'Local_X', lanemarks_coors_meter = np.array([0, 11.98, 24.05, 35.84, 48.15, 60.11])/3.2808, lateral_scale_to_meter = 1/3.2808, moment_scale_to_second = .1, normalized_lw = 3.5, ):
        """
        normalize the data to lat at zero. 
        
        CALLBACK:
        
        
            data1 = pd.read_csv(dataspath + 'zjh/' + 'ngsim/NgsimLaneKeepData_1.csv')
            data2 = pd.read_csv(dataspath + 'zjh/' + 'ngsim/NgsimLaneKeepData_2.csv')
            data3 = pd.read_csv(dataspath + 'zjh/' + 'ngsim/NgsimLaneKeepData_3.csv')
            trs_ngsim = {}
            for d in [data1,data2,data3]:
                vids = set(d.Vehicle_ID)
                for vid in vids:
                    vdata = d[d.Vehicle_ID==vid]
                    trs_ngsim[vid] = vdata[['Frame_ID', 'Local_X', 'Local_Y']]
                    
            neww = sdelat.dataanalysis.normalize_lateral_to_zero_given_dicts_of_df(trs_ngsim)

            ax= sdelat.dataanalysis.plot_trs(neww, N = 300)
                    
        
            ++++++++++++++++++++++++++++++++++++++++++++++++++
            lanes of ngsim is lanemarks_coors_meter = np.array([0, 11.98, 24.05, 35.84, 48.15, 60.11])/3.2808
            lanes of highD is lanemarks_coors_meter = np.array([8.09, 11.59, 14.5, 20.45, 23.97, 27])
        
        
        ------------------------------------
        @input: trs_dict_vids
            
            a dict. 
            
            trs_dict_vids[vid] is a dataframe. info should be include t,x,y.
        
        @OUTPUT: normalized_trs_dict_vids
        
            a dict. 
            
            normalized_trs_dict_vids[vid] is a dataframe. 
            
            
            
        
        """
        #
        normalized_trs_dict_vids = {}
        
        
        #----------------------------------
        #lanesCoors, boundary of each lane. laneid is labeled from 0~N
        lanecoors = {}
        for i in range(1, len(lanemarks_coors_meter)):
            lanecoors[i-1] = (lanemarks_coors_meter[i-1], lanemarks_coors_meter[i])
        
        #------------------------------------------
        def findlaneid(loc, lanecoors):
            #find the lane id the vehicle is in. 
            for k in lanecoors.keys():
                if loc>=min(lanecoors[k]) and loc<max(lanecoors[k]):
                    return k
            return False
        
        #---------------------------------------
        
        for vid in trs_dict_vids.keys():
            #a dataframe. 
            vdata = trs_dict_vids[vid]
            
            #-------get the lateral data in a Series, 
            lateraldata  = (vdata[lateral_col].values)*lateral_scale_to_meter
            tdata = vdata[moment_col].values
            tdata = (tdata-min(tdata))*moment_scale_to_second
            
            #
            observation = pd.Series(lateraldata, index = tdata).sort_index()
            
            #-----shift to the laterl = 0.
            laneid = findlaneid(lateraldata[0], lanecoors)
            if laneid==False:continue
            #
            lw = (max(lanecoors[laneid]) - min(lanecoors[laneid]))/2.0
            middle = (max(lanecoors[laneid]) + min(lanecoors[laneid]))/2.0
            #shift = laneid*lw + lw/2.0
            #print(laneid,shift)
            lateraldata_new = (lateraldata - middle)*normalized_lw/lw
            
            #
            if max(lateraldata_new)>lw or min(lateraldata_new)<-lw:continue
            
            #assign values
            normalized_trs_dict_vids[vid] = copy.deepcopy(trs_dict_vids[vid])
            normalized_trs_dict_vids[vid][lateral_col] = lateraldata_new
            normalized_trs_dict_vids[vid][moment_col] = tdata
            normalized_trs_dict_vids[vid].sort_values(by = moment_col, inplace=True)
        
        return normalized_trs_dict_vids
        
        pass
    
    
    @classmethod
    def Convert_trs_dict2observation(self, trs_dict,  moment_col = 'frame', lateral_col = 'y', T_horizon_sec_min  = 10, LW = False, UP = False):
        """
        @input: LW and UP
        
            used to filter the trajectory
        
        """
        
        observations = {}
        
        
        for vid in trs_dict.keys():
            ts  = trs_dict[vid][moment_col].values
            if max(ts)-min(ts)<=T_horizon_sec_min:continue
            lats = copy.deepcopy(trs_dict[vid][lateral_col].values) 
        
            tmp = pd.Series( lats, index = ts)
            
            if not isinstance(LW, bool):
                if min(tmp)<LW:
                    continue
                
            if not isinstance(UP, bool):
                if max(tmp)>UP:
                    continue
            
            observations[vid] = pd.Series( lats, index = ts)
        
        return observations
    
    
    @classmethod
    def BKP_normalize_lateral_to_zero_given_dicts_of_df(self, trs_dict_vids, moment_col = 'Frame_ID', lateral_col = 'Local_X', lanemarks_coors_meter = np.array([0, 11.98, 24.05, 35.84, 48.15, 60.11])/3.2808, lateral_scale_to_meter = 1/3.2808, moment_scale_to_second = .1, lw = 3.6, normalized_lw = 3.5):
        """
        normalize the data to lat at zero. 
        
        CALLBACK:
        
        
            data1 = pd.read_csv(dataspath + 'zjh/' + 'ngsim/NgsimLaneKeepData_1.csv')
            data2 = pd.read_csv(dataspath + 'zjh/' + 'ngsim/NgsimLaneKeepData_2.csv')
            data3 = pd.read_csv(dataspath + 'zjh/' + 'ngsim/NgsimLaneKeepData_3.csv')
            trs_ngsim = {}
            for d in [data1,data2,data3]:
                vids = set(d.Vehicle_ID)
                for vid in vids:
                    vdata = d[d.Vehicle_ID==vid]
                    trs_ngsim[vid] = vdata[['Frame_ID', 'Local_X', 'Local_Y']]
                    
            neww = sdelat.dataanalysis.normalize_lateral_to_zero_given_dicts_of_df(trs_ngsim)

            ax= sdelat.dataanalysis.plot_trs(neww, N = 300)
                    
        
        
        ------------------------------------
        @input: trs_dict_vids
            
            a dict. 
            
            trs_dict_vids[vid] is a dataframe. info should be include t,x,y.
        
        @output: normalized_trs_dict_vids
        
            a dict. 
            
            normalized_trs_dict_vids[vid] is a dataframe. 
            
        
        """
        #
        normalized_trs_dict_vids = {}
        
        
        #----------------------------------
        #lanesCoors, boundary of each lane. laneid is labeled from 0~N
        lanecoors = {}
        for i in range(1, len(lanemarks_coors_meter)):
            lanecoors[i-1] = (lanemarks_coors_meter[i-1], lanemarks_coors_meter[i])
        
        #------------------------------------------
        def findlaneid(loc, lanecoors):
            #find the lane id the vehicle is in. 
            for k in lanecoors.keys():
                if loc>=min(lanecoors[k]) and loc<max(lanecoors[k]):
                    return k
            return False
        
        #---------------------------------------
        
        for vid in trs_dict_vids.keys():
            #a dataframe. 
            vdata = trs_dict_vids[vid]
            
            #-------get the lateral data in a Series, 
            lateraldata  = (vdata[lateral_col].values)*lateral_scale_to_meter
            tdata = vdata[moment_col].values
            tdata = (tdata-min(tdata))*moment_scale_to_second
            #
            observation = pd.Series(lateraldata, index = tdata).sort_index()
            
            #-----shift to the laterl = 0.
            laneid = findlaneid(lateraldata[0], lanecoors)
            if laneid==False:continue
            #
            shift = laneid*lw + lw/2.0
            #print(laneid,shift)
            lateraldata_new = lateraldata-shift
            
            #assign values
            normalized_trs_dict_vids[vid] = copy.deepcopy(trs_dict_vids[vid])
            normalized_trs_dict_vids[vid][lateral_col] = lateraldata_new
            normalized_trs_dict_vids[vid][moment_col] = tdata
            normalized_trs_dict_vids[vid].sort_values(by = moment_col, inplace=True)
        
        return normalized_trs_dict_vids
        
        pass
    

    
    @classmethod
    def plot_observations(self, observations, N = np.inf, figsize = (5,4), ax  =False):
        """
        
        
        """
        if ax==False:
            fig,ax = plt.subplots(figsize = figsize)
        i = 1
        for vid in [np.random.choice(list(observations.keys())) for i in range(int(min(len(observations), N)))]:
            ax.plot( observations[vid].index , observations[vid].values)
            i = i+1
            if i>=N:
                ax.grid()
                return ax
        ax.grid()
        return ax
    
    @classmethod
    def plot_trs(self, trs_dict_vids, ax = False, plot_x_col = 'Local_X', plot_y_col = 'Local_Y', N = 100, figsize = (5,5)):
        """
        
        
        @input: trs_dict_vids
            
            trs_dict_vids[vid] = dataframe. columns include t,x,y
        
        """
        if ax == False:
            fig,ax = plt.subplots(figsize = figsize)
        for vid in np.random.choice(list(trs_dict_vids.keys()), min(N, len(trs_dict_vids))):
            ax.plot(trs_dict_vids[vid][plot_x_col], trs_dict_vids[vid][plot_y_col])
        
        return ax
        pass
    
    @classmethod
    def plot_distributions_dict_tmp(self, distributions, post_legend = ' sec', normalize = True, ax = False, figsize = (3,3)):
        """
        when normalize, just divide the number, 
        
        """
        if isinstance(ax,bool):
            fig,ax = plt.subplots(figsize = figsize)
        
        for k in distributions.keys():
            hs,es = distributions[k]
            if normalize:
                ax.plot(es,hs/sum(hs)/(es[-1]-es[-2]), label = str(k)+ post_legend)
                
            else:
                ax.plot(es, hs, label = str(k)+ post_legend)
                
        ax.legend()
        return ax
        pass

    @classmethod
    def plot_distributions_dict(self, distributions, post_legend = ' sec', ax = False,normalize = True, figsize = (5,3)):
        """
        
        """
        if isinstance(ax,bool):
            fig,ax = plt.subplots(figsize = figsize)
        
        for k in distributions.keys():
            hs,es = distributions[k]
            if normalize:
                ax.plot(es,hs/sum(hs)/(es[-1]-es[-2]), label = str(k)+ post_legend)
                
            else:
                ax.plot(es, hs, label = str(k)+ post_legend)
                
        ax.legend()
        return ax
        pass
    
    
    
    @classmethod
    def filter_by_lateral_init_loc(self, observations, lateral_col = 'Local_X', moment_col = 'Frame_ID', moments = np.linspace(1, 90, 5), bins = 30, bgein_interval = [-.2,.2], convert2sec = .1):
        """
        Callback method
            
            #-------------------------------------------------------------------------
            Ys,Zs = sdelat.sde_util.GenerateSimulationPaths_tanh_noise_erf_Y()
            #----------
            dis = sdelat.sde_util.marginal_from_dataframe(Zs, bins = 20, moments = [1,5, 10, 20 ,30])
            #----------
            fig,ax = plt.subplots()
            for k in dis.keys():
                hist,edges = dis[k]
                ax.plot(edges, hist, label = str(k))
            ax.legend()
        
        ----------------------------------------------------
        @input: trs_dict_vids
            
            trs_dict_vids[vid] is a data frame 
            
            THe moment columns shoule be in second. 
            
            
        @input: 
        
        """
        filtered_observations  = {}
        
        for vid in observations.keys():
            
            #
            if not min(bgein_interval)<=observations[vid].iloc[0]<=max(bgein_interval):continue
            
            filtered_observations[vid] = copy.deepcopy(observations[vid])
        return filtered_observations
        
        
        pass

    @classmethod
    def marginal_from_dict_of_observations_limitbegin(self, observations, lateral_col = 'Local_X', moment_col = 'Frame_ID', moments = np.linspace(1, 90, 5), bins = 30, bgein_interval = [-.2,.2], convert2sec = .1):
        """
        Callback method
            
            #-------------------------------------------------------------------------
            Ys,Zs = sdelat.sde_util.GenerateSimulationPaths_tanh_noise_erf_Y()
            #----------
            dis = sdelat.sde_util.marginal_from_dataframe(Zs, bins = 20, moments = [1,5, 10, 20 ,30])
            #----------
            fig,ax = plt.subplots()
            for k in dis.keys():
                hist,edges = dis[k]
                ax.plot(edges, hist, label = str(k))
            ax.legend()
        
        ----------------------------------------------------
        @input: trs_dict_vids
            
            trs_dict_vids[vid] is a data frame 
            
            THe moment columns shoule be in second. 
            
            
        @input: 
        
        """
        distributions = {}
        #ymax = data_dp.max().max()
        #ymin = data_dp.min().min()
        
        #
        samples = {moment:[] for moment in moments}
        for moment in moments:
            
            for vid in observations.keys():
                
                #insert the d2
                #   find the idx of d2 in d2s such that
                #   ts[idx] < moment < ts[idx+1]
                #   idx = np.searchsorted(d2s,d2)-1
                ts = convert2sec*(observations[vid].index)
                if max(ts)<=moment:continue
                #
                if not min(bgein_interval)<=observations[vid].iloc[0]<=max(bgein_interval):continue
                
                #find the idex that time equals moment
                idx = np.searchsorted(ts,moment)-1
                #   compute the N and the CCC
                #
                delta_d = 1.0*observations[vid].iloc[idx+1] - observations[vid].iloc[idx]
                delta_t = 1.0*ts[idx+1] - ts[idx]
                #
                data = observations[vid].iloc[idx] + delta_d/delta_t*(moment-ts[idx])
                
                samples[moment].append(data)
        
        #-----------------------------
        distributions = {}
        for moment in moments:
            hist, edges = np.histogram(samples[moment], bins = bins)
            #distributions[moment] = (hist/(edges[1]-edges[0])/(sum(hist)), edges[1:])
            distributions[int(moment*100)/100.0] = (hist, edges[1:])
        
        return distributions

    @classmethod
    def GetIdxWithin(self, data, loc = 1.75):
        """
        find the index that satisfy:
        
            data[idx]>=loc and data[idx+1]<loc
        
        """
        data = np.array(data)
        #find the idx that data_dp[col].iloc[idx]>=loc and data_dp[col].iloc[idx+1]<loc
        #note that idxs1 is obtaied by plus 1
        idxs1 = np.where(data<=loc)[0]
        idxs2 = np.where(data>loc)[0]-1
        
        #find the common idx
        commonidxs = np.intersect1d(idxs1, idxs2)
        
        if len(commonidxs)>0:
            return max(commonidxs)
        else:
            return False


    
    @classmethod
    def MoveTR_2Origin(self, tr_series, loc = 1.75):
        """
        Move the intersection point by the tr and the loc to the origin. 
        
        @input: tr_series
            
            the lateral trajectory. the index is the moments. 
        @output: successornot,new_tr
            
            successornot is a bool. 
            
            
        """
        
        #find the momentt when the trajectory reach loc
        idx = self.GetIdxWithin(tr_series.values, loc = loc)
        if idx==False:
            return False,False
        
        #
        intersected_t = tr_series.index[idx]+(tr_series.index[idx+1]-tr_series.index[idx])*(loc-tr_series.values[idx])/(tr_series.values[idx+1]-tr_series.values[idx])
        tmp = tr_series.index - intersected_t
        #tmp = tr_series.index-tr_series.index[idx]
        #print(idx)
        #newtr = pd.Series(index = tmp)
        #
        #newtr.values = tr_series.values-loc
        newtr = pd.Series(tr_series.values-loc, index = tmp)
        
        return True,newtr

    @classmethod
    def MoveTR_2Origin_omit_negative_moments(self, tr_series, loc = 1.75):
        """
        Move the intersection point by the tr and the loc to the origin. 
        
        @input: tr_series
            
            the lateral trajectory. the index is the moments. 
        @output: successornot,new_tr
            
            successornot is a bool. 
            
            
        """
        
        #find the momentt when the trajectory reach loc
        idx = self.GetIdxWithin(tr_series.values, loc = loc)
        if idx==False:
            return False,False
        
        #
        intersected_t = tr_series.index[idx]+(tr_series.index[idx+1]-tr_series.index[idx])*(loc-tr_series.values[idx])/(tr_series.values[idx+1]-tr_series.values[idx])
        tmp = tr_series.index - intersected_t
        #print(tmp.max())
        #tmp = tr_series.index-tr_series.index[idx]
        #print(idx)
        #newtr = pd.Series(index = tmp)
        #
        #newtr.values = tr_series.values-loc
        newtr = pd.Series(tr_series.values-loc, index = tmp)
        #print(idx)
        return True,newtr.iloc[idx:]


    @classmethod
    def marginal_from_dict_of_dataframe_limitbegin(self, trs_dict_vids, lateral_col = 'Local_X', moment_col = 'Frame_ID', moments = np.linspace(1, 90, 5), bins = 30, bgein_interval = [-.2,.2]):
        """
        Callback method
            
            #-------------------------------------------------------------------------
            Ys,Zs = sdelat.sde_util.GenerateSimulationPaths_tanh_noise_erf_Y()
            #----------
            dis = sdelat.sde_util.marginal_from_dataframe(Zs, bins = 20, moments = [1,5, 10, 20 ,30])
            #----------
            fig,ax = plt.subplots()
            for k in dis.keys():
                hist,edges = dis[k]
                ax.plot(edges, hist, label = str(k))
            ax.legend()
        
        ----------------------------------------------------
        @input: trs_dict_vids
            
            trs_dict_vids[vid] is a data frame 
            
            THe moment columns shoule be in second. 
            
            
        @input: 
        
        """
        distributions = {}
        #ymax = data_dp.max().max()
        #ymin = data_dp.min().min()
        
        #
        samples = {moment:[] for moment in moments}
        for moment in moments:
            for vid in trs_dict_vids.keys():
                vdata = trs_dict_vids[vid]
                #insert the d2
                #   find the idx of d2 in d2s such that
                #   ts[idx] < moment < ts[idx+1]
                #   idx = np.searchsorted(d2s,d2)-1
                ts = vdata[moment_col].values
                if max(ts)<=moment:continue
                if not min(bgein_interval)<=vdata[lateral_col].iloc[0]<=max(bgein_interval):continue
                
                idx = np.searchsorted(ts,moment)-1
                #   compute the N and the CCC
                #
                delta_d = 1.0*vdata[lateral_col].iloc[idx+1] - vdata[lateral_col].iloc[idx]
                delta_t = 1.0*ts[idx+1] - ts[idx]
                #
                data = vdata[lateral_col].iloc[idx] + delta_d/delta_t*(moment-ts[idx])
                
                samples[moment].append(data)
        
        #-----------------------------
        distributions = {}
        for moment in moments:
            hist, edges = np.histogram(samples[moment], bins = bins)
            #distributions[moment] = (hist/(edges[1]-edges[0])/(sum(hist)), edges[1:])
            distributions[moment] = (hist, edges[1:])
        
        return distributions

    @classmethod
    def marginal_distribution_times_given_dp(self,data_dp, moments = np.linspace(1, 90, 5), bins = 30):
        """
        Get the distributions along certain times. 
        
        ----------------------------
        @output: distributions
            a dict 
            
            distributions[moment] = (edges, hists), two lists with the same length. 
            sum(hist) should be 1. 
            
            
        """
        distributions = {}
        #ymax = data_dp.max().max()
        #ymin = data_dp.min().min()
        
        #
        ts = data_dp.index
        
        for moment in moments:
            
            #insert the d2
            #   find the idx of d2 in d2s such that
            #   ts[idx] < moment < ts[idx+1]
            #   idx = np.searchsorted(d2s,d2)-1
            idx = np.searchsorted(ts,moment)-1
            #   compute the N and the CCC
            delta_d = 1.0*data_dp.iloc[idx+1,:] - data_dp.iloc[idx,:]
            delta_t = 1.0*ts[idx+1] - ts[idx]
            data = data_dp.iloc[idx,:] + delta_d/delta_t*(moment-ts[idx])
            
            hist, edges = np.histogram(data, bins = bins)
            
            #distributions[moment] = (hist/(edges[1]-edges[0])/(sum(hist)), edges[1:])
            distributions[moment] = (hist, edges[1:])
            
        return distributions
        
        
        
    @classmethod
    def marginal_from_dict_of_dataframe(self, trs_dict_vids, lateral_col = 'Local_X', moment_col = 'Frame_ID', moments = np.linspace(1, 90, 5), bins = 30):
        """
        Callback method
            
            #-------------------------------------------------------------------------
            Ys,Zs = sdelat.sde_util.GenerateSimulationPaths_tanh_noise_erf_Y()
            #----------
            dis = sdelat.sde_util.marginal_from_dataframe(Zs, bins = 20, moments = [1,5, 10, 20 ,30])
            #----------
            fig,ax = plt.subplots()
            for k in dis.keys():
                hist,edges = dis[k]
                ax.plot(edges, hist, label = str(k))
            ax.legend()
        
        ----------------------------------------------------
        @input: trs_dict_vids
            
            trs_dict_vids[vid] is a data frame 
            
            THe moment columns shoule be in second. 
            
            
        @input: 
        
        """
        distributions = {}
        #ymax = data_dp.max().max()
        #ymin = data_dp.min().min()
        
        #
        samples = {moment:[] for moment in moments}
        for moment in moments:
            for vid in trs_dict_vids.keys():
                vdata = trs_dict_vids[vid]
                #insert the d2
                #   find the idx of d2 in d2s such that
                #   ts[idx] < moment < ts[idx+1]
                #   idx = np.searchsorted(d2s,d2)-1
                ts = vdata[moment_col].values
                if max(ts)<=moment:continue
                idx = np.searchsorted(ts,moment)-1
                #   compute the N and the CCC
                #
                delta_d = 1.0*vdata[lateral_col].iloc[idx+1] - vdata[lateral_col].iloc[idx]
                delta_t = 1.0*ts[idx+1] - ts[idx]
                #
                data = vdata[lateral_col].iloc[idx] + delta_d/delta_t*(moment-ts[idx])
                
                samples[moment].append(data)
        
        #-----------------------------
        distributions = {}
        for moment in moments:
            hist, edges = np.histogram(samples[moment], bins = bins)
            #distributions[moment] = (hist/(edges[1]-edges[0])/(sum(hist)), edges[1:])
            distributions[moment] = (hist, edges[1:])
        
        return distributions
    
    @classmethod
    def plot_lc_changing_with_observation(self, observation, lc_changings, lane_marks_lateral_coors = [] , alpha = .3):
        """
        
        plot the lane changing identification results. 
        
        
        
        -------------------------------------------
        @input: observation
            
            pd.Series. An observation. 
            
            observation.index is the moments. 
        
        @ininput: lc_changings
        
            a dict. key are the idx in observation that is identified as lane changing and value is either 1 or -1, indicate the lane changing to positive or negative direction. 
            
            lc_changings[idx] = rho_after_changing. 
            
            Note that there is always exist lc_changings[0]. 
            
        @input: lane_marks_lateral_coors
            
            the coordinates of the lane marks. 
        
        """
        import matplotlib.patches as patches

        fig,ax = plt.subplots()
        #
        ax.plot(range(len(observation)), observation.values)
        
        
        #
        if len(lc_changings)==1:return
        
        for idx1,idx2 in zip(sorted(lc_changings.keys())[:-1], sorted(lc_changings.keys())[1:]):
            if lc_changings[idx1]==1:
                color = np.array([1.0, 0,0])
                #color = np.random.uniform(size  = (3,))
            elif lc_changings[idx1]==-1:
                color = np.array([0, 0, 1])
            else:
                color = np.array([0, 1, 0])
            #
            xy = idx1,min(observation)
            width = idx2 - idx1
            height = max(observation) - min(observation)
            
            #
            ax.add_patch(patches.Rectangle(xy, width, height, alpha= alpha, color = color))
            
            pass
        
        return ax
        
        
        
        pass

    
    @classmethod
    def distribution_reach_coor_given_observations(self, observations, begin = 0, locs = [ 1.75, 3.5], bins = 30, convert2sec = 1.0):
        """
        obtain the distribution when trajectory reach certain locs for the first time.  
        each column is a sample path
        
        
        --------------------------------------------
        @input: observations    
            
            a dict. the key is the vehicle id. 
            
            Each value is a series. the index is the moment, usually in second. 
        
        
        
        """
        #distributions[loc] = (hist, edges)
        distributions = {}
        
        
        for loc in locs:
            #
            ts = []
            for vid in observations.keys():
                observation = observations[vid]
                moments = (observation.index.values - min(observation.index.values))*convert2sec
                data = observation.values
                
                #----find the moment that begin
                #find the idx that data_dp[col].iloc[idx]>=loc and data_dp[col].iloc[idx+1]<loc
                #note that idxs2 is obtaied by plus 1
                idxs1 = np.where(data<=begin)[0]
                idxs2 = np.where(data>begin)[0]-1
                #find the common idx
                commonidxs_begin = np.intersect1d(idxs1, idxs2)
                
                
                #find the idx that data_dp[col].iloc[idx]>=loc and data_dp[col].iloc[idx+1]<loc
                #note that idxs2 is obtaied by plus 1
                idxs1 = np.where(data<=loc)[0]
                idxs2 = np.where(data>loc)[0]-1
                #
                #find the common idx
                commonidxs_loc = np.intersect1d(idxs1, idxs2)
                
                if len(commonidxs_loc)>0 and len(commonidxs_begin)>0:
                    #print(commonidxs)
                    ts.append(-moments[max(commonidxs_begin)] + moments[min(commonidxs_loc)])
            
            hist, edges = np.histogram(ts, bins = bins)
            distributions[loc] = (hist, edges[1:])
        
        return distributions
    
    
    
    
    @classmethod
    def distribution_samples_reach_coor_given_dp(self, data_dp, locs = [ 1.75, 3.5], bins = 30):
        """
        obtain the distribution when trajectory reach certain locs for the first time.  
        each column is a sample path
        
        
        --------------------------------------------
        
        """
        #distributions[loc] = (hist, edges)
        samples = {}
        
        
        for loc in locs:
            #
            ts = []
            for col in data_dp.columns:
                
                #find the idx that data_dp[col].iloc[idx]>=loc and data_dp[col].iloc[idx+1]<loc
                #note that idxs2 is obtaied by plus 1
                idxs1 = np.where(data_dp[col].values<=loc)[0]
                idxs2 = np.where(data_dp[col].values>loc)[0]-1
                
                #find the common idx
                commonidxs = np.intersect1d(idxs1, idxs2)
                
                if len(commonidxs)>0:
                    #print(commonidxs)
                    ts.append(data_dp.index[min(commonidxs)])
            
            samples[loc] = ts
        
        return samples
    
    
    @classmethod
    def distribution_reach_coor_given_dp(self, data_dp, locs = [ 1.75, 3.5], bins = 30):
        """
        obtain the distribution when trajectory reach certain locs for the first time.  
        each column is a sample path
        
        
        --------------------------------------------
        
        """
        #distributions[loc] = (hist, edges)
        distributions = {}
        
        
        for loc in locs:
            #
            ts = []
            for col in data_dp.columns:
                
                #find the idx that data_dp[col].iloc[idx]>=loc and data_dp[col].iloc[idx+1]<loc
                #note that idxs2 is obtaied by plus 1
                idxs1 = np.where(data_dp[col].values<=loc)[0]
                idxs2 = np.where(data_dp[col].values>loc)[0]-1
                
                #find the common idx
                commonidxs = np.intersect1d(idxs1, idxs2)
                
                if len(commonidxs)>0:
                    #print(commonidxs)
                    ts.append(data_dp.index[min(commonidxs)])
            
            hist, edges = np.histogram(ts, bins = bins)
            distributions[loc] = (hist, edges[1:])
        
        return distributions
    

    
    @classmethod
    def LC_durations_from_data(self, data_pd, lw = 3.5):
        """
        
        Callback:
        
            sdelat.dataanalysis.LC_duration_from_data(data_pd, lw = 3.5)
            
        --------------------------------------------------------
        @input: data_pd
        
            the simulation data. 
            Each sample is a trajectory. 
            
            data_pd.shape = (*, N), N is the samples number. 
            
            data_pd.index is the moments. 
            
            Each row correspond to a moment. 
        @OUTPUT:
            lc_durations, a list. 
        """
        samplesize = data_pd.shape[1]
        momentN = data_pd.shape[0]
        
        #each element is a duration for a sample. 
        lc_durations = []
        for column in data_pd.columns:
            #a sample is a series
            sample = data_pd[column]
            
            #find the index that traverse the lw. 
            #tmp = (array_of_bool, )
            tmp =  np.where((np.array(sample.iloc[:-1])<=lw) & (np.array(sample.iloc[1:])>=lw))[0]
            valid1 = np.array(range(1, momentN))[tmp]
            tmp =  np.where((np.array(sample.iloc[:-1])>=lw) & (np.array(sample.iloc[1:])<=lw))[0]
            valid2 = np.array(range(1, momentN))[tmp]
            
            #index of the moment when traverse the lw. 
            idx = np.inf
            if len(valid1)>0:
                idx = min(idx, valid1[0])
            if len(valid2)>0:
                idx = min(idx, valid2[0])
            
            if idx==np.inf:continue
            
            #print(idx)
            
            lc_durations.append(data_pd.index[idx])
        
        return lc_durations
    
    
    pass



class sde_util():
    """
    All utilities of the sde methods. 
    
    
    
    """
    
    @classmethod
    def LCIMs_using_lateral_zerospeex(self, observations, convert2sec = .04, tspan = np.linspace(0.0, 20.0, 100), interval_deviation = 0, target_lats_inference = [1.75, 3.0], lw = 3.5):
        """
        
        """
        lcims = {}
        lcim_laterals = {}
        leading_ts = {}
        
        #---------------------------------
        for vid in observations.keys():
            t_index,leading_t  =  self.LCIM_using_lateral_zerospeex(observation = observations[vid])
            
            #
            if t_index==False:
                continue
                
            #
            lcims[vid] = observations[vid].index[t_index]
            lcim_laterals[vid] = observations[vid].iloc[t_index]
            leading_ts[vid] = leading_t
            
        return lcims,lcim_laterals,leading_ts
        #------------------------------
    
    
    @classmethod
    def LCIM_using_lateral_zerospeex(self, observation, convert2sec = .04, tspan = np.linspace(0.0, 20.0, 100), interval_deviation = 0, target_lats_inference = [1.75, 3.0], lw = 3.5):
        """
        
        -------------------------------------------------
        @input: observation
            
            a pandas.Series. 
            
            The index are moment and the values are lateral displacement. 
        
        
        
        @input: interval_deviation_from_middleline
        
            a scalar that describe the 
        
            
        
        
        
        @input: batch_res
        
            the inference results. a dict. 
            
            batch_res[vehicle_id] is a dict. keys inlcude:
            
                dict_keys(['k', 'sigma', 'k_path', 'sigma_path', 'changing_moments_and_rhos', 'lc_onset_lat', 'y0', \
                #'clmm_values', 'lcd_values', 'true_durations', 1.8, 3.6, 'success'])
        
        @input: observations_dict
        
            a dict. 
            
            observations_dict[vehicle_id] is a pd.Series. 
        """
        #
        lateral = np.diff(observation.values)
        
        #many possibliltyes 
        #tspan = np.where(lateral==0)[0]
        #if len(tspan)>0:
        #    return tspan
        
        zerospeed_ts = []
        for i,speed in enumerate(lateral[:-1]):
            if lateral[i]==0:
                zerospeed_ts.append(i)
                continue
            elif lateral[i]<0 and lateral[i]>0:
                zerospeed_ts.append(i)
                continue
            elif lateral[i]>0 and lateral[i]<0:
                zerospeed_ts.append(i)
                continue
            
        if len(zerospeed_ts)==0:
            return False,False
        
        #
        zerospeed_ts = np.array(zerospeed_ts)
        #Find the moment when observation traverse lw/2
        #------------------------
        tmp = observation.values - lw/2.0
        min_tmp = min(observation.values - lw/2)
        index_traverse_lane_mark = np.argmin(tmp)
        
        #------------------------
        #print(zerospeed_ts, index_traverse_lane_mark)
        filtered = zerospeed_ts[np.array(zerospeed_ts)<index_traverse_lane_mark]
        
        if len(filtered)==0:
            return False,False
        
        leading_t = observation.index[index_traverse_lane_mark]-observation.index[filtered[-1]]
        
        return filtered[-1],leading_t
        
        
        
        
        
        
        
        
        return zerospeed_ts
        
        #
        tspan_idx_negative = set(np.where(lateral<0)[0])
        tspan_idx_positive = set(np.where(lateral>0)[0])
        
        if len(tspan_idx_negative)==0 or len(tspan_idx_positive)==0:
            return False
        
        
        #find the intersection. 
        transitions0 = set(tspan_idx_negative+1).intersection(tspan_idx_positive)
        transitions1 = set(tspan_idx_positive+1).intersection(tspan_idx_negative)
        
        return list(transitions0)+list(transitions1)
        
    
    @classmethod
    def mean_variance_given_dicts_of_observations(self, observations, moments = np.linspace(1, 60, 300), convert2sec = 1.0, bgein_interval = [-0.1, 0.1]):
        """
        
        
        -------------------------------------------------------------
        @input: data_dp
            
            a pd.Dataframe. The index are the moments and the columns are the simulated paths. 
            
        -------------------------------------------------------------
        
        """
        
        #
        


        #
        samples = {moment:[] for moment in moments}
        for moment in moments:
            
            for vid in observations.keys():
                
                #insert the d2
                #   find the idx of d2 in d2s such that
                #   ts[idx] < moment < ts[idx+1]
                #   idx = np.searchsorted(d2s,d2)-1
                ts = convert2sec*(observations[vid].index)
                if max(ts)<=moment:continue
                #
                if not min(bgein_interval)<=observations[vid].iloc[0]<=max(bgein_interval):continue
                
                #find the idex that time equals moment
                idx = np.searchsorted(ts,moment)-1
                #   compute the N and the CCC
                #
                delta_d = 1.0*observations[vid].iloc[idx+1] - observations[vid].iloc[idx]
                delta_t = 1.0*ts[idx+1] - ts[idx]
                #
                data = observations[vid].iloc[idx] + delta_d/delta_t*(moment-ts[idx])
                
                samples[moment].append(data)
            
        #
        means = []
        stds = []
        varss = []
        for moment in moments:
            if len(samples[moment])>0:
                means.append(np.mean(samples[moment]))
                stds.append(np.std(samples[moment]))
                varss.append(np.var(samples[moment]))
            else:
                means.append(np.nan)
                stds.append(np.nan)
                varss.append(np.nan)
            
            
        return {'means':means, 'stds':stds, 'vars':varss}
    
    
    @classmethod
    def mean_variance_given_simulated_paths(self, data_dp, moments = np.linspace(1, 60, 300)):
        """
        
        
        -------------------------------------------------------------
        @input: data_dp
            
            a pd.Dataframe. The index are the moments and the columns are the simulated paths. 
            
        -------------------------------------------------------------
        
        """
        
        #
        ts = data_dp.index
        
        means = []
        stds = []
        varss = []
        for moment in moments:
            
            #insert the d2
            #   find the idx of d2 in d2s such that
            #   ts[idx] < moment < ts[idx+1]
            #   idx = np.searchsorted(d2s,d2)-1
            idx = np.searchsorted(ts,moment)-1
            #   compute the N and the CCC
            delta_d = 1.0*data_dp.iloc[idx+1,:] - data_dp.iloc[idx,:]
            delta_t = 1.0*ts[idx+1] - ts[idx]
            data = data_dp.iloc[idx,:] + delta_d/delta_t*(moment-ts[idx])
            #print(data.shape)
            means.append(np.mean(data.values))
            stds.append(np.std(data.values))
            varss.append(np.var(data.values))
            
        return {'means':means, 'stds':stds, 'vars':varss}
    
    
    @classmethod
    def logprobs2meanlogprob(self, logprobs0):
        """
        calculate the log of mean probability from a series of logprobs.
        
        """
        logprobs = sorted(logprobs0)
        
        digitsMultipled = np.log10(-max(logprobs))
        #sample_logs_added = [i for i in 10**(digitsMultipled)+sample_logs]
        sample_logs_added = -logprobs[-1]+logprobs
        
        
        return -10**digitsMultipled+np.log(np.mean(np.exp(sample_logs_added)))
        
    
    @classmethod
    def distribution_reach_locs(self, data_dp, locs = [1, 1.7, 3.5], bins = 30):
        """
        The distributions when the system path traverse certain locs.
        
        """
        #distributions[loc] = (hist, edges)
        distributions = {}
        #
        for loc in locs:
            #
            ts = []
            for col in data_dp.columns:
                
                #find the idx that data_dp[col].iloc[idx]>=loc and data_dp[col].iloc[idx+1]<loc
                #note that idxs2 is obtaied by plus 1
                idxs1 = np.where(data_dp[col].values<=loc)[0]
                idxs2 = np.where(data_dp[col].values>loc)[0]-1
                
                #find the common idx
                commonidxs = np.intersect1d(idxs1, idxs2)
                
                if len(commonidxs)>0:
                    #print(commonidxs)
                    ts.append(data_dp.index[min(commonidxs)])
            
            hist, edges = np.histogram(ts, bins = bins)
            distributions[loc] = (hist, edges[1:])
        
        return distributions
    
    @classmethod
    def marginal_from_dataframe(self, data_dp, moments = np.linspace(1, 90, 5), bins = 30):
        """
        Callback method
            
            #-------------------------------------------------------------------------
            Ys,Zs = sdelat.sde_util.GenerateSimulationPaths_tanh_noise_erf_Y()
            #----------
            dis = sdelat.sde_util.marginal_from_dataframe(Zs, bins = 20, moments = [1,5, 10, 20 ,30])
            #----------
            fig,ax = plt.subplots()
            for k in dis.keys():
                hist,edges = dis[k]
                ax.plot(edges, hist, label = str(k))
            ax.legend()
        
        ----------------------------------------------------
        @input: data_dp
            
            a dataframe which represent the system path. 
            
            index are the moments, each columns is a system path. 
            
            
        @input: 
        
        """
        distributions = {}
        #ymax = data_dp.max().max()
        #ymin = data_dp.min().min()
        
        #
        ts = data_dp.index
        
        for moment in moments:
            
            moment = int(moment*100)/100.0
            
            #insert the d2
            #   find the idx of d2 in d2s such that
            #   ts[idx] < moment < ts[idx+1]
            #   idx = np.searchsorted(d2s,d2)-1
            idx = np.searchsorted(ts,moment)-1
            #   compute the N and the CCC
            delta_d = 1.0*data_dp.iloc[idx+1,:] - data_dp.iloc[idx,:]
            delta_t = 1.0*ts[idx+1] - ts[idx]
            data = data_dp.iloc[idx,:] + delta_d/delta_t*(moment-ts[idx])
            
            hist, edges = np.histogram(data, bins = bins)
            
            #distributions[moment] = (hist/(edges[1]-edges[0])/(sum(hist)), edges[1:])
            distributions[moment] = (hist, edges[1:])
            
        return distributions
        
        
        
    @classmethod
    def plot_dataframe_from_simulationresult(self, data_dp, N_plotted= 100, figsize = (5,3), ax = False):
        """
        plot the simulated data. 
        --------------------------------
        @input: data_dp
            
            a data frame. index are the moments and each column is a system path. 
        
        --------------------------------
        
        """
        if isinstance(ax, bool):
            fig,ax = plt.subplots(figsize = figsize)
        
        #
        selected_columns = np.random.choice(data_dp.columns,N_plotted)
        
        for col in selected_columns:
            ax.plot(data_dp.index, data_dp[col])
        
        
        
        return ax
        
        
        
    @classmethod
    def plot_distribution_reach_locs(self, dis):
        """
        
        -----------------------------------------------------------
        
        """
        
        
        
        
        pass
    
    
    
    @classmethod
    def MonteCarloSimulationConfidenceInterval(self, alpha = 4.0,  N_paths  = 100, Z0 = np.array([.01, .01]), tspan = np.linspace(0.0, 90, 100), sde_args = {'sigma':.03, 'thetaY':.1, 'thetaZ':.01, 'lw':3.5, 'target_Y':3.5, 'v':1}, invervals = [.95, .9, .85]):
        """
        Generate the mean and the confidence interval of the trajectories. 
        
        -----------------------------------------------------------
        @input: Z0 = np.array([.01, .01])
            plot_distributions_dict
            Z0 = [lateral_inital, noise]
        
        
        @output: paths
            
            a pd. 
            
        @output: Ys,Zs
        
            
            Both are pd.dataframe. The index are moments, and each column is a system path. 
        
        
        -------------------------------------------------------------------------------
        
        """
        
        
        
        
        pass
    
    
    @classmethod
    def MonteCarloSimulationPaths_tanh_noise_erf_Y_test(self, alpha = 4.0,  N_paths  = 100, Z0 = np.array([.01, .01]), tspan = np.linspace(0.0, 90, 100), sde_args = {'sigma':.03, 'thetaY':.1, 'thetaZ':.01, 'lw':3.5, 'target_Y':3.5, 'v':1}):
        """
        Generate the simulation path for the erf based noise and tanh type noise and erf based regression.
        
        
        -----------------------------------------------------------
        @input: Z0 = np.array([.01, .01])
            plot_distributions_dict
            Z0 = [lateral_inital, noise]
        
        
        @output: paths
            
            a pd. 
            
        @output: Ys,Zs
        
            
            Both are pd.dataframe. The index are moments, and each column is a system path. 
            
        """
        import sdeint

        sigma = sde_args['sigma']
        thetaY = sde_args['thetaY']
        thetaZ = sde_args['thetaZ']
        lw = sde_args['lw']
        target_Y = sde_args['target_Y']
        v = sde_args['v']

        #
        res_Y = []
        res_Z = []
        
        #--------------------------------------
        #erf implemetation
        def exp_in_errorfun(Y,lw, alpha = alpha):
            return np.exp(-(Y+lw/2.0)**2/((lw/alpha)**2))
        #---------------------------------
        def exp_in_errorfun_im(Y,lw):
            return exp_in_errorfun(Y,lw)-exp_in_errorfun(Y,-lw)
        
        #-----------------------------------------------
        #Z = (Z[0], Z[1]), Z[0] is the lateral displacement and Z[1] is the nosie. 
        def f_Z(Z,t, sigma = sigma, lw = lw,thetaY = thetaY, target_Y = target_Y, thetaZ = thetaZ, v = v):
            #normalize2_1(Z[0]-target_Y,M)*
            tmp0 = thetaY*(exp_in_errorfun_im(Z[0]- target_Y,lw))- thetaZ*(1- Z[1]*Z[1]*4/lw/lw)*Z[1]/np.sqrt((1e-10+Z[1]*Z[1]))
            tmp1 = -2.0*Z[1]*(v**2)*(sigma**2)*(lw**2-4.0*(Z[1]**2))/2/(lw**2)
            #tmp1 = -.5*2*lanewidth/np.pi*sigma**2*np.tan(Z[1]*np.pi/lanewidth)/((np.tan(Z[1]*np.pi/lanewidth)**2+1)**2)
            #return np.array([[tmp0, 0], [0, tmp1]])
            return np.array([tmp0, tmp1])
        def g_Z(Z,t, lw = lw, sigma = sigma):
            tmp0 = 0
            tmp1 = lw*v*sigma*(lw**2 - 4.0*( Z[1]**2))/2/(lw**2)
            #tmp1 = 1.0*lanewidth/np.pi*sigma/(np.tan(Z[1]*np.pi/lanewidth)**2+1)
            return np.array([[0, 0], [0, tmp1]])
            
        #---------------------------------------------
        for iterr in range(N_paths):
            #result.shape is (len(tspan),2), '2' corresponds to Z.
            #   the first column is lateral displacement, and the 2nd columns is the noize. 
            result = sdeint.itoint(f_Z, g_Z, Z0, tspan)
            res_Y.append(result[:,0])
            res_Z.append(result[:,1])
        #-----------------------------------------
        #np.array(res_Y) shape is N_paths,len(tspan)
        #np.array(res_Z) shape is N_paths,len(tspan)
        return pd.DataFrame(res_Y, columns =  tspan).T,pd.DataFrame(res_Z, columns =  tspan).T
    
    @classmethod
    def LCD_based_on_lateral_zero_speed(self, Ys, limit  = 3.5/2):
        """
        @input: Ys
            a pd.DataFrame. 
            
            Each row is a moment and each columns is a sample path. 
            
            Ys.index gives the simulated moments and Ys.columns are the index of the simulated path. 
        
        
        """
        samples = []
        #
        for i in range(Ys.shape[1]):
            #single simulated path. 
            path = Ys.iloc[:,i]
            #
            #Find the moment when the lateral speed is zero.
            #   it is determined by speeds
            speeds = np.diff(path)/(np.diff(path.index))
            
            #find the index that the speed change sign.
            productwithsign  = speeds[:-1]*speeds[1:]
            #
            idxs = np.where(productwithsign<=0)[0]
            if len(idxs)==0:continue
            
            if path.iloc[idxs[0]+2]<=limit:continue
            
            #
            idx = idxs[0]
            samples.append((path.index[idx+1]+path.index[idx+2])/2.0)
        
        return samples

    @classmethod
    def MonteCarloSimulationPaths_tanh_noise_erf_Y_assymetric(self, alpha_l = 4.0, alpha_r = 3.0, N_paths  = 100, Z0 = np.array([.0, .0]), tspan = np.linspace(0.0, 90, 100), sde_args = {'sigma':.03, 'thetaY':.1, 'thetaZ':.01, 'lw':3.5, 'target_Y':3.5, 'v':1}):
        """
        Generate the simulation path for the erf based noise and tanh type noise and erf based regression.
        
        Difference:
        
            - MonteCarloSimulationPaths_tanh_noise_erf_Y_assymetric(), the doerf function is assymetric. which is represented by alpha_l and alpha_r, while there is only a alpha in the following formula. 
            - MonteCarloSimulationPaths_tanh_noise_erf_Y()
        -----------------------------------------------------------
        @input: Z0 = np.array([.01, .01])
            
            Z0 = [lateral_inital, noise]
        
        
        @output: paths
            
            a pd. 
            
        @output: Ys,Zs
        
            
            Both are pd.dataframe. The index are moments, and each column is a system path. 
            
        """
        import sdeint

        sigma = sde_args['sigma']
        thetaY = sde_args['thetaY']
        thetaZ = sde_args['thetaZ']
        lw = sde_args['lw']
        target_Y = sde_args['target_Y']
        v = sde_args['v']

        #
        res_Y = []
        res_Z = []
        
        #--------------------------------------
        #erf implemetation
        def exp_in_errorfun(Y,lw, alpha):
            return np.exp(-(Y+lw/2.0)**2/((lw/alpha)**2))
        #---------------------------------
        def exp_in_errorfun_im(Y,lw, alpha_l , alpha_r):
            return exp_in_errorfun(Y=Y,lw=lw, alpha = alpha_r)-exp_in_errorfun(Y=Y,lw = -lw, alpha = alpha_l)
        
        #-----------------------------------------------
        #Z = (Z[0], Z[1]), Z[0] is the lateral displacement and Z[1] is the nosie. 
        def f_Z(Z,t, sigma = sigma, lw = lw,thetaY = thetaY, target_Y = target_Y, thetaZ = thetaZ, v = v, alpha_l  = alpha_l, alpha_r = alpha_r):
            #normalize2_1(Z[0]-target_Y,M)*
            tmp0 = thetaY*(exp_in_errorfun_im(Y = Z[0] - target_Y,lw = lw, alpha_l = alpha_l, alpha_r = alpha_r)) - thetaZ* Z[1]
            tmp1 = -2.0*Z[1]*(v**2)*(sigma**2)*(lw**2-4.0*(Z[1]**2))/2/(lw**2)
            #tmp1 = -.5*2*lanewidth/np.pi*sigma**2*np.tan(Z[1]*np.pi/lanewidth)/((np.tan(Z[1]*np.pi/lanewidth)**2+1)**2)
            #return np.array([[tmp0, 0], [0, tmp1]])
            return np.array([tmp0, tmp1])
        def g_Z(Z,t, lw = lw, sigma = sigma):
            tmp0 = 0
            tmp1 = lw*v*sigma*(lw**2 - 4.0*( Z[1]**2))/2/(lw**2)
            #tmp1 = 1.0*lanewidth/np.pi*sigma/(np.tan(Z[1]*np.pi/lanewidth)**2+1)
            return np.array([[0, 0], [0, tmp1]])
            
        #---------------------------------------------
        for iterr in range(N_paths):
            #result.shape is (len(tspan),2), '2' corresponds to Z.
            #   the first column is lateral displacement, and the 2nd columns is the noize. 
            result = sdeint.itoint(f_Z, g_Z, Z0, tspan)
            res_Y.append(result[:,0])
            res_Z.append(result[:,1])
        #-----------------------------------------
        #np.array(res_Y) shape is N_paths,len(tspan)
        #np.array(res_Z) shape is N_paths,len(tspan)
        return pd.DataFrame(res_Y, columns =  tspan).T,pd.DataFrame(res_Z, columns =  tspan).T



    @classmethod
    def MonteCarloSimulationPaths_tanh_noise_erf_Y_with_lc_failiure(self, alpha = 4.0,  N_paths  = 100, Z0 = np.array([.0, .0]), tspan = np.linspace(0.0, 90, 100), sde_args = {'sigma':.03, 'thetaY':.1, 'thetaZ':.01, 'lw':3.5, 'target_Y':3.5, 'v':1}, return_lat_loc = 1.2, tolerance = 1e-2):
        """
        Generate the simulation path for the erf based noise and tanh type noise and erf based regression.
        
        
        -----------------------------------------------------------
        @input: Z0 = np.array([.01, .01])
            
            Z0 = [lateral_inital, noise]
        
        @input: return_lat_loc
        
            when the vehicle reach this lat loc, then return. 
        
        @output: paths
            
            a pd. 
            
        @output: Ys,Zs
        
            
            Both are pd.dataframe. The index are moments, and each column is a system path. 
            
        """
        import sdeint

        sigma = sde_args['sigma']
        thetaY = sde_args['thetaY']
        thetaZ = sde_args['thetaZ']
        lw = sde_args['lw']
        target_Y = sde_args['target_Y']
        v = sde_args['v']

        #
        res_Y = []
        res_Z = []
        
        #--------------------------------------
        #erf implemetation
        def exp_in_errorfun(Y,lw, alpha = alpha):
            return np.exp(-(Y+lw/2.0)**2/((lw/alpha)**2))
        #---------------------------------
        def exp_in_errorfun_im(Y,lw):
            return exp_in_errorfun(Y,lw)-exp_in_errorfun(Y,-lw)
        
        #-----------------------------------------------
        returned = False
        #Z = (Z[0], Z[1]), Z[0] is the lateral displacement and Z[1] is the nosie. 
        def f_Z(Z,t, sigma = sigma, lw = lw,thetaY = thetaY, target_Y = target_Y, thetaZ = thetaZ, v = v, Z0 = Z0):
            #normalize2_1(Z[0]-target_Y,M)*
            if Z[0] > return_lat_loc:
                target_Y1  =  -3.5
            else:
                target_Y1 = target_Y
            #
            tmp0 = thetaY*(exp_in_errorfun_im(Z[0] - target_Y1,lw)) - thetaZ* Z[1]
            tmp1 = -2.0*Z[1]*(v**2)*(sigma**2)*(lw**2-4.0*(Z[1]**2))/2/(lw**2)

            #tmp1 = -.5*2*lanewidth/np.pi*sigma**2*np.tan(Z[1]*np.pi/lanewidth)/((np.tan(Z[1]*np.pi/lanewidth)**2+1)**2)
            #return np.array([[tmp0, 0], [0, tmp1]])
            return np.array([tmp0, tmp1])
        def g_Z(Z,t, lw = lw, sigma = sigma):
            tmp0 = 0
            tmp1 = lw*v*sigma*(lw**2 - 4.0*( Z[1]**2))/2/(lw**2)
            #tmp1 = 1.0*lanewidth/np.pi*sigma/(np.tan(Z[1]*np.pi/lanewidth)**2+1)
            return np.array([[0, 0], [0, tmp1]])
            
        #---------------------------------------------
        for iterr in range(N_paths):
            #result.shape is (len(tspan),2), '2' corresponds to Z.
            #   the first column is lateral displacement, and the 2nd columns is the noize. 
            result = sdeint.itoint(f_Z, g_Z, Z0, tspan)
            res_Y.append(result[:,0])
            res_Z.append(result[:,1])
        #-----------------------------------------
        #np.array(res_Y) shape is N_paths,len(tspan)
        #np.array(res_Z) shape is N_paths,len(tspan)
        return pd.DataFrame(res_Y, columns =  tspan).T,pd.DataFrame(res_Z, columns =  tspan).T





    @classmethod
    def MonteCarloSimulationPaths_tanh_noise_erf_Y(self, alpha = 4.0,  N_paths  = 100, Z0 = np.array([.0, .0]), tspan = np.linspace(0.0, 90, 100), sde_args = {'sigma':.03, 'thetaY':.1, 'thetaZ':.01, 'lw':3.5, 'target_Y':3.5, 'v':1}):
        """
        Generate the simulation path for the erf based noise and tanh type noise and erf based regression.
        
        
        -----------------------------------------------------------
        @input: Z0 = np.array([.01, .01])
            
            Z0 = [lateral_inital, noise]
        
        
        @output: paths
            
            a pd. 
            
        @output: Ys,Zs
        
            
            Both are pd.dataframe. The index are moments, and each column is a system path. 
            
        """
        import sdeint

        sigma = sde_args['sigma']
        thetaY = sde_args['thetaY']
        thetaZ = sde_args['thetaZ']
        lw = sde_args['lw']
        target_Y = sde_args['target_Y']
        v = sde_args['v']

        #
        res_Y = []
        res_Z = []
        
        #--------------------------------------
        #erf implemetation
        def exp_in_errorfun(Y,lw, alpha = alpha):
            return np.exp(-(Y+lw/2.0)**2/((lw/alpha)**2))
        #---------------------------------
        def exp_in_errorfun_im(Y,lw):
            return exp_in_errorfun(Y,lw)-exp_in_errorfun(Y,-lw)
        
        #-----------------------------------------------
        #Z = (Z[0], Z[1]), Z[0] is the lateral displacement and Z[1] is the nosie. 
        def f_Z(Z,t, sigma = sigma, lw = lw,thetaY = thetaY, target_Y = target_Y, thetaZ = thetaZ, v = v):
            #normalize2_1(Z[0]-target_Y,M)*
            tmp0 = thetaY*(exp_in_errorfun_im(Z[0] - target_Y,lw)) - thetaZ* Z[1]
            tmp1 = -2.0*Z[1]*(v**2)*(sigma**2)*(lw**2-4.0*(Z[1]**2))/2/(lw**2)
            #tmp1 = -.5*2*lanewidth/np.pi*sigma**2*np.tan(Z[1]*np.pi/lanewidth)/((np.tan(Z[1]*np.pi/lanewidth)**2+1)**2)
            #return np.array([[tmp0, 0], [0, tmp1]])
            return np.array([tmp0, tmp1])
        def g_Z(Z,t, lw = lw, sigma = sigma):
            tmp0 = 0
            tmp1 = lw*v*sigma*(lw**2 - 4.0*( Z[1]**2))/2/(lw**2)
            #tmp1 = 1.0*lanewidth/np.pi*sigma/(np.tan(Z[1]*np.pi/lanewidth)**2+1)
            return np.array([[0, 0], [0, tmp1]])
            
        #---------------------------------------------
        for iterr in range(N_paths):
            #result.shape is (len(tspan),2), '2' corresponds to Z.
            #   the first column is lateral displacement, and the 2nd columns is the noize. 
            result = sdeint.itoint(f_Z, g_Z, Z0, tspan)
            res_Y.append(result[:,0])
            res_Z.append(result[:,1])
        #-----------------------------------------
        #np.array(res_Y) shape is N_paths,len(tspan)
        #np.array(res_Z) shape is N_paths,len(tspan)
        return pd.DataFrame(res_Y, columns =  tspan).T,pd.DataFrame(res_Z, columns =  tspan).T





    
    @classmethod
    def prob_Y_i_given_previousYandZ(self, Y_i1, Y_i, Z_i_1, deltat = .5, lw = 3.5, args_sde = {'theta':.5, 'sigma':.01},multiplier = 1e7):
        """
        @input: Y_i1, Y_i, Z_i_1
        
            Y_i1 is the Y at moment i+1
            Y_i is the Y at moment i
            Z_i_1 is the Z at moment i-1
        @OUTUPUT: prob
        
        """
        meany,sigmasquarey = self.meanY_sigmaY_given_previousYandZ(Y_i = Y_i, Z_i_1 = Z_i_1, deltat = deltat, lw = lw, args_sde = args_sde)
        
        print(Y_i1, meany,sigmasquarey)
        #print(normal(Y_i1, meany, sigmasquarey))
        return multiplier*normal(Y_i1, meany, sigmasquarey)
    
    
    @classmethod
    def logprob_givenbreakingpoints_singlelc_positive_tanh_erf(self, o, rho_old = 0, switchidx1 = 1, rho = 1, switchidx2 = 10, Zs=False, logprobs_Zs = False, convert2sec = 1.0, alpha = 4.0, lw = 3.5,  sde_args =  {'sigma':.05, 'thetaY': .8, 'thetaZ':.4, 'lw':3.5, 'target_Y':.0, 'v':1}, sde_args_lc = {'sigma':.08, 'thetaY':5.2, 'thetaZ':.14, 'lw':3.5, 'target_Y':3.5, 'v':1}, N_discretize_Z = 1e3, backlookingstep = 20, step = 20, n_mcint = 1e5):
        """
        Calculate the logprob when the two break points (as index in the observation) are given 
        
        The two switching points are given as rho_old = 0, switchidx1 = 1, rho = 1, switchidx2 = 10
        
        
        The logprob is calculated using:
        
            - OBJECTIVE_omit_smaller_number4optimize_mc_erf(o, sigma=.05, thetaY = .05, thetaZ = .05, lw = 3.5, alpha = 4.0,  n_mcint = 1e5, Zs = False, logprobs_Zs = False, convert2sec = .1, limit_observation_len = 20000)
        
        
        --------------------------------------------------
        @input: rho_old = 0, switchidx1 = 1, rho = 1, switchidx2 = 10
        
            rho_old is the initial rho.
            
            switchidx1 is the index of the first switching moment and rho is the first lc direction. 1 means turn to positive lateral direction and -1 means turn to the negative lateral direction.
            
            switchidx2 is the return to the lane keeping stage. 
            
            NOTE that switchidx1 and switchidx2 should be smaller thant the len of the observation-1. 
            
            
        
        @input: sde_arg sde_arg_lc
            
            they are the args if the sde. 
            The first is the lane keeping arg and the second is the lane changing arg. 
        
        --------------------------------------------------
        
        """
        #-------------------------------------
        import builtins
        if len(o)<=backlookingstep:
            raise ValueError('sdsdds')
        
        #------------------------------------
        if isinstance(Zs,bool):
            #
            Zs = np.random.uniform(low = -lw/2.0, high = lw/2.0, size = (len(observation),int(n_mcint)))
            #
            #meanZ and sigmaZ shape are both (len(observation)-1, n_mcint)
            meanZ,sigmaZ = sde_util.meanZ_sigmaZ_given_previousZ_erf(Zs[:-1, :])
            #logprobs_Zs.shape is (len(observation)-1, n_mcint)
            logprobs_Zs = log_normal(Zs[1:,:], meanZ,sigmaZ)
        
        #--------------------
        #THree segments: o1,o2 and o3.
        o1 = o.iloc[:switchidx1:step];o2 = o.iloc[switchidx1:switchidx2:step];o3 = o.iloc[switchidx2::step]
        #
        o1 = o1 - rho_old*lw
        o2 = o2 - (rho_old+rho)*lw
        o3 = o3 - (rho_old+rho)*lw
        #
        logprob_o1 = OBJECTIVE_omit_smaller_number4optimize_mc_erf(observation0 = o1, sigma= sde_args['sigma'], thetaY = sde_args['thetaY'], thetaZ = sde_args['thetaZ'], lw = lw, alpha = alpha,  n_mcint = n_mcint, Zs = Zs, logprobs_Zs = logprobs_Zs, convert2sec = convert2sec, limit_observation_len = 20000)
        #
        logprob_o2 = OBJECTIVE_omit_smaller_number4optimize_mc_erf(observation0 = o2, sigma= sde_args_lc['sigma'], thetaY = sde_args_lc['thetaY'], sde_args_lc = sde_arg['thetaZ'], lw = lw, alpha = alpha,  n_mcint = n_mcint, Zs = Zs, logprobs_Zs = logprobs_Zs, convert2sec = convert2sec, limit_observation_len = 20000)
        #
        logprob_o3 = OBJECTIVE_omit_smaller_number4optimize_mc_erf(observation0 = o3, sigma= sde_args['sigma'], thetaY = sde_args['thetaY'], thetaZ = sde_args['thetaZ'], lw = lw, alpha = alpha,  n_mcint = n_mcint, Zs = Zs, logprobs_Zs = logprobs_Zs, convert2sec = convert2sec, limit_observation_len = 20000)
        
        
        return logprob_o1+logprob_o2+logprob_o3
    
    
    @classmethod
    def logprobs_each_breakpoints_pair(self, o, convert2sec = 1.0, backlookingstep = 20, step = 20, rho_old = 0, switchidx1 = 1, rho = 1, switchidx2 = 10, Zs = False, alpha = 4.0, lw = 3.5,  sde_args =  {'sigma':.05, 'thetaY': .8, 'thetaZ':.4, 'lw':3.5, 'target_Y':.0, 'v':1}, sde_args_lc = {'sigma':.08, 'thetaY':5.2, 'thetaZ':.14, 'lw':3.5, 'target_Y':3.5, 'v':1}, N_discretize_Z = 1e3,n_mcint = 1e5):
        """
        
        """
        #------------------------------------
        if isinstance(Zs,bool):
            #
            Zs = np.random.uniform(low = -lw/2.0, high = lw/2.0, size = (len(o),int(n_mcint)))
            #
            #meanZ and sigmaZ shape are both (len(observation)-1, n_mcint)
            meanZ,sigmaZ = sde_util.meanZ_sigmaZ_given_previousZ_erf(Zs[:-1, :])
            #logprobs_Zs.shape is (len(observation)-1, n_mcint)
            logprobs_Zs = log_normal(Zs[1:,:], meanZ,sigmaZ)
        
        #-------------------------------------------------------------------------------
        #get the switching points. switchidx1s is a list. switchidx2s_dict is a dict. keys are in switchidx1s.
        switchidx1s,switchidx2s_dict = self.BreakingPointsGenerate(o = o, convert2sec = convert2sec, backlookingstep = backlookingstep, step = step)
        
        Niterations = len(switchidx1s)+sum([len(switchidx2s_dict[i]) for i in switchidx1s])
        
        #logprobs keys are str(switchidx1)+'_'+str(switchidx2)
        i=0
        logprobs = {}
        for switchidx1 in switchidx1s:
            for switchidx2 in switchidx2s_dict[switchidx1]:
                #
                
                res = self.logprob_givenbreakingpoints_singlelc_positive_tanh_erf( o = o, rho_old = rho_old, switchidx1 = switchidx1, rho = rho, switchidx2 = switchidx2, Zs = Zs, convert2sec = convert2sec, alpha = alpha, lw = lw,  sde_args =  sde_args, sde_args_lc = sde_args_lc, N_discretize_Z = N_discretize_Z, backlookingstep = backlookingstep, step = step)
                
                logprobs[str(switchidx1)+ '_' + str(switchidx2)] = res
                i = i+1
                print(i, Niterations)
            
        return logprobs
                
    
    
    @classmethod
    def BreakingPointsGenerate(self, o, convert2sec = 1.0, backlookingstep = 20, step = 20):
        """
        
        ------------------------------------------------
        @input: o
            the observation. It is the lateral trajectory. The index are the moments and values are the lateral position. 
            
        @input: convert2sec
        
            o.index*convert2sec unit is the seconds. 
        
        @input:  backlookingstep = 20, step = 20
            
            backlookingstep is how far the algorithm look back, and step is the stride each calculation.
        
        @OUTPUT: 
        
        ----------------------------------------------------
            
        
        """
        
        #generate break points: switchidx1 and switchidx2
        #   they must satisfy the following requirements:
        #       switchidx1 >=backlookingstep, switchidx2>=switchidx1+backlookingstep, switchidx2<=len(o)-backlookingstep
        #       then switchidx1<=switchidx2-backlookingstep<=len(o)-2*backlookingstep
        switchidx1s = range(backlookingstep, len(o)-2*backlookingstep, step)
        switchidx2s_dict = {}
        for switchidx1 in switchidx1s:
            switchidx2s_dict[switchidx1] = range(switchidx1 + backlookingstep, len(o)-2, step)
            
        return switchidx1s,switchidx2s_dict
        
    
    
    @classmethod
    def lc_flags_given_breakpoint(self, len_observations = 100, k = 10):
        """
        Givn the length of the observations and 
        ------------------------------
        @input: len_observations
        
            the length of the observations. 
        
        @input: k
            
            an integer, that indicate the index of the lane changing. 
            
            k must be greater or equal than 1. 
        
        @input: 
        
        
        
        @OUTPUT: lane_flags_positive,lane_flags_negative.
            
            both are list. 
            
            the length is len(len_observations)-1.
            
            postive means lane change to the lateral-increasing-ccordinate direction. 
        ----------------------------------------------
            
            
        """
        #lane keeping, thus it is 
        first_sgage = [0 for i in range(k)]
        
        #
        secondstage_positive = [1 for i in range(k, len(observation)-1)]
        #
        secondstage_negative = [-1 for i in range(k, len(observation)-1)]
        
        
        return first_sgage+secondstage_positive,first_sgage+secondstage_negative
        
    
    
    
    @classmethod
    def lc_logprob_given_breakpoint(self, observation, k =1, lw = 3.5,args_sde = {'theta':.5, 'sigma':.01}):
        """
        THe log prob when the lane changing point occurs at k-th 
        -------------------------------------------
        
        @input: observation
            
            a pd.Series. THe index should be time and unit is sec. 
        
        @output: logprob_positive,logprob_negative
            
            the log prob to the positive lateral coor or toward to the negative lateral coor. 
        ------------------------------------------
        
        """
        #get the lane changing flags.
        #   both are list of 0,1 or -1. 
        #   will be used in the meanY and sigmaY calculating. 
        rho_es_positive,rho_es_negative = self.GenerateLaneIdentifications(len_observation = len(observation), k = k, lw = lw)
        
        #
        
        
        
        
        
        pass


    @classmethod
    def LC_initization_identification_erf_simple(self, observation, args_sde = {'thetaY':.5, 'thetaZ':.5, 'sigma':.01}, delta_roll_sec = 1, convert2sec = .04, lw = 3.5, n_mcint = 1e5, start_rho = 0):
        """
        identify the lc initialization moment for the erf function. 
        Difference between 
            
            - self.LC_initization_identification_erf_simple()
            - self.LC_initization_identification_erf
            - self.LC_initization_identification()
        
        The first one simplify the transition point. 
        The second focus on the sde derived from the error function. 
        
        
        Basic method is as follows:
            
            at each index (say idx) in the observation, compare the following two cases:
            
                - for t in 0~index, all are lane keeping case
                - for t in 0 to index-delta_back_sec, it is lane keeping; while for delta_back_sec to current moment (idx), lane changing. 
                
            Once it is found that the lane changing prob is greater than the lane keeping, then set the 
        
        ----------------------------------------
        @input: start_rho
        
            the rho of the start trajectory. 
            
            If the begining trajetory's lateral position locates at lane the middle of which is zero, then start_rho=0. 
            
        @input: n_mcint
        
            the number of sampling points for integration. 
        
        @input: convert2sec
        
            the factor that conver the time index to seconds. 
            
            for the highD and Mirror, it is 0.04 (or there are 25 moments within one second)
            
            For ngsim, it is 0.1 (there are 10 frames within one second. )
            
            
        @input: observation
            
            pd.Series. An observation. 
            
            observation.index is the moments. 
        
        
        @input: delta_back_sec
            
        @input: lc_checking_threshold_sec
            
            continuously check the obseervation. Until
        
        @output: lc_changings
        
            a dict. 
            
            lc_changings[idx] = rho_after_changing. 
            
                
        ---------------------------------
        
        """
        
        t = time.time()
        
        def mc_single_point_sum(log_probs):
            #print(np.format_float_scientific(log_probs[0]).split('+')[1])
            #given many points from monte carlao. calculate the log of the mean of the prob. 
            #
            digitsNs = []
            for i in log_probs:
                if i<=-1:
                    #print(i, np.format_float_scientific(i).split('+'), '--')
                    #print(int(np.format_float_scientific(i).split('+')[1]))
                    digitsNs.append(int(np.format_float_scientific(i).split('+')[1]))
                else:
                    #when i<=-1, then np.format_float_scientific(-.9) return '-9.e-01'
                    digitsNs.append(0)
            #digitsNs = [int(np.format_float_scientific(i).split('+')[1]) for i in log_probs]
            digitsMIN = min(digitsNs)
            
            #the condition >=-10 is to omit the smaller number. 
            sample_logs_added = [i for i in 10**(digitsMIN)+log_probs if i>=-10]
            
            return -10**digitsMIN+np.log(np.mean(np.exp(sample_logs_added)))
            
        
        #the frames. 
        frames_delta_roll_sec = max(1, int(1.0*delta_roll_sec/convert2sec))
        
        #current_idx = 1 to len(observation)-1
        current_idx = frames_delta_roll_sec
        
        #
        Zs = np.random.uniform(low = -lw/2.0, high = lw/2.0, size = (int(n_mcint),))
        
        #lc_changings[idx] = 0, -1, or 1
        lc_changings = {0:start_rho}
        current_rho = start_rho
        
        
        while current_idx<len(observation)-1:
            
            #print(len(observation))
            #------------calculate the prob of rho=0, -1 and 1 respectively. 
            #calculate the prob of rho = 0
            #the index that the logprobs need to be caculated should be 
            #   indexs = observation.index[current_idx-frames_delta_back_sec+1:current_idx+1]
            Y_i1_es = observation.iloc[current_idx-frames_delta_roll_sec+1:current_idx+1]
            Y_i_es = observation.iloc[current_idx-frames_delta_roll_sec:current_idx]
            #
            #lane keeping
            logs_res_lanekeeping = []
            for Y_i1,Y_i in zip(Y_i1_es, Y_i_es):
                tmp = self.logprob_Y_i_given_previousYandZ_erf(Y_i1, Y_i, Z_i_1 = Zs, deltat = .5, lw = lw, args_sde = args_sde, rho = current_rho + 0)
                builtins.tmp = tmp
                logs_res_lanekeeping.append(mc_single_point_sum(tmp))
            #turn to positive
            logs_res_lc_positive = []
            for Y_i1,Y_i in zip(Y_i1_es, Y_i_es):
                tmp = self.logprob_Y_i_given_previousYandZ_erf(Y_i1, Y_i, Z_i_1 = Zs, deltat = .5, lw = lw, args_sde = args_sde, rho = current_rho + 1)
                logs_res_lc_positive.append(mc_single_point_sum(tmp))
            #turn to negative
            logs_res_lc_negative = []
            for Y_i1,Y_i in zip(Y_i1_es, Y_i_es):
                tmp = self.logprob_Y_i_given_previousYandZ_erf(Y_i1, Y_i, Z_i_1 = Zs, deltat = .5, lw = lw, args_sde = args_sde, rho = current_rho -1)
                logs_res_lc_negative.append(mc_single_point_sum(tmp))
            
            #-------------Compare lane keeping case and lane changing case. 
            #print(logs_res_lanekeeping)
            #print(sum(logs_res_lanekeeping), sum(logs_res_lc_positive), sum(logs_res_lc_negative))
            to_compares = [sum(logs_res_lanekeeping), sum(logs_res_lc_positive), sum(logs_res_lc_negative)]
            #print(sum(logs_res_lanekeeping))
            max_idx = to_compares.index(max(to_compares))
            if max_idx==0:
                
                #lane keeping, do nothing. 
                pass
            elif max_idx==1:
                #lane changing to positive direction. 
                continuous_check_counter_sec = 0
                
                print('Lane changing to positive')
                
                lc_changings[current_idx] = 1
                current_rho = current_rho + 1
                
            elif max_idx==2:
                
                print('Lane changing to Negative')
                
                lc_changings[current_idx] = -1
                current_rho = current_rho -1

            
            #move time forward. 
            current_idx = min(current_idx + frames_delta_roll_sec, len(observation)-1)
            
            print(current_idx,  time.time() - t, to_compares)
            t =  time.time()
            
        return lc_changings


    @classmethod
    def LC_initization_identification_erf(self, observation, args_sde = {'thetaY':.5, 'thetaZ':.5, 'sigma':.01}, delta_roll_sec = 1, delta_back_sec = 1, convert2sec = .04, lw = 3.5, lc_checking_threshold_sec = 15, n_mcint = 1e5, start_rho = 0):
        """
        identify the lc initialization moment for the erf function. 
        Difference between 
            - self.LC_initization_identification_erf
            - self.LC_initization_identification()
            
        The former focus on the sde derived from the error function. 
        
        
        Basic method is as follows:
            
            at each index (say idx) in the observation, compare the following two cases:
            
                - for t in 0~index, all are lane keeping case
                - for t in 0 to index-delta_back_sec, it is lane keeping; while for delta_back_sec to current moment (idx), lane changing. 
                
            Once it is found that the lane changing prob is greater than the lane keeping, then set the 
        
        ----------------------------------------
        @input: start_rho
        
            the rho of the start trajectory. 
            
            If the begining trajetory's lateral position locates at lane the middle of which is zero, then start_rho=0. 
            
        @input: n_mcint
        
            the number of sampling points for integration. 
        
        @input: convert2sec
        
            the factor that conver the time index to seconds. 
            
            for the highD and Mirror, it is 0.04 (or there are 25 moments within one second)
            
            For ngsim, it is 0.1 (there are 10 frames within one second. )
            
            
        @input: observation
            
            pd.Series. An observation. 
            
            observation.index is the moments. 
        
        
        @input: delta_back_sec
            
        @input: lc_checking_threshold_sec
            
            continuously check the obseervation. Until
        
        @output: rhos
            
            a list of 0, 1 or -1.
            
            0 means lane keeping and -1 and 1 means either turn left or turn right. 
            
            len(rhos) = len(observation)-1. 
            
            rhos[0] is the sign of the rho of observation.iloc[1]
                
        ---------------------------------
        
        """
        #
        #logprobs is a series. len(logprobs) = len(observation)-1
        logprobs = pd.Series(dtype = 'float64', index = observation.index[:])
        logprobs[observation.index[0]] = 1.0
        
        #the frames. 
        frames_delta_back_sec = max(1, int(1.0*delta_back_sec/convert2sec))
        frames_delta_roll_sec = max(1, int(1.0*delta_roll_sec/convert2sec))
        
        #current_idx = 1 to len(observation)-1
        current_idx = frames_delta_back_sec
        
        #
        Zs = np.random.uniform(low = -lw/2.0, high = lw/2.0, size = (int(n_mcint),))
        
        #
        current_rho = start_rho
        while current_idx<=len(observation)-1:
            
            #------------calculate the prob of rho=0, -1 and 1 respectively. 
            #calculate the prob of rho = 0
            #the index that the logprobs need to be caculated should be 
            #   indexs = observation.index[current_idx-frames_delta_back_sec+1:current_idx+1]
            Y_i1_es = observation.iloc[current_idx-frames_delta_back_sec+1:current_idx+1]
            Y_i_es = observation.iloc[current_idx-frames_delta_back_sec:current_idx]
            #
            logs_res_lanekeeping = []
            for Y_i1,Y_i in zip(Y_i1_es, Y_i_es):
                tmp = self.logprob_Y_i_given_previousYandZ_erf(Y_i1, Y_i, Z_i_1 = Zs, deltat = .5, lw = lw, args_sde = args_sde, rho = current_rho + 0)
                logs_res_lanekeeping.append(tmp)
            #turn to positive
            logs_res_lc_positive = []
            for Y_i1,Y_i in zip(Y_i1_es, Y_i_es):
                tmp = self.logprob_Y_i_given_previousYandZ_erf(Y_i1, Y_i, Z_i_1 = Zs, deltat = .5, lw = lw, args_sde = args_sde, rho = current_rho + 1)
                logs_res_lc_positive.append(tmp)
            #turn to negative
            logs_res_lc_negative = []
            for Y_i1,Y_i in zip(Y_i1_es, Y_i_es):
                tmp = self.logprob_Y_i_given_previousYandZ_erf(Y_i1, Y_i, Z_i_1 = Zs, deltat = .5, lw = lw, args_sde = args_sde, rho = current_rho -1)
                logs_res_lc_negative.append(tmp)
            
            #-------------Compare lane keeping case and lane changing case. 
            to_compares = [sum(logs_res_lanekeeping), sum(logs_res_lc_positive), sum(logs_res_lc_negative)]
            max_idx = to_compares.index(max(to_compares))
            if max_idx==0:
                
                #lane keeping, do nothing. 
                pass
            elif max_idx==1:
                #lane changing to positive direction. 
                continuous_check_counter_sec = 0
                
                
                
                
                
            elif max_idx==2:
                
                
                
                pass
                
            
            
            
            #move time forward. 
            current_idx = min(current_idx + frames_delta_roll_sec, len(observation)-1)
            


    @classmethod
    def LC_initization_identification(self, observation, args_sde = {'theta':.5, 'sigma':.01}):
        """
        
        
        """
        
        
        pass
    
    
    
    @classmethod
    def MonteCarloSimulationPaths_tanh_noise_parabolic_Y(self,  N_paths  = 100, Z0 = np.array([.0, .0]), tspan = np.linspace(0.0, 90, 100), sde_args = {'sigma':.03, 'thetaY':.1, 'thetaZ':.01, 'lw':3.5, 'target_Y':3.5, 'v':1, 'M':1e6, 'A':1.0}):
        """
        
        """
        import sdeint

        sigma = sde_args['sigma']
        thetaY = sde_args['thetaY']
        thetaZ = sde_args['thetaZ']
        lw = sde_args['lw']
        target_Y = sde_args['target_Y']
        v = sde_args['v']
        M = sde_args['M']
        A = sde_args['A']

        #
        res_Y = []
        res_Z = []
        
        #--------------------------
        def parabolic(y):
            return -A*y**2*(y/np.sqrt(1.0/M+y**2))
                
                
        def f_Z(Z,t, sigma = sigma, ):
            #normalize2_1(Z[0]-target_Y,M)*
            tmp0 = thetaY*parabolic(Z[0]-target_Y) - thetaZ* Z[1]
            tmp1 = -2.0*Z[1]*(v**2)*(sigma**2)*(lw**2-4.0*(Z[1]**2))/2/(lw**2)
            #tmp1 = -.5*2*lanewidth/np.pi*sigma**2*np.tan(Z[1]*np.pi/lanewidth)/((np.tan(Z[1]*np.pi/lanewidth)**2+1)**2)
            #return np.array([[tmp0, 0], [0, tmp1]])
            return np.array([tmp0, tmp1])
        def g_Z(Z,t, ):
            tmp0 = 0
            tmp1 = lw*v*sigma*(lw**2 - 4.0*( Z[1]**2))/2/(lw**2)
            #tmp1 = 1.0*lanewidth/np.pi*sigma/(np.tan(Z[1]*np.pi/lanewidth)**2+1)
            return np.array([[0, 0], [0, tmp1]])
        
        
        #---------------------------------------------
        for iterr in range(N_paths):
            #result.shape is (len(tspan),2), '2' corresponds to Z.
            #   the first column is lateral displacement, and the 2nd columns is the noize. 
            result = sdeint.itoint(f_Z, g_Z, Z0, tspan)
            res_Y.append(result[:,0])
            res_Z.append(result[:,1])
        #-----------------------------------------
        #np.array(res_Y) shape is N_paths,len(tspan)
        #np.array(res_Z) shape is N_paths,len(tspan)
        return pd.DataFrame(res_Y, columns =  tspan).T,pd.DataFrame(res_Z, columns =  tspan).T
    
    
    @classmethod
    def logprob_Y_i_given_previousYandZ_erf(self, Y_i1, Y_i, Z_i_1, deltat = .5, lw = 3.5, args_sde = {'thetaY':.5, 'thetaZ':.5, 'sigma':.01, 'alpha':4.0}, multiplier = 1.0, rho = 0):
        """
        Get the log prob. 
        
        Difference:
        
            - logprob_Y_i_given_previousYandZ(), the probability from -Y. 
            -logprob_Y_i_given_previousYandZ_erf(), the probability from error function. 
        
        
        
        Z_i_1 can be scalar or 1d array. 
        
        If it is array, then the returned is also an array with the same length. 
        -------------------------------------------
        @input: rho
            the deviation is calculated as:
            
                y-rho*lw
            
        @input: Y_i1, Y_i, Z_i_1
        
            Y_i1 is the Y at moment i+1
            Y_i is the Y at moment i
            Z_i_1 is the Z at moment i-1.
            
        @OUTUPUT: prob
        
        """
        #if Z_i_1 is scalar, then meany,sigmasquarey are also scalar
        #   if Z_i_1 is 1d array, then meany,sigmasquarey are also array with the same length. 
        meany,sigmasquarey = self.meanY_sigmaY_given_previousYandZ_erf(Y_i = Y_i, Z_i_1 = Z_i_1, deltat = deltat, lw = lw, args_sde = args_sde, rho = rho)
        
        #print(Y_i1, meany,sigmasquarey)
        #print(normal(Y_i1, meany, sigmasquarey))
        return multiplier*log_normal(Y_i1, meany, sigmasquarey)

    
    @classmethod
    def logprob_Y_i_given_previousYandZ(self, Y_i1, Y_i, Z_i_1, deltat = .5, lw = 3.5, args_sde = {'theta':.5, 'sigma':.01},multiplier = 1.0):
        """
        Z_i_1 can be scalar or 1d array. 
        
        If it is array, then the returned is also an array with the same length. 
        -------------------------------------------
        @input: Y_i1, Y_i, Z_i_1
        
            Y_i1 is the Y at moment i+1
            Y_i is the Y at moment i
            Z_i_1 is the Z at moment i-1.
            
            
            
        @OUTUPUT: prob
        
        """
        #if Z_i_1 is scalar, then meany,sigmasquarey are also scalar
        #   if Z_i_1 is 1d array, then meany,sigmasquarey are also array with the same length. 
        meany,sigmasquarey = self.meanY_sigmaY_given_previousYandZ(Y_i = Y_i, Z_i_1 = Z_i_1, deltat = deltat, lw = lw, args_sde = args_sde)
        
        #print(Y_i1, meany,sigmasquarey)
        #print(normal(Y_i1, meany, sigmasquarey))
        return multiplier*log_normal(Y_i1, meany, sigmasquarey)

        
    @classmethod
    def BKP_logprobs_given_single_series4integration(self, Z, observation, lw = 3.5,  args_sde = {'theta':.5, 'sigma':.01}, multiplier = 1.0):
        """
        The prob given single series.
        
        @input: observation
            
            a pd.Series. THe index should be time and unit is sec. 
        
        @input: Z
            a array, len(Z)=len(observation)-1
            
        """
        
        #deltat, a np.array, 1d
        deltats = np.diff(observation.index)
        
        #
        logprobs = []
        for idx in range(1, len(observation)):
            #
            deltat = deltats[idx-1]
            z = Z[idx-1]
            y = observation.iloc[idx-1]
            y1 = observation.iloc[idx] 
            
            #print(z, y, y1, deltat)
            
            #single prob
            logprob = self.logprob_Y_i_given_previousYandZ(Y_i1 = y1, Y_i = y, Z_i_1 = z, deltat = deltat, lw = lw, args_sde = args_sde, multiplier = multiplier)
            
            logprobs.append(logprob)
            
        return logprobs
        #return np.prod(probs)
        
    
    @classmethod
    def LOGPROBS_given_single_seriesmultiple_Z_4integration_erf(self, Zs, observation, lw = 3.5, args_sde = {'thetaY':.5, 'thetaZ':.5, 'sigma':.01, 'alpha':4.0}, multiplier = 1.0, convert2sec = .1):
        """
        The prob given single series.
        
        Difference :
        
            - logprobs_given_single_seriesmultiple_Z_4integration(), the dY is constructed from -Y.
            - logprobs_given_single_seriesmultiple_Z_4integration_erf(), the dY is construnted from erro function. 
            
            
            
        --------------------------------------------------------------------
        
        @input: Zs
        
            an array. shape is (len(observation)-1, N_montecarlo)
            
        @input: observation
            
            a pd.Series. THe index should be time and unit is sec. 
        
        @input: Z
            a array, len(Z)=len(observation)-1
            
        """
        
        #deltat, a np.array, 1d
        deltats = convert2sec*np.diff(observation.index)
        #
        logprobs = []
        for idx in range(1,len(observation)):
            #
            deltat = deltats[idx-1]
            z = Zs[idx, :]#array. 
            y = observation.iloc[idx-1]
            y1 = observation.iloc[idx] 
            
            #print(z, y, y1, deltat)
            
            #single prob
            logprob = self.logprob_Y_i_given_previousYandZ_erf(Y_i1 = y1, Y_i = y, Z_i_1 = z, deltat = deltat, lw = lw, args_sde = args_sde, multiplier = multiplier)
            
            logprobs.append(logprob)
        
        #shape is the same as Zs
        return np.array(logprobs)
    
    @classmethod
    def logprobs_given_single_seriesmultiple_Z_4integration_erf(self, Zs, observation, lw = 3.5, args_sde = {'thetaY':.5, 'thetaZ':.5, 'sigma':.01}, multiplier = 1.0):
        """
        The prob given single series.
        
        Difference :
        
            - logprobs_given_single_seriesmultiple_Z_4integration(), the dY is constructed from -Y.
            - logprobs_given_single_seriesmultiple_Z_4integration_erf(), the dY is construnted from erro function. 
            
            
            
        --------------------------------------------------------------------
        
        @input: Zs
        
            an array. shape is (len(observation)-1, N_montecarlo)
            
        @input: observation
            
            a pd.Series. THe index should be time and unit is sec. 
        
        @input: Z
            a array, len(Z)=len(observation)-1
            
        """
        
        #deltat, a np.array, 1d
        deltats = np.diff(observation.index)
        
        #
        logprobs = []
        for idx in range(1, len(observation)):
            #
            deltat = deltats[idx-1]
            z = Zs[idx-1, :]#array. 
            y = observation.iloc[idx-1]
            y1 = observation.iloc[idx] 
            
            #print(z, y, y1, deltat)
            
            #single prob
            logprob = self.logprob_Y_i_given_previousYandZ_erf(Y_i1 = y1, Y_i = y, Z_i_1 = z, deltat = deltat, lw = lw, args_sde = args_sde, multiplier = multiplier)
            
            logprobs.append(logprob)
        
        #shape is the same as Zs
        return np.array(logprobs)
    
    @classmethod
    def logprobs_given_single_seriesmultiple_Z_4integration(self, Zs, observation, lw = 3.5,  args_sde = {'theta':.5, 'sigma':.01}, multiplier = 1.0):
        """
        The prob given single series.
        
        @input: Zs
        
            an array. shape is (len(observation)-1, N_montecarlo)
            
        @input: observation
            
            a pd.Series. THe index should be time and unit is sec. 
        
        @input: Z
            a array, len(Z)=len(observation)-1
            
        """
        
        #deltat, a np.array, 1d
        deltats = np.diff(observation.index)
        
        #
        logprobs = []
        for idx in range(1, len(observation)):
            #
            deltat = deltats[idx-1]
            z = Zs[idx-1, :]#array. 
            y = observation.iloc[idx-1]
            y1 = observation.iloc[idx] 
            
            #print(z, y, y1, deltat)
            
            #single prob
            logprob = self.logprob_Y_i_given_previousYandZ(Y_i1 = y1, Y_i = y, Z_i_1 = z, deltat = deltat, lw = lw, args_sde = args_sde, multiplier = multiplier)
            
            logprobs.append(logprob)
        
        #shape is the same as Zs
        return np.array(logprobs)
        return logprobs
        #return np.prod(probs)

    @classmethod
    def meanZ_sigmaZ_given_previousZ_erf(self, Z_i_1, deltat = .5, lw = 3.5, args_sde = {'thetaY':.5, 'thetaZ':.5, 'sigma':.01}, rho = 0):
        """
        This method return the mean and the sigma of Y_i+1 given Y_i and Z_i_2. 
        
        The parameters include sigma, thetaY and thetaZ. 
        
        DIfference:
        
            - self.meanY_sigmaY_given_previousYandZ()
            
            - self.meanY_sigmaY_given_previousYandZ_erf(), constructed from error fucntion. 
            
        NOTE:
            
            if Y_i, Z_i_1 are both float, then the returned value is also scalar. 
        
            While if Z_i_1 is an array, the returned meanY,sigmasquareY are both array. size the same as Z_i_1.
        
        ----------------------------------------
        @input: rho
            the lane changing flag
            
            If is either 0, or 1 or -1. 
            
        @input: Y_i, Z_i_1
            
            the input methods. 
        
            Z_i_1 can be an array. Then the returned value is also array with the same length. 
        
        @input: deltat
            the time interval in sec
        
        @input: lw
            lane width in meter. 
            
        
        """
        #
        thetaY = args_sde['thetaY']
        thetaZ = args_sde['thetaZ']
        sigma = args_sde['sigma']
        
        meanZ = None
        sigmasquareZ = None
        
        #   Z_i_1*np.pi/lw 
        #       np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)**2
        #meanY = -(Z_i_1 - lw*sigma*sigma*np.tan(Z_i_1*np.pi/lw )/(np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)**2)*deltat)*theta*deltat + (Y_i - theta*Y_i*deltat + rho*deltat*theta*lw)
        #meanY = Y_i + deltat*thetaY*Transformed_erf_integrator(Y_i,lw = lw, rho = rho) - deltat*thetaZ*(Z_i_1-deltat*2.0*Z_i_1*(sigma**2)*(lw**2 - 4.0*Z_i_1**2)/(2.0*lw**2))
        meanZ = Z_i_1-deltat*2.0*Z_i_1*(sigma**2)*(lw**2 - 4.0*Z_i_1**2)/(2.0*lw*lw)

        #sigmasquareY = (sigma*lw/(np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)))**2*deltat*((theta*deltat)**2)
        #sigmasquareY = deltat*((thetaZ*deltat)**2)*((sigma*lw)**2)*((lw**2 - 4.0*Z_i_1**2)/(2.0*lw*lw))**2
        sigmasquareZ = deltat*((sigma*lw)**2)*((lw**2 - 4.0*(Z_i_1**2))/(2.0*lw*lw))**2
        
        return meanZ,sigmasquareZ


    @classmethod
    def logprobs_with_different_switchs(self, o, rho_old = 0, rho = 1, Zs = False, logprobs_Zs = False, convert2sec = 1.0, alpha = 4.0, lw = 3.5,  sde_args =  {'sigma':.01, 'thetaY': .8, 'thetaZ':.4, 'lw':3.5, 'target_Y':.0, 'v':1}, sde_args_lc = {'sigma':.01, 'thetaY':5.2, 'thetaZ':.14, 'lw':3.5, 'target_Y':3.5, 'v':1}, N_discretize_Z = 1e3, step = 10, n_mcint = 1e5):
        """
        Calculate the logprob when the two break points (as index in the observation) are given 
        
        The two switching points are given as rho_old = 0, switchidx1 = 1, rho = 1, switchidx2 = 10
        
        
        The logprob is calculated using:
        
            - OBJECTIVE_omit_smaller_number4optimize_mc_erf(o, sigma=.05, thetaY = .05, thetaZ = .05, lw = 3.5, alpha = 4.0,  n_mcint = 1e5, Zs = False, logprobs_Zs = False, convert2sec = .1, limit_observation_len = 20000)
        
        
        --------------------------------------------------
        @input: rho_old = 0, switchidx1 = 1, rho = 1, switchidx2 = 10
        
            rho_old is the initial rho.
            
            switchidx1 is the index of the first switching moment and rho is the first lc direction. 1 means turn to positive lateral direction and -1 means turn to the negative lateral direction.
            
            switchidx2 is the return to the lane keeping stage. 
            
            NOTE that switchidx1 and switchidx2 should be smaller thant the len of the observation-1. 
            
            
        
        @input: sde_arg sde_arg_lc
            
            they are the args if the sde. 
            The first is the lane keeping arg and the second is the lane changing arg. 
        
        --------------------------------------------------
        
        """
        switchidxs = range(2*step, len(o)-2*step)
        logprobs = {}
        for i in switchidxs:
            
            logprobs[i] = self.logprob_observation_single_switch(o = o, rho_old = rho_old, switchidx = i, rho = rho, Zs=Zs, logprobs_Zs = logprobs_Zs, convert2sec = convert2sec, alpha = alpha, lw = lw,  sde_args =  sde_args, sde_args_lc = sde_args_lc, N_discretize_Z = N_discretize_Z, step = step, n_mcint = n_mcint)
        
        return logprobs
        


    @classmethod
    def logprob_observation_single_switch(self, o, rho_old = 0, switchidx = 10, rho = 1, Zs=False, logprobs_Zs = False, convert2sec = 1.0, alpha = 4.0, lw = 3.5,  sde_args =  {'sigma':.05, 'thetaY': .8, 'thetaZ':.4, 'lw':3.5, 'target_Y':.0, 'v':1}, sde_args_lc = {'sigma':.08, 'thetaY':5.2, 'thetaZ':.14, 'lw':3.5, 'target_Y':3.5, 'v':1}, N_discretize_Z = 1e3, step = 25, n_mcint = 1e5):
        """
        Calculate the logprob when the two break points (as index in the observation) are given 
        
        The two switching points are given as rho_old = 0, switchidx1 = 1, rho = 1, switchidx2 = 10
        
        
        The logprob is calculated using:
        
            - OBJECTIVE_omit_smaller_number4optimize_mc_erf(o, sigma=.05, thetaY = .05, thetaZ = .05, lw = 3.5, alpha = 4.0,  n_mcint = 1e5, Zs = False, logprobs_Zs = False, convert2sec = .1, limit_observation_len = 20000)
        
        
        --------------------------------------------------
        @input: rho_old = 0, switchidx1 = 1, rho = 1, switchidx2 = 10
        
            rho_old is the initial rho.
            
            switchidx1 is the index of the first switching moment and rho is the first lc direction. 1 means turn to positive lateral direction and -1 means turn to the negative lateral direction.
            
            switchidx2 is the return to the lane keeping stage. 
            
            NOTE that switchidx1 and switchidx2 should be smaller thant the len of the observation-1. 
            
            
        
        @input: sde_arg sde_arg_lc
            
            they are the args if the sde. 
            The first is the lane keeping arg and the second is the lane changing arg. 
        
        --------------------------------------------------
        
        """
        if switchidx<=2*step:
            switchidx = 2*step
        
        #-------------------------------------
        import builtins
        if len(o)<=step*2:
            raise ValueError('sdsdds')
        
        
        #--------------------
        #THree segments: o1,o2 and o3.
        o1 = o.iloc[:switchidx:step]
        o2 = o.iloc[switchidx::step]
        #
        o1 = o1 - rho_old*lw
        o2 = o2 - (rho_old+rho)*lw
        
        #print(len(o1), len(o2))
        #
        logprob_o1 = OBJECTIVE_omit_smaller_number4optimize_mc_erf(observation0 = o1, sigma= sde_args['sigma'], thetaY = sde_args['thetaY'], thetaZ = sde_args['thetaZ'], lw = lw, alpha = alpha,  n_mcint = n_mcint, Zs = False, logprobs_Zs = False, convert2sec = convert2sec, limit_observation_len = 20000)
        #
        logprob_o2 = OBJECTIVE_omit_smaller_number4optimize_mc_erf(observation0 = o2, sigma= sde_args_lc['sigma'], thetaY = sde_args_lc['thetaY'], thetaZ = sde_args_lc['thetaZ'], lw = lw, alpha = alpha,  n_mcint = n_mcint, Zs = False, logprobs_Zs = False, convert2sec = convert2sec, limit_observation_len = 20000)
        
        
        return logprob_o1+logprob_o2
        


    @classmethod
    def lc_identification_single_transition_tanh_erf_simple(self, observation0, lookingbacksteps = 4, steplength = 25,  Zs=False,convert2sec = 1/25.0, lw = 3.5, rhos = [1, 0, -1], rho_old = 0, alpha = 3.0, sde_args =  {'sigma':.05, 'thetaY': .8, 'thetaZ':.4, 'lw':3.5, 'target_Y':.0, 'v':1}, sde_args_lc = {'sigma':.08, 'thetaY':5.2, 'thetaZ':.14, 'lw':3.5, 'target_Y':3.5, 'v':1}, N_discretize_Z = 1e3):
        """
        The lane changing moment identification. simple in the method name means that there is only one lane changing occurrence. 
        
        The tanh means that the error function is tanh and erf means the lateral dynamics is erf type. 
        
        
        Note that the observation is normalized, i.e. the trajectory begins from the lane whose middle line is defined as the y=0. 
        
        ----------------------------------------------
        
        @input: N_discretize_Z
            the number of discretization of Z. 
            
            The discretized Z are within [-lw/2, lw/2]
        
        
        @input: rho_old
            
            the rho before the lane changing. 
            
        @input: rhos
            
            a list. 
            
            either [1] or [-1] or [1, -1]. 
            
            if only [1] meahs only turn to the positive lateral coor. 
            
            only [-1] means only turn to negative lateral coor. 
            
        @input: observation
            
            a series which represents the observation. THe index are either moments or the longitudinal coordinates. 
        
        @input: lw
        
            the lane width. 
        
        @input: YC_pre and YC_new
        
            YC_pre is the previsous lane middle line coordinate (which is the lane previsous of the last lane changing. )
            
            YC_new is the middle line of the lane before the detected lane changing (the lane should be also where the trajectory begins. )
        
        """
        import builtins
        if len(observation0)<=lookingbacksteps:
            raise ValueError('sdsdds')
        
        #------------------------------------
        if isinstance(Zs,bool):
            
            Zs = np.random.uniform(low = -lw/2.0, high = lw/2.0, size = int(N_discretize_Z))
            #Zs = np.linspace(-lw/2, lw/2, int(N_discretize_Z))
            #Zs = Zs[1:-1]
        #--------------------------------------------------------------
        #RETURNED VALUES
        #key are rhos and values are logprob
        #logprobs = {rho:{i:None for i in range(2, len(observation0))} for rho in rhos}
        #logprobs[rho] = (ks, logprob), both are list. 
        logprobs = {rho:{} for rho in rhos}
        
        observation0 = observation0.sort_index()
        #
        deltats = convert2sec*np.diff(observation0.index)
        Y_i_es = observation0.iloc[1:].values
        Y_i_1_es = observation0.iloc[:-1].values
        
        
        BaseChanged = False
        #--------------------------------------------------------------
        #LOGPROBS[rho][idx] = logprob
        LOGPROBS = {rho:pd.Series(index = observation0.index[lookingbacksteps-1::steplength]) for rho in rhos}
        #
        #print(range(len(observation0))[rollingsteps-1::steplength])
        for idx in range(len(observation0))[lookingbacksteps-1::steplength]:
            #
            Yi_es_series = observation0.iloc[idx-lookingbacksteps+1:idx+1]
            #
            for rho0 in rhos:
                #
                #rho = max(0, min(1, rho_old+rho0))
                rho = rho_old+rho0
                if rho0!=0 and rho_old!=0:
                    #
                    res = self.Conditional_given_Ys(Yi_es_series, convert2sec = convert2sec, rho = rho, lw = lw, alpha = alpha, Zs = Zs, sde_args = sde_args_lc, N_discretize_Z = N_discretize_Z)
                else:
                    res = self.Conditional_given_Ys(Yi_es_series, convert2sec = convert2sec, rho = rho, lw = lw, alpha = alpha, Zs = Zs, sde_args = sde_args, N_discretize_Z = N_discretize_Z)
                
                #LOGPROBS[rho0].append(res)
                LOGPROBS[rho0][observation0.index[idx]] = res
            
            #
            #determine which one is the greatest. 
            #print(LOGPROBS)
            tmp = [LOGPROBS[rho][observation0.index[idx]] for rho in rhos]
            rho_with_max = rhos[tmp.index(max(tmp))]
            #rho_old = max(0, min(1, rho_old + rho_with_max))
            if not BaseChanged:
                rho_old = rho_old + rho_with_max
                BaseChanged = True
        
        #for rho0 in rhos:LOGPROBS[rho0]=pd.Series(LOGPROBS[rho0])
        
        return LOGPROBS
            


    @classmethod
    def find_break_points_given_logprobs_dict(self, logprobs):
        """
        Find the break index fiven the log of the probabilities. 
        -------------------------------
        @input:  logprobs
            
            a dict. 
            The keys should be -1, 0, 1 which means turn right, lane keep and turn left. 
        
            logprobs[0] is a list. 
        
        @OUTPUT: breakpoints_dict
            
            breakpoints_dict[idx] = 0, means that after the moment idx, key 0 has the greatest value. 
        
        """
        breakpoints_dict = {}
        keys = list(logprobs.keys())
         
        #the greatest value of the onset moment
        firstvalues = [logprobs[k][0] for k in keys]
        breakpoints_dict[0] = keys[firstvalues.index(max(firstvalues))]
        #
        pastmaxkey = keys[firstvalues.index(max(firstvalues))]
        #
        for i in range(1, len(logprobs[keys[0]])):
            #get the key that have greatest value. 
            values = [logprobs[k][i] for k in keys]
            maxkey = keys[values.index(max(values))]
            if maxkey!=pastmaxkey:
                breakpoints_dict[i] = maxkey
                #
                pastmaxkey = maxkey
                
        return breakpoints_dict



    @classmethod
    def lc_identification_single_transition_tanh_erf_single_steps(self, observation0, convert2sec = 1/25.0, lw = 3.5, rhos = [1, 0, -1], rho_old = 0, alpha = 3.0, sde_args = {'sigma':.04, 'thetaY':2.0, 'thetaZ':1.0, 'lw':3.5, 'target_Y':.0, 'v':1, 'alpha':3.0}, N_discretize_Z = 1e3):
        """
        The lane changing moment identification. simple in the method name means that there is only one lane changing occurrence. 
        
        The tanh means that the error function is tanh and erf means the lateral dynamics is erf type. 
        
        
        Note that the observation is normalized, i.e. the trajectory begins from the lane whose middle line is defined as the y=0. 
        
        ----------------------------------------------
        
        @input: N_discretize_Z
            the number of discretization of Z. 
            
            The discretized Z are within [-lw/2, lw/2]
        
        
        @input: rho_old
            
            the rho before the lane changing. 
            
        @input: rhos
            
            a list. 
            
            either [1] or [-1] or [1, -1]. 
            
            if only [1] meahs only turn to the positive lateral coor. 
            
            only [-1] means only turn to negative lateral coor. 
            
        @input: observation
            
            a series which represents the observation. THe index are either moments or the longitudinal coordinates. 
        
        @input: lw
        
            the lane width. 
        
        @input: YC_pre and YC_new
        
            YC_pre is the previsous lane middle line coordinate (which is the lane previsous of the last lane changing. )
            
            YC_new is the middle line of the lane before the detected lane changing (the lane should be also where the trajectory begins. )
        
        """
        import builtins
        
        #------------------------------------
        Zs = np.linspace(-lw/2, lw/2, int(N_discretize_Z))
        Zs = Zs[1:-1]
        #--------------------------------------------------------------
        #RETURNED VALUES
        #key are rhos and values are logprob
        #logprobs = {rho:{i:None for i in range(2, len(observation0))} for rho in rhos}
        #logprobs[rho] = (ks, logprob), both are list. 
        logprobs = {rho:{} for rho in rhos}
        
        observation0 = observation0.sort_index()
        #
        deltats = convert2sec*np.diff(observation0.index)
        Y_i_es = observation0.iloc[1:].values
        Y_i_1_es = observation0.iloc[:-1].values
        
        #--------------------------------------------------------------
        #LOGPROBS[rho][idx] = logprob
        LOGPROBS = {rho:[-np.inf for i in range(len(deltats))] for rho in rhos}
        #
        for idx in range(len(deltats)):
            #
            Y_i = observation0.iloc[idx+1]
            Y_i_1 = observation0.iloc[idx]
            #
            for rho0 in rhos:
                
                #
                rho = max(0, min(1, rho_old+rho0))
                
                logprobs0  = self.ConditionalY_given_previsous_Y_1(Y_i=Y_i, Y_i_1 = Y_i_1, deltat = deltats[idx], rho = rho, lw = lw, alpha = alpha, Zs = Zs, sde_args =sde_args, N_discretize_Z = N_discretize_Z)
                logprob = self.logprobs2meanlogprob(logprobs0)
                LOGPROBS[rho0][idx] = logprob
            #
            #determine which one is the greatest. 
            tmp = [LOGPROBS[rho][idx] for rho in rhos]
            rho_with_max = rhos[tmp.index(max(tmp))]
            rho_old = max(0, min(1, rho_old + rho_with_max))
            
            #print(rho_old)
            
        return LOGPROBS
            

        
        
        
        
        
        
        #########################################################
        import builtins
        
        #--------------------------------------------------------------
        #RETURNED VALUES
        #key are rhos and values are logprob
        #logprobs = {rho:{i:None for i in range(2, len(observation0))} for rho in rhos}
        #logprobs[rho] = (ks, logprob), both are list. 
        logprobs = {rho:{} for rho in rhos}
        
        observation0 = observation0.sort_index()
        #
        deltats = convert2sec*np.diff(observation0.index)
        Y_i_es = observation0.iloc[1:].values
        Y_i_1_es = observation0.iloc[:-1].values
        
        builtins.tmp = []
        
        #--------------------------------------------------------------
        #LOGPROBS[rho][YC] = (ks, logprobs)
        LOGPROBS = {rho:{YC:[] for YC in [YC_pre, YC_new]} for rho in rhos}
        for rho in rhos:
            for YC in [YC_pre, YC_new]:
                #
                YC_es = np.array([YC for i in range(1, len(observation0))])
                rhos4logprob = np.array([rho for i in range(1, len(observation0))])
                #
                ks = []
                logps = []
                #
                for k in range(0, len(observation0)-1):
                    #
                    if Y_i_1_es[k]==YC_es[k]:
                        ks.append(k)
                        logps.append(0)
                        continue
                    #
                    ks.append(k)
                    #
                    logps.append(self.logprobs(deltats = np.array([deltats[k]]), Y_i_es = np.array([Y_i_es[k]]), Y_i_1_es= np.array([Y_i_1_es[k]]), rhos = np.array([rhos4logprob[k]]), lw = lw, YC_es = np.array([YC_es[k]]), sigma = sigma, k = k))
                #
                LOGPROBS[rho][YC] = (ks, logps)
        
        #---------------------------------------------------------
        for rho in rhos:
            ks = []
            logstmp = []
            for k in range(0, len(observation0)-1):
                #
                ks.append(k)
                
                #
                tmp = sum(LOGPROBS[rho][YC_pre][1][:k]) + sum(LOGPROBS[rho][YC_new][1][k:])
                logstmp.append(tmp)
            logprobs[rho] = logstmp
        #
        return logprobs

    
    @classmethod
    def Conditional_given_Ys(self, Yi_es_series, convert2sec = .04, rho = 0, lw = 3.5, alpha = 3.0, Zs = False, sde_args = {'sigma':.01, 'thetaY':.2, 'thetaZ':.3, 'lw':3.5, 'target_Y':.0, 'v':1}, N_discretize_Z = 1e4):
        """
        Calculate the conditional distribution of Y_i (current value) given Y_i_1 (which is the previous value of Y).
        
        ---------------------------------------
        @input: N_discretize_Z
            the number of discretization of Z. 
            
            The discretized Z are within [-lw/2, lw/2]. 
        
        @input: Yi_es, deltats
        
            both are lists or .
        
        @input: Zs
            
            Zs is the discretization of Z, the noise. 
        
        @input: alpha
            
            the parameter used in erf scale. in meanY_sigmaY_given_previousYandZ_tanh_erf(alpha = alpha)
        
        @OUTPUT: LOGPROB
        
            the length should be len(Yi_es_series)-1.
            
            Each value is 
            
        -----------------------------------------
        @Step:
            - For each 
        """
        #
        if isinstance(Zs, bool):
            Zs = np.linspace(-lw/2, lw/2, int(N_discretize_Z))
            Zs = Zs[1:-1]
            
        #a list which means the log of p(Y1|Y0), .....p(Y_i|Y_i_1)
        #   it length should be len(Yi_es_series)-1
        LOGPROBS_series = []
        
        #
        deltats = convert2sec*np.diff(Yi_es_series.index)
        
        meanYs,sigmasquareYs = [],[]
        #Y_i is the newer value and Y_i_1 is the old value. 
        for Y_i,Y_i_1,deltat in zip(Yi_es_series.iloc[1:], Yi_es_series.iloc[:-1], deltats):
            #
            #meanys and sigmasquareys should be the same size as Zs
            meanys,sigmasquareys = self.meanY_sigmaY_given_previousYandZ_tanh_erf(Y_i = Y_i_1, Z_i_1 =  Zs, deltat = deltat, lw = lw, alpha = alpha,args_sde = sde_args, rho = rho)
            #
            #-----------------------logprob is a array, shape is the same as Zs
            logprob = log_normal(Y_i, meanys, np.sqrt(sigmasquareys))
            
            #
            LOGPROBS_series.append(self.logprobs2meanlogprob(logprob))
        
        #---------------
        
        
        
        return sum(LOGPROBS_series)
        
        
        #--------------------------------
        for Z in Zs:
            meany,sigmasquarey = self.meanY_sigmaY_given_previousYandZ_tanh_erf(Y_i = Y_i, Z_i_1 =  Z, deltat = deltat, lw = lw, alpha = alpha,args_sde = sde_args, rho = rho)
            logprob = log_normal(Y_i, meany, np.sqrt(sigmasquarey))
            
            logprobs.append(logprob)
            
        return logprobs

    @classmethod
    def kernel_max(self, data_1darray0, edges,percentage = 0.999, N = 20):
        """
        Given a 1d array. Find the max using kernel method. 
        
        @input: percentage
            
            only use the latter percentage data. 
        
        @input: N
            
            the discretization. 
        
        """
        data_1darray = copy.deepcopy(data_1darray0)
        
        
        # from scipy.stats import lognorm
        # data_1darray = data_1darray[-int(len(data_1darray)*percentage):]
        # s, loc, scale = lognorm.fit(data_1darray)
        # #estimated_mu = np.log(scale)
        # #estimated_sigma = s
        # return np.log(scale)
        
        
        from scipy import stats
        data_1darray = data_1darray[-int(len(data_1darray)*percentage):]
        kernel = stats.gaussian_kde(data_1darray)
        #
        evalutpoints = np.linspace(min(edges), max(edges), N)
        kernel_values = kernel(evalutpoints)
        return evalutpoints[kernel_values.argmax()]
        
        pass
    
    @classmethod
    def ConditionalY_given_previsous_Y_1(self, Y_i=1.0, Y_i_1 = 1.0, deltat = .5, rho = 0, lw = 3.5, alpha = 3.0, Zs = False, sde_args = {'sigma':.01, 'thetaY':.2, 'thetaZ':.3, 'lw':3.5, 'target_Y':.0, 'v':1}, N_discretize_Z = 1e4):
        """
        Calculate the conditional distribution of Y_i (current value) given Y_i_1 (which is the previous value of Y).
        
        ---------------------------------------
        @input: N_discretize_Z
            the number of discretization of Z. 
            
            The discretized Z are within [-lw/2, lw/2]. 
        
        @input: Y_i and Y_i_1
        
            Y_i is the current value
            Y_i_1 is the previous value. 
        
        @input: Zs
            
            Zs is the discretization of Z, the noise. 
        
        @input: alpha
            
            the parameter used in erf scale. in meanY_sigmaY_given_previousYandZ_tanh_erf(alpha = alpha)
        
        -----------------------------------------
        @Step:
            - For each 
        """
        #
        if isinstance(Zs, bool):
            Zs = np.linspace(-lw/2, lw/2, int(N_discretize_Z))
            Zs = Zs[1:-1]
        logprobs = []
        
        meanYs,sigmasquareYs = [],[]
        for Z in Zs:
            meany,sigmasquarey = self.meanY_sigmaY_given_previousYandZ_tanh_erf(Y_i = Y_i, Z_i_1 =  Z, deltat = deltat, lw = lw, alpha = alpha,args_sde = sde_args, rho = rho)
            logprob = log_normal(Y_i, meany, np.sqrt(sigmasquarey))
            logprobs.append(logprob)
        
        #a list. its shape is the same as Zs
        return logprobs


    @classmethod
    def meanZ_sigmaZ_given_previousZ_tanh(self, Z_i_1, deltat = .5, lw = 3.5, sde_args = {'sigma':.01, 'thetaY':.2, 'thetaZ':.3, 'lw':3.5, 'target_Y':.0, 'v':1}):
        """
        for tahn type 
        """
        #
        thetaY = args_sde['thetaY']
        thetaZ = args_sde['thetaZ']
        sigma = args_sde['sigma']
        v = args_sde['v']
        
        meanZ = None
        sigmasquareZ = None
        
        #   Z_i_1*np.pi/lw 
        #       np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)**2
        #meanY = -(Z_i_1 - lw*sigma*sigma*np.tan(Z_i_1*np.pi/lw )/(np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)**2)*deltat)*theta*deltat + (Y_i - theta*Y_i*deltat + rho*deltat*theta*lw)
        #meanY = Y_i + deltat*thetaY*Transformed_erf_integrator(Y_i,lw = lw, rho = rho) - deltat*thetaZ*(Z_i_1-deltat*2.0*Z_i_1*(sigma**2)*(lw**2 - 4.0*Z_i_1**2)/(2.0*lw**2))
        meanZ = Z_i_1-deltat*2.0*Z_i_1*(v**2)*(sigma**2)*(lw**2 - 4.0*Z_i_1**2)/(2.0*lw*lw)

        #sigmasquareY = (sigma*lw/(np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)))**2*deltat*((theta*deltat)**2)
        #sigmasquareY = deltat*((thetaZ*deltat)**2)*((sigma*lw)**2)*((lw**2 - 4.0*Z_i_1**2)/(2.0*lw*lw))**2
        sigmasquareZ = deltat*((sigma*lw*v)**2)*((lw**2 - 4.0*(Z_i_1**2))/(2.0*lw*lw))**2
        
        return meanZ,sigmasquareZ


    @classmethod
    def meanY_sigmaY_given_previousYandZ_tanh_erf(self, Y_i, Z_i_1, deltat = .5, lw = 3.5, alpha = 3.5,args_sde = {'thetaY':.5, 'thetaZ':.5, 'sigma':.01}, rho = 0):
        """
        This method return the mean and the sigma of Y_i+1 given Y_i and Z_i_2. 
        
        The parameters include sigma, thetaY and thetaZ. 
        
        DIfference:
        
            - self.meanY_sigmaY_given_previousYandZ()
            
            - self.meanY_sigmaY_given_previousYandZ_erf(), constructed from error fucntion. 
            
        NOTE:
            
            if Y_i, Z_i_1 are both float, then the returned value is also scalar. 
        
            While if Z_i_1 is an array, the returned meanY,sigmasquareY are both array. size the same as Z_i_1.
        
        ----------------------------------------
        @input: rho
            the lane changing flag
            
            If is either 0, or 1 or -1. 
            
        @input: Y_i, Z_i_1
            
            the input methods. 
        
            Z_i_1 can be an array. Then the returned value is also array with the same length. 
        
        @input: deltat
            the time interval in sec
        
        @input: lw
            lane width in meter. 
            
        
        """
        #
        thetaY = args_sde['thetaY']
        thetaZ = args_sde['thetaZ']
        sigma = args_sde['sigma']
        v = args_sde['v']
        
        meanY = None
        sigmasquareY = None
        
        #   Z_i_1*np.pi/lw 
        #       np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)**2
        #meanY = -(Z_i_1 - lw*sigma*sigma*np.tan(Z_i_1*np.pi/lw )/(np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)**2)*deltat)*theta*deltat + (Y_i - theta*Y_i*deltat + rho*deltat*theta*lw)
        meanY = Y_i + deltat*thetaY*Transformed_erf_integrator(Y_i,lw = lw, rho = rho, alpha = alpha) - deltat*thetaZ*(Z_i_1-deltat*2.0*Z_i_1*(sigma**2)*(v**2)*(lw**2 - 4.0*Z_i_1**2)/(2.0*lw**2))

        #sigmasquareY = (sigma*lw/(np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)))**2*deltat*((theta*deltat)**2)
        sigmasquareY = deltat*((thetaZ*deltat)**2)*((sigma*lw*v)**2)*((lw**2 - 4.0*Z_i_1**2)/(2.0*lw*lw))**2

        return meanY,sigmasquareY

    @classmethod
    def meanY_sigmaY_given_previousYandZ_erf(self, Y_i, Z_i_1, deltat = .5, lw = 3.5, args_sde = {'thetaY':.5, 'thetaZ':.5, 'sigma':.01, 'alpha':4.0}, rho = 0):
        """
        This method return the mean and the sigma of Y_i+1 given Y_i and Z_i_2. 
        
        The parameters include sigma, thetaY and thetaZ. 
        
        DIfference:
        
            - self.meanY_sigmaY_given_previousYandZ()
            
            - self.meanY_sigmaY_given_previousYandZ_erf(), constructed from error fucntion. 
            
        NOTE:
            
            if Y_i, Z_i_1 are both float, then the returned value is also scalar. 
        
            While if Z_i_1 is an array, the returned meanY,sigmasquareY are both array. size the same as Z_i_1.
        
        ----------------------------------------
        @input: rho
            the lane changing flag
            
            If is either 0, or 1 or -1. 
            
        @input: Y_i, Z_i_1
            
            the input methods. 
        
            Z_i_1 can be an array. Then the returned value is also array with the same length. 
        
        @input: deltat
            the time interval in sec
        
        @input: lw
            lane width in meter. 
            
        
        """
        #
        thetaY = args_sde['thetaY']
        thetaZ = args_sde['thetaZ']
        sigma = args_sde['sigma']
        alpha = args_sde['alpha']
        
        meanY = None
        sigmasquareY = None
        
        #   Z_i_1*np.pi/lw 
        #       np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)**2
        #meanY = -(Z_i_1 - lw*sigma*sigma*np.tan(Z_i_1*np.pi/lw )/(np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)**2)*deltat)*theta*deltat + (Y_i - theta*Y_i*deltat + rho*deltat*theta*lw)
        meanY = Y_i + deltat*thetaY*Transformed_erf_integrator(y = Y_i,alpha = alpha,lw = lw, rho = rho) - deltat*thetaZ*(Z_i_1-deltat*2.0*Z_i_1*(sigma**2)*(lw**2 - 4.0*Z_i_1**2)/(2.0*lw**2))

        #sigmasquareY = (sigma*lw/(np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)))**2*deltat*((theta*deltat)**2)
        sigmasquareY = deltat*((thetaZ*deltat)**2)*((sigma*lw)**2)*((lw**2 - 4.0*Z_i_1**2)/(2.0*lw*lw))**2

        return meanY,sigmasquareY

    @classmethod
    def meanY_sigmaY_given_previousYandZ(self, Y_i, Z_i_1, deltat = .5, lw = 3.5, args_sde = {'theta':.5, 'sigma':.01}, rho = 0):
        """
        THis method return the mean and the sigma of Y_i+1 given Y_i and Z_i_2. 
        
        NOTE:
            
            if Y_i, Z_i_1 are both float, then the returned value is also scalar. 
        
            While if Z_i_1 is an array, the returned meanY,sigmasquareY are both array. size the same as Z_i_1.
        
        ----------------------------------------
        @input: rho
            the lane changing flag
            
            If is either 0, or 1 or -1. 
            
        @input: Y_i, Z_i_1
            
            the input methods. 
        
            Z_i_1 can be an array. Then the returned value is also array with the same length. 
        
        @input: deltat
            the time interval in sec
        
        @input: lw
            lane width in meter. 
            
        
        """
        #
        theta = args_sde['theta']
        sigma = args_sde['sigma']
        
        meanY = None
        sigmasquareY = None
        
        #   Z_i_1*np.pi/lw 
        #       np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)**2
        meanY = -(Z_i_1 - lw*sigma*sigma*np.tan(Z_i_1*np.pi/lw )/(np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)**2)*deltat)*theta*deltat + (Y_i - theta*Y_i*deltat + rho*deltat*theta*lw)
        
        sigmasquareY = (sigma*lw/(np.pi*(np.tan(Z_i_1*np.pi/lw )**2 + 1)))**2*deltat*((theta*deltat)**2)

        return meanY,sigmasquareY
    
    @classmethod
    def prob_given_observationseries(self, observation_series, args_sde = {'theta':.5, 'sigma':.01}, multiplier = 1e10):
        """
        Compute the probability given that the transition moment occurs at certain interval, i.e. transition_t
        ---------------------------------------------------------------------
        @input: theta and sigma
            
            the parameters of the SDE. 
        
        @input: observation_series
            a series. 
            
            The index are moments or the longitudinal coordinates. 
            
            The values are the lateral displacement 
        
        @input: transition_t
            
            the transition moment. 
        
        @input: multiplier
            
            to prevent the overflow of the probability value. 
        
        @OUTUPUT: log_prob
            
             a float. THe log of the probability
        
        
        """
        #
        moments = sorted(observation_series.index)
        delta_ts = np.diff(moments)
        
        #
        
        
        pass
    
    
    pass



class LC_identification():
    """
    find the transition moment of the lane changing process using the SDE routine. 
    
    """
    

    
    pass



class SDECalirbration():
    """
    calibration of the parameters in the sde. 
    
    The parameters that need to be calibrated include \theta and \rho
    """
    
    

    
    
    pass

def sample_Z(observation, lw = 3.5):
    """
    @input: lw is the lane width. 
    
    
    return a np.array 1d. 
    
    length is len(observation)-1
    
    #np.random.uniform(low=0.0, high=1.0, size=None)
    """
    
    
    return np.random.uniform(lw = -lw/2.0, high =  lw/2.0, size = (len(observation)-1 ,))
    

#import mcint
def objective4optimize(theta, sigma, observation, lw = 3.5,  n_mcint = 1e2):
    """
    @input: sigma
        0.0<sigma<=1.0
        
    @input: theta
        
        0<theta<10
        
    @input: observation
        
        a pd.Series. THe index should be time and unit is sec. 
        value unit is meter. 
    
    @input: n_mcint
        points for the 
    """
    
    
    #logprob_given_single_series4integration(self, Z, observation, lw = 3.5,  args_sde = {'theta':.5, 'sigma':.01})
    #result, error = mcint.integrate(integrand, sampler(), measure=math.pi/4)
    result, error = mcint.integrate(sde_util.logprob_given_single_series4integration(observation = observation, lw = lw, args_sde = {'theta':theta, 'sigma':sigma}), sample_Z(observation = observation), n = n_mcint)
    return result,error
    pass


def prepare_Zs_logprobs_Zs(observation, lw = 3.5, n_mcint = 1e5,):
    """
    
    """
    
    #
    Zs = np.random.uniform(low = -lw/2.0, high = lw/2.0, size = (len(observation),int(n_mcint)))
    #
    #meanZ and sigmaZ shape are both (len(observation)-1, n_mcint)
    meanZ,sigmaZ = sde_util.meanZ_sigmaZ_given_previousZ_erf(Zs[:-1, :])
    #logprobs_Zs.shape is (len(observation)-1, n_mcint)
    logprobs_Zs = log_normal(Zs[1:,:], meanZ,sigmaZ)
    
    return Zs,logprobs_Zs



def EquilibriumPoints_tanh_erf_SDE_guarantee_stability(lw = 3.5, sde_args = {'theta_Z':1.0, 'theta_Y':.4, 'sigma':.3, 'alpha':7.0}, N_discretize = 1000, tolerance = 1e-2):
    """
    Find the equilibrium points in the stability analysis. 
    ---------------------------------------------------------------
    @input: lw
    
        the lane width, unit is meter. 
        
    @input： tolerance
        
        the tolrance of the error. 
    
    @OUTPUT: success,equilibriumpoints
        
        success is a bool. True means there is equilibrium points while False means there is no equilibrium points. 
        
        equilibriumpoints is a list of Y, which makes Transformed_erf_integrator() equals either theta_Z/theta_Y*lw/2.0 or -theta_Z/theta_Y*lw/2.0. 
        
    """
    equilibriumpoints = []
    
    
    theta_Z = sde_args['theta_Z']
    theta_Y = sde_args['theta_Y']
    sigma = sde_args['sigma']
    alpha = sde_args['alpha']
    
    #-------------------------------
    Ys = np.linspace(-lw/2.0, lw/2.0, N_discretize)
    
    #-------------------------------
    #differences=0 is the solution of equilibrium point. 
    #  There are two equilibrium points of Z: lw/2 and -lw/2. 
    #first case: lw/2
    differences =  np.array([Transformed_erf_integrator(y=y, lw = lw, alpha = alpha)-theta_Z/theta_Y*lw/2.0 for y in Ys])
    #if 
    if sum(np.abs(differences)<tolerance)>0:
        idxs = np.where(differences<=tolerance)[0]
        equilibriumpoints.extend([Ys[i] for i in idxs])
    #second case: -lw/2
    differences =  np.array([Transformed_erf_integrator(y=y, lw = lw, alpha = alpha)+theta_Z/theta_Y*lw/2.0 for y in Ys])
    #if 
    if sum(np.abs(differences)<tolerance)>0:
        idxs = np.where(differences<=tolerance)[0]
        equilibriumpoints.extend([Ys[i] for i in idxs])
    
    if len(equilibriumpoints)==0:
        return False,[]
    else:
        return True,equilibriumpoints
    

def PSO_fun_by_pyswarm_calllback(particles, observation0 = 0, sigma=.05, theta_Y = 5.0, theta_Z = .05, alpha = 1.0, lw = 3.5,  n_mcint = 1e5, Zs = False, logprobs_Zs = False, M_penelty = 1e100, convert2sec = .1, limit_observation_len = 100):
    """
    
    PSO callbacK:
    
        
        from pyswarms.single.global_best import GlobalBestPSO

        # instatiate the optimizer
        x_max = 10 * np.ones(2)
        x_min = -1 * x_max
        bounds = (x_min, x_max)
        options = {'c1': 0.5, 'c2': 0.3, 'w': 0.9}
        optimizer = GlobalBestPSO(n_particles=10, dimensions=2, options=options, bounds=bounds)

        # now run the optimization, pass a=1 and b=100 as a tuple assigned to args

        cost, pos = optimizer.optimize(rosenbrock_with_args, 1000, a=1, b=100, c=0)
    ------------------------------------------------
    @input: particles
    
        shape is particles_N, dim.
        
        shoule be (N, 3)
        
        
    """
    res_iter = []
    for dim in range(particles.shape[0]):
        sigma = particles[dim, 0]
        theta_Y = particles[dim, 1]
        theta_Z = particles[dim, 2]
        res_iter.append(PSO_objective_maximize_tanh_erf(observation0 = observation0, sigma = sigma, theta_Y = theta_Y, theta_Z = theta_Z, alpha = alpha, lw = lw,  n_mcint = n_mcint, Zs = Zs, logprobs_Zs = logprobs_Zs, M_penelty = M_penelty, convert2sec = convert2sec, limit_observation_len = limit_observation_len))
        
    return np.array(res_iter)
        


def PSO_objective_maximize_tanh_erf(observation0, sigma=.05, theta_Y = 5.0, theta_Z = .05, alpha = 1.0, lw = 3.5,  n_mcint = 1e5, Zs = False, logprobs_Zs = False, M_penelty = 1e100, convert2sec = .1, limit_observation_len = 100):
    """
    the objective function for the pso optimiztion. 
    
    The objective is to maximize the value. 
    
    ---------------------------------------------------------------
    @input: sigma
        0.0<sigma<=.1
        
    @input: theta
        
        0<thetaY<6
        0<thetaZ<.3
        
        
    @input: observation
        
        a pd.Series. THe index should be time and unit is sec. 
        
        value unit is meter. 
    
    @input: n_mcint
        points for the 
    
    @iput: M_penelty
        
        penelty for stability condition violations. 
        
        Lyapunov function should be negative. If it is positive, then the penelty is -M_penelty.
    
    
    @OUTPUT: objective_value
    
        it is the sum of logprob of observation and penelty (for stablity. )
    
    
    """
    #if theta_Z>theta_Y/20.0:
    #    return -M_penelty
    
    #-----------First term: the true objective value
    #--------
    obj_value = OBJECTIVE_omit_smaller_number4optimize_mc_erf(observation0 = observation0, sigma= sigma, thetaY = theta_Y, thetaZ = theta_Z, alpha = alpha, lw = lw,  n_mcint = n_mcint, Zs = Zs, logprobs_Zs = logprobs_Zs, convert2sec = convert2sec, limit_observation_len = limit_observation_len)
    #-----------Second term: the penelty
    
    penelty = lyapunov_function_tanh_erf_SDE_guarantee_stability(lw = lw, sde_args = {'theta_Z':theta_Z, 'theta_Y':theta_Y, 'sigma':sigma, 'alpha':alpha}, M = M_penelty)
    #
    return obj_value+penelty


def lyapunov_function_tanh_erf_SDE_guarantee_stability(lw = 3.5, sde_args = {'theta_Z':1.0, 'theta_Y':.4, 'sigma':.3, 'alpha':7.0}, M = 1e50, N_discretize_equilibrium_numerical = 1000, tolerance = 1e-2):
    """
    the penelty of the PSO objective. The objective of the method is to guarantee the stability. 
    
    The returned value is the lyapunov function value. IT SHOULD BE SMALLER OR EQUAL THAN 0, to guarantee the stability. 
    
    As the requirement of the stability is lyapunov<=0, thus if it is positive, return -M*lyapunov.
    
    ------------------------------------------------------
    @input: 
        
        theta_Z = 1.0, sigma = .4, alpha = 7.0
    
    @OUTUT:
        
        lyfunctionvalue.
        
        It may be either -M (M is one input arg), or the calculated value. 
        
        Returning M means that the system is not stable given the input args of sde. 
        
    """
    Ys_Zs_meshgrid = np.meshgrid(np.linspace(-lw/2.0, lw/2.0, 100), np.linspace(-lw/2.0, lw/2.0, 100))
    #discretize the Y (lateral coordinate) and Z (lateral noise).
    #Ys = np.linspace(-lw/2.0, lw/2.0, 100);Zs = np.linspace(-lw/2.0, lw/2.0, 100)
    #Ys_grid and Zs_grid are both 2d array. 
    Ys_grid,Zs_grid = Ys_Zs_meshgrid
    
    
    #
    theta_Z = sde_args['theta_Z']
    theta_Y = sde_args['theta_Y']
    sigma = sde_args['sigma']
    alpha = sde_args['alpha']
    
    #--------------Find all the equilibrimp poinsts
    has_equilibrium_points,equilibriumpoints = EquilibriumPoints_tanh_erf_SDE_guarantee_stability(lw = lw, sde_args = sde_args, N_discretize = N_discretize_equilibrium_numerical, tolerance = tolerance)
    #
    if not has_equilibrium_points:
        return -M
    
    stable_or_not = True
    for equilibrium_y in equilibriumpoints:
        term1 = 2*(Ys_grid-equilibrium_y)
        term2 = theta_Y*Transformed_erf_integrator(Ys_grid, lw = lw, alpha = alpha)-theta_Z*Zs_grid
        term3  = -(2*Zs_grid-lw)*Zs_grid*sigma*sigma*(lw*lw - 4*np.multiply(Zs_grid,Zs_grid))/(lw*lw) + sigma*sigma*(lw*lw - 4*np.multiply(Zs_grid,Zs_grid))/4/(lw*lw)
        
        res = np.multiply(term1,term2) + term3
        #if all res negative, then it is stable. 
        
        if res.max()>0:
            stable_or_not = False
            break
        
    #
    if stable_or_not:
        return 0.0
    else:
        return -M


def OBJECTIVE_omit_smaller_number4optimize_mc_erf(observation0, sigma=.08, thetaY = .05, thetaZ = .05, lw = 3.5, alpha = 4.0,  n_mcint = 1e5, Zs = False, logprobs_Zs = False, convert2sec = .1, limit_observation_len = 200):
    """
    The objective function for single objervatio. 
    DIfference between:
        - objective_omit_smaller_number4optimize_mc_erf()
        - objective_omit_smaller_number4optimize_mc()
    
    
    There are three parameters that need to be calibrated: sigma, thetaY, thetaZ. 
    
        sigma domain (0, .5]
        thetaY domain (0, 4)
        thetaZ domain (0, 2)
    
    
    -----------------------------------------------
    @input: alpha
        
        the arg used in the transformed error function. 
        
    @input: limit_observation_len
    
        if the observation is too long, then just use the last part of the observation. 
        
        This arg is used to limit the observations in log prob. 
        
    @input: sigma
        0.0<sigma<=1.0
        
    @input: theta
        
        0<theta<10
        
    @input: observation
        
        a pd.Series. THe index should be time and unit is sec. 
        value unit is meter. 
    
    @input: n_mcint
        points for the 
    
    @OUTPUT: digitalsN,averageProb
        
        digitalsN is an integer. It means the number of zeros. 
        
        averageProb is a prob. 
        
    """
    #-------------------------------------
    if len(observation0)>=limit_observation_len:
        observation = observation0.iloc[-limit_observation_len:]
    else:
        observation = observation0
    
    #------------random sample the Z
    #   Zs shape is (len(observation)-1, n_mcint)
    if isinstance(Zs, bool):
        #
        Zs = np.random.uniform(low = -lw/2.0, high = lw/2.0, size = (len(observation),int(n_mcint)))
        #
        #meanZ and sigmaZ shape are both (len(observation)-1, n_mcint)
        meanZ,sigmaZ = sde_util.meanZ_sigmaZ_given_previousZ_erf(Zs[:-1, :])
        #logprobs_Zs.shape is (len(observation)-1, n_mcint)
        logprobs_Zs = log_normal(Zs[1:,:], meanZ,sigmaZ)
    
    #------------res.shape is the same as Zs, i.e. (len(observation)-1, n_mcint)
    #res.sum(0).shape is  n_mcint
    logprobs_Ys = sde_util.LOGPROBS_given_single_seriesmultiple_Z_4integration_erf(Zs, observation, args_sde = {'thetaY':thetaY, 'thetaZ':thetaZ,  'sigma':sigma, 'alpha':alpha}, convert2sec = convert2sec)
    
    #return res
    
    #-------------The reason for the concat is that p(y) = p(z_i | z_i_1) p( y_i | z_i), thus the log of Ys and Zs need to be concated. 
    #shape is (2*(len(observation)-1), n_mcint)
    #print(Zs.shape, logprobs_Ys.shape, len(observation), type(logprobs_Ys))
    res = np.concatenate((logprobs_Zs, logprobs_Ys))
    #print(sample_logs.shape)
    
    sample_logs = np.array(sorted(res.sum(0)))
    #if sum(np.isinf([np.inf, -np.inf]))>0:
    #    pass
    
    digitsMultipled = np.log10(-max(sample_logs))
    #sample_logs_added = [i for i in 10**(digitsMultipled)+sample_logs]
    sample_logs_added = -sample_logs[-1]+sample_logs
    
    
    return -10**digitsMultipled+np.log10(np.mean(np.exp(sample_logs_added)))
    # #########################################################
    # #digitsNs is the length of zeros in the scientific notation. 
    # #   np.format_float_scientific() return a str
    # try:
        # digitsNs = []
        # prefix_values = []
        # for i in sample_logs:
            # if i>-1:
                # digitsNs.append(0)
                # prefix_values.append(0)
            # else:
                # digitsNs.append(int(np.format_float_scientific(i).split('+')[1]))
                # prefix_values.append(float(np.format_float_scientific(i).split('+')[0][:-1]))
    # except:
        # pickled_data = {'thetaY':thetaY, 'thetaZ':thetaZ,  'sigma':sigma},res,observation,Zs
        # pickle.dump(pickled_data, open('pickled_data.pickle', 'wb'))
        # raise ValueError('dafasfasadfdafasfasadfdafasfasadfdafasfasadfdafasfasadfdafasfasadf')
    # #digitsNs = [int(np.format_float_scientific(i).split('+')[1]) for i in sample_logs]
    
    # #
    # digitsMIN = min(digitsNs)
    # prefix_digitsMIN =  prefix_values[digitsNs.index(digitsMIN)]
    
    # #
    # #the condition >=-10 is to omit the smaller number. 
    # #sample_logs_added = [i for i in 10**(digitsMIN)+sample_logs if i>=-10]
    # sample_logs_added = [i for i in 10**(digitsMIN)+sample_logs]
    
    
    # #
    # return -10**digitsMIN+np.log(np.mean(np.exp(sample_logs_added)))


def objective_omit_smaller_number4optimize_mc_erf(observation, sigma=.05, thetaY = .05, thetaZ = .05, lw = 3.5,  n_mcint = 1e5):
    """
    The objective function for single objervatio. 
    DIfference between:
        - objective_omit_smaller_number4optimize_mc_erf()
        - objective_omit_smaller_number4optimize_mc()
    
    
    There are three parameters that need to be calibrated: sigma, thetaY, thetaZ. 
    
        sigma domain (0, 1]
        thetaY domain (0, 2)
        thetaZ domain (0, 2)
    
    
    
    -----------------------------------------------
    @input: sigma
        0.0<sigma<=1.0
        
    @input: theta
        
        0<theta<10
        
    @input: observation
        
        a pd.Series. THe index should be time and unit is sec. 
        value unit is meter. 
    
    @input: n_mcint
        points for the 
    
    @OUTPUT: digitalsN,averageProb
        
        digitalsN is an integer. It means the number of zeros. 
        
        averageProb is a prob. 
        
    """
    
    
    #------------random sample the Z
    #   Zs shape is (len(observation)-1, n_mcint)
    Zs = np.random.uniform(low = -lw/2.0, high = lw/2.0, size = (len(observation)-1,int(n_mcint)))
    
    #------------res.shape is the same as Zs, i.e. (len(observation)-1, n_mcint)
    #res.sum(0).shape is  n_mcint
    res = sde_util.logprobs_given_single_seriesmultiple_Z_4integration_erf(Zs, observation, args_sde = {'thetaY':thetaY, 'thetaZ':thetaZ,  'sigma':sigma})
    
    #the value of the log(p) for certain Z. 
    #shape is n_mcint
    sample_logs = np.array(sorted(res.sum(0)))
    #if sum(np.isinf([np.inf, -np.inf]))>0:
    #    pass
    
    #digitsNs is the length of zeros in the scientific notation. 
    #   np.format_float_scientific() return a str
    try:
        digitsNs = []
        for i in sample_logs:
            if i>-1:
                digitsNs.append(1)
            else:
                digitsNs.append(int(np.format_float_scientific(i).split('+')[1]))
    except:
        pickled_data = {'thetaY':thetaY, 'thetaZ':thetaZ,  'sigma':sigma},res,observation,Zs
        pickle.dump(pickled_data, open('pickled_data.pickle', 'wb'))
        raise ValueError('dafasfasadfdafasfasadfdafasfasadfdafasfasadfdafasfasadfdafasfasadf')
    #digitsNs = [int(np.format_float_scientific(i).split('+')[1]) for i in sample_logs]
    digitsMIN = min(digitsNs)
    
    #the condition >=-10 is to omit the smaller number. 
    sample_logs_added = [i for i in 10**(digitsMIN)+sample_logs if i>=-10]
    
    #
    return -10**digitsMIN+np.log(np.mean(np.exp(sample_logs_added)))

def DifferentiateBatch(observations_dict, limit_observation_len= 100, convert2sec = 1.0, n_mcint = 1e6):
    """
    @input: observations_dict
    
        a dict. The keys are 
    @input: limit_observation_len
    
        If the observation length is too long, just use part of it to compute the logprob. 
    
    
    @OUTPUT: logav_versus_logprobhdv
    
        a list. 
        
        Each element is calculated as logprob_av/logprob_hdv
        
        logprob_av is the log of the probability that it is a av. 
        
        logprob_hdv is the log of the probability that it is a hdv.  
        
    """
    logav_versus_logprobhdv = []
    for vid in sorted(observations_dict.keys()):
        #
        res = Differentiate(observation0  = observations_dict[vid], limit_observation_len= limit_observation_len, convert2sec = convert2sec, n_mcint = n_mcint)
        #
        logav_versus_logprobhdv.append(res)
    #
    return np.array(logav_versus_logprobhdv)
    #

def Differentiate(observation0, limit_observation_len= 100, convert2sec = 0.02, n_mcint =  1e5):
    """
        
    @input: observation
        
        a pd.Series. THe index should be time and unit is sec. 
        value unit is meter. 
    
    @input: limit_observation_len
    
        
        
    """
    #
    logprobcav= OBJECTIVE_omit_smaller_number4optimize_mc_erf(observation0 = observation0, limit_observation_len= limit_observation_len,  sigma=.05, thetaY = 3.35, thetaZ = .17, convert2sec = convert2sec, n_mcint = n_mcint)
    #
    logprobhdv = OBJECTIVE_omit_smaller_number4optimize_mc_erf(observation0 = observation0, limit_observation_len= limit_observation_len, sigma=.05,thetaY = 1.2, thetaZ = .14, convert2sec = convert2sec, n_mcint = n_mcint)
    #
    return logprobcav/logprobhdv


#import mcint
def OBJECTIVE_omit_smaller_number4optimize_mc(observation0, theta = .1, sigma=.05, lw = 3.5,  n_mcint = 1e5, Zs = False):
    """
    As the prob is small, the smaller ones are omitted. 
    
    -----------------------------------------------
    @input: sigma
        0.0<sigma<=1.0
        
    @input: theta
        
        0<theta<10
        
    @input: observation
        
        a pd.Series. THe index should be time and unit is sec. 
        value unit is meter. 
    
    @input: n_mcint
        points for the 
    
    @OUTPUT: digitalsN,averageProb
        
        digitalsN is an integer. It means the number of zeros. 
        
        averageProb is a prob. 
        
    """
    #-------------------------------------
    if len(observation0)>=200:
        observation = observation0.iloc[-200:]
    else:
        observation = observation0
    
    #------------random sample the Z
    #   Zs shape is (len(observation)-1, n_mcint)
    if Zs==False:
        #
        Zs = np.random.uniform(low = -lw/2.0, high = lw/2.0, size = (len(observation),int(n_mcint)))
        #
        #meanZ and sigmaZ shape are both (len(observation)-1, n_mcint)
        meanZ,sigmaZ = sde_util.meanZ_sigmaZ_given_previousZ_erf(Zs[:-1, :])
        #logprobs_Zs.shape is (len(observation)-1, n_mcint)
        logprobs_Zs = log_normal(Zs[1:,:], meanZ,sigmaZ)
    
    #------------res.shape is the same as Zs, i.e. (len(observation)-1, n_mcint)
    #res.sum(0).shape is  n_mcint
    logprobs_Ys = sde_util.logprobs_given_single_seriesmultiple_Z_4integration(Zs, observation, args_sde = {'theta':theta, 'sigma':sigma})
    
    #-------------The reason for the concat is that p(y) = p(z_i | z_i_1) p( y_i | z_i), thus the log of Ys and Zs need to be concated. 
    #shape is (2*(len(observation)-1), n_mcint)
    res = np.concatenate((logprobs_Zs, logprobs_Ys))

    sample_logs = np.array(sorted(res.sum(0)))
    #if sum(np.isinf([np.inf, -np.inf]))>0:
    #    pass
    
    digitsMultipled = np.log10(-max(sample_logs))
    #sample_logs_added = [i for i in 10**(digitsMultipled)+sample_logs]
    sample_logs_added = -sample_logs[-1]+sample_logs
    
    return -10**digitsMultipled+np.log(np.mean(np.exp(sample_logs_added)))
    
    
    
    
    
    
    
    
    #digitsNs is the length of zeros in the scientific notation. 
    #   np.format_float_scientific() return a str
    digitsNs = [int(np.format_float_scientific(i).split('+')[1]) for i in sample_logs]
    digitsMIN = min(digitsNs)
    
    #the condition >=-10 is to omit the smaller number. 
    sample_logs_added = [i for i in 10**(-digitsMIN)*sample_logs if i>=-10]
    
    
    #
    return -10**digitsMIN+np.log(np.mean(np.exp(sample_logs_added)))
    
    return digitsMIN,np.mean(np.exp(sample_logs_added))
    return sample_logs
    #-------------------
    return multiplier+res.sum(0)
    return np.mean(np.exp(multiplier+res.sum(0)))
    
    #logprob_given_single_series4integration(self, Z, observation, lw = 3.5,  args_sde = {'theta':.5, 'sigma':.01})
    #result, error = mcint.integrate(integrand, sampler(), measure=math.pi/4)
    result, error = mcint.integrate(sde_util.logprob_given_single_series4integration(observation = observation, lw = lw, args_sde = {'theta':theta, 'sigma':sigma}), sample_Z(observation = observation), n = n_mcint)
    return result,error
    pass


#import mcint
def objective_omit_smaller_number4optimize_mc(observation, theta = .1, sigma=.05, lw = 3.5,  n_mcint = 1e5, Zs = False):
    """
    As the prob is small, the smaller ones are omitted. 
    
    -----------------------------------------------
    @input: sigma
        0.0<sigma<=1.0
        
    @input: theta
        
        0<theta<10
        
    @input: observation
        
        a pd.Series. THe index should be time and unit is sec. 
        value unit is meter. 
    
    @input: n_mcint
        points for the 
    
    @OUTPUT: digitalsN,averageProb
        
        digitalsN is an integer. It means the number of zeros. 
        
        averageProb is a prob. 
        
    """
    
    
    #------------random sample the Z
    #   Zs shape is (len(observation)-1, n_mcint)
    if Zs==False:
        Zs = np.random.uniform(low = -lw/2.0, high = lw/2.0, size = (len(observation)-1,int(n_mcint)))
    
    #------------res.shape is the same as Zs, i.e. (len(observation)-1, n_mcint)
    #res.sum(0).shape is  n_mcint
    res = sde_util.logprobs_given_single_seriesmultiple_Z_4integration(Zs, observation, args_sde = {'theta':theta, 'sigma':sigma})
    
    #the value of the log(p) for certain Z. 
    sample_logs = np.array(sorted(res.sum(0)))
    
    #digitsNs is the length of zeros in the scientific notation. 
    #   np.format_float_scientific() return a str
    digitsNs = [int(np.format_float_scientific(i).split('+')[1]) for i in sample_logs]
    digitsMIN = min(digitsNs)
    
    #the condition >=-10 is to omit the smaller number. 
    sample_logs_added = [i for i in 10**(-digitsMIN)*sample_logs if i>=-10]
    
    
    #
    return -10**digitsMIN+np.log(np.mean(np.exp(sample_logs_added)))
    
    return digitsMIN,np.mean(np.exp(sample_logs_added))
    return sample_logs
    #-------------------
    return multiplier+res.sum(0)
    return np.mean(np.exp(multiplier+res.sum(0)))
    
    #logprob_given_single_series4integration(self, Z, observation, lw = 3.5,  args_sde = {'theta':.5, 'sigma':.01})
    #result, error = mcint.integrate(integrand, sampler(), measure=math.pi/4)
    result, error = mcint.integrate(sde_util.logprob_given_single_series4integration(observation = observation, lw = lw, args_sde = {'theta':theta, 'sigma':sigma}), sample_Z(observation = observation), n = n_mcint)
    return result,error
    pass


def normal(x, mean, sigma):
    """
    calculate the normal distribution value at x, given the mean and sigma. 
    """
    sigma = abs(sigma)
    return 1.0/(np.sqrt(2.0*np.pi)* sigma)*np.exp(-(x-mean)**2/(2.0*sigma*sigma))


def log_normal(x, mean, sigma):
    """
    calculate the normal distribution value at x, given the mean and sigma. 
    """
    sigma = abs(sigma)
    
    tmp1 = -np.sqrt(2.0*np.pi)* sigma
    tmp2 = -(x-mean)**2/(2.0*sigma*sigma)
    
    return tmp1+tmp2
    
    return 1.0/(np.sqrt(2.0*np.pi)* sigma)*np.exp(-(x-mean)**2/(2.0*sigma*sigma))













