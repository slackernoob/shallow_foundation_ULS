import math
import csv

'''
UNDRAINED DESIGN BEARING CAPACITY:
Rd / A' = 5.14 * cu * sc * ic * bc + q
'''

# To get help, do help(FUNCTION)
# E.g. help(UndrainedCalc().q))

class UndrainedCalc:
    '''
    A class for Undrained calculations.
    '''
    
    def spt_N_value(self, D):
        '''
        Calculate SPT-N value using linear relationship obtained from borehole data

        Parameters:    
            D (float): depth of foundation base (in m)

        Returns:
            float: SPT-N value
        '''
        if D <= 3.225: # Lower Limit
            return 9
        elif D <= 12.225: 
            return 4.1317 * D - 0.1076
        else:
            return 50 # Upper Limit

    def q(self, D, gamma=20):
        '''
        Calculate total overburden pressure at level of foundation base

        Parameters:    
            gamma (float): bulk unit weight of soil, assuming constant 20kN/m^3
            D (float): depth of foundation base (in m)

        Returns:
            float: total overburden pressure
        '''
        return gamma * D
    
    def ic(self, H, A_prime, cu):
        '''
        Calculate inclination of loading due to horizontal load H

        Parameters:    
            H (float): bulk unit weight of soil, assuming constant 20kN/m^3
            A_prime (float): effective area (in m^2)
            cu (float): undrained shear strength (in kN/m^2)

        Returns:
            float: inclination of loading
        '''
        assert A_prime > 0 
        assert cu > 0
        return 0.5 * (1 + (1 - H / ((A_prime * cu)) ** 0.5) )

    def sc(self, B_prime, L_prime):
        '''
        Calculate shape of foundation (rectangular/square)

        Parameters:
            B_prime (float): effective width (in m)
            L_prime (float): effective length (in m)

        Returns:
            float: shape of foundation
        '''
        assert B_prime > 0
        assert L_prime > 0
        return 1 + 0.2 * (B_prime / L_prime)
    
    def bc(self, alpha):
        '''
        Calculate inclination of foundation base

        Parameters:
            alpha (float): angle of inclination (in degrees)

        Returns:
            float: inclination of foundation base
        '''
        return 1 - (2 * alpha) / (2 + math.pi)

    def dc(self, D, B_prime):
        '''
        Calculate depth of embedment

        Parameters:
            D (float): depth of foundation base (in m)
            B_prime (float): effective width (in m)

        Returns:
            float: depth of embedment
        '''
        assert D > 0
        assert B_prime > 0
        res = 1 + 0.35 * (D / B_prime)
        return res if res <= 1.7 else 1.7
    
    def volume(self, B, L):
        '''
        Calculate volume of foundation base

        Parameters:
            B (float): width of foundation base (in m)
            L (float): length of foundation base (in m)

        Returns:
            float: volume of foundation base
        '''
        return B * L * 0.5
    
    def price(self, D, volume):
        '''
        Calculate cost of concrete for foundation base

        Parameters:
            D (float): depth of foundation base (in m)
            volume (float): volume of foundation base (in m^3)

        Returns:
            float: volume of foundation base
        '''
        base_price = volume * 100000
        if D <= 1:
            return base_price
        elif D <= 2:
            return base_price * 1.05
        elif D <= 3:
            return base_price * 1.1
        else:
            return base_price * 1.15

'''
DRAINED DESIGN BEARING CAPACITY:
Rd / A' = c' * Nc * sc * ic * bc * dc + q' * Nq * sq * iq * bq * dq + 0.5 gamma' * B'
'''

class DrainedCalc:
    def q_prime(self, D, density, water_table_depth):
        '''
        Calculate design effective overburden pressure at level of foundation base

        Parameters:    
            D (float): depth of foundation base (in m)
            density (float): bulk density of soil (in kN/m^3)
            water_table_depth (float): depth of water table (in m)

        Returns:
            float: total overburden pressure
        '''
        if D <= water_table_depth:
            return D * density
        else:
            return water_table_depth * density + (D - water_table_depth) * (density - 10)
    
    def eb(self, H, V, h=0): # H refers to force acting in direction of B, h is distance from H to bottom left corner
        '''
        Calculate eccentricity when force acting in direction of width

        Parameters:    
            H (float): horizontal force acting in direction of width (in kN)
            V (float): vertical force (in kN)
            h (float): moment arm of H (in m)

        Returns:
            float: eccentricity (in m)
        '''
        return H * h / V
    
    def el(self, H, V, h=0): # H refers to force acting in direction of L, h is distance from H to bottom left corner
        '''
        Calculate eccentricity when force acting in direction of length

        Parameters:    
            H (float): horizonal force acting in direction of length (in kN)
            V (float): vertical force (in kN)
            h (float): moment arm of H (in m)

        Returns:
            float: eccentricity (in m)
        '''
        return H * h / V

    def B_prime(self, B, eb):
        '''
        Calculate effective width

        Parameters:    
            B (float): width of foundation base (in m)
            eb (float): eccentricity (in m)

        Returns:
            float: effective width (in m)
        '''
        return 2 * (B/2 - eb)
    
    def L_prime(self, L, el):
        '''
        Calculate effective length

        Parameters:    
            L (float): length of foundation base (in m)
            el (float): eccentricity (in m)

        Returns:
            float: effective length (in m)
        '''
        return 2 * (L/2 - el)
    
    def A_prime(self, B_prime, L_prime):
        '''
        Calculate effective area

        Parameters:    
            B_prime (float): effective width (in m)
            L_prime (float): effective length (in m)

        Returns:
            float: effective area (in m^2)
        '''
        return B_prime * L_prime
    
    def Wb_prime(self, D, density, water_table_depth, thickness, B, L):
        '''
        Calculate weight of backfill material

        Parameters:    
            D (float): depth of foundation base (in m)
            density (float): bulk density of soil (in kN/m^3)
            water_table_depth (float): depth of water table (in m)
            B (float): width of foundation base (in m)
            L (float): length of foundation base (in m)

        Returns:
            float: weight of backfill material (in kN)
        '''
        if D <= thickness:
            return 0
        else:
            if (D - thickness) <= water_table_depth:
                return B * L * (D - thickness) * density
            else:
                return B * L * (D - thickness - water_table_depth) * (density - 10) + B * L * (water_table_depth) * density
    
    def Wf(self, B, L, thickness, unit_weight):
        '''
        Calculate weight of foundation base

        Parameters:    
            B (float): width of foundation base (in m)
            L (float): length of foundation base (in m)
            thickness (float): thickness of foundation base (in m)
            unit_weight (float): unit weight of foundation base (in kN/m^3)

        Returns:
            float: weight of foundation base (in kN)
        '''
        return B * L * thickness * unit_weight
    
    def U1(self, D, B, L, water_table_depth, thickness):
        '''
        Calculate force due to water pressure on top of foundation base

        Parameters:    
            D (float): depth of foundation base (in m)
            B (float): width of foundation base (in m)
            L (float): length of foundation base (in m)
            water_table_depth (float): depth of water table (in m)
            thickness (float): thickness of foundation base (in m)

        Returns:
            float: force due to water pressure on top of foundation base (in kN)
        '''
        if (D - thickness) <= water_table_depth:
            return 0
        else:
            return (D - thickness - water_table_depth) * B * L * 10
    
    def U2(self, D, B, L, water_table_depth):
        '''
        Calculate force due to water pressure on underside of foundation

        Parameters:    
            D (float): depth of foundation base (in m)
            B (float): width of foundation base (in m)
            L (float): length of foundation base (in m)
            water_table_depth (float): depth of water table (in m)

        Returns:
            float: force due to water pressure on underside of foundation (in kN)
        '''
        if D <= water_table_depth:
            return 0
        else:
            return (D - water_table_depth) * 10 * B * L
        
    def eff_density(self, D, density, water_table_depth):
        '''
        Calculate effective density of soil

        Parameters:    
            D (float): depth of foundation base (in m)
            density (float): bulk density of soil (in kN/m^3)
            water_table_depth (float): depth of water table (in m)

        Returns:
            float: effective density of soil (in kN/m^3)
        '''
        if D <= water_table_depth:
            return density
        else:
            return density - 10
        
    def Nq(self, phi):
        '''
        Calculate bearing resistance factor Nq

        Parameters:    
            phi (float): drained friction angle (in degrees)

        Returns:
            float: bearing resistance factor Nq
        '''
        return math.exp(math.pi * math.tan(math.radians(phi))) * (math.tan(math.radians(45 + phi/2)))**2
    
    def Ng(self, phi):
        '''
        Calculate bearing resistance factor Ng

        Parameters:    
            phi (float): drained friction angle (in degrees)

        Returns:
            float: bearing resistance factor Ng
        '''
        return 2 * (self.Nq(phi) - 1) * math.tan(math.radians(phi))

    def Nc(self, phi):
        '''
        Calculate bearing resistance factor Nc

        Parameters:    
            phi (float): drained friction angle (in degrees)

        Returns:
            float: bearing resistance factor Nc
        '''
        return (self.Nq(phi) -1) / (math.tan(math.radians(phi)))
    
    def Sq(self, phi, B_prime, L_prime):
        '''
        Calculate shape foundation factor Sq

        Parameters:    
            phi (float): drained friction angle (in degrees)
            B_prime (float): effective width (in m)
            L_prime (float): effective length (in m)

        Returns:
            float: shape foundation factor Sq
        '''
        return 1 + (B_prime/L_prime) * math.sin(math.radians(phi))
    
    def Sg(self, B_prime, L_prime):
        '''
        Calculate shape foundation factor Sg

        Parameters:    
            phi (float): drained friction angle (in degrees)
            B_prime (float): effective width (in m)
            L_prime (float): effective length (in m)

        Returns:
            float: shape foundation factor Sg
        '''
        return 1 - 0.3 * (B_prime/L_prime)
    
    def Sc(self, Sq, Nq):
        '''
        Calculate shape foundation factor Sc

        Parameters:    
            Sq (float): shape foundation factor Sq
            Nq (float): bearing resistance factor Nq

        Returns:
            float: shape foundation factor Sc
        '''
        return (Sq * Nq - 1) / (Nq - 1)

    def iq(self, H, V, A_prime, c_prime, phi, B_prime, L_prime, direction=2): # direction 0 means H acts in direction of B', 1 means direction of L'
        '''
        Calculate load inclination factor iq

        Parameters:    
            H (float): horizonal force acting in direction of length (in kN)
            V (float): vertical force (in kN)
            A_prime (float): effective area (in m^2)
            c_prime (float): effective cohesion (in kPa)
            phi (float): drained friction angle (in degrees)
            B_prime (float): effective width (in m)
            L_prime (float): effective length (in m)
            direction (int): direction that H is acting in
            
        Returns:
            float: load inclination factor iq
        '''
        if H == 0:
            return 1
        else:
            assert direction == 1 or direction == 0
            if direction == 0:
                m = (2 + (B_prime/L_prime)) / (1 + (B_prime/L_prime))
                return (1 - H / (V + A_prime * c_prime * (math.tan(math.radians(phi)))**-1)) ** m
            else:
                m = (2 + (L_prime/B_prime)) / (1 + (L_prime/B_prime))
                return (1 - H / (V + A_prime * c_prime * (math.tan(math.radians(phi)))**-1)) ** m
            
    def ig(self, H, V, A_prime, c_prime, phi, B_prime, L_prime, direction=2): # direction 0 means H acts in direction of B', 1 means direction of L'
        '''
        Calculate load inclination factor ig

        Parameters:    
            H (float): horizonal force acting in direction of length (in kN)
            V (float): vertical force (in kN)
            A_prime (float): effective area (in m^2)
            c_prime (float): effective cohesion (in kPa)
            phi (float): drained friction angle (in degrees)
            B_prime (float): effective width (in m)
            L_prime (float): effective length (in m)
            direction (int): direction that H is acting in
            
        Returns:
            float: load inclination factor ig
        '''
        if H == 0:
            return 1
        else:
            assert direction == 1 or direction == 0
            if direction == 0:
                m = (2 + (B_prime/L_prime)) / (1 + (B_prime/L_prime))
                return (1 - H / (V + A_prime * c_prime * (math.tan(math.radians(phi)))**-1)) ** (m+1)
            else:
                m = (2 + (L_prime/B_prime)) / (1 + (L_prime/B_prime))
                return (1 - H / (V + A_prime * c_prime * (math.tan(math.radians(phi)))**-1)) ** (m+1)
    
    def ic(self, iq, Nc, phi):
        '''
        Calculate load inclination factor ic

        Parameters:    
            iq : load inclination factor iq
            Nc : bearing resistance factor Nc
            phi (float): drained friction angle (in degrees)
            
        Returns:
            float: load inclination factor ic
        '''
        return iq - (1 - iq) / (Nc * math.tan(math.radians(phi)))

    def bq(self, alpha, phi):
        '''
        Calculate base inclination factor bq

        Parameters:    
            alpha (float): angle of inclination (in degrees)
            phi (float): drained friction angle (in degrees)
            
        Returns:
            float: base inclination factor bq
        '''
        return (1 - alpha * math.tan(math.radians(phi))) ** 2
    
    def bg(self, alpha, phi):
        '''
        Calculate base inclination factor bg

        Parameters:    
            alpha (float): angle of inclination (in degrees)
            phi (float): drained friction angle (in degrees)
            
        Returns:
            float: base inclination factor bg
        '''
        return (1 - alpha * math.tan(math.radians(phi))) ** 2
    
    def bc(self, bq, Nc, phi):
        '''
        Calculate base inclination factor bc

        Parameters:    
            bq : base inclination factor bq
            Nc (float): bearing resistance factor Nc
            phi (float): drained friction angle (in degrees)
            
        Returns:
            float: base inclination factor bc
        '''
        return bq - (1 - bq) / (Nc * math.tan(math.radians(phi)))
    
    def dq(self, D, B_prime):
        '''
        Calculate embedment depth factor dq

        Parameters:    
        D (float): depth of foundation base (in m)
        B_prime (float): effective width (in m)
            
        Returns:
            float: embedment depth factor dq
        '''
        res = 1 + 0.35 * (D / B_prime)
        return res if res <= 1.7 else 1.7
    
    def dc(self, D, B_prime):
        '''
        Calculate embedment depth factor dc

        Parameters:    
        D (float): depth of foundation base (in m)
        B_prime (float): effective width (in m)
            
        Returns:
            float: embedment depth factor dc
        '''
        res = 1 + 0.35 * (D / B_prime)
        return res if res <= 1.7 else 1.7
    
    def dg(self):
        '''
        Return embedment depth factor dq

        Parameters:    
        -nil-

        Returns:
            float: embedment depth factor dg
        '''
        return 1.0
    
    def volume(self, B, L):
        '''
        Calculate volume of foundation base

        Parameters:
            B (float): width of foundation base (in m)
            L (float): length of foundation base (in m)

        Returns:
            float: volume of foundation base
        '''
        return B * L * 0.5
    
    def price(self, D, volume):
        '''
        Calculate cost of concrete for foundation base

        Parameters:
            D (float): depth of foundation base (in m)
            volume (float): volume of foundation base (in m^3)

        Returns:
            float: volume of foundation base
        '''
        base_price = volume * 100000
        if D <= 1:
            return base_price
        elif D <= 2:
            return base_price * 1.05
        elif D <= 3:
            return base_price * 1.1
        else:
            return base_price * 1.15


def iterate(D_min, D_max, D_step, B_min, B_max, B_step, L_min, L_max, L_step):
    results_undrained = [ ["B", "L", "D", "C1 Load", "R1 Resistance", "C2 Load", "R2 Resistance", "Volume", "Price", ] ]
    results_drained = [ ["B", "L", "D", "C1 Load", "R1 Resistance", "C2 Load", "R2 Resistance", "Volume", "Price", ] ]
    vertical_variable_load = 120
    vertical_permanent_load = 250
    horizontal_variable_load = 0
    horizontal_permanent_load = 0
    unit_weight_concrete = 2400*10/1000
    bulk_density_soil = 20
    thickness_foundation = 0.5
    alpha = 0
    c_prime = 0
    water_table_depth = 2 
    phi = 35 # in degrees
    
    ## UNDRAINED CALCULATIONS
    UC = UndrainedCalc()
    d = D_min
    while d <= D_max:
        b = B_max
        while b >= B_min:
            l = L_max
            while l >= L_min and l >= b:

                w_foundation = b * l * thickness_foundation * unit_weight_concrete # Weight of foundation
                w_backfill = (d - thickness_foundation) * b * l * bulk_density_soil # Weight of backfill
                q = UC.q(d) # Shared q for both combi 1 and 2

                spt_N = UC.spt_N_value(d)
                cu = spt_N * 5

                # COMBI 1
                c1_v1 = 1.35 * (vertical_permanent_load + w_foundation + w_backfill) + 1.5 * vertical_variable_load
                c1_h1 = 1.35 * horizontal_permanent_load + 1.5 * horizontal_variable_load
                c1_e = c1_h1 * thickness_foundation / c1_v1
                c1_B_prime = b - 2*c1_e
                c1_L_prime = l
                c1_A_prime = c1_B_prime * c1_L_prime

                c1_cu = cu
                c1_sc = UC.sc(c1_B_prime, c1_L_prime)
                c1_ic = UC.ic(c1_h1, c1_A_prime, c1_cu)
                c1_bc = UC.bc(alpha)
                c1_dc = UC.dc(d, c1_B_prime)

                c1_r1 = ((2+math.pi) * c1_cu * c1_sc * c1_ic * c1_bc * c1_dc + q) * c1_A_prime
                c1_satisfied = (c1_r1 >= c1_v1)

                # COMBI 2
                c2_v2 = 1.0 * (vertical_permanent_load + w_foundation + w_backfill) + 1.3 * vertical_variable_load
                c2_h2 = 1.35 * horizontal_permanent_load + 1.5 * horizontal_variable_load
                c2_e = c2_h2 * thickness_foundation / c2_v2
                c2_B_prime = b - 2*c2_e
                c2_L_prime = l
                c2_A_prime = c2_B_prime * c2_L_prime

                c2_cu = cu / 1.4
                c2_sc = UC.sc(c2_B_prime, c2_L_prime)
                c2_ic = UC.ic(c2_h2, c2_A_prime, c2_cu)
                c2_bc = UC.bc(alpha)
                c2_dc = UC.dc(d, c2_B_prime)

                c2_r2 = ((2+math.pi) * c2_cu * c2_sc * c2_ic * c2_bc * c2_dc + q) * c2_A_prime
                c2_satisfied = (c2_r2 >= c2_v2)

                if (c1_satisfied and c2_satisfied):
                    volume = UC.volume(b, l)
                    price = UC.price(d, volume)
                    
                    results_undrained.append([b, l, d, c1_v1, c1_r1, c2_v2, c2_r2, volume, price])
                    l -= L_step
                else:
                    break
            b -= B_step
        d += D_step

    ## DRAINED CALCULATIONS
    DC = DrainedCalc()
    d = D_min
    while d <= D_max:
        b = B_max
        while b >= B_min:
            l = L_max
            while l >= L_min and l >= b:
                
                ## DRAINED CALCULATIONS
                ws = vertical_permanent_load
                wb_prime = DC.Wb_prime(d, bulk_density_soil, water_table_depth, thickness_foundation, b, l)
                wf = DC.Wf(b, l, thickness_foundation, unit_weight_concrete)
                u1 = DC.U1(d, b, l, water_table_depth, thickness_foundation)
                u2 = DC.U2(d, b, l, water_table_depth)
                H = 0
                q_prime = DC.q_prime(d, bulk_density_soil, water_table_depth)
                eff_density = DC.eff_density(d, bulk_density_soil, water_table_depth)

                # COMBI 1
                v1 = 1.35 * (ws + wb_prime + wf + u1 - u2) + 1.5 * vertical_variable_load # Fd 1
                c1_phi = phi
                c1_eb = DC.eb(H, v1)
                c1_B_prime = DC.B_prime(b, c1_eb)
                c1_el = DC.el(H, v1)
                c1_L_prime = l
                c1_A_prime = DC.A_prime(c1_B_prime, c1_L_prime)
                c1_c_prime = c_prime

                c1_nq = DC.Nq(c1_phi)
                c1_ng = DC.Ng(c1_phi)
                c1_nc = DC.Nc(c1_phi)

                c1_sq = DC.Sq(c1_phi, c1_B_prime, c1_L_prime)
                c1_sg = DC.Sg(c1_B_prime, c1_L_prime)
                c1_sc = DC.Sc(c1_sq, c1_nq)

                c1_iq = DC.iq(H, v1, c1_A_prime, c1_c_prime, c1_phi, c1_B_prime, c1_L_prime)
                c1_ig = DC.ig(H, v1, c1_A_prime, c1_c_prime, c1_phi, c1_B_prime, c1_L_prime)
                c1_ic = DC.ic(c1_iq, c1_nc, c1_phi)

                c1_bq = DC.bq(alpha, c1_phi)
                c1_bg = DC.bg(alpha, c1_phi)
                c1_bc = DC.bc(c1_bq, c1_nc, c1_phi)

                c1_dq = DC.dq(d, c1_B_prime)
                c1_dc = DC.dc(d, c1_B_prime)
                c1_dg = DC.dg()

                c1_r = c1_A_prime * (c1_c_prime * c1_nc * c1_sc * c1_ic * c1_bc * c1_dc 
                                     + q_prime * c1_nq * c1_sq * c1_iq * c1_bq * c1_dq
                                     + 0.5 * eff_density * c1_B_prime * c1_ng * c1_sg * c1_ig * c1_bg * c1_dg)
                c1_satisfied = (c1_r >= v1)

                # COMBI 2
                v2 = 1.0 * (ws + wb_prime + wf + u1 - u2) + 1.3 * vertical_variable_load # Fd 2
                c2_phi = math.degrees(math.atan( math.tan(math.radians(phi))/1.25 ))
                c2_eb = DC.eb(H, v2)
                c2_B_prime = DC.B_prime(b, c2_eb)
                c2_el = DC.el(H, v2)
                c2_L_prime = l
                c2_A_prime = DC.A_prime(c2_B_prime, c2_L_prime)
                c2_c_prime = c_prime / 1.25

                c2_nq = DC.Nq(c2_phi)
                c2_ng = DC.Ng(c2_phi)
                c2_nc = DC.Nc(c2_phi)

                c2_sq = DC.Sq(c2_phi, c2_B_prime, c2_L_prime)
                c2_sg = DC.Sg(c2_B_prime, c2_L_prime)
                c2_sc = DC.Sc(c2_sq, c2_nq)

                c2_iq = DC.iq(H, v2, c2_A_prime, c2_c_prime, c2_phi, c2_B_prime, c2_L_prime)
                c2_ig = DC.ig(H, v2, c2_A_prime, c2_c_prime, c2_phi, c2_B_prime, c2_L_prime)
                c2_ic = DC.ic(c2_iq, c2_nc, c2_phi)

                c2_bq = DC.bq(alpha, c2_phi)
                c2_bg = DC.bg(alpha, c2_phi)
                c2_bc = DC.bc(c2_bq, c2_nc, c2_phi)

                c2_dq = DC.dq(d, c2_B_prime)
                c2_dc = DC.dc(d, c2_B_prime)
                c2_dg = DC.dg()

                c2_r = c2_A_prime * (c2_c_prime * c2_nc * c2_sc * c2_ic * c2_bc * c2_dc 
                                     + q_prime * c2_nq * c2_sq * c2_iq * c2_bq * c2_dq
                                     + 0.5 * eff_density * c2_B_prime * c2_ng * c2_sg * c2_ig * c2_bg * c2_dg)
                c2_satisfied = (c2_r >= v2)

                if (c1_satisfied and c2_satisfied):
                # if True:
                    volume = DC.volume(b, l)
                    price = DC.price(d, volume)
                    # print([c2_nq, c2_sq, c2_iq, c2_bq, c2_dq])
                    results_drained.append([b, l, d, v1, c1_r, c2_v2, c2_r, volume, price])
                    l -= L_step
                else:
                    break
            b -= B_step
        d += D_step

    # with open("pad_foundation_undrained.csv","w",newline="") as f:
    with open("pad_foundation_undrained_SPECIFIC.csv","w",newline="") as f:
        writer=csv.writer(f)
        writer.writerows(results_undrained)

    # with open("pad_foundation_drained.csv","w",newline="") as f:
    with open("pad_foundation_drained_SPECIFIC.csv","w",newline="") as f:
        writer=csv.writer(f)
        writer.writerows(results_drained)

    print("Done")
                
'''
First iteration
'''
# iterate(D_min=0.5, D_max=4.0, D_step=0.05, B_min=0.1, B_max=4, B_step=0.01, L_min=0.1, L_max=4, L_step=0.05)

'''
Second iteration
'''
# iterate(D_min=3.8, D_max=4.0, D_step=0.05, B_min=0.1, B_max=1, B_step=0.01, L_min=0.1, L_max=1, L_step=0.01)

