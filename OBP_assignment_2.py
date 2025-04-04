import math
import matplotlib.pyplot as plt
import numpy as np
import time
import matplotlib.ticker as ticker
import streamlit as st

def plot_system(system):
    """
    Plots the PMF and CCDF of the system side-by-side.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))  # 1 row, 2 columns

    # PMF (Bar Chart)
    plot_system_pdf(system, axes[0])

    # CCDF
    plot_system_ccdf(system, axes[1])

    fig.suptitle(f"Reliability of a K = {K} out of N = {system.N} system, \n with {system.S} repairman and {system.type}", fontsize=12)  # Overall figure title
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]) # Adjust layout to prevent overlap
    plt.show()

def plot_system_pdf(system, ax):
    
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.bar(list(range(K)), system.stat_dist[:K], color = 'tab:blue', label = 'System is down')
    ax.bar(list(range(K, system.N + 1)), system.stat_dist[K:], color='tab:green', label='System is up')
    
    ax.set_title("Probability Mass Function (PMF)")
    ax.set_xlabel("Number of Working Components")
    ax.set_ylabel("P(X=x)")

    ax.legend()
    
def plot_system_ccdf(system, ax): 
    
    ax.grid(True, linestyle='--', alpha=0.7)

    # Calculate CCDF
    system.ccdf = get_ccdf_inclusive(system)

    ax.plot(list(range(system.N + 1)), system.ccdf, marker='o', linestyle='-')
    ax.fill_between(list(range(K, system.N + 1)), system.ccdf[K:system.N+1], color='tab:green', label='System is up')
    # ax.set_xlim(system.N+1-0.5, - 0.5)

    # Add annotation line.
    frac_up = sum(system.stat_dist[K:])
    ax.hlines(frac_up, 0, system.N, color='red', linestyle='dashed', label=f'P(System Up) = {frac_up:.3f}')

    ax.set_xlabel("Number of Working Components")
    ax.set_ylabel("P(X $\geq$ x)")
    ax.set_title("Reliability Function: Probability of X or More Working Components")
    ax.legend()

def get_fraction_up(system):
    system.stat_dist = system.get_stationary_distribution()
    frac_up = 0

    for i in range(K, system.N + 1):
        
        frac_up += system.stat_dist[i]

    return frac_up

def get_ccdf_inclusive(system):
    """
    Calculates the Complementary Cumulative Distribution Function (CCDF)
    representing P(X >= x).
    """
    stat_dist_ccdf = []
    tot = 0
    stat_dist = system.stat_dist[::-1] #reverse the array.
    for element in stat_dist:
        tot += element
        stat_dist_ccdf.append(tot)

    return stat_dist_ccdf[::-1] #reverse the array again to return it to the original order.

class WarmStandby():
    def __init__(self, s, n):
        self.type = "Warm standby" 
        self.S = s
        self.N = n       
    
    def get_pi0(self):
        """Compute π0 using instance variables, that all machines are down."""

        pi_0 = 0
        for i in range(self.N - self.S + 1):
            pi_0 += (self.S*RHO)**i / math.factorial(i)

        for i in range(self.N - self.S + 1, self.N + 1):
            num = (self.S**(self.N - self.S)) * math.factorial(self.S)
            den = (math.factorial(self.N - i)) * (math.factorial(i))

            pi_0 += num/den * RHO**i

        return 1/pi_0 
    
    def get_stationary_distribution(self):
        '''Calculates the stationary distribution'''
        stat_dist = []
        self.pi_0 = self.get_pi0()
        pi_tot = 0
                
        for i in range(self.N - self.S + 1):
            pi_i = (self.S*RHO)**i / math.factorial(i) * self.pi_0
            stat_dist.append(pi_i)

        for i in range(self.N - self.S + 1, self.N + 1):
            num = (self.S**(self.N - self.S)) * math.factorial(self.S)
            den = (math.factorial(self.N - i)) * (math.factorial(i))

            pi_i = num/den * RHO**i * self.pi_0
            stat_dist.append(pi_i)

        if abs(sum(stat_dist) - 1) < 0.0001:
            return stat_dist
        else:
            print("Total law of probability not satisfied") 

class ColdStandby():
    def __init__(self, s, n):
        self.type = "Cold standby"    
        self.S = s
        self.N = n     

    def get_pi_down(self):
        pi_down = 0

        if K > self.N + 1 - self.S:
            for i in range(K-1, self.N + 1):
                num = math.factorial(K-1) * math.factorial(self.N - K + 1)
                den = math.factorial(i) * math.factorial(self.N - i)
                pi_down += (num/den) * RHO**(i - K + 1)
        
            return 1/pi_down
        
        for i in range(K-1, self.N - self.S + 2):
            num = math.factorial(K-1)
            den = math.factorial(i)
            pi_down += (num/den) * (self.S * RHO)**(i-K+1)

        for i in range( self.N - self.S + 2, self.N + 1):
            num = math.factorial(K-1) * math.factorial(self.S)
            den = math.factorial(i) * math.factorial(self.N - 1)
            
            s_pow = self.S**(self.N - self.S - K + 1)
            rho_pow = (RHO)**(i-K+1)

            pi_down += (num/den) * s_pow * rho_pow

        return 1/pi_down

    def get_stationary_distribution(self):
        pi_down = self.get_pi_down()
        
        stat_dist = [0]*(K-1)
        pi_tot = 0    

        if K > self.N + 1  - self.S:
            for i in range(K-1, self.N + 1):
                num = math.factorial(K-1) * math.factorial(self.N - K + 1)
                den = math.factorial(i) * math.factorial(self.N - i)
                pi_i = (num/den) * RHO**(i - K + 1) * pi_down
                
                stat_dist.append(pi_i)
                
            
        else:
            for i in range(K-1, self.N - self.S + 2):
                num = math.factorial(K-1)
                den = math.factorial(i)
                pi_i = (num/den) * (self.S * RHO)**(i-K+1) * pi_down
                
                stat_dist.append(pi_i)
                

            for i in range( self.N - self.S + 2, self.N + 1):
                num = math.factorial(K-1) * math.factorial(self.S)
                den = math.factorial(i) * math.factorial(self.N - 1)
                
                s_pow = self.S**(self.N - self.S - K + 1)
                rho_pow = (RHO)**(i-K+1)

                pi_i = (num/den) * s_pow * rho_pow * pi_down
                stat_dist.append(pi_i)
                

        if abs(sum(stat_dist) - 1) < 0.001:            
              return stat_dist
        else:
            print(sum(stat_dist))
            print(self.S, self.N)
            print("Total law of probability not satisfied")
            

class SystemOptimizer:
    def __init__(self, type, cost_repairman, cost_component, cost_downtime):
        self.type = type
        self.cost_repairman = cost_repairman
        self.cost_component = cost_component
        self.cost_downtime = cost_downtime

    def get_downtime(self, s, n):
        if self.type == "Warm":
            system = WarmStandby(s, n)
        if self.type == "Cold":
            system = ColdStandby(s, n)
        return 1-get_fraction_up(system)

    def calculate_costs(self, s, n):
        if s > n:
            return float('inf')
        
        downtime = self.get_downtime(s, n)

        return self.cost_repairman*s + self.cost_component*n + downtime*self.cost_downtime

    def optimize(self, max_downtime = None, exploration_constant = None):
        """
        Optimizes the number of repairmen (s) and components (n) to minimize total cost,
        optionally subject to a maximum downtime constraint.

        Args:
            max_downtime (float, optional): The maximum allowable downtime. If provided,
                combinations of s and n that exceed this downtime are discarded.
                If None (default), no downtime constraint is applied.
            exploration_constant (int, optional): Controls the number of iterations
                the algorithm explores without improvement. Higher values increase exploration but also
                computation time. If None (default), it defaults to 2 * K, where
                K is a class attribute.

        Returns:
            self (object): The optimizer object with optimal s and n values, and the
                        generated cost surface.
            """

        self.n_range = None
        self.optimal_n = None
        self.optimal_s = None

        if exploration_constant == None:
            exploration_constant = 2*K

        self.min_cost = float('inf')
        n = K-1
        
        counter = 0
        cost_surface_rows = []

        current_row_min = None
        previous_row_min = float('inf')

        ## keep adding rows, until the minimal rowcost of 'exploration_constant' consecutive new rows is higher than the previous
        while counter < exploration_constant:        
            #dynamically increase the search space to find the optimal   
            n += 1
            cost_row = []
            
            # s cannot be larger than n, more repairman than components does not make sense
            for s in range(1, n):
                # if maxdowntime constraint is added, downtime higher than treshold punishes the costs
                if max_downtime != None and self.get_downtime(s, n) > max_downtime:
                    cost = float('inf')      

                else:
                    cost = self.calculate_costs(s, n)

                    if cost < self.min_cost:
                        # update class variables
                        self.min_cost = cost
                        self.optimal_n = n
                        self.optimal_s = s
                        counter = 0  

                cost_row.append(cost)

            current_row_min = min(cost_row)
            # current row has lower cost than previous
            if (current_row_min < previous_row_min):           
                counter = 0

            # if a minimum has been found, count the added rows without improvement
            if self.optimal_n != None:
                counter += 1

            previous_row_min = current_row_min
            cost_surface_rows.append(cost_row)

        self.n_range = list(range(K, n + 1))
        self.s_range = list(range(1, n + 1))

        # Determine matrix dimensions
        num_rows = len(cost_surface_rows)  # Number of rows
        num_cols = len(cost_surface_rows[-1])  # Width determined by the last row

        # Create a rectangular matrix filled with zeros
        self.cost_surface = np.full((num_rows, num_cols), np.nan)

        # Fill the matrix with data
        for i, row in enumerate(cost_surface_rows):
            self.cost_surface[i, :len(row)] = row  # Copy existing values

        return self

    def plot_loss_surface_2d(self):
        """Plots the loss surface as a 2D heatmap with improved visualization."""
        fig, ax = plt.subplots(figsize=(8, 6))

        # Calculate adjusted extent
        s_min = min(self.s_range) - 0.5
        s_max = max(self.s_range) + 0.5
        n_min = min(self.n_range) - 0.5
        n_max = max(self.n_range) + 0.5

        # Use a better colormap (e.g., 'RdYlGn_r' or 'Greens_r')
        im = ax.imshow(self.cost_surface,
                        extent=[s_min, s_max, n_min, n_max], # changed n_max and n_min
                        aspect='auto', cmap='RdYlGn_r', origin='lower') # added origin='lower'

        # Add colorbar with proper label
        cbar = fig.colorbar(im, ax=ax, label="Total Cost")

        # Set axis labels and title
        ax.set_xlabel(f"Number of Repairmen (S)\nOptimum at (S, N) = ({self.optimal_s}, {self.optimal_n}), Cost = {self.calculate_costs(self.optimal_s, self.optimal_n):.2f}")
        ax.set_ylabel("Number of Components (N)")
        ax.set_title(f"Cost Heatmap in {self.type} standby with K={K}")
        
        # Add a marker for the optimal point
        ax.plot(self.optimal_s + 0.5, self.optimal_n, 'ro', markersize=8, label='Optimum')

        # Add a legend
        ax.legend()
        ax.set_xlim(0.5, max(self.s_range))
        # Improve tick label appearance
        ax.tick_params(axis='both', which='major', labelsize=10)
            # Force integer ticks on both axes
        ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

        plt.tight_layout()
        plt.show()

def main():
    st.title("OBP assignment 2")
    st.subheader("a) Reliability of K-out-of-N system:")

    global K
    global FAILURE_RATE
    global REPAIR_RATE
    global RHO

    # Default constants
    FAILURE_RATE = st.number_input("Failure rate", min_value=0.0, value=3.0)
    REPAIR_RATE = st.number_input("Repair rate", min_value=0.0, value=2.0)
    RHO = REPAIR_RATE / FAILURE_RATE

    # Input widgets
    K = st.number_input("Number of compontens for the system to function", min_value=1, value=3)
    s = st.number_input("Number of repairmen", min_value=1, value=3)
    n = st.number_input("Number of components", min_value=K, value=7)

    type = st.radio("Choose System Type", options=["Cold", "Warm"])

    # Create the appropriate system
    if type == "Warm":
        system = WarmStandby(s, n)
    else:
        system = ColdStandby(s, n)

    # Calculate and display fraction up
    system.frac_up = get_fraction_up(system)
    st.subheader(f"Probability system is up: {system.frac_up:.3f}")
    
    plot_system(system)
    st.pyplot(plt)
    
    #part b
    st.subheader("b) Optimization of the number of repairmen and components:")

    cost_repairman = st.number_input("Cost per repairmen", min_value=0.0, value=10.0)
    cost_component = st.number_input("Cost per component", min_value=0.0, value=50.0)
    cost_downtime = st.number_input("Downtime cost", min_value=0.0, value=1000.0) 

    optimizer = SystemOptimizer(type, cost_repairman, cost_component, cost_downtime)
    optimizer.optimize(max_downtime=None, exploration_constant=None)

    st.write(f"Optimal number of components: **{optimizer.optimal_n}**")
    st.write(f"Optimal number of repairmen: **{optimizer.optimal_s}**")
    st.write(f"Minimum total cost: **{optimizer.min_cost:.2f}**")

    if type == "Warm":
        opt_system = WarmStandby(optimizer.optimal_s, optimizer.optimal_n)
    else:
        opt_system = ColdStandby(optimizer.optimal_s, optimizer.optimal_n)

    opt_system.frac_down = 1 - get_fraction_up(opt_system)
    st.write(f" Fraction system is down: **{optimizer.get_downtime(optimizer.optimal_s, optimizer.optimal_n):.4f}**")

    # Show plots
    optimizer.plot_loss_surface_2d()
    st.pyplot(plt)

    plot_system(opt_system)
    st.pyplot(plt)
    
    

if __name__ == "__main__":
    main()    