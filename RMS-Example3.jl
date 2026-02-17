### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 79bde88e-0788-4d67-975c-5d47b72846ff
using PlutoUI, Plots, PlutoUI, LaTeXStrings, Printf, PlutoTeachingTools

# ╔═╡ b8fd235e-072b-4697-97a3-b537aba895f2
md"""# Model of the vehicle
The auxiliary code is hidden by default. To see the code, put the cursor on top of the cell and press the eye to show or hide the code
"""


# ╔═╡ 5516b798-bf15-4cc3-a190-dca31aea2304
PlutoUI.TableOfContents()

# ╔═╡ b68dcc2e-d48f-4286-9baf-b7d28a17bd2b
function medusa_model(yaw,u,v,r,a_s,a_p)
    # Dynamics do not depend on x and y, so the only input state is yaw
    
	# Thrusters
    K=3.6e-3; 
    a_max=100.0; 

    a_s=min(a_max,max(-a_max,a_s)); # Saturation
    a_p=min(a_max,max(-a_max,a_p));

    Fs=K * abs(a_s)* a_s;
    Fp=K * abs(a_p)* a_p;

    # Dynamics (or Kinetics)

    # Mass, inertia and added masses
    m_u = 37.0; # kg
    m_v = 47.0; # kg
    m_r = 4.64; # kg*m*m
    m_uv = m_u - m_v;

    X_u = -0.2; # kg/s Linear damping
    Y_v = -55.1; # kg/s
    N_r = -4.14; # kg*m*m/s

    X_uu = -25.0; # kg/m Quadratic damping 
    Y_vv = -101.3; # kg/m 
    N_rr = -6.23; # kg*m*m 

	# Drag
    D_u = -X_u - X_uu*abs(u);
    D_v = -Y_v - Y_vv*abs(v);
    D_r = -N_r - N_rr*abs(r);
	
	# Force
    d=0.25;# m Distance from the actuator to the symmetry axis
    tau_u=Fs+Fp; # Total force N
    tau_r=d*(Fp-Fs); # Torque N*m
	
    # Final Dynamics
    du = (tau_u + m_v*v*r - D_u*u)/m_u;
    dv = (-m_u*u*r - D_v*v)/m_v;
    dr = (tau_r + m_uv*u*v - D_r*r)/m_r;

    # Kinematics
    dx = u*cos(yaw) - v*sin(yaw);
    dy = u*sin(yaw) + v*cos(yaw);
    dyaw = r;

    return (dx,dy,dyaw,du,dv,dr)
end

# ╔═╡ d6d4ee7b-7a67-4ae6-9d72-4c21837b8182
function paint_vehicle!(p,x,y,yaw)
	formx=[  -1,    0, 1, 0, -1, -1]
	formy=[-0.5, -0.5, 0, +0.5, +0.5, -0.5 ]
	
	X=x .+ formx*cos(yaw)-formy*sin(yaw)
	Y=y .+ formx*sin(yaw)+formy*cos(yaw)
	vehiculo=Shape(Y,X) # NED
	p=plot!(p,vehiculo,label="Vehicle")
end

# ╔═╡ ae3961d4-be42-4268-b60f-653fe5842adc
set_point=[0.0, 0.0] # x_r, y_r

# ╔═╡ 91139377-20b7-428c-ae62-8be53d5eb685
Tf=1600.0

# ╔═╡ 25e7d4b0-51aa-4064-ad7b-07593d75c628
dt=0.1;

# ╔═╡ 68042e13-9df3-484c-afb5-dd7cc549eb0d
## Energy consumption:

# ╔═╡ 4da9ca61-abc6-4b59-93d6-8bb8c83d8aa6
md"""# Enhancements of the control law
In the last example, we arrived at a pretty good controller that has nice properties:

- It is simple (good for implementation and debugging)

- It is continuous (good for the actuators)

- It *seems* to work from any initial conditions


The problems that we face are: 

- We don't have formal warranties that the vehicle will reach and remain in a recovery area from a set of initial conditions (it usually does, but we are not sure by design)

- As the actuators are active all the time, the control law could drain the batteries fast (bad for a rescue scenario)"""




# ╔═╡ baec733f-9474-465f-90a5-401c76e0c0fd
md""" 
## Initial and final set

In order to have formal warranties of stability, we will define two sets:

$\Omega =\{\chi| x^2+y^2 <R^2, u=v=r=0\}$

$\Omega_r =\{\chi| (x-x_r)^2+(y-y_r)^2 < R_f^2\}$

The set $\Omega$ defines the possible initial conditions considered in the optimisation problem (a circle of radius $R$), and $\Omega_r$ a neighbourhood of the recovery point that we will **enforce** the trajectory to end (a circle of radius $R_f$).

In this example, we will use $R=100m$ and $R_f=3m$

Notice that, for the optimisation problem to be feasible, $R_f$ must be chosen appropiately, in our case we choose $R_f=3m$ because the reference controller (with parameters $a_1=30$, $a_2=60$ and $\Delta \psi =-63\pi/180$) is able to drive the vehicle to this set.

"""


# ╔═╡ 87cc3dca-6b63-4e90-9b90-8dbf3c9574ed
begin
	R=100
	Rf=3
end

# ╔═╡ ed5233f6-ff79-4086-9c51-b5e90eff7648
md"""## Energy consumption

To compute the energy consumption, we integrate the instantaneous power: 

$E=\int_{0}^{T} P(t) \; \mathrm{d}t \approx \Delta t\sum_{i=0}^{T/\Delta t} P(t_i)=J_E$

In order to compute the power, we need a model of the propeller that returns the torque that is needed to apply in any given condition. This model also needs to compute the hydrodynamic drag coefficient for all possible angles of attack $\beta$:  



"""

# ╔═╡ 31b5d97b-96c7-48a0-87a3-450a8170ef7e
function four_quadrant_model(β)
# "Four-Quadrant Model" adapted from de A. Häusler PhD Thesis.
	
	# Model Parameters
	OL=-1.6157*pi/180;
	OD=1.9309*pi/180;
	CDmin=0.0273;
	CDmax=1.0383;
	CLmax=0.5749;
	ζ=atan(1.35/pi); # Pitch angle of the propeller 
	
	# Angle of attack at the propeller blade
	γ=ζ-β;
	
	# Drag and lift coefficients CL y CD
	CLL=CLmax*sin(2*(γ-OL));
	CDL=(CDmax-CDmin)*(1-cos(2*(γ-OD)))/2 + CDmin;
	
	CT=CLL*cos(β)-CDL*sin(β);
	CQ=(-0.7/2)*(-CLL*sin(β)-CDL*cos(β));

	return (CT,CQ)
end

# ╔═╡ 49dd9b51-29a3-4851-bd12-8d386a2c6612
function propeller_model(wp,u,r) 
# "Four-Quadrant Model" adapted from de A. Häusler PhD Thesis.
	ly=0.15; # 'y' component of the distance from the centre of the propeller to the centre of mass of the vehicle in the body frame
	dp=0.076; # Propeller diameter
	Rp=dp/2; # Propeller radius
	rho=1023; # Density of water
	
	vap=u-ly.*r; # Advance velocity of the propeller
	vpp=0.7*Rp.*wp; # Lateral velocity 
	β=atan(vap,vpp); # Advance angle of the propeller
	_,CQ=four_quadrant_model(β);
	Qp=rho/2*CQ*(vap^2+vpp^2)*pi*Rp^2*dp;

	return Qp
end

# ╔═╡ ab53825c-7bdf-40ed-8126-b24780c84140
md"Now that we have the torque, we can compute the power consumption using an electric model of the motor:"

# ╔═╡ fa590ceb-e8b4-4e0f-956a-1b0ec9b62d0f
function instant_power_compsumption(ap,u,r)
# Adapted from de A. Häusler PhD Thesis. 
	nmax=71.6; # Max rps of the actuator
	np=nmax*ap/100; # Turning rate rps
	wp=2*pi*np; # Angular velocity rad/s.
 
	Qp=propeller_model(wp,u,r); # Torque
	Kt=0.2599022
	b=1e-6
	Ip=(b*wp+Qp)/Kt # Actuator current
	Ra=0.66
	Ke=Kt
	Pp=26; # Internal idle power W

	# Finally
	V=Ra*Ip +Ke*np;
	P_motor=V*Ip;
	P_motor=max(P_motor,0*P_motor); # Do not regenerate power by breaking
	P=P_motor+Pp;
	return P
end


# ╔═╡ cee61b43-23d9-420d-afee-c29c24d1dd08
md"And finally, we are ready to compute the energy consumption"

# ╔═╡ c8acb0b9-5c22-4bf5-b181-c943313f9da1
function J_E(dt,ap,u,r)
	E=0.0
	for i in 1:length(ap)
		E+= instant_power_compsumption(ap[i],u[i],r[i])
	end
	return dt*E
end

# ╔═╡ 80fa9ca7-a435-4f5e-9d36-d0a3deb2a3b4
md"""## Final cost function
Now we combine both the ISE and the energy consumption to get:

$J=(1-\lambda) \frac{J_{ISE}}{J_{ISE_{ref}}} +  \lambda \frac{ J_{E} }{J_{E_{ref}} }$

Where $J_{ISE_{ref}}$ and $J_{E_{ref}}$ are the ISE and energy costs of the reference controller starting from initial condition $\boldsymbol{\chi_{0_w}}=[70.26, 70.26, 5.14, 0, 0, 0]^T$ and $0<\lambda<1$ is a weighting factor that tunes the importance of ISE and energy in the cost function.

Notice that the normalisation has two purposes:

- It makes the cost dimensionless (we cannot mix Jules with $m^2s$)
- It makes the cost easy to interpret: $J=1$ means equal to the reference controller in the worst case for any initial condition in $\Omega$ (which is attained in $\boldsymbol{\chi_{0_w}}$). Thus, the optimal controller must have $J<1$ for all initial conditions in $\Omega$.

"""

# ╔═╡ 71d22223-7bec-4409-862b-eca7854368b8
x0w=[70.26, 70.26, 5.14, 0, 0, 0]

# ╔═╡ e160754f-6894-462f-8550-209f818b485f
function position_control(x,y,ψ,u,v,r,t,θ,sp)
	a1,a2,Δψ=θ # Get the values of the parameters
	a0=(a1+a2)/2
	Δa=(a2-a1)/2
	xr,yr=sp
    a_s=0.0 # BROKEN ACTUATOR
	ψ_r=atan(yr-y,xr-x)
	a_p=a0+Δa*sign(sin(ψ-ψ_r+Δψ))
    return (a_s,a_p)
end

# ╔═╡ f819631b-81b4-46c7-ac8c-0ec210a0aed1
function J_ISE(dt,x,y,sp)
	xr,yr=sp
	return dt*sum( (x.-xr).^2 + (y.-yr).^2 )
end

# ╔═╡ 826badb3-f4b0-4049-87fb-273498117dce
md"""# Minimax optimisation problem
We want the vehicle to:
- Optimise the cost $J$
- Converge to the set $\Omega_r$ from **all** initial condition $\chi_0 \in \Omega$

Thus, we perform a minimax optimisation


$\boldsymbol{\theta}^* = \underset{\boldsymbol{\theta}}{\operatorname{arg\,min}}\ \underset{\boldsymbol{\chi_0} \in \Omega}{max}\ J(\boldsymbol{\theta},\boldsymbol{\chi_0})$

Subject to the dynamics of the system and:

$\boldsymbol{\chi}(t) \in \Omega_r \ \ for \ \ t_f<t<T_f$

Where $t_f$ is the transitory time and $T_f$ is the final recovery time.

In this this example, we will consider $t_f=1280s$ and $T_f=1600s$.

"""

# ╔═╡ 5e225716-e396-42ed-a415-99c45a522bfa
tf=1280

# ╔═╡ c8717d42-4423-4d70-b8ea-8c7f48a2184b
function restriction(x,y)
	n=Int(tf/dt)
	d=sqrt.(x[n:end].^2+y[n:end].^2)
    if maximum(d)<Rf
		return true
	else 
		return false
	end
end


# ╔═╡ bbc464ef-852f-4dc3-ae75-f48516327bb2
md"""
Again, it is important to choose values that are feasible, in this case these values are feasible but the reference controller *is not able to satisfy the restriction*, as we can see in the figure, at $t=t_f$ the vehicle is at more than $R_f$ distance of the target point. Thus, the optimal controller must improve the reference one.

"""

# ╔═╡ 8829a36c-266e-40a7-9098-dc7ad3c65d64
md"As the figure shows, the trajectory does not enter $\Omega_r$ in time $t_f$:"

# ╔═╡ 5a6ad0a6-9402-4b8c-997f-69ebbb2bf73a
md"""## Time to play
We have optimised the cost for several values of $\lambda$ and the best results are:
"""

# ╔═╡ 17c6b93f-3820-4957-b9a4-055e6092c212
begin

lambdas=collect(0.05:0.05:0.95)


op_coeficients=[100.532867769320	100.532867769320	100.532867769320	101.326505518243	101.326505518243	99.4563792572124	99.4563792572124	98.8535908909393	99.4563792572124	99.4563792572124	98.2916539818171	97.3563442909694	97.3563442909694	97.3563442909694	95.9467159031433	95.9467159031433	95.9467159031433	94.7466959224330	94.7466959224330;
-29.0876917035200	-29.0876917035200	-29.0876917035200	-27.4160760510659	-27.4160760510659	-29.2451262165943	-29.2451262165943	-29.2451262165943	-29.2451262165943	-29.2451262165943	-29.2451262165943	-29.2451262165943	-29.2451262165943	-29.2451262165943	-29.2451262165943	-29.2451262165943	-29.2451262165943	-29.2451262165943	-29.2451262165943;
12.6961906899468	14.0528338936849	14.3191592709080	14.6600521598577	14.6600521598577	14.9484490548879	15.7774913291685	15.2951433520933	17.9574121613783	17.9880129522896	18.3293373548552	17.9880129522896	17.9880129522896	17.9880129522896	17.9880129522896	17.9880129522896	17.9880129522896	17.9880129522896	17.9880129522896;
-8.38275153739373	-8.54520323905089	-8.63465695427911	-12.5747913567929	-12.5747913567929	-9.13155266786018	-9.31694659382359	-9.05550056832865	-9.38666094190833	-9.77631806123682	-10.2096930810532	-10.2096930810532	-10.2096930810532	-10.2096930810532	-10.2096930810532	-11.9520546019131	-11.9520546019131	-13.8933147490217	-13.8933147490217;
4.05404613900951	4.05404613900951	4.05404613900951	7.63544455550241	7.63544455550241	5.04324498151246	4.35884481732756	5.04324498151246	4.47962131688959	5.04324498151246	3.44845724494511	4.22800360946867	4.22800360946867	4.22800360946867	7.19975783588264	7.19975783588264	7.19975783588264	7.19975783588264	7.19975783588264;
-3.67869927767402	-3.67869927767402	-3.67869927767402	0.390370134729956	0.390370134729956	-2.61512354855936	-2.61512354855936	-1.88972522321389	-2.61512354855936	-2.61512354855936	-1.78609689102168	-1.78609689102168	-1.78609689102168	-1.78609689102168	-1.78609689102168	-1.78609689102168	-1.78609689102168	-0.988158733141668	-0.988158733141668;
-5.23953236146430	-6.24589003403479	-6.82375090301668	-12.0427367669335	-12.0427367669335	-9.49236546636493	-9.49236546636493	-10.4033893827602	-11.3656175837566	-11.6793392244289	-11.4422401642715	-12.7394681267236	-12.7463563101752	-14.2468534005232	-14.8355106528975	-15.8928105108896	-15.8350556116205	-16.5267490903696	-16.5267490903696;
3.48028303203501	3.48028303203501	3.48028303203501	1.67664637793198	-0.324064944549125	3.48028303203501	3.48028303203501	2.21099941787525	2.62974038581252	2.09700902467407	2.31140345718952	2.40159691623454	2.31140345718952	2.40159691623454	2.40159691623454	2.81512541713905	2.53263269635654	2.53263269635654	2.53263269635654;
-6.55885411619219	-7.44943277489752	-6.55885411619219	-4.73038197844138	-5.07189697262884	-6.55885411619219	-6.55885411619219	-6.68098687161938	-6.68098687161938	-6.68098687161938	-6.68098687161938	-6.68098687161938	-6.68098687161938	-5.39908235913398	-5.39908235913398	-5.39908235913398	-5.39908235913398	-4.73038197844138	-4.73038197844138;
2.61993846576766	2.61993846576766	2.61993846576766	3.42831334268964	3.42831334268964	2.10551875029399	2.44000172497672	2.90729850265861	2.90729850265861	3.42831334268964	3.05898636747777	3.42831334268964	3.42831334268964	3.42831334268964	3.42831334268964	3.42831334268964	3.42831334268964	5.39402057739657	5.39402057739657;
3.14463787294054	4.39780108861901	3.85120862220606	4.20142579291427	5.75012757280287	5.68356697871918	6.43013229869785	6.77675191154508	6.77675191154508	7.62996942009213	8.26590182805132	8.87272355942139	8.87272355942139	9.09184672641997	7.34025738474273	8.33607458340484	8.26590182805132	8.87272355942139	8.87272355942139]
	
end

# ╔═╡ 8db37db2-116f-4d71-ba01-495d8a584a98
md"""
Now, change $\lambda$ to see the trade-off between converging fast to the recovery point and reducing the energy consumption:
"""

# ╔═╡ b212918a-0325-4691-b464-22ab12aa30cb
@bind vy_current Slider(0:0.001:0.07, default=0, show_value=true)

# ╔═╡ 25a4d818-4cde-46b2-8af3-34f63647558a
function simulator(X0,Tf,dt,controller,θ,sp)
	N=Int(ceil(Tf/dt))
        
    x=zeros(1,N+1);
    y=zeros(1,N+1); 
    u=zeros(1,N+1); 
    v=zeros(1,N+1); 
    yaw=zeros(1,N+1); 
    r=zeros(1,N+1);
    a_s=zeros(1,N+1);
	a_p=zeros(1,N+1);
    t=0:dt:Tf;
    # Initial values
    
    x[1]=X0[1];
    y[1]=X0[2];
    yaw[1]=X0[3];
    u[1]=X0[4];
    v[1]=X0[5];
    r[1]=X0[6]; 
    
    # Euler method (not very sophisticated)
    for i=1:N 
		# Compute the control
        a_s[i],a_p[i]=controller(x[i],y[i],yaw[i],u[i],v[i],r[i],t[i],θ,sp);
		# Apply to the vehicle
        (dx,dy,dyaw,du,dv,dr)=medusa_model(yaw[i],u[i],v[i],r[i],a_s[i],a_p[i]);     
        # Update the states
        x[i+1]=x[i]+dx*dt;
        y[i+1]=y[i]+dy*dt  + vy_current*dt; # Disturbance
        yaw[i+1]=yaw[i]+dyaw*dt;
        u[i+1]=u[i]+du*dt;
        v[i+1]=v[i]+dv*dt;
        r[i+1]=r[i]+dr*dt;
    end
	a_s[N+1]=a_s[N]
	a_p[N+1]=a_p[N]

    return (x=x,y=y,yaw=yaw,u=u,v=v,r=r,t=t,a_s=a_s,a_p=a_p)
end

# ╔═╡ 126f1b41-6b18-4359-9945-ec18daf32eee
begin
	x_base,y_base,_,u_base,_,r_base,_,_,ap_base=simulator(x0w,Tf,dt,position_control,[30 60 -63*pi/180],set_point)
	Jref_ISE=J_ISE(dt,x_base,y_base,set_point)
	Jref_E=J_E(dt,ap_base,u_base,r_base)
	
	(Jref_ISE=Jref_ISE,Jref_E=Jref_E)
end

# ╔═╡ 2dc7ec04-faea-4f6a-99bb-f226b7075cba
function J(dt,x,y,sp,u,r,ap,λ)
	JE=J_E(dt,ap,u,r)
	JISE=J_ISE(dt,x,y,sp)
	return (1-λ)JISE/Jref_ISE + λ*JE/Jref_E
end

# ╔═╡ 60aaed79-479c-4d4a-bfcf-03ec77b71ad6
restriction(x_base,y_base) # The reference controller is not fast enough

# ╔═╡ 65efcc7e-d127-4a71-9412-3a774b46b507
function select_coef(λ)
	# Select the appropriate values for this lambda
	index=Int(round(λ/0.05))
	return op_coeficients[:,index]
end

# ╔═╡ 64cb1add-744c-42fc-90ca-ef6817bc2fc0
function fourier_control(x,y,ψ,u,v,r,t,θ,sp)
	xr,yr=sp
    a_s=0.0 # BROKEN ACTUATOR
	ψ_r=atan(yr-y,xr-x)
	
	α=ψ - ψ_r
    # Fourier series of the relative angle alpha
    a_cos=0.0;
    a_sin=0.0;
    for j=2:2:(length(θ)-1)
        a_cos+= θ[j]*cos((j/2)*α);
        a_sin+= θ[j+1]*sin((j/2)*α);
    end
    a=θ[1]/2 + a_cos + a_sin;
    
	a_max=60.0; # Saturation
    a_p=min(a_max,max(-a_max,a));
	
    return (a_s,a_p)
end

# ╔═╡ faa44634-49c8-43dc-8600-77b93e6ce15b
begin
	X0=x0w
    #X0=[-20.0,-20.0,0.0,0.0,0.0,0.0,0.0]# You could change if you wish
end

# ╔═╡ b959045e-ce5b-45aa-a9f0-540a18c1c2de
# If you want to save the last figure, uncomment this line:
#savefig("figs.pdf")

# ╔═╡ 1038d5ce-8aeb-4d8f-bbb0-0514fc34e7b1
@bind λ Slider(0.05:0.05:0.95 ; default=0.5, show_value=true)

# ╔═╡ e11b958e-fb9b-4811-9922-089c3e7aacae
function show_figure(now,dt,Tf,control,θ,sp)
	
	x,y,yaw,u,v,r,t,a_s,a_p=simulator(X0,Tf,dt,control,θ,sp)

	
	N=length(t)
	n=Int(ceil((N-1)*now/Tf)+1) # Index of the actual time 

	# Trajectory in 2D plane
	p1=plot(y'[1:(n-1)],x'[1:(n-1)],aspect_ratio=:equal,label=L"Trajectory\ for\ t<t_f", legend=:bottomright)
	plot!(p1,  y'[n:end],x'[n:end],aspect_ratio=:equal, label=L"""Trajectory\ for\ t > t_f""")
	xlabel!(L"y") # NED
	ylabel!(L"x")

    p1=paint_vehicle!(p1,x[n],y[n],yaw[n])
	p1=scatter!(p1,[sp[1]],[sp[2]],label="Setpoint")

    circlex=Rf*sin.(0:0.01:2*pi)
	circley=Rf*cos.(0:0.01:2*pi)

    plot!(p1, circley,circlex,color=:green, label=L"""\Omega_r""")

    # Evolution of the state and control actions
	p2=plot(t,u',xlabel=L"t",ylabel=L"u",label=:none)
	p2=scatter!(p2,[t[n]],[u[n]],label=:none)
	
	p3=plot(t,v',xlabel=L"t",ylabel=L"v",label=:none)
	p3=scatter!(p3,[t[n]],[v[n]],label=:none)
	
	p4=plot(t,r',xlabel=L"t",ylabel=L"r",label=:none)
	p4=scatter!(p4,[t[n]],[r[n]],label=:none)
	
	p5=plot(t,a_s',xlabel=L"t",ylabel=L"a_s",label=:none)
	p5=scatter!(p5,[t[n]],[a_s[n]],label=:none)
	
	p6=plot(t,a_p',xlabel=L"t",ylabel=L"a_p",label=:none)
	p6=scatter!(p6,[t[n]],[a_p[n]],label=:none)	
	
	l = @layout [a{0.5h} ; b c d; e f]

	cost=J(dt,x,y,sp,u,r,a_p,λ)
	Jtext= @sprintf("%5.3f",cost)
	titulo=L"J=%$Jtext"
    figs=plot(p1, p2, p3, p4, p5, p6, layout = l, plot_title=titulo)
	plot!(size=(600,800))
    return figs
end


# ╔═╡ 6d3e390b-7677-4494-a76b-adeaf21f31ba
show_figure(tf,dt,Tf,position_control,[30 60 -63*pi/180],set_point)

# ╔═╡ f607f9ed-6bb3-4eff-b69f-6d32e04b43bc
coef_op=select_coef(λ)

# ╔═╡ c9d19fe0-c247-4eb5-a290-5f9f2b7ea24b
show_figure(tf,dt,Tf,fourier_control,coef_op,set_point)

# ╔═╡ a7f24817-692e-4915-801e-38631feeee48
begin
	Js_ISE=[]
	Js_E=[]
	λs=0.05:0.05:0.95
	current_J_e=current_J_ise=0
	
	for lambda in λs
		coef_op=select_coef(lambda)
		x,y,_,u,_,r,_,_,ap=simulator(x0w,Tf,dt,fourier_control,coef_op,set_point)
		J_ise=J_ISE(dt,x,y,set_point)
		J_e=J_E(dt,ap,u,r)
		push!(Js_ISE,J_ise/Jref_ISE)
		push!(Js_E,J_e/Jref_E)
		if lambda==λ
			current_J_ise=J_ise/Jref_ISE
			current_J_e=J_e/Jref_E
		end
			
	end
	
	p1=plot(Js_E,Js_ISE, xlabel=L"\frac{J_E}{J_{E_{ref}}}", ylabel=L"\frac{J_{ISE}}{J_{ISE_{ref}}}", label="Pareto front")	

	scatter!(p1,[current_J_e],[current_J_ise],label=L"""\lambda=%$λ""")



	p2=plot(λs, Js_E, label=L"J_E/J_{E_{ref}}",xlabel=L"\lambda",color=:blue)
	plot!(p2,λs, Js_ISE, label=L"J_{ISE}/J_{ISE_{ref}}")
	plot!(p2,λs, (1 .- λs).*Js_ISE.+λs.*Js_E, label=L"J")	


	scatter!(p2,[λ], [current_J_e], label="", color=:blue)
	scatter!(p2,[λ], [current_J_ise], label="",color=:red)
	scatter!(p2,[λ], [(1-λ)*current_J_ise + λ*current_J_e], label="",color=:green)

	plot(p1,p2)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[compat]
LaTeXStrings = "~1.4.0"
Plots = "~1.38.16"
PlutoTeachingTools = "~0.4.7"
PlutoUI = "~0.7.51"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.3"
manifest_format = "2.0"
project_hash = "eb3a1ec1fd19ff7589f21787a7ee9db864416bd6"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fde3bf89aead2e723284a8ff9cdf5b551ed700e8"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.5+0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

    [deps.ColorTypes.weakdeps]
    StyledStrings = "f489334b-da3d-4c2e-b8f0-e476e12c162b"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e357641bb3e0638d353c4b29ea0e40ea644066a6"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.3"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "01ba9d15e9eae375dc1eb9589df76b3572acd3f2"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "b7bfd56fa66616138dfe5237da4dc13bbd83c67f"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "27442171f28c952804dede8ff72828a96f2bfc1f"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.72.10"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "025d171a2847f616becc0f84c8dc62fe18f0f6dd"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.72.10+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "5e6fe50ae7f23d171f44e311c2960294aaa0beb5"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.19"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "d1a86724f81bcd184a38fd284ce183ec067d71a0"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "1.0.0"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "0ee181ec08df7d7c911901ea38baf16f755114dc"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "1.0.0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "b3ad4a0255688dcb895a52fafbaae3023b588a90"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.4.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "97bbca976196f2a1eb9607131cb108c69ec3f8a6"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d0205286d9eceadc518742860bf23f703779a3d6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.3+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c9cbeda6aceffc52d8a0017e71db27c7a7c0beaf"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Preferences", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "9f8675a55b37a70aa23177ec110f6e3f4dd68466"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.38.17"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "Latexify", "Markdown", "PlutoUI"]
git-tree-sha1 = "90b41ced6bacd8c01bd05da8aed35c5458891749"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.4.7"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Downloads", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "3ac7038a98ef6977d44adeadc73cc6f596c08109"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.79"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "522f093a29b31a93e34eaea17ba055d850edea28"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.1"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "9297459be9e338e546f5c4bedb59b3b5674da7f1"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.6.2"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "311349fd1c93a31f783f977a71e8b062a57d4101"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.13"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "6258d453843c466d84c17a58732dda5deeb8d3af"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.24.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    ForwardDiffExt = "ForwardDiff"
    InverseFunctionsUnitfulExt = "InverseFunctions"
    PrintfExt = "Printf"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"
    Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "af305cc62419f9bd61b6644d19170a4d258c7967"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.7.0"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "9cce64c0fdd1960b597ba7ecda2950b5ed957438"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.2+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "b5899b25d17bf1889d25906fb9deed5da0c15b3b"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.12+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a4c0ee07ad36bf8bbce1c3bb52d21fb1e0b987fb"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.7+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a5bc75478d323358a90dc36766f3c99ba7feb024"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.6+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "aff463c82a773cb86061bce8d53a0d976854923e"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.5+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "e3150c7400c41e207012b41659591f083f3ef795"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.3+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "6ab498eaf50e0495f89e7a5b582816e2efb95f64"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.54+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╟─b8fd235e-072b-4697-97a3-b537aba895f2
# ╟─5516b798-bf15-4cc3-a190-dca31aea2304
# ╟─b68dcc2e-d48f-4286-9baf-b7d28a17bd2b
# ╟─25a4d818-4cde-46b2-8af3-34f63647558a
# ╟─d6d4ee7b-7a67-4ae6-9d72-4c21837b8182
# ╟─e11b958e-fb9b-4811-9922-089c3e7aacae
# ╠═79bde88e-0788-4d67-975c-5d47b72846ff
# ╠═ae3961d4-be42-4268-b60f-653fe5842adc
# ╠═91139377-20b7-428c-ae62-8be53d5eb685
# ╠═25e7d4b0-51aa-4064-ad7b-07593d75c628
# ╠═68042e13-9df3-484c-afb5-dd7cc549eb0d
# ╟─4da9ca61-abc6-4b59-93d6-8bb8c83d8aa6
# ╟─baec733f-9474-465f-90a5-401c76e0c0fd
# ╠═87cc3dca-6b63-4e90-9b90-8dbf3c9574ed
# ╟─ed5233f6-ff79-4086-9c51-b5e90eff7648
# ╠═31b5d97b-96c7-48a0-87a3-450a8170ef7e
# ╠═49dd9b51-29a3-4851-bd12-8d386a2c6612
# ╟─ab53825c-7bdf-40ed-8126-b24780c84140
# ╠═fa590ceb-e8b4-4e0f-956a-1b0ec9b62d0f
# ╟─cee61b43-23d9-420d-afee-c29c24d1dd08
# ╠═c8acb0b9-5c22-4bf5-b181-c943313f9da1
# ╟─80fa9ca7-a435-4f5e-9d36-d0a3deb2a3b4
# ╠═71d22223-7bec-4409-862b-eca7854368b8
# ╟─e160754f-6894-462f-8550-209f818b485f
# ╠═126f1b41-6b18-4359-9945-ec18daf32eee
# ╟─f819631b-81b4-46c7-ac8c-0ec210a0aed1
# ╠═2dc7ec04-faea-4f6a-99bb-f226b7075cba
# ╟─826badb3-f4b0-4049-87fb-273498117dce
# ╠═5e225716-e396-42ed-a415-99c45a522bfa
# ╠═c8717d42-4423-4d70-b8ea-8c7f48a2184b
# ╟─bbc464ef-852f-4dc3-ae75-f48516327bb2
# ╠═60aaed79-479c-4d4a-bfcf-03ec77b71ad6
# ╟─8829a36c-266e-40a7-9098-dc7ad3c65d64
# ╠═6d3e390b-7677-4494-a76b-adeaf21f31ba
# ╟─5a6ad0a6-9402-4b8c-997f-69ebbb2bf73a
# ╟─17c6b93f-3820-4957-b9a4-055e6092c212
# ╟─8db37db2-116f-4d71-ba01-495d8a584a98
# ╠═b212918a-0325-4691-b464-22ab12aa30cb
# ╠═65efcc7e-d127-4a71-9412-3a774b46b507
# ╠═f607f9ed-6bb3-4eff-b69f-6d32e04b43bc
# ╟─64cb1add-744c-42fc-90ca-ef6817bc2fc0
# ╠═faa44634-49c8-43dc-8600-77b93e6ce15b
# ╠═c9d19fe0-c247-4eb5-a290-5f9f2b7ea24b
# ╠═b959045e-ce5b-45aa-a9f0-540a18c1c2de
# ╟─1038d5ce-8aeb-4d8f-bbb0-0514fc34e7b1
# ╟─a7f24817-692e-4915-801e-38631feeee48
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
