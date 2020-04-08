#=
main:
- Julia version: 1.3.1
- Author: abd_alazez_ahmed
- Date: 2020-04-02
=#
using Plots
pyplot()
thetaRange=range(-pi,stop=pi,length=100)
println("ϵₐ = ")
ep_a = 1200
println("ϵᵦ = ")
ep_b = 600
println("ϵᵪ = ")
ep_c = -700
println("θₐ = ")
theta_a = -60
println("θᵦ = ")
theta_ab = 60
println("θᵪ = ")
theta_ac = 120
ep = inv([sind(theta_ab)^2 0.5*sind(2*theta_ab); sind(theta_ac)^2 0.5*sind(2*theta_ac)])*[ep_b - ep_a*cosd(theta_ab)^2 ; ep_c - ep_a*cosd(theta_ab)^2]
final_ep = [cosd(theta_a)^2 sind(theta_a)^2 sind(2*theta_a);sind(theta_a)^2 cosd(theta_a)^2 -sind(2*theta_a);-0.5*sind(2*theta_a) 0.5*sind(2*theta_a) cosd(2*theta_a)]*[ep_a;ep[1];ep[2]/2]
println(ep)
final_ep[3]*=2
println(final_ep)
E=2e11
v=0.3
sigMat=E/(1-v^2)*[1 v 0;v 1 0;0 0 (1-v)/2]*final_ep
sigma_x = sigMat[1]
sigma_y = sigMat[2]
tau = sigMat[3]
r = sqrt((sigma_x-sigma_y)^2 /4 + tau^2)
plot((sigma_x+sigma_y)/2 .+ r*cos.(thetaRange), r.*sin.(thetaRange),aspect_ratio=1,framestyle = :origin)
plot!([sigma_x,sigma_x],[0,-tau])
hline!([-tau])
plot!([sigma_y,(sigma_x+sigma_y)/2 + 1.1*r],[-tau,0.1*r*(tau)/((sigma_x+sigma_y)/2 + r-sigma_y)])
plot!([sigma_y,(sigma_x+sigma_y)/2],[-tau,-r])
savefig("Stress.svg")
r = sqrt((final_ep[1]-final_ep[2])^2 /4 + final_ep[3]^2)
plot((final_ep[1]+final_ep[2])/2 .+ r*cos.(thetaRange), r.*sin.(thetaRange),aspect_ratio=1,framestyle = :origin)
plot!([final_ep[1],final_ep[1]],[0,-final_ep[3]])
hline!([-final_ep[3]])
plot!([final_ep[2],(final_ep[1]+final_ep[2])/2 + 1.1*r],[-final_ep[3],0.1*r*(final_ep[3])/((final_ep[1]+final_ep[2])/2 + r-final_ep[2])])
plot!([final_ep[2],(final_ep[1]+final_ep[2])/2],[-final_ep[3],-r])
savefig("Strain.svg")
