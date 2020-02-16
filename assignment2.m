%% Assignment 2
% By: Anton Campbell, 100975168
%% Question 1 - Simple Cases
% Two simple cases of 2-D electrostatic potential were created. The y
% dimension length was set to be $W=60$. The x dimension length was
% set to be $L=50$. The conductivity was set at $\sigma=1$ everywhere.
%
% For case 1, the top and bottom we set as free by setting 
% $\frac{dV}{dy}=0$. The left boundary ($x=0$) was set to $V=V_0=1V$.
% The right boundary ($x=L$) was set to $V=0$.
% 
% For case 2, the top and bottom were set to $V=0$. The left and right were
% set to $V=V_0=1V$. The analytical solution was computed using:
%
% $$V\left ( x,y \right
% )=\frac{4V_{0}}{\pi}\sum_{n=1,3,5,...}^{\infty}\frac{1}{n}\frac{cosh\left
% ( \frac{n\pi x}{a} \right )}{cosh\left ( \frac{n\pi b}{a} \right )}sin\left ( \frac{n\pi y}{a} \right )$$
%
% where $a=W$ and $b=L/2$.






%%
% 
% *(a)*
box_potential_1a()
%%
%
% Figure 1 shows the electrostatic potential for case 1 that was solved
% using the finite difference method. The voltage changes linearly between
% the left boundary and the right boundary. The voltage is uniform in the
% y direction.
% 


%%
% 
% *(b)*
box_potential_1b()
%%
% 
% Figure 2 shows the electrostatic potential for case 2 that was solved
% using the finite difference method. Figure 3 shows the electrostatic
% potential for case 2 that was solved using the analytical series
% solution. The voltage map for both solutions has a saddle shape.
% As seen, both solutions are nearly identical. The corners of the
% analytical solution are slightly low.
%
% The analytical solution was performed for a sum going up to $n=115$. At
% this point the analytical solution is approaching the numerical
% solution. The numerical solution using using the finite difference
% method requires more memory. The numerical solution gives a distinct
% solution for the given mesh spacing.
% The analytical solution can require iteration.
%




%% Question 2 - Bottle Neck
% Case 1 was from Question 1 was modified by adding a "bottle-neck". The y
% dimension length was set to be $W=100$. The x dimension length was
% set to be $L=200$.
% Two rectangular "boxes" with low conductivity of $\sigma =10^{-2}$ were
% created with heights of $W_b =40$ and widths of $L_b =40$.
% The conductivity was set at $\sigma=1$ outside the boxes. 
% 

%%
% 
% *(a)*
box_potential_2a()
%%
%
% The relative current through the left contact is $0.23853$ which is the
% same as the current through the right contact.
%
% Figure 4 shows the conductivity map. The conductivity is $\sigma
% =10^{-2}$ inside the "boxes" and $\sigma=1$ outside the boxes.
%
% Figure 5 shows the voltage map with equipotential lines. The
% equipotential lines have closer spacing at the 
% "bottle-neck". Thus, the greatest voltage drop is across the 
% "bottle-neck". The low conductivity near the "bottle-neck" equates to a 
% higher resistance. Since current should be conserved, voltage drop 
% is higher for higher resistance according to the equation $V=IR$.
%
% Figure 6 shows a map of electric field vectors. The vectors are all
% pointing to the right. The magnitude of the electric field is much higher
% inside the boxes because the voltage change is much quicker as seen in 
% Figure 5. 
%
% Figure 7 shows a map of the current density vectors. The vectors are all
% pointing to the right. The current density
% was computed as $\vec{J}=\sigma \vec{E}$. The current density is much
% lower inside the "boxes" which makes sense given their low conductivity.

%%
% 
% *(b)*
box_potential_2b()
%%
% 
% Figure 8 shows the relative current on both contacts as a function of the
% mesh density. The number of mesh divisions in both the x and y directions 
% was changed using multiplying factor. The relative sizes of the
% "bottle-neck" was kept constant.
%
% As the mesh density was increased, the current exponetially approached a
% value, which is the ideal current. As mesh density increases, the
% simulation more accurately represents a real situation. 

%%
% 
% *(c)*
box_potential_2c()
%%
% 
% Figure 9 shows the current as a function of the box heights, $W_b$. 
% As the box heights increase, the width of the higher conductive path
% through the bottle neck decreases. This causes the resistance of the 
% "bottle-neck" to increase. Since the voltages on the left and right
% contacts are fixed, an increased resistance decreases the current
% according to $V=IR$. Thus, increasing the box heights causes the relative
% current to decrease.
% 
% Figure 10 shows the current as a function of the box widths, $L_b$. 
% As the box widths increase, the length of the "bottle-neck" path
% increases. This causes the resistance of the 
% "bottle-neck" to increase. Since the voltages on the left and right
% contacts are fixed, an increased resistance decreases the current
% according to $V=IR$. Thus, increasing the box widths causes the relative
% current to decrease.

%%
% 
% *(d)*
box_potential_2d()
%%
% 
% Figure 11 shows the current as a function of the conductivity in the
% "boxes". As the conductivity in the boxes increases, the overall
% resistance decreases. Since the voltages on the left and right
% contacts are fixed, an decreased resistance increases the current
% according to $V=IR$. Thus, increasing the box conductivity causes the
% relative current to increase.
% 