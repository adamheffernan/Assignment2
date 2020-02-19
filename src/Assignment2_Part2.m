function Assignment2_Part2(part)
%% Function for Part 2 of Assignment 2 4700. Adam Heffernan 100977570
% Set Part to 1 for 2A, 2 for 2B, 3 for 2C and 4 for 2D 

%Variables
len=30;
width=20;%Height
sigma_map=ones(width,len);
F=zeros(1,len*width);
G=zeros(len*width,len*width);
sigma_map(((3/4)*(width)):width,(len/3):(2*(len/3)))=1e-2;
sigma_map(1:(width/4),(len/3):(2*(len/3)))=1e-2;

%% Part 2A
if part == 1 
for i=1:len
    for j=1:width
        Updater=j+(i-1)*width;
        
        if(i==1)%left
           
            G(Updater,Updater)=1;
            F(Updater)=1;
            
        elseif(i==len)%right
            
            G(Updater,Updater)=1;
            F(Updater)=0;
            
        elseif(j==1)%bottom
 
            Nyu=j+1+(i-1)*width;
            Nxu=j+(i)*width;
            Nxd=j+(i-1-1)*width;

            Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
            Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
            Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
            
            
            G(Updater,Updater)=-(Ryu+Rxu+Rxd);
            G(Updater,Nyu)=Ryu;
            G(Updater,Nxu)=Rxu;
            G(Updater,Nxd)=Rxd;
            
            
        elseif(j==width)%top

            Nyd=j-1+(i-1)*width;
            Nxu=j+(i)*width;
            Nxd=j+(i-1-1)*width;

            Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
            Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
            Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
            
            
            G(Updater,Updater)=-(Ryd+Rxu+Rxd);
            G(Updater,Nyd)=Ryd;
            G(Updater,Nxu)=Rxu;
            G(Updater,Nxd)=Rxd;
            
        else%middle
            Nyu=j+1+(i-1)*width;
            Nyd=j-1+(i-1)*width;
            Nxu=j+(i)*width;
            Nxd=j+(i-1-1)*width;

            Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
            Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
            Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
            Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
            
            
            G(Updater,Updater)=-(Ryu+Ryd+Rxu+Rxd);
            G(Updater,Nyu)=Ryu;
            G(Updater,Nyd)=Ryd;
            G(Updater,Nxu)=Rxu;
            G(Updater,Nxd)=Rxd;
            
        end        
    end
end

F=F';

V=G\F;
V_2D = zeros(width,len);
for i=1:len
    for j=1:width
       Updater=j+(i-1)*width;
        V_2D(j,i)=V(Updater);
    end
end    

[Ex,Ey] = gradient(-V_2D);


figure(4)
contourf(V_2D,30)
title('Voltage Contour Map')
xlabel('X Coordinate')
ylabel('Y Coordinate')
cb3=colorbar;
cb3.Label.String = 'V (V)';

figure(5)
quiver(Ex,Ey)
title('Electric Field Vector Map')
xlabel('X Coordinate')
ylabel('Y Coordinate')
xlim([1 len])
ylim([1 width])
cb3=colorbar;
cb3.Label.String = 'V (V)';

figure(6)
surf(sigma_map,'EdgeColor','flat')
title('Conductivity Map')
xlabel('X Coordinate')
ylabel('Y Coordinate')
view(0,90)
xlim([1 len])
ylim([1 width])
cb3=colorbar;
cb3.Label.String = 'V (V)';

figure(7)
quiver(Ex.*sigma_map,Ey.*sigma_map)
title('Current Density Map')
xlabel('X Coordinate')
ylabel('Y Coordinate')
xlim([1 len])
ylim([1 width])
cb3=colorbar;
cb3.Label.String = 'V (V)';

% figure(19)
% surf(E_Magnitude)
% title('Electric Field Magnitude Map')
% xlabel('X Coordinate')
% ylabel('Y Coordinate')
% view(0,90)

%% Part 2B     

elseif part == 2 

    
for i=1:len
    for j=1:width
        Updater=j+(i-1)*width;
        
        if(i==1)%left
           
            G(Updater,Updater)=1;
            F(Updater)=1;
            
        elseif(i==len)%right
            
            G(Updater,Updater)=1;
            F(Updater)=0;
            
        elseif(j==1)%bottom
 
            Nyu=j+1+(i-1)*width;
            Nxu=j+(i)*width;
            Nxd=j+(i-1-1)*width;

            Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
            Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
            Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
            
            
            G(Updater,Updater)=-(Ryu+Rxu+Rxd);
            G(Updater,Nyu)=Ryu;
            G(Updater,Nxu)=Rxu;
            G(Updater,Nxd)=Rxd;
            
            
        elseif(j==width)%top

            Nyd=j-1+(i-1)*width;
            Nxu=j+(i)*width;
            Nxd=j+(i-1-1)*width;

            Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
            Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
            Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
            
            
            G(Updater,Updater)=-(Ryd+Rxu+Rxd);
            G(Updater,Nyd)=Ryd;
            G(Updater,Nxu)=Rxu;
            G(Updater,Nxd)=Rxd;
            
        else%middle
            Nyu=j+1+(i-1)*width;
            Nyd=j-1+(i-1)*width;
            Nxu=j+(i)*width;
            Nxd=j+(i-1-1)*width;

            Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
            Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
            Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
            Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
            
            
            G(Updater,Updater)=-(Ryu+Ryd+Rxu+Rxd);
            G(Updater,Nyu)=Ryu;
            G(Updater,Nyd)=Ryd;
            G(Updater,Nxu)=Rxu;
            G(Updater,Nxd)=Rxd;
            
        end        
    end
end

F=F';

V=G\F;
V_2D = zeros(width,len);
for i=1:len
    for j=1:width
       Updater=j+(i-1)*width;
        V_2D(j,i)=V(Updater);
    end
end    

[Ex,Ey] = gradient(-V_2D);
E_Magnitude = sqrt((Ex.^2)+(Ey.^2));    
%% Part 2C  
elseif part ==3
for i=1:len
    for j=1:width
        Updater=j+(i-1)*width;
        
        if(i==1)%left
           
            G(Updater,Updater)=1;
            F(Updater)=1;
            
        elseif(i==len)%right
            
            G(Updater,Updater)=1;
            F(Updater)=0;
            
        elseif(j==1)%bottom
 
            Nyu=j+1+(i-1)*width;
            Nxu=j+(i)*width;
            Nxd=j+(i-1-1)*width;

            Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
            Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
            Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
            
            
            G(Updater,Updater)=-(Ryu+Rxu+Rxd);
            G(Updater,Nyu)=Ryu;
            G(Updater,Nxu)=Rxu;
            G(Updater,Nxd)=Rxd;
            
            
        elseif(j==width)%top

            Nyd=j-1+(i-1)*width;
            Nxu=j+(i)*width;
            Nxd=j+(i-1-1)*width;

            Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
            Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
            Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
            
            
            G(Updater,Updater)=-(Ryd+Rxu+Rxd);
            G(Updater,Nyd)=Ryd;
            G(Updater,Nxu)=Rxu;
            G(Updater,Nxd)=Rxd;
            
        else%middle
            Nyu=j+1+(i-1)*width;
            Nyd=j-1+(i-1)*width;
            Nxu=j+(i)*width;
            Nxd=j+(i-1-1)*width;

            Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
            Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
            Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
            Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
            
            
            G(Updater,Updater)=-(Ryu+Ryd+Rxu+Rxd);
            G(Updater,Nyu)=Ryu;
            G(Updater,Nyd)=Ryd;
            G(Updater,Nxu)=Rxu;
            G(Updater,Nxd)=Rxd;
            
        end        
    end
end

F=F';

V=G\F;
V_2D = zeros(width,len);
for i=1:len
    for j=1:width
       Updater=j+(i-1)*width;
        V_2D(j,i)=V(Updater);
    end
end    

[Ex,Ey] = gradient(-V_2D);
E_Magnitude = sqrt((Ex.^2)+(Ey.^2));       
%%Part 2D 
elseif part ==4
for i=1:len
    for j=1:width
        Updater=j+(i-1)*width;
        
        if(i==1)%left
           
            G(Updater,Updater)=1;
            F(Updater)=1;
            
        elseif(i==len)%right
            
            G(Updater,Updater)=1;
            F(Updater)=0;
            
        elseif(j==1)%bottom
 
            Nyu=j+1+(i-1)*width;
            Nxu=j+(i)*width;
            Nxd=j+(i-1-1)*width;

            Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
            Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
            Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
            
            
            G(Updater,Updater)=-(Ryu+Rxu+Rxd);
            G(Updater,Nyu)=Ryu;
            G(Updater,Nxu)=Rxu;
            G(Updater,Nxd)=Rxd;
            
            
        elseif(j==width)%top

            Nyd=j-1+(i-1)*width;
            Nxu=j+(i)*width;
            Nxd=j+(i-1-1)*width;

            Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
            Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
            Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
            
            
            G(Updater,Updater)=-(Ryd+Rxu+Rxd);
            G(Updater,Nyd)=Ryd;
            G(Updater,Nxu)=Rxu;
            G(Updater,Nxd)=Rxd;
            
        else%middle
            Nyu=j+1+(i-1)*width;
            Nyd=j-1+(i-1)*width;
            Nxu=j+(i)*width;
            Nxd=j+(i-1-1)*width;

            Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
            Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
            Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
            Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
            
            
            G(Updater,Updater)=-(Ryu+Ryd+Rxu+Rxd);
            G(Updater,Nyu)=Ryu;
            G(Updater,Nyd)=Ryd;
            G(Updater,Nxu)=Rxu;
            G(Updater,Nxd)=Rxd;
            
        end        
    end
end

F=F';

V=G\F;
V_2D = zeros(width,len);
for i=1:len
    for j=1:width
       Updater=j+(i-1)*width;
        V_2D(j,i)=V(Updater);
    end
end

[Ex,Ey] = gradient(-V_2D);
E_Magnitude = sqrt((Ex.^2)+(Ey.^2));       
else 
    return
end
end