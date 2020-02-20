function Assignment2_Part2(part)
%% Function for Part 2 of Assignment 2 4700. Adam Heffernan 100977570
% Set Part to 1 for 2A, 2 for 2B, 3 for 2C and 4 for 2D

%% Part 2A
if part == 1
    %Variables
    clear
    clc
    len=30;
    width=20;%Height
    sigma_map=ones(width,len);
    F=zeros(1,len*width);
    G=zeros(len*width,len*width);
    sigma_map(((3/4)*(width)):width,(len/3):(2*(len/3)))=1e-2;
    sigma_map(1:(width/4),(len/3):(2*(len/3)))=1e-2;
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
    c=colorbar;
    c.Label.String = 'V (V)';
    
    figure(5)
    quiver(Ex,Ey)
    title('Electric Field Vector Map')
    xlabel('X Coordinate')
    ylabel('Y Coordinate')
    xlim([1 len])
    ylim([1 width])
    c=colorbar;
    c.Label.String = 'V (V)';
    
    figure(6)
    surf(sigma_map,'EdgeColor','flat')
    title('Conductivity Map')
    xlabel('X Coordinate')
    ylabel('Y Coordinate')
    view(0,90)
    xlim([1 len])
    ylim([1 width])
    c=colorbar;
    c.Label.String = 'V (V)';
    
    figure(7)
    quiver(Ex.*sigma_map,Ey.*sigma_map)
    title('Current Density Map')
    xlabel('X Coordinate')
    ylabel('Y Coordinate')
    xlim([1 len])
    ylim([1 width])
    c=colorbar;
    c.Label.String = 'V (V)';
    
    
    %% Part 2B Multiple different mesh sizes 
    
elseif part == 2
  
    for multiplier=1:5 
        clear cMap
        clear B
        clear G
        
        %Variables
        len=30*multiplier;
        width=20*multiplier;
        F=zeros(1,len*width);
        G=zeros(len*width,len*width);
        %Set conductivity map
        sigma_map=ones(width,len);
        sigma_map(((3/4)*(width)):width,(len/3):(2*(len/3)))=1e-2;
        sigma_map(1:(width/4),(len/3):(2*(len/3)))=1e-2;
        
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;
                
                if(i==1)%left
                    
                    G(updater,updater)=1;
                    F(updater)=1;
                    
                elseif(i==len)%right
                    
                    %G(N,:)=0;
                    G(updater,updater)=1;
                    F(updater)=0;
                    
                elseif(j==1)%bottom
                    Nyu=j+1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ryu+Rxu+Rxd);
                    G(updater,Nyu)=Ryu;
                    G(updater,Nxu)=Rxu;
                    G(updater,Nxd)=Rxd;
                       
                elseif(j==width)%top
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ryd+Rxu+Rxd);
                    G(updater,Nyd)=Ryd;
                    G(updater,Nxu)=Rxu;
                    G(updater,Nxd)=Rxd;
                    
                else%middle
                    Nyu=j+1+(i-1)*width;
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ryu+Ryd+Rxu+Rxd);
                    G(updater,Nyu)=Ryu;
                    G(updater,Nyd)=Ryd;
                    G(updater,Nxu)=Rxu;
                    G(updater,Nxd)=Rxd;
                    
                end
            end
        end
        
        F = F';
        V=G\F;
        V_2D = zeros(width, len);
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;
                V_2D(j,i)=V(updater);
            end
        end
        
        [Ex,Ey]=gradient(-V_2D);
        %Emag=sqrt(Ex.^2+Ey.^2);
        leftContact(multiplier)=sum(sigma_map(:,1).*Ex(:,1));
        rightContact(multiplier)=sum(sigma_map(:,len).*Ex(:,len));
        
    end
    
    figure(8);
    hold on;
    plot(leftContact,'LineStyle','-','Marker','+','Color','blue');
    plot(rightContact,'LineStyle','-.','Marker','o','Color','red');
    legend('Left Contact','Right Contact');
    hold on;
    title('Current v Mesh Size');
    ylabel('Relative Current');
    xlabel('Mesh size multiplier');

    
    %% Part 2C
elseif part ==3
 %Variables 
    len=30;
    width=20;
    F=zeros(1,len*width);
    G=zeros(len*width,len*width);
    leftContact=zeros(1,50);
    rightContact=zeros(1,50);
    
    for varrying_width=1:width
        %Set conductivity map
        sigma_map=ones(width,len);
        sigma_map(1:varrying_width,(len/2):(2*(len/3)))=1e-2;
        sigma_map((width-varrying_width+1):width,(len/3):(2*(len/3)))=1e-2;
        
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;
                
                if(i==1)%left
                    
                    G(updater,updater)=1;
                    F(updater)=1;
                    
                elseif(i==len)%right
                    
                    G(updater,updater)=1;
                    F(updater)=0;
                    
                elseif(j==1)%bottom
                    Nyu=j+1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ryu+Rxu+Rxd);
                    G(updater,Nyu)=Ryu;
                    G(updater,Nxu)=Rxu;
                    G(updater,Nxd)=Rxd;
                      
                elseif(j==width)%top
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ryd+Rxu+Rxd);
                    G(updater,Nyd)=Ryd;
                    G(updater,Nxu)=Rxu;
                    G(updater,Nxd)=Rxd;
                    
                else%middle
                    Nyu=j+1+(i-1)*width;
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ryu+Ryd+Rxu+Rxd);
                    G(updater,Nyu)=Ryu;
                    G(updater,Nyd)=Ryd;
                    G(updater,Nxu)=Rxu;
                    G(updater,Nxd)=Rxd;
                    
                end
            end
        end
        
        V=G\F';
        V_2D = zeros(width,len);
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;            
                V_2D(j,i)=V(updater);        
            end
        end
        
        [Ex,Ey]=gradient(-V_2D);        
        leftContact(varrying_width)=sum(sigma_map(:,1).*Ex(:,1));
        rightContact(varrying_width)=sum(sigma_map(:,len).*Ex(:,len));
        
    end
    figure(9)
    hold on;
    plot(leftContact,'LineStyle','-','Color','blue');
    plot(rightContact,'LineStyle','-.','Color','red');
    legend('Left Contact','Right Contact');
    xlim([0 20])
    hold off;
    title('Current v Varrying Width')
    xlabel('Width of Boxes')
    ylabel('Relative Current')
    
  
    
    for varrying_length=1:len
        %Set conductivity map
        sigma_map=ones(width,len);
        sigma_map(((3/4)*(width)):width,1:varrying_length)=1e-2;
        sigma_map(1:(width/4),(len-varrying_length+1):len)=1e-2;
        
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;
                
                if(i==1)%left
                    
                    G(updater,updater)=1;
                    F(updater)=1;
                    
                elseif(i==len)%right
                    
                    G(updater,updater)=1;
                    F(updater)=0;
                    
                elseif(j==1)%bottom
                    Nyu=j+1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ryu+Rxu+Rxd);
                    G(updater,Nyu)=Ryu;
                    G(updater,Nxu)=Rxu;
                    G(updater,Nxd)=Rxd;
                      
                elseif(j==width)%top
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ryd+Rxu+Rxd);
                    G(updater,Nyd)=Ryd;
                    G(updater,Nxu)=Rxu;
                    G(updater,Nxd)=Rxd;
                    
                else%middle
                    Nyu=j+1+(i-1)*width;
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ryu+Ryd+Rxu+Rxd);
                    G(updater,Nyu)=Ryu;
                    G(updater,Nyd)=Ryd;
                    G(updater,Nxu)=Rxu;
                    G(updater,Nxd)=Rxd;
                    
                end
            end
        end
        
        V=G\F';
        V_2D = zeros(width,len);
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;            
                V_2D(j,i)=V(updater);        
            end
        end
        
        [Ex,Ey]=gradient(-V_2D);        
        leftContact(varrying_width)=sum(sigma_map(:,1).*Ex(:,1));
        rightContact(varrying_width)=sum(sigma_map(:,len).*Ex(:,len));
        
    end
    figure(10)
    hold on;
    plot(leftContact,'LineStyle','-','Color','blue');
    plot(rightContact,'LineStyle','-.','Color','red');
    legend('Left Contact','Right Contact');
    xlim([0 30])
    hold off;
    title('Current v Varrying Length')
    xlabel('Length of boxes')
    ylabel('Relative Current')
    

    %%Part 2D
elseif part ==4
    %Variables 
    len=30;
    width=20;
    F=zeros(1,len*width);
    G=zeros(len*width,len*width);
    leftContact=zeros(1,50);
    rightContact=zeros(1,50);
    
    for multiplier=1:50
        %Set conductivity map
        sigma_map=ones(width,len);
        sigma_map(((3/4)*(width)):width,(len/3):(2*(len/3)))=1e-2*multiplier;
        sigma_map(1:(width/4),(len/3):(2*(len/3)))=1e-2*multiplier;
        
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;
                
                if(i==1)%left
                    
                    G(updater,updater)=1;
                    F(updater)=1;
                    
                elseif(i==len)%right
                    
                    G(updater,updater)=1;
                    F(updater)=0;
                    
                elseif(j==1)%bottom
                    Nyu=j+1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ryu+Rxu+Rxd);
                    G(updater,Nyu)=Ryu;
                    G(updater,Nxu)=Rxu;
                    G(updater,Nxd)=Rxd;
                      
                elseif(j==width)%top
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ryd+Rxu+Rxd);
                    G(updater,Nyd)=Ryd;
                    G(updater,Nxu)=Rxu;
                    G(updater,Nxd)=Rxd;
                    
                else%middle
                    Nyu=j+1+(i-1)*width;
                    Nyd=j-1+(i-1)*width;
                    Nxu=j+(i-1+1)*width;
                    Nxd=j+(i-1-1)*width;
                    
                    Ryu=(sigma_map(j,i)+sigma_map(j+1,i))/2.0;
                    Ryd=(sigma_map(j,i)+sigma_map(j-1,i))/2.0;
                    Rxu=(sigma_map(j,i)+sigma_map(j,i+1))/2.0;
                    Rxd=(sigma_map(j,i)+sigma_map(j,i-1))/2.0;
                    
                    G(updater,updater)=-(Ryu+Ryd+Rxu+Rxd);
                    G(updater,Nyu)=Ryu;
                    G(updater,Nyd)=Ryd;
                    G(updater,Nxu)=Rxu;
                    G(updater,Nxd)=Rxd;
                    
                end
            end
        end
        
        V=G\F';
        V_2D = zeros(width,len);
        for i=1:len
            for j=1:width
                updater=j+(i-1)*width;            
                V_2D(j,i)=V(updater);        
            end
        end
        
        [Ex,Ey]=gradient(-V_2D);        
        leftContact(multiplier)=sum(sigma_map(:,1).*Ex(:,1));
        rightContact(multiplier)=sum(sigma_map(:,len).*Ex(:,len));
        
    end
    figure(11)
    hold on;
    plot(1e-2*(1:50),leftContact,'LineStyle','-','Color','blue');
    plot(1e-2*(1:50),rightContact,'LineStyle','-.','Color','red');
    legend('Left Contact','Right Contact');
    hold off;
    title('Current v Conductivity')
    xlabel('Conductivity')
    ylabel('Relative Current')
    
    
end